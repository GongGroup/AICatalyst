import copy
import re
import string
from collections import defaultdict
from pathlib import Path

import numpy as np

from AICatalysis.common.error import ParseError
from AICatalysis.common.species import TransMetal, Solvent, Reagent, Time, Temperature, Gas, Ligand, Base, Additive, \
    Oxidant
from AICatalysis.common.utils import get_tokens, flatten


class FileIO(object):
    def __init__(self, file):
        self.file = file

    @property
    def strings(self):
        with open(self.file, "r", encoding="utf-8") as f:
            content = f.read()

        _strings = [''.join(x for x in line if x.isprintable()) for line in content.splitlines()]

        return _strings


class CSVReader(FileIO):
    def __init__(self, file):
        super(CSVReader, self).__init__(file)

    def parse(self):

        # split tables
        start, end = [], []
        for index, line in enumerate(self.strings):
            if line.startswith('Table'):
                start.append(index)
            elif line.startswith('http'):
                end.append(index)

        for s, e in zip(start, end):
            table = self.strings[s:e + 1]
            caption = table[0]
            doi = table[-1]

            # obtain table head
            if table[1].lower().startswith('entry'):
                head = table[1]
            else:
                raise ParseError("The table has not `entry` field, please check!!")
            checked_AllCol = self._parse_head(head.split(","))
            if checked_AllCol['yield'] is None:
                print("Warning: `yield` feature not exist, continue")
                continue

            # locate table footnote
            foot_start = -1
            for index, line in enumerate(table):
                if line.startswith('[a]') or line.startswith('a'):
                    foot_start = index
                    break

            # parse body && footnotes
            body = self._parse_body(table[2:foot_start], checked_AllCol)  # type -> list(dict)
            multi_foots, base_i, footnotes = self._parse_footnote(table[foot_start:-1])
            records = self._merge_bf(body, base_i, footnotes)
            print()
            pass

    @staticmethod
    def _parse_head(columns):
        AllCol = {'entry': None,
                  'metal': None,
                  'ligand': None,
                  'gas': None,
                  'solvent': None,
                  'reagent': None,
                  'time': None,
                  'temperature': None,
                  'yield': None,
                  'base': None,
                  'additive': None,
                  'oxidant': None}

        for key in AllCol.keys():
            for index, col in enumerate(columns):
                if key in col.lower():
                    AllCol[key] = index
                if "temp" in col and key == 'temperature':
                    AllCol[key] = index
                if ("catalyst" in col or "Palladium" == col or "[Pd]" == col) and key == 'metal':
                    AllCol[key] = index

        return AllCol

    @staticmethod
    def _parse_body(lines, AllCol):

        def unify_value(value):
            patten = re.compile("([0-9]+)([a-z])")
            if patten.search(value) is not None:
                symbol = patten.search(value).groups()
                return symbol[0] + "[" + symbol[1] + "]"
            else:
                return value

        lines_np = np.array([line.split(",") for line in lines])
        search_tuple = [(key, index) for key, index in AllCol.items() if index is not None]
        body = []
        for line in lines_np:
            record = {key: unify_value(line[index]) for key, index in search_tuple if index < len(line)}
            body.append(record)

        return body

    @staticmethod
    def _parse_footnote(lines):

        def sub_parse_species(class_name, item):
            species = class_name(item)
            species.parse()
            return species.formula, species.content

        def unify_ref(lines):
            unify_lines = []
            for l in lines:
                if re.search(r'\[[a-z]+]', l) is not None:
                    unify_lines.append(l)
                else:
                    patten1 = re.compile(r'^([a-z]{1})\s')
                    patten2 = re.compile(r'\s([a-z]{1})\s')
                    l = patten1.sub(r'[\1] ', l)
                    l = patten2.sub(r'[\1] ', l)
                    unify_lines.append(l)
            return unify_lines

        base_i, footnotes, multi_foots = None, [], []
        join_lines = [' '.join(lines)]
        unify_lines = unify_ref(join_lines)
        tokens = get_tokens(unify_lines)

        for line in tokens:
            line_split_index = []
            for index, token in enumerate(line):
                if token.string == '[' and line[index + 2].string == ']':
                    line_split_index.append(token.start[1])
            else:
                line_split_index.append(line[index - 2].end[1])
            # split every footnote into list, e.g., ['[a]xxxx', '[b]xxxx', ...]
            multi_foots = [line[1].line[s:e] for s, e in zip(line_split_index[:-1], line_split_index[1:])]
            print("\n".join(multi_foots))
            for index, single_foot in enumerate(multi_foots):
                ReaCon = defaultdict(list)
                single_foot = re.sub(r'\[[a-z]+\]', r'', single_foot)  # remove [a] in the first
                cond = re.split(r': |,|;', single_foot)
                if 'condition' in single_foot.lower():
                    base_i = index
                for item in cond:
                    if TransMetal.is_or_not(item):
                        ReaCon['metal'].append(sub_parse_species(TransMetal, item))
                    elif Time.is_or_not(item):
                        species = Time(item)
                        ReaCon['time'].append(species.name)
                    elif Ligand.is_or_not(item):
                        ReaCon['ligand'].append(sub_parse_species(Ligand, item))
                    elif Solvent.is_or_not(item):
                        ReaCon['solvent'].append(sub_parse_species(Solvent, item))
                    elif not TransMetal.is_or_not(item) and Reagent.is_or_not(item):
                        ReaCon['reagent'].append(sub_parse_species(Reagent, item))
                    elif Gas.is_or_not(item):
                        species = Gas(item)
                        ReaCon['gas'].append(species.name)
                    elif Base.is_or_not(item):
                        ReaCon['base'].append(sub_parse_species(Base, item))
                    elif Oxidant.is_or_not(item):
                        ReaCon['oxidant'].append(sub_parse_species(Oxidant, item))
                    elif Temperature.is_or_not(item):
                        species = Temperature(item)
                        ReaCon['temperature'].append(species.name)
                    elif Additive.is_or_not(item):
                        ReaCon['additive'].append(sub_parse_species(Additive, item))
                    elif "reaction" in item.lower():
                        continue
                    else:
                        # print(f"Can't recognize `{item}`")
                        pass
                footnotes.append(ReaCon)

        return multi_foots, base_i, footnotes

    @staticmethod
    def _merge_bf(body, base_i, footnotes):

        def parse_ref(ref: str):
            symbol = re.search("\[([a-z]+)]", ref).groups()[0]
            if symbol in string.ascii_lowercase:
                return string.ascii_lowercase.index(symbol)

        base_cond = footnotes[base_i] if base_i is not None else None
        features = ['metal', 'ligand', 'gas', 'solvent', 'reagent', 'time', 'temperature', 'yield', 'base', 'additive',
                    'oxidant']
        species_class = {'ligand': Ligand, 'solvent': Solvent, 'metal': TransMetal, 'base': Base, 'additive': Additive,
                         'oxidant': Oxidant}
        records = []
        for item in body:
            patten = re.compile("(\[[a-z]+])")
            indicator = [patten.search(value).groups() if patten.search(value) is not None else None
                         for value in item.values()]
            indicator = [item for item in flatten(indicator) if item is not None]
            symbols = [parse_ref(ii) for ii in indicator]
            other_cond = [footnotes[ii] for ii in symbols] if len(symbols) else []
            temp_record = {}
            for fea in features:
                if base_cond is not None and base_cond.get(fea, None) is not None:  # base condition
                    temp_record[fea] = copy.deepcopy(base_cond[fea])  # memory view may fail
                if item.get(fea, None) is not None:  # sub base-condition with body content
                    if temp_record.get(fea, None) is not None and fea in ['metal', 'ligand', 'solvent', 'base',
                                                                          'additive', 'oxidant']:
                        ll = species_class[fea](item[fea])  # body content
                        ll.parse()
                        if len(temp_record[fea]) == 1:
                            if ll.formula is not None:
                                temp_record[fea][0] = (ll.formula, temp_record[fea][0][1])
                            if ll.content is not None:
                                temp_record[fea][0] = (temp_record[fea][0][0], ll.content)
                        else:
                            print("Multi-Match found, please check <merge_bf>")
                    else:
                        temp_record[fea] = item[fea]
            else:  # expand footnotes-ref
                if len(other_cond):
                    for oc in other_cond:
                        for key in oc.keys():
                            if key == "ligand":
                                if oc[key][0][0] == "Ligand":
                                    temp_record[key][0] = (temp_record[key][0][0], oc[key][0][1])
                                else:
                                    temp_record[key][0] = oc[key][0]
                            else:
                                temp_record[key] = oc[key]
            records.append(temp_record)

        return records


if __name__ == '__main__':
    csv_dir = "."
    files = [file for file in Path(csv_dir).iterdir() if file.suffix == ".csv"]
    file = files[7]
    print(file)
    csvreader = CSVReader(file)
    results = csvreader.parse()
    pass

# --*--exclude--*--
# 02215d7edfcd500d46ee0fc005b9422a.csv
# 0558c776679c8f021e1f74c348648d45.csv
# 065f47cd7e8a3482f0a620c2fd3d7a87.csv
