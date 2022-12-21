import copy
import re
import string
from collections import defaultdict
from pathlib import Path

import numpy as np

from AICatalysis.common.error import ParseError
from AICatalysis.common.species import CarbonylCatalyst, Solvent, Reagent, Time, Temperature, Gas, Ligand, Base, Acid, \
    Additive, Oxidant
from AICatalysis.common.utils import get_tokens, flatten, is_number, is_ratio

Features = ['catalyst', 'reagent', 'acid', 'base', 'additive', 'oxidant', 'ligand', 'solvent', 'gas', 'time', 'yield',
            'temperature']
SpeciesClass = {'catalyst': CarbonylCatalyst, 'ligand': Ligand, 'solvent': Solvent, 'gas': Gas, 'acid': Acid,
                'base': Base, 'reagent': Reagent, 'additive': Additive, 'oxidant': Oxidant}
GlobalExclude = ["reaction", "1a (", "yield"]


class FileIO(object):
    def __init__(self, file):
        self.file = file

    @property
    def strings(self):
        with open(self.file, "r", encoding="utf-8") as f:
            content = f.read()

        _strings = [''.join(x for x in line if x.isprintable()) for line in content.splitlines()]

        return _strings


class TableTransformer(FileIO):
    def __init__(self, file):
        super(TableTransformer, self).__init__(file)
        self.records = []

    def write(self):
        print(";".join(Features))
        for table in self.records:
            for record in table:
                for fea in Features:
                    print(record.get(fea, ''), end="; ")
                print()

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
            # print(doi)

            # obtain table head
            if table[1].lower().startswith('entry'):
                head = table[1]
            else:
                raise ParseError("The table has not `entry` field, please check!!")

            # parse head
            checked_AllCol = self._parse_head(head.split(","))
            if checked_AllCol['yield'] is None:
                # print("Warning: `yield` feature not exist, continue")
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
            self.records.append(records)

    @staticmethod
    def _parse_head(columns):
        AllCol = {fea: (None, None) for fea in Features}
        AllCol.update({'entry': None})

        for key in AllCol.keys():
            for index, col in enumerate(columns):
                unit = re.findall("\((.*?)\)", col)  # unit use (); ref use []
                COD1 = key in col.lower()
                COD2 = key == 'temperature' and ("temp" in col.lower() or "T " in col)
                COD3 = key == 'catalyst' and ("catalyst" in col.lower() or "Palladium" == col or "Pd" in col or
                                              "LTA" in col or "cat." in col.lower())
                COD4 = key == 'acid' and "HCO2H" in col
                COD5 = key == 'gas' and ("PCO [MPa]" == col or "CO" in col)
                COD6 = key == 'reagent' and ("carbonyl" in col)
                COD7 = key == 'ligand' and ('PR3' in col or 'L' == col)
                COD8 = key == 'time' and 't ' in col

                if COD1 or COD2 or COD3 or COD4 or COD5 or COD6 or COD7 or COD8:
                    AllCol[key] = (index, unit[-1]) if len(unit) else (index, '')

        return AllCol

    @staticmethod
    def _parse_body(lines, AllCol):

        def unify_value(value):
            """
            Add [] at left and right sides of footnote ref
            """
            patten = re.compile("(.*[0-9]+)([a-z])")
            if patten.search(value) is not None:
                symbol = patten.search(value).groups()
                return symbol[0] + "[" + symbol[1] + "]"
            else:
                return value

        lines_np = np.array([line.split(",") for line in lines])
        search_tuple = [(key, index, unit) for key, (index, unit) in AllCol.items() if index is not None]
        body = []
        for line in lines_np:
            record = {key: (unify_value(line[index]), unit) for key, index, unit in search_tuple if index < len(line)}
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
            # print("\n".join(multi_foots))
            for index, single_foot in enumerate(multi_foots):
                ReaCon = defaultdict(tuple)
                single_foot = re.sub(r'\[[a-z]+\]', r'', single_foot)  # remove [a] in the first
                cond = re.split(r': |,|;', single_foot)
                if 'condition' in single_foot.lower():
                    base_i = index
                for item in cond:
                    if Time.is_or_not(item) and ReaCon.get('time', None) is None:
                        species = Time(item)
                        ReaCon['time'] = species.name
                        continue
                    if Temperature.is_or_not(item) and ReaCon.get('temperature', None) is None:
                        species = Temperature(item)
                        ReaCon['temperature'] = species.name
                        continue

                    # catch species in footnotes => (formula, content)
                    findit_flag = False
                    for key, class_ in SpeciesClass.items():
                        if class_.is_or_not(item) and ReaCon.get(key, None) is None:
                            ReaCon[key] = (sub_parse_species(class_, item))
                            findit_flag = True
                            break
                    if findit_flag:
                        continue

                    # special indicator, continue
                    for ex in GlobalExclude:
                        if ex in item.lower():
                            break
                    else:
                        # print(f"Warning: Can't recognize `{item}`")
                        pass

                footnotes.append(ReaCon)

        return multi_foots, base_i, footnotes

    @staticmethod
    def _merge_bf(body, base_i, footnotes):

        def parse_ref(ref: str):
            symbol = re.search("\[([a-z]+)]", ref).groups()[0]
            if symbol in string.ascii_lowercase:
                return string.ascii_lowercase.index(symbol)

        def parse_species(item, fea):  # TODO: need merge
            patten1 = re.compile("(.*?)\s\(([0-9]\.?[0-9]?)\)")  # e.g., PPh3 (3)
            patten2 = re.compile("(.*)/([0-9]\.?[0-9]?)")  # e.g., Pd(PPh3)2Cl2/5
            patten3 = re.compile("(.*)\(([0-9]+\.?[0-9]?)\)")  # e.g., AgOAc(1.5)
            patten4 = re.compile("(.*?)\s\((–)\)")  # e.g., PdCl2 (–)
            if is_number(item[fea][0]):  # e.g.,item[fea] = 5.0 => acid (5.0)
                ll = SpeciesClass[fea](f"{fea} ({item[fea][0]} {item[fea][1]})")
            elif is_ratio(item[fea][0]):
                ll = SpeciesClass[fea](f"{fea} ({item[fea][0]})")
            elif patten1.search(item[fea][0]) is not None:
                match = patten1.search(item[fea][0])
                ll = SpeciesClass[fea](f"{match.groups()[0]} ({match.groups()[1] + item[fea][1]})")
            elif patten2.search(item[fea][0]) is not None:
                match = patten2.search(item[fea][0])
                ll = SpeciesClass[fea](f"{match.groups()[0]} ({match.groups()[1] + item[fea][1]})")
            elif patten3.search(item[fea][0]) is not None:
                match = patten3.search(item[fea][0])
                ll = SpeciesClass[fea](f"{match.groups()[0]} ({match.groups()[1] + item[fea][1]})")
            elif patten4.search(item[fea][0]) is not None:
                match = patten4.search(item[fea][0])
                ll = SpeciesClass[fea](f"{match.groups()[0]}")
            elif item[fea][0] == "–":
                ll = SpeciesClass[fea](item[fea][0])
            else:
                ll = SpeciesClass[fea](''.join(item[fea]))  # body content
            ll.parse()

            return ll

        base_cond = footnotes[base_i] if base_i is not None else None

        records = []
        for body_item in body:
            patten = re.compile("(\[[a-z]+])")
            indicator = [patten.findall(value) if patten.search(value) is not None else None
                         for value, unit in body_item.values()]
            indicator = [item for item in flatten(indicator) if item is not None]
            symbols = [parse_ref(ii) for ii in indicator]
            other_cond = [footnotes[ii] for ii in symbols] if len(symbols) else []
            temp_record = {}
            for fea in Features:
                # First: fill with base condition
                if base_cond is not None and base_cond.get(fea, None) is not None:
                    temp_record[fea] = copy.deepcopy(base_cond[fea])  # memory view may fail
                # Second: modify with the body item
                if body_item.get(fea, None) is not None:
                    # fea in body and base-cond
                    if temp_record.get(fea, None) is not None and fea in SpeciesClass.keys():
                        ll = parse_species(body_item, fea)
                        ll.formula = None if ll.formula.strip() == fea else ll.formula
                        if ll.formula is not None:
                            if ll.formula != "–":
                                temp_record[fea] = (ll.formula, temp_record[fea][1])
                            else:  # if formula == "–" => (–, None)
                                temp_record[fea] = (ll.formula, None)
                        if ll.content is not None:
                            temp_record[fea] = (temp_record[fea][0], ll.content)

                    # fea in body but not in base cond
                    else:
                        if fea in ['catalyst', 'oxidant', 'gas', 'base', 'reagent']:
                            ll = parse_species(body_item, fea)
                            temp_record[fea] = (ll.formula, ll.content)
                        else:
                            temp_record[fea] = ''.join(body_item[fea]) if body_item[fea][0] != '–' else body_item[fea][
                                0]
            # Three: Expand with other-conditions
            else:
                if len(other_cond):
                    for oc in other_cond:
                        for fea in oc.keys():
                            # fea in other-cond and temp_record
                            if temp_record.get(fea, None) is not None:
                                if fea == "ligand":
                                    if oc[fea][0].strip().lower() == fea:
                                        temp_record[fea] = (temp_record[fea][0], oc[fea][1])
                                    else:
                                        if temp_record[fea][0] != "–":
                                            temp_record[fea] = oc[fea]
                                elif fea == "additive":
                                    if fea in temp_record.keys():
                                        temp_record[fea] += oc[fea]
                                    else:
                                        temp_record[fea] = oc[fea]
                                else:
                                    if temp_record[fea][0][0] != "–":
                                        temp_record[fea] = oc[fea]
                            # fea in other-cond but not in temp_record
                            else:
                                temp_record[fea] = oc[fea]

            records.append(temp_record)

        return records

# if __name__ == '__main__':
# csv_dir = "tcsv"
# files = [file for file in Path(csv_dir).iterdir() if file.suffix == ".csv"]
# file = files[0]
# print(file)
# csvreader = TableTransformer(file)
# csvreader.parse()
# csvreader.write()
# pass

# --*--exclude--*--
# 02215d7edfcd500d46ee0fc005b9422a.csv
# 0558c776679c8f021e1f74c348648d45.csv
# 065f47cd7e8a3482f0a620c2fd3d7a87.csv
# 155a9f614f863ec1f3fc1b094f7569e5.csv
# 17a6965ceeefdf8e491c794ad96643c6.csv
