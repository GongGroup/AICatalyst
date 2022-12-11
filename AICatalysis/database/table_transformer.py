import re
from collections import defaultdict
from pathlib import Path

import numpy as np

from AICatalysis.common.error import ParseError
from AICatalysis.common.species import TransMetal, Solvent, Metal, Time, Temperature, Gas, Ligand
from AICatalysis.common.utils import get_tokens


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

            # locate table footnote
            foot_start = -1
            for index, line in enumerate(table):
                if line.startswith('[a]'):
                    foot_start = index
                    break

            # parse body && footnotes
            body = self._parse_body(table[2:foot_start], checked_AllCol)  # type -> list(dict)
            footnotes = self._parse_footnote(table[foot_start:-1])
            records = self._merge_bf(body, footnotes[0])
            pass

    @staticmethod
    def _parse_head(columns):
        AllCol = {'metal': None,
                  'ligand': None,
                  'gas': None,
                  'solvent': None,
                  'reagent': None,
                  'time': None,
                  'temperature': None,
                  'yield': None}

        for key in AllCol.keys():
            for index, col in enumerate(columns):
                if key in col.lower():
                    AllCol[key] = index

        return AllCol

    @staticmethod
    def _parse_body(lines, AllCol):
        lines_np = np.array([line.split(",") for line in lines])
        search_tuple = [(key, index) for key, index in AllCol.items() if index is not None]
        body = []
        for line in lines_np:
            record = {key: line[index] for key, index in search_tuple}
            body.append(record)

        return body

    @staticmethod
    def _parse_footnote(lines):

        def sub_parse_species(class_name, item):
            species = class_name(item)
            species.parse()
            return species.formula, species.content

        footnotes = []
        tokens = get_tokens(lines)

        for line in tokens:
            line_split_index = []
            for index, token in enumerate(line):
                if token.string == '[' and line[index + 2].string == ']':
                    line_split_index.append(token.start[1])

            # split every footnote into list, e.g., ['[a]xxxx', '[b]xxxx', ...]
            multi_foots = [line[1].line[s:e] for s, e in zip(line_split_index[:-1], line_split_index[1:])]
            ReaCon = defaultdict(list)
            if 'condition' in multi_foots[0]:
                bcon = re.split(r'[:,]', multi_foots[0])
                print("Start analyse the `base reaction condition`...")
                for item in bcon:
                    if TransMetal.is_or_not(item):
                        ReaCon['metal'].append(sub_parse_species(TransMetal, item))
                    elif not TransMetal.is_or_not(item) and Metal.is_or_not(item):
                        ReaCon['reagent'].append(sub_parse_species(Metal, item))
                    elif Ligand.is_or_not(item):
                        ReaCon['ligand'].append(sub_parse_species(Ligand, item))
                    elif Solvent.is_or_not(item):
                        ReaCon['solvent'].append(sub_parse_species(Solvent, item))
                    elif Gas.is_or_not(item):
                        species = Gas(item)
                        ReaCon['gas'].append(species.name)
                    elif Time.is_or_not(item):
                        species = Time(item)
                        ReaCon['time'].append(species.name)
                    elif Temperature.is_or_not(item):
                        species = Temperature(item)
                        ReaCon['temperature'].append(species.name)
                    elif "reaction" in item.lower():
                        continue
                    else:
                        print(f"Can't recognize `{item}`")
                print()
            else:
                print("Can't find `base reaction condition` in footnotes!")

            footnotes.append(ReaCon)

        return footnotes

    @staticmethod
    def _merge_bf(body, footnote):
        features = ['metal', 'ligand', 'gas', 'solvent', 'reagent', 'time', 'temperature', 'yield']
        records = []
        for item in body:
            temp_record = {}
            for fea in features:
                if item.get(fea, None) is not None:
                    temp_record[fea] = item[fea]
                elif footnote.get(fea, None) is not None:
                    temp_record[fea] = footnote[fea]
            records.append(temp_record)

        return records


if __name__ == '__main__':
    csv_dir = "."
    files = [file for file in Path(csv_dir).iterdir() if file.suffix == ".csv"]
    csvreader = CSVReader(files[0])
    results = csvreader.parse()
    pass
