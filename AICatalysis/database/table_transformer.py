import re
from collections import defaultdict
from pathlib import Path

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
            columns = head.split(",")

            # obtain table footnote
            foot_start = -1
            for index, line in enumerate(table):
                if line.startswith('[a]'):
                    foot_start = index
                    break
            footnotes = self._parse_footnote(table[foot_start:-1])

            body = table[2:foot_start]
            pass

    @staticmethod
    def _parse_footnote(lines):
        footnotes = []
        tokens = get_tokens(lines)

        for line in tokens:
            line_split_index = []
            for index, token in enumerate(line):
                if token.string == '[' and line[index + 2].string == ']':
                    line_split_index.append(token.start[1])

            # split every footnote into list, e.g., ['[a]xxxx', '[b]xxxx', ...]
            line_split = [line[1].line[s:e] for s, e in zip(line_split_index[:-1], line_split_index[1:])]
            ReaCon = defaultdict(list)
            if 'condition' in line_split[0]:
                bcon = re.split(r'[:,]', line_split[0])
                print("Start analyse the `base reaction condition`...")
                for item in bcon:
                    if TransMetal.is_or_not(item):
                        species = TransMetal(item)
                        species.parse()
                        ReaCon['TM'].append((species.formula, species.content))
                        print("TM: " + species.formula, species.content)
                    elif not TransMetal.is_or_not(item) and Metal.is_or_not(item):
                        species = Metal(item)
                        species.parse()
                        ReaCon['M'].append((species.formula, species.content))
                        print("M: " + species.formula, species.content)
                    elif Ligand.is_or_not(item):
                        species = Ligand(item)
                        species.parse()
                        ReaCon['L'].append((species.formula, species.content))
                        print("L: " + species.formula, species.content)
                    elif Solvent.is_or_not(item):
                        species = Solvent(item)
                        species.parse()
                        ReaCon['Sol'].append((species.formula, species.content))
                        print("Sol: " + species.formula, species.content)
                    elif Gas.is_or_not(item):
                        species = Gas(item)
                        ReaCon['Gas'].append(species.name)
                        print("Gas: " + species.name)
                    elif Time.is_or_not(item):
                        species = Time(item)
                        ReaCon['Time'].append(species.name)
                        print("Time: " + species.name)
                    elif Temperature.is_or_not(item):
                        species = Temperature(item)
                        ReaCon['Temp'].append(species.name)
                        print("Temperature: " + species.name)
                    elif "reaction" in item.lower():
                        continue
                    else:
                        print(f"Can't recognize `{item}`")
                print()
            else:
                print("Can't find `base reaction condition` in footnotes!")

            footnotes.append(ReaCon)

        return footnotes


if __name__ == '__main__':
    csv_dir = "."
    files = [file for file in Path(csv_dir).iterdir() if file.suffix == ".csv"]
    csvreader = CSVReader(files[0])
    results = csvreader.parse()
    pass
