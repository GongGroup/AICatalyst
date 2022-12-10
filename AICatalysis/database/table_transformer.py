import re
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

    def _parse_footnote(self, lines):
        footnotes = []
        tokens = get_tokens(lines)

        for line in tokens:
            line_split_index = []
            for index, token in enumerate(line):
                if token.string == '[' and line[index + 2].string == ']':
                    line_split_index.append(token.start[1])

            # split every footnote into list, e.g., ['[a]xxxx', '[b]xxxx', ...]
            line_split = [line[1].line[s:e] for s, e in zip(line_split_index[:-1], line_split_index[1:])]
            if 'condition' in line_split[0]:
                bcon = re.split(r':|,', line_split[0])
                print("Start analyse the `base reaction condition`...")
                for item in bcon:
                    if TransMetal.is_or_not(item):
                        print("TM: " + TransMetal(item).name)
                    elif not TransMetal.is_or_not(item) and Metal.is_or_not(item):
                        print("M: " + Metal(item).name)
                    elif Ligand.is_or_not(item):
                        print("L: " + Ligand(item).name)
                    elif Solvent.is_or_not(item):
                        print("Sol: " + Solvent(item).name)
                    elif Gas.is_or_not(item):
                        print("Gas: " + Gas(item).name)
                    elif Time.is_or_not(item):
                        print("Time: " + Time(item).name)
                    elif Temperature.is_or_not(item):
                        print("Temperature: " + Temperature(item).name)
                    elif "reaction" in item.lower():
                        continue
                    else:
                        print(f"Can't recognize `{item}`")
                print()
            else:
                print("Can't find `base reaction condition` in footnotes!")

        return footnotes


if __name__ == '__main__':
    csv_dir = "."
    files = [file for file in Path(csv_dir).iterdir() if file.suffix == ".csv"]
    csvreader = CSVReader(files[0])
    results = csvreader.parse()
    pass
