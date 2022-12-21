from collections import defaultdict
from pathlib import Path

from AICatalysis.common.species import ChemFormula
from AICatalysis.database.table_transformer import TableTransformer

Features = ['catalyst', 'reagent', 'acid', 'base', 'additive', 'oxidant', 'ligand', 'solvent', 'gas', 'time', 'yield',
            'temperature']

if __name__ == '__main__':
    csv_dir = "tcsv"
    files = [file for file in Path(csv_dir).iterdir() if file.suffix == ".csv"]
    file = files[20]
    csvreader = TableTransformer(file)
    csvreader.parse()
    records = csvreader.records

    total_merged = []
    # print("catalyst;ligand;species;solvent;gas;temperature;time;yield")
    for table in records:
        for item in table:
            merged_item = defaultdict(list)
            for fea, value in item.items():
                if fea in Features[:6]:
                    if ChemFormula.is_or_not(value[0]):
                        species = ChemFormula(value[0])
                        ions = species.split()
                        if not isinstance(ions, tuple):
                            ions = (ions, '')
                            print(f"{value[0]} can't split")
                        if fea == 'catalyst':
                            merged_item[fea] = (*ions, value[1])
                        else:
                            merged_item['species'].append((*ions, value[1]))
                    else:
                        print(f"{value[0]} is not formula")
                        if fea == 'catalyst':
                            merged_item[fea] = (value[0], '', value[1])
                        else:
                            merged_item['species'].append((value[0], '', value[1]))
                else:
                    merged_item[fea] = value
            # print(merged_item.get('catalyst', None), end=";")
            # print(merged_item.get('ligand', None), end=";")
            # print(merged_item.get('species', None), end=";")
            # print(merged_item.get('solvent', None), end=";")
            # print(merged_item.get('gas', None), end=";")
            # print(merged_item.get('temperature', None), end=";")
            # print(merged_item.get('time', None), end=";")
            # print(merged_item.get('yield', None), end=";")
            # print()
            total_merged.append(merged_item)
    pass
