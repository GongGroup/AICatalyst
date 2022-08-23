import json
import os
import shutil
from functools import partial
from pathlib import Path

import execjs
import requests
from rdkit import Chem, RDConfig
from rdkit.Chem import FragmentCatalog

from common.fio import JsonIO, ftemp, fcopy
from common.logger import logger

# Net const
OPSINRoot = 'https://opsin.ch.cam.ac.uk/'
headers = {
    "user-agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) "
                  "Chrome/104.0.5112.81 Safari/537.36 Edg/104.0.1293.47"}

# Directory const
ChemDir = Path("../chemical")

# File const
FJSFunc = "func.js"
FChemical = ChemDir / "chemical.json"
FFormula = ChemDir / "formula.json"
FOpsinRecord = ChemDir / "opsin_record.json"
FIChemRecord = ChemDir / "ichem_record.json"
FNameName = ChemDir / "name_name.json"


class OPSINCrawler(object):
    def __init__(self, file=None):
        self.io = JsonIO
        self.encode = partial(execjs.compile(OPSINCrawler.load_js(FJSFunc)).call, "encode")

        self.file = file
        self.chemicals = self.io.read(self.file) if self.file is not None else None

    def get_records(self):
        """
        According to the chemical.json [list] obtain the smiles, inchi, inchi-key and so on
        """

        if not Path(FOpsinRecord).exists():
            self.io.write([], FOpsinRecord)

        records = self.io.read(FOpsinRecord)
        records_name = [item['name'] for item in records]

        for index, chemical in enumerate(self.chemicals):
            if chemical in records_name:
                logger.debug(f"`{chemical}` has been stored in database, continue")
                continue

            url = OPSINRoot + "opsin/" + self.encode(chemical)
            html = requests.get(url, headers=headers)
            if html.status_code == 200:
                content = json.loads(html.content)
                content['name'] = chemical
                records.append(content)
                logger.info(f"`{chemical}` successfully store in database")
            else:
                content = {'name': chemical, 'status': 'FAILURE'}
                records.append(content)
                logger.warning(f"`{chemical}` : Response != 200, please check")

        shutil.copy(FOpsinRecord, fcopy(FOpsinRecord))
        self.io.write(records, ftemp(FOpsinRecord))
        shutil.move(ftemp(FOpsinRecord), FOpsinRecord)

        logger.info("Successfully update the opsin_record.json")

    def update_records(self):
        formatted_func = lambda x: x.replace('<i>', '').replace('</i>', '').replace('<sup>', '').replace('</sup>', '')
        chemicals = [item for item in JsonIO.read(FNameName) if item["new_name"] != "FAILURE"]
        formatted_chemicals = [{'old_name': item['old_name'], 'new_name': formatted_func(item['new_name'])}
                               for item in chemicals]
        formatted_chemicals_name = [item['old_name'] for item in formatted_chemicals]

        records = [item for item in self.io.read(FOpsinRecord) if item not in formatted_chemicals_name]
        for item in formatted_chemicals:
            url = OPSINRoot + "opsin/" + self.encode(item['new_name'])
            html = requests.get(url, headers=headers)
            if html.status_code == 200:
                content = json.loads(html.content)
                content.update({'name': item['old_name'], 'new_name': item['new_name']})
                records.append(content)
                logger.info(f"`{item['old_name']}` successfully store in database")
            else:
                content = {'name': item['old_name'], 'new_name': item['new_name'], 'status': 'FAILURE'}
                records.append(content)
                logger.warning(f"`{item['old_name']}` : Response != 200, please check")

        self.io.write(records, ftemp(FOpsinRecord))
        shutil.move(ftemp(FOpsinRecord), FOpsinRecord)

        logger.info("Successfully update the opsin_record.json")

    def update_formulas(self):
        formulas = {key: value for key, value in self.io.read(FFormula).items() if value != "FAILURE"}

        records = []
        for item in self.io.read(FOpsinRecord):
            if item['name'] in formulas.keys() and item.get('formula', None) is None:
                item['formula'] = formulas[item['name']]
                logger.info(f"update `{item['name']}` formula successfully")
            records.append(item)

        shutil.copy(FOpsinRecord, fcopy(FOpsinRecord))
        self.io.write(records, ftemp(FOpsinRecord))
        shutil.move(ftemp(FOpsinRecord), FOpsinRecord)

        logger.info("Successfully update the opsin_record.json")

    def update_groups(self):
        fname = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
        fparams = FragmentCatalog.FragCatParams(1, 6, fname)

        def get_group(smiles):
            m = Chem.MolFromSmiles(smiles)
            fcat = FragmentCatalog.FragCatalog(fparams)
            fcgen = FragmentCatalog.FragCatGenerator()
            if m is None:
                return None
            fcgen.AddFragsFromMol(m, fcat)

            group_count = fcat.GetNumEntries()
            group_id = [list(fcat.GetEntryFuncGroupIds(i)) for i in range(group_count)]
            group_found = set(sum(group_id, []))
            single_group = []
            for group in group_found:
                func_group = fparams.GetFuncGroup(group)
                single_group.append(func_group.GetProp('_Name'))
            return single_group

        records = []
        for item in self.io.read(FOpsinRecord):
            if item.get('group', None) is None and item.get('smiles', None) is not None:
                smiles = item['smiles']
                group = get_group(smiles)
                if group is not None:
                    item['group'] = group
            records.append(item)

        shutil.copy(FOpsinRecord, fcopy(FOpsinRecord))
        self.io.write(records, ftemp(FOpsinRecord))
        shutil.move(ftemp(FOpsinRecord), FOpsinRecord)

        logger.info("Successfully update the opsin_record.json")

    @staticmethod
    def get_failures():
        """
        chemicals in opsin site failure but success in ichem
        """

        result = JsonIO.read(FOpsinRecord)
        failure = [item['name'] for item in result if item['status'] == "FAILURE"]

        ichem = [key for key, value in JsonIO.read(FIChemRecord).items() if value.startswith("Ok.")]

        return [item for item in failure if item in ichem]

    @staticmethod
    def load_js(file):
        """
        Js func library

        Args:
            file: .js file

        Returns:
            content: .js file in string format

        """
        with open(file, 'r', encoding='UTF-8') as f:
            content = f.read()

        return content


if __name__ == '__main__':
    # get records
    # opsin = OPSINCrawler(FChemical)
    # opsin.get_records()

    # get failures
    # opsin = OPSINCrawler(FChemical)
    # chemicals = opsin.get_failures()

    # update records
    # opsin = OPSINCrawler()
    # opsin.update_records()

    # update formulas
    # opsin = OPSINCrawler()
    # opsin.update_formulas()

    # update groups
    # opsin = OPSINCrawler()
    # opsin.update_groups()

    records = JsonIO.read(FOpsinRecord)
    print()
