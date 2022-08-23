import random
import shutil
import time
from pathlib import Path

import requests
from pyquery import PyQuery

from common.fio import JsonIO, ftemp
from common.logger import logger

# Net const
IChemRoot = "http://www.ichemistry.cn/structure.asp"
IChemStructure = "http://www.ichemistry.cn/ketcher/Name2Structure/?action=mol&input="
IChemName = "http://www.ichemistry.cn/ketcher/Name2Structure/?action=ename&input="
IChemFormula = "http://search.ichemistry.cn/?keys="
headers = {
    "user-agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) "
                  "Chrome/104.0.5112.81 Safari/537.36 Edg/104.0.1293.47"}

# Directory const
ChemDir = Path("../chemical")

# file const
FChemical = ChemDir / "chemical.json"
FFormula = ChemDir / "formula.json"
FRecord = ChemDir / "ichem_record.json"  # success structure in sdf
FNameSmile = ChemDir / "name_smile.json"  # name~smile mapping
FNameName = ChemDir / "name_name.json"  # name~name mapping


class IChemCrawler(object):
    def __init__(self, file=None):
        self.io = JsonIO
        self.file = file
        self.chemicals = self.io.read(self.file) if self.file is not None else None

    def get_structure(self):
        """
        According to the chemical.json [list] crawler the sdf structure
        """
        records = self.io.read(FRecord)

        for chemical in self.chemicals:
            if chemical in records.keys():
                logger.debug(f"`{chemical}` has been searched, continue")
                continue
            time.sleep(random.random() + 0.5)
            html = requests.get(IChemStructure + chemical, headers=headers)
            if html.status_code == 200:
                content = str(html.content, encoding='raw_unicode_escape')
                if not content.startswith("Error"):
                    records.update({chemical: content})
                    logger.info(f"`{chemical}` store in database")
                else:
                    records.update({chemical: "Error"})
                    logger.warning(f"`{chemical}` not in this site")
            else:
                records.update({chemical: "Error"})
                logger.warning(f"`{chemical}` search error")

        # update record.json
        self.io.write(records, ftemp(FRecord))
        shutil.move(ftemp(FRecord), FRecord)

        logger.info("Update record.json successfully")

    @staticmethod
    def get_name():
        """
        According to the smile obtain the name
        """
        data = JsonIO.read(FNameSmile)

        if not Path(FNameName).exists():
            JsonIO.write([], FNameName)

        records = JsonIO.read(FNameName)
        searched = [item['old_name'] for item in records]
        get_kv = lambda x: list(x.items())[0]
        for index, item in enumerate(data):
            name, smiles = get_kv(item)
            if name in searched:
                continue
            time.sleep(random.random() + 0.5)
            response = requests.get(IChemName + smiles, headers=headers)
            text = response.text
            if text.startswith("Ok."):
                new_name = text.split("Ok.")[1]
                logger.info(f"{name} => {new_name}")
            else:
                new_name = 'FAILURE'
                logger.warning(f"transform error")
            searched.append({"old_name": name, "smiles": smiles, "new_name": new_name})

            if index % 20 == 0 or index == len(data) - 1:
                JsonIO.write(searched, ftemp(FNameName))

        shutil.move(ftemp(FNameName), FNameName)

    def get_formula(self):

        formulas = self.io.read(FFormula)
        formulas_name = list(formulas.keys())

        for index, chemical in enumerate(self.chemicals):
            if chemical in formulas_name:
                continue
            response = requests.get(IChemFormula + chemical, headers=headers)
            response.encoding = 'gbk'
            doc = PyQuery(response.text)
            records = doc.find('#container-right tr:nth-child(n+2) td:nth-child(5)')
            if len(records):
                formula = [record.text() for record in records.items()]
                formulas.update({chemical: formula[0]})
                logger.info(f"Transform [{chemical}] => [{formula[0]}] successfully")
            else:
                formulas.update({chemical: "FAILURE"})
                logger.info(f"Transform [{chemical}] error")

            if index % 20 == 0 or index == len(self.chemicals) - 1:
                JsonIO.write(formulas, ftemp(FFormula))
                shutil.move(ftemp(FFormula), FFormula)


if __name__ == '__main__':
    # get_structure
    # ichem = IChemCrawler(file=FChemical)
    # ichem.get_structure()

    # get new_name
    # IChemCrawler.get_name()

    # get formula
    ichem = IChemCrawler(file=FChemical)
    ichem.get_formula()

    print()
