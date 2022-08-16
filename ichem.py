import json
import random
import shutil
import time
from pathlib import Path

import requests
from requests_html import HTMLSession

from fio import JsonIO, temp
from logger import logger

# Net const
IChemRoot = "http://www.ichemistry.cn/structure.asp"
IChemStructure = "http://www.ichemistry.cn/ketcher/Name2Structure/?action=mol&input="
IChemName = "http://www.ichemistry.cn/ketcher/Name2Structure/?action=ename&input="
headers = {
    "user-agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) "
                  "Chrome/104.0.5112.81 Safari/537.36 Edg/104.0.1293.47"}

# Directory const
ChemDir = Path("./chemical")

# file const
FRecord = ChemDir / "record.json"  # success structure in sdf
FError = ChemDir / "error"  # error chemical
FNameSmile = ChemDir / "name_smile.json"  # name~smile mapping
FNameName = ChemDir / "name_name.json"  # name~name mapping


class IChemCrawler(object):
    def __init__(self, file):
        self.io = JsonIO
        self.file = file
        self.chemicals = self.io.read(self.file)

    def get_htmls(self):
        """
        According to the chemical.json [list] crawler the sdf structure
        """
        records = self.io.read(FRecord)

        with open(FError, "r", encoding="utf-8") as f:
            exclude = [line.rstrip() for line in f.readlines()]

        for chemical in self.chemicals:
            if chemical in records.keys() or chemical in exclude:
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
                    exclude.append(chemical)
                    logger.warning(f"`{chemical}` not in this site")
            else:
                exclude.append(chemical)
                logger.warning(f"`{chemical}` search error")

        # update record.json
        self.io.write(records, temp(FRecord))
        shutil.move(temp(FRecord), FRecord)

        # update error
        with open(FError, "w", encoding="utf-8") as f:
            f.write('\n'.join(exclude))

        logger.info("Update record.json successfully")

    def smile_name(self):
        """
        According to the smile obtain the name
        """
        data = self.io.read(FNameSmile)

        if not Path(FNameName).exists():
            self.io.write([], FNameName)

        records = self.io.read(FNameName)
        searched = [item['old_name'] for item in records]
        for index, item in enumerate(data):
            if item['old_name'] in searched:
                continue
            time.sleep(random.random() + 0.5)
            response = requests.get(IChemName + item['smiles'], headers=headers)
            text = response.text
            if text.startswith("Ok."):
                new_name = text.split("Ok.")[1]
                item['new_name'] = new_name
                logger.info(f"{item['old_name']} => {item['new_name']}")
            else:
                item['new_name'] = 'FAILURE'
                logger.warning(f"transform error")
            searched.append(item)

            if index % 20 == 0 or index == len(data) - 1:
                self.io.write(searched, temp(FNameName))

        shutil.move(temp(FNameName), FNameName)

    def name_smile(self, names):
        """
        According to the name obtain the smile

        Args:
            names: chemicals' name, list

        """
        url = 'http://127.0.0.1:5500/structure.html?input='  # local net, need structure.html && ketcher && ketcher.js
        session = HTMLSession()

        results = self.io.read(FNameSmile)
        searched = [list(item.keys())[0] for item in results]
        for name in names:
            if name in searched:
                continue

            response = session.get(url=url + name)
            response.html.render(wait=2, sleep=2)
            smiles = response.html.find('#iename')[0].text
            logger.info(f"{name} : {smiles}")
            results.append({name: smiles})

        try:
            self.io.write(results, temp(FNameSmile))
        except:
            logger.error("Error!")
        else:
            shutil.move(temp(FNameSmile), FNameSmile)


if __name__ == '__main__':
    ichem = IChemCrawler(ChemDir / "chemical.json")
    # ichem.get_htmls()
    # failures = OPSINCrawler.failure()
    # IChemCrawler.name_smile(failures)
    # IChemCrawler.smile_name()
    # with open("chemical/chemical_new.json", "r", encoding='utf-8') as f:
    #     data = json.load(f)
    #
    # new_data = [item['new_name'] for item in data]
    #
    # with open("chemical/chemical_new.json", "w", encoding='utf-8') as f:
    #     json.dump(new_data, f)
    #
    print()
