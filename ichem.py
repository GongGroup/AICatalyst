import json
import random
import shutil
import time
from pathlib import Path
from requests_html import HTMLSession
import requests
from pyquery import PyQuery
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as ec

from driver import ChromeDriver
from logger import logger
from opsin import OPSINCrawler

ChemDir = Path("./chemical")

IChemRoot = "http://www.ichemistry.cn/structure.asp"
IChemStructure = "http://www.ichemistry.cn/ketcher/Name2Structure/?action=mol&input="
IChemName = "http://www.ichemistry.cn/ketcher/Name2Structure/?action=ename&input="
headers = {
    "user-agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) "
                  "Chrome/104.0.5112.81 Safari/537.36 Edg/104.0.1293.47"}


class IChemCrawler(object):
    def __init__(self, file):
        self.file = file
        self.chemicals = None

        self.load()

    def load(self):
        with open(self.file, "r") as f:
            self.chemicals = json.load(f)

    def get_htmls(self):
        with open(ChemDir / "record.json", "r", encoding="utf-8") as f:
            records = json.load(f)

        with open(ChemDir / "error", "r", encoding="utf-8") as f:
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
        try:
            with open(ChemDir / "_record.json", "w", encoding="utf-8") as f:
                json.dump(records, f)
        except Exception as error:
            raise error
        else:
            shutil.move(ChemDir / "_record.json", ChemDir / "record.json")

        # update error
        try:
            with open(ChemDir / "_error", "w", encoding="utf-8") as f:
                f.write('\n'.join(exclude))
        except Exception as error:
            raise error
        else:
            shutil.move(ChemDir / "_error", ChemDir / "error")

        logger.info("Update record.json successfully")

    @staticmethod
    def smile_name():

        with open("name_smile.json", "r", encoding='utf-8') as f:
            data = json.load(f)

        if not Path("name_name.json").exists():
            with open("name_name.json", "w", encoding="utf-8") as f:
                json.dump([], f, indent=2)

        with open("name_name.json", "r", encoding="utf-8") as f:
            records = json.load(f)

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
                with open("_name_name.json", "w", encoding="utf-8") as f:
                    json.dump(searched, f, indent=2)

        shutil.move("_name_name.json", "name_name.json")

    @staticmethod
    def name_smile(names):
        url = 'http://127.0.0.1:5500/structure.html?input='

        with open('name_smile.json', 'r', encoding='utf-8') as f:
            result = json.load(f)

        searched = [list(item.keys())[0] for item in result]
        session = HTMLSession()
        for name in names:
            if name in searched:
                continue

            response = session.get(url=url + name)
            response.html.render(wait=2, sleep=2)
            smiles = response.html.find('#iename')[0].text
            logger.info(f"{name} : {smiles}")
            result.append({name: smiles})

        try:
            with open('_smile.json', 'w', encoding='utf-8') as f:
                json.dump(result, f, indent=2)
        except:
            logger.error("Error!")
        else:
            shutil.move('_smile.json', 'name_smile.json')


if __name__ == '__main__':
    # ichem = IChemCrawler(ChemDir / "chemical.json")
    # ichem.get_htmls()
    # failures = OPSINCrawler.failure()
    # IChemCrawler.name_smile(failures)
    # IChemCrawler.smile_name()
    with open("chemical_new.json", "r", encoding='utf-8') as f:
        data = json.load(f)

    new_data = [item['new_name'] for item in data]

    with open("chemical_new.json", "w", encoding='utf-8') as f:
        json.dump(new_data, f)

    print()
