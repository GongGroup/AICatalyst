import json
import random
import shutil
import time
from pathlib import Path

import requests

from logger import logger

ChemDir = Path("./chemical")

IChemRoot = "http://www.ichemistry.cn/structure.asp"
IChemSearch = "http://www.ichemistry.cn/ketcher/Name2Structure/?action=mol&input="
headers = {
    "user-agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/104.0.5112.81 Safari/537.36 Edg/104.0.1293.47"}


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
            html = requests.get(IChemSearch + chemical, headers=headers)
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


if __name__ == '__main__':
    ichem = IChemCrawler(ChemDir / "chemical.json")
    ichem.get_htmls()
    print()
