import json
import shutil
from functools import partial
from pathlib import Path

import execjs
import requests

from logger import logger

ChemDir = Path("./chemical")
OPSINRoot = 'https://opsin.ch.cam.ac.uk/'

headers = {
    "user-agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) "
                  "Chrome/104.0.5112.81 Safari/537.36 Edg/104.0.1293.47"}


class OPSINCrawler(object):
    def __init__(self, file):
        self.file = file
        self.chemicals = None
        self.encode = partial(execjs.compile(OPSINCrawler.js_from_file("func.js")).call, "encode")
        self.load()

    def load(self):
        with open(self.file, "r") as f:
            self.chemicals = json.load(f)

    def get_htmls(self):
        with open(ChemDir / "opsin_record.json", "r", encoding="utf-8") as f:
            records = json.load(f)

        records_name = [item['name'] for item in records]

        for chemical in self.chemicals:
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

        # update opsin_record.json
        try:
            with open(ChemDir / "_record.json", "w", encoding="utf-8") as f:
                json.dump(records, f, indent=2)
        except Exception as error:
            raise error
        else:
            shutil.move(ChemDir / "_record.json", ChemDir / "opsin_record.json")
        logger.info("Successfully update the opsin_record.json")

    @staticmethod
    def js_from_file(file):
        with open(file, 'r', encoding='UTF-8') as f:
            content = f.read()

        return content


if __name__ == '__main__':
    opsin = OPSINCrawler("chemical/chemical.json")
    opsin.get_htmls()
    # result = opsin.encode.call("encode", "2,4,6-trinitrotoluene")
    # execjs.eval("encodeURIComponent", "2,4,6-trinitrotoluene")
    print()
