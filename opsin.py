import copy
import json
import shutil
from functools import partial
from pathlib import Path

import execjs
import requests

from fio import JsonIO, temp
from logger import logger

# Net const
OPSINRoot = 'https://opsin.ch.cam.ac.uk/'
headers = {
    "user-agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) "
                  "Chrome/104.0.5112.81 Safari/537.36 Edg/104.0.1293.47"}

# Directory const
ChemDir = Path("./chemical")

# File const
FJSFunc = "func.js"
FRecord = ChemDir / "opsin_record.json"
FIChemRecord = ChemDir / "record.json"


class OPSINCrawler(object):
    def __init__(self, file):
        self.io = JsonIO
        self.file = file
        self.chemicals = self.io.read(self.file)
        self.encode = partial(execjs.compile(OPSINCrawler.load_js(FJSFunc)).call, "encode")

    def get_htmls(self):
        """
        According to the chemical.json [list] obtain the smiles, inchi, inchi-key and so on
        """

        if not Path(FRecord).exists():
            self.io.write([], FRecord)

        records = self.io.read(FRecord)
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

            if index % 20 == 0 or index == len(self.chemicals) - 1:
                self.io.write(records, temp(FRecord))
                shutil.move(temp(FRecord), FRecord)

        logger.info("Successfully update the opsin_record.json")

    def failure(self):
        """
        chemicals in opsin site failure but success in ichem
        """

        result = self.io.read(FRecord)
        failure = [item['name'] for item in result if len(item) == 2]

        ichem = list(self.io.read(FIChemRecord).keys())

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
    # opsin = OPSINCrawler("chemical_new.json")
    # opsin.get_htmls()

    with open("chemical/name_name.json", "r", encoding='utf-8') as f:  # old_name ~ new_name mapping
        name2name = json.load(f)

    with open("chemical/opsin_record_old.json", "r", encoding='utf-8') as f:  # last opsin_record
        opsin_record = json.load(f)

    with open("chemical/opsin_record.json", "r", encoding='utf-8') as f:  # newly succeed opsin_record
        opsin_record_new = json.load(f)

    new_records = []

    # transform the list to dict for the subsequent quickly mapping
    name2name_dict = {item['old_name']: item['new_name'] for item in name2name}
    opsin_record_new_dict = {item['name']: {key: item[key] for key in item.keys() if key != "name"} for item in
                             opsin_record_new}
    for item in copy.deepcopy(opsin_record):
        if item['name'] in list(name2name_dict.keys()):
            # obtain new_name from old_name
            new_name = name2name_dict[item['name']].replace('<i>', '').replace('</i>', '')
            if new_name != "FAILURE":
                for key in opsin_record_new_dict[new_name].keys():
                    item[key] = opsin_record_new_dict[new_name][key]  # substitute old_info with new_info
                    item['new_name'] = new_name  # add `new_name` key
            else:
                item['new_name'] = "FAILURE"
        else:
            item['new_name'] = item['name']
        new_records.append(item)

    new_records_sort = [{item2[0]: item2[1] for item2 in sorted(item.items(), key=lambda x: x[0])} for item in
                        new_records]

    with open("chemical/opsin_record_new.json", "w", encoding='utf-8') as f:
        json.dump(new_records_sort, f, indent=2)
