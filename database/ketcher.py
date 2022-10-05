import shutil

from requests_html import HTMLSession

from common.constant import ChemDir
from common.file import JsonIO, ftemp
from common.logger import logger
from database.opsin import OPSINCrawler

# Net constant
URL = 'http://127.0.0.1:5500/structure.html?input='

# file const
FNameSmile = ChemDir / "name_smile.json"  # name~smile mapping


class Ketcher(object):
    def __init__(self):
        self.io = JsonIO

    def get_smile(self, names):
        """
        According to the name obtain the smile

        Args:
            names: chemicals' name, list

        """
        # local net, need structure.html && ketcher && ketcher.js
        session = HTMLSession()

        results = self.io.read(FNameSmile)
        searched = [list(item.keys())[0] for item in results]
        for name in names:
            if name in searched:
                continue

            response = session.get(url=URL + name)
            response.html.render(wait=2, sleep=2)
            smiles = response.html.find('#iename')[0].text
            logger.info(f"{name} : {smiles}")
            results.append({name: smiles})

        try:
            self.io.write(results, ftemp(FNameSmile))
        except:
            logger.error("Error!")
        else:
            shutil.move(ftemp(FNameSmile), FNameSmile)


if __name__ == '__main__':
    chemicals = OPSINCrawler.get_failures()
    ketcher = Ketcher()
    ketcher.get_smile(chemicals)
