from database.parse_reaxys import ParseReaxys
from database.extract_chemical import ExtractChemical
from database.opsin import OPSINCrawler
from database.chemical_search import SucReaxys

class ChemDatabase:
    def __init__(self):
        pass

    def sequential(self, *args):
        self._sequential = args

    def work(self):
        for method_name in self._sequential:
            tmp_name = method_name()
            tmp_name.forward()


if __name__ == '__main__':
    test = ChemDatabase()
    test.sequential(ParseReaxys, ExtractChemical, OPSINCrawler, SucReaxys)
    #test.sequential(OPSINCrawler, SucReaxys)
    test.work()