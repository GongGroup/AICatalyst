# input file
from common.constant import FReaxys
# output file
from common.constant import FChemical

from common.fio import JsonIO
from common.logger import logger
import pandas as pd
from functools import reduce


class ExtractChemical:
    def __init__(self, features=None):
        if features:
            self.features = features
        else:
            self.features = ['reagent', 'solvent', 'reactant', 'product']
        self.reaxys = pd.DataFrame.from_records(JsonIO.read(FReaxys))
        self.chemicals = []

    def forward(self):
        for feature in self.features:
            self.chemicals.extend(self._set_chemical(feature))
        JsonIO.write(self.chemicals, FChemical)
        logger.info("Successfully extract chemistry")

    def _set_chemical(self, feature):
        return list(set(reduce(lambda a, b: a + b, self.reaxys[feature].to_list())))


if __name__ == '__main__':
    tmp = ExtractChemical()
    tmp.forward()
