from pathlib import Path
from typing import Iterable

import numpy as np
from rdkit import Chem

from fio import JsonIO

# Directory constant
ChemDir = Path("./chemical")

# file constant
FReaxys = ChemDir / "opsin_reaxys.json"
FRecord = ChemDir / "opsin_record_new.json"

# Variable constant
ChemInfo = {item['name']: {key: value for key, value in item.items() if key != 'name'} for item in JsonIO.read(FRecord)}


def list_float(values):
    formatted_values = []
    for value in values:
        if ' - ' in value:
            formatted_values.append(list(map(float, [item for item in value.split(' - ') if len(item)])))
        else:
            formatted_values.append(float(value))
    return formatted_values


def list_smiles(values, fingerprint=True):
    if fingerprint:
        fps = []
        for value in values:
            smiles = ChemInfo[value].get('smiles', None)
            try:
                fps.append(np.array(Chem.RDKFingerprint(Chem.MolFromSmiles(smiles)).ToList()))
            except IndexError:
                fps.append(smiles)
        return fps
    else:
        return [ChemInfo[value].get('smiles', None) for value in values]


class ReaxysRecord(object):
    _total_attr = ["reactant", "product", "SKW", "reaction type", "title", "text", "PRO", "productivity", "reagent",
                   "solvent", "time", "temperature", "pressure", "doi", "journal", "vol", "year", "page", "ISSN",
                   "citation"]

    _trans_func = {
        "productivity": list_float,
        "time": list_float,
        "temperature": list_float,
        "pressure": list_float,
        "reactant": list_smiles,
        "product": list_smiles,
        "PRO": list_smiles,
        "reagent": list_smiles,
    }

    def __init__(self, record: dict):
        self._record = record

        self.initialize()

    def initialize(self):
        for key, value in self._record.items():
            if key == "yield":
                key = "productivity"
            if key in self._trans_func.keys():
                setattr(self, key, self._trans_func[key](value))
            else:
                setattr(self, key, value)


class ReaxysRecords(list):
    def __init__(self, seq: Iterable):
        super(ReaxysRecords, self).__init__(seq)


class Reaxys(object):
    def __init__(self, file):
        self.io = JsonIO
        self.file = file
        self._records = None

    def flatten(self):
        var_records = []

        original_records = self.io.read(self.file)

        for reaction in original_records:
            header = reaction[0]
            for record in reaction[1:]:
                record.update(header)
                var_records.append(record)
        self._records = var_records

    @property
    def records(self):
        if self._records is None:
            self.flatten()
        seq = [ReaxysRecord(item) for item in self._records]

        return ReaxysRecords(seq)

    def output_chemicals(self, name=ChemDir / "chemical.json"):
        chemicals = []
        for record in self.records:
            chemicals.append([record.PRO, record.product, record.reactant, record.reagent])
        total_chemicals = list(set(sum(sum(chemicals, []), [])))

        self.io.write(total_chemicals, name)


if __name__ == '__main__':
    reaxys = Reaxys(FReaxys)
    records = reaxys.records
    for record in records:
        for attr in ["reactant", "product", "PRO", "reagent"]:
            for item in getattr(record, attr):
                if not isinstance(item, np.ndarray):
                    print(item)
    # ms = [Chem.MolFromSmiles(record.PRO[0]) for record in records]
    # fps = np.array([Chem.RDKFingerprint(x).ToList() for x in ms])
    # sum_fps = np.sum(fps, axis=1)
    print()
