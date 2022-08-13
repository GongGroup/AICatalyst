import json
from pathlib import Path
from typing import Iterable

ChemDir = Path("./chemical")


def list_float(values):
    if isinstance(values, list):
        if len(values):
            try:
                return [float(value) for value in values]
            except ValueError:
                return values
        else:
            return values
    elif isinstance(values, str):
        return float(values)
    else:
        raise TypeError(f"{type(values)} can't transform to float")


class ReaxysRecord(object):
    _total_attr = ["reactant", "product", "SKW", "reaction type", "title", "text", "PRO", "productivity", "reagent",
                   "solvent", "time", "temperature", "pressure", "doi", "journal", "vol", "year", "page", "ISSN",
                   "citation"]

    _trans_func = {
        "productivity": list_float,
        "time": list_float,
        "temperature": list_float,
        "pressure": list_float,
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
        self.file = file
        self._records = None

    def flatten(self):
        records = []
        with open(self.file, "r", encoding="utf-8") as f:
            original_records = json.load(f)

        for reaction in original_records:
            header = reaction[0]
            for record in reaction[1:]:
                record.update(header)
                records.append(record)
        self._records = records

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
        with open(name, "w", encoding='utf-8') as f:
            json.dump(total_chemicals, f)


if __name__ == '__main__':
    reaxys = Reaxys("reaxys_json.json")
    reaxys.output_chemicals()
    print()
