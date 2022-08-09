import json


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


class ReaxysRecords(object):
    def __init__(self, records_file: str):
        self.records_file = records_file
        self.records = None

        self._records = None

        self.flatten()
        self.register()

    def flatten(self):
        records = []
        with open(self.records_file, "r", encoding="utf-8") as f:
            original_records = json.load(f)

        for reaction in original_records:
            header = reaction[0]
            for record in reaction[1:]:
                record.update(header)
                records.append(record)
        self._records = records

    def register(self):
        self.records = [ReaxysRecord(item) for item in self._records]


if __name__ == '__main__':
    reaxys = ReaxysRecords("reaxys_json.json")
    print()
