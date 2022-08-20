from collections import Counter, defaultdict
from typing import Iterable

import numpy as np
from rdkit import Chem, DataStructs
from rdkit import RDLogger

from constant import FOpsinRecord, FChemical, FReaxys, FReactions
from fio import JsonIO, CatalystJsonIO
from species import MCatalyst

# Variable constant
ChemInfo = {item['name']: {key: value for key, value in item.items() if key != 'name'}
            for item in JsonIO.read(FOpsinRecord)}

# close the rdkit warning
RDLogger.DisableLog('rdApp.warning')


def to_int(values):
    return [int(value) for value in values]


def to_float(values):
    formatted_values = []
    for value in values:
        if ' - ' in value:
            formatted_values.append(list(map(float, [item for item in value.split(' - ') if len(item)])))
        else:
            formatted_values.append(float(value))
    return formatted_values


def to_fingerprint(values):
    fps = []
    for value in values:
        smiles = ChemInfo[value].get('smiles', None)
        if smiles is None:
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        try:
            fps.append(np.array(Chem.RDKFingerprint(mol).ToList()))
        except IndexError:
            fps.append(np.array(list(Chem.RDKFingerprint(mol).ToBitString())).astype(int))
    return fps


class ReaxysRecord(object):
    _total_attr = ["reactant", "product", "SKW", "reaction type", "title", "text", "PRO", "productivity", "reagent",
                   "solvent", "time", "temperature", "pressure", "doi", "journal", "vol", "year", "page", "ISSN",
                   "citation"]

    _trans_func = {
        "citation": to_int,
        "productivity": to_float,
        "time": to_float,
        "temperature": to_float,
        "pressure": to_float,
        "reactant": to_fingerprint,
        "product": to_fingerprint,
        "PRO": to_fingerprint,
        "reagent": to_fingerprint,
        "solvent": to_fingerprint,
    }

    def __init__(self, record: dict, transform=True):
        self._record = record

        self.initialize()
        if transform:
            self.transform()

    def initialize(self):
        for key, value in self._record.items():
            if key == "yield":
                key = "productivity"
            setattr(self, key, value)

    def transform(self):
        for key in self._record.keys():
            if key in self._trans_func.keys():
                setattr(self, key, self._trans_func[key](getattr(self, key)))

    def to_dict(self):
        return {attr: getattr(self, attr) for attr in ReaxysRecord._total_attr}


class ReaxysRecords(list):
    def __init__(self, seq: Iterable):
        super(ReaxysRecords, self).__init__(seq)
        self.io = JsonIO

    def output_chemicals(self, name=FChemical):
        chemicals = []
        for record in self:
            chemicals.append([record.PRO, record.product, record.reactant, record.reagent, record.solvent])
        total_chemicals = list(set(sum(sum(chemicals, []), [])))

        self.io.write(total_chemicals, name)


class ReaxysReaction(object):
    _total_attr = ["reactant", "reagent", "solvent", "product", "PRO", "productivity"]

    def __init__(self, record: dict):
        self._record = record
        self.reagent = []
        self.productivity = []

        self.reactant = []
        self.catalyst = []
        self.product = []

        self.initialize()

    def initialize(self):
        for key, values in self._record.items():
            for value in values:
                if key in ["reactant", "reagent", "solvent"]:
                    try:
                        catalyst = MCatalyst(value)
                    except ValueError:
                        self.reagent.append(value)
                    else:
                        self.catalyst.append(catalyst)
                elif key == "yield":
                    self.productivity.append(value)
                elif key in ["PRO", "product"]:
                    self.product.append(value)
                else:
                    continue

        self.product = list(set(self.product))
        self.reagent = list(set(self.reagent))
        self.productivity = list(set(self.productivity))
        # search reactant to product from fingerprint mapping
        for product in self.product:
            try:
                pro_smiles = ChemInfo[product]['smiles']
            except KeyError:
                continue
            pro_fp = Chem.RDKFingerprint(Chem.MolFromSmiles(pro_smiles))
            similarity = {}

            for reactant in self.reagent:
                try:
                    rea_smiles = ChemInfo[reactant]['smiles']
                except KeyError:
                    continue
                rea_fp = Chem.RDKFingerprint(Chem.MolFromSmiles(rea_smiles))
                similarity.update({reactant: DataStructs.FingerprintSimilarity(pro_fp, rea_fp)})

            self.reactant.append(similarity)

        self.reactant = [{key: value for key, value in item.items() if value == max(item.values())}
                         for item in self.reactant]

        new_reactant = []
        for item in self.reactant:
            reactant = list(item.keys())[0]
            group = ChemInfo[reactant]['group']
            new_reactant.append({reactant: {'group': group, 'similarity': item[reactant]}})
        self.reactant = new_reactant

        new_product = []
        for item in self.product:
            try:
                group = ChemInfo[item]['group']
            except KeyError:
                continue
            new_product.append({item: {'group': group}})

        self.product = new_product


class ReaxysReactions(list):
    def __init__(self, seq: Iterable):
        super(ReaxysReactions, self).__init__(seq)

    def to_json(self):
        reactions = []
        attrs = ['catalyst', 'product', 'productivity', 'reactant', 'reagent']
        for item in self:
            reactions.append({attr: getattr(item, attr) for attr in attrs})
        CatalystJsonIO.write(reactions, FReactions)


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
    def reactions(self):
        if self._records is None:
            self.flatten()
        seq = [ReaxysReaction(item) for item in self._records]

        return ReaxysReactions(seq)

    def records(self, transform=True):
        if self._records is None:
            self.flatten()
        seq = [ReaxysRecord(item, transform) for item in self._records]

        return ReaxysRecords(seq)

    def get_catalyst_reaction_mapping(self):
        reactions = self.reactions
        catalysts = Counter([item.name for reaction in reactions for item in reaction.catalyst])
        catalysts_sort = sorted(catalysts.items(), key=lambda x: x[1], reverse=True)
        catalyst_reaction_mapping = defaultdict(list)
        for key, value in catalysts_sort:
            for reaction in reactions:
                if key in [item.name for item in reaction.catalyst]:
                    for i, j in zip(reaction.reactant, reaction.product):
                        catalyst_reaction_mapping[key].append(
                            (list(i.values())[0]['group'], list(j.values())[0]['group']))
        return catalyst_reaction_mapping


if __name__ == '__main__':
    # get records
    # reaxys = Reaxys(FReaxys)
    # records = reaxys.records()

    # get reactions
    reaxys = Reaxys(FReaxys)
    reaxys.reactions.to_json()
    # catalyst_reaction_mapping = reaxys.get_catalyst_reaction_mapping()

    # for reaction in reactions:
    #     if 'palladium diacetate' in [item.name for item in reaction.catalyst]:
    #         test.append(reaction.product)
    # print(reactions)
    # shuffle_index = list(range(len(records)))
    # random.shuffle(shuffle_index)
    # length = len(records)
    # train_size = int(length * 0.8)
    # train_records = [records[i].to_dict() for i in shuffle_index[:train_size]]
    # test_records = [records[i].to_dict() for i in shuffle_index[train_size:]]
    # np.save("train_set.npy", train_records, allow_pickle=True)
    # np.save("test_set.npy", test_records, allow_pickle=True)
    # train_records = np.load("train_set_v1.npy", allow_pickle=True)
    # test_records = np.load("test_set_v1.npy", allow_pickle=True)

    print()
