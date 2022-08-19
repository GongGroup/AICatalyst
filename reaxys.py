import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np
from rdkit import Chem, DataStructs, RDConfig
from rdkit.Chem import FragmentCatalog

from fio import JsonIO

# Directory constant
from species import MCatalyst

ChemDir = Path("./chemical")

# file constant
FReaxys = ChemDir / "opsin_reaxys.json"
FRecord = ChemDir / "opsin_record_new.json"

# Variable constant
ChemInfo = {item['name']: {key: value for key, value in item.items() if key != 'name'} for item in JsonIO.read(FRecord)}


def list_int(values):
    return [int(value) for value in values]


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
                fps.append(np.array(list(Chem.RDKFingerprint(Chem.MolFromSmiles(smiles)).ToBitString())).astype(int))
        return fps
    else:
        return [ChemInfo[value].get('smiles', None) for value in values]


class ReaxysRecord(object):
    _total_attr = ["reactant", "product", "SKW", "reaction type", "title", "text", "PRO", "productivity", "reagent",
                   "solvent", "time", "temperature", "pressure", "doi", "journal", "vol", "year", "page", "ISSN",
                   "citation"]

    _trans_func = {
        "citation": list_int,
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

    def to_dict(self):
        return {attr: getattr(self, attr) for attr in ReaxysRecord._total_attr}


class ReaxysRecords(list):
    def __init__(self, seq: Iterable):
        super(ReaxysRecords, self).__init__(seq)


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
            pro_smiles = ChemInfo[product]['smiles']
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

        # add group info
        fName = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
        fparams = FragmentCatalog.FragCatParams(1, 6, fName)

        def group_info(name):
            smiles = ChemInfo[name]['smiles']
            m = Chem.MolFromSmiles(smiles)
            fcat = FragmentCatalog.FragCatalog(fparams)
            fcgen = FragmentCatalog.FragCatGenerator()
            fcgen.AddFragsFromMol(m, fcat)

            group_count = fcat.GetNumEntries()
            group_id = [list(fcat.GetEntryFuncGroupIds(i)) for i in range(group_count)]
            group_found = set(sum(group_id, []))
            single_group = []
            for group in group_found:
                func_group = fparams.GetFuncGroup(group)
                single_group.append(func_group.GetProp('_Name'))
            return single_group

        new_reactant = []
        for item in self.reactant:
            reactant = list(item.keys())[0]
            single_group = group_info(reactant)
            new_reactant.append({reactant: {'group': single_group, 'similarity': item[reactant]}})
        self.reactant = new_reactant

        new_product = []
        for item in self.product:
            single_group = group_info(item)
            new_product.append({item: {'group': single_group}})

        self.product = new_product


class ReaxysReactions(list):
    def __init__(self, seq: Iterable):
        super(ReaxysReactions, self).__init__(seq)


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
    reactions = reaxys.reactions
    # for reaction in reactions:
    #     print(reaction.reactant, [item.name for item in reaction.catalyst], reaction.product)
    catalysts = Counter([item.name for reaction in reactions for item in reaction.catalyst])
    catalysts_sort = sorted(catalysts.items(), key=lambda x: x[1], reverse=True)
    catalyst_reaction_mapping = defaultdict(list)
    for key, value in catalysts_sort:
        for reaction in reactions:
            if key in [item.name for item in reaction.catalyst]:
                for i, j in zip(reaction.reactant, reaction.product):
                    # if (list(i.keys())[0], j) not in catalyst_reaction_mapping[key]:
                    catalyst_reaction_mapping[key].append((list(i.values())[0]['group'], list(j.values())[0]['group']))

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
    # smi = 'C1(=CC=C(C=C1)I)C'
    # smi = 'CC1=CC=C(C=O)C=C1'
    # mol = Chem.MolFromSmiles(smi)
    #
    # for _ in range(10):
    #     smi = Chem.MolToSmiles(mol, doRandom=True)
    #     print(smi)

    print()
