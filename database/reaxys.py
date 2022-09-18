from collections import Counter, defaultdict, namedtuple
from itertools import combinations
from typing import Iterable

import numpy as np
from rdkit import Chem, DataStructs
from rdkit import RDLogger

from common.constant import FOpsinRecord, FChemical, FReactions, FReaxysYield, DrawDir, FReaxys
from common.draw import draw_svg
from common.fio import JsonIO, CatalystJsonIO, md5
from common.species import MCatalyst

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

            effective_reagent = {}
            for reactant in self.reagent:
                try:
                    rea_smiles = ChemInfo[reactant]['smiles']
                except KeyError:
                    continue
                else:
                    rea_fp = Chem.RDKFingerprint(Chem.MolFromSmiles(rea_smiles))
                    effective_reagent.update({reactant: rea_fp})

            candidate_reagent = [{i: effective_reagent[i] for i in item} for num in range(len(effective_reagent))
                                 for item in combinations(effective_reagent, num + 1)]

            reactant_p = []
            for item in candidate_reagent:
                sum_fp = Chem.RDKFingerprint(Chem.MolFromSmiles('[H]'))
                for value in item.values():
                    sum_fp |= value

                similarity = {"reactant": list(item.keys()),
                              "similarity": DataStructs.FingerprintSimilarity(pro_fp, sum_fp)}

                reactant_p.append(similarity)

            max_similarity = 0.
            max_reactant = None
            for item in reactant_p:
                if item['similarity'] >= max_similarity:
                    max_similarity = item['similarity']
                    max_reactant = item['reactant']
            self.reactant.append({"reactant": max_reactant, "similarity": max_similarity})
        # self.reactant = [{key: value for key, value in item.items() if value == max(item.values())}
        #                  for item in self.reactant]

        new_reactant = []
        for item in self.reactant:
            reactants = item['reactant']
            group = sum([ChemInfo[reactant]['group'] for reactant in reactants], [])
            new_reactant.append({"reactant": item['reactant'], 'group': group, 'similarity': item['similarity']})
        self.reactant = new_reactant

        new_product = []
        for item in self.product:
            try:
                group = ChemInfo[item]['group']
            except KeyError:
                continue
            new_product.append({"product": item, 'group': group})
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
        self.file = file
        self._records = None

    def flatten(self):
        var_records = []

        original_records = JsonIO.read(self.file)

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

    @staticmethod
    def get_reaction_mapping():
        reactions = JsonIO.read(FReactions)
        catalysts = Counter([item for reaction in reactions for item in reaction['catalyst']])
        catalysts_sort = sorted(catalysts.items(), key=lambda x: x[1], reverse=True)
        mapping_count = {}
        mapping_type = defaultdict(list)
        Mapping = namedtuple("Mapping", ("reactant", "product"))
        for key, value in catalysts_sort:
            for reaction in reactions:
                if key in [item for item in reaction['catalyst']]:
                    for i, j in zip(reaction['reactant'], reaction['product']):
                        i_group = i['group']
                        j_group = j['group']

                        for item in i_group:
                            if item in j_group:
                                i_group.remove(item)
                                j_group.remove(item)

                        mapping = Mapping(tuple(i_group), tuple(j_group))
                        mapping_type[key].append(mapping)
            mapping_count.update({key: Counter(mapping_type[key])})
        return mapping_type, mapping_count


if __name__ == '__main__':
    # get records
    reaxys = Reaxys(FReaxysYield)
    records = reaxys.records(transform=False)

    # get reactions
    reactions = reaxys.reactions

    # reaxys.reactions.to_json()
    product_type = defaultdict(list)
    for reaction in reactions:
        for p in reaction.product:
            product_type[p['product']].append(reaction)
    product_type_sorted = {key: value
                           for key, value in sorted(product_type.items(), key=lambda x: len(x[1]), reverse=True)}

    new_product_type = {key: value for key, value in product_type_sorted.items() if len(value) > 3}

    # name = 'toluene'
    # special_catalyst = [reaction.catalyst for reaction in new_product_type[name]]
    # special_records = [record for record in records if (name) in record.product and len(record.doi)]
    # special_product = new_product_type[name]
    # for key, value in new_product_type.items():
    #     mol = Chem.MolFromSmiles(ChemInfo[key]['smiles'])
    #     draw_svg(mol, file=DrawDir / f"{md5(key)}.svg")

    reaction_mapping = Reaxys.get_reaction_mapping()
    # for key, value in reaction_mapping[1].items():
    #     if len(value):
    #         for item in value.most_common(1):
    #             if item[1] >= 10:
    #                 print(key, value.most_common(1))

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
    print()
