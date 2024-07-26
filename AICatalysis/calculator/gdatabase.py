import os
import re
import json
import pandas as pd
import numpy as np
from pathlib import Path
import shutil

from AICatalysis.calculator import descriptor
from AICatalysis.common.constant import *
from AICatalysis.common.file import CSVIO

from AICatalysis.calculator.gaussian import GJFFile
from AICatalysis.calculator.smiles import SmilesFile


class GaussianInDir:
    def __init__(self, dir_name, parent_path=GaussianInDataDir):
        self.parent_path = parent_path
        self.dir_name = dir_name
        self.dir_path = os.path.join(self.parent_path, self.dir_name)
        self.des_path = os.path.join(DescriptorDataDir, dir_name)

    def updatedir(self, smiles_file=None):
        if smiles_file is None:
            smiles_file = DefaultUpdateSmilesFile
        else:
            smiles_file = SmilesDir / smiles_file
        old_smiles_file = SmilesFile(ReactSmilesFile)
        update_smiles_file = SmilesFile(smiles_file)

        repeat_index = [index for index, smiles in enumerate(update_smiles_file.smiles.values)
                        if smiles in old_smiles_file.smiles.values]

        copy_name = []
        for index in repeat_index:
            tmp_smile = old_smiles_file[old_smiles_file.smiles == update_smiles_file.smiles[index]][0]
            copy_name.append(str(tmp_smile[0]) + '-' + tmp_smile[2])

        repeat_name = [str(id) + '-' + str(name) for id, name in zip(update_smiles_file.type_id[repeat_index],
                                                                update_smiles_file.name[repeat_index])]
        gjf_index = [index for index in range(len(update_smiles_file)) if index not in repeat_index]

        if self.dir_name == 'reaction_compound':
            old_smiles_file.concat(update_smiles_file)
            old_smiles_file.save()

        self._create_gjfs(update_smiles_file, gjf_index)
        self.update_descriptor_dir(copy_name, repeat_name)

    def _create_gjfs(self, smiles_file, index):
        if self.dir_name == 'reaction_compound':
            smiles_name = smiles_file.file_name[index]
        else:
            smiles_name = smiles_file.name[index]
        for smiles, name in zip(smiles_file.smiles[index], smiles_name):
            self.create_gjf(smiles, name)

    def create_gjf(self, smiles, name):
        gjf_file = GJFFile()
        gjf_file.write(smiles, self.dir_path, name)

    def update_descriptor_dir(self, copy_name, repeat_name):
        for copy_des, repeat_des in zip(copy_name, repeat_name):
            if os.path.exists(os.path.join(self.des_path, copy_des + '.json')):
                shutil.copyfile(os.path.join(self.des_path, copy_des + '.json'),
                            os.path.join(self.des_path, repeat_des + '.json'))
            else:
                CSVIO.write([copy_des, repeat_des], DescriptorDir_ / "hold_list", delimiter=' ')


class DescriptorDir:
    def __init__(self, root_path=DescriptorDataDir):
        self.root_path = root_path

    def makedir(self, dir_name):
        dir_path = os.path.join(self.root_path, dir_name)
        os.makedirs(dir_path)

    def makedirs(self, name_list):
        if len(name_list) < 1:
            exit('No dir')
        elif len(name_list) == 1:
            self.makedir(name_list)
        else:
            for dir_name in name_list:
                self.makedir(str(dir_name))

    def extractor(self, dir_name):
        pass


class GaussianOutDir:
    root_path = GaussianOutDataDir
    def __init__(self, chemical: str):
        self.out_path = self.root_path / chemical / 'out'
        self.fchk_path = self.root_path / chemical / 'fchk'
        #self.del_prefix()

    @staticmethod
    def out_file(chemical):
        out_files = [os.path.join(f"../../database/struct_data/gaussian_data/output/out/{chemical}/", item)
                     for item in os.listdir(f"../../database/struct_data/gaussian_data/output/out/{chemical}/")]
        return out_files

    @staticmethod
    def fchk_file(chemical):
        fchk_files = [os.path.join(f"../../database/struct_data/gaussian_data/output/fchk/{chemical}/", item)
                     for item in os.listdir(f"../../database/struct_data/gaussian_data/output/fchk/{chemical}/")]
        return fchk_files


    def extractor_liter(self, chemical, test=False):
        out_file_path = os.path.join(self.root_path, 'out', chemical)
        out_file_list = [os.path.join(out_file_path, out_file) for out_file in os.listdir(out_file_path)
                         if out_file.endswith('.out')]
        json_root_path = os.path.join(DescriptorDataDir, chemical)

        # smiles_table = pd.read_csv(ReactSmilesFile, delimiter=' ')
        # smiles_table.set_index(['index'], inplace=True)

        smiles_table = SmilesFile.read(ReactSmilesFile)

        # compound_descriptor_list = [descriptor.Descriptor(compound).print() for compound in out_file_list]

        for compound in out_file_list:
            json_path = Path(json_root_path) / (Path(compound).stem + '.json')
            if os.path.exists(json_path):
                continue

            try:
                smiles_type_id = int(Path(compound).stem.split('-')[0])
                smiles_name = Path(compound).stem.split('-')[1]
                smiles = smiles_table.loc[(smiles_table['type_id'] == smiles_type_id) &
                                          (smiles_table['out_name'] == smiles_name), 'smiles'].values[0]
                compound_descriptor = descriptor.Descriptor(compound, smiles=smiles).print(test)
            except:
                compound_descriptor = descriptor.Descriptor(compound).print(test)
            dict_json = json.dumps(compound_descriptor)
            with open(json_path, 'w+') as f:
                f.write(dict_json)


    def del_prefix(self):
        for liter_file in self.liter_dir:
            out_file_path = os.path.join(self.root_path, liter_file, 'out')
            fchk_file_path = os.path.join(self.root_path, liter_file, 'fchk')
            out_file_list = [os.path.join(out_file_path, out_file) for out_file in os.listdir(out_file_path)]
            fchk_file_list = [os.path.join(fchk_file_path, fchk_file) for fchk_file in os.listdir(fchk_file_path)]
            prefix = liter_file + '-'
            for out_file, fchk_file in zip(out_file_list, fchk_file_list):
                if os.path.basename(out_file).startswith(prefix):
                    os.rename(out_file, re.sub(prefix, '', out_file))
                if os.path.basename(fchk_file).startswith(prefix):
                    os.rename(fchk_file, re.sub(prefix, '', fchk_file))

    @property
    def liter_dir(self):
        liter_dir = os.listdir(self.root_path)
        return liter_dir


class ReactionDatabase():
    def __init__(self, root_path=ReactDataDir, react_data=ReactionFile,
                 des_data=TotalDesDir):
        self.root_path = root_path
        self.des_data = des_data
        self.database, self.not_find = self._concat(react_data)
        self.output()

    def _concat(self, react_data):
        react_data = self._pre_transform(react_data)

        concat_data = pd.DataFrame()
        not_find_file = set()
        for _, sample in react_data.iterrows():
            concat_sample, not_find_file = self._concat_chemical(sample, not_find_file)
            concat_data = pd.concat([concat_data, concat_sample], axis=1)
        concat_data = concat_data.T
        concat_data.index = range(len(concat_data))

        return concat_data, not_find_file

    def _pre_transform(self, react_data):
        # fill paper id
        react_data = pd.read_csv(react_data)
        react_data.iloc[:, 0].fillna(method='ffill', inplace=True)
        react_data['paper_id'] = react_data['paper_id'].astype('int')

        return react_data

    def _concat_chemical(self, sample, not_find_file, pca_data=False):
        concat_sample = sample.copy()
        for chemical in ChemicalList:
            concat_sample.drop(chemical, inplace=True)

        reference_index = str(sample[0])

        for chemical in ChemicalList[:3]:
            if sample[chemical] is np.nan:
                continue
            sub_pro = reference_index + '-' + sample[chemical]
            print(sub_pro)
            sub_pro_json = self.des_data / 'reaction_compound' / (sub_pro + '.json')
            if not os.path.exists(sub_pro_json):
                import shutil
                smiles_file = SmilesFile(ReactSmilesFile)
                sub_pro_smiles = smiles_file[smiles_file.file_name == sub_pro, 'smiles'][0]
                candidate_json = [self.des_data / 'reaction_compound' / (smiles_list[3] + '.json')
                                  for smiles_list in smiles_file[smiles_file.smiles == sub_pro_smiles]
                                  if smiles_list[3] != sub_pro]
                for cur_json in candidate_json:
                    if os.path.exists(cur_json):
                        shutil.copyfile(cur_json, sub_pro_json)
                        break
                else:
                    not_find_file.add(sub_pro_json.parent.name + '/' + sub_pro_json.name)
                    continue
            reagent = pd.read_json(sub_pro_json, orient='index')
            index = [chemical + '_' + index for index in reagent.index]
            reagent.index = index
            concat_sample = pd.concat([concat_sample, reagent.iloc[:, 0]])

        if pca_data:
            reagent_path_list = [self.des_data / chemical / ('nan.json') if sample[chemical] is np.nan
                                else self.des_data / chemical / (sample[chemical] + '.json')
                                for chemical in ChemicalList[3:]]
        else:
            reagent_path_list = [self.des_data / chemical / (sample[chemical] + '.json')
                                 for chemical in ChemicalList[3:]
                                 if sample[chemical] is not np.nan]

        for reagent_json in reagent_path_list:
            if os.path.exists(reagent_json):
                reagent = pd.read_json(reagent_json, orient='index')
                chemical_name = reagent_json.parent.name
                index = [chemical_name + '_' + str(index) for index in reagent.index]
                reagent.index = index
                concat_sample = pd.concat([concat_sample, reagent.iloc[:, 0]])
            else:
                not_find_file.add(reagent_json.parent.name + '/' + reagent_json.name)

        return concat_sample, not_find_file

    def _post_transform(self, concat_data):
        concat_data.fillna(0, inplace=True)
        concat_data = concat_data.loc[:, ~(concat_data == 0).all(axis=0)].copy()
        concat_data.drop(concat_data.columns[concat_data.dtypes == object], axis=1, inplace=True)

        columns_name = list(concat_data.columns)
        columns_name[0] = 'reference_id'
        concat_data.columns = columns_name

        temp = concat_data['yield']
        concat_data.drop(labels=['yield'], axis=1, inplace=True)
        concat_data.insert(concat_data.shape[-1], 'yield', temp)

        concat_data.iloc[:,1:-1] = (concat_data.iloc[:,1:-1] - concat_data.iloc[:,1:-1].mean()) \
                                   / concat_data.iloc[:,1:-1].std()

        return concat_data

    def output(self):
        self.database.to_csv(ModelData)
        with open(NotFindFile, 'w', encoding='utf-8') as f:
            for n_file in sorted(list(self.not_find)):
                print(n_file, file=f)


if __name__ == '__main__':
    pass
