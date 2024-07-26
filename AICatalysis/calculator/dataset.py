import pandas as pd
from sklearn.model_selection import train_test_split

from AICatalysis.common.constant import *
from AICatalysis.calculator.functool import ss_split


class MetaDataSet:
    reaction_type_list = ['paper_id', 'reaction_type', 'r1', 'r2', 'r', 'cycle']
    chemical_list = ['delta_sub1', 'delta_sub2', 'Pd_precursor', 'ligand_P',
                     'base1', 'base2', 'add', 'oxidant','solvent1', 'solvent2','carbon_source']
    reaction_list = ['delta_sub1', 'delta_sub2']
    reagent_list = ['Pd_precursor', 'ligand_P', 'base1', 'base2',
                    'add', 'oxidant','solvent1', 'solvent2','carbon_source']
    label_list = ['yield']

    def combine(self, list1: list):
        combine1_data = pd.concat([getattr(self, feature) for feature in list1], axis=1)
        list2 = [feature for feature in self.chemical_list if feature not in list1]
        combine2_data = pd.concat([getattr(self, feature) for feature in list2], axis=1)
        return combine1_data, combine2_data

    def all_split(self):
        return [getattr(self, chemical) for chemical in self.chemical_list]

    def _split(self, low=2, step=5):
        import numpy as np
        train_index = np.array([], dtype=int)
        unique_paper_ids = self.data['paper_id'].unique()

        for paper_id in unique_paper_ids:
            sample_index = self.data[self.data['paper_id'] == paper_id].sort_values(
                by='yield').index
            delete_index = [sample_index[0]] if sample_index.size < low else np.delete(sample_index,
                                                            np.arange(low, sample_index.size, step))
            train_index = np.concatenate((train_index, delete_index))
        np.random.seed(42)
        np.random.shuffle(train_index)

        test_index = np.delete(self.data.index, train_index)
        train_set = SubDataSet(self.data.loc[train_index, :])
        test_set = SubDataSet(self.data.loc[test_index, :])
        return train_set, test_set

    def split_by_label(self, test_size, label=None, random_state=42):
        if label in self.data.columns:
            if (self.data[label].value_counts() == 1).any():
                tmp_feature = label + "_cut"
                self.data[tmp_feature] = pd.cut(self.data[label], bins=10)
                train_set, test_set = ss_split(self.data, tmp_feature,
                                                     test_size=test_size, random_state=random_state, n_splits=1)
                self.data.drop(tmp_feature, axis=1, inplace=True)
                train_set.drop(tmp_feature, axis=1, inplace=True)
                test_set.drop(tmp_feature, axis=1, inplace=True)
            else:
                train_set, test_set = ss_split(self.data, label,
                                                     test_size=test_size, random_state=random_state, n_splits=1)
        elif label is None:
            train_set, test_set = train_test_split(self.data, test_size=test_size, random_state=random_state)
        else:
            exit(f'{label} is not find in feature')
        return SubDataSet(train_set), SubDataSet(test_set)

    @property  # X for input
    def x(self):
        try:
            return self.data.drop(self.reaction_type_list + self.label_list, axis=1)
        except:
            old_reaction_type_list = ['paper_id']
            return self.data.drop(old_reaction_type_list + self.label_list, axis=1)

    @property  # label
    def y(self):
        return self.yield_

    @property
    def paper_id(self):
        return self.data['paper_id']

    @property
    def yield_(self):
        return self.data['yield']

    @property
    def chemical(self):
        chemical_columns = [column for column in self._data.columns
                            for chemical in self.chemical_list if column.startswith(chemical)]
        return self._data[chemical_columns]

    @property
    def condition(self):
        condition_columns = self._data.columns.difference(self.chemical.columns)
        if 'paper_id' in self.data.columns:
            return self.data[condition_columns].drop(self.reaction_type_list, axis=1)
        else:
            tmp_col = self.reaction_type_list
            tmp_col.pop(0)
            return self.data[condition_columns].drop(tmp_col, axis=1)

    @property
    def reaction(self):
        reaction_columns = [column for column in self._data.columns
                            for reaction in self.reaction_list if column.startswith(reaction)]
        return self._data[reaction_columns]

    @property
    def reagent(self):
        reagent_columns = [column for column in self._data.columns
                            for reagent in self.reagent_list if column.startswith(reagent)]
        return self._data[reagent_columns]

    def __getattr__(self, name):
        if name in self.chemical_list:
            tmp_column = [column for column in self._data.columns
                          if column.startswith(name)]
            return self.data[tmp_column]
        elif hasattr(self.data, name):
            return getattr(self.data, name)
        else:
            raise AttributeError(f'\'{self.__class__.__name__}\' object has no attribute \'{name}\'')


class DataSet(MetaDataSet):
    def __init__(self, dir_name='dataset', file=ModelData):
        from pathlib import Path
        try:
            self.data = self.read_file(file)
        except AttributeError:
            self.data = self.read_file(Path(file))
        self._data = self.data.drop(['yield'], axis=1)
        self.dir_path = file.parent / dir_name

    def read_file(self, file):
        if file.suffix == '.csv':
            data = pd.read_csv(file)
        else:
            raise NotImplementedError(f'not support {file.suffix} file')
        return data

    def make_dir(self):
        if not os.path.lexists(self.dir_path):
            os.makedirs(self.dir_path)

    def drop(self, threshold):
        self.data.drop(self.data.index[self.data['yield'] <= threshold], inplace=True)
        self.data.reset_index(inplace=True)
        self._data = self.data.drop(['yield'], axis=1)

    def drop0(self):
        self.data.drop(self.data.index[self.data['yield'] == 0], inplace=True)
        self.data.reset_index(inplace=True)
        self._data = self.data.drop(['yield'], axis=1)


    def split(self, train_size, test_size, label=None, random_state=42):
        if label in self.data.columns:
            if (self.data[label].value_counts() == 1).any():
                tmp_feature = label + "_cut"
                self.data[tmp_feature] = pd.cut(self.data[label], bins=10)
                train_valid_set, test_set = ss_split(self.data, tmp_feature,
                                                     test_size=test_size, random_state=random_state, n_splits=1)
                test_size = 1 - train_size / (1 - test_size)
                train_set, valid_set = ss_split(train_valid_set, tmp_feature,
                                                test_size=test_size, random_state=random_state, n_splits=1)
                train_set, valid_set, test_set = [data.drop(tmp_feature, axis=1, inplace=True)
                                                  for data in [train_set, valid_set, test_set]]
            else:
                train_valid_set, test_set = ss_split(self.data, label,
                                                     test_size=test_size, random_state=random_state, n_splits=1)
                test_size = 1 - train_size / (1 - test_size)
                train_set, valid_set = ss_split(train_valid_set, label,
                                                test_size=test_size, random_state=random_state, n_splits=1)
        elif label is None:
            train_valid_set, test_set = train_test_split(self.data, test_size=test_size, random_state=random_state)
            train_size = train_size / (1 - test_size)
            train_set, valid_set = train_test_split(train_valid_set, train_size=train_size, random_state=random_state)
        else:
            exit(f'{label} is not find in feature')
        return SubDataSet(train_set), SubDataSet(valid_set), SubDataSet(test_set)


class SubDataSet(MetaDataSet):
    def __init__(self, data):
        self.data = data.copy()
        if 'yield' in self.data.columns:
            self._data = self.data.drop(['yield'], axis=1)
        else:
            self._data = self.data



if __name__ == '__main__':
    pass
