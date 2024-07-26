import pandas as pd


class SmilesFile:
    header = ['type_id', 'smiles', 'name']

    def __init__(self, file):
        self.file = file
        self.data = self.read(self.file)

    @classmethod
    def read(cls, file):
        data = pd.read_csv(file, delimiter=' ', names=cls.header)
        data['file_name'] =[str(id) + '-' + str(name)
                            for id, name in zip(data.type_id, data.name)]
        return data

    def concat(self, smiles):
        if isinstance(smiles, pd.DataFrame):
            self.data = pd.concat([self.data, smiles])
            self.data.drop_duplicates(inplace=True)
        elif isinstance(smiles, SmilesFile):
            self.data = pd.concat([self.data, smiles.data])
            self.data.drop_duplicates(inplace=True)
        else:
            raise NotImplementedError

    def save(self, file=None):
        self.data.drop('file_name', axis=1, inplace=True)
        if file is None:
            self.data.to_csv(self.file, sep=' ', header=0, index=0)
        else:
            self.data.to_csv(file, sep=' ', header=0, index=0)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        return self.data.loc[item].values

    def __setitem__(self, item, value):
        if len(value) != len(self.header):
            pass
        else:
            self.data.loc[item] = value

    def __delitem__(self, item):
        self.data.drop(item, axis=0)

    def __getattr__(self, name):
        if hasattr(self.data, name):
            return getattr(self.data, name)
        else:
            raise AttributeError(f'\'{self.__class__.__name__}\' object has no attribute \'{name}\'')