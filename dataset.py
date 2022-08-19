from pathlib import Path

import numpy as np
from torch.utils.data import Dataset

from reaxys import Reaxys

# Directory constant
ChemDir = Path("./chemical")

# file constant
FReaxys = ChemDir / "opsin_reaxys.json"


class ModelDataset(Dataset):
    def __init__(self, file):
        self.records = np.load(file, allow_pickle=True)

    def __getitem__(self, index):
        pass

    def __len__(self):
        pass


if __name__ == '__main__':
    dataset = ModelDataset("train_set_v1.npy")
    print()
