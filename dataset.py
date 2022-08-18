from pathlib import Path

from torch.utils.data import Dataset

from reaxys import Reaxys

# Directory constant
ChemDir = Path("./chemical")

# file constant
FReaxys = ChemDir / "opsin_reaxys.json"


class ModelDataset(Dataset):
    def __init__(self, records):
        self.records = records

    def __getitem__(self, index):
        pass

    def __len__(self):
        pass


if __name__ == '__main__':
    reaxys = Reaxys(FReaxys)
    dataset = ModelDataset(reaxys.records)
    print()
