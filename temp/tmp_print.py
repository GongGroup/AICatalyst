import json
from common.constant import FSucReaxys

import pandas as pd
from collections import Counter
import rdkit
from rdkit import DataStructs
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import PandasTools

with open(FSucReaxys, 'r') as f:
    doc = json.load(f)

