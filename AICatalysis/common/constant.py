import inspect
import os
import time
from pathlib import Path

PeriodicTable = [
    'H', 'D', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'P', 'S', 'Cl', 'Ar', 'Si',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
    'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
]

ChemicalList = ['sub1', 'sub2', 'product', 'Pd_precursor', 'ligand_P',
                 'base1', 'base2', 'add', 'oxidant', 'solvent1', 'solvent2', 'carbon_source',]
ReactionList = ['sub1', 'sub2', 'product']
ReagentList = ['Pd_precursor', 'ligand_P', 'base1', 'base2',
                'add', 'oxidant','solvent1', 'solvent2','metal_carbonyl']

DATE = time.strftime("%Y-%m-%d", time.localtime())

BLACK = "\033[1;30m"
RED = "\033[1;31m"
GREEN = "\033[1;32m"
YELLOW = "\033[1;33m"
BLUE = "\033[1;34m"
MAGENTA = "\033[1;35m"
CYAN = "\033[1;36m"
WHITE = "\033[1;37m"
BOLD = "\033[1m"
RESET = "\033[0m"

CurrentDir = os.path.dirname(os.path.abspath(os.path.realpath(inspect.getfile(inspect.currentframe()))))
SourceDir = os.path.dirname(CurrentDir)
RootDir = os.path.dirname(SourceDir)

# Directory constant here
Calculator = Path("../calculator")
DatabaseDir = Path("../../database")

ModelDataDir = DatabaseDir / "model_data"
ReactDataDir = DatabaseDir / "reaction_data"
StructDataDir = DatabaseDir / "struct_data"

DescriptorDir_ = StructDataDir / "descriptor_data"
TotalDesDir = DescriptorDir_ / "total_des"
RdkitDesDir = DescriptorDir_ / "rdkit_des"
GaussDesDir = DescriptorDir_ / "gaussian_des"
MultiwinDesDir = DescriptorDir_ / "multiwin_des"
DescriptorDataDir = DescriptorDir_ / "descriptor_database"
PCATotalDesDataDir = DescriptorDir_ / "pca_total_des"
PCARdkitDesDataDir = DescriptorDir_ / "pca_rdkit_des"
PCAMulDesDataDir = DescriptorDir_ / "pca_multiwin_des"

GaussianDataDir = StructDataDir / "gaussian_data"

GaussianInDataDir = GaussianDataDir / "input"
GaussianOutDataDir = GaussianDataDir / "output"

SmilesDir = ReactDataDir / "smiles_dir"

# file constant
ReactionFile = ReactDataDir / "reaction_dataset.csv"
NotFindFile = ReactDataDir / "not_find_file_list"

ModelData = ModelDataDir / "model_dataset.csv"

ReactSmilesFile = ReactDataDir / "smiles"
DefaultUpdateSmilesFile = SmilesDir / "smiles"

EAconfig = ModelDir / "ea_config.json"

Elements = Calculator / "element.yaml"

