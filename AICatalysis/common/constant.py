import inspect
import os
import time
from pathlib import Path

from AICatalysis.common.file import JsonIO, YamlIO, QMIO

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
LogDir = os.path.join(CurrentDir, f"{RootDir}/logs")

# Directory constant here
ChemDir = Path("../../chemical")
DrawDir = Path("../../draw")
Calculator = Path("../calculator")

# file constant
FChemical = ChemDir / "chemical.json"
FFormula = ChemDir / "formula.json"
FReaxys = ChemDir / "reaxys.json"
FReactions = ChemDir / "reactions.json"
FReaxysYield = ChemDir / "opsin_reaxys.json"
FOpsinRecord = ChemDir / "opsin_record.json"
FRXconfig = ChemDir / "RX_config.json"
FRXDconfig = ChemDir / "RXD_config.json"
FReaxysXML = ChemDir / "reaxys_xml.xml"
FSucReaxys = ChemDir / "suc_reaxys.json"
Elements = Calculator / "element.yaml"
Angle = Calculator / "angle.dat"
FQM1 = Calculator / "QM1.dat"
FQM2 = Calculator / "QM2.dat"

# Variable constant
ChemInfo = {item['name']: {key: value for key, value in item.items() if key != 'name'}
            for item in JsonIO.read(FOpsinRecord)}
ElementInfo = YamlIO.read(Elements)
QM1 = QMIO.read1(FQM1)
QM2 = QMIO.read2(FQM2)
