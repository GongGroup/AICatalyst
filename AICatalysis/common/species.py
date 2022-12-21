import logging
import re
import shutil

from AICatalysis.common.constant import FFormula, FChemical
from AICatalysis.common.descriptor import ReagentDescriptor, FormulaDescriptor, SolDescriptor, TimeDescriptor, \
    LigandDescriptor, GasDescriptor, AdditiveDescriptor, OxidantDescriptor
from AICatalysis.common.file import JsonIO, ftemp
from AICatalysis.database.ichem import IChemCrawler

logger = logging.getLogger(__name__)

MetalElement = [
    'Li', 'Be',
    'Na', 'Mg', 'Al',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',
    'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po',
    'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv',
    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
]

MetalElementName = [
    'Lithium', 'Beryllium',
    'Sodium', 'Magnesium', 'Aluminium',
    'Potassium', 'Calcium', 'Scandium', 'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron', 'Cobalt', 'Nickel',
    'Copper', 'Zinc', 'Gallium', 'Germanium', 'Arsenic',

    'Rubidium', 'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum', 'Technetium', 'Ruthenium', 'Rhodium',
    'Palladium', 'Silver', 'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',

    'Cesium', 'Barium', 'Lanthanum', 'Cerium', 'Praseodymium', 'Neodymium', 'Promethium', 'Samarium', 'Europium',
    'Gadolinium', 'Terbium', 'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium', 'Lutetium', 'Hafnium',
    'Tantalum', 'Tungsten', 'Rhenium', 'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury', 'Thallium', 'Lead',
    'Bismuth', 'Polonium',

    'Francium', 'Radium', 'Actinium', 'Thorium', 'Protactinium', 'Uranium', 'Neptunium', 'Plutonium', 'Americium',
    'Curium', 'Berkelium', 'Californium', 'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium', 'Lawrencium',

    'Rutherfordium', 'Dubnium', 'Seaborgium', 'Bohrium', 'Hassium', 'Meitnerium', 'Darmstadtium', 'Roentgenium',
    'Copernicium', 'Nihonium', 'Flerovium', 'Moscovium', 'Livermorium',
]

CarbonylCatalystType = ['–', 'Pd', 'catalyst', 'metal', 'NaI', "TBAI", "THAI", "KI", "LTA", "Pb(NO3)2", "CuCl2·2H2O",
                        "NiCl2·DME", "NiCl2", "NiBr2", 'Sulfur', "Cu(NO3)2", "CuCO3", "Cu(OTf)2", "CuI", "CuCl", "Cu2O"]

ReagentType = ["reagent", 'Co2(CO)8', 'Fe2(CO)9', 'Mo(CO)6', 'Cr(CO)6', 'carbonyl', 'dioxane', 'Rh4(CO)12',
               "NaBPh4"] + MetalElement + MetalElementName + [item.lower() for item in MetalElementName]

LigandType = ['–', "ligand", "Ligand", "L1", 'L2', 'L3', 'L4', 'L5', 'L6', "Ph3P", "Bu(Ad)2P", "Cy3P", "(o-tolyl)3P",
              "XPhos", 'XantPhos', 'xantphos', 'dcpp', 'Sphos', "dppb", "dppe", "dppp", "BINAP", "Xantphos", "dppf",
              "DPEphos", 'PPh3', "P(o-tolyl)3", "P(p-anisyl)3", "PCy3", "BuPAd2", "P(o-Tol)3", "Xphos", "PtBu3",
              "tBu2P(2-biphenyl)", "tBu2P(2;4;6-iPr3C6H2)"]

SolventType = ["mL", 'solvent', "ACN", "DMF", "NMP", "MeCN", "DMSO", "Toluene", "THF", 'PEG-400', "Ph", "DCE", "PhCH3",
               "4-BQ", 'DME', "dioxane", "toluene", 'Dioxane', "Glycol", "H2O", "CH2Cl2", "1_4-Dioxane", "Benzene",
               "Xylene", "o-xylene", "CH3CN", '1_4-dioxane', 'BMImBF4', 'BMImCl', 'BMImPF6', 'Methanol']

BaseType = ['–', "base", "Na2CO3", "Cs2CO3", "TEA", "DIEA", "DBU", "Et3N", "DIPEA", "TMEDA", "NEt3", "K2CO3", "K3PO4",
            "K2HPO4", "NaH2PO4", "KF", "NaOAc", "AgOAc", "nBuONa", "CsF", "dbu", "dabco", "B1", "B2", "B3", "KOAc",
            "Na3PO4", "Li2CO3", "NaHCO3", "KHCO3"]

AdditiveType = ["–", "additive", "MI", "TBAB", "TBAC", "TBAI", "NaI", "Bu4NI", "I2", "KI", "K2CO3", "KOH", "t-BuOK",
                "K3PO4", "AgOAc", "Ag2O", "K2S2O8", "Mn(OAc)2", 'NaHCO3', 'NaBF4', 'KOAc', 'NaOAc']

OxidantType = ["oxidant", "BQ", "Cu(OAc)2", "AgOAc", "BzOOBz", "TBHP", "Oxone", "CuBr2"]

AcidType = ["acid", "HCOOH", "HCO2H"]


class ChemFormula(object):
    mapping = JsonIO.read(FFormula)
    name = FormulaDescriptor('name')

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return f"<ChemFormula [{self.name}]>"

    @staticmethod
    def is_or_not(name):
        try:
            ChemFormula(name)
        except ValueError:
            return False
        else:
            return True

    def split(self):
        patten1 = re.compile("([A-Z][a-z])\(([A-Za-z]{3,4})\)[0-9]$")  # Pd(OAc)2, Pd(OPiv)2
        patten2 = re.compile("^(:?[A-Z][a-z]|K)([A-Za-z]{1,3})")  # NaI, NaOAc, KI
        patten3 = re.compile("(:?[A-Z][a-z]|K)[0-9]([A-Za-z]{2}[0-9])")  # Na2CO3, K2CO3
        patten4 = re.compile("([A-Z][a-z])[0-9]\(([A-Za-z]{3})\)[0-9]")  # Pd3(dba)2
        patten5 = re.compile("([A-Z][a-z])([A-Za-z]{2})[0-9]")  # PdCl2
        patten6 = re.compile("([A-Z][a-z])\(([A-Za-z]{3}[0-9])\)[0-9]")  # Pd(PPh3)4
        patten7 = re.compile("([A-Z][a-z][0-9]N)([A-Za-z]{1})")  # Bu4NI
        patten8 = re.compile("(H)(COOH)")  # HCOOH
        patten9 = re.compile("(K)[0-9](S2O8)")  # K2S2O8
        patten10 = re.compile("(:?t-BuO|nBuO)(:?K|Na)")  # t-BuOK
        if re.search(patten1, self.name) is not None:
            match = re.search(patten1, self.name)
            return match.groups()
        elif re.search(patten2, self.name) is not None:
            match = re.search(patten2, self.name)
            return match.groups()
        elif re.search(patten3, self.name) is not None:
            match = re.search(patten3, self.name)
            return match.groups()
        elif re.search(patten4, self.name) is not None:
            match = re.search(patten4, self.name)
            return match.groups()
        elif re.search(patten5, self.name) is not None:
            match = re.search(patten5, self.name)
            return match.groups()
        elif re.search(patten6, self.name) is not None:
            match = re.search(patten6, self.name)
            return match.groups()
        elif re.search(patten7, self.name) is not None:
            match = re.search(patten7, self.name)
            return match.groups()
        elif re.search(patten8, self.name) is not None:
            match = re.search(patten8, self.name)
            return match.groups()
        elif re.search(patten9, self.name) is not None:
            match = re.search(patten9, self.name)
            return match.groups()
        elif re.search(patten10, self.name) is not None:
            match = re.search(patten10, self.name)
            return match.groups()[::-1]
        else:
            return self.name

    @staticmethod
    def new(name):
        if name in list(ChemFormula.mapping.keys()):
            return ChemFormula(ChemFormula.mapping[name])
        else:
            formula = IChemCrawler().get_formula(name)
            if formula is not None:
                ChemFormula.mapping[name] = formula
                JsonIO.write(ChemFormula.mapping, ftemp(FFormula))
                shutil.move(ftemp(FFormula), FFormula)
                logger.info(f"{name}~{formula} mapping stored in database")
                return ChemFormula(formula)
            else:
                raise RuntimeError("Create ChemFormula instance error")


class Reagent(object):
    name = ReagentDescriptor('name', ReagentType)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            Reagent(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "mol" in self.name or "equiv" in self.name:
            match = re.search(r'(.*)\s?\((.*?\s?[0-9]+\.?[0-9]*\s?(mol)?(mmol)?(equiv)?.*)\)', self.name)
            if match is None:
                self.formula = self.name
            else:
                self.formula, self.content = match.groups()[0:2]
        else:
            self.formula = self.name


class Acid(object):
    name = ReagentDescriptor('name', AcidType)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            Acid(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "mol" in self.name or "equiv" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.[0-9]+\s?(mol)?(mmol)?(equiv)?.*)\)', self.name)
            self.formula, self.content = match.groups()[0:2]
        elif ":" in self.name:
            match = re.search(r'(.*)\s?\(([0-9](:[0-9])+)\)', self.name)
            if match is None:
                self.formula = self.name
            else:
                self.formula, self.content = match.groups()[0:2]
        else:
            self.formula = self.name


class Base(object):
    name = ReagentDescriptor('name', BaseType)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            Base(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "mol" in self.name or "equiv" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.?[0-9]*\s?(mol)?(mmol)?(equiv)?.*)\)', self.name)
            self.formula, self.content = match.groups()[0:2]
        else:
            self.formula = self.name


class Additive(object):
    name = AdditiveDescriptor('name', AdditiveType)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            Additive(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "mol" in self.name or "equiv" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.?[0-9]*\s?(mol)?(mmol)?(equiv)?.*)\)', self.name)
            self.formula, self.content = match.groups()[0:2]
        else:
            self.formula = self.name


class Oxidant(object):
    name = OxidantDescriptor('name', OxidantType)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            Oxidant(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "mol" in self.name or "equiv" in self.name:
            match = re.search(r'(.*)\s\(([0-9]*\.?[0-9]+\s?(mol)?(mmol)?(equiv\.?)?)\)', self.name)
            self.formula, self.content = match.groups()[0:2]
        else:
            self.formula = self.name


class CarbonylCatalyst(object):
    name = ReagentDescriptor('name', CarbonylCatalystType)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            CarbonylCatalyst(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "mol" in self.name or "%" in self.name or "equiv" in self.name:
            match1 = re.search(r'(.*)\s\(([0-9]+\.*[0-9]*\s?m?(mol)?%?(equiv)?.*)\)', self.name)
            if match1 is not None:
                self.formula, self.content = match1.groups()[:2]
            else:
                self.formula = self.name
        else:
            self.formula = self.name


class Ligand(object):
    name = LigandDescriptor('name', LigandType)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            Ligand(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "mol" in self.name:
            match = re.search(r'(.*)\s\((.*[0-9]+\.?[0-9]*\s?m?μ?mol.*)\)', self.name)
            self.formula, self.content = match.groups()
        else:
            self.formula = self.name


class Solvent(object):
    name = SolDescriptor('name', SolventType)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            Solvent(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "ml" in self.name or "mL" in self.name or "g" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.*[0-9]*\s?(m[lL])?g?.*?)\)', self.name)
            self.formula, self.content = match.groups()[:2]
        else:
            self.formula = self.name


class Time(object):
    name = TimeDescriptor('name')

    def __init__(self, name):
        self.name = name

    @staticmethod
    def is_or_not(name):
        try:
            Time(name)
        except ValueError:
            return False
        else:
            return True


class Temperature(object):
    name = SolDescriptor('name', ['RT', "°C", "room", 'r.t.'])

    def __init__(self, name):
        self.name = name

    @staticmethod
    def is_or_not(name):
        try:
            Temperature(name)
        except ValueError:
            return False
        else:
            return True


class Gas(object):
    name = GasDescriptor('name', ['N2', 'CO', 'gas', 'ambient', 'argon'])

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            Gas(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "MPa" in self.name or ":" in self.name or "atm" in self.name or "bar" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.?:?[0-9]*\s?(MPa)?(atm)?(bar)?.*?)\)', self.name)
            if match is not None:
                self.formula, self.content = match.groups()[0:2]
            else:
                self.formula = self.name
        else:
            self.formula = self.name


if __name__ == '__main__':
    chemicals = JsonIO.read(FChemical)
    catalysts = {chemical: Reagent.is_or_not(chemical) for chemical in chemicals[2500:]}
    # for chemical in chemicals:
    #     try:
    #         formula = ChemFormula(chemical)
    #     except ValueError:
    #         print(chemical)
    #         pass
    #     else:
    #         # print(formula.name)
    #         pass
    # ChemFormula("(C6H11)3P · HBF4")
    print()
