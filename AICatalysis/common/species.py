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

TransMetalElement = ['–', 'Pd', 'catalyst', 'metal']

LigandType = ["1", "ligand", "Ligand", "Ph3P", "Bu(Ad)2P", "Cy3P", "(o-tolyl)3P", "XPhos", "dppb", "dppe", "dppp",
              "BINAP", "Xantphos", "dppf", "DPEphos", '–', 'PPh3', "P(o-tolyl)3", "P(p-anisyl)3", "PCy3", "BuPAd2",
              "P(o-Tol)3"]

SolventType = ["1", "mL", "dioxane", "toluene", "ACN", "DMF", "NMP", "MeCN", "DMSO", "Toluene", "THF",
               "Benzene", "Xylene", "DCE", "PhCH3", "4-BQ", "o-xylene", "CH3CN", "PhCl"]

BaseType = ["base", "Na2CO3", "Cs2CO3", "TEA", "DIEA", "DBU", "Et3N", "DIPEA", "TMEDA", "NEt3", "K2CO3", "K3PO4"]

AdditiveType = ["", "additive", "MI", "TBAB", "TBAC", "TBAI", ]

OxidantType = ["oxidant"]

AcidType = ["acid", "HCOOH", "HCO2H"]


class ChemFormula(object):
    mapping = JsonIO.read(FFormula)
    name = FormulaDescriptor('name')

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return f"<ChemFormula [{self.name}]>"

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
    name = ReagentDescriptor('name', MetalElement + MetalElementName + [item.lower() for item in MetalElementName])

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    def __repr__(self):
        return f"<MCatalyst [{self.name}]>"

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
            match = re.search(r'(.*)\s\(([0-9]+\.?[0-9]*\s?(mol)?(mmol)?(equiv)?.*)\)', self.name)
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
            match = re.search(r'(.*)\s\(([0-9]+\.[0-9]+\s?(mol)?(mmol)?(equiv)?)\)', self.name)
            self.formula, self.content = match.groups()[0:2]
        else:
            self.formula = self.name


class TransMetal(object):
    name = ReagentDescriptor('name', TransMetalElement)

    def __init__(self, name):
        self.name = name
        self.formula = None
        self.content = None

    @staticmethod
    def is_or_not(name):
        try:
            TransMetal(name)
        except ValueError:
            return False
        else:
            return True

    def parse(self):
        if "mol" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.*[0-9]*\s?m?mol.*)\)', self.name)
            self.formula, self.content = match.groups()
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
        patten2 = re.compile("(.*)\s(\([0-9]+\))")
        if "mol" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.*[0-9]*\s?m?mol.*)\)', self.name)
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
        if "mL" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.*[0-9]*\s?mL)\)', self.name)
            self.formula, self.content = match.groups()
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
    name = SolDescriptor('name', ['RT', "°C"])

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
    name = GasDescriptor('name', ['N2', 'CO', 'gas'])

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
        if "MPa" in self.name:
            match = re.search(r'(.*)\s\(([0-9]+\.*[0-9]*\s?MPa)\)', self.name)
            self.formula, self.content = match.groups()
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
