import logging
import shutil

from AICatalysis.common.constant import FFormula, FChemical
from AICatalysis.common.descriptor import CatalystDescriptor, FormulaDescriptor
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

TransMetalElement = ['Pd']


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


class Metal(object):
    name = CatalystDescriptor('name', MetalElement + MetalElementName + [item.lower() for item in MetalElementName])

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return f"<MCatalyst [{self.name}]>"

    @staticmethod
    def is_metal(name):
        try:
            Metal(name)
        except ValueError:
            return False
        else:
            return True


class TransMetal(object):
    name = CatalystDescriptor('name', TransMetalElement)

    def __init__(self, name):
        self.name = name

    @staticmethod
    def is_transmetal(name):
        try:
            TransMetal(name)
        except ValueError:
            return False
        else:
            return True


if __name__ == '__main__':
    chemicals = JsonIO.read(FChemical)
    catalysts = {chemical: Metal.is_metal(chemical) for chemical in chemicals[2500:]}
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
