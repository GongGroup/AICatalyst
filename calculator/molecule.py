from calculator.rbase import RMolecule
from common.fio import JsonIO
from common.logger import logger


class Ligand(RMolecule):
    _LigandInfo = JsonIO.read("ligand.json")

    def __init__(self, name=None, smiles=None, rmol=None):
        super(Ligand, self).__init__(rmol)
        self.name = name
        self.smiles = smiles

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.name}>"

    @staticmethod
    def from_strings(symbol, name=None, addH=True):
        if name is None:
            name = symbol

        smiles = Ligand._LigandInfo[symbol]
        rmol = RMolecule.from_smiles(smiles, addH)

        return Ligand(name, smiles, rmol)

    def get_unsaturated_atoms(self):
        unsaturated_atoms = []
        for atom in self.atoms:
            if atom.is_unsaturated:
                unsaturated_atoms.append(atom)

        if len(unsaturated_atoms) > 1:
            logger.warning(f"Unsaturated atoms more than 1, may failed")
        elif not len(unsaturated_atoms):
            logger.warning(f"Unsaturated atoms equal to 0, may be a saturated molecule")

        return unsaturated_atoms


class MCenter(RMolecule):
    def __init__(self, name=None, rmol=None):
        super(MCenter, self).__init__(rmol)
        self.name = name

    @staticmethod
    def from_strings(symbol, name=None, addH=False):
        if name is None:
            name = symbol

        smiles = symbol
        if "[" not in symbol:
            smiles = f"[{symbol}]"

        rmol = RMolecule.from_smiles(smiles, addH)

        return MCenter(name, rmol)


class Molecule(object):
    def __init__(self, center, ligands):
        self.center = center
        self.ligands = ligands

        self.rearrange()

    def __repr__(self):
        return f"{self.center}{self.ligands}"

    def rearrange(self):
        angle = 360.0 / len(self.ligands)
        print()


if __name__ == '__main__':
    OAc = Ligand.from_strings("OAc")
    Pd = MCenter.from_strings("Pd")
    uatoms = OAc.get_unsaturated_atoms()
    mol = Molecule(Pd, [OAc] * 2)

    print()
