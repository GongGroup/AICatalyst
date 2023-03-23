import subprocess
from collections import Counter

from AICatalysis.calculator.gaussian import OUTFile
from AICatalysis.calculator.rbase import RMolecule
from AICatalysis.common.utils import flatten


class Descriptor(object):

    def __init__(self, name):
        self.name = name
        self._out_gaussian = OUTFile(name).read()
        self._rmol = self._convert_rdkit()
        self._convert_rdkit()

    def _convert_rdkit(self):
        process = subprocess.Popen(f"obabel -ig16 {self.name} -osdf", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        content = process.stdout.read().decode(encoding='utf-8')
        lines = content.splitlines()[:-4]
        with open("temp.mol", "w") as f:
            f.write("\n".join(lines))
        rmol, _ = RMolecule._from_mol_file("temp.mol")

        return RMolecule(rmol)

    @property
    def Weight(self):
        """
        molecular weight and average atomic weight

        """
        return self._rmol.weight

    @property
    def TAtoms(self):
        """
        total number of atoms in the molecule

        """
        return len(self._out_gaussian.input_atoms)

    @property
    def Atoms(self):
        """
        absolute numbers of atoms of certain chemical identity (C, H, O, N, F, etc.) in the molecule

        """
        symbols = [atom[0] for atom in self._out_gaussian.input_atoms]
        return Counter(symbols)

    @property
    def TBonds(self):
        """
        total number of bonds in the molecule

        """
        _bonds = self._rmol.bonds

        return len(_bonds)

    @property
    def Bonds(self):
        """
        absolute numbers of single, double, triple, aromatic or other bonds in the molecule

        """

        _bonds_type = [_bond['bond_type'].name for _bond in self._rmol.bonds]

        return Counter(_bonds_type)

    @property
    def Groups(self):
        """
        absolute and relative numbers of certain chemical groups and functionalities in the molecule

        """
        _groups = self._rmol.groups

        return list(Counter(flatten(_groups)).keys())

    @property
    def Rings(self):
        """
        total number of rings, number of rings divided by the total number of atoms

        """

        return len(self._rmol.rings)

    @property
    def AromaticRings(self):
        """
        total and relative number of 6-atoms aromatic rings

        """

        return len(self._rmol.aromatic_rings)


if __name__ == '__main__':
    m_descriptor = Descriptor("../database/chemical-gjf/Ac2O.out")
    pass
