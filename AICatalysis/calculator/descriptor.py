import subprocess
from collections import Counter

from AICatalysis.calculator.gaussian import OUTFile
from AICatalysis.calculator.rbase import RMolecule


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
    def Groups(self):
        """
        absolute and relative numbers of certain chemical groups and functionalities in the molecule

        """
        return


if __name__ == '__main__':
    m_descriptor = Descriptor("../database/chemical-gjf/Ac2O.out")
    pass
