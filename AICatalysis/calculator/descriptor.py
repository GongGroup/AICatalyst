import math
import subprocess
from collections import Counter

from rdkit.Chem import GraphDescriptors

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

        return RMolecule(rmol, remove_H=False)

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

    @property
    def WienerIndex(self):
        """
        Wiener Index (non-hydrogen atoms)

            W = \frac{1}{2} \sum_{(i,j)}^{N_{SA}}d_{ij}

        """
        rmol = RMolecule(self._rmol._rmol, remove_H=True)

        res = 0
        for i in range(rmol._rmol.GetNumAtoms()):
            for j in range(i + 1, rmol._rmol.GetNumAtoms()):
                res += rmol.distance_matrix[i][j]

        return res

    @property
    def ChiIndex(self):
        """
        Chi0 = \sum(\delta _i)^{-1/2}

        Chi0v = \sum(\delta _i^v)^{-1/2}, \delta _i^v = Z_i^v-b_i
                Z_i^v is number of valence electrons, b_i is number of H atoms bonded to i

        Chi0n is similar to Chi0v, but use nVal instead of valence

        Chi1 = \sum(\delta _i \delta _j)^{-1/2}

        References:
            https://doi.org/10.1002/9780470125793.ch9
        """

        rmol = RMolecule(self._rmol._rmol, remove_H=True)

        _chi_index = {'chi0': GraphDescriptors.Chi0(rmol._rmol),
                      'chi0v': GraphDescriptors.Chi0v(rmol._rmol),
                      'chi0n': GraphDescriptors.Chi0n(rmol._rmol),
                      'chi1': GraphDescriptors.Chi1(rmol._rmol),
                      'chi1v': GraphDescriptors.Chi1v(rmol._rmol),
                      'chi1n': GraphDescriptors.Chi1n(rmol._rmol),
                      'chi2n': GraphDescriptors.Chi2n(rmol._rmol),
                      'chi2v': GraphDescriptors.Chi2v(rmol._rmol),
                      'chi3n': GraphDescriptors.Chi3n(rmol._rmol),
                      'chi3v': GraphDescriptors.Chi3v(rmol._rmol),
                      'chi4n': GraphDescriptors.Chi4n(rmol._rmol),  # GraphDescriptors.ChiNn_(rmol._rmol, n)
                      'chi4v': GraphDescriptors.Chi4v(rmol._rmol)}  # GraphDescriptors.ChiNv_(rmol._rmol, n)

        return _chi_index

    @property
    def BalabanIndex(self):
        """
        J=\frac{m}{\gamma +1}\sum_{(i,j)\in E(G)}(D_iD_j)^{-1/2}

        """
        rmol = RMolecule(self._rmol._rmol, remove_H=True)

        return GraphDescriptors.BalabanJ(rmol._rmol)

    @property
    def KierShapeIndex(self):
        """
        Kappa1 = (A+\alpha)(A+\alpha-1)^2/(^1P+\alpha)^2
        Kappa2 = (A+\alpha-1)(A+\alpha-2)^2/(^2P+\alpha)^2
        Kappa3 = (A+\alpha-1)(A+\alpha-3)^2/(^3P+\alpha)^2      (A is odd)
        Kappa3 = (A+\alpha-2)(A+\alpha-3)^2/(^3P+\alpha)^2      (A is even)

        """

        rmol = RMolecule(self._rmol._rmol, remove_H=True)

        _kier_shape_index = {'kappa1': GraphDescriptors.Kappa1(rmol._rmol),
                             'kappa2': GraphDescriptors.Kappa2(rmol._rmol),
                             'kappa3': GraphDescriptors.Kappa3(rmol._rmol)}

        return _kier_shape_index

    @property
    def KierFlexIndex(self):
        """
        \Phi =\frac{^1\kappa ^2\kappa}{N_{SA}}

        """

        rmol = RMolecule(self._rmol._rmol, remove_H=True)
        N_SA = rmol.num_atoms

        return self.KierShapeIndex['kappa1'] * self.KierShapeIndex['kappa2'] / N_SA

    @property
    def IC(self):
        """
        Mean information content index

        IC = -\sum_{i}p_i\cdot  log_2\,p_i

        References:
            https://doi.org/10.1002/jps.2600730403

        """
        N = self._rmol.num_atoms
        p = [value / N for value in self._rmol.coordination_info.values()]
        _ic = -sum([item * math.log2(item) for item in p])

        return _ic


if __name__ == '__main__':
    m_descriptor = Descriptor("../database/chemical-gjf/Ac2O.out")
    pass