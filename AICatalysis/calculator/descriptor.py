import math
import os
import subprocess
from collections import Counter
from pathlib import Path

from rdkit.Chem import GraphDescriptors

from AICatalysis.calculator.gaussian import OUTFile
from AICatalysis.calculator.rbase import RMolecule
from AICatalysis.common.utils import flatten


def charge_selector(charge):
    def inner(func):
        def wrapper(self, *args, **kargs):
            if charge == "mulliken":
                self._partial_charge = self.MullikenCharge
            elif charge == "gasteiger":
                self._partial_charge = self.GasteigerCharge

            return func(self, *args, **kargs)

        return wrapper

    return inner


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
        _rmol, _ = RMolecule._from_mol_file("temp.mol")
        rmol = RMolecule(_rmol, remove_H=False)

        mulliken_charge = self._out_gaussian.mulliken_charge
        for atom, charge in zip(rmol.atoms, mulliken_charge):
            atom.mulliken_charge = charge

        rmol.compute_gasteiger_charge()

        return rmol

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
    def InfoContentIndex(self):
        """
        information content index

        IC = -\sum_{i}p_i\cdot  log_2\,p_i
        SIC = IC / log_2\,N
        CIC = log_2\,N -IC
        BIC = IC / log_2\,q

        References:
            https://doi.org/10.1002/jps.2600730403

        """
        N = self._rmol.num_atoms  # Number of atoms
        q = len(self._rmol.bonds)  # Number of edges
        p = [value / N for value in self._rmol.coordination_info.values()]  # number of type-i atom / N (coordination)

        ic = -sum([item * math.log2(item) for item in p])  # Mean information content index
        sic = ic / math.log2(N)  # Structural information content index
        cic = math.log2(N) - ic  # Complementary information content index
        bic = ic / math.log2(q)  # Bonding information content index

        _info_content = {"ic": ic,
                         "sic": sic,
                         "cic": cic,
                         "bic": bic}

        return _info_content

    @property
    def TopologicalElectronicIndes(self):
        """
        T^E = \sum_{(i,j)}^{N_{SA}}\frac{|q_i-q_j|}{r_{ij}^2} (i = 1..N, j=i+1..N)

        References:
            https://doi.org/10.1016/S0021-9673(01)86894-1
        """

        _topo_elect = 0.
        for i in range(self._rmol.num_atoms):
            for j in range(i + 1, self._rmol.num_atoms):
                _topo_elect += math.fabs(self._rmol.atoms[i].mulliken_charge - self._rmol.atoms[j].mulliken_charge) / \
                               (self._rmol.distance_matrix_3d[i][j]) ** 2

        return _topo_elect

    @property
    def WfnSurfVolume(self):
        """
        Use Gaussian + Multiwfn method to calculate the molecule surface area && Volume

        References:
            http://sobereva.com/487, http://sobereva.com/102

        """
        cal_log = "log"
        _surf_area = None  # (unit: Angstrom^2)
        _volume = None  # (unit: Angstrom^3)

        out_file = Path(self.name)
        wfn_file = out_file.parent / (out_file.stem + ".wfn")
        if not Path(wfn_file).exists():
            return _surf_area

        os.system(f"bash multiwfn.sh {wfn_file.as_posix()} {cal_log}")

        with open(cal_log, "r") as f:
            _content = f.readlines()
        os.remove(cal_log)

        for line in _content:
            if line.startswith(' Volume:'):
                _volume = float(line.split()[4])

            if line.startswith(' Overall surface area'):
                _surf_area = float(line.split()[6])
                break

        return {'surface area': _surf_area,
                'volume': _volume}

    @property
    def SolAccMolSurfArea(self):
        """
        Solvent-accessible molecular surface area

        """
        return self._rmol.FreeSASA

    @property
    def Volume(self):
        return self._rmol.volume

    @property
    def GravitationalIndex(self):
        """
        G = \sum_{(i,j)}^{N_{SA}}\frac{m_im_j}{r_{ij}^2} (i = 1..N, j=i+1..N)

        """

        _grav = 0.
        for i in range(self._rmol.num_atoms):
            for j in range(i + 1, self._rmol.num_atoms):
                _grav += self._rmol.atoms[i].mass * self._rmol.atoms[j].mass / (
                    self._rmol.distance_matrix_3d[i][j]) ** 2

        return _grav

    @property
    def PMI(self):
        """
        Principal moments of inertia of a molecule

        I_k = \sum_{i}m_ir_{ik}^2

        """
        return {"PMI1": self._rmol.PMI1,
                "PMI2": self._rmol.PMI2,
                "PMI3": self._rmol.PMI3}

    @property
    def ShadowArea(self):
        pass
        return

    @property
    def GasteigerCharge(self):
        return [atom.gasteiger_charge for atom in self._rmol.atoms]

    @property
    def MullikenCharge(self):
        return [atom.mulliken_charge for atom in self._rmol.atoms]

    @property
    @charge_selector(charge="mulliken")
    def PolarityParam(self):

        Q_max, Q_min = max(self._partial_charge), min(self._partial_charge)
        Q_max_arg, Q_min_arg = self._partial_charge.index(Q_max), self._partial_charge.index(Q_min)
        R_mm = self._rmol.distance_matrix_3d[Q_max_arg][Q_min_arg]

        P = Q_max - Q_min
        P1 = P / R_mm
        P2 = P1 / R_mm

        return {"P": P, "P1": P1, "P2": P2}

    @property
    def DipoleMoment(self):
        return self._out_gaussian.dipole_moment


if __name__ == '__main__':
    m_descriptor = Descriptor("../database/chemical-gjf/CH3OH.out")
    print(m_descriptor.PolarityParam)
    pass
