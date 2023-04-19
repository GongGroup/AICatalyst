import os
from collections import Counter

import numpy as np
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem, MolSurf, rdFreeSASA, Descriptors3D, rdPartialCharges
from rdkit.Chem import FragmentCatalog
from rdkit.Chem import RDConfig

from AICatalysis.common.constant import ElementInfo
from AICatalysis.common.error import FileFormatError

# close the rdkit warning
RDLogger.DisableLog('rdApp.warning')


class RAtom(object):
    searched = []

    def __init__(self, ratom, rposition, mol):
        self._ratom = ratom
        self._rposition = rposition
        self._mol = mol

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.symbol}{self.explicit_valence} : {self.position}>"

    def __sub__(self, other):
        return np.array(self._rposition - other._rposition)

    @property
    def is_unsaturated(self):
        try:
            if self.explicit_valence < ElementInfo[f'Element {self.symbol}']['valence']:
                return True
        except KeyError:
            print(self.symbol)
        return False

    @property
    def order(self):
        return self._ratom.GetIdx()

    @property
    def mass(self):
        return self._ratom.GetMass()

    @property
    def symbol(self):
        return self._ratom.GetSymbol()

    @property
    def position(self):
        return self._rposition.x, self._rposition.y, self._rposition.z

    @property
    def atom_map_num(self):
        return self._ratom.GetAtomMapNum()

    @property
    def atomic_num(self):
        return self._ratom.GetAtomicNum()

    @property
    def formal_charge(self):
        return self._ratom.GetFormalCharge()

    @property
    def mulliken_charge(self):
        try:
            return self._ratom.GetDoubleProp("MullikenCharge")
        except KeyError:
            return

    @mulliken_charge.setter
    def mulliken_charge(self, value: float):
        self._ratom.SetDoubleProp("MullikenCharge", value)

    @property
    def gasteiger_charge(self):
        return float(self._ratom.GetProp("_GasteigerCharge"))

    @property
    def _bonds(self):
        return self._ratom.GetBonds()

    @property
    def degree(self):
        return self._ratom.GetDegree()

    @property
    def total_degree(self):
        return self._ratom.GetTotalDegree()

    @property
    def hybridization(self):
        return self._ratom.GetHybridization()

    @property
    def neighbors(self):
        return [self._mol.atoms[atom.GetIdx()] for atom in self._ratom.GetNeighbors()]

    @property
    def coordination_info(self):
        """
        Used for calculating IC descriptor

        """
        return [atom.total_valence for atom in self.neighbors]

    @property
    def connected(self):
        """
        Obtain the all connected atoms, which can be seen as one ligand

        Returns:
            connection (list): store the atoms' orders with the connection

        """
        connection = []
        self._connected(connection)

        return connection

    def _connected(self, connection):
        for neighbor in self.neighbors:
            if neighbor.order not in connection:
                connection.append(neighbor.order)
                neighbor._connected(connection)

    def get_connected_without_idx(self, idx):
        """
        Obtain the connected atoms without one index

        Returns:
            connection (list): store the atoms' orders with the connection (include self)

        """
        connection = []
        self._get_connected_without_idx(connection, idx)

        return connection

    def _get_connected_without_idx(self, connection, idx):
        for neighbor in self.neighbors:
            if neighbor.order not in connection and neighbor.order not in idx:
                connection.append(neighbor.order)
                neighbor._get_connected_without_idx(connection, idx)

    @property
    def explicit_valence(self):
        return self._ratom.GetExplicitValence()

    @property
    def implicit_valence(self):
        return self._ratom.GetImplicitValence()

    @property
    def total_valence(self):
        return self._ratom.GetTotalValence()

    @staticmethod
    def from_sp(symbol, position):
        return RAtom(Chem.Atom(symbol), position, mol=None)


class RMolecule(object):

    def __init__(self, rmol, remove_H=False):
        if not remove_H:
            self._rmol = rmol
        else:
            self._rmol = Chem.RemoveHs(rmol)

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.rsmiles}>"

    @property
    def rsmiles(self):
        return Chem.MolToSmiles(self._rmol)

    @property
    def weight(self):
        return Chem.rdMolDescriptors.CalcExactMolWt(self._rmol)

    @property
    def atoms(self):
        return [RAtom(atom, self._rmol.GetConformer().GetAtomPosition(index), self) for index, atom in
                enumerate(self._rmol.GetAtoms())]

    @property
    def positions(self):
        return np.array([atom.position for atom in self.atoms])

    @property
    def mass_center(self):
        return np.mean(self.positions, axis=0)

    @property
    def num_atoms(self):
        return self._rmol.GetNumAtoms()

    @property
    def rfrags(self):
        return [rmol for rmol in Chem.GetMolFrags(self._rmol, asMols=True)]

    @property
    def distance_matrix_3d(self):
        return Chem.Get3DDistanceMatrix(self._rmol)

    @property
    def distance_matrix(self):
        return Chem.GetDistanceMatrix(self._rmol)

    @property
    def adjacency_matrix(self):
        return Chem.rdmolops.GetAdjacencyMatrix(self._rmol)

    @property
    def chiral_centers(self):
        return Chem.FindMolChiralCenters(self._rmol)

    @property
    def bonds(self):

        return [{"idx": bond.GetIdx(),
                 "bond_type": bond.GetBondType(),
                 "bond_type_as_double": bond.GetBondTypeAsDouble(),
                 "aromatic": bond.GetIsAromatic(),
                 "conjugated": bond.GetIsConjugated(),
                 "in_ring": bond.IsInRing(),
                 "degree": (self.atoms[bond.GetBeginAtomIdx()].degree * self.atoms[
                     bond.GetEndAtomIdx()].degree) ** -0.5,
                 "begin": bond.GetBeginAtomIdx(),
                 "end": bond.GetEndAtomIdx()} for bond in self._rmol.GetBonds()]

    @property
    def rings(self):
        ssr = Chem.GetSymmSSSR(self._rmol)
        return [list(ring) for ring in ssr]

    @property
    def aromatic_rings(self):
        ring_info = self._rmol.GetRingInfo()
        atoms_in_rings = ring_info.AtomRings()

        _aromatic_rings = []
        for ring in atoms_in_rings:
            aromatic_atom_in_ring = 0
            for atom_id in ring:
                atom = self._rmol.GetAtomWithIdx(atom_id)
                if atom.GetIsAromatic():
                    aromatic_atom_in_ring += 1
            if aromatic_atom_in_ring == len(ring):
                _aromatic_rings.append(ring)

        return _aromatic_rings

    @property
    def groups(self):
        fName = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
        fparams = FragmentCatalog.FragCatParams(1, 6, fName)
        fcat = FragmentCatalog.FragCatalog(fparams)
        fcgen = FragmentCatalog.FragCatGenerator()
        fcgen.AddFragsFromMol(self._rmol, fcat)
        mol_frags = [fcat.GetEntryDescription(index) for index in range(fcat.GetNumEntries())]

        _groups = []
        for index in range(len(mol_frags)):
            mol_frag_groups = list(fcat.GetEntryFuncGroupIds(index))
            _in_groups = []
            for fg in mol_frag_groups:
                funcgroup = fparams.GetFuncGroup(fg)
                _in_groups.append(funcgroup.GetProp('_Name'))
            _groups.append(_in_groups)

        return _groups

    @property
    def rotate_bonds(self):
        _rotate_bonds = []
        for bond in self.bonds:
            bond_type = bond['bond_type']
            conjugated = bond['conjugated']
            begin_atom = self.atoms[bond['begin']]
            end_atom = self.atoms[bond['end']]
            if bond_type == Chem.rdchem.BondType.SINGLE and not conjugated and begin_atom.symbol != "H" and end_atom.symbol != "H":
                _rotate_bonds.append(bond)
        return _rotate_bonds

    @property
    def fragments(self):
        """
        Cut molecule based on rotate bonds

        Returns:
            frags (list[tuple, list]): tuple represent the cut-bonds, while the list represent the fragment
        """
        frags = []
        for bond in self.rotate_bonds:
            begin_idx = bond['begin']
            end_idx = bond['end']
            frags.append(((begin_idx, end_idx), sorted(self.atoms[begin_idx].get_connected_without_idx([end_idx])),))
            frags.append(((begin_idx, end_idx), sorted(self.atoms[end_idx].get_connected_without_idx([begin_idx])),))
        return frags

    @property
    def coordination_info(self):
        """
        Used for calculating IC descriptor

        """
        _coord = []
        for atom in self.atoms:
            _coord.append(tuple(sorted(atom.coordination_info)))

        return Counter(tuple(_coord))

    @property
    def TPSA(self):
        return MolSurf.TPSA(self._rmol)

    @property
    def FreeSASA(self):
        radii = rdFreeSASA.classifyAtoms(self._rmol)
        return rdFreeSASA.CalcSASA(self._rmol, radii)

    @property
    def volume(self):
        return Chem.AllChem.ComputeMolVolume(self._rmol)

    @property
    def PMI1(self):
        return Descriptors3D.PMI1(self._rmol)

    @property
    def PMI2(self):
        return Descriptors3D.PMI2(self._rmol)

    @property
    def PMI3(self):
        return Descriptors3D.PMI3(self._rmol)

    def compute_gasteiger_charge(self):
        rdPartialCharges.ComputeGasteigerCharges(self._rmol)

    @staticmethod
    def _from_smiles(smiles, addHs=True):
        rmol = Chem.MolFromSmiles(smiles)
        if addHs:
            rmol = AllChem.AddHs(rmol)
        AllChem.EmbedMolecule(rmol)
        AllChem.MMFFOptimizeMolecule(rmol)
        return rmol

    @staticmethod
    def _from_mol_file(file, removeHs=False):
        rmol = Chem.MolFromMolFile(file, removeHs=removeHs)
        if rmol is None:
            raise FileFormatError(f"The format of {file} is not correct")
        return rmol, Chem.MolToSmiles(rmol)


if __name__ == '__main__':
    # smiles = "CC([O-])=O"
    # smiles = '[H]C([H])([H])C(=O)OC(=O)C([H])([H])[H]'
    smiles = 'C(C)(C)O'
    rmol = RMolecule._from_smiles(smiles)
    rmol = RMolecule(rmol, remove_H=False)
    print()
