import numpy as np
from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import AllChem

from common.constant import ElementInfo
from common.error import FileFormatError

# close the rdkit warning
RDLogger.DisableLog('rdApp.warning')


class RAtom(object):

    def __init__(self, ratom, rposition, mol):
        self._ratom = ratom
        self._rposition = rposition
        self._mol = mol

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.symbol}{self.explicit_valence} : {self.position}>"

    @property
    def is_unsaturated(self):
        if self.explicit_valence < ElementInfo[f'Element {self.symbol}']['valence']:
            return True
        return False

    @property
    def order(self):
        return self._ratom.GetIdx()

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

    def __init__(self, rmol):
        self._rmol = rmol

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.rsmiles}>"

    @property
    def rsmiles(self):
        return Chem.MolToSmiles(self._rmol)

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
    smiles = "CC([O-])=O"
    rmol = RMolecule._from_smiles(smiles)
    rmol = RMolecule(rmol)
    print()
