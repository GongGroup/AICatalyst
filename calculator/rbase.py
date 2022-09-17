from rdkit import Chem
from rdkit.Chem import AllChem

from common.constant import ElementInfo


class RAtom(object):

    def __init__(self, ratom, rposition):
        self._ratom = ratom
        self._rposition = rposition

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.symbol}{self.total_valence} : {self.position}>"

    @property
    def is_unsaturated(self):
        if self.total_valence < ElementInfo[f'Element {self.symbol}']['valence']:
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
    def explicit_valence(self):
        return self._ratom.GetExplicitValence()

    @property
    def implicit_valence(self):
        return self._ratom.GetImplicitValence()

    @property
    def total_valence(self):
        return self._ratom.GetTotalValence()


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
        return [RAtom(atom, self._rmol.GetConformer().GetAtomPosition(index)) for index, atom in
                enumerate(self._rmol.GetAtoms())]

    @property
    def num_atoms(self):
        return self._rmol.GetNumAtoms()

    @staticmethod
    def from_smiles(smiles, addH=True):
        rmol = Chem.MolFromSmiles(smiles)
        if addH:
            rmol = AllChem.AddHs(rmol)
        AllChem.EmbedMolecule(rmol)
        AllChem.MMFFOptimizeMolecule(rmol)
        return rmol


if __name__ == '__main__':
    smiles = "CC([O-])=O"
    rmol = RMolecule.from_smiles(smiles)
    rmol = RMolecule(rmol)
    print()
