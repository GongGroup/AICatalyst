from rdkit import Chem
from rdkit.Chem import AllChem


class RAtom(object):

    def __init__(self, atom):
        self._atom = atom

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.symbol}{self.total_degree}>"

    @property
    def symbol(self):
        return self._atom.GetSymbol()

    @property
    def atom_map_num(self):
        return self._atom.GetAtomMapNum()

    @property
    def atomic_num(self):
        return self._atom.GetAtomicNum()

    @property
    def formal_charge(self):
        return self._atom.GetFormalCharge()

    @property
    def bonds(self):
        return self._atom.GetBonds()

    @property
    def degree(self):
        return self._atom.GetDegree()

    @property
    def total_degree(self):
        return self._atom.GetTotalDegree()

    @property
    def total_valence(self):
        return self._atom.GetTotalValence()


class RMolecule(object):

    def __init__(self, mol):
        self._mol = mol

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.smiles}>"

    @property
    def smiles(self):
        return Chem.MolToSmiles(self._mol)

    @property
    def atoms(self):
        atoms = []
        for atom in self._mol.GetAtoms():
            atoms.append(RAtom(atom))
        return atoms

    @property
    def num_atoms(self):
        return self._mol.GetNumAtoms()

    @staticmethod
    def from_smiles(smiles):
        mol = Chem.MolFromSmiles(smiles)
        mol = AllChem.AddHs(mol)
        return RMolecule(mol)


if __name__ == '__main__':
    smiles = "CC([O-])=O"
    mol = RMolecule.from_smiles(smiles)
    print()
