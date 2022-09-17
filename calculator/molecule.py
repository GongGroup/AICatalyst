import subprocess
from typing import Optional

from calculator.base import Atom, Atoms
from calculator.rbase import RMolecule
from common.fio import JsonIO


class Ligand(object):
    _LigandInfo = JsonIO.read("ligand.json")

    def __init__(self, name=None, atoms: Optional[Atoms] = None, smiles=None):
        self.name = name
        self.atoms = atoms
        self.smiles = smiles

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.name}>"

    @property
    def rdkit_mol(self):
        if self.smiles is not None:
            return RMolecule.from_smiles(self.smiles)
        return None

    @staticmethod
    def from_strings(symbol, name=None):
        smiles = Ligand._LigandInfo[symbol]
        process = subprocess.Popen(f"obabel -:{smiles} --gen3d -oxyz", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        content = process.stdout.read().decode(encoding='utf-8')
        lines = content.splitlines()

        atoms = []
        for line in lines:
            items = line.split()
            if len(items) == 4:
                coords = list(map(float, items[1:]))
                atoms.append(Atom(items[0], cart_coord=coords))
        atoms = Atoms.from_list(atoms)
        if name is None:
            name = symbol
        return Ligand(name, atoms, smiles)


class Molecule(object):
    def __init__(self, center, ligands):
        self.center = center
        self.ligands = ligands


if __name__ == '__main__':
    # LigandInfo = JsonIO.read("ligand.json")
    # smiles = LigandInfo['OAc']
    # process = subprocess.Popen(f"obabel -:{smiles} --gen3d -oxyz", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # content = process.stdout.read().decode(encoding='utf-8')
    OAc = Ligand.from_strings("OAc")
    mol = Molecule("Pd", [OAc] * 2)
    # lines = content.splitlines()
    # print(lines)
    print()
