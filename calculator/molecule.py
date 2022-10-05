import os
import shutil
from pathlib import Path

import numpy as np

from calculator.matrix import r_matrix
from calculator.rbase import RMolecule
from common.constant import QM1, QM2
from common.error import StructureError
from common.file import JsonIO
from common.logger import logger
from common.species import Metal


class Ligand(RMolecule):
    _LigandInfo = JsonIO.read("ligand.json")

    def __init__(self, name=None, smiles=None, rmol=None):
        super(Ligand, self).__init__(rmol)
        self.name = name
        self.smiles = smiles

    def __repr__(self):
        return f"<{self.__class__.__name__} : {self.name}>"

    @staticmethod
    def from_strings(symbol, name=None, addHs=True):
        if name is None:
            name = symbol

        smiles = Ligand._LigandInfo[symbol]
        rmol = RMolecule._from_smiles(smiles, addHs)

        return Ligand(name, smiles, rmol)

    @staticmethod
    def from_file(file, format="mol"):
        if format != "mol":
            raise NotImplementedError(f"{format} file is not supported now")

        rmol, smiles = RMolecule._from_mol_file(file)
        return Ligand(Path(file).stem, smiles, rmol)

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
    def from_strings(symbol, name=None, addHs=False):
        if name is None:
            name = symbol

        smiles = symbol
        if "[" not in symbol:
            smiles = f"[{symbol}]"

        rmol = RMolecule._from_smiles(smiles, addHs)

        return MCenter(name, rmol)


class Molecule(object):
    def __init__(self, center, ligands, gfnff=True):
        self.center = center
        self.ligands = ligands
        self.optimized_position = None
        self.gfnff = gfnff

        self.optimize()

    def __repr__(self):
        return f"{self.center}{self.ligands}"

    def optimize(self):

        def translate(ligand_uatoms, ligand_positions):
            """
            Translate Operation

            Args:
                ligand_uatoms (list[RAtom]): unsaturated <RAtom object> in each ligand
                ligand_positions (list[np.array]): positions before translate

            Returns:
                ligand_positions (list[np.array]): positions after translate
            """

            for single_uatoms, single_positions in zip(ligand_uatoms, ligand_positions):
                anchor_atom_type = f"{single_uatoms.symbol}_{single_uatoms.hybridization}"
                translate_vector = np.array(QM1[anchor_atom_type]) - np.array(single_uatoms.position)
                single_positions += translate_vector
            return ligand_positions

        def rotate(ligand_uatoms, ligand_positions):
            """
            Rotate Operation

            Args:
                ligand_uatoms (list[RAtom]): unsaturated <RAtom object> in each ligand
                ligand_positions (list[np.array]): positions before rotate

            Returns:
                ligand_positions (list[np.array]): positions after rotate
            """

            ligand_rotation_positions = []
            for single_uatoms, single_positions in zip(ligand_uatoms, ligand_positions):
                anchor_atom_type = f"{single_uatoms.symbol}_{single_uatoms.hybridization}"
                neighbors = single_uatoms.neighbors  # get neighbors
                for neighbor in neighbors:
                    neighbor_atom_type = f"{neighbor.symbol}_{neighbor.hybridization}"
                    if QM2.get(anchor_atom_type, None) is not None and \
                            QM2[anchor_atom_type].get(neighbor_atom_type, None) is not None:
                        rotation_before = np.array(neighbor.position) - np.array(single_uatoms.position)
                        rotation_after = np.array(QM2[anchor_atom_type][neighbor_atom_type])
                        rotation_matrix = r_matrix(rotation_before, rotation_after)
                        single_positions = np.dot(single_positions, rotation_matrix)
                        ligand_rotation_positions.append(single_positions)
                        break
                else:
                    raise TypeError("The atom pairs is not exist in QM2 file")
            ligand_positions = ligand_rotation_positions
            return ligand_positions

        def symmetry(ligand_positions, directions):
            """
            Symmetry Operation

            Args:
                ligand_positions: (list[np.array]): positions before symmetry
                directions (np.array): symmetry directions depends on CN

            Returns:
                ligand_positions: (list[np.array]): positions after symmetry
            """

            for single_position, single_direction in zip(ligand_positions, directions):
                single_position *= single_direction
            return ligand_positions

        if len(self.ligands) == 2:
            directions = np.array([(1, 1, 1), (-1, -1, -1)])
        else:
            raise NotImplementedError(f"Number of ligands equal to {len(self.ligands)} is not supported now")

        # ligand Atom list
        ligand_atoms = [ligand.atoms for ligand in self.ligands]  # atoms in each ligand
        ligand_uatoms = sum([ligand.get_unsaturated_atoms() for ligand in self.ligands], [])  # unsaturated_atoms

        # ligand atom position && symbol array
        ligand_positions = [np.array([atom.position for atom in ligand]) for ligand in ligand_atoms]
        ligand_symbols = [[atom.symbol for atom in ligand] for ligand in ligand_atoms]

        if len(ligand_uatoms) != len(self.ligands):
            raise RuntimeError(
                f"Anchor atoms num({len(ligand_uatoms)}) is not equal ligands num({len(self.ligands)})")

        ligand_positions = translate(ligand_uatoms, ligand_positions)  # translate operation
        ligand_positions = rotate(ligand_uatoms, ligand_positions)  # rotation operation
        ligand_positions = symmetry(ligand_positions, directions)  # symmetry operation

        # add symbols
        self.optimized_position = []
        center_position = np.array(self.center.atoms[0].position)
        for index in range(len(self.ligands)):
            for symbol, position in zip(ligand_symbols[index], ligand_positions[index]):
                self.optimized_position.append((symbol, position))

        self.optimized_position.append((self.center.atoms[0].symbol, center_position))

    def write_to_xyz(self, name="molecule.xyz"):
        with open("temp.xyz", "w") as f:
            f.write(f"{len(self.optimized_position)} \n")
            f.write(f"\n")
            for item in self.optimized_position:
                f.write(f"{item[0]}\t {'    '.join(item[1].astype(str))} \n")

        if self.gfnff:
            logger.info("Use xtb as the post-opt backend")
            os.system(f"bash xtb.sh")
            shutil.move("xtbopt.xyz", name)
        else:
            shutil.move("temp.xyz", name)

    @staticmethod
    def from_file(file, format="mol"):
        if format != "mol":
            raise NotImplementedError(f"{format} file is not supported now")

        rmol, _ = RMolecule._from_mol_file(file)
        rmol = RMolecule(rmol)
        atoms = rmol.atoms
        metals = [atom for atom in atoms if Metal.is_metal(atom.symbol)]
        if len(metals) > 1:
            raise StructureError(f"The num of metal elements is `{len(metals)}`, should be `1`")
        center = MCenter.from_strings(metals[0].symbol)

        uatoms = [atom for atom in atoms if atom.is_unsaturated]
        fragments = [RMolecule(item) for item in rmol.mol_frags]
        return RMolecule(rmol)


if __name__ == '__main__':
    # I = Ligand.from_strings("I")
    Cl = Ligand.from_strings("OAc")
    # CH2Cl = Ligand.from_file("CH2Cl.mol")
    center = MCenter.from_strings("Pd")
    mol = Molecule(center, [Cl, Cl], gfnff=True)
    mol.write_to_xyz()

    # mol = Molecule.from_file("PdOAc.mol")

    print()
