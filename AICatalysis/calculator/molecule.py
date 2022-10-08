import logging
import math
import os
import shutil
from pathlib import Path

import numpy as np

from AICatalysis.calculator.matrix import r_matrix
from AICatalysis.calculator.rbase import RMolecule
from AICatalysis.common.constant import QM1, QM2
from AICatalysis.common.error import StructureError
from AICatalysis.common.file import JsonIO
from AICatalysis.common.species import Metal
from AICatalysis.common.utils import get_combinations

logger = logging.getLogger(__name__)


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

        try:
            smiles = Ligand._LigandInfo[symbol]
        except KeyError:
            logger.error(f"`{symbol}` can not recognized as a valid ligand")
            exit(1)

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
                try:
                    translate_vector = np.array(QM1[anchor_atom_type]) - np.array(single_uatoms.position)
                except KeyError:
                    logger.error(f"atom_type: `{anchor_atom_type}` is not exist for translate")
                    exit(1)
                else:
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
                if len(neighbors) == 0:
                    unit_matrix = np.eye(3)
                    single_positions = np.dot(single_positions, unit_matrix)
                    ligand_rotation_positions.append(single_positions)
                    continue
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
                    raise TypeError("No atom pairs exist for rotate operation")
            ligand_positions = ligand_rotation_positions
            return ligand_positions

        def symmetry(ligand_uatoms_order, ligand_positions, unit_directions):
            """
            Symmetry Operation

            Args:
                ligand_positions (list[np.array]): positions before symmetry
                unit_directions (np.array): symmetry directions depends on CN (unit vector)

            Returns:
                ligand_positions (list[np.array]): positions after symmetry
            """
            final_positions = []
            distances = [np.sum(ligand[order] ** 2) ** 0.5
                         for ligand, order in zip(ligand_positions, ligand_uatoms_order)]
            directions = [distances[i] * unit_directions[i] for i in range(len(distances))]
            for single_position, vector_after, order in zip(ligand_positions, directions, ligand_uatoms_order):
                vector_before = single_position[order]
                rotation_matrix = r_matrix(vector_before, vector_after)
                single_position = np.dot(single_position, rotation_matrix)
                final_positions.append(single_position)
            return np.array(final_positions, dtype=object)

        if len(self.ligands) == 2:
            unit_directions = np.array([(0, 1, 0), (0, -1, 0)])
        elif len(self.ligands) == 4:
            unit_directions = np.array(
                [(0.26629634, 0.59716506, 0.75662418), (-0.78360941, -0.56477488, 0.2588158),
                 (-0.26830478, 0.58738621, -0.7635378), (0.77721406, -0.57470191, -0.25623428)])
        else:
            raise NotImplementedError(f"Number of ligands equal to {len(self.ligands)} is not supported now")

        def rotate_again(ligand_uatoms_order, ligand_positions, max_count=500, min_dist=1.5):
            """
            Rotate Operation around axis

            Args:
                ligand_uatoms_order (list): record the order of anchor atoms for each ligand
                ligand_positions (list[np.array]): positions before rotate

            Returns:
                ligand_positions (list[np.array]): positions after rotate_again
            """
            vector_axis = [ligand[order] for ligand, order in zip(ligand_positions, ligand_uatoms_order)]
            overlap = check_overlap(ligand_positions, min_dist=min_dist)
            count = 0
            while len(overlap):
                for pair in overlap:
                    angle = 0.
                    while angle <= 2 * math.pi:
                        logger.debug(f"angle = {angle}")
                        rotation_matrix = r_matrix(angle, vector_axis[pair[1]])
                        ligand_positions[pair[1]] = np.dot(ligand_positions[pair[1]], rotation_matrix)
                        if len(check_overlap(ligand_positions, min_dist=min_dist)):
                            angle += 0.03
                        else:
                            break
                    if angle > 2 * math.pi:
                        logger.warning(f"Step {count}: Rotation have performed one cycle, still overlap")
                        count += 1
                if count >= 10:
                    min_dist -= 0.01
                    logger.warning(f"Step {count}: --> Failed exceed 10 times, decrease min_dist to {min_dist}")
                if count >= max_count:
                    logger.error(f"Step {count}: --> Exceed max iteration ({max_count}), please check!")
                    exit(1)
                overlap = check_overlap(ligand_positions, min_dist=min_dist)
            return ligand_positions

        def check_overlap(ligand_positions, min_dist=1.5):
            length = len(ligand_positions)
            combinations = get_combinations(range(length), range(length))
            overlap = []
            for item in combinations:
                position_1 = ligand_positions[item[0]]
                position_2 = ligand_positions[item[1]]
                for index, atom_position in enumerate(position_2):
                    distance = np.sum((position_1 - atom_position) ** 2, axis=1) ** 0.5
                    min_distance = np.min(distance)
                    if min_distance <= min_dist:
                        overlap.append(item)
                        logger.debug(
                            f"The minium distance of ligand_{item[0]} and ligand_{item[1]} is {min_distance} (Atom {index})")
            return set(overlap)

        # ligand Atom list
        ligand_atoms = [ligand.atoms for ligand in self.ligands]  # atoms in each ligand
        ligand_uatoms = sum([ligand.get_unsaturated_atoms() for ligand in self.ligands], [])  # unsaturated_atoms
        ligand_uatoms_order = [uatom.order for uatom in ligand_uatoms]

        # ligand atom position && symbol array
        ligand_positions = [np.array([atom.position for atom in ligand]) for ligand in ligand_atoms]
        ligand_symbols = [[atom.symbol for atom in ligand] for ligand in ligand_atoms]

        if len(ligand_uatoms) != len(self.ligands):
            raise RuntimeError(
                f"Anchor atoms num({len(ligand_uatoms)}) is not equal ligands num({len(self.ligands)})")

        ligand_positions = translate(ligand_uatoms, ligand_positions)  # translate operation
        ligand_positions = rotate(ligand_uatoms, ligand_positions)  # rotation operation
        ligand_positions = symmetry(ligand_uatoms_order, ligand_positions, unit_directions)  # symmetry operation
        ligand_positions = rotate_again(ligand_uatoms_order, ligand_positions)

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
        logger.info(f"molecule has been writen to `{name}`")

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
    ligand = Ligand.from_strings("PPh3")
    ligand2 = Ligand.from_strings("Cl")
    center = MCenter.from_strings("Pd")
    mol = Molecule(center, [ligand, ligand, ligand, ligand], gfnff=False)
    mol.write_to_xyz()

    print()
