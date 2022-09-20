import os
import shutil
from pathlib import Path

import numpy as np

from calculator.matrix import r_matrix
from calculator.optimizer import Optimizer
from calculator.rbase import RMolecule, RAtom
from common.constant import ElementInfo, QM1, QM2  # , FFAngle
from common.file import JsonIO
from common.logger import logger


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

    def rearrange_old(self):
        if len(self.ligands) == 2:
            matrix = np.array([(0, 1, 0), (0, -1, 0)])
        else:
            raise NotImplementedError(f"Number of ligands equal to {len(self.ligands)} is not supported now")

        default_bonds = ElementInfo[f'Element {self.center.name}']['default_bonds']
        ligand_atoms = [ligand.atoms for ligand in self.ligands]
        ligand_uatoms = [ligand.get_unsaturated_atoms() for ligand in self.ligands]

        # obtain the anchor atoms of each ligand_position
        anchor_atoms = [(uatom.symbol, uatom.position, uatom.order) for ligand in ligand_uatoms for uatom in ligand]
        anchor_atoms_2 = [uatom for ligand in ligand_uatoms for uatom in ligand]

        if len(anchor_atoms) != len(self.ligands):
            raise RuntimeError(f"Anchor atoms num({len(anchor_atoms)}) is not equal ligands num({len(self.ligands)})")

        # calculate the delta vectors to the target position of the anchor atoms
        target_positions = []
        orders = []
        for anchor_atom, vector in zip(anchor_atoms, matrix):
            target_positions.append(default_bonds[f'Element {anchor_atom[0]}'] * vector)
            orders.append(anchor_atom[2])

        # search vector after rotated
        vector_rotated = []
        for ligand_anchor_atom, direction in zip(anchor_atoms_2, matrix):
            atom_type = f"{ligand_anchor_atom.symbol}_{ligand_anchor_atom.hybridization}"
            for parameter in QM1:
                if atom_type == parameter.atom_type:
                    vector_rotated.append(
                        np.array([parameter[1] * direction[1], parameter[2] * direction[1], parameter[3]]))

        # initial translate
        ligand_atoms_position = [np.array([atom.position for atom in ligand]) for ligand in ligand_atoms]
        vector_before = []
        for ligand_anchor_atom in anchor_atoms_2:
            vector_before.append(np.array(ligand_anchor_atom.position) - np.array(self.center.atoms[0].position))

        rotated_matrices = []
        for va, vb in zip(vector_before, vector_rotated):
            rotated_matrices.append(r_matrix(vb, va))

        rotated_atoms_position = []
        for atoms_position, rotated_matrix in zip(ligand_atoms_position, rotated_matrices):
            rotated_atoms_position.append(np.dot(atoms_position, rotated_matrix))

        deltas = [vector_rotated[index] - rotated_atoms_position[index][orders[index]] for index in
                  range(len(self.ligands))]

        rotated_atoms_position = [item + delta for item, delta in zip(rotated_atoms_position, deltas)]

        # calculation target angles
        # angle_parameters = []
        # for ligand_anchor_atom in anchor_atoms_2:  # anchor atom of each ligand
        #     atom2_type = f"{ligand_anchor_atom.symbol}_{ligand_anchor_atom.hybridization}"
        #     ligand_angle_parameters = []
        #     for neighbor_atom in ligand_anchor_atom.neighbors:  # each neighbour atom of ligand_anchor_atom
        #         atom3_type = f"{neighbor_atom.symbol}_{neighbor_atom.hybridization}"
        #         for parameter in FFAngle:  # each force field parameter
        #             if atom2_type == parameter.atom_2 and atom3_type == parameter.atom_3 \
        #                     and parameter.atom_1 == "*" or self.center.atoms[0].symbol == parameter.atom_1:
        #                 ligand_angle_parameters.append(
        #                     (self.center, ligand_anchor_atom.order, neighbor_atom.order, parameter.value))
        #
        #     angle_parameters.append(ligand_angle_parameters)

        self.optimized_position = []
        # for ligand_atom, target_position, order in zip(ligand_atoms, target_positions, orders):
        #     center_position = np.array(self.center.atoms[0].position)
        #     ligand_position = np.array([atom.position for atom in ligand_atom])
        #     ligand_symbol = [atom.symbol for atom in ligand_atom]
        #     optimizer = Optimizer(center=center_position, ligand=ligand_position)
        #     center_position, ligand_position = optimizer.optimize(target_position, order)
        #     for symbol, position in zip(ligand_symbol, ligand_position):
        #         self.optimized_position.append((symbol, position))
        center_position = np.array(self.center.atoms[0].position)
        for index in range(len(self.ligands)):
            for ligand_atom, rearrange_position in zip(ligand_atoms[index], rotated_atoms_position[index]):
                symbol = ligand_atom.symbol
                position = rearrange_position
                #   ligand_position = np.array([atom.position for atom in ligand_atom])
                #     ligand_symbol = [atom.symbol for atom in ligand_atom]
                #     optimizer = Optimizer(center=center_position, ligand=ligand_position)
                #     center_position, ligand_position = optimizer.optimize(target_position, order)
                #     for symbol, position in zip(ligand_symbol, ligand_position):
                self.optimized_position.append((symbol, position))

        self.optimized_position.append((self.center.atoms[0].symbol, center_position))

    def optimize(self):
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

        # translate operation
        for single_uatoms, single_positions in zip(ligand_uatoms, ligand_positions):
            anchor_atom_type = f"{single_uatoms.symbol}_{single_uatoms.hybridization}"
            translate_vector = np.array(QM1[anchor_atom_type]) - np.array(single_uatoms.position)
            single_positions += translate_vector

        # rotation operation
        ligand_rotation_positions = []
        for single_uatoms, single_positions in zip(ligand_uatoms, ligand_positions):
            anchor_atom_type = f"{single_uatoms.symbol}_{single_uatoms.hybridization}"
            neighbors = single_uatoms.neighbors
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

        # symmetry operation
        for single_position, single_direction in zip(ligand_positions, directions):
            single_position *= single_direction

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
        return RMolecule(rmol)


if __name__ == '__main__':
    # I = Ligand.from_strings("I")
    Cl = Ligand.from_strings("OAc")
    # CH2Cl = Ligand.from_file("CH2Cl.mol")
    center = MCenter.from_strings("Pd")
    mol = Molecule(center, [Cl, Cl], gfnff=True)
    mol.write_to_xyz()

    # mol = Molecule.from_file("OAc.mol")

    print()
