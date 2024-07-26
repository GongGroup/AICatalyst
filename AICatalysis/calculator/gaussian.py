import subprocess
import re
import os
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

import numpy as np


class Gaussian(object):
    pass


class GJFFile(object):
    def __init__(self, keyword="opt freq=noraman nmr pop=nboread b3lyp/def2tzvp int(grid=ultrafine)"):
        self.keyword = keyword

    def read(self):
        pass

    def write(self, smiles: str, file_path, name):

        sdf_path = os.path.join(file_path, name + '.sdf')
        gjf_path = os.path.join(file_path, name + '.gjf')

        mol = Chem.MolFromSmiles(smiles)
        mol = AllChem.AddHs(mol)
        AllChem.EmbedMultipleConfs(mol, numConfs=50)
        AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=500)
        Chem.MolToMolFile(mol, sdf_path)
        os.system(f'obabel {sdf_path} -isdf -ogjf -O {gjf_path}')
        os.remove(sdf_path)

        with open(gjf_path, 'r', encoding='utf-8') as f_gjf:
            text = f_gjf.readlines()

            text[0] = '%nprocshared=48\n'
            text[1] = '%mem=20GB\n'
            text[2] = f'%chk={name}.chk\n'
            text[3] = f'# {self.keyword}\n'
            text[4] = '\n'
            text.insert(5, Path(sdf_path).stem + '\n')
            text.insert(6, '\n')
            text.append('$nbo bndidx $end\n')

        with open(gjf_path, 'w', encoding='utf-8') as f_gjf:
            f_gjf.writelines(text)


class OUTFile(object):
    def __init__(self, name="result.out"):
        self.name = name
        self.input_atoms = None
        self.mulliken_charge = None
        self.dipole_moment = None
        self.homo, self.lumo = None, None
        self.homo_index, self.lumo_index = None, None
        self.energy = None

    def read(self):
        with open(self.name, "r") as f:
            self._strings = f.readlines()

        # read input_atoms
        lns, lne = None, None
        for index, line in enumerate(self._strings):
            if line.startswith(" Symbolic Z-matrix"):
                lns = index
            elif line.startswith(" GradGrad") or line.startswith(" Add virtual"):
                lne = index
            if lns is not None and lne is not None:
                break

        _format = lambda x: [str(x[0]), float(x[1]), float(x[2]), float(x[3])]
        atom_lines = [l.split() for l in self._strings[lns + 2:lne - 2]]
        self.input_atoms = list(map(_format, atom_lines))

        # read Mulliken charges
        mulliken_index = [index for index, line in enumerate(self._strings)
                          if line.startswith(" Mulliken charges:") or
                          line.startswith(" Mulliken charges and spin densities:")]
        mulliken_charge = self._strings[mulliken_index[-1] + 2:mulliken_index[-1] + 2 + len(self.input_atoms)]
        self.mulliken_charge = [float(line.split()[2]) for line in mulliken_charge]

        # read Dipole moment
        dipole_index = [index for index, line in enumerate(self._strings) if line.startswith(" Dipole moment")][-1] + 1
        dipole_moment = list(map(float, self._strings[dipole_index].split()[1::2]))
        self.dipole_moment = {"X": dipole_moment[0], "Y": dipole_moment[1], "Z": dipole_moment[2],
                              "Tot": dipole_moment[3]}

        # read HOMO/LUMO
        try:
            orb_s_index = [index for index, line in enumerate(self._strings)
                           if line.startswith(" The electronic state")][-1] + 1
        except:
            orb_s_index = [index for index, line in enumerate(self._strings)
                           if line.startswith(" Unable to determine electronic state:")][-1] + 1
        orb_e_index = [index for index, line in enumerate(self._strings) if
                       line.startswith("          Condensed to atoms")][-1]
        orbital_occ = list(map(float,
                               [item for line in self._strings[orb_s_index:orb_e_index] if "occ." in line
                                for item in re.sub('-', ' ', line.split("--")[1]).split()]))
        orbital_virt = list(map(float,
                                [item for line in self._strings[orb_s_index:orb_e_index] if "virt." in line
                                 for item in line.split("--")[1].split()]))
        self.homo, self.homo_index = -orbital_occ[-1], len(orbital_occ) - 1
        self.lumo, self.lumo_index = orbital_virt[0], len(orbital_occ)

        # read energy
        out_s_index = [index for index, line in enumerate(self._strings) if
                       line.startswith(' ----------------------------------------------------------------------')][-1]
        out_e_index = [index for index, line in enumerate(self._strings) if line.endswith('@\n')][-1]
        energy = "".join(self._strings[out_s_index:out_e_index + 1]).replace('\n', '').replace(' ', '').split('\\')
        keywords = {'HF', 'ZeroPoint', 'Thermal', 'ETot', 'HTot', 'GTot'}
        energy = {line.split("=")[0]: float(line.split("=")[1]) for line in energy if line.split("=")[0] in keywords}
        self.energy = energy

        return self


class FCHKFile(object):
    def __init__(self, name="result.fchk"):
        self.name = name
        self.coeff = None

    def read(self):
        with open(self.name, "r") as f:
            self._strings = f.readlines()

        # read Alpha MO coefficients
        basis_num, lns, lne = None, -1, -1
        for index, line in enumerate(self._strings):
            if line.startswith("Number of basis functions"):
                basis_num = int(line.split()[-1])
            elif line.startswith("Alpha MO coefficients"):
                lns = index
            elif line.startswith("Beta MO coefficients"):
                lne = index
                break
            elif line.startswith("Orthonormal basis"):
                lne = index
                break
        coeff = list(map(float, [item for line in self._strings[lns + 1:lne] for item in line.split()]))
        coeff = np.array(coeff).reshape((basis_num, -1))
        self.coeff = coeff

        return self


if __name__ == '__main__':
    pass