import json
import os
import re
import subprocess
from functools import wraps
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from AICatalysis.calculator.gaussian import OUTFile, FCHKFile
from AICatalysis.calculator.rbase import RMolecule
from AICatalysis.calculator.smiles import SmilesFile
from AICatalysis.common.constant import TotalDesDir, RdkitDesDir, MultiwinDesDir, DescriptorDataDir, \
    ReactSmilesFile, PCAMulDesDataDir
from AICatalysis.common.file import JsonIO
from AICatalysis.common.utils import float_


def print_status(func):
    @wraps(func)
    def wrapper(self, *args, **kwargs):
        result = func(self, *args, **kwargs)
        print(f"{Path(self.out_name).stem} has saved\n")
        return result

    return wrapper


class Descriptor(object):
    des_type = None
    def __init__(self, out_name):
        self.out_name = out_name

    def calc_descriptor(self):
        pass

    @print_status
    def write(self):
        json_path = self.des_type / self.out_name.split('/')[-2] / (Path(self.out_name).stem + '.json')
        if not os.path.exists(json_path):
            JsonIO.write(self.calc_descriptor(), json_path)


class RdkitDescriptor(Descriptor):
    des_type = RdkitDesDir
    def __init__(self, out_name, smiles=None):
        super().__init__(out_name)
        self.smiles = smiles

    def _convert_rdkit(self):
        if self.smiles is None:
            process = subprocess.Popen(f"obabel -ig16 {self.out_name} -osdf", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            content = process.stdout.read().decode(encoding='utf-8')
            lines = content.splitlines()[:-4]
            with open("temp.mol", "w") as f:
                f.write("\n".join(lines))
            try:
                _rmol, _ = RMolecule._from_mol_file("temp.mol")
            except:
                return None
        else:
            _rmol = RMolecule._from_smiles(self.smiles)
        rmol = RMolecule(_rmol, remove_H=False)

        return rmol

    def calc_descriptor(self):
        self._rmol = self._convert_rdkit()
        if self._rmol is None:
            return None
        nms = [x[0] for x in Descriptors._descList]
        calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
        return dict(zip(nms, calc.CalcDescriptors(self._rmol._rmol)))


class SmilesDescriptor(Descriptor):
    def __init__(self, out_name, smiles=None):
        super().__init__(out_name)
        self.smiles = smiles


    def calc_descriptor(self):
        if self.smiles is None:
            process = subprocess.Popen(f"obabel -ig16 {self.out_name} -osdf", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            content = process.stdout.read().decode(encoding='utf-8')
            lines = content.splitlines()[:-4]
            with open("temp.mol", "w") as f:
                f.write("\n".join(lines))
            try:
                _rmol, _ = RMolecule._from_mol_file("temp.mol")
            except:
                return None
            smiles = Chem.MolToSmiles(_rmol._rmol)
        else:
            smiles = self.smiles

        return smiles

    def write(self):
        return self.calc_descriptor()

class GaussianDescriptor(Descriptor):

    def __init__(self, out_name):
        super().__init__(out_name)
        self._out_gaussian = OUTFile(out_name).read()

class MultiwinDescriptor(Descriptor):
    mod = ["WfnSurfESP", "WfnLength", "WfnSurfALIE", "OrbRelated", "EnergyRelated"]
    des_type = MultiwinDesDir

    def __init__(self, out_name):
        super().__init__(out_name)
        self.fchk_name = re.sub(r'\.out', '.fchk', re.sub('out/', 'fchk/',self.out_name))


    def calc_descriptor(self):
        self._out_gaussian = OUTFile(self.out_name).read()
        des = {}
        for mod_fun in self.mod:
            des.update(getattr(self, mod_fun)())
        return  des

    def WfnSurfESP(self):
        """
        Use Gaussian + Multiwfn method to calculate the molecule surface area && Volume

        References:
            http://sobereva.com/487, http://sobereva.com/102, http://sobereva.com/159

        """
        cal_log = "log"
        wfn_config = "config"
        _esp_min, _esp_max = None, None
        _volume = None  # (unit: Angstrom^3)
        _density = None #(unit: g/cm^3)
        _surf_area, _pos_surf_area, _neg_surf_area = None, None, None  # (unit: Angstrom^2)
        _elec_stat_pot, _pos_elec_stat_pot, _neg_elec_stat_pot = None, None, None # (unit: a. u.)
        _esp_variance, _pos_esp_variance, _neg_esp_variance = None, None, None # (unit: a. u. ^2)
        _balance_charge = None
        _mol_polarity_index = None # (unit: eV)
        _nonpolar_surf_area, _polar_surf_area  = None, None # (unit: Angstrom ^2)
        _nonpolar_surf_area_percent, _polar_surf_area_percent = None, None

        out_file = Path(self.out_name)
        fchk_file = self.fchk_name
        if not Path(fchk_file).exists():
            return None

        with open(wfn_config, "w") as f:
            f.writelines('12\n0\n-1\n-1\nq')

        os.system(f"Multiwfn.exe {fchk_file} < {wfn_config} > {cal_log}")

        with open(cal_log, "r") as f:
            _content = f.readlines()

        # os.remove(cal_log)
        # os.remove(wfn_config)

        for line in _content:
            if line.startswith(' Global surface minimum:'):
                _esp_min = float(line.split()[3])
            elif line.startswith(' Global surface maximum:'):
                _esp_max = float(line.split()[3])
            elif line.startswith(' Volume:'):
                try:
                    _volume = float(line.split()[4])
                except:
                    _volume = format(line.split()[3].split('(')[1])
            elif line.startswith(' Estimated density '):
                _density = float(line.split()[8])
            elif line.startswith(' Overall surface area'):
                try:
                    _surf_area = float(line.split()[6])
                except:
                    _surf_area = float(line.split()[5].split('(')[1])
            elif line.startswith(' Positive surface area:'):
                _pos_surf_area = float(line.split()[6])
            elif line.startswith(' Negative surface area:'):
                _neg_surf_area = float(line.split()[6])
            elif line.startswith(' Overall average value:'):
                _elec_stat_pot = float(line.split()[3])
            elif line.startswith(' Positive average value:'):
                _pos_elec_stat_pot = float(line.split()[3])
            elif line.startswith(' Negative average value:'):
                _neg_elec_stat_pot = float(line.split()[3])
            elif line.startswith(' Overall variance (sigma^2_tot):'):
                _esp_variance = float(line.split()[3])
            elif line.startswith(' Positive variance:'):
                _pos_esp_variance = float(line.split()[2])
            elif line.startswith(' Negative variance:'):
                _neg_esp_variance = float(line.split()[2])
            elif line.startswith(' Balance of charges (nu):'):
                _balance_charge = float(line.split()[4])
            elif line.startswith(' Molecular polarity index (MPI):'):
                _mol_polarity_index = float(line.split()[4])
            elif line.startswith(' Nonpolar surface area'):
                _nonpolar_surf_area = float(line.split()[7])
                try:
                    _nonpolar_surf_area_percent = float(line.split()[10])
                except:
                    _nonpolar_surf_area_percent = 100.0
            elif line.startswith(' Polar surface area'):
                _polar_surf_area = float(line.split()[7])
                _polar_surf_area_percent = float(line.split()[10])
                break

        return {'esp_min': _esp_min, 'esp_max': _esp_max,
                'volume': _volume, 'density': _density,
                'surf_area': _surf_area, 'pos_surf_area': _pos_surf_area, 'neg_surf_area': _neg_surf_area,
                'elec_stat_pot': _elec_stat_pot, 'pos_elec_stat_pot': _pos_elec_stat_pot, 'neg_elec_stat_pot': _neg_elec_stat_pot,
                'esp_variance': _esp_variance, 'pos_esp_variance': _pos_esp_variance, 'neg_esp_variance': _neg_esp_variance,
                'balance_charge': _balance_charge,
                'mol_polarity_index': _mol_polarity_index,
                'nonpolar_surf_area': _nonpolar_surf_area, 'polar_surf_area': _polar_surf_area,
                'nonpolar_surf_area_percent': _nonpolar_surf_area_percent, 'polar_surf_area_percent': _polar_surf_area_percent}

    def WfnLength(self):
        """
            References:
            http://sobereva.com/426, http://sobereva.com/190

        """
        cal_log = "log"
        wfn_config = "config"
        _length_x, _length_y, _length_z = None, None, None
        _radius = None

        out_file = Path(self.out_name)
        fchk_file = self.fchk_name

        if not Path(fchk_file).exists():
            return None

        with open(wfn_config, "w") as f:
            f.writelines("100\n21\nsize\n0\nq\n0\nq")

        os.system(f"Multiwfn.exe {fchk_file} < {wfn_config} > {cal_log}")

        with open(cal_log, "r") as f:
            _content = f.readlines()
        os.remove(cal_log)

        for line in _content:
            if line.startswith(' Radius of the system:'):
                _radius = float(line.split()[4])
            elif line.startswith(' Length of the three sides:'):
                _length_x, _length_y, _length_z = [float(x) for x in line.split()[5:8]]
                break

        return {'radius': _radius,
                'length_x': _length_x,
                'length_y': _length_y,
                'length_z': _length_z}

    def WfnSurfALIE(self):
        cal_log = "log"
        wfn_config = "config"
        _ALIE = None # (unit: a. u.)

        out_file = Path(self.out_name)
        fchk_file = self.fchk_name

        if not Path(fchk_file).exists():
            return None

        with open(wfn_config, "w") as f:
            f.writelines("12\n2\n2\n0\n-1\n-1\nq")

        os.system(f"Multiwfn.exe {fchk_file} < {wfn_config} > {cal_log}")

        with open(cal_log, "r") as f:
            _content = f.readlines()
        os.remove(cal_log)

        for line in _content:
            if line.startswith(' Average value:'):
                _ALIE = float(line.split()[2])
                break

        return {'ALIE': _ALIE}

    def OrbRelated(self):
        homo = self._out_gaussian.homo
        lumo = self._out_gaussian.lumo

        # calculate Homo, Lumo
        out_file = Path(self.out_name)
        fchk_file = self.fchk_name
        if not Path(fchk_file).exists():
            return {"HOMO": self._out_gaussian.homo,
                    "LUMO": self._out_gaussian.lumo,
                    "AbsHardness": (self._out_gaussian.lumo - self._out_gaussian.homo) / 2}

        # calculate Fukui-related
        fchk = FCHKFile(fchk_file).read()

        coeff_homo = fchk.coeff[self._out_gaussian.homo_index]
        coeff_lumo = fchk.coeff[self._out_gaussian.lumo_index]
        Fukui_nucleophilic = np.sum(coeff_homo ** 2) / (1 - homo)
        Fukui_electrophilic = np.sum(coeff_lumo ** 2) / (lumo + 10)
        Fukui_one_electron = np.sum(np.kron(coeff_homo, coeff_lumo)) / (lumo - homo)

        # calculate Mulliken Bond Order
        cal_log = "log"
        wfn_config = "config"

        with open(wfn_config, "w") as f:
            f.writelines("9\n4\n0\n0\nq")

        os.system(f"Multiwfn.exe {fchk_file} < {wfn_config} > {cal_log}")

        with open(cal_log, "r") as f:
            _content = f.readlines()
        os.remove(cal_log)

        lns, lne = -1, 0
        for index, line in enumerate(_content):
            if line.startswith(' Bond orders with absolute value'):
                lns = index
            elif line.startswith(' If outputting bond order matrix'):
                lne = index
                break

        bond_order = float_([line.split()[-1] for line in _content[lns + 1:lne - 1]])
        avg_bond_order = sum(bond_order) / len(bond_order)

        # calculate Total and Free Valence
        with open(wfn_config, "w") as f:
            f.writelines("9\n1\n0\n0\nq")

        os.system(f"Multiwfn.exe {fchk_file} < {wfn_config} > {cal_log}")

        with open(cal_log, "r") as f:
            _content = f.readlines()
        os.remove(cal_log)

        lns, lne = -1, 0
        for index, line in enumerate(_content):
            if line.startswith(' Total valences and free valences defined by Mayer'):
                lns = index
            elif line.startswith(' If outputting bond order matrix'):
                lne = index
                break
        valence = [(line.split()[-2], line.split()[-1]) for line in _content[lns + 1:lne - 1]]
        tot_valence, free_valence = list(map(list, zip(*valence)))
        avg_tot_valence, avg_free_valence = sum(float_(tot_valence)) / len(tot_valence), \
                                            sum(float_(free_valence)) / len(free_valence)

        return {"HOMO": self._out_gaussian.homo,
                "LUMO": self._out_gaussian.lumo,
                "AbsHardness": (self._out_gaussian.lumo - self._out_gaussian.homo) / 2,
                "FukuiNucleophilic": Fukui_nucleophilic,
                "FukuiElectrophilic": Fukui_electrophilic,
                "FukuiOneElectron": Fukui_one_electron,
                "AvgMullikenBondOrder": avg_bond_order,
                "TotValence": avg_tot_valence,
                "FreeValence": avg_free_valence}


    def EnergyRelated(self):
        """
        total number of atoms in the molecule

        """
        return {"TAtom": len(self._out_gaussian.input_atoms)} | self._out_gaussian.energy

class TotalDescriptor(Descriptor):
    des_type = TotalDesDir

    def __init__(self, out_name):
        super().__init__(out_name)

    def calc_descriptor(self):
        rdkit_des_path = RdkitDesDir / self.out_name.split('/')[-2] / (Path(self.out_name).stem + '.json')
        #如果没有对应的描述符文件，就生成一个
        if not os.path.exists(rdkit_des_path):
            reaction_smiles = SmilesFile(ReactSmilesFile)
            if Path(self.out_name).stem in reaction_smiles.file_name.values:
                smiles = reaction_smiles.smiles[reaction_smiles.file_name == Path(self.out_name).stem].values[0]
                rdkit_descriptor = RdkitDescriptor(self.out_name, smiles=smiles)
            else:
                rdkit_descriptor = RdkitDescriptor(self.out_name)
            rdkit_descriptor.write()
        multiwin_des_path = MultiwinDesDir / self.out_name.split('/')[-2] / (Path(self.out_name).stem + '.json')
        # 如果没有对应的描述符文件，就生成一个
        if not os.path.exists(multiwin_des_path):
            multiwin_descriptor = MultiwinDescriptor(self.out_name)
            multiwin_descriptor.write()
        try:
            des = JsonIO.read(rdkit_des_path) | JsonIO.read(multiwin_des_path)
        except TypeError:
            #其中一类描述符为空
            des = JsonIO.read(rdkit_des_path) if JsonIO.read(rdkit_des_path) else JsonIO.read(multiwin_des_path)
        return des

class DescriptorDatabase:
    des_type = None
    def __init__(self, chemical_type):
        self.chemical_des_dir = self.des_type / chemical_type
        self.data = self.read_json()
        self.pca_data = self.pca()

    def read_json(self):
        data_dict = {}

        for filename in os.listdir(self.chemical_des_dir):
            if filename.endswith('.json'):
                with open(os.path.join(self.chemical_des_dir, filename)) as f:
                    data_dict[os.path.splitext(filename)[0]] = json.load(f)

        try:
            df = pd.DataFrame.from_dict(data_dict, orient='index')
        except:
            df = pd.DataFrame(data_dict).transpose()
        df.loc['nan'] = np.zeros(df.shape[1])
        df.fillna(0, inplace=True)
        return df

    def pca(self, thresh=0.9):
        from sklearn.preprocessing import StandardScaler
        from sklearn.decomposition import PCA

        scaler = StandardScaler()
        self.data.iloc[:, :] = scaler.fit_transform(self.data.iloc[:, :])
        pca = PCA()
        pca.fit(self.data.values)
        cumsum = np.cumsum(pca.explained_variance_ratio_)
        d = np.argmax(cumsum >= thresh) + 1
        pca = PCA(n_components=2)
        return pd.DataFrame(pca.fit_transform(self.data), index=self.data.index)


    def concat(self):
        des_list = []
        des_files = [os.path.join(self.chemical_des_dir, file) for file in os.listdir(self.chemical_des_dir)]
        for des_file in des_files:
            des_list.append(JsonIO.read(des_file))
        return des_list

    def update(self):
        pass

    def write_csv(self):
        write_file = DescriptorDataDir / (self.chemical_des_dir.stem + '_' + self.des_type.stem + '.csv')
        self.pca_data.to_csv(write_file)

    def write_json(self):
        for index, row in self.pca_data.iterrows():
            json_data = row.to_json()
            with open(PCAMulDesDataDir / self.chemical_des_dir.stem / f"{index}.json", "w") as file:
                file.write(json_data)


class RdkitDescriptorDatabase(DescriptorDatabase):
    des_type = RdkitDesDir

    def __init__(self, chemical_type):
        super().__init__(chemical_type)


class MultiwinDescriptorDatabase(DescriptorDatabase):
    des_type = MultiwinDesDir

    def __init__(self, chemical_type):
        super().__init__(chemical_type)


class TotalDescriptorDatabase(DescriptorDatabase):
    des_type = TotalDesDir

    def __init__(self, chemical_type):
        super().__init__(chemical_type)


if __name__ == '__main__':
    pass