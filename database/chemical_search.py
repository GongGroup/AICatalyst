import os
import pprint
from collections import Counter

from IPython.display import display
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

#import common
from common.constant import FOpsinRecord
from common.constant import FSucReaxys
from common.constant import FReaxys
from common.constant import FChemical
from common.fio import JsonIO
from common.logger import logger


class FeatureError(Exception):
    pass


class SucReaxys:
    Features = ['reagent', 'solvent', 'product', 'reactant']


    def __init__(self, chemical=None, feature=None):
        self._index_list = []
        self.chemical = chemical
        self.feature = feature
        self._pd = pd.DataFrame()
        self.ChemInfo = {item['name']: {key: value for key, value in item.items() if key != 'name'}
                        for item in JsonIO.read(FOpsinRecord)}
        self._initialize()

    def _initialize(self):
        if not os.path.exists(FSucReaxys):
            self._write_suc_reaxys()
        self.reaxys_doc = self._pd.from_records(JsonIO.read(FSucReaxys))

    def forward(self):
        self._write_suc_reaxys()
        logger.info("Successfully generate suc_reaxys")

    def _write_suc_reaxys(self):
        self.all_chemistry = JsonIO.read(FChemical)
        self.origin_reaxys = JsonIO.read(FReaxys)
        reaxys_list = []
        for index, reaction in enumerate(self.origin_reaxys):
            if not reaction['yield']:
                continue
            tmp_chemistry_list = [self.ChemInfo[tmp_chemistry]['status'] == 'SUCCESS'
                                  for feature in self.Features
                                  for tmp_chemistry in reaction[feature]]
            if sum(tmp_chemistry_list) == len(tmp_chemistry_list):
                reaxys_list.append(reaction)
        JsonIO.write(reaxys_list, FSucReaxys)

    def search(self, chemical=None, feature=None, show_reagent=False,
               save_img=False, show_inline=True):
        if chemical:
            self.chemical = chemical
        if feature:
            self.feature = feature
        status = self.ChemInfo[self.chemical]['status']
        if not self.feature in self.Features:
            raise FeatureError('Feature illegal!')
        if status == 'FAILURE':
            print("{} is not found".format(self.chemical))
        else:
            self._extract_reaxys()
            if not self._index_list:
                print("{} is not found in {}".format(self.chemical, self.feature))
            elif save_img:
                self._save_png(show_reagent)
            else:
                for img in self._show_reaction(show_reagent, show_inline):
                    if img:
                        try:
                            img.show()
                        except:
                            display(img)
                    else:
                        pass
                    exit = input('Press any key to show next, q to quit:')
                    if exit.startswith('q'):
                        break

    def _extract_reaxys(self):
        self._index_list = []
        for index, chemical_list in zip(self.reaxys_doc.index, self.reaxys_doc[self.feature]):
            if self.chemical in chemical_list:
                self._index_list.append(index)

    def Name2Smart(self, chemical_list):
        if len(chemical_list) > 0:
            chemical_mol = [Chem.MolFromSmiles(self.ChemInfo[rea]['smiles']) for rea in chemical_list]
            return '.'.join([Chem.MolToSmarts(rea) for rea in chemical_mol])
        else:
            return None

    def _show_reaction(self, show_reagent, show_inline):
        for i, index in enumerate(self._index_list):
            reactants = self.reaxys_doc.loc[index, 'reactant']
            products = self.reaxys_doc.loc[index, 'product']
            try:
                reactants_smart = self.Name2Smart(reactants)
                products_smart = self.Name2Smart(products)
                if show_reagent:
                    solvents = self.reaxys_doc.loc[index, 'solvent']
                    reagents = self.reaxys_doc.loc[index, 'reagent']
                    reagents.extend(solvents)
                    sol_rea_smart = self.Name2Smart(reagents)
                    reaction_smart = '>'.join([reactants_smart, sol_rea_smart, products_smart])
                else:
                    reaction_smart = '>>'.join([reactants_smart, products_smart])
                reaction = AllChem.ReactionFromSmarts(reaction_smart)
                if show_inline:
                    img = reaction
                else:
                    img = Draw.ReactionToImage(reaction)
                print('Find {} reactions, {} of {}:'.format(
                    len(self._index_list), i + 1, len(self._index_list)))
                yield img
                print('='*100)
            except:
                print('Find {} reactions, {} of {} can\'t show'.format(
                    len(self._index_list), i + 1, len(self._index_list)))
                yield None
                print('='*100)

    def _save_png(self, show_reagent):
        img_count = 1
        dir_name = '../' + '_'.join([self.chemical, self.feature])
        if not os.path.exists(dir_name):
            os.mkdir(dir_name)
        for img in self._show_reaction(show_reagent, show_inline=False):
            if img:
                filename = os.path.join(dir_name, '_'.join([self.chemical, self.feature, str(img_count)]))
                img.save(filename + '.png')
                img_count += 1


    def counter(self, feature=None, most_common=20):
        if feature:
            self.feature = feature
        chem_list = []
        for doc_feature in self.reaxys_doc[self.feature]:
            if not doc_feature is None:
                chem_list.extend(doc_feature)
        chem_counter = Counter(chem_list)
        pprint.pprint(chem_counter.most_common(most_common))

if __name__ == '__main__':
    reaction = SucReaxys()
    #reaction.search('ibuprofen','product', save_img=True, show_inline=False)
    reaction.counter('product')