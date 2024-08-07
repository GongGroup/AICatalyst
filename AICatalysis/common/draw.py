from rdkit import Chem
from rdkit.Chem import Draw, AllChem

from AICatalysis.common.constant import ChemInfo


def draw_png(mol, width=800, height=400, file="figure.png"):
    d = Draw.MolDraw2DCairo(width, height)
    Draw.PrepareMolForDrawing(mol)
    d.DrawMolecule(mol)
    d.FinishDrawing()
    d.WriteDrawingText(file)


def draw_svg(mol, width=800, height=400, file="figure.svg"):
    d = Draw.MolDraw2DSVG(width, height)
    Draw.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    with open(file, 'w') as f:
        f.write(svg)


def draw_reaction(rea: list[str], reagent=None, pro=None, template=None, file="reaction.png"):
    def kekulize(smarts: str):
        odd = True
        new_smarts = []
        for char in smarts:
            if char == ":":
                if odd:
                    char = "="
                    odd = False
                else:
                    char = "-"
                    odd = True
            new_smarts.append(char)
        return "".join(new_smarts)

    reactants = [Chem.MolFromSmiles(ChemInfo[reactant]['smiles']) for reactant in rea]
    reactant_smarts = [Chem.MolToSmarts(reactant) for reactant in reactants]
    if pro is not None:
        products = [Chem.MolFromSmiles(ChemInfo[product]['smiles']) for product in pro]
        product_smarts = [Chem.MolToSmarts(product) for product in products]
    elif template is not None:
        reaction_template = AllChem.ReactionFromSmarts(template)

        products = reaction_template.RunReactants(reactants)
        product_smarts = [Chem.MolToSmarts(product) for product in products[0]]
    else:
        raise TypeError(f"<pro> and <template> should have one")
    if reagent is not None:
        reagents = [Chem.MolFromSmiles(ChemInfo[rea]['smiles']) for rea in reagent]
        reagent_smarts = [Chem.MolToSmarts(reagent) for reagent in reagents]
    else:
        reagent_smarts = []
    reaction_smarts = f"{'.'.join(reactant_smarts)} >{'.'.join(reagent_smarts)}> {'.'.join(product_smarts)}"
    reaction_smarts = kekulize(reaction_smarts)
    reaction = AllChem.ReactionFromSmarts(reaction_smarts)
    image = Draw.ReactionToImage(reaction)
    image.save(file)


if __name__ == '__main__':
    pass
