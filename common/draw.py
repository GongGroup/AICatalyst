from rdkit.Chem import Draw


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
