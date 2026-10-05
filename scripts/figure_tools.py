"""Small helpers for publication figures (scripts/paper_figure.py).

- `isotopologue_png`: a structure drawing of one 13C isotopologue (RDKit, optional): the labelled carbon and its
  protons highlighted in the isotopologue's colour, the other carbon-bound protons (the coupled ones) in grey,
  heteroatom-bound protons left plain (decoupled by fast exchange). Returns PNG bytes with a transparent
  background, or None when RDKit is not installed.
- `PALETTE`: isotopologue colours in fixed order (checked with the dataviz palette validator: lightness band,
  chroma, colour-vision separation and contrast pass on a white surface).
"""
from typing import Optional, Sequence

PALETTE = ("#c0469e", "#2f8f5b", "#2a78d6", "#eb6834")
INK = "#1d1d1f"
MUTED = "#6b6b70"
GRID = "#dcdcd8"
SIMULATION = "#7b2f8f"


def _rgba(hex_color: str, alpha: float = 1.0):
    h = hex_color.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)) + (alpha,)


def isotopologue_png(smiles: str, labelled_atom: int, color: str, size: Sequence[int] = (520, 440),
                     exchanged: bool = True, coupled_alpha: float = 0.22, label_alpha: float = 0.32,
                     bond_note: Optional[str] = None, return_coords: bool = False, mark_bond: bool = False):
    """Structure of the molecule `smiles` (heavy-atom indices as in the SMILES) with atom `labelled_atom` as 13C.
    exchanged=True: protons on N/O/S are drawn plain (decoupled by fast exchange); False: grey like the others.
    bond_note: text written by RDKit on one bond from the 13C to its protons (small); return_coords=True returns
    (png, {"c": (x, y), "h": (x, y)}) with the pixel positions of the 13C and of one of its protons, for a label
    drawn by the caller instead; mark_bond=True draws that 13C-H bond thick in the colour."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError:
        return None
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    label = mol.GetAtomWithIdx(int(labelled_atom))
    label.SetIsotope(13)
    group = {label.GetIdx()} | {n.GetIdx() for n in label.GetNeighbors() if n.GetAtomicNum() == 1}
    atoms, colors, radii = [], {}, {}
    for i in sorted(group):
        atoms.append(i)
        colors[i] = _rgba(color, label_alpha)
        radii[i] = 0.42
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 1 or a.GetIdx() in group:
            continue
        on = a.GetNeighbors()[0].GetAtomicNum()
        if on == 6 or not exchanged:
            atoms.append(a.GetIdx())
            colors[a.GetIdx()] = (0.55, 0.55, 0.58, coupled_alpha)
            radii[a.GetIdx()] = 0.36
    if bond_note:
        hs = [n.GetIdx() for n in label.GetNeighbors() if n.GetAtomicNum() == 1]
        if hs:
            mol.GetBondBetweenAtoms(label.GetIdx(), hs[0]).SetProp("bondNote", bond_note)
    bonds, bond_colors = [], {}
    if mark_bond:
        hs = [n.GetIdx() for n in label.GetNeighbors() if n.GetAtomicNum() == 1]
        if hs:
            # the proton lowest in the drawing (largest y after the flip): the label goes below the molecule
            conf = mol.GetConformer()
            h = min(hs, key=lambda i: conf.GetAtomPosition(i).y)
            b = mol.GetBondBetweenAtoms(label.GetIdx(), h).GetIdx()
            bonds, bond_colors = [b], {b: _rgba(color, 0.9)}
            mol.SetProp("_marked_h", str(h))
    d = rdMolDraw2D.MolDraw2DCairo(int(size[0]), int(size[1]))
    o = d.drawOptions()
    o.clearBackground = False
    o.bondLineWidth = 3
    o.padding = 0.06
    o.minFontSize = 30
    o.maxFontSize = 40
    o.baseFontSize = 1.0
    o.isotopeLabels = True
    o.annotationFontScale = 0.62
    o.setAtomPalette({-1: _rgba(INK)[:3], 1: _rgba(INK)[:3], 7: _rgba(MUTED)[:3], 8: _rgba(MUTED)[:3]})
    o.highlightBondWidthMultiplier = 4
    d.DrawMolecule(mol, highlightAtoms=atoms, highlightAtomColors=colors, highlightAtomRadii=radii,
                   highlightBonds=bonds, highlightBondColors=bond_colors)
    d.FinishDrawing()
    if not return_coords:
        return d.GetDrawingText()
    hs = [n.GetIdx() for n in label.GetNeighbors() if n.GetAtomicNum() == 1]
    pos = lambda i: tuple(float(v) for v in (d.GetDrawCoords(i).x, d.GetDrawCoords(i).y))
    marked = int(mol.GetProp("_marked_h")) if mol.HasProp("_marked_h") else (hs[0] if hs else label.GetIdx())
    return d.GetDrawingText(), {"c": pos(label.GetIdx()), "h": pos(marked), "hs": [pos(i) for i in hs],
                                "size": (int(size[0]), int(size[1]))}
