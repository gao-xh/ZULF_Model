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
                     exchanged: bool = True, coupled_alpha: float = 0.22, label_alpha: float = 0.32) -> Optional[bytes]:
    """Structure of the molecule `smiles` (heavy-atom indices as in the SMILES) with atom `labelled_atom` as 13C.
    exchanged=True: protons on N/O/S are drawn plain (decoupled by fast exchange); False: grey like the others."""
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
    d = rdMolDraw2D.MolDraw2DCairo(int(size[0]), int(size[1]))
    o = d.drawOptions()
    o.clearBackground = False
    o.bondLineWidth = 3
    o.padding = 0.06
    o.minFontSize = 30
    o.maxFontSize = 40
    o.baseFontSize = 1.0
    o.isotopeLabels = True
    o.setAtomPalette({-1: _rgba(INK)[:3], 1: _rgba(INK)[:3], 7: _rgba(MUTED)[:3], 8: _rgba(MUTED)[:3]})
    d.DrawMolecule(mol, highlightAtoms=atoms, highlightAtomColors=colors, highlightAtomRadii=radii,
                   highlightBonds=[], highlightBondColors={})
    d.FinishDrawing()
    return d.GetDrawingText()
