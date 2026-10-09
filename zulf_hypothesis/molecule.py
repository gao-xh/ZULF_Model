"""Molecules and structure specifications: the interface between a chemical structure (RDKit molecule, SMILES or
mol-file text; later a molecule sketcher) and the structure specifications the simulation and the fits use.

    structure_from_molecule(mol)   molecule -> {"compound", "chain": {"groups", "bonds", "sym"}, "one_bond",
                                   "molecule": {"smiles", "sites"}}
    molecule_from_structure(spec)  structure specification -> (RDKit molecule, {site label: atom index}) for drawing

A chain specification lists every heavy atom as a site (element and the number of non-exchangeable protons on it)
and the bonds between them; zulf_hypothesis.motifs._chain then assigns the couplings by bond distance with
saturated (sp3) defaults, and the 1J of every protonated site comes from "one_bond". structure_from_molecule
fills "one_bond" with starting guesses from the hybridisation (to be refined or fitted) and finds the symmetry
of the molecule (its graph automorphisms) so that equivalent sites share couplings. The "molecule" entry keeps the
SMILES and the site-to-atom map, so a drawing shows the real bond orders; the fits ignore it.

Limits: the chain couplings are the saturated ones, so an aromatic or unsaturated molecule converts with a warning
(use a registered motif such as "benzene ring" for those); protons on O and S are exchangeable and left out unless
keep_exchangeable is set.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

from .structure_spec import fragment_from_spec

# 1J(X,H) starting guesses (Hz) by element and hybridisation; 15N couplings are negative (negative gamma), the
# chain factory applies the sign itself
ONE_BOND_GUESS = {("C", "SP3"): 125.0, ("C", "SP2"): 160.0, ("C", "SP"): 250.0, ("N", "SP3"): 75.0,
                  ("N", "SP2"): 90.0}
EXCHANGEABLE_ON = ("O", "S")
MAX_SYMMETRY_GENERATORS = 48


def _rdkit():
    try:
        from rdkit import Chem
    except ImportError as exc:                    # pragma: no cover - rdkit is in the zulf environment
        raise ImportError("RDKit is needed for molecules (conda environment zulf)") from exc
    return Chem


def to_molecule(molecule):
    """An RDKit molecule from an RDKit molecule, a SMILES string or mol-file (V2000/V3000) text."""
    Chem = _rdkit()
    if isinstance(molecule, Chem.Mol):
        return Chem.Mol(molecule)
    text = str(molecule)
    if "M  END" in text or "V2000" in text or "V3000" in text:
        mol = Chem.MolFromMolBlock(text)          # not stripped: the first (title) line may be empty
    else:
        text = text.strip()
        mol = Chem.MolFromSmiles(text)
    if mol is None:
        raise ValueError(f"not a valid SMILES or mol block: {text[:60]!r}")
    return mol


def structure_from_molecule(molecule, one_bond: Optional[Dict[str, float]] = None, compound: Optional[str] = None,
                            keep_exchangeable: bool = False) -> dict:
    """Structure specification (chain) of a molecule. Sites are the heavy atoms, labelled element + running
    number in atom order (C1, C2, O1, ...); one_bond overrides the 1J guesses by site label. The result is
    checked by building its fragment. Returns the specification with "notes" (warnings) when there are any."""
    Chem = _rdkit()
    mol = Chem.RemoveHs(to_molecule(molecule))
    if mol.GetNumAtoms() == 0:
        raise ValueError("the molecule has no heavy atoms")
    counts: Dict[str, int] = {}
    labels: List[str] = []
    for atom in mol.GetAtoms():
        el = atom.GetSymbol()
        counts[el] = counts.get(el, 0) + 1
        labels.append(f"{el}{counts[el]}")
    groups, guesses, notes = [], {}, []
    for atom, label in zip(mol.GetAtoms(), labels):
        el = atom.GetSymbol()
        n_h = atom.GetTotalNumHs()
        if n_h and el in EXCHANGEABLE_ON and not keep_exchangeable:
            n_h = 0                                # exchangeable: no resolved couplings in solution
        groups.append([label, el, int(n_h)])
        if n_h:
            hyb = str(atom.GetHybridization()).split(".")[-1]
            guess = ONE_BOND_GUESS.get((el, hyb))
            if guess is not None:
                guesses[label] = guess
            elif el in ("C", "N"):
                notes.append(f"no 1J guess for {label} ({el}, {hyb}): set one_bond[{label!r}]")
    bonds = [[labels[b.GetBeginAtomIdx()], labels[b.GetEndAtomIdx()]] for b in mol.GetBonds()]
    if any(b.GetBondType() != Chem.BondType.SINGLE for b in mol.GetBonds()):
        notes.append("unsaturated or aromatic bonds: the chain uses saturated (sp3) 2J/3J defaults; a registered "
                     "motif (e.g. 'benzene ring') models aromatic couplings")
    sym = []
    for match in mol.GetSubstructMatches(mol, uniquify=False, useChirality=False, maxMatches=1000):
        perm = {labels[i]: labels[j] for i, j in enumerate(match) if i != j}
        if perm and perm not in sym:
            sym.append(perm)
        if len(sym) >= MAX_SYMMETRY_GENERATORS:
            break
    order = {l: i for i, l in enumerate(labels)}
    for perm in sym:                               # 1J guesses only on one site of each symmetry orbit
        for a, b in perm.items():
            if a in guesses and b in guesses and order[b] < order[a]:
                guesses.pop(a)
    sym = [_with_protons(p) for p in sym]
    smiles = Chem.MolToSmiles(mol)
    spec = {"compound": compound or smiles,
            "chain": {"groups": groups, "bonds": bonds, "sym": sym},
            "one_bond": {**guesses, **(one_bond or {})},
            "molecule": {"smiles": Chem.MolToSmiles(mol, canonical=False), "sites": {l: i for i, l in
                                                                                      enumerate(labels)}}}
    fragment_from_spec(spec)                       # raises if the specification is inconsistent
    if notes:
        spec["notes"] = notes
    return spec


def molecule_record(spec_or_molecule) -> dict:
    """The "molecule" entry ({"smiles", "sites"}) of a structure specification (its recorded one, else its
    skeleton) or of a molecule (SMILES, mol-file text or RDKit molecule; sites labelled as in
    structure_from_molecule). A spin-system specification carries it for drawing only."""
    Chem = _rdkit()
    if isinstance(spec_or_molecule, dict):
        if spec_or_molecule.get("molecule"):
            return dict(spec_or_molecule["molecule"])
        mol, sites = molecule_from_structure(spec_or_molecule)
        order = sorted(sites, key=sites.get)
        return {"smiles": Chem.MolToSmiles(mol, canonical=False), "sites": {l: i for i, l in enumerate(order)}}
    return dict(structure_from_molecule(spec_or_molecule)["molecule"])


def _with_protons(perm: Dict[str, str]) -> Dict[str, str]:
    """A site permutation also maps the proton groups on the sites (H<site> -> H<image>)."""
    return {**perm, **{f"H{a}": f"H{b}" for a, b in perm.items()}}


def molecule_from_structure(spec: dict) -> Tuple[object, Dict[str, int]]:
    """RDKit molecule and {site label: atom index} of a structure specification, for drawing. With a "molecule"
    entry (structure_from_molecule) the recorded SMILES (real bond orders); otherwise the heavy-atom skeleton of
    the fragment (motif or chain) with single bonds and the proton counts of its groups."""
    Chem = _rdkit()
    rec = spec.get("molecule")
    if rec and rec.get("smiles"):
        mol = Chem.MolFromSmiles(rec["smiles"])
        if mol is not None:
            return mol, {k: int(v) for k, v in rec.get("sites", {}).items()}
    if "spin_system" in spec:
        raise ValueError("a typed-in spin system has no molecule (draw it as a spin network)")
    frag = fragment_from_spec(spec)
    rw = Chem.RWMol()
    index = {}
    for site in frag.sites:
        atom = Chem.Atom(site.element)
        atom.SetNoImplicit(True)
        index[site.label] = rw.AddAtom(atom)
    for group in frag.protons:
        atom = rw.GetAtomWithIdx(index[group.site])
        atom.SetNumExplicitHs(atom.GetNumExplicitHs() + group.size)
    for a, b in frag.bonds:
        if rw.GetBondBetweenAtoms(index[a], index[b]) is None:
            rw.AddBond(index[a], index[b], Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    return mol, index
