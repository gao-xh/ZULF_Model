"""Molecule <-> structure specification (zulf_hypothesis.molecule). References: the hand-written chain
specifications and registered motifs of the same compounds (independently written), compared through the
transitions they produce."""
import unittest

import numpy as np

try:
    import rdkit  # noqa: F401
    HAVE_RDKIT = True
except ImportError:                                   # pragma: no cover
    HAVE_RDKIT = False


def lines_components(spec):
    from zulf_studio.session import StudioSession
    return [c["label"] for c in StudioSession(spec, log_file=False).components()]


def lines(spec):
    from zulf_studio.session import StudioSession
    s = StudioSession(spec, log_file=False)
    # relative intensities: the overall scale depends on which isotopologue a specification lists first
    return sorted((round(r["frequency_hz"], 6), round(r["relative"], 6)) for r in s.lines(1e-4))


@unittest.skipUnless(HAVE_RDKIT, "RDKit not installed")
class MoleculeTests(unittest.TestCase):
    def test_ethanol_smiles_gives_the_hand_written_ethanol(self):
        from zulf_hypothesis.molecule import structure_from_molecule
        hand = {"compound": "ethanol", "chain": {"groups": [["C1", "C", 3], ["C2", "C", 2]],
                                                 "bonds": [["C1", "C2"]]}, "one_bond": {"C1": 126.0, "C2": 141.0}}
        spec = structure_from_molecule("CCO", one_bond={"C1": 126.0, "C2": 141.0})
        self.assertEqual(spec["chain"]["groups"], [["C1", "C", 3], ["C2", "C", 2], ["O1", "O", 0]])  # OH dropped
        self.assertEqual(lines(spec), lines(hand))

    def test_propane_symmetry_is_found_and_matches_the_motif(self):
        from zulf_hypothesis.molecule import structure_from_molecule
        spec = structure_from_molecule("CCC", one_bond={"C1": 126.0, "C2": 128.0})
        self.assertEqual(spec["one_bond"], {"C1": 126.0, "C2": 128.0})   # one guess per symmetry orbit
        self.assertTrue(any(p.get("C1") == "C3" and p.get("HC1") == "HC3" for p in spec["chain"]["sym"]))
        # the motif labels the CH2 C1 and the methyls C2, C3; SMILES order gives CH3 C1, CH2 C2, CH3 C3
        motif = {"motif": "CH2(CH3)2", "one_bond": {"C1": 128.0, "C2": 126.0}}
        self.assertEqual(lines(spec), lines(motif))

    def test_symmetry_with_unprotonated_sites(self):
        from zulf_hypothesis.molecule import structure_from_molecule
        glycerol = structure_from_molecule("OCC(O)CO")                   # O-H dropped: the O sites have no H
        self.assertTrue(any(p.get("C1") == "C3" and p.get("HC1") == "HC3" and "HO1" not in p
                            for p in glycerol["chain"]["sym"]))
        self.assertEqual([c for c in lines_components(glycerol)], ["13C@C1 (x2)", "13C@C2"])

    def test_one_bond_guesses_follow_the_neighbours(self):
        from zulf_hypothesis.molecule import structure_from_molecule
        # literature 1J(C,H): ethane 125, methylamine 133, methanol 141, bromomethane 152, isopropylamine CH 133.5
        # / CH3 124.4 (docs/analysis/2026-10-03_isopropylamine_complex-fit.md)
        self.assertEqual(structure_from_molecule("CC")["one_bond"], {"C1": 125.0})
        self.assertEqual(structure_from_molecule("CN")["one_bond"]["C1"], 133.0)
        self.assertEqual(structure_from_molecule("CO")["one_bond"]["C1"], 141.0)
        self.assertEqual(structure_from_molecule("CBr")["one_bond"]["C1"], 152.0)
        ipa = structure_from_molecule("CC(C)N")["one_bond"]
        self.assertEqual((ipa["C1"], ipa["C2"]), (125.0, 133.0))
        self.assertEqual(structure_from_molecule("CC(=O)O")["one_bond"]["C1"], 129.0)

    def test_mol_block_input_unsaturation_warning_and_errors(self):
        from rdkit import Chem
        from zulf_hypothesis.molecule import structure_from_molecule
        block = Chem.MolToMolBlock(Chem.MolFromSmiles("CCN"))
        spec = structure_from_molecule(block)
        self.assertEqual([g[0] for g in spec["chain"]["groups"]], ["C1", "C2", "N1"])
        self.assertEqual(spec["chain"]["groups"][2][2], 2)              # N-H kept (slow exchange possible)
        self.assertIn("notes", structure_from_molecule("CC#N"))         # triple bond: saturated defaults warned
        with self.assertRaises(ValueError):
            structure_from_molecule("C1CC")                              # not a valid SMILES

    def test_molecule_from_structure_for_drawing(self):
        from zulf_hypothesis.molecule import molecule_from_structure, structure_from_molecule
        mol, sites = molecule_from_structure({"motif": "isopropyl"})
        self.assertEqual(mol.GetNumAtoms(), 3)
        self.assertEqual(sorted(a.GetTotalNumHs() for a in mol.GetAtoms()), [1, 3, 3])
        self.assertEqual(set(sites), {"C1", "C2", "C3"})
        mol, sites = molecule_from_structure(structure_from_molecule("CC=O"))   # recorded SMILES: real bonds
        self.assertEqual(str(mol.GetBondWithIdx(1).GetBondType()), "DOUBLE")
        self.assertEqual(sites, {"C1": 0, "C2": 1, "O1": 2})


if __name__ == "__main__":
    unittest.main()
