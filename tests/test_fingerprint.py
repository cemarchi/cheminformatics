import unittest
from src.fingerprint import Fingerprint

class TestFingerprint(unittest.TestCase):
    
    def test_create_atom_pair_returns_atom_pair_fingerprint(self):
        fingerprint = Fingerprint()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'

        atom_pair_fingerprint = fingerprint.create_atom_pair(phenylephrine_cononical_smiles)
        
        self.assertIsNotNone(atom_pair_fingerprint)
        self.assertEqual(atom_pair_fingerprint.GetLength(), 2048)

    def test_create_atom_pair_returns_none(self):
        fingerprint = Fingerprint()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        atom_pair_fingerprint = fingerprint.create_atom_pair(invalid_smiles)

        self.assertIsNone(atom_pair_fingerprint)

    def test_create_maccs_returns_maccs_fingerprint(self):
        fingerprint = Fingerprint()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'

        maccs_fingerprint = fingerprint.create_maccs(phenylephrine_cononical_smiles)
        
        self.assertIsNotNone(maccs_fingerprint)
        self.assertEqual(maccs_fingerprint.GetNumBits(), 167)

    def test_create_maccs_returns_none(self):
        fingerprint = Fingerprint()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        maccs_fingerprint = fingerprint.create_maccs(invalid_smiles)
        
        self.assertIsNone(maccs_fingerprint)

    def test_create_extended_connectivity_returns_extended_connectivity_fingerprint(self):
        fingerprint = Fingerprint()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'

        extended_connectivity_fingerprint = fingerprint.create_extended_connectivity(phenylephrine_cononical_smiles)
        
        self.assertIsNotNone(extended_connectivity_fingerprint)
        self.assertEqual(extended_connectivity_fingerprint.GetNumBits(), 2048)

    def test_extended_connectivity_returns_none(self):
        fingerprint = Fingerprint()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        extended_connectivity_fingerprint = fingerprint.create_extended_connectivity(invalid_smiles)

        self.assertIsNone(extended_connectivity_fingerprint)

    def test_create_feature_based_morgan_returns_feature_based_morgan_fingerprint(self):
        fingerprint = Fingerprint()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'

        feature_based_morgan_fingerprint = fingerprint.create_feature_based_morgan(phenylephrine_cononical_smiles)
        
        self.assertIsNotNone(feature_based_morgan_fingerprint)
        self.assertEqual(feature_based_morgan_fingerprint.GetNumBits(), 2048)

    def test_create_feature_based_morgan_returns_none(self):
        fingerprint = Fingerprint()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        feature_based_morgan_fingerprint = fingerprint.create_feature_based_morgan(invalid_smiles)

        self.assertIsNone(feature_based_morgan_fingerprint)

    def test_create_topological_torsion_returns_topological_torsion_fingerprint(self):
        fingerprint = Fingerprint()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'

        topological_torsion_fingerprint = fingerprint.create_topological_torsion(phenylephrine_cononical_smiles)
        
        self.assertIsNotNone(topological_torsion_fingerprint)
        self.assertEqual(topological_torsion_fingerprint.GetLength(), 2048)

    def test_create_topological_torsion_returns_none(self):
        fingerprint = Fingerprint()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        topological_torsion_fingerprint = fingerprint.create_topological_torsion(invalid_smiles)

        self.assertIsNone(topological_torsion_morgan_fingerprint)

    def test_create_avalon_count_returns_avalon_count_fingerprint(self):
        fingerprint = Fingerprint()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'

        avalon_count_fingerprint = fingerprint.create_avalon_count(phenylephrine_cononical_smiles)
        
        self.assertIsNotNone(avalon_count_fingerprint)
        self.assertEqual(avalon_count_fingerprint.GetLength(), 512)

    def test_create_avalon_count_returns_none(self):
        fingerprint = Fingerprint()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        avalon_count_fingerprint = fingerprint.create_avalon_count(invalid_smiles)

        self.assertIsNone(avalon_count_fingerprint)

if __name__ == '__main__':
    unittest.main()