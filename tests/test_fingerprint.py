import unittest
from src.fingerprint import Fingerprint

class TestFingerprint(unittest.TestCase):
    
    def test_create_atom_pair_returns_atom_pair_fingerprint(self):
        fingerprint = Fingerprint()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'

        atom_pair_fingerprint = fingerprint.create_atom_pair(phenylephrine_cononical_smiles)
        
        self.assertIsNotNone(atom_pair_fingerprint)

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

    def test_create_maccs_returns_none(self):
        fingerprint = Fingerprint()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        maccs_fingerprint = fingerprint.create_maccs(invalid_smiles)

        self.assertIsNone(maccs_fingerprint)

    def test_create_circular_morgan_returns_circular_morgan_fingerprint(self):
        fingerprint = Fingerprint()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'

        morgan_fingerprint = fingerprint.create_morgan(phenylephrine_cononical_smiles)
        
        self.assertIsNotNone(morgan_fingerprint)

if __name__ == '__main__':
    unittest.main()