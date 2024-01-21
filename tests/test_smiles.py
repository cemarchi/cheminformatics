import unittest
from src.smiles import Smiles

class TestSmiles(unittest.TestCase):

    def test_convert_smiles_to_canonical_smiles_returns_canonical_smiles(self):
        smiles = Smiles()
        phenylephrine_smiles = 'CNC[C@H](O)C1=CC(O)=CC=C1'

        phenylephrine_canonical_smiles = smiles.convert_smiles_to_canonical_smiles(phenylephrine_smiles)

        self.assertTrue(phenylephrine_canonical_smiles)
        self.assertEqual(phenylephrine_canonical_smiles, 'CNC[C@H](O)c1cccc(O)c1')

    def test_convert_smiles_to_canonical_smiles_returns_none(self):
        smiles = Smiles()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        phenylephrine_canonical_smiles = smiles.convert_smiles_to_canonical_smiles(invalid_smiles)

        self.assertIsNone(phenylephrine_canonical_smiles)

if __name__ == '__main__':
    unittest.main()