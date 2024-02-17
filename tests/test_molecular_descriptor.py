import unittest
from mordred import descriptors

from src.molecular_descriptor import MolecularDescriptor

class TestMolecularDescriptor(unittest.TestCase):

    def test_calculate_molecular_descriptors_returns_molecular_descriptors(self):
        molecular_descriptor = MolecularDescriptor()
        phenylephrine_smiles = 'CNC[C@H](O)C1=CC(O)=CC=C1'

        phenylephrine_descriptors = molecular_descriptor.calculate_molecular_descriptors(phenylephrine_smiles)

        self.assertTrue(phenylephrine_descriptors)
        self.assertTrue(len(phenylephrine_descriptors), 1826)

    def test_calculate_molecular_descriptors_returns_none(self):
        molecular_descriptor = MolecularDescriptor()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        phenylephrine_descriptors = molecular_descriptor.calculate_molecular_descriptors(invalid_smiles)

        self.assertIsNone(phenylephrine_descriptors)

    def test_calculate_2d_molecular_descriptors_returns_molecular_descriptors(self):
        molecular_descriptor = MolecularDescriptor()
        phenylephrine_smiles = 'CNC[C@H](O)C1=CC(O)=CC=C1'

        phenylephrine_descriptors = molecular_descriptor.calculate_2d_molecular_descriptors(phenylephrine_smiles)

        self.assertTrue(phenylephrine_descriptors)
        self.assertTrue(len(phenylephrine_descriptors), 1613)

    def test_calculate_2d_molecular_descriptors_returns_none(self):
        molecular_descriptor = MolecularDescriptor()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        phenylephrine_descriptors = molecular_descriptor.calculate_2d_molecular_descriptors(invalid_smiles)

        self.assertIsNone(phenylephrine_descriptors)

    def test_calculate_3d_molecular_descriptors_returns_molecular_descriptors(self):
        molecular_descriptor = MolecularDescriptor()
        phenylephrine_smiles = 'CNC[C@H](O)C1=CC(O)=CC=C1'

        phenylephrine_descriptors = molecular_descriptor.calculate_3d_molecular_descriptors(phenylephrine_smiles)

        self.assertTrue(phenylephrine_descriptors)
        self.assertTrue(len(phenylephrine_descriptors), 1613)

    def test_calculate_3d_molecular_descriptors_returns_none(self):
        molecular_descriptor = MolecularDescriptor()
        invalid_smiles = 'CNC[C@H](O)C1SxCC(O)!CC=C1'

        phenylephrine_descriptors = molecular_descriptor.calculate_3d_molecular_descriptors(invalid_smiles)

        self.assertIsNone(phenylephrine_descriptors)

if __name__ == '__main__':
    unittest.main()