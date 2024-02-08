import unittest
from parameterized import parameterized, param

from src.molecular_fingerprint_similarity import MolecularFingerprintSimilarity

class TestMolecularFingerprintSimilarity(unittest.TestCase):
    
    def test_calculate_atom_pair_similarity_returns_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_atom_pair_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.75)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])        
    def test_calculate_atom_pair_similarity_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()       
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_atom_pair_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)

    def test_calculate_maccs_similarity_returns_very_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_maccs_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.80)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])  
    def test_calculate_maccs_similarity_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
            
        similarity_ratio = molecular_fingerprint_similarity.calculate_maccs_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)

    def test_calculate_extended_connectivity2_similarity_returns_very_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_extended_connectivity2_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.80)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])  
    def test_calculate_extended_connectivity2_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
            
        similarity_ratio = molecular_fingerprint_similarity.calculate_extended_connectivity2_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)

    def test_calculate_extended_connectivity4_similarity_returns_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_extended_connectivity4_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.70)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])  
    def test_calculate_extended_connectivity4_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
            
        similarity_ratio = molecular_fingerprint_similarity.calculate_extended_connectivity4_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)

    def test_calculate_extended_connectivity6_similarity_returns_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_extended_connectivity6_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.70)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])  
    def test_calculate_extended_connectivity6_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
            
        similarity_ratio = molecular_fingerprint_similarity.calculate_extended_connectivity6_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)

    def test_calculate_feature_based_morgan_similarity_returns_very_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_feature_based_morgan_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.85)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])  
    def test_calculate_feature_based_morgan_similarity_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
            
        similarity_ratio = molecular_fingerprint_similarity.calculate_feature_based_morgan_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)

    def test_calculate_topological_torsion_similarity_returns_highly_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_topological_torsion_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.90)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])  
    def test_calculate_topological_torsion_similarity_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
            
        similarity_ratio = molecular_fingerprint_similarity.calculate_topological_torsion_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)
    
    def test_calculate_avalon_count_similarity_returns_highly_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_avalon_count_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.95)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])  
    def test_calculate_avalon_count_similarity_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
            
        similarity_ratio = molecular_fingerprint_similarity.calculate_avalon_count_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)

    def test_calculate_avalon_bit_similarity_returns_highly_similar(self):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
        phenylephrine_cononical_smiles = 'CNC[C@H](O)c1cccc(O)c1'
        etilefrine_cononical_smiles = 'CCNCC(O)c1cccc(O)c1'
    
        similarity_ratio = molecular_fingerprint_similarity.calculate_avalon_bit_similarity(phenylephrine_cononical_smiles, etilefrine_cononical_smiles)
        
        self.assertIsNotNone(similarity_ratio)
        self.assertGreater(similarity_ratio, 0.90)

    @parameterized.expand([
        ('CNC[C@H](O)C1SxCC(O)!CC=C1', 'CCNCC(O)c1cccc(O)c1'),
        ('CCNCC(O)c1cccc(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1'),
        ('CC#CC(O)c1c@!c!(O)c1', 'CNC[C@H](O)C1SxCC(O)!CC=C1')])  
    def test_calculate_avalon_bit_similarity_returns_none(self, smiles1, smiles2):
        molecular_fingerprint_similarity = MolecularFingerprintSimilarity()
            
        similarity_ratio = molecular_fingerprint_similarity.calculate_avalon_bit_similarity(smiles1, smiles2)
        
        self.assertIsNone(similarity_ratio)

if __name__ == '__main__':
    unittest.main()