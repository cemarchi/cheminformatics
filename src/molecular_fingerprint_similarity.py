from rdkit.Chem import DataStructs

from src.molecular_fingerprint import MolecularFingerprint

class MolecularFingerprintSimilarity:
    def __init__(self):
        self.__molecular_fingerprint = MolecularFingerprint()

    def calculate_atom_pair_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_atom_pair(smiles1)
        fingerprint2 = self.__molecular_fingerprint.create_atom_pair(smiles2)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.DiceSimilarity(fingerprint1, fingerprint2)

    def calculate_maccs_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_maccs(smiles1)
        fingerprint2 = self.__molecular_fingerprint.create_maccs(smiles2)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.TanimotoSimilarity(fingerprint1, fingerprint2)

    def calculate_extended_connectivity2_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_extended_connectivity(smiles1, radius=2)
        fingerprint2 = self.__molecular_fingerprint.create_extended_connectivity(smiles2, radius=2)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.DiceSimilarity(fingerprint1, fingerprint2)

    def calculate_extended_connectivity4_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_extended_connectivity(smiles1, radius=4)
        fingerprint2 = self.__molecular_fingerprint.create_extended_connectivity(smiles2, radius=4)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.DiceSimilarity(fingerprint1, fingerprint2)

    def calculate_extended_connectivity6_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_extended_connectivity(smiles1, radius=6)
        fingerprint2 = self.__molecular_fingerprint.create_extended_connectivity(smiles2, radius=6)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.DiceSimilarity(fingerprint1, fingerprint2)


    def calculate_feature_based_morgan_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_feature_based_morgan(smiles1)
        fingerprint2 = self.__molecular_fingerprint.create_feature_based_morgan(smiles2)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.DiceSimilarity(fingerprint1, fingerprint2)

    def calculate_topological_torsion_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_topological_torsion(smiles1)
        fingerprint2 = self.__molecular_fingerprint.create_topological_torsion(smiles2)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.DiceSimilarity(fingerprint1, fingerprint2)

    def calculate_avalon_count_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_avalon_count(smiles1)
        fingerprint2 = self.__molecular_fingerprint.create_avalon_count(smiles2)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.DiceSimilarity(fingerprint1, fingerprint2)

    def calculate_avalon_bit_similarity(self, smiles1, smiles2):
        fingerprint1 = self.__molecular_fingerprint.create_avalon_bit(smiles1)
        fingerprint2 = self.__molecular_fingerprint.create_avalon_bit(smiles2)

        if not fingerprint1 or not fingerprint2:
            return None

        return DataStructs.DiceSimilarity(fingerprint1, fingerprint2)