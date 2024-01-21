from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys


class Fingerprint:

    def create_atom_pair(self, smiles, nBits=2048):
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:             
            return AllChem.GetHashedAtomPairFingerprint(mol, nBits=nBits)
            
        return None

    def create_maccs(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:         
            return MACCSkeys.GenMACCSKeys(mol)

        return None    

    def create_morgan(self, smiles, radius=2, nBits=2048):
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:         
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits)
            
        return None

    def create_feature_based_morgan(self, smiles, radius=2, nBits=2048):
        mol = Chem.MolFromSmiles(smiles)
        
        if mol:         
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nBits, useFeatures=True)
            
        return None

