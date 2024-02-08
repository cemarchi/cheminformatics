from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys
from rdkit.Avalon import pyAvalonTools


class MolecularFingerprint:

    def create_atom_pair(self, smiles, nBits=2048):
        molecule = Chem.MolFromSmiles(smiles)
        
        if molecule:             
            return AllChem.GetHashedAtomPairFingerprint(molecule, nBits=nBits)
            
        return None

    def create_maccs(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        
        if molecule:         
            return MACCSkeys.GenMACCSKeys(molecule)

        return None    

    def create_extended_connectivity(self, smiles, radius=2, nBits=2048):
        molecule = Chem.MolFromSmiles(smiles)
        
        if molecule:         
            return AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=nBits)
            
        return None

    def create_feature_based_morgan(self, smiles, radius=2, nBits=2048):
        molecule = Chem.MolFromSmiles(smiles)
        
        if molecule:         
            return AllChem.GetMorganFingerprintAsBitVect(molecule, radius, nBits=nBits, useFeatures=True)
            
        return None

    def create_topological_torsion(self, smiles, nBits=2048):
        molecule = Chem.MolFromSmiles(smiles)
        
        if molecule:         
            return AllChem.GetHashedTopologicalTorsionFingerprint(molecule, nBits=nBits)
            
        return None

    def create_avalon_count(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)

        if molecule:
            return pyAvalonTools.GetAvalonCountFP(molecule)

        return None

    def create_avalon_bit(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)

        if molecule:
            return pyAvalonTools.GetAvalonFP(molecule)
        
        return None


