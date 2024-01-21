from rdkit import Chem


class Smiles:
    def convert_smiles_to_canonical_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
    
        if mol:            
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
            return canonical_smiles
        
        return None
  
    