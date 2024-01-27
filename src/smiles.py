from rdkit import Chem


class Smiles:
    def convert_smiles_to_canonical_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
    
        if mol:            
            return Chem.MolToSmiles(mol, isomericSmiles=True)
        
        return None
  
    