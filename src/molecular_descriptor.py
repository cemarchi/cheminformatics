from rdkit import Chem
from mordred import Calculator, descriptors


class MolecularDescriptor:
    def calculate_molecular_descriptors(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        
        if not molecule:         
            return None

        calculator = Calculator(descriptors, ignore_3D=False)

        return calculator(molecule).asdict()

    def calculate_2d_molecular_descriptors(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        
        if not molecule:         
            return None
        
        calculator = Calculator(descriptors, ignore_3D=True)
        descriptor_values = calculator(molecule).asdict()

        return {str(descriptor): descriptor_values[str(descriptor)] for descriptor in calculator.descriptors if not descriptor.require_3D}
        
    def calculate_3d_molecular_descriptors(self, smiles):
        molecule = Chem.MolFromSmiles(smiles)
        
        if not molecule:         
            return None
        
        calculator = Calculator(descriptors, ignore_3D=False)
        descriptor_values = calculator(molecule).asdict()
        
        return {str(descriptor): descriptor_values[str(descriptor)] for descriptor in calculator.descriptors if descriptor.require_3D}
