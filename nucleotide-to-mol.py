from openbabel import pybel

def smiles_to_3d(smiles, output_file, output_format):
    """Convert SMILES to 3D structure file (PDB/XYZ)."""
    try:
        mol = pybel.readstring("smi", smiles)
        mol.addh()  # Add hydrogen atoms
        mol.make3D(forcefield="mmff94", steps=200)
        
        # Ensure output format is correctly passed
        mol.write(output_format, output_file, overwrite=True)
        print(f"Generated: {output_file}")

    except Exception as e:
        print(f"Error processing {output_file}: {str(e)}")

# Dictionary of nucleobases with their SMILES representations
"""
nucleobases = {
    'adenine': 'C1=NC2=NC=NC(=C2N1)N',
    'uracil': 'C1=CNC(=O)NC1=O',
    'thymine': 'CC1=CNC(=O)NC1=O',
    'guanine': 'C1=NC2=C(N1)C(=O)N=C(N2)N',
    'cytosine': 'C1=C(NC(=O)N=C1)N'
}
"""

nucleotides = {
    'adeninosine': 'C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)N',
    'uridine': 'C1=CN(C(=O)NC1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O',
    'guanosine': 'C1=NC2=C(N1[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O)O)N=C(NC2=O)N',
    'cytidine': 'C1=CN(C(=O)N=C1N)[C@H]2[C@@H]([C@@H]([C@H](O2)CO)O)O'
}

# Convert and save each nucleobase in PDB and XYZ format
for name, smiles in nucleotides.items():
    smiles_to_3d(smiles, f"{name}.pdb", "pdb")  # Save in PDB format
    smiles_to_3d(smiles, f"{name}.xyz", "xyz")  # Save in XYZ format