from openbabel import pybel # type: ignore

# SMILES string
smiles = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)OC[C@@H]4[C@H]([C@H]([C@@H](O4)N5C=NC6=C5N=C(NC6=O)N)O)O)O)N"

# Step 1: Convert SMILES to 3D structure with MMFF94 force field
mol = pybel.readstring("smi", smiles)
mol.make3D(forcefield="gaff", steps=500)  # More steps for better structure
mol.addh()
# Step 2: Save as PDB file
mol.write("pdb", "molecule-gaff.pdb", overwrite=True)

print("Generated molecule.pdb with gaff force field")