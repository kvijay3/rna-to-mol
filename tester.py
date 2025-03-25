from openbabel import pybel
from simtk.openmm import app
import openmm as mm
import simtk.unit as unit

# Step 1: Convert SMILES to PDB using Open Babel
def smiles_to_3d(smiles, output_file):
    """Convert SMILES to a 3D PDB file with MMFF94."""
    try:
        mol = pybel.readstring("smi", smiles)
        mol.make3D(forcefield="mmff94", steps=500)  # More steps for better structure
        mol.write("pdb", output_file, overwrite=True)
        print(f"Generated: {output_file}")
    except Exception as e:
        print(f"Error processing {output_file}: {str(e)}")

# Step 2: Fix PDB atom names to match OpenMM expectations
def fix_pdb_atom_names(input_pdb, output_pdb):
    with open(input_pdb, "r") as infile, open(output_pdb, "w") as outfile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                modified_line = line[:12] + line[12:16].replace("*", "'") + line[16:]
                outfile.write(modified_line)
            else:
                outfile.write(line)

# Define RNA fragment
AG = "C1=NC(=C2C(=N1)N(C=N2)C3C(C(C(O3)CO)OP(=O)(O)OCC4C(C(C(O4)N5C=NC6=C5N=C(NC6=O)N)O)O)O)N"
smiles_to_3d(AG, "AG-mmff94.pdb")

# Fix atom names
fix_pdb_atom_names("AG-mmff94.pdb", "AG-mmff94-e.pdb")

# Step 3: Load the PDB and apply AMBER OL3 force field in OpenMM
pdb = app.PDBFile("AG-mmff94-e.pdb")

# ✅ Use RNA OL3 force field for accurate RNA parameterization
forcefield = app.ForceField("amber14/RNA.OL3.xml", "amber14/tip3p.xml")

# Ensure hydrogen atoms are added correctly
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield)

# Create the system with no cutoff (important for RNA stacking)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds)

# Set up a Langevin integrator
integrator = mm.LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picosecond, 
                                   0.002 * unit.picoseconds)

# Set up the simulation
platform = mm.Platform.getPlatformByName("CPU")
simulation = app.Simulation(modeller.topology, system, integrator, platform)

# Initialize positions
simulation.context.setPositions(modeller.positions)

# Minimize energy (important for correcting structure)
simulation.minimizeEnergy()

# Save the optimized structure
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open("AG-amber-ol3.pdb", "w"))

print("Optimized structure saved as AG-amber-ol3.pdb")