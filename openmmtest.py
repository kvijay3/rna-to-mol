from openbabel import pybel
from openff.toolkit import Molecule
from openmm import app, unit
from openmmforcefields.generators import GAFFTemplateGenerator
import openmm as mm

# Corrected SMILES with balanced parentheses
PUBCHEM_SMILES = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)OC[C@@H]4[C@H]([C@H]([C@@H](O4)N5C=NC6=C5N=C(NC6=O)N)O)O)O)N"

def generate_3d_structure(smiles, output_file):
    """Generate 3D structure with explicit stereochemistry"""
    try:
        mol = pybel.readstring("smi", smiles)
        mol.addh()
        mol.make3D(forcefield="mmff94", steps=500)
        mol.localopt()
        mol.write("pdb", output_file, overwrite=True)
        print(f"Generated 3D structure: {output_file}")
    except Exception as e:
        print(f"Structure generation failed: {str(e)}")
        raise

def load_molecule(pdb_file):
    """Load molecule using modern OpenFF methods"""
    try:
        # Use from_file with explicit format specification
        molecule = Molecule.from_file(pdb_file, file_format="pdb", allow_undefined_stereo=True)
        print("Successfully loaded molecule from PDB")
        return molecule
    except Exception as e:
        print(f"Molecule loading failed: {str(e)}")
        raise

def parameterize_with_gaff(molecule, pdb_file):
    """Create AMBER/GAFF parameterized system"""
    try:
        gaff_generator = GAFFTemplateGenerator(molecules=molecule)
        forcefield = app.ForceField("amber14-all.xml")
        forcefield.registerTemplateGenerator(gaff_generator.generator)

        pdb = app.PDBFile(pdb_file)
        
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.NoCutoff,
            constraints=app.HBonds,
            rigidWater=True,
            hydrogenMass=1.5*unit.amu
        )
        print("Created parameterized system")
        return pdb, system
    except Exception as e:
        print(f"Parameterization failed: {str(e)}")
        raise

def optimize_geometry(pdb, system):
    """Hybrid minimization protocol"""
    try:
        integrator = mm.LangevinMiddleIntegrator(
            300*unit.kelvin,
            1.0/unit.picosecond,
            0.002*unit.picoseconds
        )
        simulation = app.Simulation(
            pdb.topology,
            system,
            integrator,
            mm.Platform.getPlatformByName("CPU")
        )
        simulation.context.setPositions(pdb.positions)

        print("\nInitial energy:", simulation.context.getState(getEnergy=True).getPotentialEnergy())
        for tol in [100, 20, 5]:
            simulation.minimizeEnergy(tolerance=tol*unit.kilojoule_per_mole)
            print(f"Energy after {tol} kJ/mol:", 
                  simulation.context.getState(getEnergy=True).getPotentialEnergy())

        return simulation
    except Exception as e:
        print(f"Optimization failed: {str(e)}")
        raise

def main():
    generate_3d_structure(PUBCHEM_SMILES, "AG-mmff94.pdb")
    molecule = load_molecule("AG-mmff94.pdb")
    pdb, system = parameterize_with_gaff(molecule, "AG-mmff94.pdb")
    simulation = optimize_geometry(pdb, system)
    
    state = simulation.context.getState(getPositions=True)
    app.PDBFile.writeFile(
        simulation.topology,
        state.getPositions(),
        open("AG-amber-optimized.pdb", "w")
    )
    print("\nOptimized structure saved to AG-amber-optimized.pdb")

if __name__ == "__main__":
    main()