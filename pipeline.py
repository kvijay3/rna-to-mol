from openbabel import pybel
from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np

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

def prepare_openmm_structure(input_pdb, output_pdb):
    """Prepare PDB for OpenMM simulation"""
    # Load RNA force field
    forcefield = ForceField('amber14-RNA.xml')
    
    # Load PDB
    pdb = PDBFile(input_pdb)
    
    # Rename residues explicitly
    topology = pdb.topology
    for residue in topology.residues():
        # Explicitly set residue names
        if residue.name == 'A':
            residue.name = 'DA5'  # 5' adenosine
        elif residue.name == 'G':
            residue.name = 'DG3'  # 3' guanosine
    
    # Create modeller
    modeller = Modeller(topology, pdb.positions)
    
    # Create system
    system = forcefield.createSystem(
        modeller.topology, 
        nonbondedMethod=NoCutoff,
        constraints=HBonds
    )
    
    # Define integrator
    integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
    
    # Create simulation
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    
    # Minimize energy
    simulation.minimizeEnergy()
    
    # Save processed PDB
    with open(output_pdb, "w") as pdbfile:
        PDBFile.writeFile(
            simulation.topology, 
            simulation.context.getState(getPositions=True).getPositions(), 
            pdbfile
        )
    
    print(f"Processed PDB saved as {output_pdb}")

# SMILES for AG
AG = "C1=NC(=C2C(=N1)N(C=N2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO)OP(=O)(O)OC[C@@H]4[C@H]([C@H]([C@@H](O4)N5C=NC6=C5N=C(NC6=O)N)O)O)O)N"

# Main workflow
generate_3d_structure(AG, 'initial_structure.pdb')
prepare_openmm_structure('initial_structure.pdb', 'final_structure.pdb')