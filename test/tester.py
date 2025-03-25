from openmm.app import PDBFile, Modeller, ForceField
from openmm.app import Simulation, NoCutoff
from openmm import LangevinIntegrator, Platform
import openmm.unit as unit

# ✅ Ensure you've converted MOL2 to PDB first!
pdb = PDBFile("output.pdb")

# ✅ Load correct AMBER force field for RNA
forcefield = ForceField("amber14/rna.xml", "amber14/tip3p.xml")

# ✅ Create modeller & add hydrogens
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.0)

# ✅ Create system without cutoffs (important for RNA base stacking)
system = forcefield.createSystem(modeller.topology, 
                                 nonbondedMethod=NoCutoff,
                                 constraints=None)

# ✅ Langevin integrator for stability
integrator = LangevinIntegrator(300 * unit.kelvin, 
                                1.0 / unit.picosecond, 
                                0.002 * unit.picoseconds)

# ✅ Run on CPU (or change to "CUDA" for GPU)
platform = Platform.getPlatformByName("CPU")

simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)

# ✅ Minimize energy
simulation.min
