import os
import requests

def download_file(url, filename):
    """Download a file from a given URL"""
    try:
        response = requests.get(url, allow_redirects=True)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        with open(filename, 'wb') as file:
            file.write(response.content)
        
        print(f"Successfully downloaded {filename}")
        return True
    except Exception as e:
        print(f"Error downloading {url}: {e}")
        return False

def create_rna_forcefield():
    """Create RNA-specific forcefield XML manually"""
    rna_forcefield_content = '''<?xml version="1.0" encoding="UTF-8"?>
<ForceField>
    <Residues>
        <Residue name="RA">
            <PatchedAtom name="P" type="P" charge="-0.4960"/>
            <PatchedAtom name="O1P" type="OP" charge="-0.4960"/>
            <PatchedAtom name="O2P" type="OP" charge="-0.4960"/>
            <PatchedAtom name="O5'" type="OS" charge="-0.4989"/>
            <PatchedAtom name="C5'" type="CT" charge="0.0066"/>
            <PatchedAtom name="C4'" type="CT" charge="0.1065"/>
            <PatchedAtom name="O4'" type="OS" charge="-0.3548"/>
            <PatchedAtom name="C1'" type="CT" charge="0.0970"/>
            <PatchedAtom name="N9" type="N*" charge="-0.0251"/>
            <PatchedAtom name="C8" type="CB" charge="0.2681"/>
            <PatchedAtom name="N7" type="NC" charge="-0.5453"/>
            <PatchedAtom name="C5" type="CB" charge="0.1245"/>
            <PatchedAtom name="C6" type="CA" charge="0.4960"/>
            <PatchedAtom name="N6" type="N2" charge="-0.9568"/>
            <PatchedAtom name="N1" type="NC" charge="-0.7614"/>
            <PatchedAtom name="C2" type="CA" charge="0.5716"/>
            <PatchedAtom name="N3" type="NC" charge="-0.5838"/>
            <PatchedAtom name="C4" type="CB" charge="0.3053"/>
        </Residue>
        <Residue name="RG">
            <PatchedAtom name="P" type="P" charge="-0.4960"/>
            <PatchedAtom name="O1P" type="OP" charge="-0.4960"/>
            <PatchedAtom name="O2P" type="OP" charge="-0.4960"/>
            <PatchedAtom name="O5'" type="OS" charge="-0.4989"/>
            <PatchedAtom name="C5'" type="CT" charge="0.0066"/>
            <PatchedAtom name="C4'" type="CT" charge="0.1065"/>
            <PatchedAtom name="O4'" type="OS" charge="-0.3548"/>
            <PatchedAtom name="C1'" type="CT" charge="0.0970"/>
            <PatchedAtom name="N9" type="N*" charge="-0.0251"/>
            <PatchedAtom name="C8" type="CB" charge="0.2681"/>
            <PatchedAtom name="N7" type="NC" charge="-0.5453"/>
            <PatchedAtom name="C5" type="CB" charge="0.1245"/>
            <PatchedAtom name="C6" type="CA" charge="0.4960"/>
            <PatchedAtom name="O6" type="O" charge="-0.5226"/>
            <PatchedAtom name="N1" type="NC" charge="-0.7614"/>
            <PatchedAtom name="C2" type="CA" charge="0.5716"/>
            <PatchedAtom name="N2" type="N2" charge="-0.9568"/>
            <PatchedAtom name="N3" type="NC" charge="-0.5838"/>
            <PatchedAtom name="C4" type="CB" charge="0.3053"/>
        </Residue>
    </Residues>
    
    <HarmonicBondForce>
        <!-- Add RNA-specific bond parameters here if needed -->
    </HarmonicBondForce>
    
    <HarmonicAngleForce>
        <!-- Add RNA-specific angle parameters here if needed -->
    </HarmonicAngleForce>
</ForceField>
'''
    
    # Ensure directory exists
    os.makedirs(os.path.expanduser('~/openmm/data'), exist_ok=True)
    
    # Write forcefield XML
    forcefield_path = os.path.expanduser('~/openmm/data/amber14-RNA.xml')
    with open(forcefield_path, 'w') as f:
        f.write(rna_forcefield_content)
    
    print(f"Created RNA forcefield XML at {forcefield_path}")

# Run the forcefield creation
create_rna_forcefield()