from Bio.PDB import PDBParser
import numpy as np
import os

# ------------ INPUTS ----------------
input_pdb = r"D:\PhD_Thesis\LmrR_EVB\structures\6i8n.pdb"

# Residues to measure distance from MPO N1
target_residues = [
    ('E', 7), ('A', 11), ('L', 18), ('N', 19), ('K', 22),
    ('N', 88), ('M', 89), ('A', 92), ('F', 93), ('S', 95),
    ('S', 97), ('D', 100)
]

# Chain where your protein residues are located
chain_id = "A"

# MPO atom details
ligand_resname = "MPO"  # Ligand residue name in PDB
ligand_atom_name = "N1"  # Atom in MPO to measure from
# -----------------------------------

# Load PDB structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", input_pdb)

# Function to find MPO N1 atom coordinates
def get_mpo_n1_coord(structure):
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ligand_resname:
                    for atom in residue:
                        if atom.get_name().strip() == ligand_atom_name:
                            return atom.coord
    return None

# Get MPO N1 coordinate
mpo_coord = get_mpo_n1_coord(structure)
if mpo_coord is None:
    raise ValueError(f"Could not find atom {ligand_atom_name} in residue {ligand_resname} in the PDB!")

# Function to calculate Euclidean distance
def distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

# Collect distances
distances = []

for model in structure:
    for chain in model:
        if chain.id == chain_id:
            for residue in chain:
                resname = residue.get_resname().strip()
                resnum = residue.get_id()[1]
                for target_name, target_num in target_residues:
                    if resnum == target_num:
                        # Use CA atom for the residue
                        if "CA" in residue:
                            ca_coord = residue["CA"].coord
                            dist = distance(mpo_coord, ca_coord)
                            distances.append((f"{resname}{resnum}", round(dist, 3)))

# Save results to a text file
output_file = os.path.join(os.path.dirname(input_pdb), "residue_distances.txt")
with open(output_file, "w") as f:
    f.write("Residue\tDistance_from_MPO_N1(Å)\n")
    for res, dist in distances:
        f.write(f"{res}\t{dist}\n")

# Print to terminal
print("Distances from MPO N1:")
for res, dist in distances:
    print(f"{res} -> {dist} Å")

print(f"\nSaved to: {output_file}")
