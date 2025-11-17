from Bio.PDB import PDBParser, PDBIO, Selection
import os

# ----------- INPUTS -----------
input_pdb = r"D:\PhD_Thesis\LmrR_EVB\structures\6i8n.pdb"
output_dir = r"D:\PhD_Thesis\LmrR_EVB\prep_structures\alanine_mutants"

# Residues you want to mutate
residues_to_mutate = [
    ('E', 7), ('A', 11), ('L', 18), ('N', 19), ('K', 22),
    ('N', 88), ('M', 89), ('A', 92), ('F', 93), ('S', 95),
    ('S', 97), ('D', 100)
]
# --------------------------------

# Make sure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load PDB structure
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", input_pdb)

# Function to mutate to ALA by simply changing the residue name
def mutate_residue_to_alanine(structure, chain_id, residue_number):
    for model in structure:
        for chain in model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.get_id()[1] == residue_number:
                        residue.resname = "ALA"
                        return structure
    return None

# Generate mutants
for original_res, position in residues_to_mutate:
    # Reload structure fresh each time
    structure = parser.get_structure("protein", input_pdb)

    # Assuming chain is 'A'
    chain_id = "A"

    # Mutate residue
    mutated_structure = mutate_residue_to_alanine(structure, chain_id, position)

    if mutated_structure:
        output_path = os.path.join(output_dir, f"{original_res}{position}A_mutant.pdb")
        io = PDBIO()
        io.set_structure(mutated_structure)
        io.save(output_path)
        print(f"Saved mutant: {output_path}")
    else:
        print(f"Residue {original_res}{position} not found in chain {chain_id}!")
