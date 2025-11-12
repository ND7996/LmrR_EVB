#!/usr/bin/env python3
"""
Script to create alanine mutations of specific residues in a PDB file.
Creates separate PDB files with single mutations.
"""

import os
import sys
from pathlib import Path

def read_pdb_file(filename):
    """Read PDB file and return lines."""
    try:
        with open(filename, 'r') as f:
            return f.readlines()
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found!")
        sys.exit(1)

def get_residue_name(pdb_lines, residue_number):
    """Get the original residue name for a given residue number."""
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM')):
            res_num = int(line[22:26].strip())
            if res_num == residue_number:
                return line[17:20].strip()
    return "UNK"

def mutate_residue(pdb_lines, target_residue, target_aa):
    """
    Mutate a specific residue to target amino acid.
    For alanine: keeps only backbone atoms (N, CA, C, O) and CB
    For leucine: keeps backbone + CB, CG, CD1, CD2
    For glutamate: keeps backbone + CB, CG, CD, OE1, OE2
    For tryptophan: keeps backbone + CB, CG, CD1, CD2, NE1, CE2, CE3, CZ2, CZ3, CH2
    """
    mutated_lines = []
    
    # Define atoms to keep for each amino acid
    allowed_atoms = {
        'ALA': ['N', 'CA', 'C', 'O', 'CB'],
        'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
        'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2'],
        'TRP': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
    }
    
    for line in pdb_lines:
        if line.startswith(('ATOM', 'HETATM')):
            res_num = int(line[22:26].strip())
            atom_name = line[12:16].strip()
            
            if res_num == target_residue:
                # Keep only allowed atoms for the target amino acid
                if atom_name in allowed_atoms[target_aa]:
                    # Change residue name
                    new_line = line[:17] + target_aa + line[20:]
                    mutated_lines.append(new_line)
                # Skip all other sidechain atoms
            else:
                # Keep all atoms for other residues
                mutated_lines.append(line)
        else:
            # Keep non-ATOM lines (headers, etc.)
            mutated_lines.append(line)
    
    return mutated_lines

def write_pdb_file(lines, filename):
    """Write PDB lines to file."""
    with open(filename, 'w') as f:
        f.writelines(lines)

def create_mutation_files(input_pdb, mutation_dict):
    """Create separate PDB files for each mutation."""
    
    # Check if input file exists
    if not os.path.exists(input_pdb):
        print(f"Error: Input PDB file '{input_pdb}' not found!")
        return
    
    # Read the original PDB file
    original_lines = read_pdb_file(input_pdb)
    
    # Create output directory
    output_dir = "mutations"
    Path(output_dir).mkdir(exist_ok=True)
    
    print(f"Creating mutations from {input_pdb}")
    print(f"Target mutations: {mutation_dict}")
    print(f"Output directory: {output_dir}/")
    print("-" * 50)
    
    # Create mutation for each residue
    for residue_num, target_aa in mutation_dict.items():
        # Get original residue name
        original_res = get_residue_name(original_lines, residue_num)
        
        # Create mutation
        mutated_lines = mutate_residue(original_lines, residue_num, target_aa)
        
        # Generate output filename
        base_name = os.path.splitext(os.path.basename(input_pdb))[0]
        aa_code = {'ALA': 'A', 'LEU': 'L', 'GLU': 'E', 'TRP': 'W'}[target_aa]
        output_filename = f"{output_dir}/{base_name}_{original_res}{residue_num}{aa_code}.pdb"
        
        # Write mutated PDB
        write_pdb_file(mutated_lines, output_filename)
        
        print(f"Created: {output_filename} ({original_res}{residue_num} → {target_aa})")
    
    print("-" * 50)
    print(f"Successfully created {len(mutation_dict)} mutation files in {output_dir}/")

def main():
    # Configuration - Define all mutations
    input_pdb_file = "WT.pdb"  # Change this to your input file name
    
    # Mutation dictionary: residue_number -> target_amino_acid
    mutations = {
        # Alanine mutations
        96: 'ALA', 94: 'ALA', 18: 'ALA', 17: 'ALA', 21: 'ALA',
        87: 'ALA', 99: 'ALA', 7: 'ALA', 88: 'ALA', 92: 'ALA', 95: 'ALA',
        11: 'LEU',  # A11L
        91: 'GLU'   # A91E
    }
    
    # Check command line arguments
    if len(sys.argv) > 1:
        input_pdb_file = sys.argv[1]
    
    # Create mutation files
    create_mutation_files(input_pdb_file, mutations)

if __name__ == "__main__":
    main()
