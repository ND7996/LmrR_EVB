#!/usr/bin/env python3
"""
Script to create alanine mutations of specific residues in a PDB file.
Creates separate PDB files with single mutations.
OPLS force field compatible - cleans PDB and handles histidine properly.
"""

import os
import sys
from pathlib import Path

def read_pdb_file(filename):
    """Read PDB file and return clean lines."""
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        # Filter out git conflict markers and other non-PDB lines
        clean_lines = []
        seen_atoms = set()  # Track atom lines to remove duplicates
        
        for line in lines:
            # Skip git conflict markers
            if line.startswith(('<<<<<<', '======', '>>>>>>')):
                continue
            
            # For ATOM lines, check for duplicates
            if line.startswith('ATOM  '):
                # Create a unique key from atom number and coordinates
                try:
                    atom_num = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    res_num = int(line[22:26].strip())
                    x = line[30:38].strip()
                    y = line[38:46].strip()
                    z = line[46:54].strip()
                    atom_key = (atom_num, atom_name, res_num, x, y, z)
                    
                    if atom_key in seen_atoms:
                        continue  # Skip duplicate
                    seen_atoms.add(atom_key)
                    clean_lines.append(line)
                except (ValueError, IndexError):
                    # If parsing fails, keep the line anyway
                    clean_lines.append(line)
            # Keep other standard PDB records
            elif line.startswith(('HETATM', 'TER   ', 'END   ', 
                              'HEADER', 'TITLE ', 'COMPND', 'SOURCE', 
                              'KEYWDS', 'EXPDTA', 'AUTHOR', 'REVDAT',
                              'REMARK', 'SEQRES', 'CRYST1', 'MODEL ', 
                              'ENDMDL', 'CONECT')):
                clean_lines.append(line)
        
        return clean_lines
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found!")
        sys.exit(1)

def get_residue_name(pdb_lines, residue_number):
    """Get the original residue name for a given residue number."""
    for line in pdb_lines:
        if line.startswith('ATOM  '):
            res_num = int(line[22:26].strip())
            if res_num == residue_number:
                return line[17:20].strip()
    return "UNK"

def fix_histidine_atoms(line):
    """
    Don't filter histidine atoms - keep them all as they are in the WT.
    Just pass through the line unchanged.
    """
    return line

def mutate_residue(pdb_lines, target_residue, target_aa):
    """
    Mutate a specific residue to target amino acid.
    For alanine: keeps only backbone atoms (N, CA, C, O) and CB
    For leucine: keeps backbone + CB, CG, CD1, CD2
    For glutamate: keeps backbone + CB, CG, CD, OE1, OE2
    Converts all HIS to HID but keeps all atoms unchanged.
    """
    mutated_lines = []
    
    # Define atoms to keep for each amino acid
    allowed_atoms = {
        'ALA': ['N', 'CA', 'C', 'O', 'CB'],
        'LEU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2'],
        'GLU': ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2']
    }
    
    for line in pdb_lines:
        if line.startswith('ATOM  '):
            res_num = int(line[22:26].strip())
            res_name = line[17:20].strip()
            atom_name = line[12:16].strip()
            
            # Convert HIS to HID for all residues (but don't filter atoms)
            if res_name == 'HIS':
                line = line[:17] + 'HID' + line[20:]
            
            if res_num == target_residue:
                # Keep only allowed atoms for the target amino acid
                if atom_name in allowed_atoms[target_aa]:
                    # Change residue name
                    new_line = line[:17] + target_aa + line[20:]
                    mutated_lines.append(new_line)
                # Skip all other sidechain atoms
            else:
                # Keep all atoms for other residues unchanged
                mutated_lines.append(line)
        elif line.startswith(('TER   ', 'END   ')):
            # Keep TER and END records
            mutated_lines.append(line)
        # Skip all other record types for clean output
    
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
    print(f"Reading and cleaning {input_pdb}...")
    original_lines = read_pdb_file(input_pdb)
    print(f"Cleaned PDB contains {len([l for l in original_lines if l.startswith('ATOM')])} ATOM records")
    
    # Check for HID residues and their atom counts
    hid_residues = {}
    for line in original_lines:
        if line.startswith('ATOM  '):
            res_name = line[17:20].strip()
            res_num = int(line[22:26].strip())
            if res_name in ['HIS', 'HID', 'HIE', 'HIP']:
                if res_num not in hid_residues:
                    hid_residues[res_num] = {'name': res_name, 'atoms': []}
                atom_name = line[12:16].strip()
                hid_residues[res_num]['atoms'].append(atom_name)
    
    if hid_residues:
        print(f"\nHistidine residues found:")
        for res_num, info in sorted(hid_residues.items()):
            print(f"  Residue {res_num} ({info['name']}): {len(info['atoms'])} atoms - {', '.join(info['atoms'])}")
    
    # Create output directory
    output_dir = "mutations"
    Path(output_dir).mkdir(exist_ok=True)
    
    print(f"\nCreating mutations from {input_pdb}")
    print(f"Target mutations: {len(mutation_dict)} total")
    print(f"Output directory: {output_dir}/")
    print("-" * 50)
    
    # Create mutation for each residue
    for residue_num, target_aa in sorted(mutation_dict.items()):
        # Get original residue name
        original_res = get_residue_name(original_lines, residue_num)
        
        # Create mutation
        mutated_lines = mutate_residue(original_lines, residue_num, target_aa)
        
        # Generate output filename
        base_name = os.path.splitext(os.path.basename(input_pdb))[0]
        aa_code = {'ALA': 'A', 'LEU': 'L', 'GLU': 'E'}[target_aa]
        output_filename = f"{output_dir}/{base_name}_{original_res}{residue_num}{aa_code}.pdb"
        
        # Write mutated PDB
        write_pdb_file(mutated_lines, output_filename)
        
        atom_count = len([l for l in mutated_lines if l.startswith('ATOM')])
        print(f"Created: {output_filename} ({original_res}{residue_num} → {target_aa}, {atom_count} atoms)")
    
    print("-" * 50)
    print(f"Successfully created {len(mutation_dict)} mutation files in {output_dir}/")
    print("\nNote: All HIS residues converted to HID (delta-protonated) for OPLS")

def main():
    # Configuration - Define all mutations
    input_pdb_file = "/home/hp/nayanika/github/LmrR_EVB/structures/LMRR_WT2.pdb"
    
    # Mutation dictionary: residue_number -> target_amino_acid
    mutations = {
        # Original alanine mutations
        96: 'ALA', 94: 'ALA', 18: 'ALA', 17: 'ALA', 21: 'ALA',
        87: 'ALA', 99: 'ALA', 7: 'ALA', 88: 'ALA', 92: 'ALA',
        # Additional alanine mutations
        14: 'ALA',   # ASN14A
        9: 'ALA',    # LEU9A
        10: 'ALA',   # ARG10A
        100: 'ALA',  # LYS100A
        103: 'ALA',  # GLU103A
        15: 'ALA',   # ILE15A
        98: 'ALA',   # VAL98A
        95: 'ALA',   # W95A
        # Other mutations
        11: 'LEU',   # A11L
        91: 'GLU'    # A91E
    }
    
    # Check command line arguments
    if len(sys.argv) > 1:
        input_pdb_file = sys.argv[1]
    
    # Create mutation files
    create_mutation_files(input_pdb_file, mutations)

if __name__ == "__main__":
    main()