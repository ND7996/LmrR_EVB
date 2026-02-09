#!/usr/bin/env python3
"""
Convert Schrödinger ffld_server output to .prm and .lib files
Usage: python ffld_to_prm.py <ffld_log_file> <pdb_file>
"""

import sys
import os
import re

def parse_pdb(pdb_file):
    """Parse PDB file to get atom information"""
    atoms = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('HETATM') or line.startswith('ATOM'):
                # Standard PDB format
                try:
                    atom_num = int(line[6:11].strip())
                    atom_name = line[12:16].strip()
                    res_name = line[17:20].strip()
                    res_num = int(line[22:26].strip())
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    element = line[76:78].strip()
                    
                    atoms.append({
                        'num': atom_num,
                        'name': atom_name,
                        'res_name': res_name,
                        'res_num': res_num,
                        'x': x, 'y': y, 'z': z,
                        'element': element
                    })
                except (ValueError, IndexError):
                    continue
    return atoms

def parse_ffld_log(ffld_file):
    """Parse ffld_server output log file"""
    with open(ffld_file, 'r') as f:
        content = f.read()
    
    print("=== PARSING FFLD LOG ===")
    
    # Extract molecule name
    mol_name_match = re.search(r'molecule:\s*(\S+)', content, re.IGNORECASE)
    mol_name = mol_name_match.group(1) if mol_name_match else 'UNK'
    print(f"Molecule name: {mol_name}")
    
    # Extract total charge
    charge_match = re.search(r'Total charge.*?([-+]?\d*\.?\d+)', content, re.IGNORECASE)
    total_charge = float(charge_match.group(1)) if charge_match else 0.0
    print(f"Total charge: {total_charge}")
    
    # Extract atom types and charges from the detailed table
    atom_charges = {}
    atom_types = {}
    atom_names = {}
    
    # Try multiple patterns to find the atom table
    patterns = [
        r'atom\s+type\s+vdw\s+symbol\s+charge.*?\n(.*?)\n-{70,}',
        r'atom.*?type.*?vdw.*?symbol.*?charge.*?\n(.*?)(?:\n\n|\n\s*[A-Z]|\Z)',
        r'Atom.*?Type.*?Vdw.*?Symbol.*?Charge.*?\n(.*?)(?:\n\n|\n\s*[A-Z]|\Z)',
        r'atom.*?opls.*?charge.*?\n(.*?)(?:\n\n|\n\s*[A-Z]|\Z)'
    ]
    
    atom_table_found = False
    for pattern in patterns:
        atom_table_match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
        if atom_table_match:
            print("Found atom parameter table")
            atom_lines = atom_table_match.group(1).strip().split('\n')
            for line in atom_lines:
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        atom_name = parts[0]
                        opls_type = int(parts[1])
                        charge = float(parts[4])
                        
                        # Map atom names to numbers (O1 -> 1, C2 -> 2, etc.)
                        atom_num_match = re.search(r'(\d+)$', atom_name)
                        if atom_num_match:
                            atom_num = int(atom_num_match.group(1))
                            atom_charges[atom_num] = charge
                            atom_types[atom_num] = opls_type
                            atom_names[atom_num] = atom_name
                            print(f"Table atom: {atom_num}, {atom_name}, opls_{opls_type}, charge: {charge}")
                    except (ValueError, IndexError):
                        continue
            atom_table_found = True
            break
    
    if not atom_table_found:
        print("Warning: Could not find atom parameter table")
        print("Trying to extract charges from energy components...")
        
        # Try to find charges in the energy breakdown
        charge_pattern = r'(\w+\d+)\s+[-\d\.]+\s+[-\d\.]+\s+([-\d\.]+)'
        charge_matches = re.findall(charge_pattern, content)
        for atom_name, charge_str in charge_matches:
            try:
                charge = float(charge_str)
                atom_num_match = re.search(r'(\d+)$', atom_name)
                if atom_num_match:
                    atom_num = int(atom_num_match.group(1))
                    atom_charges[atom_num] = charge
                    atom_names[atom_num] = atom_name
                    print(f"Charge atom: {atom_num}, {atom_name}, charge: {charge}")
            except ValueError:
                continue
    
    # Parse bond information from various sections
    bonds = []
    bond_patterns = [
        r'Stretch.*?k.*?r0.*?\n(.*?)(?=\n\n|\n [A-Z]|\Z)',
        r'Bond.*?Force.*?Length.*?\n(.*?)(?=\n\n|\n [A-Z]|\Z)',
        r'(\w+\d+)\s+(\w+\d+)\s+([\d\.]+)\s+([\d\.]+)'
    ]
    
    bond_found = False
    for pattern in bond_patterns:
        bond_section_match = re.search(pattern, content, re.IGNORECASE | re.DOTALL)
        if bond_section_match:
            print("Found bond section")
            if pattern == bond_patterns[2]:  # Direct bond pattern
                bond_matches = re.findall(pattern, content)
                for match in bond_matches:
                    try:
                        atom1_name, atom2_name, force_const_str, length_str = match
                        force_const = float(force_const_str)
                        length = float(length_str)
                        
                        # Get atom numbers from names
                        atom1_num = int(re.search(r'(\d+)$', atom1_name).group(1))
                        atom2_num = int(re.search(r'(\d+)$', atom2_name).group(1))
                        
                        bonds.append({
                            'atom1': atom_names.get(atom1_num, atom1_name),
                            'atom2': atom_names.get(atom2_num, atom2_name),
                            'atom1_num': atom1_num,
                            'atom2_num': atom2_num,
                            'length': length,
                            'force_const': force_const
                        })
                        print(f"Bond: {atom1_name}-{atom2_name}, length: {length}, force_const: {force_const}")
                    except (ValueError, AttributeError):
                        continue
            else:
                bond_lines = bond_section_match.group(1).strip().split('\n')
                for line in bond_lines:
                    parts = line.split()
                    if len(parts) >= 5:
                        try:
                            atom1_name = parts[0]
                            atom2_name = parts[1]
                            force_const = float(parts[2])
                            length = float(parts[3])
                            
                            # Get atom numbers from names
                            atom1_num = int(re.search(r'(\d+)$', atom1_name).group(1))
                            atom2_num = int(re.search(r'(\d+)$', atom2_name).group(1))
                            
                            bonds.append({
                                'atom1': atom_names.get(atom1_num, atom1_name),
                                'atom2': atom_names.get(atom2_num, atom2_name),
                                'atom1_num': atom1_num,
                                'atom2_num': atom2_num,
                                'length': length,
                                'force_const': force_const
                            })
                            print(f"Bond: {atom1_name}-{atom2_name}, length: {length}, force_const: {force_const}")
                        except (ValueError, IndexError, AttributeError):
                            continue
            bond_found = True
            break
    
    if not bond_found:
        print("Warning: Could not find bond section")
    
    return mol_name, total_charge, atom_charges, atom_types, bonds, atom_names

def write_prm_file(mol_name, atom_types, bonds, output_file):
    """Write .prm file with OPLS parameters"""
    with open(output_file, 'w') as f:
        f.write(f"! OPLS parameters for {mol_name}\n")
        f.write(f"! Generated from ffld_server output\n\n")
        
        # Write bond parameters
        f.write("BONDS\n")
        if bonds:
            for bond in bonds:
                type1 = atom_types.get(bond['atom1_num'], 0)
                type2 = atom_types.get(bond['atom2_num'], 0)
                kb = bond['force_const'] * 100.0  # Convert to appropriate units
                r0 = bond['length']
                f.write(f"{type1:>4d} {type2:>4d} {kb:8.2f} {r0:8.4f}  ! {bond['atom1']}-{bond['atom2']}\n")
        else:
            f.write("! No bond parameters found\n")
        
        f.write("\nANGLES\n")
        f.write("! Angle parameters need manual addition\n")
        
        f.write("\nDIHEDRALS\n")
        f.write("! Dihedral parameters need manual addition\n")
        
        f.write("\nIMPROPER\n")
        f.write("! Improper parameters need manual addition\n")
        
        f.write("\nNONBONDED\n")
        unique_types = set(atom_types.values())
        if unique_types:
            for opls_type in sorted(unique_types):
                # Default OPLS parameters
                epsilon = 0.066  # kcal/mol
                sigma = 3.5      # Å
                f.write(f"{opls_type:>4d}  0.0000  {epsilon:8.4f} {sigma:8.4f}  ! OPLS type {opls_type}\n")
        else:
            f.write("! No atom types found\n")

def write_lib_file(mol_name, atoms, atom_charges, atom_types, atom_names, bonds, output_file):
    """Write .lib file with residue definition"""
    with open(output_file, 'w') as f:
        f.write(f"! Library file for {mol_name}\n")
        f.write(f"! Generated from ffld_server output\n\n")
        
        # Calculate total charge
        total_charge = sum(atom_charges.values()) if atom_charges else 0.0
        
        f.write(f"RESI {mol_name}    {total_charge:6.3f}\n")
        f.write("GROUP\n")
        
        # Write atom definitions - match PDB order with ffld charges
        for atom in atoms:
            atom_num = atom['num']
            opls_type = atom_types.get(atom_num, 0)
            charge = atom_charges.get(atom_num, 0.0)
            
            f.write(f"ATOM {atom['name']:<4s} {opls_type:<4d} {charge:8.4f}\n")
        
        # Add bond connectivity based on parsed bonds
        f.write("\nBOND\n")
        if bonds:
            for bond in bonds:
                f.write(f"BOND {bond['atom1']} {bond['atom2']}\n")
        else:
            f.write("! No bond information found - add manually based on CONECT records\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python ffld_to_prm.py <ffld_log_file> <pdb_file>")
        sys.exit(1)
    
    ffld_file = sys.argv[1]
    pdb_file = sys.argv[2]
    
    if not os.path.exists(ffld_file):
        print(f"Error: FFLD file {ffld_file} not found")
        sys.exit(1)
    
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file {pdb_file} not found")
        sys.exit(1)
    
    # Parse input files
    print("Parsing PDB file...")
    atoms = parse_pdb(pdb_file)
    
    print("\nParsing FFLD file...")
    mol_name, total_charge, atom_charges, atom_types, bonds, atom_names = parse_ffld_log(ffld_file)
    
    # Generate output filenames
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    prm_file = f"{base_name}.prm"
    lib_file = f"{base_name}.lib"
    
    # Write output files
    print(f"\nWriting {prm_file}...")
    write_prm_file(mol_name, atom_types, bonds, prm_file)
    
    print(f"Writing {lib_file}...")
    write_lib_file(mol_name, atoms, atom_charges, atom_types, atom_names, bonds, lib_file)
    
    print(f"\nGenerated {prm_file} and {lib_file}")
    print(f"Molecule: {mol_name}")
    print(f"Total charge: {total_charge}")
    print(f"Number of atoms: {len(atoms)}")
    print(f"Number of atom types found: {len(set(atom_types.values()))}")
    print(f"Number of bonds found: {len(bonds)}")
    
    # Show charges that will be written to lib file
    print("\nCharges to be written to .lib file:")
    for atom in atoms:
        charge = atom_charges.get(atom['num'], 0.0)
        print(f"  {atom['num']} {atom['name']}: {charge:8.4f}")

if __name__ == "__main__":
    main()