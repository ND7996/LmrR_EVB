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
                print(f"Parsed atom: {atom_num}, name: {atom_name}, element: {element}")
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
    
    # Extract atom types and charges
    atom_charges = {}
    atom_types = {}
    
    # Look for atom assignment section
    # Common pattern: "Atom  Type   Charge" or similar header
    lines = content.split('\n')
    in_atom_section = False
    
    for i, line in enumerate(lines):
        line_lower = line.lower()
        
        # Look for start of atom section
        if any(keyword in line_lower for keyword in ['atom', 'type', 'charge', 'opls']):
            if 'atom' in line_lower and ('type' in line_lower or 'charge' in line_lower):
                in_atom_section = True
                print(f"Found atom section header: {line.strip()}")
                continue
        
        # Process atom data lines
        if in_atom_section and line.strip():
            parts = line.split()
            if len(parts) >= 4:
                try:
                    atom_num = int(parts[0])
                    atom_name = parts[1]
                    
                    # Look for OPLS type and charge
                    charge = None
                    opls_type = None
                    
                    for part in parts[2:]:
                        if 'opls_' in part:
                            opls_type = int(part.split('_')[1])
                        elif '.' in part or '-' in part or '+' in part:
                            try:
                                charge = float(part)
                            except ValueError:
                                pass
                    
                    if opls_type is not None and charge is not None:
                        atom_charges[atom_num] = charge
                        atom_types[atom_name] = opls_type
                        print(f"Found atom: {atom_num}, {atom_name}, opls_{opls_type}, charge: {charge}")
                        
                except (ValueError, IndexError):
                    continue
            else:
                # If we hit a line that doesn't look like atom data, maybe we're done
                if len(parts) < 2 or not parts[0].isdigit():
                    in_atom_section = False
    
    # If we didn't find atoms using the section method, try pattern matching
    if not atom_charges:
        print("Trying pattern matching for atoms...")
        patterns = [
            r'^\s*(\d+)\s+(\S+)\s+opls_(\d+)\s+([-+]?\d*\.?\d+)',
            r'Atom\s+(\d+).*?opls_(\d+).*?([-+]?\d*\.?\d+)',
            r'(\d+)\s+(\S+)\s+([-+]?\d*\.?\d+)\s+opls_(\d+)',
        ]
        
        for pattern in patterns:
            matches = re.findall(pattern, content, re.MULTILINE | re.IGNORECASE)
            for match in matches:
                if len(match) == 4:
                    atom_num = int(match[0])
                    atom_name = match[1]
                    opls_type = int(match[2] if match[2].isdigit() else match[3])
                    charge = float(match[3] if '.' in match[3] or 'e' in match[3].lower() else match[2])
                    
                    atom_charges[atom_num] = charge
                    atom_types[atom_name] = opls_type
                    print(f"Pattern found atom: {atom_num}, {atom_name}, opls_{opls_type}, charge: {charge}")
    
    # Parse bond information
    bonds = []
    bond_patterns = [
        r'BOND\s+(\S+)\s+(\S+)\s+([\d.]+)\s+([\d.]+)',
        r'(\S+)-(\S+)\s+[\d.]+\s+([\d.]+)\s+([\d.]+)',
    ]
    
    for pattern in bond_patterns:
        matches = re.findall(pattern, content)
        for match in matches:
            if len(match) == 4:
                atom1, atom2, length, force_const = match
                bonds.append({
                    'atom1': atom1,
                    'atom2': atom2,
                    'length': float(length),
                    'force_const': float(force_const)
                })
                print(f"Found bond: {atom1}-{atom2}, length: {length}, force_const: {force_const}")
    
    return mol_name, total_charge, atom_charges, atom_types, bonds

def write_prm_file(mol_name, atom_types, bonds, output_file):
    """Write .prm file with OPLS parameters"""
    with open(output_file, 'w') as f:
        f.write(f"! OPLS parameters for {mol_name}\n")
        f.write(f"! Generated from ffld_server output\n\n")
        
        # Write bond parameters
        if bonds:
            f.write("BONDS\n")
            for bond in bonds:
                type1 = atom_types.get(bond['atom1'], 0)
                type2 = atom_types.get(bond['atom2'], 0)
                kb = bond['force_const'] * 100.0  # Convert to appropriate units
                r0 = bond['length']
                f.write(f"{type1:>4d} {type2:>4d} {kb:8.2f} {r0:8.4f}  ! {bond['atom1']}-{bond['atom2']}\n")
        else:
            f.write("BONDS\n")
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

def write_lib_file(mol_name, atoms, atom_charges, atom_types, output_file):
    """Write .lib file with residue definition"""
    with open(output_file, 'w') as f:
        f.write(f"! Library file for {mol_name}\n")
        f.write(f"! Generated from ffld_server output\n\n")
        
        # Calculate total charge
        total_charge = sum(atom_charges.values()) if atom_charges else 0.0
        
        f.write(f"RESI {mol_name}    {total_charge:6.3f}\n")
        f.write("GROUP\n")
        
        # Write atom definitions
        for atom in atoms:
            atom_name = atom['name']
            opls_type = atom_types.get(atom_name, 0)
            charge = atom_charges.get(atom['num'], 0.0)
            
            f.write(f"ATOM {atom_name:<4s} {opls_type:<4d} {charge:8.4f}\n")
        
        f.write("\n")
        f.write("! BOND connectivity - add manually based on CONECT records\n")
        f.write("! Example: BOND C1 C2\n")

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
    mol_name, total_charge, atom_charges, atom_types, bonds = parse_ffld_log(ffld_file)
    
    # Generate output filenames
    base_name = os.path.splitext(os.path.basename(pdb_file))[0]
    prm_file = f"{base_name}.prm"
    lib_file = f"{base_name}.lib"
    
    # Write output files
    print(f"\nWriting {prm_file}...")
    write_prm_file(mol_name, atom_types, bonds, prm_file)
    
    print(f"Writing {lib_file}...")
    write_lib_file(mol_name, atoms, atom_charges, atom_types, lib_file)
    
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