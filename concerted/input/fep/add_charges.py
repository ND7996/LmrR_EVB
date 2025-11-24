#!/usr/bin/env python3
"""
Script to calculate total charges for state 1 and state 2 from Q map and charge parameters
"""

# Charge dictionaries from your provided data
ala_charges = {
    'N': -0.500, 'HN': 0.300, 'CA': 0.140, 'HA': 0.060,
    'CB': -0.180, 'HB1': 0.060, 'HB2': 0.060, 'HB3': 0.060,
    'C': 0.500, 'O': -0.500
}

imi_charges = {
    'O1': -0.530000, 'O2': -0.440000, 'N1': -0.943900, 'N2': -0.139000,
    'C1': -0.005000, 'C2': 0.163900, 'C3': -0.115000, 'C4': -0.115000,
    'C5': -0.115000, 'C6': 0.520000, 'C7': -0.115000, 'C8': -0.115000,
    'C9': 0.122500, 'C10': 0.263400, 'C11': -0.115000, 'C12': -0.115000,
    'C13': -0.120000, 'C14': -0.120000, 'C15': -0.180000,
    'H1': 0.450000, 'H2': 0.360000, 'H3': 0.360000, 'H4': 0.480000,
    'H5': 0.060000, 'H6': 0.060000, 'H7': 0.060000, 'H8': 0.115000,
    'H9': 0.115000, 'H10': 0.115000, 'H11': 0.115000, 'H12': 0.273100,
    'H13': 0.115000, 'H14': 0.115000, 'H15': 0.060000, 'H16': 0.060000,
    'H17': 0.060000, 'H18': 0.060000, 'H19': 0.060000, 'H20': 0.060000,
    'H21': 0.060000
}

paf_charges = {
    'O1': -0.530000, 'O2': -0.440000, 'N1': -0.943900, 'N2': -0.900000,
    'C1': -0.005000, 'C2': 0.163900, 'C3': -0.115000, 'C4': -0.115000,
    'C5': -0.115000, 'C6': 0.520000, 'C7': -0.115000, 'C8': -0.115000,
    'C9': 0.180000,
    'H1': 0.450000, 'H2': 0.360000, 'H3': 0.360000, 'H4': 0.360000,
    'H5': 0.360000, 'H6': 0.060000, 'H7': 0.060000, 'H8': 0.060000,
    'H9': 0.115000, 'H10': 0.115000, 'H11': 0.115000, 'H12': 0.115000
}

ind_charges = {
    'N1': -0.478600, 'C1': 0.286000, 'C2': 0.123300, 'C3': -0.256700,
    'C4': -0.354100, 'C5': -0.224300, 'C6': -0.029300, 'C7': -0.088800,
    'C8': -0.222200, 'C9': -0.065000,
    'H1': 0.357700, 'H2': 0.137700, 'H3': 0.179600, 'H4': 0.169900,
    'H5': 0.135900, 'H6': 0.148900, 'H7': 0.060000, 'H8': 0.060000,
    'H9': 0.060000
}

prd_charges = {
    'N1': -0.478600, 'C1': 0.286000, 'C2': 0.123300, 'C3': -0.174900,
    'C4': -0.354100, 'C5': -0.224300, 'C6': -0.029300, 'C7': -0.088800,
    'C8': -0.222200, 'C9': -0.065000, 'C10': -0.004100, 'C11': -0.120000,
    'C12': 0.145000, 'O1': -0.683000, 'C13': -0.120000, 'C14': -0.120000,
    'C15': -0.180000,
    'H1': 0.357700, 'H2': 0.179600, 'H3': 0.169900, 'H4': 0.135900,
    'H5': 0.148900, 'H6': 0.060000, 'H7': 0.060000, 'H8': 0.060000,
    'H9': 0.060000, 'H10': 0.060000, 'H11': 0.060000, 'H12': 0.060000,
    'H13': 0.060000, 'H14': 0.418000, 'H15': 0.060000, 'H16': 0.060000,
    'H17': 0.060000, 'H18': 0.060000, 'H19': 0.060000, 'H20': 0.060000,
    'H21': 0.060000
}

ho2_charges = {
    'O': -0.834, 'H1': 0.417, 'H2': 0.417
}

def get_charge(resname_atom):
    """
    Extract charge for a given resname.atom combination
    """
    if '.' in resname_atom:
        resname, atom = resname_atom.split('.')
    else:
        # If no resname, try to find atom in any dictionary
        resname, atom = '', resname_atom
    
    charge_dicts = {
        'ALA': ala_charges, 'IMI': imi_charges, 'PAF': paf_charges,
        'IND': ind_charges, 'PRD': prd_charges, 'HO2': ho2_charges
    }
    
    # Try exact match first
    for dict_name, charge_dict in charge_dicts.items():
        if atom in charge_dict:
            return charge_dict[atom]
    
    # If not found, try without resname prefix
    simple_atom = atom.replace('IMI.', '').replace('PAF.', '').replace('IND.', '').replace('PRD.', '').replace('HO2.', '')
    for charge_dict in charge_dicts.values():
        if simple_atom in charge_dict:
            return charge_dict[simple_atom]
    
    print(f"Warning: Charge not found for {resname_atom}")
    return 0.0

def calculate_charges(qmap_file):
    """
    Calculate total charges for state 1 and state 2
    """
    state1_total = 0.0
    state2_total = 0.0
    state1_atoms = []
    state2_atoms = []
    
    print("Calculating charges from Q map...")
    print("-" * 80)
    
    with open(qmap_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
                
            parts = line.split()
            if len(parts) < 4:
                continue
                
            # Skip header lines
            if 'PDB ID' in line or 'STATE 1' in line or 'LIB ID' in line:
                continue
                
            # Handle both 'q' and 'n' entries
            if parts[0] in ['q', 'n']:
                state1_atom = parts[2]  # STATE 1 atom name
                state2_atom = parts[3]  # STATE 2 atom name
                
                charge1 = get_charge(state1_atom)
                charge2 = get_charge(state2_atom)
                
                state1_total += charge1
                state2_total += charge2
                state1_atoms.append((state1_atom, charge1))
                state2_atoms.append((state2_atom, charge2))
                
                print(f"{line_num:3d}: {state1_atom:15s} {charge1:8.4f}  ->  {state2_atom:15s} {charge2:8.4f}")
    
    print("-" * 80)
    print(f"State 1 total charge: {state1_total:.6f}")
    print(f"State 2 total charge: {state2_total:.6f}")
    print(f"State 1 integer check: {state1_total} - round({state1_total}) = {state1_total - round(state1_total):.6f}")
    print(f"State 2 integer check: {state2_total} - round({state2_total}) = {state2_total - round(state2_total):.6f}")
    
    # Detailed breakdown
    print("\nDetailed breakdown:")
    print("\nState 1 atoms and charges:")
    for atom, charge in state1_atoms:
        print(f"  {atom:15s}: {charge:8.4f}")
    
    print("\nState 2 atoms and charges:")
    for atom, charge in state2_atoms:
        print(f"  {atom:15s}: {charge:8.4f}")
    
    return state1_total, state2_total

def main():
    qmap_file = "FC_concerted.qmap"  # Change this to your Q map file name
    
    try:
        state1_charge, state2_charge = calculate_charges(qmap_file)
        
        print(f"\nSummary:")
        print(f"State 1 (IMI + IND + HO2 + ALA): {state1_charge:.6f}")
        print(f"State 2 (PAF + PRD + ALA): {state2_charge:.6f}")
        
        # Check if charges are integers
        if abs(state1_charge - round(state1_charge)) > 0.001:
            print(f"ERROR: State 1 charge is not integer! Difference: {state1_charge - round(state1_charge):.6f}")
        
        if abs(state2_charge - round(state2_charge)) > 0.001:
            print(f"ERROR: State 2 charge is not integer! Difference: {state2_charge - round(state2_charge):.6f}")
            
    except FileNotFoundError:
        print(f"Error: Q map file '{qmap_file}' not found.")
        print("Please make sure the file exists in the current directory.")
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()