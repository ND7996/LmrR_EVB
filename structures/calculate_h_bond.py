from pymol import cmd
import numpy as np

# Load your PDB file first
pdb_file = "/home/hp/nayanika/github/LmrR_EVB/structures/WT_2.pdb"
cmd.load(pdb_file, "my_structure")

print(f"Loaded structure from: {pdb_file}")

# Select residues
residues = [7,9,10,11,14,16,17,18,21,87,88,91,92,94,95,96,98,99,100,103]
cmd.select("my_residues", "resi " + "+".join(map(str, residues)))

# Visualization setup
cmd.show("sticks", "my_residues")
cmd.color("cyan", "my_residues")

# Output directory
output_dir = "/home/hp/nayanika/github/LmrR_EVB/structures/"

def find_hbonds_simple(selection1, selection2, cutoff=3.5):
    """
    Simple H-bond detection using only distance criteria
    """
    hbonds = []
    
    # Get all potential H-bond atoms from both selections
    atoms1 = []
    cmd.iterate(f"({selection1}) and (elem N,O)",
                "atoms1.append((model, chain, resi, resn, name, index, elem))",
                space={'atoms1': atoms1})
    
    atoms2 = []
    cmd.iterate(f"({selection2}) and (elem N,O)", 
                "atoms2.append((model, chain, resi, resn, name, index, elem))",
                space={'atoms2': atoms2})
    
    print(f"Found {len(atoms1)} atoms in selection1, {len(atoms2)} atoms in selection2")
    
    for atom1 in atoms1:
        model1, chain1, resi1, resn1, name1, idx1, elem1 = atom1
        sel1 = f"index {idx1}"
        
        for atom2 in atoms2:
            model2, chain2, resi2, resn2, name2, idx2, elem2 = atom2
            
            if idx1 == idx2:
                continue
                
            # Skip same residue if desired
            if chain1 == chain2 and resi1 == resi2:
                continue
                
            sel2 = f"index {idx2}"
            
            try:
                dist = cmd.get_distance(sel1, sel2)
                if 0 < dist <= cutoff:
                    hbonds.append((atom1, atom2, dist))
            except:
                continue
                
    return hbonds

def find_hbonds_pymol_builtin(selection1, selection2, cutoff=3.5):
    """
    Use PyMOL's built-in hbond command
    """
    hbonds = []
    
    # Create a temporary object for H-bonds
    temp_name = "temp_hbonds"
    cmd.delete(temp_name)
    
    # Use PyMOL's hbond command
    cmd.distance(temp_name, selection1, selection2, cutoff, mode=2)
    
    # Get the distance objects created
    dist_count = cmd.count_discrete(temp_name)
    print(f"PyMOL built-in found {dist_count} H-bond distance objects")
    
    if dist_count > 0:
        # Iterate through distance objects to get details
        for i in range(1, dist_count + 1):
            try:
                # Get the distance measurement
                dist_obj = f"{temp_name}_{i:04d}"
                dist = cmd.get_distance(dist_obj)
                
                # Get the atoms involved
                atoms = cmd.get_model(dist_obj).atom
                if len(atoms) >= 2:
                    atom1 = atoms[0]
                    atom2 = atoms[1]
                    
                    donor_info = (atom1.model, atom1.chain, str(atom1.resi), atom1.resn, atom1.name, atom1.index, atom1.elem)
                    acceptor_info = (atom2.model, atom2.chain, str(atom2.resi), atom2.resn, atom2.name, atom2.index, atom2.elem)
                    
                    hbonds.append((donor_info, acceptor_info, dist))
                    
            except Exception as e:
                continue
                
        # Clean up
        cmd.delete(temp_name)
    
    return hbonds

def find_hbonds_robust(selection1, selection2, cutoff=3.5):
    """
    Robust H-bond detection that tries multiple methods
    """
    print("Trying simple distance-based method...")
    hbonds = find_hbonds_simple(selection1, selection2, cutoff)
    
    if len(hbonds) == 0:
        print("Simple method found no H-bonds, trying PyMOL built-in method...")
        hbonds = find_hbonds_pymol_builtin(selection1, selection2, cutoff)
    
    return hbonds

# Calculate H-bonds
print("\nCalculating hydrogen bonds...")
print("="*90)

# Use robust method
print("Finding internal H-bonds (between selected residues)...")
internal_hbonds = find_hbonds_robust("my_residues", "my_residues", 3.5)

print("Finding all H-bonds (selected residues to all atoms)...")
all_hbonds = find_hbonds_robust("my_residues", "all", 3.5)

# Remove internal from all to get external only
external_hbonds = [hb for hb in all_hbonds if hb not in internal_hbonds]

print(f"Found {len(internal_hbonds)} internal H-bonds")
print(f"Found {len(all_hbonds)} total H-bonds")
print(f"Found {len(external_hbonds)} external H-bonds")

# If still no H-bonds, try a more permissive approach
if len(all_hbonds) == 0:
    print("\nTrying more permissive H-bond detection...")
    print("Using larger cutoff (4.0 Å) and including potential donors/acceptors...")
    
    # Try with larger cutoff
    internal_hbonds = find_hbonds_simple("my_residues", "my_residues", 4.0)
    all_hbonds = find_hbonds_simple("my_residues", "all", 4.0)
    external_hbonds = [hb for hb in all_hbonds if hb not in internal_hbonds]
    
    print(f"With 4.0 Å cutoff: {len(internal_hbonds)} internal, {len(external_hbonds)} external H-bonds")

# Visualize H-bonds
print("\nCreating visualizations...")

# Create distance objects for visualization
if internal_hbonds:
    for i, (atom1, atom2, dist) in enumerate(internal_hbonds):
        sel1 = f"index {atom1[5]}"
        sel2 = f"index {atom2[5]}"
        cmd.distance(f"internal_hb_{i}", sel1, sel2)
        cmd.set("dash_color", "red", f"internal_hb_{i}")

if external_hbonds:
    max_external = min(50, len(external_hbonds))
    for i, (atom1, atom2, dist) in enumerate(external_hbonds[:max_external]):
        sel1 = f"index {atom1[5]}"
        sel2 = f"index {atom2[5]}"
        cmd.distance(f"external_hb_{i}", sel1, sel2)
        cmd.set("dash_color", "yellow", f"external_hb_{i}")

cmd.set("dash_width", 2)
cmd.hide("labels")

# Zoom and center
cmd.zoom("my_residues", buffer=5)
cmd.center("my_residues")

# Create detailed report
print("\nGenerating report...")
report = open(output_dir + "hbonds_detailed_report.txt", "w")
report.write("="*90 + "\n")
report.write("HYDROGEN BOND ANALYSIS - DETAILED REPORT\n")
report.write("="*90 + "\n\n")

report.write(f"INPUT STRUCTURE: {pdb_file}\n")
report.write(f"STRUCTURE NAME: my_structure\n\n")

report.write("SELECTED RESIDUES:\n")
report.write("-"*90 + "\n")
report.write(str(residues) + "\n\n")

report.write("PARAMETERS:\n")
report.write("-"*90 + "\n")
report.write("H-bond distance cutoff: 3.5 Angstroms\n")
report.write("Donor/Acceptor atoms: N, O (Nitrogen, Oxygen)\n")
report.write("Note: Intra-residue H-bonds are excluded\n\n")

# Internal H-bonds
report.write("="*90 + "\n")
report.write("H-BONDS BETWEEN SELECTED RESIDUES (INTERNAL):\n")
report.write("="*90 + "\n")
report.write(f"Total internal H-bonds found: {len(internal_hbonds)}\n\n")

if internal_hbonds:
    report.write(f"{'No.':<5} {'Residue 1':<25} {'Residue 2':<25} {'Distance (Å)':<12}\n")
    report.write("-"*90 + "\n")
    for i, (atom1, atom2, dist) in enumerate(sorted(internal_hbonds, key=lambda x: x[2]), 1):
        res1_str = f"{atom1[3]}{atom1[2]}:{atom1[4]} (Ch.{atom1[1]})"
        res2_str = f"{atom2[3]}{atom2[2]}:{atom2[4]} (Ch.{atom2[1]})"
        report.write(f"{i:<5} {res1_str:<25} {res2_str:<25} {dist:<12.2f}\n")
else:
    report.write("No internal H-bonds detected.\n")

# External H-bonds
report.write("\n" + "="*90 + "\n")
report.write("H-BONDS WITH EXTERNAL ATOMS:\n")
report.write("="*90 + "\n")
report.write(f"Total external H-bonds found: {len(external_hbonds)}\n\n")

if external_hbonds:
    report.write(f"{'No.':<5} {'Selected Residue':<25} {'Partner Atom':<25} {'Distance (Å)':<12}\n")
    report.write("-"*90 + "\n")
    for i, (atom1, atom2, dist) in enumerate(sorted(external_hbonds, key=lambda x: x[2])[:100], 1):
        res1_str = f"{atom1[3]}{atom1[2]}:{atom1[4]} (Ch.{atom1[1]})"
        res2_str = f"{atom2[3]}{atom2[2]}:{atom2[4]} (Ch.{atom2[1]})"
        report.write(f"{i:<5} {res1_str:<25} {res2_str:<25} {dist:<12.2f}\n")
    
    if len(external_hbonds) > 100:
        report.write(f"\n... and {len(external_hbonds)-100} more external H-bonds\n")
else:
    report.write("No external H-bonds detected.\n")

# Summary by residue
report.write("\n" + "="*90 + "\n")
report.write("SUMMARY BY RESIDUE:\n")
report.write("="*90 + "\n\n")

residue_hbonds = {}
for hbond in all_hbonds:
    atom1, atom2, dist = hbond
    res_key = f"{atom1[3]}{atom1[2]} (Chain {atom1[1]})"
    if res_key not in residue_hbonds:
        residue_hbonds[res_key] = {'internal': 0, 'external': 0}
    
    # Check if it's internal or external
    is_internal = any(str(atom2[2]) == str(r) for r in residues)
    if is_internal:
        residue_hbonds[res_key]['internal'] += 1
    else:
        residue_hbonds[res_key]['external'] += 1

if residue_hbonds:
    report.write(f"{'Residue':<20} {'Internal H-bonds':<20} {'External H-bonds':<20} {'Total':<10}\n")
    report.write("-"*90 + "\n")
    for res, counts in sorted(residue_hbonds.items()):
        total = counts['internal'] + counts['external']
        report.write(f"{res:<20} {counts['internal']:<20} {counts['external']:<20} {total:<10}\n")
else:
    report.write("No H-bond data available for residue summary.\n")

# Additional diagnostics
report.write("\n" + "="*90 + "\n")
report.write("DIAGNOSTIC INFORMATION:\n")
report.write("="*90 + "\n")

# Count atoms in selection
num_atoms = cmd.count_atoms("my_residues")
num_n_o_atoms = cmd.count_atoms("my_residues and (elem N,O)")

report.write(f"Total atoms in selected residues: {num_atoms}\n")
report.write(f"N and O atoms in selected residues: {num_n_o_atoms}\n")
report.write(f"Total atoms in structure: {cmd.count_atoms('all')}\n")

report.write("\n" + "="*90 + "\n")
report.write("POSSIBLE REASONS FOR NO H-BONDS:\n")
report.write("-"*90 + "\n")
report.write("1. Structure may not have hydrogen atoms\n")
report.write("2. Distance cutoff may be too small\n")
report.write("3. No suitable donor-acceptor pairs in the selection\n")
report.write("4. Residues may be too far apart\n")
report.write("5. Protein may be in a conformation with few H-bonds\n")

report.write("\n" + "="*90 + "\n")
report.write("FILES TO BE SAVED:\n")
report.write("-"*90 + "\n")
report.write("- hbonds_session.pse (Complete PyMOL session)\n")
report.write("- my_residues.pdb (Selected residues structure)\n")
report.write("- hbonds_view.png (Full visualization)\n")
report.write("- hbonds_detailed_report.txt (This report)\n")
report.write("="*90 + "\n")

report.close()

# Save files
print("\nSaving files...")
cmd.save(output_dir + "hbonds_session.pse")
cmd.save(output_dir + "my_residues.pdb", "my_residues")
cmd.ray(3000, 3000)
cmd.png(output_dir + "hbonds_view.png", dpi=300)

# Print summary
print("\n" + "="*90)
print("ANALYSIS COMPLETE!")
print("="*90)
print(f"Input structure: {pdb_file}")
print(f"Files saved to: {output_dir}")
print(f"Internal H-bonds (between selected residues): {len(internal_hbonds)}")
print(f"External H-bonds (selected to other atoms): {len(external_hbonds)}")
print(f"Total H-bonds: {len(all_hbonds)}")
print("="*90)

if len(all_hbonds) == 0:
    print("\nDIAGNOSTIC INFORMATION:")
    print(f"Total atoms in selection: {num_atoms}")
    print(f"N and O atoms in selection: {num_n_o_atoms}")
    print("\nPossible reasons for no H-bonds:")
    print("1. Structure may not have hydrogen atoms")
    print("2. Try: 'show spheres, my_residues and (elem N,O)' to see potential H-bond atoms")
    print("3. Check if residues are close enough for H-bonds")

print("\nVisualization:")
print(f"  - Red dashes: Internal H-bonds ({len(internal_hbonds)} bonds)")
print(f"  - Yellow dashes: External H-bonds (first {min(50, len(external_hbonds))} of {len(external_hbonds)} shown)")
print("="*90)