from pymol import cmd
import os

# Configuration
residues = [7, 9, 10, 11, 14, 16, 17, 18, 21, 87, 88, 91, 92, 94, 95, 96, 98, 99, 100, 103]
output_dir = "D:/PhD_Thesis/LmrR_EVB/structures/"
hbond_cutoff = 3.5  # Angstroms

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

print("="*90)
print("STARTING H-BOND ANALYSIS")
print("="*90)

# Select residues
cmd.select("my_residues", "resi " + "+".join(map(str, residues)))
print(f"Selected {len(residues)} residues")

# Calculate H-bonds using PyMOL's distance function
cmd.distance("hbonds_all", "(my_residues and (elem N,O))", "(all and (elem N,O))", 
             hbond_cutoff, mode=2)
cmd.distance("hbonds_internal", "(my_residues and (elem N,O))", 
             "(my_residues and (elem N,O))", hbond_cutoff, mode=2)

# Visualization settings
cmd.show("sticks", "my_residues")
cmd.color("cyan", "my_residues")
cmd.set("dash_color", "yellow", "hbonds_all")
cmd.set("dash_color", "red", "hbonds_internal")
cmd.set("dash_width", 3)
cmd.hide("labels")
cmd.zoom("my_residues", buffer=5)

print("Visualization prepared")


def find_hbonds_detailed(selection1, selection2, cutoff=3.5):
    """
    Find hydrogen bonds between two selections.
    Returns list of tuples: (donor_info, acceptor_info, distance)
    """
    hbonds = []
    
    # Get all N and O atoms from both selections
    atoms1 = []
    cmd.iterate(f"({selection1}) and (elem N,O)",
                "atoms1.append((chain, resi, resn, name, elem))",
                space={'atoms1': atoms1})
    
    atoms2 = []
    cmd.iterate(f"({selection2}) and (elem N,O)",
                "atoms2.append((chain, resi, resn, name, elem))",
                space={'atoms2': atoms2})
    
    print(f"  Scanning {len(atoms1)} atoms in selection1 vs {len(atoms2)} atoms in selection2...")
    
    # Check all atom pairs
    for atom1 in atoms1:
        sel1 = f"chain {atom1[0]} and resi {atom1[1]} and name {atom1[3]}"
        
        for atom2 in atoms2:
            # Skip same atom
            if atom1 == atom2:
                continue
            
            # Skip same residue for internal (optional, remove if you want intra-residue bonds)
            if selection1 == selection2 and atom1[1] == atom2[1]:
                continue
                
            sel2 = f"chain {atom2[0]} and resi {atom2[1]} and name {atom2[3]}"
            
            # Measure distance
            try:
                dist = cmd.distance("temp_dist", sel1, sel2)
                cmd.delete("temp_dist")
                
                if 0 < dist <= cutoff:
                    hbonds.append((atom1, atom2, dist))
            except:
                continue
    
    return hbonds


# Find internal H-bonds (within selected residues)
print("\nCalculating internal H-bonds...")
internal_hbonds = find_hbonds_detailed("my_residues", "my_residues", hbond_cutoff)
print(f"Found {len(internal_hbonds)} internal H-bonds")

# Find all H-bonds (selected residues with all atoms)
print("\nCalculating all H-bonds...")
all_hbonds = find_hbonds_detailed("my_residues", "all", hbond_cutoff)
print(f"Found {len(all_hbonds)} total H-bonds")

# Separate external H-bonds (all - internal)
external_hbonds = [hb for hb in all_hbonds if hb not in internal_hbonds]
print(f"Found {len(external_hbonds)} external H-bonds")


# Generate detailed text report
report_file = os.path.join(output_dir, "hbonds_detailed_report.txt")
with open(report_file, "w") as report:
    report.write("="*100 + "\n")
    report.write("HYDROGEN BOND ANALYSIS - DETAILED REPORT\n")
    report.write("="*100 + "\n\n")
    
    # Parameters section
    report.write("ANALYSIS PARAMETERS:\n")
    report.write("-"*100 + "\n")
    report.write(f"Selected Residues: {residues}\n")
    report.write(f"Total residues analyzed: {len(residues)}\n")
    report.write(f"H-bond distance cutoff: {hbond_cutoff} Å\n")
    report.write(f"Donor/Acceptor atoms: Nitrogen (N) and Oxygen (O)\n\n")
    
    # Summary section
    report.write("SUMMARY:\n")
    report.write("-"*100 + "\n")
    report.write(f"Internal H-bonds (within selected residues): {len(internal_hbonds)}\n")
    report.write(f"External H-bonds (with other residues): {len(external_hbonds)}\n")
    report.write(f"Total H-bonds: {len(all_hbonds)}\n\n")
    
    # Internal H-bonds section
    report.write("="*100 + "\n")
    report.write("INTERNAL HYDROGEN BONDS (WITHIN SELECTED RESIDUES)\n")
    report.write("="*100 + "\n\n")
    
    if internal_hbonds:
        report.write(f"{'No.':<6} {'Donor Residue':<20} {'Donor Atom':<15} "
                    f"{'Acceptor Residue':<20} {'Acceptor Atom':<15} {'Distance (Å)':<12}\n")
        report.write("-"*100 + "\n")
        
        for i, (donor, acceptor, dist) in enumerate(sorted(internal_hbonds, key=lambda x: x[2]), 1):
            donor_res = f"{donor[2]}{donor[1]} (Chain {donor[0]})"
            donor_atom = f"{donor[3]} ({donor[4]})"
            acceptor_res = f"{acceptor[2]}{acceptor[1]} (Chain {acceptor[0]})"
            acceptor_atom = f"{acceptor[3]} ({acceptor[4]})"
            
            report.write(f"{i:<6} {donor_res:<20} {donor_atom:<15} "
                        f"{acceptor_res:<20} {acceptor_atom:<15} {dist:<12.2f}\n")
    else:
        report.write("No internal H-bonds detected.\n")
    
    # External H-bonds section
    report.write("\n" + "="*100 + "\n")
    report.write("EXTERNAL HYDROGEN BONDS (WITH OTHER RESIDUES)\n")
    report.write("="*100 + "\n\n")
    
    if external_hbonds:
        report.write(f"{'No.':<6} {'Selected Residue':<20} {'Atom':<15} "
                    f"{'Partner Residue':<20} {'Atom':<15} {'Distance (Å)':<12}\n")
        report.write("-"*100 + "\n")
        
        for i, (atom1, atom2, dist) in enumerate(sorted(external_hbonds, key=lambda x: x[2]), 1):
            res1 = f"{atom1[2]}{atom1[1]} (Chain {atom1[0]})"
            atom1_name = f"{atom1[3]} ({atom1[4]})"
            res2 = f"{atom2[2]}{atom2[1]} (Chain {atom2[0]})"
            atom2_name = f"{atom2[3]} ({atom2[4]})"
            
            report.write(f"{i:<6} {res1:<20} {atom1_name:<15} "
                        f"{res2:<20} {atom2_name:<15} {dist:<12.2f}\n")
            
            # Limit output for very large lists
            if i >= 200:
                report.write(f"\n... and {len(external_hbonds)-200} more external H-bonds (truncated for readability)\n")
                break
    else:
        report.write("No external H-bonds detected.\n")
    
    # Files saved section
    report.write("\n" + "="*100 + "\n")
    report.write("OUTPUT FILES:\n")
    report.write("-"*100 + "\n")
    report.write(f"1. {report_file} - This detailed report\n")
    report.write(f"2. {os.path.join(output_dir, 'hbonds_session.pse')} - PyMOL session file\n")
    report.write(f"3. {os.path.join(output_dir, 'my_residues.pdb')} - Selected residues structure\n")
    report.write(f"4. {os.path.join(output_dir, 'hbonds_view.png')} - Visualization image\n")
    report.write("="*100 + "\n")

print(f"\nDetailed report saved to: {report_file}")

# Save PyMOL session and files
cmd.save(os.path.join(output_dir, "hbonds_session.pse"))
cmd.save(os.path.join(output_dir, "my_residues.pdb"), "my_residues")

# Save visualization
print("\nGenerating high-quality image...")
cmd.ray(3000, 3000)
cmd.png(os.path.join(output_dir, "hbonds_view.png"), dpi=300)

# Print summary to console
print("\n" + "="*90)
print("ANALYSIS COMPLETE!")
print("="*90)
print(f"Output directory: {output_dir}")
print(f"Internal H-bonds: {len(internal_hbonds)}")
print(f"External H-bonds: {len(external_hbonds)}")
print(f"Total H-bonds: {len(all_hbonds)}")
print("="*90)
print("\nAll files saved successfully!")