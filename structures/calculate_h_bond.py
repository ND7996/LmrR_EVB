from pymol import cmd

# Select residues
residues = [7,9,10,11,14,16,17,18,21,87,88,91,92,94,95,96,98,99,100,103]
cmd.select("my_residues", "resi " + "+".join(map(str, residues)))

# Calculate H-bonds
cmd.distance("hbonds_all", "(my_residues and (elem N,O))", "(all and (elem N,O))", 3.5, mode=2)
cmd.distance("hbonds_internal", "(my_residues and (elem N,O))", "(my_residues and (elem N,O))", 3.5, mode=2)

# Visualization
cmd.show("sticks", "my_residues")
cmd.color("cyan", "my_residues")
cmd.set("dash_color", "yellow", "hbonds_all")
cmd.set("dash_color", "red", "hbonds_internal")
cmd.set("dash_width", 3)
cmd.hide("labels")

# Zoom to fit everything
cmd.zoom("all", buffer=10)
cmd.center("my_residues")

# Output directory
output_dir = "D:/PhD_Thesis/LmrR_EVB/structures/"

# Save files
cmd.save(output_dir + "hbonds_session.pse")
cmd.save(output_dir + "my_residues.pdb", "my_residues")

# Save full image
cmd.ray(3000, 3000)
cmd.png(output_dir + "hbonds_view.png", dpi=300)

# Get H-bond information properly
print("Extracting H-bond pairs...")

# Function to get distance measurements
def get_hbond_pairs(obj_name):
    pairs = []
    cmd.iterate_state(1, obj_name, 
        "pairs.append((chain, resi, resn, name, ID))", 
        space={'pairs': pairs})
    return pairs

# Get internal H-bonds
internal_atoms = get_hbond_pairs("hbonds_internal")

# Get all H-bonds
all_atoms = get_hbond_pairs("hbonds_all")

# Alternative method: Find H-bonds manually
def find_hbonds_manual(selection1, selection2, cutoff=3.5):
    hbonds = []
    
    # Get donor atoms (N with H, O with H)
    donors = []
    cmd.iterate(f"({selection1}) and (elem N,O)",
                "donors.append((chain, resi, resn, name))",
                space={'donors': donors})
    
    # Get acceptor atoms (N, O)
    acceptors = []
    cmd.iterate(f"({selection2}) and (elem N,O)",
                "acceptors.append((chain, resi, resn, name))",
                space={'acceptors': acceptors})
    
    # Check distances
    for d in donors:
        d_sel = f"chain {d[0]} and resi {d[1]} and name {d[3]}"
        for a in acceptors:
            a_sel = f"chain {a[0]} and resi {a[1]} and name {a[3]}"
            
            # Skip same atom
            if d == a:
                continue
                
            dist = cmd.get_distance(d_sel, a_sel)
            if dist > 0 and dist <= cutoff:
                hbonds.append((d, a, dist))
    
    return hbonds

# Find internal H-bonds
print("Calculating internal H-bonds...")
internal_hbonds = find_hbonds_manual("my_residues", "my_residues", 3.5)

# Find all H-bonds
print("Calculating all H-bonds...")
all_hbonds = find_hbonds_manual("my_residues", "all", 3.5)

# Create detailed report
report = open(output_dir + "hbonds_detailed_report.txt", "w")
report.write("="*90 + "\n")
report.write("HYDROGEN BOND ANALYSIS - DETAILED REPORT\n")
report.write("="*90 + "\n\n")

report.write("SELECTED RESIDUES:\n")
report.write("-"*90 + "\n")
report.write(str(residues) + "\n\n")

report.write("PARAMETERS:\n")
report.write("-"*90 + "\n")
report.write("H-bond distance cutoff: 3.5 Angstroms\n")
report.write("Donor/Acceptor atoms: N, O (Nitrogen, Oxygen)\n\n")

# Internal H-bonds
report.write("="*90 + "\n")
report.write("H-BONDS BETWEEN SELECTED RESIDUES (INTERNAL):\n")
report.write("="*90 + "\n")
report.write(f"Total internal H-bonds found: {len(internal_hbonds)}\n\n")

if internal_hbonds:
    report.write(f"{'No.':<5} {'Donor':<30} {'Acceptor':<30} {'Distance (Å)':<12}\n")
    report.write("-"*90 + "\n")
    for i, (donor, acceptor, dist) in enumerate(internal_hbonds, 1):
        donor_str = f"{donor[2]}{donor[1]}:{donor[3]} (Chain {donor[0]})"
        acceptor_str = f"{acceptor[2]}{acceptor[1]}:{acceptor[3]} (Chain {acceptor[0]})"
        report.write(f"{i:<5} {donor_str:<30} {acceptor_str:<30} {dist:<12.2f}\n")
else:
    report.write("No internal H-bonds detected.\n")

report.write("\n" + "="*90 + "\n")
report.write("H-BONDS WITH ALL OTHER ATOMS:\n")
report.write("="*90 + "\n")
report.write(f"Total H-bonds found: {len(all_hbonds)}\n\n")

if all_hbonds:
    # Show first 100
    report.write(f"{'No.':<5} {'Donor':<30} {'Acceptor':<30} {'Distance (Å)':<12}\n")
    report.write("-"*90 + "\n")
    for i, (donor, acceptor, dist) in enumerate(all_hbonds[:100], 1):
        donor_str = f"{donor[2]}{donor[1]}:{donor[3]} (Chain {donor[0]})"
        acceptor_str = f"{acceptor[2]}{acceptor[1]}:{acceptor[3]} (Chain {acceptor[0]})"
        report.write(f"{i:<5} {donor_str:<30} {acceptor_str:<30} {dist:<12.2f}\n")
    
    if len(all_hbonds) > 100:
        report.write(f"\n... and {len(all_hbonds)-100} more H-bonds\n")

report.write("\n" + "="*90 + "\n")
report.write("FILES SAVED:\n")
report.write("-"*90 + "\n")
report.write("- hbonds_session.pse (Complete PyMOL session)\n")
report.write("- my_residues.pdb (Selected residues structure)\n")
report.write("- hbonds_view.png (Full visualization)\n")
report.write("- hbonds_detailed_report.txt (This report)\n")
report.write("="*90 + "\n")

report.close()

print("\n" + "="*90)
print("ANALYSIS COMPLETE!")
print("="*90)
print(f"Files saved to: {output_dir}")
print(f"Internal H-bonds: {len(internal_hbonds)}")
print(f"Total H-bonds: {len(all_hbonds)}")
print("="*90)