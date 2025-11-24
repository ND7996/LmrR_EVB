from pymol import cmd

def save_nearby_atoms_csv():
    origin = cmd.get_model('resn IMI and name N2').atom[0]
    ox, oy, oz = origin.coord
    cmd.select("nearby_atoms", "all within 15 of (resn IMI and name N2)")
    atoms = cmd.get_model("nearby_atoms and not (resn IMI and name N2)")
    with open("nearby_atoms.csv", "w") as f:
        f.write("chain,resn,resi,name,distance\n")
        for atom in atoms.atom:
            dx, dy, dz = atom.coord
            dist = ((ox - dx) ** 2 + (oy - dy) ** 2 + (oz - dz) ** 2) ** 0.5
            f.write(f"{atom.chain},{atom.resn},{atom.resi},{atom.name},{dist:.3f}\n")
    print("CSV saved as nearby_atoms.csv")

cmd.extend("save_nearby_atoms_csv", save_nearby_atoms_csv)