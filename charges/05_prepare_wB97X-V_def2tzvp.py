
#!/usr/bin/env python3
from pathlib import Path
import csv, re, shutil

BASE_DIR = Path(r"D:\PhD_Thesis\LmrR_EVB")
CHARGES_DIR = BASE_DIR / "charges"
OUT_DIR = CHARGES_DIR / "qm_reference_for_evb"

NPROCS = 8
MAXCORE_MB = 3000
CHARGE = 0
MULTIPLICITY = 1

STRUCTURE_SETS = {
    "RS1_1_to_TS1_2": ("step_1_1_RS.pdb", "step_1_1_TS.pdb", "step_1_1_PS.pdb"),
    "TS1_2_to_PS1_2b": ("step_1_2_RS.pdb", "step_1_2_TS.pdb", "step_1_2_PS.pdb"),
    "RS1_2b_to_PS1_3": ("step_1_3_RS.pdb", "step_1_3_TS.pdb", "step_1_3_PS.pdb"),
    "RS2_1_to_TS2_1a": ("step_2_1_RS.pdb", "step_2_1_TS.pdb", "step_2_1_PS.pdb"),
    "TS2_1a_to_PS2_2": ("step_2_2a_RS.pdb", "step_2_2b_TS.pdb", "step_2_2b_PS.pdb"),
}

def mkdir(p):
    p.mkdir(parents=True, exist_ok=True)

def parse_pdb(path):
    atoms=[]
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ","HETATM")):
            continue
        element = line[76:78].strip()
        if not element:
            letters = "".join(c for c in line[12:16] if c.isalpha())
            element = ("Cl" if letters[:2].lower()=="cl" else
                       "Br" if letters[:2].lower()=="br" else
                       (letters[:1] or "X"))
        atoms.append({
            "serial": int(line[6:11]),
            "name": line[12:16].strip(),
            "resname": line[17:20].strip(),
            "chain": line[21:22].strip(),
            "resid": line[22:26].strip(),
            "element": element.capitalize(),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
        })
    if not atoms:
        raise ValueError(f"No atoms found in {path}")
    return atoms

def write_xyz(pdb, xyz, map_csv):
    atoms = parse_pdb(pdb)
    mkdir(xyz.parent)
    with xyz.open("w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n{pdb.name}\n")
        for a in atoms:
            f.write(f"{a['element']:2s} {a['x']:14.8f} {a['y']:14.8f} {a['z']:14.8f}\n")
    with map_csv.open("w", newline="", encoding="utf-8") as f:
        w=csv.writer(f)
        w.writerow(["xyz_index","pdb_serial","resid","resname","atom_name","element"])
        for i,a in enumerate(atoms,1):
            w.writerow([i,a["serial"],a["resid"],a["resname"],a["name"],a["element"]])

def find_optimized_geometry(step, state):
    candidates = [
        OUT_DIR / "B3LYP_def2SVP" / step / state / "geometry" / f"{step}_{state}_B3LYP_def2SVP_opt.xyz",
        OUT_DIR / "B3LYP_def2SVP" / step / state / "geometry" / f"{step}_{state}_DFT_opt.xyz",
        OUT_DIR / "B3LYP_def2SVP" / step / state / "geometry" / f"{step}_{state}_opt.xyz",
        CHARGES_DIR / "fast_lmrr_evb_qm" / "pdb_preserving_xtb_qm_inputs" / step / f"{step}_{state}_xtbopt.xyz",
        CHARGES_DIR / "fast_lmrr_evb_qm" / "pdb_preserving_xtb_qm_inputs" / step / f"{step}_{state}.xyz",
    ]
    for p in candidates:
        if p.exists():
            return p
    return None

def write_input(path, route, geometry):
    mkdir(path.parent)
    text = f"""! {route}

%pal
   nprocs {NPROCS}
end

%maxcore {MAXCORE_MB}

%output
   Print[P_Mulliken] 1
end

* xyzfile {CHARGE} {MULTIPLICITY} {geometry}
"""
    path.write_text(text, encoding="utf-8")

METHOD_DIR = OUT_DIR / "wB97X-V_def2TZVP"

def main():
    for label in STRUCTURE_SETS:
        for state in ("RS","TS","PS"):
            geometry = find_optimized_geometry(label, state)
            if geometry is None:
                print(f"WARNING: no optimized geometry found for {label}/{state}")
                continue
            state_dir = METHOD_DIR / label / state
            input_file = state_dir / "input" / f"{label}_{state}_wB97X-V_def2TZVP_SP.inp"
            write_input(input_file, "wB97X-V def2-TZVP SP TightSCF", geometry)
            (state_dir/"README.txt").write_text(
                "omegaB97X-V range-separated hybrid with VV10 nonlocal correlation / def2-TZVP. No D3BJ is added.\n"
                "Single-point benchmark on the optimized reference geometry.\n",
                encoding="utf-8")
    print(f"Prepared wB97X-V_def2TZVP jobs in {METHOD_DIR}")

if __name__ == "__main__":
    main()
