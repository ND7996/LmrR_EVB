#!/usr/bin/env python3
r"""
Fast LmrR EVB/QM preparation workflow.

This is meant to replace the slow "run everything in VeloxChem for all steps"
workflow.  By default it only prepares files; it does not launch expensive jobs.

Main ideas
----------
1. Preserve your PDB atom order and atom numbering by converting the existing
   capped PDB files to XYZ and writing a serial-to-XYZ-index map.
2. Prepare xTB jobs for fast pre-optimization of RS/TS/PS structures.
3. Optionally run xTB one step at a time if `xtb` is installed.
4. Export ORCA/Gaussian input files for final QM optimization/frequency jobs.
5. Optionally run the VeloxChem SMILES TS guesser in FAST mode only
   (no relaxed QM scan, no py3Dmol viewer).

Important
---------
- PDB serial numbers are not always the same as XYZ line indices, because the
  capped step-1 PDBs skip some serials. The generated `*_pdb_to_xyz_map.csv`
  files are therefore essential.
- The final barriers still require optimized RS/TS/PS energies and frequencies.
  This script gives fast starting points and final-QM input files.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


BASE_DIR = Path(r"D:\PhD_Thesis\LmrR_EVB")
CHARGES_DIR = BASE_DIR / "charges"
# Keep this workflow outside the `results` folder by default.
# The user wants code and generated fast-job files directly under `charges`.
RESULTS_DIR = CHARGES_DIR / "fast_lmrr_evb_qm"

HARTREE_TO_KJMOL = 2625.499638


@dataclass(frozen=True)
class PDBAtom:
    serial: int
    name: str
    resname: str
    chain: str
    resid: str
    x: float
    y: float
    z: float
    element: str

    @property
    def res_atom(self) -> str:
        return f"{self.resname}.{self.name}"

    @property
    def pdb_id(self) -> str:
        return f"{self.resid}.{self.name}"


@dataclass(frozen=True)
class BondChange:
    kind: str
    serial_i: int
    serial_j: int
    note: str


@dataclass(frozen=True)
class StructureSet:
    label: str
    rs_pdb: str
    ts_pdb: str
    ps_pdb: str
    charge: int
    multiplicity: int
    bond_changes: tuple[BondChange, ...]
    comment: str


STRUCTURE_SETS: dict[str, StructureSet] = {
    "RS1_1_to_TS1_2": StructureSet(
        label="RS1_1_to_TS1_2",
        rs_pdb="step_1_1_RS.pdb",
        ts_pdb="step_1_1_TS.pdb",
        ps_pdb="step_1_1_PS.pdb",
        charge=0,
        multiplicity=1,
        bond_changes=(
            BondChange("formed", 1, 10, "PAF.N2 attacks ENL.C10; N-C bond formation"),
            BondChange("weakened", 10, 16, "ENL.C10=O1 carbonyl weakens"),
        ),
        comment="Step 1.1 addition. Use TS as a starting TS guess; PS is the tetrahedral state.",
    ),
    "TS1_2_to_PS1_2b": StructureSet(
        label="TS1_2_to_PS1_2b",
        rs_pdb="step_1_2_RS.pdb",
        ts_pdb="step_1_2_TS.pdb",
        ps_pdb="step_1_2_PS.pdb",
        charge=0,
        multiplicity=1,
        bond_changes=(
            BondChange("broken", 1, 29, "PAF N-H breaks"),
            BondChange("formed", 27, 29, "W1 accepts H29"),
            BondChange("broken", 60, 61, "W2 O-H breaks"),
            BondChange("formed", 16, 61, "ENL.O1 accepts H61"),
        ),
        comment="Two-water proton redistribution.",
    ),
    "RS1_2b_to_PS1_3": StructureSet(
        label="RS1_2b_to_PS1_3",
        rs_pdb="step_1_3_RS.pdb",
        ts_pdb="step_1_3_TS.pdb",
        ps_pdb="step_1_3_PS.pdb",
        charge=0,
        multiplicity=1,
        bond_changes=(
            BondChange("broken", 27, 29, "W1 hydronium O-H breaks"),
            BondChange("formed", 16, 29, "carbinolamine O accepts H29"),
            BondChange("broken", 10, 16, "C-O leaving-water bond breaks"),
            BondChange("double-bond formed", 1, 10, "iminium N=C forms"),
        ),
        comment="Dehydration to iminium.",
    ),
    "RS2_1_to_TS2_1a": StructureSet(
        label="RS2_1_to_TS2_1a",
        rs_pdb="step_2_1_RS.pdb",
        ts_pdb="step_2_1_TS.pdb",
        ps_pdb="step_2_1_PS.pdb",
        charge=0,
        multiplicity=1,
        bond_changes=(
            BondChange("formed", 12, 20, "new C-C sigma bond ENL.C12--IND.C3"),
            BondChange("weakened", 1, 10, "iminium N=C weakens"),
            BondChange("strengthened", 10, 11, "C10-C11 bond strengthens"),
        ),
        comment="Friedel-Crafts C-C bond formation.",
    ),
    "TS2_1a_to_PS2_2": StructureSet(
        label="TS2_1a_to_PS2_2",
        # For the combined 2.2 tautomerization, use 2.2a RS as the early state,
        # 2.2b TS as the final-proton-transfer TS guess, and 2.2b PS as PS2.2.
        rs_pdb="step_2_2a_RS.pdb",
        ts_pdb="step_2_2b_TS.pdb",
        ps_pdb="step_2_2b_PS.pdb",
        charge=0,
        multiplicity=1,
        bond_changes=(
            BondChange("broken", 20, 50, "IND.C3-H50 breaks"),
            BondChange("formed", 60, 50, "W2 accepts H50"),
            BondChange("broken", 27, 58, "W1 O-H58 breaks"),
            BondChange("formed", 11, 58, "ENL.C11 accepts H58"),
            BondChange("double-bond formed", 1, 10, "N=C/enamine bond order changes"),
        ),
        comment="Combined 2.2 tautomerization. If you want separate 2.2a/2.2b, split this entry.",
    ),
}


STATE_SMILES = {
    "RS1.1": "CC(=O)NC(Cc1ccc(N)cc1)C(=O)NC.CCC/C=C/C=O",
    "TS1.2": "CC(=O)NC(Cc1ccc([NH2+][CH]([O-])/C=C/CCC)cc1)C(=O)NC",
    "PS1.2b": "CC(=O)NC(Cc1ccc(N[CH](O)/C=C/CCC)cc1)C(=O)NC",
    "RS1.2b": "CC(=O)NC(Cc1ccc(N[CH](O)/C=C/CCC)cc1)C(=O)NC",
    "PS1.3": "CC(=O)NC(Cc1ccc([NH+]=[CH]/C=C/CCC)cc1)C(=O)NC.O",
    "RS2.1": "CC(=O)NC(Cc1ccc([NH+]=[CH]/C=C/CCC)cc1)C(=O)NC.c1ccc2[nH]ccc2c1",
    "TS2.1a": "CC(=O)NC(Cc1ccc(N/C=C/C(CCC)C2C=[NH+]c3ccccc23)cc1)C(=O)NC",
    "PS2.2": "CC(=O)NC(Cc1ccc([NH+]=[CH]CC(CCC)c2c[nH]c3ccccc23)cc1)C(=O)NC",
    "H3O+": "[OH3+]",
    "H2O": "O",
}


SMILES_TRANSITIONS = {
    "RS1_1_to_TS1_2": (["RS1.1"], ["TS1.2"]),
    "TS1_2_to_PS1_2b": (["TS1.2"], ["PS1.2b"]),
    "RS1_2b_to_PS1_3": (["RS1.2b", "H3O+"], ["PS1.3", "H2O"]),
    "RS2_1_to_TS2_1a": (["RS2.1"], ["TS2.1a"]),
    "TS2_1a_to_PS2_2": (["TS2.1a"], ["PS2.2"]),
}


def selected_steps(step_args: Iterable[str] | None) -> list[str]:
    if not step_args:
        return list(STRUCTURE_SETS)
    out: list[str] = []
    for raw in step_args:
        for item in raw.split(","):
            item = item.strip()
            if not item:
                continue
            if item not in STRUCTURE_SETS:
                raise SystemExit(f"Unknown step {item!r}. Available: {', '.join(STRUCTURE_SETS)}")
            out.append(item)
    return out


def infer_element(atom_name: str, explicit: str = "") -> str:
    if explicit:
        return explicit.strip().capitalize()
    letters = "".join(ch for ch in atom_name if ch.isalpha())
    if not letters:
        return "X"
    two = letters[:2].capitalize()
    if two in {"Cl", "Br"}:
        return two
    return letters[0].upper()


def parse_pdb(path: Path) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atoms.append(
            PDBAtom(
                serial=int(line[6:11]),
                name=line[12:16].strip(),
                resname=line[17:20].strip(),
                chain=line[21:22].strip(),
                resid=line[22:26].strip(),
                x=float(line[30:38]),
                y=float(line[38:46]),
                z=float(line[46:54]),
                element=infer_element(line[12:16].strip(), line[76:78].strip()),
            )
        )
    if not atoms:
        raise ValueError(f"No atoms found in {path}")
    return atoms


def write_xyz_from_pdb(pdb_path: Path, xyz_path: Path, map_csv_path: Path) -> dict[int, int]:
    atoms = parse_pdb(pdb_path)
    xyz_path.parent.mkdir(parents=True, exist_ok=True)
    serial_to_xyz_index: dict[int, int] = {}
    with xyz_path.open("w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\n")
        f.write(f"from {pdb_path.name}; PDB serial map in {map_csv_path.name}\n")
        for xyz_index, atom in enumerate(atoms, start=1):
            serial_to_xyz_index[atom.serial] = xyz_index
            f.write(f"{atom.element:<2s} {atom.x:14.8f} {atom.y:14.8f} {atom.z:14.8f}\n")

    with map_csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "xyz_index",
                "pdb_serial",
                "pdb_id",
                "res_atom",
                "chain",
                "resid",
                "resname",
                "atom_name",
                "element",
            ],
        )
        writer.writeheader()
        for xyz_index, atom in enumerate(atoms, start=1):
            writer.writerow(
                {
                    "xyz_index": xyz_index,
                    "pdb_serial": atom.serial,
                    "pdb_id": atom.pdb_id,
                    "res_atom": atom.res_atom,
                    "chain": atom.chain,
                    "resid": atom.resid,
                    "resname": atom.resname,
                    "atom_name": atom.name,
                    "element": atom.element,
                }
            )
    return serial_to_xyz_index


def read_xyz(path: Path) -> tuple[list[str], list[tuple[float, float, float]], str]:
    lines = path.read_text(encoding="utf-8").splitlines()
    n = int(lines[0].strip())
    comment = lines[1] if len(lines) > 1 else ""
    labels: list[str] = []
    coords: list[tuple[float, float, float]] = []
    for line in lines[2 : 2 + n]:
        parts = line.split()
        labels.append(parts[0])
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return labels, coords, comment


def distance_from_xyz(xyz_path: Path, idx_i: int, idx_j: int) -> float:
    _, coords, _ = read_xyz(xyz_path)
    xi, yi, zi = coords[idx_i - 1]
    xj, yj, zj = coords[idx_j - 1]
    return math.sqrt((xi - xj) ** 2 + (yi - yj) ** 2 + (zi - zj) ** 2)


def write_xtb_constraint(path: Path, constraints: list[tuple[int, int, float]]) -> None:
    lines = ["$constrain"]
    for i, j, dist in constraints:
        lines.append(f"  distance: {i}, {j}, {dist:.4f}")
    lines.append("$end")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_xtb_runner(
    run_path: Path,
    xyz_path: Path,
    charge: int,
    multiplicity: int,
    constraint_file: Path | None,
    output_prefix: str,
) -> None:
    uhf = max(0, multiplicity - 1)
    cmd = [
        "xtb",
        xyz_path.name,
        "--gfn",
        "2",
        "--opt",
        "tight",
        "--chrg",
        str(charge),
        "--uhf",
        str(uhf),
    ]
    if constraint_file is not None:
        cmd.extend(["--input", constraint_file.name])
    command = " ".join(cmd)
    run_path.write_text(
        "\n".join(
            [
                "$ErrorActionPreference = 'Stop'",
                f"Set-Location -LiteralPath '{run_path.parent}'",
                command + f" *> {output_prefix}.xtb.log",
                "if (Test-Path xtbopt.xyz) {",
                f"    Copy-Item -LiteralPath xtbopt.xyz -Destination {output_prefix}_xtbopt.xyz -Force",
                "}",
                "if (Test-Path xtbopt.log) {",
                f"    Copy-Item -LiteralPath xtbopt.log -Destination {output_prefix}_xtbopt.log -Force",
                "}",
                "if (Test-Path xtbrestart) {",
                f"    Copy-Item -LiteralPath xtbrestart -Destination {output_prefix}.xtbrestart -Force",
                "}",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_orca_input(
    inp_path: Path,
    xyz_path: Path,
    charge: int,
    multiplicity: int,
    job_type: str,
    nprocs: int = 8,
    memory_mb: int = 3000,
) -> None:
    labels, coords, _ = read_xyz(xyz_path)
    if job_type == "ts":
        route = "! B3LYP D3BJ def2-SVP OptTS NumFreq TightSCF"
    elif job_type == "freq":
        route = "! B3LYP D3BJ def2-SVP Freq TightSCF"
    else:
        route = "! B3LYP D3BJ def2-SVP Opt Freq TightSCF"
    lines = [
        route,
        f"%pal nprocs {nprocs} end",
        f"%maxcore {memory_mb}",
        "",
        f"* xyz {charge} {multiplicity}",
    ]
    for label, (x, y, z) in zip(labels, coords):
        lines.append(f"  {label:<2s} {x:14.8f} {y:14.8f} {z:14.8f}")
    lines.append("*")
    lines.append("")
    inp_path.write_text("\n".join(lines), encoding="utf-8")


def write_gaussian_input(
    gjf_path: Path,
    xyz_path: Path,
    charge: int,
    multiplicity: int,
    job_type: str,
    nprocs: int = 8,
    memory_gb: int = 16,
) -> None:
    labels, coords, _ = read_xyz(xyz_path)
    opt_keyword = "Opt=(TS,CalcFC,NoEigenTest) Freq" if job_type == "ts" else "Opt Freq"
    lines = [
        f"%nprocshared={nprocs}",
        f"%mem={memory_gb}GB",
        f"#p B3LYP/def2SVP EmpiricalDispersion=GD3BJ {opt_keyword}",
        "",
        f"{xyz_path.stem} {job_type}",
        "",
        f"{charge} {multiplicity}",
    ]
    for label, (x, y, z) in zip(labels, coords):
        lines.append(f"{label:<2s} {x:14.8f} {y:14.8f} {z:14.8f}")
    lines.extend(["", ""])
    gjf_path.write_text("\n".join(lines), encoding="utf-8")


def prepare_pdb_preserving_jobs(steps: list[str], overwrite: bool = True) -> None:
    prep_dir = RESULTS_DIR / "pdb_preserving_xtb_qm_inputs"
    prep_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for step in steps:
        cfg = STRUCTURE_SETS[step]
        step_dir = prep_dir / step
        step_dir.mkdir(parents=True, exist_ok=True)
        serial_maps: dict[str, dict[int, int]] = {}

        for role, pdb_name, qm_job_type in [
            ("RS", cfg.rs_pdb, "min"),
            ("TS", cfg.ts_pdb, "ts"),
            ("PS", cfg.ps_pdb, "min"),
        ]:
            pdb_path = CHARGES_DIR / pdb_name
            if not pdb_path.exists():
                raise FileNotFoundError(pdb_path)
            xyz_path = step_dir / f"{step}_{role}.xyz"
            map_csv = step_dir / f"{step}_{role}_pdb_to_xyz_map.csv"
            if overwrite or not xyz_path.exists():
                serial_maps[role] = write_xyz_from_pdb(pdb_path, xyz_path, map_csv)
            else:
                # Read existing map.
                serial_maps[role] = {}
                with map_csv.open(newline="", encoding="utf-8") as f:
                    for row in csv.DictReader(f):
                        serial_maps[role][int(row["pdb_serial"])] = int(row["xyz_index"])

            write_orca_input(step_dir / f"{step}_{role}.inp", xyz_path, cfg.charge, cfg.multiplicity, qm_job_type)
            write_gaussian_input(step_dir / f"{step}_{role}.gjf", xyz_path, cfg.charge, cfg.multiplicity, qm_job_type)
            write_xtb_runner(
                step_dir / f"run_xtb_{role}.ps1",
                xyz_path,
                cfg.charge,
                cfg.multiplicity,
                constraint_file=None,
                output_prefix=f"{step}_{role}",
            )

            summary_rows.append(
                {
                    "step": step,
                    "role": role,
                    "pdb": str(pdb_path),
                    "xyz": str(xyz_path),
                    "map_csv": str(map_csv),
                    "orca_input": str(step_dir / f"{step}_{role}.inp"),
                    "gaussian_input": str(step_dir / f"{step}_{role}.gjf"),
                    "xtb_runner": str(step_dir / f"run_xtb_{role}.ps1"),
                    "charge": cfg.charge,
                    "multiplicity": cfg.multiplicity,
                    "comment": cfg.comment,
                }
            )

        # Write constrained xTB TS preoptimization for the TS structure.
        ts_xyz = step_dir / f"{step}_TS.xyz"
        ts_map = serial_maps["TS"]
        constraints: list[tuple[int, int, float]] = []
        for change in cfg.bond_changes:
            if change.kind in {"formed", "broken", "double-bond formed"}:
                if change.serial_i in ts_map and change.serial_j in ts_map:
                    i_xyz = ts_map[change.serial_i]
                    j_xyz = ts_map[change.serial_j]
                    dist = distance_from_xyz(ts_xyz, i_xyz, j_xyz)
                    constraints.append((i_xyz, j_xyz, dist))
        if constraints:
            constraint_file = step_dir / f"{step}_TS_xtb_constraints.inp"
            write_xtb_constraint(constraint_file, constraints)
            write_xtb_runner(
                step_dir / "run_xtb_TS_constrained.ps1",
                ts_xyz,
                cfg.charge,
                cfg.multiplicity,
                constraint_file=constraint_file,
                output_prefix=f"{step}_TS_constrained",
            )

        # Human-readable bond map for this step.
        atoms = parse_pdb(CHARGES_DIR / cfg.ts_pdb)
        by_serial = {a.serial: a for a in atoms}
        bond_lines = [
            f"# {step}",
            f"# {cfg.comment}",
            "# PDB serial -> atom identity; XYZ index is taken from TS map",
        ]
        for change in cfg.bond_changes:
            ai = by_serial[change.serial_i]
            aj = by_serial[change.serial_j]
            bond_lines.append(
                f"{change.kind:18s} PDB {change.serial_i:>3}-{change.serial_j:<3} "
                f"XYZ {ts_map[change.serial_i]:>3}-{ts_map[change.serial_j]:<3} "
                f"{ai.res_atom} -- {aj.res_atom} ; {change.note}"
            )
        (step_dir / f"{step}_bond_map.txt").write_text("\n".join(bond_lines) + "\n", encoding="utf-8")

    summary_csv = prep_dir / "prepared_fast_jobs_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "step",
                "role",
                "pdb",
                "xyz",
                "map_csv",
                "orca_input",
                "gaussian_input",
                "xtb_runner",
                "charge",
                "multiplicity",
                "comment",
            ],
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Prepared PDB-preserving xTB/QM inputs in: {prep_dir}")
    print(f"Summary: {summary_csv}")


def run_xtb_jobs(steps: list[str], roles: list[str]) -> None:
    xtb = shutil.which("xtb")
    if not xtb:
        raise SystemExit(
            "Could not find xtb on PATH. Install xTB or use the prepared run_xtb_*.ps1 files on a machine with xTB."
        )
    prep_dir = RESULTS_DIR / "pdb_preserving_xtb_qm_inputs"
    for step in steps:
        step_dir = prep_dir / step
        for role in roles:
            runner = step_dir / f"run_xtb_{role}.ps1"
            if role == "TS_constrained":
                runner = step_dir / "run_xtb_TS_constrained.ps1"
            if not runner.exists():
                print(f"Skipping missing runner: {runner}")
                continue
            print(f"Running {runner}")
            subprocess.run(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(runner)], check=True)


def _safe_float(value):
    try:
        return float(value)
    except Exception:
        return None


def _scan_items(result):
    scan = result.get("scan", {}) if isinstance(result, dict) else {}
    if not isinstance(scan, dict):
        return []
    items = []
    for lam, records in scan.items():
        if isinstance(records, dict):
            records = [records]
        for idx, rec in enumerate(records or []):
            if isinstance(rec, dict):
                items.append((lam, idx, rec))
    return items


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_veloxchem_fast(steps: list[str]) -> None:
    try:
        import veloxchem as vlx
    except Exception as exc:
        raise SystemExit(f"VeloxChem is not importable in this Python environment: {exc}") from exc

    out_dir = RESULTS_DIR / "veloxchem_fast_smiles_guesses"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_rows = []

    for step in steps:
        reactant_labels, product_labels = SMILES_TRANSITIONS[step]
        reactants = [vlx.Molecule.read_smiles(STATE_SMILES[label]) for label in reactant_labels]
        products = [vlx.Molecule.read_smiles(STATE_SMILES[label]) for label in product_labels]

        guesser = vlx.TransitionStateGuesser()
        # The key speed choices: no QM rescoring/relaxed QM scan and no viewer.
        for attr, value in [
            ("scf_scan", False),
            ("do_qm_scan", False),
            ("qm_scan", False),
            ("optimize_ts", False),
            ("do_ts_optimization", False),
            ("mute_ff_build", True),
        ]:
            if hasattr(guesser, attr):
                setattr(guesser, attr, value)

        print(f"\n=== VeloxChem FAST guess: {step} ===")
        result = guesser.find_transition_state(reactants, products)
        step_dir = out_dir / step
        step_dir.mkdir(parents=True, exist_ok=True)

        for key in ["max_mm_xyz", "max_qm_xyz"]:
            if isinstance(result, dict) and key in result:
                write_text(step_dir / f"{step}_{key}.xyz", result[key])

        scan_rows = []
        for lam, idx, rec in _scan_items(result):
            xyz = rec.get("qm_xyz") or rec.get("mm_xyz") or rec.get("xyz")
            if xyz:
                write_text(step_dir / f"lambda_{float(lam):.2f}_conf_{idx}.xyz", xyz)
            scan_rows.append(
                {
                    "step": step,
                    "lambda": float(lam),
                    "conformer_index": idx,
                    "mm_energy": rec.get("mm_energy", ""),
                    "qm_energy": rec.get("qm_energy", ""),
                    "has_xyz": bool(xyz),
                }
            )
        if scan_rows:
            write_csv(
                step_dir / f"{step}_veloxchem_fast_scan.csv",
                scan_rows,
                ["step", "lambda", "conformer_index", "mm_energy", "qm_energy", "has_xyz"],
            )

        summary_rows.append(
            {
                "step": step,
                "reactants": " + ".join(reactant_labels),
                "products": " + ".join(product_labels),
                "result_keys": ", ".join(result.keys()) if isinstance(result, dict) else "",
                "max_mm_lambda": result.get("max_mm_lambda", "") if isinstance(result, dict) else "",
                "max_qm_lambda": result.get("max_qm_lambda", "") if isinstance(result, dict) else "",
                "scan_points": len(scan_rows),
                "folder": str(step_dir),
            }
        )

    write_csv(
        out_dir / "veloxchem_fast_summary.csv",
        summary_rows,
        ["step", "reactants", "products", "result_keys", "max_mm_lambda", "max_qm_lambda", "scan_points", "folder"],
    )
    print(f"\nVeloxChem fast outputs: {out_dir}")


def write_readme() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    readme = RESULTS_DIR / "README_fast_workflow.md"
    readme.write_text(
        """# Fast LmrR EVB/QM workflow

This folder is generated by `fast_lmrr_evb_qm_workflow.py`.

Recommended fast route:

1. Prepare PDB-preserving XYZ, xTB runners, and final-QM input files:

   ```powershell
   python fast_lmrr_evb_qm_workflow.py --prepare
   ```

2. Run xTB only for one step/role at a time, for example:

   ```powershell
   python fast_lmrr_evb_qm_workflow.py --run-xtb --steps RS2_1_to_TS2_1a --roles RS,TS,PS
   ```

   Or run the generated `run_xtb_*.ps1` files manually inside each step folder.

3. Use the generated ORCA `.inp` or Gaussian `.gjf` files for final B3LYP-D/def2-SVP
   optimization/frequency calculations.

4. Use VeloxChem only for quick SMILES TS guesses:

   ```powershell
   python fast_lmrr_evb_qm_workflow.py --veloxchem-fast --steps RS2_1_to_TS2_1a
   ```

Do not run all relaxed QM scans for all five steps in one notebook. That is the slow
path that can consume days without useful intermediate stopping points.
""",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Fast LmrR EVB/QM preparation workflow.")
    parser.add_argument("--steps", nargs="*", help="Step labels, comma-separated or space-separated.")
    parser.add_argument("--prepare", action="store_true", help="Prepare PDB-preserving XYZ/xTB/QM input files.")
    parser.add_argument("--run-xtb", action="store_true", help="Run prepared xTB jobs. Requires xtb on PATH.")
    parser.add_argument(
        "--roles",
        default="RS,TS,PS",
        help="Roles for --run-xtb: RS,TS,PS,TS_constrained. Default: RS,TS,PS.",
    )
    parser.add_argument("--veloxchem-fast", action="store_true", help="Run fast VeloxChem SMILES TS guesses only.")
    parser.add_argument("--list-steps", action="store_true", help="Print available step labels and exit.")
    args = parser.parse_args()

    if args.list_steps:
        print("\n".join(STRUCTURE_SETS))
        return

    steps = selected_steps(args.steps)
    write_readme()

    if not args.prepare and not args.run_xtb and not args.veloxchem_fast:
        print("No action requested. Preparing files by default.")
        args.prepare = True

    if args.prepare:
        prepare_pdb_preserving_jobs(steps)
    if args.run_xtb:
        roles = [x.strip() for x in args.roles.split(",") if x.strip()]
        run_xtb_jobs(steps, roles)
    if args.veloxchem_fast:
        run_veloxchem_fast(steps)


if __name__ == "__main__":
    main()
