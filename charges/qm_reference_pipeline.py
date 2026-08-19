#!/usr/bin/env python3
"""
LmrR QM reference -> EVB preparation pipeline.

Workflow:

    PDB
      |
      v
    XYZ + atom mapping
      |
      v
    DFT input generation
      |
      v
    DFT calculation
      |
      +--> energy
      +--> optimized geometry
      +--> Mulliken charges
      +--> frequencies
      |
      v
    QM reference dataset
      |
      v
    QM DeltaE barrier
      |
      v
    EVB preparation / later comparison

IMPORTANT:
    QM DeltaE^‡ is NOT the same physical quantity as EVB DeltaG^‡.

The QM stage provides the electronic/reference-state information.
The later EVB stage provides the ensemble free energy in protein/solvent.

Expected source layout:

    D:/PhD_Thesis/LmrR_EVB/
        charges/
            step_1_1_RS.pdb
            step_1_1_TS.pdb
            step_1_1_PS.pdb
            step_1_2_RS.pdb
            step_1_2_TS.pdb
            step_1_2_PS.pdb
            ...

Output:

    charges/
        qm_reference_for_evb/
            step_1_1/
                RS/
                    input/
                    geometry/
                    energy/
                    charges/
                    frequencies/
                TS/
                    ...
                PS/
                    ...

            step_1_2/
                ...

            summaries/
                qm_state_summary.csv
                qm_barriers_for_evb_comparison.csv
                qm_mulliken_charges_long.csv
                qm_frequencies_long.csv
                qm_mapping_validation.csv

            evb_reference/
                qm_vs_evb_comparison_template.csv

Commands:

    PREPARE:
        python qm_reference_pipeline.py --prepare

    COLLECT:
        python qm_reference_pipeline.py --collect

    VALIDATE:
        python qm_reference_pipeline.py --validate

    ALL:
        python qm_reference_pipeline.py --all

If no argument is supplied, ONLY --prepare is run.
This prevents the script from attempting to collect QM output files
before the QM calculations have been completed.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


# =============================================================================
# USER CONFIGURATION
# =============================================================================

BASE_DIR = Path(r"D:\PhD_Thesis\LmrR_EVB")

CHARGES_DIR = BASE_DIR / "charges"

OUT_DIR = CHARGES_DIR / "qm_reference_for_evb"

HARTREE_TO_KJMOL = 2625.499638

TEMPERATURE_K = 298.15

# ORCA settings
NPROCS = 8
MAXCORE_MB = 3000

# DFT reference method
DFT_FUNCTIONAL = "PBE0"
DFT_BASIS = "def2-SVP"

# Higher-level benchmark
DFT_BENCHMARK_BASIS = "def2-TZVP"

# =============================================================================
# DATA CLASSES
# =============================================================================


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


# =============================================================================
# REACTION / STRUCTURE DEFINITIONS
# =============================================================================

STRUCTURE_SETS: dict[str, StructureSet] = {

    "step_1_1": StructureSet(
        "step_1_1",
        "step_1_1_RS.pdb",
        "step_1_1_TS.pdb",
        "step_1_1_PS.pdb",
        0,
        1,

        (
            BondChange(
                "formed",
                1,
                10,
                "PAF.N2 attacks ENL.C10; N-C formation"
            ),

            BondChange(
                "weakened",
                10,
                16,
                "ENL.C10=O1 weakens"
            ),
        ),

        "Step 1.1 addition; TS is the TS starting structure.",
    ),

    "step_1_2": StructureSet(
        "step_1_2",
        "step_1_2_RS.pdb",
        "step_1_2_TS.pdb",
        "step_1_2_PS.pdb",
        0,
        1,

        (
            BondChange(
                "broken",
                1,
                29,
                "PAF N-H breaks"
            ),

            BondChange(
                "formed",
                27,
                29,
                "W1 accepts H29"
            ),

            BondChange(
                "broken",
                60,
                61,
                "W2 O-H breaks"
            ),

            BondChange(
                "formed",
                16,
                61,
                "ENL.O1 accepts H61"
            ),
        ),

        "Two-water proton redistribution.",
    ),

    "step_1_3": StructureSet(
        "step_1_3",
        "step_1_3_RS.pdb",
        "step_1_3_TS.pdb",
        "step_1_3_PS.pdb",
        0,
        1,

        (
            BondChange(
                "broken",
                27,
                29,
                "W1 hydronium O-H breaks"
            ),

            BondChange(
                "formed",
                16,
                29,
                "carbinolamine O accepts H29"
            ),

            BondChange(
                "broken",
                10,
                16,
                "C-O leaving-water bond breaks"
            ),

            BondChange(
                "double-bond formed",
                1,
                10,
                "iminium N=C forms"
            ),
        ),

        "Dehydration to iminium.",
    ),

    "step_2_1": StructureSet(
        "step_2_1",
        "step_2_1_RS.pdb",
        "step_2_1_TS.pdb",
        "step_2_1_PS.pdb",
        0,
        1,

        (
            BondChange(
                "formed",
                12,
                20,
                "new C-C bond ENL.C12--IND.C3"
            ),

            BondChange(
                "weakened",
                1,
                10,
                "iminium N=C weakens"
            ),

            BondChange(
                "strengthened",
                10,
                11,
                "C10-C11 strengthens"
            ),
        ),

        "Friedel-Crafts C-C bond formation.",
    ),

    "step_2_2": StructureSet(
        "step_2_2",
        "step_2_2a_RS.pdb",
        "step_2_2b_TS.pdb",
        "step_2_2b_PS.pdb",
        0,
        1,

        (
            BondChange(
                "broken",
                20,
                50,
                "IND.C3-H50 breaks"
            ),

            BondChange(
                "formed",
                60,
                50,
                "W2 accepts H50"
            ),

            BondChange(
                "broken",
                27,
                58,
                "W1 O-H58 breaks"
            ),

            BondChange(
                "formed",
                11,
                58,
                "ENL.C11 accepts H58"
            ),
        ),

        "Combined 2.2 tautomerization; split into 2.2a/2.2b if required later.",
    ),
}


ROLES = ("RS", "TS", "PS")


# =============================================================================
# GENERAL HELPERS
# =============================================================================


def selected_steps(values: Iterable[str] | None) -> list[str]:
    """
    Select reaction steps.

    Example:

        --steps step_1_1 step_1_2

    or:

        --steps step_1_1,step_1_2
    """

    if not values:
        return list(STRUCTURE_SETS)

    out = []

    for value in values:

        for item in value.split(","):

            item = item.strip()

            if not item:
                continue

            if item not in STRUCTURE_SETS:
                raise SystemExit(
                    f"Unknown step '{item}'. "
                    f"Available: {', '.join(STRUCTURE_SETS)}"
                )

            out.append(item)

    # Remove duplicates while preserving order
    return list(dict.fromkeys(out))


def ensure_directory(path: Path) -> None:
    """Create directory if it does not exist."""
    path.mkdir(parents=True, exist_ok=True)


def write_csv(
    path: Path,
    rows: list[dict],
    fieldnames: list[str]
) -> None:
    """
    Write dictionaries to CSV.

    This function fixes the previous:
        NameError: name 'write_csv' is not defined
    """

    ensure_directory(path.parent)

    with path.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames,
            extrasaction="ignore"
        )

        writer.writeheader()

        for row in rows:
            writer.writerow(row)


# =============================================================================
# PDB READING
# =============================================================================


def pdb_atoms(path: Path) -> list[dict]:

    atoms = []

    text = path.read_text(
        encoding="utf-8",
        errors="replace"
    )

    for line in text.splitlines():

        if not line.startswith(("ATOM  ", "HETATM")):
            continue

        element = (
            line[76:78].strip()
            or re.sub(
                r"[^A-Za-z]",
                "",
                line[12:16]
            ).strip()[:1]
        )

        try:

            atom = {
                "serial": int(line[6:11]),
                "atom_name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resid": line[22:26].strip(),
                "element": element.capitalize(),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            }

        except ValueError as exc:

            raise ValueError(
                f"Could not parse PDB atom line:\n{line}\n"
                f"File: {path}"
            ) from exc

        atoms.append(atom)

    if not atoms:
        raise ValueError(
            f"No ATOM/HETATM records found in:\n{path}"
        )

    return atoms


# =============================================================================
# PDB -> XYZ
# =============================================================================


def write_xyz(
    pdb: Path,
    xyz: Path,
    mapping: Path
) -> None:

    atoms = pdb_atoms(pdb)

    ensure_directory(xyz.parent)

    # ---------------------------------------------------------
    # XYZ
    # ---------------------------------------------------------

    with xyz.open(
        "w",
        encoding="utf-8"
    ) as f:

        f.write(
            f"{len(atoms)}\n"
            f"{pdb.name}\n"
        )

        for atom in atoms:

            f.write(
                f"{atom['element']:2s} "
                f"{atom['x']:14.8f} "
                f"{atom['y']:14.8f} "
                f"{atom['z']:14.8f}\n"
            )

    # ---------------------------------------------------------
    # PDB -> XYZ mapping
    # ---------------------------------------------------------

    ensure_directory(mapping.parent)

    fieldnames = list(atoms[0].keys()) + ["xyz_index"]

    with mapping.open(
        "w",
        newline="",
        encoding="utf-8"
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=fieldnames
        )

        writer.writeheader()

        for i, atom in enumerate(atoms, 1):

            row = dict(atom)

            row["xyz_index"] = i

            writer.writerow(row)


# =============================================================================
# XYZ READING
# =============================================================================


def read_xyz(path: Path):

    lines = path.read_text(
        encoding="utf-8"
    ).splitlines()

    if len(lines) < 2:
        raise ValueError(
            f"Invalid XYZ file: {path}"
        )

    try:
        n = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(
            f"First line of XYZ is not an atom count: {path}"
        ) from exc

    data = []

    coordinate_lines = lines[2:2 + n]

    if len(coordinate_lines) != n:
        raise ValueError(
            f"XYZ atom count mismatch in {path}: "
            f"header says {n}, "
            f"found {len(coordinate_lines)}"
        )

    for line in coordinate_lines:

        parts = line.split()

        if len(parts) < 4:
            continue

        data.append(
            (
                parts[0],
                float(parts[1]),
                float(parts[2]),
                float(parts[3]),
            )
        )

    return data


# =============================================================================
# DIRECTORY CREATION
# =============================================================================


def create_state_directories(
    step: str,
    role: str
) -> Path:

    role_dir = OUT_DIR / step / role

    # IMPORTANT:
    # All directories are created before ANY files are written.
    for subdir in (
        "input",
        "geometry",
        "energy",
        "charges",
        "frequencies",
    ):

        ensure_directory(
            role_dir / subdir
        )

    return role_dir


# =============================================================================
# DFT INPUT GENERATION
# =============================================================================


def write_dft_inputs(
    step: str,
    role: str,
    xyz: Path,
    role_dir: Path,
    charge: int,
    mult: int
) -> None:

    # ---------------------------------------------------------
    # FIX FOR YOUR FileNotFoundError
    # ---------------------------------------------------------

    input_dir = role_dir / "input"

    geometry_dir = role_dir / "geometry"

    energy_dir = role_dir / "energy"

    charges_dir = role_dir / "charges"

    frequencies_dir = role_dir / "frequencies"

    for directory in (
        input_dir,
        geometry_dir,
        energy_dir,
        charges_dir,
        frequencies_dir,
    ):

        ensure_directory(directory)

    # ---------------------------------------------------------
    # Read XYZ
    # ---------------------------------------------------------

    data = read_xyz(xyz)

    coordinates = []

    for element, x, y, z in data:

        coordinates.append(
            f"  {element:<2} "
            f"{x:14.8f} "
            f"{y:14.8f} "
            f"{z:14.8f}"
        )

    # ---------------------------------------------------------
    # RS / PS
    #
    # Optimize minimum + frequency
    #
    # TS
    #
    # Optimize transition state + frequency
    # ---------------------------------------------------------

    if role == "TS":

        optimization = "OptTS NumFreq"

    else:

        optimization = "Opt NumFreq"

    # ---------------------------------------------------------
    # ORCA input
    # ---------------------------------------------------------

    orca_input = [

        f"! {DFT_FUNCTIONAL} "
        f"D3BJ "
        f"{DFT_BASIS} "
        f"{optimization} "
        f"TightSCF",

        "",

        "%pal",

        f"  nprocs {NPROCS}",

        "end",

        "",

        f"%maxcore {MAXCORE_MB}",

        "",

        f"* xyz {charge} {mult}",

        *coordinates,

        "*",

        "",
    ]

    dft_input = (
        input_dir /
        f"{step}_{role}_DFT.inp"
    )

    dft_input.write_text(
        "\n".join(orca_input),
        encoding="utf-8"
    )

    print(
        f"  Created DFT input:\n"
        f"    {dft_input}"
    )

    # ---------------------------------------------------------
    # TZVP benchmark
    #
    # This uses the optimized geometry generated after the
    # main DFT calculation.
    # ---------------------------------------------------------

    optimized_xyz_name = (
        f"{step}_{role}_DFT_opt.xyz"
    )

    benchmark_input = [

        f"! {DFT_FUNCTIONAL} "
        f"D3BJ "
        f"{DFT_BENCHMARK_BASIS} "
        f"SP "
        f"TightSCF",

        "",

        "%pal",

        f"  nprocs {NPROCS}",

        "end",

        "",

        f"%maxcore {MAXCORE_MB}",

        "",

        (
            f"* xyzfile "
            f"{charge} "
            f"{mult} "
            f"../geometry/{optimized_xyz_name}"
        ),

        "",
    ]

    benchmark_file = (
        input_dir /
        f"{step}_{role}_DFT_"
        f"{DFT_BENCHMARK_BASIS}_SP.inp"
    )

    benchmark_file.write_text(
        "\n".join(benchmark_input),
        encoding="utf-8"
    )

    print(
        f"  Created TZVP benchmark input:\n"
        f"    {benchmark_file}"
    )


# =============================================================================
# PREPARE
# =============================================================================


def prepare(
    steps: list[str]
) -> None:

    ensure_directory(OUT_DIR)

    print()
    print("=" * 72)
    print("PREPARING QM REFERENCE DATASET")
    print("=" * 72)

    print()
    print(f"Source:")
    print(f"    {CHARGES_DIR}")

    print()
    print(f"Output:")
    print(f"    {OUT_DIR}")

    for step in steps:

        cfg = STRUCTURE_SETS[step]

        print()
        print("-" * 72)
        print(f"STEP: {step}")
        print("-" * 72)

        step_dir = OUT_DIR / step

        ensure_directory(step_dir)

        for role, pdb_name in (
            ("RS", cfg.rs_pdb),
            ("TS", cfg.ts_pdb),
            ("PS", cfg.ps_pdb),
        ):

            print()
            print(f"[{step}] {role}")

            # -------------------------------------------------
            # Create complete directory tree FIRST
            # -------------------------------------------------

            role_dir = create_state_directories(
                step,
                role
            )

            print(
                f"  Directory: {role_dir}"
            )

            # -------------------------------------------------
            # Source PDB
            # -------------------------------------------------

            src = CHARGES_DIR / pdb_name

            if not src.exists():

                print(
                    "  WARNING: source PDB does not exist:"
                )

                print(
                    f"    {src}"
                )

                print(
                    "  Skipping this state."
                )

                continue

            print(
                f"  PDB:\n"
                f"    {src}"
            )

            # -------------------------------------------------
            # Initial XYZ
            # -------------------------------------------------

            xyz = (
                role_dir /
                "geometry" /
                f"{step}_{role}_initial.xyz"
            )

            mapping = (
                role_dir /
                "geometry" /
                f"{step}_{role}_pdb_to_xyz.csv"
            )

            write_xyz(
                src,
                xyz,
                mapping
            )

            print(
                f"  XYZ:\n"
                f"    {xyz}"
            )

            print(
                f"  Mapping:\n"
                f"    {mapping}"
            )

            # -------------------------------------------------
            # DFT inputs
            # -------------------------------------------------

            write_dft_inputs(
                step=step,
                role=role,
                xyz=xyz,
                role_dir=role_dir,
                charge=cfg.charge,
                mult=cfg.multiplicity,
            )

            # -------------------------------------------------
            # State README
            # -------------------------------------------------

            readme = (
                role_dir /
                "README.md"
            )

            readme.write_text(
                f"""# {step} / {role}

## Starting structure

`{pdb_name}`

## QM input

Main DFT:

`input/{step}_{role}_DFT.inp`

Higher-level benchmark:

`input/{step}_{role}_DFT_{DFT_BENCHMARK_BASIS}_SP.inp`

## Expected outputs

After running the main DFT job, the pipeline expects:

- QM energy
- optimized geometry
- Mulliken charges
- vibrational frequencies

The optimized geometry should preferably be stored as:

`geometry/{step}_{role}_DFT_opt.xyz`

The DFT output file can be placed in:

`input/`

or elsewhere underneath this state directory.

The collector searches recursively for:

- `.out`
- `.log`
- `.txt`

## EVB relevance

This state is one of the QM reference states used before EVB parameterization.

The QM energy is an electronic/reference quantity.

It should not be interpreted directly as the final EVB activation free energy.
""",
                encoding="utf-8"
            )

        # -----------------------------------------------------
        # Metadata
        # -----------------------------------------------------

        write_step_metadata(step)

    # ---------------------------------------------------------
    # Master README
    # ---------------------------------------------------------

    write_master_readme()

    print()
    print("=" * 72)
    print("QM PREPARATION COMPLETE")
    print("=" * 72)

    print()
    print(
        f"QM reference directory:\n"
        f"    {OUT_DIR}"
    )


# =============================================================================
# STEP METADATA
# =============================================================================


def write_step_metadata(
    step: str
) -> None:

    cfg = STRUCTURE_SETS[step]

    step_dir = OUT_DIR / step

    metadata = {

        "step": step,

        "comment": cfg.comment,

        "charge": cfg.charge,

        "multiplicity": cfg.multiplicity,

        "states": {

            "RS": cfg.rs_pdb,

            "TS": cfg.ts_pdb,

            "PS": cfg.ps_pdb,
        },

        "bond_changes": [

            {
                "kind": b.kind,
                "serial_i": b.serial_i,
                "serial_j": b.serial_j,
                "note": b.note,
            }

            for b in cfg.bond_changes
        ],

        "qm_reference": {

            "functional": DFT_FUNCTIONAL,

            "basis": DFT_BASIS,

            "benchmark_basis": DFT_BENCHMARK_BASIS,

            "temperature_K": TEMPERATURE_K,
        },
    }

    output = (
        step_dir /
        "step_metadata.json"
    )

    output.write_text(
        json.dumps(
            metadata,
            indent=2
        ),
        encoding="utf-8"
    )


# =============================================================================
# MASTER README
# =============================================================================


def write_master_readme() -> None:

    ensure_directory(OUT_DIR)

    text = []

    text.append(
        "# LmrR QM reference dataset for EVB\n"
    )

    text.append(
        "This directory contains the QM reference layer "
        "used before EVB parameterization.\n"
    )

    text.append(
        "## Workflow\n\n"
        "PDB -> XYZ -> DFT -> energy/geometry/charges/frequencies "
        "-> QM reference -> EVB calibration -> EVB/MD\n"
    )

    text.append(
        "## State organization\n\n"
        "Each reaction step contains:\n\n"
        "- RS\n"
        "- TS\n"
        "- PS\n\n"
        "Each state contains:\n\n"
        "- input/\n"
        "- geometry/\n"
        "- energy/\n"
        "- charges/\n"
        "- frequencies/\n"
    )

    text.append(
        "## QM reference quantities\n\n"
        "- electronic energy\n"
        "- optimized geometry\n"
        "- Mulliken charges\n"
        "- frequencies\n"
        "- QM activation energy\n"
        "- QM reaction energy\n"
    )

    text.append(
        "## Important interpretation\n\n"
        "QM DeltaE^‡ is an electronic/reference-model quantity. "
        "EVB DeltaG^‡ is an ensemble free-energy quantity "
        "obtained after sampling the protein/solvent environment. "
        "They should not be treated as identical observables.\n"
    )

    text.append(
        "## Commands\n\n"
        "```powershell\n"
        "python qm_reference_pipeline.py --prepare\n"
        "python qm_reference_pipeline.py --collect\n"
        "python qm_reference_pipeline.py --validate\n"
        "```\n"
    )

    output = (
        OUT_DIR /
        "README.md"
    )

    output.write_text(
        "\n".join(text),
        encoding="utf-8"
    )


# =============================================================================
# ENERGY EXTRACTION
# =============================================================================


def extract_energy(
    text: str
):

    patterns = [

        # ORCA
        r"FINAL SINGLE POINT ENERGY\s+"
        r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",

        # Generic
        r"TOTAL ENERGY\s*=\s*"
        r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",

        r"total energy\s*[:=]\s*"
        r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)",

        r"E\s*=\s*"
        r"(-?\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?)"
        r"\s*Eh",
    ]

    found = []

    for pattern in patterns:

        matches = re.findall(
            pattern,
            text,
            flags=re.IGNORECASE
        )

        found.extend(matches)

    if not found:
        return None

    return float(found[-1])


# =============================================================================
# MULLIKEN CHARGES
# =============================================================================


def extract_orca_mulliken(
    text: str
):

    pattern = (
        r"MULLIKEN ATOMIC CHARGES"
        r".*?"
        r"(?:\n\s*Sum of atomic charges"
        r"|\n\s*\*{3,}"
        r"|\Z)"
    )

    match = re.search(
        pattern,
        text,
        flags=re.IGNORECASE | re.DOTALL
    )

    if not match:
        return []

    charges = []

    for line in match.group(0).splitlines():

        # Typical ORCA format:
        #
        # 0 C :  0.123
        #
        # or
        #
        # 1 C  0.123

        numbers = re.findall(
            r"[-+]?\d+(?:\.\d+)?",
            line
        )

        if not numbers:
            continue

        # Avoid capturing atom number only
        if len(numbers) >= 2:

            try:

                charges.append(
                    float(numbers[-1])
                )

            except ValueError:
                pass

    return charges


# =============================================================================
# xTB CHARGES
# =============================================================================


def extract_xtb_mulliken(
    text: str
):

    match = re.search(
        r"Mulliken charges:\s*\n"
        r"(.*?)(?:\n\s*\n|\Z)",
        text,
        flags=re.IGNORECASE | re.DOTALL
    )

    if not match:
        return []

    values = []

    for line in match.group(1).splitlines():

        parts = line.split()

        if len(parts) < 2:
            continue

        try:

            values.append(
                float(parts[-1])
            )

        except ValueError:
            continue

    return values


# =============================================================================
# FREQUENCY EXTRACTION
# =============================================================================


def extract_frequencies(
    text: str
):

    values = []

    # ORCA normal-mode blocks often contain
    # frequency values in cm^-1.

    for line in text.splitlines():

        if (
            "cm" not in line.lower()
            and "freq" not in line.lower()
        ):
            continue

        numbers = re.findall(
            r"[-+]?\d+(?:\.\d+)?",
            line
        )

        if not numbers:
            continue

        # Avoid lines where the only number is an index
        for value in numbers:

            try:

                values.append(
                    float(value)
                )

            except ValueError:
                pass

    return values


# =============================================================================
# FILE SEARCH
# =============================================================================


def find_outputs(
    role_dir: Path,
    step: str,
    role: str
):

    candidates = []

    for extension in (
        "*.out",
        "*.log",
        "*.txt",
    ):

        candidates.extend(
            role_dir.rglob(extension)
        )

    results = []

    for path in candidates:

        if not path.is_file():
            continue

        if path.name.lower() == "readme.md":
            continue

        # Avoid parsing generated CSV/JSON/etc.
        results.append(path)

    return sorted(
        set(results)
    )


# =============================================================================
# COPY HELPER
# =============================================================================


def copy_if_exists(
    src: Path,
    dst: Path
) -> bool:

    if not src.exists():
        return False

    ensure_directory(
        dst.parent
    )

    shutil.copy2(
        src,
        dst
    )

    return True


# =============================================================================
# COLLECT QM RESULTS
# =============================================================================


def collect(
    steps: list[str]
) -> None:

    summary_rows = []

    charge_rows = []

    frequency_rows = []

    print()
    print("=" * 72)
    print("COLLECTING QM RESULTS")
    print("=" * 72)

    for step in steps:

        for role in ROLES:

            role_dir = (
                OUT_DIR /
                step /
                role
            )

            if not role_dir.exists():

                print(
                    f"WARNING: missing directory "
                    f"{role_dir}"
                )

                continue

            print()
            print(
                f"[{step}] {role}"
            )

            outputs = find_outputs(
                role_dir,
                step,
                role
            )

            energies = []

            charges = []

            frequencies = []

            # -------------------------------------------------
            # Parse outputs
            # -------------------------------------------------

            for output in outputs:

                try:

                    text = output.read_text(
                        encoding="utf-8",
                        errors="replace"
                    )

                except Exception as exc:

                    print(
                        f"  Could not read {output}: "
                        f"{exc}"
                    )

                    continue

                # Energy
                energy = extract_energy(text)

                if energy is not None:

                    energies.append(
                        (
                            output,
                            energy
                        )
                    )

                # Charges
                q = (
                    extract_orca_mulliken(text)
                    or extract_xtb_mulliken(text)
                )

                if q:

                    charges.append(
                        (
                            output,
                            q
                        )
                    )

                # Frequencies
                freq = extract_frequencies(text)

                if freq:

                    frequencies.append(
                        (
                            output,
                            freq
                        )
                    )

            # -------------------------------------------------
            # Energy
            # -------------------------------------------------

            if energies:

                source, energy = energies[-1]

                energy_file = (
                    role_dir /
                    "energy" /
                    f"{step}_{role}_energy.json"
                )

                energy_file.write_text(
                    json.dumps(
                        {
                            "hartree": energy,

                            "energy_kJ_mol_absolute":
                                energy *
                                HARTREE_TO_KJMOL,

                            "source":
                                str(source),
                        },
                        indent=2
                    ),
                    encoding="utf-8"
                )

                print(
                    f"  Energy: "
                    f"{energy:.12f} Eh"
                )

            else:

                source = None

                energy = None

                print(
                    "  WARNING: no QM energy found"
                )

            # -------------------------------------------------
            # Charges
            # -------------------------------------------------

            if charges:

                charge_source, charge_values = (
                    charges[-1]
                )

                charge_file = (
                    role_dir /
                    "charges" /
                    f"{step}_{role}_mulliken.csv"
                )

                charge_rows_for_state = []

                for index, value in enumerate(
                    charge_values,
                    start=1
                ):

                    charge_rows_for_state.append(
                        {
                            "step": step,
                            "state": role,
                            "atom_index": index,
                            "mulliken_charge": value,
                            "source": str(
                                charge_source
                            ),
                        }
                    )

                write_csv(
                    charge_file,
                    charge_rows_for_state,
                    [
                        "step",
                        "state",
                        "atom_index",
                        "mulliken_charge",
                        "source",
                    ]
                )

                charge_rows.extend(
                    charge_rows_for_state
                )

                print(
                    f"  Mulliken charges: "
                    f"{len(charge_values)} atoms"
                )

            else:

                print(
                    "  WARNING: no Mulliken charges found"
                )

            # -------------------------------------------------
            # Frequencies
            # -------------------------------------------------

            if frequencies:

                freq_source, freq_values = (
                    frequencies[-1]
                )

                frequency_rows_for_state = []

                for index, value in enumerate(
                    freq_values,
                    start=1
                ):

                    frequency_rows_for_state.append(
                        {
                            "step": step,
                            "state": role,
                            "frequency_index": index,
                            "frequency_cm-1": value,
                            "source": str(
                                freq_source
                            ),
                        }
                    )

                frequency_file = (
                    role_dir /
                    "frequencies" /
                    f"{step}_{role}_frequencies.csv"
                )

                write_csv(
                    frequency_file,
                    frequency_rows_for_state,
                    [
                        "step",
                        "state",
                        "frequency_index",
                        "frequency_cm-1",
                        "source",
                    ]
                )

                frequency_rows.extend(
                    frequency_rows_for_state
                )

                print(
                    f"  Frequencies: "
                    f"{len(freq_values)} values"
                )

            else:

                print(
                    "  WARNING: no frequencies found"
                )

            # -------------------------------------------------
            # Geometry
            # -------------------------------------------------

            geometry_dir = (
                role_dir /
                "geometry"
            )

            geometry_candidates = [

                geometry_dir /
                f"{step}_{role}_DFT_opt.xyz",

                geometry_dir /
                f"{step}_{role}_opt.xyz",

                role_dir /
                f"{step}_{role}_DFT_opt.xyz",

                role_dir /
                "xtbopt.xyz",
            ]

            geometry_found = False

            for geometry in geometry_candidates:

                if geometry.exists():

                    final_geometry = (
                        geometry_dir /
                        f"{step}_{role}_final.xyz"
                    )

                    copy_if_exists(
                        geometry,
                        final_geometry
                    )

                    print(
                        f"  Geometry:\n"
                        f"    {final_geometry}"
                    )

                    geometry_found = True

                    break

            if not geometry_found:

                print(
                    "  WARNING: optimized geometry "
                    "not found"
                )

            # -------------------------------------------------
            # Summary
            # -------------------------------------------------

            summary_rows.append(
                {
                    "step": step,

                    "state": role,

                    "energy_hartree":
                        energy
                        if energy is not None
                        else "",

                    "energy_kJ_mol_absolute":
                        (
                            energy *
                            HARTREE_TO_KJMOL
                        )
                        if energy is not None
                        else "",

                    "energy_source":
                        str(source)
                        if source is not None
                        else "",

                    "mulliken_atoms":
                        len(
                            charges[-1][1]
                        )
                        if charges
                        else 0,

                    "frequency_values":
                        len(
                            frequencies[-1][1]
                        )
                        if frequencies
                        else 0,

                    "optimized_geometry_found":
                        geometry_found,
                }
            )

    # =========================================================
    # WRITE SUMMARY FILES
    # =========================================================

    summaries = (
        OUT_DIR /
        "summaries"
    )

    ensure_directory(
        summaries
    )

    # ---------------------------------------------------------
    # State summary
    # ---------------------------------------------------------

    write_csv(
        summaries /
        "qm_state_summary.csv",

        summary_rows,

        [
            "step",
            "state",
            "energy_hartree",
            "energy_kJ_mol_absolute",
            "energy_source",
            "mulliken_atoms",
            "frequency_values",
            "optimized_geometry_found",
        ]
    )

    # ---------------------------------------------------------
    # Long charge table
    # ---------------------------------------------------------

    write_csv(
        summaries /
        "qm_mulliken_charges_long.csv",

        charge_rows,

        [
            "step",
            "state",
            "atom_index",
            "mulliken_charge",
            "source",
        ]
    )

    # ---------------------------------------------------------
    # Long frequency table
    # ---------------------------------------------------------

    write_csv(
        summaries /
        "qm_frequencies_long.csv",

        frequency_rows,

        [
            "step",
            "state",
            "frequency_index",
            "frequency_cm-1",
            "source",
        ]
    )

    # ---------------------------------------------------------
    # Barrier table
    # ---------------------------------------------------------

    make_barrier_summary(
        steps
    )

    print()
    print("=" * 72)
    print("QM COLLECTION COMPLETE")
    print("=" * 72)


# =============================================================================
# QM BARRIER SUMMARY
# =============================================================================


def make_barrier_summary(
    steps: list[str]
) -> None:

    rows = []

    for step in steps:

        values = {}

        for role in ROLES:

            energy_file = (
                OUT_DIR /
                step /
                role /
                "energy" /
                f"{step}_{role}_energy.json"
            )

            if not energy_file.exists():
                continue

            try:

                data = json.loads(
                    energy_file.read_text(
                        encoding="utf-8"
                    )
                )

                values[role] = data["hartree"]

            except Exception:

                continue

        # -----------------------------------------------------
        # QM activation energy
        #
        # Delta E‡ = E_TS - E_RS
        # -----------------------------------------------------

        if (
            "RS" in values
            and
            "TS" in values
        ):

            delta_e_activation = (
                values["TS"]
                -
                values["RS"]
            ) * HARTREE_TO_KJMOL

        else:

            delta_e_activation = ""

        # -----------------------------------------------------
        # QM reaction energy
        #
        # Delta E_rxn = E_PS - E_RS
        # -----------------------------------------------------

        if (
            "RS" in values
            and
            "PS" in values
        ):

            delta_e_reaction = (
                values["PS"]
                -
                values["RS"]
            ) * HARTREE_TO_KJMOL

        else:

            delta_e_reaction = ""

        rows.append(
            {
                "step": step,

                "RS_Eh":
                    values.get(
                        "RS",
                        ""
                    ),

                "TS_Eh":
                    values.get(
                        "TS",
                        ""
                    ),

                "PS_Eh":
                    values.get(
                        "PS",
                        ""
                    ),

                "QM_DeltaE_activation_kJ_mol":
                    delta_e_activation,

                "QM_DeltaE_reaction_kJ_mol":
                    delta_e_reaction,
            }
        )

    output = (
        OUT_DIR /
        "summaries" /
        "qm_barriers_for_evb_comparison.csv"
    )

    write_csv(
        output,
        rows,
        [
            "step",
            "RS_Eh",
            "TS_Eh",
            "PS_Eh",
            "QM_DeltaE_activation_kJ_mol",
            "QM_DeltaE_reaction_kJ_mol",
        ]
    )


# =============================================================================
# VALIDATION
# =============================================================================


def validate(
    steps: list[str]
) -> None:

    print()
    print("=" * 72)
    print("VALIDATING RS / TS / PS MAPPING")
    print("=" * 72)

    rows = []

    for step in steps:

        paths = {}

        for role in ROLES:

            geometry_dir = (
                OUT_DIR /
                step /
                role /
                "geometry"
            )

            final_geometry = (
                geometry_dir /
                f"{step}_{role}_final.xyz"
            )

            initial_geometry = (
                geometry_dir /
                f"{step}_{role}_initial.xyz"
            )

            if final_geometry.exists():

                paths[role] = final_geometry

            elif initial_geometry.exists():

                paths[role] = initial_geometry

        counts = {}

        for role, path in paths.items():

            try:

                counts[role] = len(
                    read_xyz(path)
                )

            except Exception:

                counts[role] = ""

        valid_counts = [
            value
            for value in counts.values()
            if isinstance(value, int)
        ]

        same_atom_count = (
            len(valid_counts) == 3
            and
            len(set(valid_counts)) == 1
        )

        rows.append(
            {
                "step": step,

                "RS_atoms":
                    counts.get(
                        "RS",
                        ""
                    ),

                "TS_atoms":
                    counts.get(
                        "TS",
                        ""
                    ),

                "PS_atoms":
                    counts.get(
                        "PS",
                        ""
                    ),

                "same_atom_count":
                    same_atom_count,
            }
        )

        print()
        print(
            f"{step}:"
        )

        print(
            f"  RS atoms = "
            f"{counts.get('RS', 'missing')}"
        )

        print(
            f"  TS atoms = "
            f"{counts.get('TS', 'missing')}"
        )

        print(
            f"  PS atoms = "
            f"{counts.get('PS', 'missing')}"
        )

        print(
            f"  Same atom count = "
            f"{same_atom_count}"
        )

    output = (
        OUT_DIR /
        "summaries" /
        "qm_mapping_validation.csv"
    )

    write_csv(
        output,
        rows,
        [
            "step",
            "RS_atoms",
            "TS_atoms",
            "PS_atoms",
            "same_atom_count",
        ]
    )


# =============================================================================
# EVB COMPARISON TEMPLATE
# =============================================================================


def make_evb_comparison_template(
    steps: list[str]
) -> None:

    output = (
        OUT_DIR /
        "evb_reference" /
        "qm_vs_evb_comparison_template.csv"
    )

    ensure_directory(
        output.parent
    )

    rows = []

    for step in steps:

        rows.append(
            {
                "step": step,

                "QM_DeltaE_activation_kJ_mol":
                    "",

                "QM_DeltaE_reaction_kJ_mol":
                    "",

                "QM_DeltaG_activation_kJ_mol":
                    "",

                "EVB_reference_DeltaE_activation_kJ_mol":
                    "",

                "EVB_DeltaG_activation_kJ_mol":
                    "",

                "EVB_DeltaG_reaction_kJ_mol":
                    "",

                "EVB_minus_QM_activation_kJ_mol":
                    "",

                "notes":
                    "",
            }
        )

    write_csv(
        output,
        rows,
        [
            "step",
            "QM_DeltaE_activation_kJ_mol",
            "QM_DeltaE_reaction_kJ_mol",
            "QM_DeltaG_activation_kJ_mol",
            "EVB_reference_DeltaE_activation_kJ_mol",
            "EVB_DeltaG_activation_kJ_mol",
            "EVB_DeltaG_reaction_kJ_mol",
            "EVB_minus_QM_activation_kJ_mol",
            "notes",
        ]
    )


# =============================================================================
# MAIN
# =============================================================================


def main():

    parser = argparse.ArgumentParser(
        description=(
            "LmrR QM reference -> EVB preparation pipeline"
        )
    )

    parser.add_argument(
        "--steps",
        nargs="*",
        help=(
            "Steps to process. "
            "Example: step_1_1 step_1_2"
        )
    )

    parser.add_argument(
        "--prepare",
        action="store_true",
        help=(
            "Create directories, PDB->XYZ files "
            "and DFT inputs"
        )
    )

    parser.add_argument(
        "--collect",
        action="store_true",
        help=(
            "Collect completed QM outputs"
        )
    )

    parser.add_argument(
        "--validate",
        action="store_true",
        help=(
            "Validate RS/TS/PS atom counts"
        )
    )

    parser.add_argument(
        "--all",
        action="store_true",
        help=(
            "Run prepare + collect + validate"
        )
    )

    args = parser.parse_args()

    steps = selected_steps(
        args.steps
    )

    # =========================================================
    # IMPORTANT
    #
    # If no command is specified:
    #
    #     python qm_reference_pipeline.py
    #
    # ONLY PREPARE.
    #
    # This avoids attempting to collect output files that
    # do not exist yet.
    # =========================================================

    if not (
        args.prepare
        or
        args.collect
        or
        args.validate
        or
        args.all
    ):

        args.prepare = True

    # =========================================================
    # --ALL
    # =========================================================

    if args.all:

        args.prepare = True
        args.collect = True
        args.validate = True

    # =========================================================
    # PREPARE
    # =========================================================

    if args.prepare:

        prepare(
            steps
        )

    # =========================================================
    # COLLECT
    # =========================================================

    if args.collect:

        collect(
            steps
        )

    # =========================================================
    # VALIDATE
    # =========================================================

    if args.validate:

        validate(
            steps
        )

    # =========================================================
    # EVB COMPARISON TEMPLATE
    # =========================================================

    make_evb_comparison_template(
        steps
    )

    print()
    print("=" * 72)
    print("DONE")
    print("=" * 72)

    print()
    print(
        "QM reference dataset:"
    )

    print(
        f"    {OUT_DIR}"
    )

    print()
    print(
        "Important summary files:"
    )

    print(
        f"    {OUT_DIR / 'summaries' / 'qm_state_summary.csv'}"
    )

    print(
        f"    {OUT_DIR / 'summaries' / 'qm_barriers_for_evb_comparison.csv'}"
    )

    print(
        f"    {OUT_DIR / 'summaries' / 'qm_mulliken_charges_long.csv'}"
    )

    print(
        f"    {OUT_DIR / 'summaries' / 'qm_mapping_validation.csv'}"
    )

    print()
    print(
        "EVB comparison template:"
    )

    print(
        f"    {OUT_DIR / 'evb_reference' / 'qm_vs_evb_comparison_template.csv'}"
    )


# =============================================================================
# ENTRY POINT
# =============================================================================


if __name__ == "__main__":

    main()