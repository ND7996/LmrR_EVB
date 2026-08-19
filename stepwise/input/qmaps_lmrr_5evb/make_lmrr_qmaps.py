#!/usr/bin/env python3
r"""
Create LmrR Friedel-Crafts EVB qmap files for the five slide-state transitions.

The script reads the current capped QM-cluster PDB files from:
    D:\PhD_Thesis\LmrR_EVB\charges

and writes qmap files plus bond-change summaries into:
    D:\PhD_Thesis\LmrR_EVB\stepwise\input\qmaps_lmrr_5evb

Two qmap styles are generated:

1) qmaps_pdb_resname
   STATE columns use the residue names that are already in the PDB files
   (PAF, ENL, WAT, WA2, IND, ACE, NME). Use this if your .lib files keep
   those names and the FEP file handles the bond/charge changes.

2) qmaps_state_alias
   STATE columns use short state/component aliases, e.g. R11P.N2, T12E.C10.
   Use this if you create separate EVB library residues for each valence state.
   In that case, the residue block names in your .lib files must match these
   aliases, or you should edit COMPONENT_SUFFIX / state_code below.

All atoms are written as q atoms. No n boundary atoms are included because these
PDB clusters are already capped as ACE-pAF-NME and are not covalently connected
to a larger protein fragment in the qmap PDBs.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


BASE_DIR = Path(r"D:\PhD_Thesis\LmrR_EVB")
PDB_DIR = BASE_DIR / "charges"
OUT_DIR = BASE_DIR / "stepwise" / "input" / "qmaps_lmrr_5evb"


COMPONENT_SUFFIX = {
    "PAF": "P",
    "ENL": "E",
    "WAT": "W",
    "WA2": "2",
    "IND": "I",
    "ACE": "A",
    "NME": "N",
}


@dataclass(frozen=True)
class Atom:
    serial: int
    atom: str
    resname: str
    chain: str
    resid: str
    element: str

    @property
    def pdb_id(self) -> str:
        """GPX-style PDB ID column: residue-number.atom-name."""
        return f"{self.resid}.{self.atom}"

    @property
    def res_atom(self) -> str:
        return f"{self.resname}.{self.atom}"

    @property
    def comment_id(self) -> str:
        chain = self.chain if self.chain else "-"
        return f"serial {self.serial} {self.resname} {chain}{self.resid}.{self.atom}"


@dataclass(frozen=True)
class BondChange:
    kind: str
    a: int
    b: int
    note: str


@dataclass(frozen=True)
class Transition:
    name: str
    state1: str
    state1_code: str
    pdb1: str
    state2: str
    state2_code: str
    pdb2: str
    description: str
    bond_changes: tuple[BondChange, ...]


TRANSITIONS: tuple[Transition, ...] = (
    Transition(
        name="RS1_1_to_TS1_2",
        state1="RS1.1",
        state1_code="R11",
        pdb1="step_1_1_RS.pdb",
        state2="TS1.2",
        state2_code="T12",
        pdb2="step_1_1_PS.pdb",
        description=(
            "Nucleophilic addition of the pAF amine to the enal carbonyl carbon; "
            "this gives the tetrahedral ammonium/alkoxide state used as TS1.2."
        ),
        bond_changes=(
            BondChange("formed", 1, 10, "pAF amino N to enal carbonyl C; N-C bond formation"),
            BondChange("weakened", 10, 16, "enal C=O pi bond weakens as alkoxide forms"),
        ),
    ),
    Transition(
        name="TS1_2_to_PS1_2b",
        state1="TS1.2",
        state1_code="T12",
        pdb1="step_1_2_RS.pdb",
        state2="PS1.2b",
        state2_code="P12",
        pdb2="step_1_2_PS.pdb",
        description=(
            "Two-water proton redistribution from ammonium/alkoxide to neutral "
            "carbinolamine plus W1 hydronium / W2 hydroxide."
        ),
        bond_changes=(
            BondChange("broken", 1, 29, "pAF N-H proton leaves ammonium"),
            BondChange("formed", 27, 29, "W1 accepts H29 from pAF nitrogen"),
            BondChange("broken", 60, 61, "W2 loses H61"),
            BondChange("formed", 16, 61, "enal alkoxide O accepts H61 to give carbinolamine OH"),
        ),
    ),
    Transition(
        name="RS1_2b_to_PS1_3",
        state1="RS1.2b",
        state1_code="P12",
        pdb1="step_1_3_RS.pdb",
        state2="PS1.3",
        state2_code="P13",
        pdb2="step_1_3_PS.pdb",
        description=(
            "Carbinolamine dehydration to iminium; W1 hydronium protonates the OH "
            "leaving group and the C-N bond becomes iminium-like."
        ),
        bond_changes=(
            BondChange("broken", 27, 29, "W1 hydronium O-H bond breaks"),
            BondChange("formed", 16, 29, "carbinolamine O accepts H29, making leaving water"),
            BondChange("broken", 10, 16, "C-O bond to leaving water breaks"),
            BondChange("double-bond formed", 1, 10, "pAF N/enal C iminium double bond forms"),
        ),
    ),
    Transition(
        name="RS2_1_to_TS2_1a",
        state1="RS2.1",
        state1_code="R21",
        pdb1="step_2_1_RS.pdb",
        state2="TS2.1a",
        state2_code="T21",
        pdb2="step_2_1_PS.pdb",
        description=(
            "Friedel-Crafts C-C sigma-bond formation between the iminium/enal chain "
            "and indole. The step_2_1_PS structure is used as the TS2.1a valence state."
        ),
        bond_changes=(
            BondChange("formed", 12, 20, "enal beta/gamma carbon to indole C3; new C-C sigma bond"),
            BondChange("weakened", 1, 10, "iminium N=C bond weakens during conjugate addition"),
            BondChange("strengthened", 10, 11, "enal chain C10-C11 bond strengthens/enamine character changes"),
        ),
    ),
    Transition(
        name="TS2_1a_to_PS2_2",
        state1="TS2.1a",
        state1_code="T21",
        pdb1="step_2_2a_RS.pdb",
        state2="PS2.2",
        state2_code="P22",
        pdb2="step_2_2b_PS.pdb",
        description=(
            "Combined iminium/enamine tautomerization to the final alkylated product. "
            "This combines the 2.2a indolium deprotonation and 2.2b alpha-carbon protonation."
        ),
        bond_changes=(
            BondChange("broken", 20, 50, "indolium C3-H bond breaks"),
            BondChange("formed", 60, 50, "W2 accepts indole proton H50"),
            BondChange("broken", 27, 58, "W1 O-H bond breaks"),
            BondChange("formed", 11, 58, "enal alpha carbon accepts H58"),
            BondChange("double-bond formed", 1, 10, "N=C iminium/enamine bond order changes toward product state"),
        ),
    ),
)


def parse_pdb_atoms(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            atoms.append(
                Atom(
                    serial=int(line[6:11]),
                    atom=line[12:16].strip(),
                    resname=line[17:20].strip(),
                    chain=line[21:22].strip(),
                    resid=line[22:26].strip(),
                    element=line[76:78].strip(),
                )
            )
    return atoms


def atom_signature(atom: Atom) -> tuple[int, str, str, str]:
    """Identity that must stay fixed between EVB endpoints."""
    return (atom.serial, atom.atom, atom.resname, atom.resid)


def validate_atom_identity(atoms1: list[Atom], atoms2: list[Atom], transition: Transition) -> None:
    if len(atoms1) != len(atoms2):
        raise ValueError(
            f"{transition.name}: atom count mismatch: {len(atoms1)} in {transition.pdb1}, "
            f"{len(atoms2)} in {transition.pdb2}"
        )
    mismatches = []
    for a1, a2 in zip(atoms1, atoms2):
        if atom_signature(a1) != atom_signature(a2):
            mismatches.append((a1, a2))
    if mismatches:
        lines = [f"{transition.name}: atom identity/order mismatch:"]
        for a1, a2 in mismatches[:20]:
            lines.append(f"  {a1.comment_id} != {a2.comment_id}")
        if len(mismatches) > 20:
            lines.append(f"  ... and {len(mismatches) - 20} more")
        raise ValueError("\n".join(lines))


def state_alias(state_code: str, resname: str) -> str:
    suffix = COMPONENT_SUFFIX.get(resname, resname[:1].upper())
    return f"{state_code}{suffix}"


def lib_id(atom: Atom, transition: Transition, side: int, mode: str) -> str:
    if mode == "pdb":
        return atom.res_atom
    if mode == "state":
        code = transition.state1_code if side == 1 else transition.state2_code
        return f"{state_alias(code, atom.resname)}.{atom.atom}"
    raise ValueError(f"Unknown mode: {mode}")


def format_qmap(transition: Transition, mode: str) -> str:
    pdb1 = PDB_DIR / transition.pdb1
    pdb2 = PDB_DIR / transition.pdb2
    atoms1 = parse_pdb_atoms(pdb1)
    atoms2 = parse_pdb_atoms(pdb2)
    validate_atom_identity(atoms1, atoms2, transition)

    lines: list[str] = []
    lines.append(f"# {transition.name}")
    lines.append(f"# {transition.state1} ({transition.pdb1})  ->  {transition.state2} ({transition.pdb2})")
    lines.append(f"# {transition.description}")
    lines.append("#")
    lines.append("# q atom (q)         PDB ID          LIB ID (resname.name)")
    lines.append(f"# or neighbour (n)?  (resid.atom)    {transition.state1:<14} {transition.state2:<14}")
    lines.append("#")
    if mode == "state":
        lines.append("# STATE-ALIAS MODE: make sure your .lib residue names match these aliases.")
        lines.append("# Component aliases:")
        for resname in sorted({a.resname for a in atoms1}):
            lines.append(
                f"#   {resname}: {state_alias(transition.state1_code, resname)} "
                f"-> {state_alias(transition.state2_code, resname)}"
            )
        lines.append("#")
    else:
        lines.append("# PDB-RESNAME MODE: STATE columns keep the PDB residue names.")
        lines.append("# If you use separate EVB libraries for each valence state, use qmaps_state_alias instead.")
        lines.append("#")

    for atom1, atom2 in zip(atoms1, atoms2):
        lines.append(
            f"q                   {atom1.pdb_id:<15} "
            f"{lib_id(atom1, transition, 1, mode):<16} "
            f"{lib_id(atom2, transition, 2, mode):<16} "
            f"# {atom1.comment_id}"
        )

    lines.append("")
    lines.append("# Bond changes:")
    serial_map = {atom.serial: atom for atom in atoms1}
    for change in transition.bond_changes:
        a = serial_map[change.a]
        b = serial_map[change.b]
        lines.append(
            f"# {change.kind}: {change.a}-{change.b} "
            f"({a.res_atom} -- {b.res_atom}); {change.note}"
        )
    lines.append("")
    return "\n".join(lines)


def write_qmaps(modes: Iterable[str]) -> list[Path]:
    written: list[Path] = []
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for mode in modes:
        mode_dir = OUT_DIR / f"qmaps_{'pdb_resname' if mode == 'pdb' else 'state_alias'}"
        mode_dir.mkdir(parents=True, exist_ok=True)
        for transition in TRANSITIONS:
            qmap_text = format_qmap(transition, mode)
            path = mode_dir / f"{transition.name}.qmap"
            path.write_text(qmap_text, encoding="utf-8")
            written.append(path)
    return written


def write_bond_summary() -> tuple[Path, Path]:
    rows = []
    md: list[str] = []
    md.append("# LmrR five-EVB qmap bond-change summary")
    md.append("")
    md.append("These bond changes are taken from the current corrected capped PDB states.")
    md.append("")

    for transition in TRANSITIONS:
        atoms = parse_pdb_atoms(PDB_DIR / transition.pdb1)
        serial_map = {atom.serial: atom for atom in atoms}
        md.append(f"## {transition.name}")
        md.append("")
        md.append(f"{transition.state1} (`{transition.pdb1}`) -> {transition.state2} (`{transition.pdb2}`)")
        md.append("")
        md.append(transition.description)
        md.append("")
        for change in transition.bond_changes:
            a = serial_map[change.a]
            b = serial_map[change.b]
            rows.append(
                {
                    "transition": transition.name,
                    "state1": transition.state1,
                    "state2": transition.state2,
                    "change": change.kind,
                    "atom_i_serial": change.a,
                    "atom_i": a.res_atom,
                    "atom_j_serial": change.b,
                    "atom_j": b.res_atom,
                    "note": change.note,
                }
            )
            md.append(f"- {change.kind}: `{change.a}-{change.b}` = `{a.res_atom}` -- `{b.res_atom}`; {change.note}")
        md.append("")

    csv_path = OUT_DIR / "bond_changes_summary.csv"
    md_path = OUT_DIR / "bond_changes_summary.md"
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "transition",
                "state1",
                "state2",
                "change",
                "atom_i_serial",
                "atom_i",
                "atom_j_serial",
                "atom_j",
                "note",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    md_path.write_text("\n".join(md), encoding="utf-8")
    return csv_path, md_path


def write_readme() -> Path:
    readme = OUT_DIR / "README_qmaps.md"
    readme.write_text(
        """# LmrR five-EVB qmap files

Generated qmaps for the five slide-state EVB calculations:

1. `RS1_1_to_TS1_2`
2. `TS1_2_to_PS1_2b`
3. `RS1_2b_to_PS1_3`
4. `RS2_1_to_TS2_1a`
5. `TS2_1a_to_PS2_2`

Two qmap styles are provided:

- `qmaps_pdb_resname/`: STATE columns keep the current PDB residue names
  (`PAF`, `ENL`, `WAT`, `WA2`, `IND`, `ACE`, `NME`). This is easiest if your
  `.lib` files use the same residue names.
- `qmaps_state_alias/`: STATE columns use state-specific component aliases such
  as `R11P.N2` and `T12E.C10`. This is safer when every EVB valence state has
  its own library residue name. If you use these, make sure your `.lib` residue
  block names match the aliases or edit `make_lmrr_qmaps.py`.

All atoms are marked as `q`, because these are capped QM-cluster PDBs and there
are no covalent protein boundary atoms in the qmap PDBs. If you later insert
these into the full LmrR protein, boundary atoms may need to be added as `n`.

Important: qmap files map atoms to library identities; they do not by themselves
define the EVB bond/charge changes. The corresponding `.fep` files still need to
encode the bond formation/breaking and state charges.
""",
        encoding="utf-8",
    )
    return readme


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate LmrR five-EVB qmap files.")
    parser.add_argument(
        "--mode",
        choices=["all", "pdb", "state"],
        default="all",
        help="Which qmap style to write. Default: all.",
    )
    args = parser.parse_args()

    modes = ["pdb", "state"] if args.mode == "all" else [args.mode]
    written = write_qmaps(modes)
    csv_path, md_path = write_bond_summary()
    readme = write_readme()

    print(f"Wrote {len(written)} qmap files to {OUT_DIR}")
    for path in written:
        print(f"  {path}")
    print(f"Wrote {csv_path}")
    print(f"Wrote {md_path}")
    print(f"Wrote {readme}")


if __name__ == "__main__":
    main()
