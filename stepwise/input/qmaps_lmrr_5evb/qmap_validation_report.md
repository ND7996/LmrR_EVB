# LmrR qmap validation report

Generated: 2026-08-14T14:20:34

Validation rules:

- Each qmap PDB ID (`resid.atom`, GPX/Q style) must match exactly one atom in endpoint 1 and endpoint 2 PDBs.
- Atom serial, residue name, residue number, and atom name must remain identical between endpoint PDBs.
- Each bond-change serial number must resolve to the intended atom in both endpoint PDBs.

## RS1_1_to_TS1_2

RS1.1 `step_1_1_RS.pdb` -> TS1.2 `step_1_1_PS.pdb`

- Atom counts: 57 -> 57
- Endpoint atom identity/order mismatches: 0
- qmaps_pdb_resname: 57 q lines checked, 0 PDB-ID errors
- qmaps_state_alias: 57 q lines checked, 0 PDB-ID errors

| Change | Serial pair | Endpoint 1 atoms | Endpoint 2 atoms | Status |
|---|---:|---|---|---|
| formed | `1-10` | `PAF.N2 -- ENL.C10` | `PAF.N2 -- ENL.C10` | OK |
| weakened | `10-16` | `ENL.C10 -- ENL.O1` | `ENL.C10 -- ENL.O1` | OK |

## TS1_2_to_PS1_2b

TS1.2 `step_1_2_RS.pdb` -> PS1.2b `step_1_2_PS.pdb`

- Atom counts: 57 -> 57
- Endpoint atom identity/order mismatches: 0
- qmaps_pdb_resname: 57 q lines checked, 0 PDB-ID errors
- qmaps_state_alias: 57 q lines checked, 0 PDB-ID errors

| Change | Serial pair | Endpoint 1 atoms | Endpoint 2 atoms | Status |
|---|---:|---|---|---|
| broken | `1-29` | `PAF.N2 -- PAF.H22` | `PAF.N2 -- PAF.H22` | OK |
| formed | `27-29` | `WAT.O3 -- PAF.H22` | `WAT.O3 -- PAF.H22` | OK |
| broken | `60-61` | `WA2.O3 -- WA2.H22` | `WA2.O3 -- WA2.H22` | OK |
| formed | `16-61` | `ENL.O1 -- WA2.H22` | `ENL.O1 -- WA2.H22` | OK |

## RS1_2b_to_PS1_3

RS1.2b `step_1_3_RS.pdb` -> PS1.3 `step_1_3_PS.pdb`

- Atom counts: 57 -> 57
- Endpoint atom identity/order mismatches: 0
- qmaps_pdb_resname: 57 q lines checked, 0 PDB-ID errors
- qmaps_state_alias: 57 q lines checked, 0 PDB-ID errors

| Change | Serial pair | Endpoint 1 atoms | Endpoint 2 atoms | Status |
|---|---:|---|---|---|
| broken | `27-29` | `WAT.O3 -- PAF.H22` | `WAT.O3 -- PAF.H22` | OK |
| formed | `16-29` | `ENL.O1 -- PAF.H22` | `ENL.O1 -- PAF.H22` | OK |
| broken | `10-16` | `ENL.C10 -- ENL.O1` | `ENL.C10 -- ENL.O1` | OK |
| double-bond formed | `1-10` | `PAF.N2 -- ENL.C10` | `PAF.N2 -- ENL.C10` | OK |

## RS2_1_to_TS2_1a

RS2.1 `step_2_1_RS.pdb` -> TS2.1a `step_2_1_PS.pdb`

- Atom counts: 76 -> 76
- Endpoint atom identity/order mismatches: 0
- qmaps_pdb_resname: 76 q lines checked, 0 PDB-ID errors
- qmaps_state_alias: 76 q lines checked, 0 PDB-ID errors

| Change | Serial pair | Endpoint 1 atoms | Endpoint 2 atoms | Status |
|---|---:|---|---|---|
| formed | `12-20` | `ENL.C12 -- IND.C3` | `ENL.C12 -- IND.C3` | OK |
| weakened | `1-10` | `PAF.N2 -- ENL.C10` | `PAF.N2 -- ENL.C10` | OK |
| strengthened | `10-11` | `ENL.C10 -- ENL.C11` | `ENL.C10 -- ENL.C11` | OK |

## TS2_1a_to_PS2_2

TS2.1a `step_2_2a_RS.pdb` -> PS2.2 `step_2_2b_PS.pdb`

- Atom counts: 76 -> 76
- Endpoint atom identity/order mismatches: 0
- qmaps_pdb_resname: 76 q lines checked, 0 PDB-ID errors
- qmaps_state_alias: 76 q lines checked, 0 PDB-ID errors

| Change | Serial pair | Endpoint 1 atoms | Endpoint 2 atoms | Status |
|---|---:|---|---|---|
| broken | `20-50` | `IND.C3 -- IND.H2` | `IND.C3 -- IND.H2` | OK |
| formed | `60-50` | `WA2.O3 -- IND.H2` | `WA2.O3 -- IND.H2` | OK |
| broken | `27-58` | `WAT.O3 -- WAT.H22` | `WAT.O3 -- WAT.H22` | OK |
| formed | `11-58` | `ENL.C11 -- WAT.H22` | `ENL.C11 -- WAT.H22` | OK |
| double-bond formed | `1-10` | `PAF.N2 -- ENL.C10` | `PAF.N2 -- ENL.C10` | OK |

## Final result

- Total validation errors: 0
- PASS: all qmap PDB IDs and bond-change atom serials match the endpoint PDBs.