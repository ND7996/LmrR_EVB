# LmrR five-EVB qmap files

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
