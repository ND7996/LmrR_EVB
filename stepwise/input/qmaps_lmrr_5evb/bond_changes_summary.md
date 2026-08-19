# LmrR five-EVB qmap bond-change summary

These bond changes are taken from the current corrected capped PDB states.

## RS1_1_to_TS1_2

RS1.1 (`step_1_1_RS.pdb`) -> TS1.2 (`step_1_1_PS.pdb`)

Nucleophilic addition of the pAF amine to the enal carbonyl carbon; this gives the tetrahedral ammonium/alkoxide state used as TS1.2.

- formed: `1-10` = `PAF.N2` -- `ENL.C10`; pAF amino N to enal carbonyl C; N-C bond formation
- weakened: `10-16` = `ENL.C10` -- `ENL.O1`; enal C=O pi bond weakens as alkoxide forms

## TS1_2_to_PS1_2b

TS1.2 (`step_1_2_RS.pdb`) -> PS1.2b (`step_1_2_PS.pdb`)

Two-water proton redistribution from ammonium/alkoxide to neutral carbinolamine plus W1 hydronium / W2 hydroxide.

- broken: `1-29` = `PAF.N2` -- `PAF.H22`; pAF N-H proton leaves ammonium
- formed: `27-29` = `WAT.O3` -- `PAF.H22`; W1 accepts H29 from pAF nitrogen
- broken: `60-61` = `WA2.O3` -- `WA2.H22`; W2 loses H61
- formed: `16-61` = `ENL.O1` -- `WA2.H22`; enal alkoxide O accepts H61 to give carbinolamine OH

## RS1_2b_to_PS1_3

RS1.2b (`step_1_3_RS.pdb`) -> PS1.3 (`step_1_3_PS.pdb`)

Carbinolamine dehydration to iminium; W1 hydronium protonates the OH leaving group and the C-N bond becomes iminium-like.

- broken: `27-29` = `WAT.O3` -- `PAF.H22`; W1 hydronium O-H bond breaks
- formed: `16-29` = `ENL.O1` -- `PAF.H22`; carbinolamine O accepts H29, making leaving water
- broken: `10-16` = `ENL.C10` -- `ENL.O1`; C-O bond to leaving water breaks
- double-bond formed: `1-10` = `PAF.N2` -- `ENL.C10`; pAF N/enal C iminium double bond forms

## RS2_1_to_TS2_1a

RS2.1 (`step_2_1_RS.pdb`) -> TS2.1a (`step_2_1_PS.pdb`)

Friedel-Crafts C-C sigma-bond formation between the iminium/enal chain and indole. The step_2_1_PS structure is used as the TS2.1a valence state.

- formed: `12-20` = `ENL.C12` -- `IND.C3`; enal beta/gamma carbon to indole C3; new C-C sigma bond
- weakened: `1-10` = `PAF.N2` -- `ENL.C10`; iminium N=C bond weakens during conjugate addition
- strengthened: `10-11` = `ENL.C10` -- `ENL.C11`; enal chain C10-C11 bond strengthens/enamine character changes

## TS2_1a_to_PS2_2

TS2.1a (`step_2_2a_RS.pdb`) -> PS2.2 (`step_2_2b_PS.pdb`)

Combined iminium/enamine tautomerization to the final alkylated product. This combines the 2.2a indolium deprotonation and 2.2b alpha-carbon protonation.

- broken: `20-50` = `IND.C3` -- `IND.H2`; indolium C3-H bond breaks
- formed: `60-50` = `WA2.O3` -- `IND.H2`; W2 accepts indole proton H50
- broken: `27-58` = `WAT.O3` -- `WAT.H22`; W1 O-H bond breaks
- formed: `11-58` = `ENL.C11` -- `WAT.H22`; enal alpha carbon accepts H58
- double-bond formed: `1-10` = `PAF.N2` -- `ENL.C10`; N=C iminium/enamine bond order changes toward product state
