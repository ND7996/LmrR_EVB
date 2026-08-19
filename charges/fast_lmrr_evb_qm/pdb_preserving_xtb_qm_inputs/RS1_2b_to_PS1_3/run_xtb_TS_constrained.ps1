$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath 'D:\PhD_Thesis\LmrR_EVB\charges\fast_lmrr_evb_qm\pdb_preserving_xtb_qm_inputs\RS1_2b_to_PS1_3'
xtb RS1_2b_to_PS1_3_TS.xyz --gfn 2 --opt tight --chrg 0 --uhf 0 --input RS1_2b_to_PS1_3_TS_xtb_constraints.inp *> RS1_2b_to_PS1_3_TS_constrained.xtb.log
if (Test-Path xtbopt.xyz) {
    Copy-Item -LiteralPath xtbopt.xyz -Destination RS1_2b_to_PS1_3_TS_constrained_xtbopt.xyz -Force
}
if (Test-Path xtbopt.log) {
    Copy-Item -LiteralPath xtbopt.log -Destination RS1_2b_to_PS1_3_TS_constrained_xtbopt.log -Force
}
if (Test-Path xtbrestart) {
    Copy-Item -LiteralPath xtbrestart -Destination RS1_2b_to_PS1_3_TS_constrained.xtbrestart -Force
}
