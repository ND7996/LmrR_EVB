$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath 'D:\PhD_Thesis\LmrR_EVB\charges\fast_lmrr_evb_qm\pdb_preserving_xtb_qm_inputs\RS2_1_to_TS2_1a'
xtb RS2_1_to_TS2_1a_TS.xyz --gfn 2 --opt tight --chrg 0 --uhf 0 --input RS2_1_to_TS2_1a_TS_xtb_constraints.inp *> RS2_1_to_TS2_1a_TS_constrained.xtb.log
if (Test-Path xtbopt.xyz) {
    Copy-Item -LiteralPath xtbopt.xyz -Destination RS2_1_to_TS2_1a_TS_constrained_xtbopt.xyz -Force
}
if (Test-Path xtbopt.log) {
    Copy-Item -LiteralPath xtbopt.log -Destination RS2_1_to_TS2_1a_TS_constrained_xtbopt.log -Force
}
if (Test-Path xtbrestart) {
    Copy-Item -LiteralPath xtbrestart -Destination RS2_1_to_TS2_1a_TS_constrained.xtbrestart -Force
}
