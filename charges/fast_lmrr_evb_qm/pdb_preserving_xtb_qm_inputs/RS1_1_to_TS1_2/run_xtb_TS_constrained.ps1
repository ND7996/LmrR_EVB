$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath 'D:\PhD_Thesis\LmrR_EVB\charges\fast_lmrr_evb_qm\pdb_preserving_xtb_qm_inputs\RS1_1_to_TS1_2'
xtb RS1_1_to_TS1_2_TS.xyz --gfn 2 --opt tight --chrg 0 --uhf 0 --input RS1_1_to_TS1_2_TS_xtb_constraints.inp *> RS1_1_to_TS1_2_TS_constrained.xtb.log
if (Test-Path xtbopt.xyz) {
    Copy-Item -LiteralPath xtbopt.xyz -Destination RS1_1_to_TS1_2_TS_constrained_xtbopt.xyz -Force
}
if (Test-Path xtbopt.log) {
    Copy-Item -LiteralPath xtbopt.log -Destination RS1_1_to_TS1_2_TS_constrained_xtbopt.log -Force
}
if (Test-Path xtbrestart) {
    Copy-Item -LiteralPath xtbrestart -Destination RS1_1_to_TS1_2_TS_constrained.xtbrestart -Force
}
