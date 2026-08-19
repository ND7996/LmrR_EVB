$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath 'D:\PhD_Thesis\LmrR_EVB\charges\fast_lmrr_evb_qm\pdb_preserving_xtb_qm_inputs\RS1_1_to_TS1_2'
xtb RS1_1_to_TS1_2_PS.xyz --gfn 2 --opt tight --chrg 0 --uhf 0 *> RS1_1_to_TS1_2_PS.xtb.log
if (Test-Path xtbopt.xyz) {
    Copy-Item -LiteralPath xtbopt.xyz -Destination RS1_1_to_TS1_2_PS_xtbopt.xyz -Force
}
if (Test-Path xtbopt.log) {
    Copy-Item -LiteralPath xtbopt.log -Destination RS1_1_to_TS1_2_PS_xtbopt.log -Force
}
if (Test-Path xtbrestart) {
    Copy-Item -LiteralPath xtbrestart -Destination RS1_1_to_TS1_2_PS.xtbrestart -Force
}
