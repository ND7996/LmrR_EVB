$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath 'D:\PhD_Thesis\LmrR_EVB\charges\fast_lmrr_evb_qm\pdb_preserving_xtb_qm_inputs\TS2_1a_to_PS2_2'
xtb TS2_1a_to_PS2_2_PS.xyz --gfn 2 --opt tight --chrg 0 --uhf 0 *> TS2_1a_to_PS2_2_PS.xtb.log
if (Test-Path xtbopt.xyz) {
    Copy-Item -LiteralPath xtbopt.xyz -Destination TS2_1a_to_PS2_2_PS_xtbopt.xyz -Force
}
if (Test-Path xtbopt.log) {
    Copy-Item -LiteralPath xtbopt.log -Destination TS2_1a_to_PS2_2_PS_xtbopt.log -Force
}
if (Test-Path xtbrestart) {
    Copy-Item -LiteralPath xtbrestart -Destination TS2_1a_to_PS2_2_PS.xtbrestart -Force
}
