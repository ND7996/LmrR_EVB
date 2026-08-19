$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath 'D:\PhD_Thesis\LmrR_EVB\charges\fast_lmrr_evb_qm\pdb_preserving_xtb_qm_inputs\TS1_2_to_PS1_2b'
xtb TS1_2_to_PS1_2b_TS.xyz --gfn 2 --opt tight --chrg 0 --uhf 0 *> TS1_2_to_PS1_2b_TS.xtb.log
if (Test-Path xtbopt.xyz) {
    Copy-Item -LiteralPath xtbopt.xyz -Destination TS1_2_to_PS1_2b_TS_xtbopt.xyz -Force
}
if (Test-Path xtbopt.log) {
    Copy-Item -LiteralPath xtbopt.log -Destination TS1_2_to_PS1_2b_TS_xtbopt.log -Force
}
if (Test-Path xtbrestart) {
    Copy-Item -LiteralPath xtbrestart -Destination TS1_2_to_PS1_2b_TS.xtbrestart -Force
}
