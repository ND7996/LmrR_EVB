$ErrorActionPreference = 'Stop'
Set-Location -LiteralPath 'D:\PhD_Thesis\LmrR_EVB\charges\fast_lmrr_evb_qm\pdb_preserving_xtb_qm_inputs\RS1_2b_to_PS1_3'
xtb RS1_2b_to_PS1_3_RS.xyz --gfn 2 --opt tight --chrg 0 --uhf 0 *> RS1_2b_to_PS1_3_RS.xtb.log
if (Test-Path xtbopt.xyz) {
    Copy-Item -LiteralPath xtbopt.xyz -Destination RS1_2b_to_PS1_3_RS_xtbopt.xyz -Force
}
if (Test-Path xtbopt.log) {
    Copy-Item -LiteralPath xtbopt.log -Destination RS1_2b_to_PS1_3_RS_xtbopt.log -Force
}
if (Test-Path xtbrestart) {
    Copy-Item -LiteralPath xtbrestart -Destination RS1_2b_to_PS1_3_RS.xtbrestart -Force
}
