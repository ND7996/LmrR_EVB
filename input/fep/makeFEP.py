from Qpyl.qmakefep import make_fep

fepstring = make_fep(qmap_file="FC_concerted.qmap", 
                     ignore_errors=True,
                     pdb_file="/home/hp/nayanika/github/LMRR_vaccum/WT_solvated.pdb", 
                     forcefield="oplsaa",
                     parm_files=["/home/hp/nayanika/github/LmrR_EVB/parameters/qoplsaa_all.prm"],
                     lib_files=["/home/hp/nayanika/github/LmrR_EVB/parameters/qoplsaa.lib", "/home/hp/nayanika/github/LmrR_EVB/parameters/LMRR.lib"])

open("LMRR_WT4.fep", "w").write(fepstring)
