from pathlib import Path
import pandas as pd
import sys
import multiprocessing

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import refine_hs



if __name__ == "__main__":
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/169_N8/55/output_349/pdbs/253.pdb")), m, IMP.atom.AllPDBSelector())

    n_step = 20
    refine_hs(hs=hs, n_step=20, log_file=Path(Path.home(), "xray/tmp/log.csv"))

    IMP.atom.write_multimodel_pdb(hs, str(Path(Path.home(), "xray/tmp/253_ref_{}.pdb".format(n_step))))

