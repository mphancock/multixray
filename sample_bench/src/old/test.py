from pathlib import Path
import sys
import os
import shutil

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/sim_anneal/src")))
sys.path.append(str(Path(Path.home(), "xray/src")))
import simulated_annealing


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")
    data_dir = Path(xray_dir, "data")
    output_dir = Path(os.getcwd(), "../test")

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    Path(output_dir).mkdir(parents=False)
    Path(output_dir, "pdbs").mkdir(parents=False)

    ref_pdb_file = Path(data_dir, "pdbs/4n7f/4n7f_heavy.pdb")
    # pdb_file = Path(data_dir, "pdbs/4n7f/4n7f_heavy.pdb")
    pdb_file=Path(Path.home(), "xray/md_xray/4n7f_100k/output_0/pdbs/md_50.pdb")

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)

    # cif_file = Path(data_dir, "reflections/4n7f/4n7f_heavy_b_20_err_10.cif")
    cif_file = Path(data_dir, "reflections/4n7f/4n7f.cif")
    sg_symbol = "P 41 2 2"
    uc_dim = (68.411, 68.411, 37.248, 90.0, 90.0, 90.0)

    steps = 100
    t_step = 2
    w_xray = 10000
    pdb_freq = 1
    n_cycles = 3

    simulated_annealing.sim_anneal(
        output_dir=output_dir,
        h=h,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        cif_file=cif_file,
        w_xray=w_xray,
        t_step=t_step,
        n_cycles=n_cycles,
        n_steps=steps,
        temp_schedule=[3000,1000,300],
        d_min_schedule=[4,2,None],
        pdb_freq=pdb_freq,
        h_0=h_0
    )