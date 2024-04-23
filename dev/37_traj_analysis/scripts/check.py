from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import os

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
from get_stat_df_simple import get_stat_df
sys.path.append(str(Path(Path.home(), "xray/src")))
# from refine import refine, pool_refine
from utility import pool_read_pdb
sys.path.append(str(Path(Path.home(), "xray/sample_bench/scripts/analysis_exp")))
from score_analysis_exp import field_id_mapper, map_std_field_to_field


if __name__ == "__main__":
    exp_name = "180_wxray_bench"
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name)
    analysis_dir = Path("/wynton/home/sali/mhancock/xray/sample_bench/data/7mhf", exp_name)
    analysis_dir.mkdir(exist_ok=True)
    n_states = [1, 2, 4, 8, 16]
    # job_ids = list(range(63))
    job_ids = list(range(20))
    w_xrays = [0.03125, 0.0625, .125, .25, .5, 1, 2, 4, 8, 16, 32]
    # cif_df = pd.read_csv(Path("/wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv"), index_col=0)
    cif_df = pd.read_csv(Path("/wynton/home/sali/mhancock/xray/dev/29_synthetic_native_3/data/cifs/csvs/7mhf_30.csv"), index_col=0)

    w_xray_df = pd.DataFrame()

    # h_rts =["08:00:00", "12:00:00", "16:00:00", "24:00:00", "56:00:00"]
    h_rt = "24:00:00"

    for state_id in range(len(n_states)):
        n_state = n_states[state_id]
        # h_rt = h_rts[state_id]
        for job_id in job_ids:
            for w_xray_id in range(len(w_xrays)):
                w_xray = w_xrays[w_xray_id]

                job_name = "N{}_J{}_W{}".format(state_id, job_id, w_xray_id)
                job_dir = Path(exp_dir, job_name)
                out_dirs = [Path(job_dir, "output_{}".format(i)) for i in range(5)]

                avg_n_pdb_files = 0
                for out_dir in out_dirs:
                    redo = False
                    if not out_dir.exists():
                        print("{} does not exist".format(out_dir))
                        redo = True
                    else:
                        pdb_files = out_dir.glob("pdbs/*")
                        n_pdb_files = len(list(pdb_files))

                        if n_pdb_files < 250:
                            print("{} does not contain enough ({})".format(out_dir, n_pdb_files))
                            redo = True

                    params = "--input_csv /wynton/home/sali/mhancock/xray/dev/35_cif_combos/data/7mhf.csv --job_id {} --w_xray {} --n_state {} --init_weights rand --sa {{step3000,T300,dofA,pdb1,w1,res0}} --steps 2".format(job_id, w_xray, n_state)
                    offset = str(out_dir).split("_")[-1]

                    if redo:
                        bash_str = "qsub -N w{} -l h_rt={} -l mem_free=1G -l scratch=1G -t 1-1 $HOME/xray/sample_bench/scripts/sample/run_slave.sh {} {} '{}' {}".format(job_name, h_rt, job_name, job_dir, params, offset)
                        # print(bash_str)
                        # os.system(bash_str)

                    avg_n_pdb_files += n_pdb_files

                if len(out_dirs) > 0:
                    avg_n_pdb_files /= len(out_dirs)

                print(state_id, job_id, w_xray_id, len(out_dirs), avg_n_pdb_files)

