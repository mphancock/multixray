from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing
import pickle
import argparse
import time
import sys

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df
import sample_bench
sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


if __name__ == "__main__":
    # job_dir = Path("/wynton/group/sali/mhancock/xray/dev/12_sulfur_only/data/out/46_synth_sa_3/9417275")
    job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/46_synth_sa_3/9417275")
    out_dirs = list(job_dir.glob("output*"))
    out_dirs = out_dirs

    max_frames = list()
    log_files = list()
    for out_dir in out_dirs:
        log_files.append(Path(out_dir, "log.csv"))

        pdb_dir = Path(out_dir, "pdbs")
        n_pdbs = len(list(pdb_dir.glob("*.pdb")))
        max_frame = (n_pdbs-1)*10
        max_frames.append(max_frame)

    field = "rmsd"
    log_lookup_df = get_stat_df.get_stat_df(
        log_file_groups=[[log_file] for log_file in log_files],
        fields=[field],
        stats=["min"],
        N=1,
        offset=10,
        equil=100,
        max_frames=max_frames
    )
    print(log_lookup_df.head())

    print(log_lookup_df["{}_min_0".format(field)].min())
    print(log_lookup_df["{}_min_0".format(field)].mean())

    for i in range(5):
        print(log_lookup_df["{}_min_0".format(field)].iloc[i])
        print(log_lookup_df["{}_min_0_log".format(field)].iloc[i])
        print(log_lookup_df["{}_min_0_id".format(field)].iloc[i])