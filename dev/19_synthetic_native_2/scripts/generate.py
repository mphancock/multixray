from pathlib import Path
import sys
import pandas as pd
import random
random.seed(0)

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/decoys/scripts")))
import create_decoy_set_from_sample
sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


if __name__ == "__main__":
    # sample_job_dirs = list()
    # sample_job_dirs.append(Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/77_sn_realistic/817051/"))

    # log_files = list()
    # for sample_job_dir in sample_job_dirs:
    #     for out_dir in sample_job_dir.glob("output_*"):
    #         log_files.append(Path(out_dir, "log.csv"))

    # random_df = create_decoy_set_from_sample.get_random_sample_df(
    #     log_files=log_files,
    #     N=10000,
    #     equil=1
    # )

    # random_df.to_csv(Path(Path.home(), "xray/dev/19_synthetic_native_realistic/data/random_df.csv"))
    random_df = pd.read_csv(Path(Path.home(), "xray/dev/19_synthetic_native_realistic/data/random_df.csv"))
    n_repeats = 50

    n_state = 4
    decoy_meta_file = Path(Path.home(), "xray/dev/19_synthetic_native_realistic/data/decoys/{}x_decoy_df.csv".format(n_state))
    cols = list()
    cols.append("decoy_file")
    for i in range(n_state):
        cols.append("pdb_{}".format(i))
        cols.append("id_{}".format(i))
        cols.append("w_{}".format(i))

    decoy_df = pd.DataFrame(index=list(range(n_repeats)), columns=cols)
    rand_pdb_files = random_df.sample(n_state*n_repeats)["pdb"].values

    for i in range(n_repeats):
        decoy_file  = Path(Path.home(), "xray/dev/19_synthetic_native_realistic/data/pdbs/{}_state/{}.pdb".format(n_state, i))
        decoy_df.loc[i, "decoy_file"] = str(decoy_file)

        occs = [random.random() for j in range(n_state)]
        occs = [occ/sum(occs) for occ in occs]

        for j in range(n_state):
            decoy_df.loc[i, "pdb_{}".format(j)] = rand_pdb_files[i*n_state+j]
            decoy_df.loc[i, "id_{}".format(j)] = 0
            decoy_df.loc[i, "w_{}".format(j)] = occs[j]

    decoy_df.to_csv(decoy_meta_file)









