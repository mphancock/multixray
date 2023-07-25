from pathlib import Path
import sys
import pandas as pd
import random
random.seed(0)

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/decoys/scripts")))
import create_decoy_set_from_sample


if __name__ == "__main__":
    # sample_job_dirs = list()
    # sample_job_dirs.append(Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/53_100/9520043"))
    # sample_job_dirs.append(Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/54_1000/9520046"))

    # log_files = list()
    # for sample_job_dir in sample_job_dirs:
    #     for out_dir in sample_job_dir.glob("output_*"):
    #         log_files.append(Path(out_dir, "log.csv"))

    # random_df = create_decoy_set_from_sample.get_random_sample_df(
    #     log_files=log_files,
    #     N=10000,
    #     equil=1
    # )

    # random_df.to_csv(Path(Path.home(), "xray/dev/17_synthetic_native/data/random_df.csv"))
    random_df = pd.read_csv(Path(Path.home(), "xray/dev/17_synthetic_native/data/random_df.csv"))
    rmsd_ranges = [(0,.5),(.5,1),(1,2),(2,3)]
    n_repeats = 10

    n_state = 1
    decoy_meta_file = Path(Path.home(), "xray/dev/17_synthetic_native/data/{}x_decoy_df.csv".format(n_state))
    cols = list()
    cols.append("decoy_file")
    for i in range(n_state):
        cols.append("pdb_{}".format(i))
        cols.append("id_{}".format(i))
        cols.append("w_{}".format(i))

    n_decoys = len(rmsd_ranges)*n_repeats
    decoy_df = pd.DataFrame(index=list(range(n_decoys)), columns=cols)

    decoy_id = 0
    for min_rmsd, max_rmsd in rmsd_ranges:
        n_struct = n_state*n_repeats
        max_rmsd_df = random_df[(random_df["rmsd"] < max_rmsd) & (random_df["rmsd"] > min_rmsd)]
        rand_pdb_files = max_rmsd_df.sample(n_struct)["pdb"].values

        for i in range(n_repeats):
            decoy_file  = Path(Path.home(), "xray/dev/17_synthetic_native/data/pdbs/{}_state/{}.pdb".format(n_state, decoy_id))
            decoy_df.loc[decoy_id, "decoy_file"] = str(decoy_file)

            occs = [random.random() for j in range(n_state)]
            occs = [occ/sum(occs) for occ in occs]

            for j in range(n_state):
                decoy_df.loc[decoy_id, "pdb_{}".format(j)] = rand_pdb_files[i*n_state+j]
                decoy_df.loc[decoy_id, "id_{}".format(j)] = 0
                decoy_df.loc[decoy_id, "w_{}".format(j)] = occs[j]

            decoy_id = decoy_id+1

    decoy_df.to_csv(decoy_meta_file)








