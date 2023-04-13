from pathlib import Path
import random


def get_n_random_files(
        pdb_dir,
        N
):
    pdb_files = set(pdb_dir.glob("*"))
    rand_subset = random.sample(pdb_files, N)

    return rand_subset


if __name__ == "__main__":
    bench_dir = Path(Path.home(), "xray/benchmark")
    md_dir = Path(Path.home(), "xray/md_xray")

    decoy_files = list()
    for job_name in ["4n7f_ff_100K", "4n7f_ff_300K"]:
        for i in range(10):
            pdb_dir = Path(md_dir, job_name, "output_{}/pdbs".format(i))
            decoy_files.extend(
                get_n_random_files(
                    pdb_dir=pdb_dir,
                    N=50
                )
            )

    decoy_f = open(Path(bench_dir, "data/4n7f_decoy_single.txt"), "w")
    for decoy_file in decoy_files:
        decoy_f.write(str(decoy_file) + "\n")
    decoy_f.close()
