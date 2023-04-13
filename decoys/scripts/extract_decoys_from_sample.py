from pathlib import Path
import random
import pandas as pd
random.seed(0)
import IMP
import IMP.atom
import sys
import multiprocessing

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import extract_pdbs
import merge_pdbs


def create_decoy_dataset(
        decoy_pdb_dir,
        pdb_dirs,
        n_decoys,
        n_struct
):
    pdb_sample = extract_pdbs.get_n_random_files_from_pdb_dirs(
        pdb_dirs=pdb_dirs,
        n=n_decoys*n_struct,
        equil=100
    )

    decoy_df = pd.DataFrame(index=list(range(n_decoys)), columns=range(n_struct))

    pool_params = list()
    for i in range(n_decoys):
        out_file = Path(decoy_pdb_dir, "{}.pdb".format(i))
        decoy_pdb_files = pdb_sample[i*n_struct:(i+1)*n_struct]
        for j in range(n_struct):
            decoy_df.loc[i,j]=decoy_pdb_files[j]

        params_dict = dict()
        params_dict["merge_pdb_file"] = out_file
        params_dict["pdb_files"] = decoy_pdb_files
        params_dict["occs"] = [1/n_struct]*n_struct
        params_dict["n"] = -1
        params_dict["order"] = False
        params_dict["id"] = i

        pool_params.append(params_dict)

        # merge_pdbs.write_merge_pdb_file(
        #     merge_pdb_file=out_file,
        #     pdb_files=decoy_pdb_files,
        #     occs=[1/n_struct]*n_struct,
        #     n=-1,
        #     order=False
        # )
    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        merge_pdbs.write_merge_pdf_file_pool,
        pool_params
    )

    for id in pool_results:
        print(id)

    return decoy_df


if __name__ == "__main__":
    job_lookup = list()
    job_lookup.append(("24", "2406616"))

    for job_pair in job_lookup:
        print(job_pair)
        job_name, job_num = job_pair

        n_struct = 2
        n_decoys = 1000
        decoy_name = "3ca7_N_1000_x2"

        decoy_pdb_dir = Path("/wynton/group/sali/mhancock/xray/decoys/data/{}/{}".format(job_name, decoy_name))
        decoy_pdb_dir.mkdir(exist_ok=True, parents=True)

        decoy_meta_dir = Path(Path.home(), "xray/decoys/data/{}".format(job_name))
        decoy_meta_dir.mkdir(exist_ok=True, parents=True)
        decoy_meta_file = Path(decoy_meta_dir, "{}.csv".format(decoy_name))

        # job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/{}/{}".format(job_name, job_num))
        job_dir = Path(Path.home(), "xray/sample_bench/out/old/md_out/3ca7_300")

        # Get and save the dataframe containing all the decoy entries (each decoy entry containing 1 or more random pdb files and an equal number of weights).
        # pdb_dirs = [Path(out_dir, "pdbs") for out_dir in job_dir.glob("output_*")]
        pdb_dirs = list()
        for out_dir in job_dir.glob("output_*"):
            if Path(out_dir, "pdbs").exists():
                pdb_dirs.append(Path(out_dir, "pdbs"))

        # pdb_dirs = pdb_dirs[:10]
        decoy_df = create_decoy_dataset(
                decoy_pdb_dir=decoy_pdb_dir,
                pdb_dirs=pdb_dirs,
                n_decoys=n_decoys,
                n_struct=n_struct
        )
        decoy_df.to_csv(decoy_meta_file)

