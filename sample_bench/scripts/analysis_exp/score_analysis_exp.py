from pathlib import Path
import sys
import pandas as pd
import numpy as np
import argparse

sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import get_stat_df_simple


"""
This function takes the job id which is between 0 and 62 and returns the cif file id for that job corresponding to the field id. Recall that the fields in a log_df are just indexed 0, ..., n_cond - 1 and therefore the field ids will correspond to different cif files. Also recall that we need to adjust the job_id because the binary strings are job_id + 1.
"""
def field_id_mapper(
        job_id,
        field_id
):
    mapper_df = pd.read_csv(Path(Path.home(), "xray/dev/35_cif_combos/data/7mhf.csv"))

    try:
        cif_file = Path(mapper_df.loc[job_id, "cifs"].split(",")[field_id])
        return cif_file.stem
    except IndexError:
        return None


def map_std_field_to_field(
        std_field,
        job_id
):
    mapper_df = pd.read_csv(Path(Path.home(), "xray/dev/35_cif_combos/data/7mhf.csv"))
    cif_files = [Path(cif_file) for cif_file in mapper_df.iloc[job_id]["cifs"].split(",")]

    field_id = -1
    for i in range(len(cif_files)):
        if cif_files[i].stem == std_field.split("_")[-1]:
            field_id = i
            break

    if std_field[0] == "x":
        field_name = "xray"
    elif std_field[0] == "r":
        field_name = "r_free"
    elif std_field[0] == "w":
        field_name = "w_{}".format(std_field.split("_")[1])

    if field_id >= 0:
        return "{}_{}".format(field_name, field_id)
    else:
        return None



def get_joint_field(
        n_cond
):
    field = ""
    for j in range(n_cond):
        field = field + "xray_{}".format(j)
        if j < n_cond-1:
            field = field + "+"

    return field

def get_n_cond(
        i
):
    mapper_df = pd.read_csv(Path(Path.home(), "xray/dev/35_cif_combos/data/7mhf.csv"))
    n_cond = len(mapper_df.iloc[i]["cifs"].split(","))
    return n_cond


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exp_name")
    args = parser.parse_args()

    # params = [("165_J1_i", 1), ("164_J2_ik", 2), ("162_J3_ijk", 3), ("163_J4_fijk", 4)]
    # n_cond = args.n_cond
    # job_name = args.job_name

    exp_name = args.exp_name
    n_state = int(exp_name.split("_")[-1][1:])

    analysis_dir = Path(Path.home(), "xray/sample_bench/data/7mhf", exp_name)
    analysis_dir.mkdir(exist_ok=True)

    stat_df_names = ["joint"]
    cif_names = ["7mhf", "7mhg", "7mhh", "7mhi", "7mhj", "7mhk"]
    for cif_name in cif_names:
        stat_df_names.append("r_free_{}".format(cif_name))

    # for std_field in ["joint"]:
    for std_field in stat_df_names:
        all_std_fields = list()
        if std_field == "joint":
            all_std_fields.append("joint")

        for cif_name in cif_names:
            all_std_fields.append("r_free_{}".format(cif_name))
            all_std_fields.append("xray_{}".format(cif_name))

        all_std_fields.extend(["pdb", "ff"])

        # Add all possible weights for all possible conditions
        for j in range(6):
            for i in range(n_state):
                all_std_fields.append("w_{}_{}".format(i, cif_names[j]))

        all_std_bonus_fields = all_std_fields.copy()
        all_std_bonus_fields.remove(std_field)

        score_file = Path(analysis_dir, "{}.csv".format(std_field))
        score_df = pd.DataFrame(index=range(63), columns=all_std_fields)

        for job_id in range(63):
        # for job_id in [57]:
            print(job_id, std_field)

            n_cond = get_n_cond(job_id)
            std_field_to_field_dict = dict()

            std_bonus_fields = all_std_bonus_fields.copy()

            for std_field_tmp in all_std_fields:
                if std_field_tmp == "joint":
                    std_field_to_field_dict[std_field_tmp] = get_joint_field(n_cond)
                elif std_field_tmp in ["pdb", "ff"]:
                    std_field_to_field_dict[std_field_tmp] = std_field_tmp
                else:
                    tmp_field = map_std_field_to_field(std_field_tmp, job_id)
                    std_field_to_field_dict[std_field_tmp] = tmp_field

                    # Remove invalid fields from std_bonus_fields
                    if not tmp_field and std_field_tmp in std_bonus_fields:
                        std_bonus_fields.remove(std_field_tmp)

        #     break
        # break

            # This occurs if the std_field has no corresponding field in the log_dfs. For example job_id 1 only has r_free_7mhg and no r_free_7mhf field.
            if not std_field_to_field_dict[std_field]:
                continue

            print("dict: ", std_field_to_field_dict)
            print("std_field: ", std_field)
            print("std_bonus_fields: ", std_bonus_fields)

            valid_std_fields = [std_field]
            valid_std_fields.extend(std_bonus_fields)

            bonus_fields = [std_field_to_field_dict[std_field_tmp] for std_field_tmp in std_bonus_fields]
            field = std_field_to_field_dict[std_field]

            job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf", exp_name, str(job_id))
            log_files = [Path(out_dir, "log.csv") for out_dir in job_dir.glob("output_*")]

            stat_df = get_stat_df_simple.get_stat_df(
                log_files=log_files,
                field=field,
                N=1,
                bonus_fields=bonus_fields,
                equil=1,
                pdb_only=True
            )

            for std_field_tmp in valid_std_fields:
                score_df.loc[job_id, std_field_tmp] = stat_df.iloc[0][std_field_to_field_dict[std_field_tmp]]

        score_df.to_csv(score_file)
