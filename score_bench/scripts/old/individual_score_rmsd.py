from pathlib import Path
from score_bench_1ejg import score_pool_f


if __name__ == "__main__":
    job_name = "4n7f_heavy"
    job_dir = Path(Path.home(), "xray/score_rmsd_benchmark", job_name)

    xray_dir = Path(Path.home(), "xray")
    pdb_file = Path(xray_dir, "data/pdbs/4n7f_heavy.pdb")
    cif_file = Path(xray_dir, "data/reflections/4n7f_heavy.cif")

    pts_file = Path(job_dir, "plots/native_pts.txt")

    pts_f = open(pts_file, "w")

    for score_f in ["ml", "ls", "rf", "ff"]:
        param_dict = dict()
        param_dict["decoy_file"] = pdb_file
        param_dict["decoy_num"] = 0
        param_dict["ref_file"] = pdb_file
        param_dict["cif_file"] = cif_file
        param_dict["uc_dim"] = (68.411, 68.411, 37.248, 90.00, 90.00, 90.00)
        param_dict["sg"] = "P 41 2 2"
        param_dict["score_f"] = score_f

        decoy_num, score, rmsd = score_pool_f(
            params=param_dict
        )

        print(score, rmsd)
        pts_f.write("pdb_file: {}\n".format(pdb_file))
        pts_f.write("ref_file: {}\n".format(pdb_file))
        pts_f.write("cif_file: {}\n".format(cif_file))
        pts_f.write("score, rmsd: {}, {}\n".format(score, rmsd))
        pts_f.write("\n\n")

    pts_f.close()
