from pathlib import Path
import sys

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    job_name = sys.argv[1]
    xray_dir = Path(Path.home(), "xray")
    decoy_dir = Path(xray_dir, "decoys/data/decoy_sets")
    bench_dir = Path(xray_dir, "score_bench")
    job_dir = Path(bench_dir, "out/{}".format(job_name))

    ref_dir = Path(xray_dir, "data/pdbs")
    cif_dir = Path(xray_dir, "data/reflections")

    params_dict = dict()
    params_dict["decoy"], params_dict["native"], params_dict["res"], params_dict["uc_dim"], params_dict["sg"], params_dict["w_xray"], params_dict["score_f"] = list(), list(), list(), list(), list(), list(), list()

    # Decoy is the directory containing the decoy files.
    params_dict["decoy"].append(Path(decoy_dir, "1ejg/1ejg_N_1000_x1"))

    # Native is a pair (native_pdb_file, cif_file).
    params_dict["native"].append((Path(ref_dir, "1ejg/1ejg_heavy.pdb"), Path(cif_dir, "1ejg/1ejg_heavy.cif")))
    params_dict["res"].append(5)

    params_dict["uc_dim"].append((40.824, 18.498, 22.371, 90.00, 90.47, 90.00))
    params_dict["sg"].append("P 1 21 1")

    # w_xray is the weight of the X-ray term in scoring.
    params_dict["w_xray"].append(10000)

    params_dict["score_f"].append("tot")

    combos = True
    params = score_rmsd.get_params(
        params_dict,
        combos
    )

    print(params)

    params_file = Path(job_dir, "params.txt")
    params_file.unlink(missing_ok=True)
    for i in range(len(params)):
        param_list = params[i]
        rmsds_file = Path(job_dir, "rmsds_{}.p".format(i))
        scores_file = Path(job_dir, "scores_{}.p".format(i))

        score_rmsd.score_vs_rmsd(
            params_file=params_file,
            decoys_dir=param_list[0],
            native=param_list[1],
            res=param_list[2],
            uc_dim=param_list[3],
            sg=param_list[4],
            w_xray=param_list[5],
            score_f=param_list[6],
            rmsds_file=rmsds_file,
            scores_file=scores_file
        )