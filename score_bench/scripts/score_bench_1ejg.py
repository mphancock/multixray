from pathlib import Path
import sys

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    job_name = sys.argv[1]
    job_id = int(sys.argv[2])
    xray_dir = Path(Path.home(), "xray")
    decoy_dir = Path(xray_dir, "decoys/data/decoy_sets")
    bench_dir = Path(xray_dir, "score_bench")
    job_dir = Path(bench_dir, "out/1ejg/{}".format(job_name))

    ref_dir = Path(xray_dir, "data/pdbs")
    cif_dir = Path(xray_dir, "data/reflections")

    # Decoy is the directory containing the decoy files.
    decoy_dir = Path(decoy_dir, "1ejg/{}".format(sys.argv[3]))

    # Native is a pair (native_pdb_file, cif_file).
    native_pdb_file = Path(ref_dir, "1ejg/1ejg_heavy.pdb")
    native_cif_file = Path(cif_dir, "1ejg/{}".format(sys.argv[4]))

    if sys.argv[5]:
        min_res = float(sys.argv[5])
    else:
        min_res = None

    uc_dim = (40.824, 18.498, 22.371, 90.00, 90.47, 90.00)
    sg = "P 1 21 1"

    # w_xray is the weight of the X-ray term in scoring.
    w_xray = int(sys.argv[6])

    score_f = "tot"

    params_file = Path(job_dir, "params.txt")
    if job_id == 0:
        params_file.unlink(missing_ok=True)

    rmsds_file = Path(job_dir, "rmsds_{}.p".format(job_id))
    scores_file = Path(job_dir, "scores_{}.p".format(job_id))

    score_rmsd.score_vs_rmsd(
        params_file=params_file,
        decoys_dir=decoy_dir,
        native_pdb_file=native_pdb_file,
        native_cif_file=native_cif_file,
        min_res=min_res,
        uc_dim=uc_dim,
        sg=sg,
        w_xray=w_xray,
        score_f=score_f,
        rmsds_file=rmsds_file,
        scores_file=scores_file,
        add_native=True
    )