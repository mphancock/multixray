from pathlib import Path
import sys
import argparse

sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
import score_rmsd


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir")
    parser.add_argument("--cif_file")
    parser.add_argument("--min_res", type=float)
    parser.add_argument("--w_xray", type=float)
    parser.add_argument("--score_file")
    parser.add_argument("--include_native", type=int)
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()
    print(args.pdb_dir)
    print(args.cif_file)
    print(args.min_res)
    print(args.w_xray)
    print(args.score_file)
    print(args.include_native)
    print(args.test)

    xray_dir = Path(Path.home(), "xray")
    bench_dir = Path(xray_dir, "score_bench")

    ref_dir = Path(xray_dir, "data/pdbs")

    # Decoy is the directory containing the decoy files.
    pdb_files = list(Path(args.pdb_dir).glob("*.pdb"))

    native_pdb_file = Path(xray_dir, "data/pdbs/3ca7/3ca7_clean.pdb")
    native_cif_file = Path(args.cif_file)
    flags_file = native_cif_file

    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg = "C 1 2 1"

    score_fs = list()
    score_fs.append("r_all")
    score_fs.append("r_work")
    score_fs.append("r_free")
    score_fs.append("total")
    score_fs.append("rmsd")

    # params_file = Path(job_dir, "params.txt")
    # if job_id == 0:
    #     params_file.unlink(missing_ok=True)

    score_rmsd.score_vs_rmsd(
        params_file=None,
        pdb_files=pdb_files,
        native_pdb_file=native_pdb_file,
        native_cif_file=native_cif_file,
        flags_file=flags_file,
        min_res=args.min_res,
        uc_dim=uc_dim,
        sg=sg,
        w_xray=args.w_xray,
        score_fs=score_fs,
        scores_file=args.score_file,
        add_native=args.include_native,
        test=False
    )