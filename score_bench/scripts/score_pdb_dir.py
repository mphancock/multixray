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
    parser.add_argument("--uc_dim")
    parser.add_argument("--sg_symbol")
    parser.add_argument("--ref_pdb_file")
    parser.add_argument("--param_file")
    parser.add_argument("--score_file")
    parser.add_argument("--add_native", action="store_true")
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()
    print(args.pdb_dir)
    print(args.cif_file)
    print(args.min_res)
    print(args.score_file)
    print(args.add_native)
    print(args.param_file)
    print(args.test)

    xray_dir = Path(Path.home(), "xray")
    bench_dir = Path(xray_dir, "score_bench")

    ref_dir = Path(xray_dir, "data/pdbs")

    # Decoy is the directory containing the decoy files.
    pdb_files = list(Path(args.pdb_dir).glob("*.pdb"))

    if args.test:
        pdb_files = pdb_files[:1]

    native_pdb_file = Path(args.ref_pdb_file)

    if args.add_native:
        pdb_files.append(native_pdb_file)

    native_cif_file = Path(args.cif_file)
    flags_file = native_cif_file

    uc_dim = tuple([float(args.uc_dim.split(" ")[i]) for i in range(6)])
    sg_symbol = args.sg_symbol

    score_fs = list()
    score_fs.append("ml")
    score_fs.append("rmsd")

    score_rmsd.score_vs_rmsd(
        params_file=args.param_file,
        pdb_files=pdb_files,
        native_pdb_file=native_pdb_file,
        native_cif_file=native_cif_file,
        flags_file=flags_file,
        min_res=args.min_res,
        uc_dim=uc_dim,
        sg=sg_symbol,
        score_fs=score_fs,
        scores_file=args.score_file
    )