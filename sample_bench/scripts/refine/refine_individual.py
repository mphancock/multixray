import sys
from pathlib import Path

sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import read_pdb_and_refine_to_max_ff


if __name__ == "__main__":
    pool_params = dict()
    pool_params["pdb_file"] = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/180_wxray_bench_ref_10000/n0_j1_w5/output_2/pdbs/66.pdb")
    pool_params["out_pdb_file"] = None
    pool_params["max_ff"] = 10000
    pool_params["log_file"] = None

    read_pdb_and_refine_to_max_ff(pool_params)
