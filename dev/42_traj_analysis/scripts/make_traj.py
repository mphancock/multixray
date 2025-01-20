from pathlib import Path

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
from merge_pdbs import write_merge_pdb_file


if __name__ == "__main__":
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/275_visualization")

    traj_dir = Path("../data/275/trajs")
    traj_dir.mkdir(exist_ok=True)

    for i in [0,5]:
        job_dir = Path(exp_dir, str(i))
    # for job_dir in exp_dir.glob("*"):
        job_id = job_dir.stem
        pdb_dir = Path(job_dir, "output_0/pdbs")

        print(pdb_dir)

        pdb_files = list(pdb_dir.glob("*.pdb"))
        # print(pdb_files)

        out_file = Path(traj_dir, "traj_{}.pdb".format(job_id))
        write_merge_pdb_file(
            merge_pdb_file=out_file,
            pdb_files=pdb_files,
            occs=None,
            n=-1,
            order=True,
            state=0
        )
