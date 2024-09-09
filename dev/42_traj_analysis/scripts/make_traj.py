from pathlib import Path

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
from merge_pdbs import write_merge_pdb_file


if __name__ == "__main__":
    traj_dir = Path("../data/204/trajs")
    traj_dir.mkdir(exist_ok=True)

    for i in range(13):
        pdb_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/204_inflection/{}/output_0/pdbs".format(i))
        pdb_files = list(pdb_dir.glob("*.pdb"))
        print(pdb_files)

        out_file = Path(traj_dir, "traj_{}.pdb".format(i))
        write_merge_pdb_file(
            merge_pdb_file=out_file,
            pdb_files=pdb_files,
            occs=None,
            n=-1,
            order=True,
            state=0
        )
