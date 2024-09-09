from pathlib import Path
import shutil


if __name__ == "__main__":
    for i in range(1,22):
        job_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/201_lower/{}".format(i))
        dest_dir = Path(Path.home(), "xray/dev/40_higher_energy/data/201/{}".format(i))
        dest_dir.mkdir(exist_ok=True)
        for j in range(10):
            out_dir = Path(job_dir, "output_{}".format(j))
            pdb_dir = Path(out_dir, "pdbs")
            n_pdb = len([pdb_file for pdb_file in pdb_dir.glob("*.pdb")])
            pdb_file = Path(pdb_dir, "{}.pdb".format(n_pdb-1))
            dest_file = Path(dest_dir, "{}.pdb".format(j))

            shutil.copy(pdb_file, dest_file)

