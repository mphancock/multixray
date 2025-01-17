import sys
sys.path.append("../../../src")
from merge_pdbs import write_merge_pdb_file



if __name__ == "__main__":
    pdb_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/259_sb_temp/16/output_0/pdbs")
    pdb_files = list(pdb_dir.glob("*.pdb"))
    pdb_files = sorted(pdb_files, key=lambda x: int(Path(x).stem))

    pdb_files = pdb_files[::10]

    out_file = Path(Path.home(), "xray/tmp/traj.pdb")
    write_merge_pdb_file(
        merge_pdb_file=out_file,
        pdb_files=pdb_files,
        occs=None,
        n=-1,
        order=True,
        state=0
    )
