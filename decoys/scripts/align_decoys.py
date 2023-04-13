from pathlib import Path
import sys

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import align_imp


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")
    job_dir = Path(xray_dir, "decoys/md_out/3ca7_300")

    out_dirs = job_dir.glob("output*")

    ref_pdb_file = Path(xray_dir, "data/pdbs/3ca7/3ca7_clean.pdb")
    m_ref = IMP.Model()
    h_ref = IMP.atom.read_pdb(str(ref_pdb_file), m_ref, IMP.atom.AllPDBSelector())

    for out_dir in out_dirs:
        print(out_dir)
        out_id = out_dir.stem.split("_")[1]

        for pdb_file in Path(out_dir, "pdbs").glob("*.pdb"):
            pdb_id = int(pdb_file.stem)

            m = IMP.Model()
            h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

            align_imp.align_h(
                h=h,
                h_ref=h_ref
            )

            save_pdb_file = Path(job_dir, "align/{}_{}.pdb".format(out_id, pdb_id))
            IMP.atom.write_pdb(
                mhd=h,
                out=str(save_pdb_file)
            )