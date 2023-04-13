from pathlib import Path
import sys
import shutil

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/decoys/src")))
import gen_decoys


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")

    T = 300
    job_dir = Path(xray_dir, "decoys/md_out/3ca7_{}".format(T))

    for i in range(50):
        if i < 25:
            n_step = 200
            pdb_write_freq = 2
        else:
            n_step = 2000
            pdb_write_freq = 20

        output_dir = Path(job_dir, "output_{}".format(i))
        print(output_dir)

        if output_dir.exists():
            shutil.rmtree(path=output_dir)

        output_dir.mkdir(
            exist_ok=True
        )
        Path(output_dir, "pdbs").mkdir()

        pdb_file = Path(xray_dir, "data/pdbs/3ca7/3ca7_clean.pdb")

        m = IMP.Model()
        h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.ATOMPDBSelector())

        m_0 = IMP.Model()
        h_0 = IMP.atom.read_pdb(str(pdb_file), m_0, IMP.atom.ATOMPDBSelector())

        print(len(m.get_particle_indexes()), len(m_0.get_particle_indexes()))

        gen_decoys.gen_decoys(
            output_dir=output_dir,
            h=h,
            h_0=h_0,
            T=T,
            t_step=2,
            n_step=n_step,
            pdb_write_freq=pdb_write_freq
        )