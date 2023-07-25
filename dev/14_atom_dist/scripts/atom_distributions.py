from pathlib import Path
import sys
import multiprocessing

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import average_structure
import align_imp


if __name__ == "__main__":
    hs = list()
    ms = list()
    pdb_files = list()
    for i in range(10):
        pdb_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/51_synth_sa_5/9485465/output_{}/pdbs".format(i))
        pdb_files.extend(pdb_dir.glob("*.pdb"))

    # pdb_files = list(pdb_dir.glob("*.pdb"))
    for i in range(len(pdb_files)):
        pdb_file = pdb_files[i]
        print(pdb_file)
        m = IMP.Model()
        h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        hs.append(h)
        ms.append(m)

    coord_avg_dict = average_structure.get_coord_avg_dict(
        hs=hs
    )

    ref_pdb_file = str(Path(Path.home(), "xray/data/pdbs/7mhf/7mhf_refine.pdb"))
    ref_m = IMP.Model()
    ref_h = IMP.atom.read_pdb(ref_pdb_file, ref_m, IMP.atom.AllPDBSelector())

    avg_m = IMP.Model()
    avg_h = IMP.atom.read_pdb(ref_pdb_file, avg_m, IMP.atom.AllPDBSelector())

    for pid in coord_avg_dict.keys():
        avg_coords = coord_avg_dict[pid]

        xyz = IMP.core.XYZ(avg_m, pid)
        xyz.set_coordinates(avg_coords)

    print(coord_avg_dict.values())

    rmsd = IMP.atom.get_rmsd(IMP.atom.Selection(ref_h), IMP.atom.Selection(avg_h))
    print(rmsd)



