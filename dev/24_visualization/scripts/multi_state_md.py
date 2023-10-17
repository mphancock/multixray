from pathlib import Path
import shutil
import math

import IMP
import IMP.atom


if __name__ == "__main__":
    out_file = Path(Path.home(), "xray/tmp/md_tmp.pdb")
    all_hs = list()
    all_ms = list()

    pdb_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/60_7mhj_2/1018128/output_639/pdbs")

    for i in range(0,801):
        print(i)
        m = IMP.Model()

        pdb_file = Path(pdb_dir, "{}.pdb".format(i))

        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        h_0 = hs[0]
        for i in range(1, len(hs)):
            h_0.add_child(IMP.atom.get_root(hs[i]).get_children()[0])

        hchains = list()
        for pid in m.get_particle_indexes():
            if IMP.atom.Chain.get_is_setup(m, pid):
                hchain = IMP.atom.Chain(m, pid)
                # print(hchain.get_id())
                hchains.append(hchain)

        for j in range(len(hchains)):
            hchain = hchains[j]
            chain_id = chr(ord('A') + j)
            hchain.set_id(chain_id)
            # print(hchain.get_id())

        all_hs.append(h_0)
        all_ms.append(m)

    IMP.atom.write_multimodel_pdb(all_hs, str(out_file))