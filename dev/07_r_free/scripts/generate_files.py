from pathlib import Path

import IMP
import IMP.atom


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhf.pdb")

    for i in range(16):
        h20 = i // 8 % 2
        occ = i // 4 % 2
        h = i // 2 % 2
        b = i // 1 % 2
        print(h20, occ, h, b)

        file_name = "7mhf_h20_{}_occ_{}_h_{}_b_{}".format(h20, occ, h, b)

        new_pdb_file = Path(Path.home(), "xray/dev/07_r_free/data/{}.pdb".format(file_name))

        sels = [IMP.atom.NonWaterNonHydrogenPDBSelector()]

        if h20:
            sels.append(IMP.atom.WaterPDBSelector())
        if h:
            sels.append(IMP.atom.HydrogenPDBSelector())

        sel_all = sels[0]
        for sel in sels:
            sel_all = IMP.atom.OrPDBSelector(sel_all, sel)

        m = IMP.Model()
        h = IMP.atom.read_pdb(str(pdb_file), m, sel_all)
        residues = list(range(307))
        residues.extend(list(range(601,677)))
        sel_res = IMP.atom.Selection(h, residue_indexes=residues)
        pids = sel_res.get_selected_particle_indexes()

        if not occ:
            for pid in pids:
                atom = IMP.atom.Atom(m, pid)
                atom.set_occupancy(1)

        if not b:
            for pid in pids:
                atom = IMP.atom.Atom(m, pid)
                atom.set_temperature_factor(0)

        IMP.atom.write_pdb(sel_res, str(new_pdb_file))