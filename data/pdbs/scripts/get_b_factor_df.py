from pathlib import Path
import pandas as pd

import IMP
import IMP.atom


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7.pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    b_factor_file = Path(Path.home(), "xray/data/pdbs/3ca7/etc/b_factor.csv")
    b_factor_df = pd.DataFrame(index=list(range(48,98)),columns=["b_factor"])
    for res_id in range(48,98):
        s = IMP.atom.Selection(h, residue_index=res_id)
        pids = s.get_selected_particle_indexes()

        avg_b_factor = 0
        for pid in pids:
            at = IMP.atom.Atom(m, pid)
            avg_b_factor = avg_b_factor + at.get_temperature_factor()

        avg_b_factor = avg_b_factor / len(pids)

        b_factor_df.loc[res_id, "b_factor"] = avg_b_factor
        print(res_id, avg_b_factor)

    b_factor_df.to_csv(b_factor_file)