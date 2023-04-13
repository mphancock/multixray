from pathlib import Path
import sys

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import cctbx_scores
import xray_restraint


if __name__ == "__main__":
    xray_dir = Path(Path.home(), "xray")
    pdb_file = Path(Path.home(), "xray/decoys/data/decoy_sets/3ca7/3ca7_N_1000_x1/315.pdb")
    cif_file = Path(xray_dir, "data/reflections/3ca7/3ca7_clean.cif")

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()

    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg_symbol = "C 1 2 1"
    d_min = None

    r_xray = xray_restraint.XtalRestraint(
        m=m,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=cif_file,
        d_min=d_min,
        d_max=None,
        scale=True,
        target="ml"
    )

    print(r_xray.evaluate(True))
