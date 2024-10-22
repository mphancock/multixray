from pathlib import Path

import IMP
import IMP.atom

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
from cctbx_score import get_score
from miller_ops import get_miller_array, clean_miller_array

if __name__ == "__main__":
    pdb_file = Path("../data/pdbs/7mhf.pdb")

    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    pids = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

    print(len(hs))

    f_obs_file = Path(Path.home(), "xray/data/cifs/3k0m/3k0m.cif")

    f_obs = get_miller_array(
            f_obs_file=f_obs_file,
            label="_refln.intensity_meas"
    )
    print(f_obs.indices()[0], f_obs.data()[0])

    f_obs = clean_miller_array(f_obs)
    status_array = get_miller_array(
        f_obs_file=f_obs_file,
        label="_refln.status"
    )
    flags = status_array.customized_copy(data=status_array.data()=="f")
    f_obs, flags = f_obs.common_sets(other=flags)

    # d_min = 2.0
    # f_obs = miller_ops.filter_f_obs_resolution(
    #         f_obs=f_obs,
    #         d_max=None,
    #         d_min=d_min
    # )
    # flags = miller_ops.filter_f_obs_resolution(
    #         f_obs=flags,
    #         d_max=None,
    #         d_min=d_min
    # )

    print(f_obs.indices()[0], f_obs.data()[0])

    for h in hs:
        score_dict = get_score(
            hs=[h],
            occs=[1],
            pids=pids,
            f_obs=f_obs,
            r_free_flags=flags,
            target="ml",
            ab_file=None,
            update_scale=True,
            update_k1=False,
            u_aniso_file=pdb_file,
            delta=None
        )
        print(score_dict["r_free"], score_dict["r_work"], score_dict["score"], score_dict["ff"])
