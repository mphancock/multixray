from pathlib import Path
import pickle
import sys

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from cctbx_score import get_score
import miller_ops


if __name__ == "__main__":
    m = IMP.Model()
    sel = IMP.atom.AllPDBSelector()

    pdb_file = Path(Path.home(), "xray/tmp/output_0/pdbs/10.pdb")

    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, sel)

    pids = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

    print(len(pids))

    file = open(Path(Path.home(), "xray/tmp/coords_dict.pkl"), 'rb')

    # dump information to that file
    coords_dict = pickle.load(file)
    for pid in pids:
        xyz = IMP.core.XYZR(m, pid)
        xyz.set_coordinates(coords_dict[pid])

    cif_file = Path(Path.home(), "xray/dev/29_synthetic_native_3/data/cifs/7mhf_30/0/0.cif")
    f_obs_array = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.F_meas_au"
    )
    f_obs_array = miller_ops.clean_miller_array(f_obs_array)

    # Set flags from file.
    status_array = miller_ops.get_miller_array(
        f_obs_file=cif_file,
        label="_refln.status"
    )
    flags_array = status_array.customized_copy(data=status_array.data()=="f")
    f_obs_array, flags_array = f_obs_array.common_sets(other=flags_array)

    d_min = 1.5
    f_obs = miller_ops.filter_f_obs_resolution(
            f_obs=f_obs_array,
            d_max=None,
            d_min=d_min
    )
    flags = miller_ops.filter_f_obs_resolution(
            f_obs=flags_array,
            d_max=None,
            d_min=d_min
    )


    results_dict = get_score(
        hs=hs,
        occs=[1],
        pids=pids,
        f_obs=f_obs,
        r_free_flags=flags,
        target="ml",
        update_scale=True,
        update_k1=True,
        u_aniso_file=None,
        delta=None
    )
    print(results_dict["r_free"])