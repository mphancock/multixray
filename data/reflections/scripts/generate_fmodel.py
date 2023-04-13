from pathlib import Path
import sys

import mmtbx.f_model
import cctbx.crystal
import cctbx.xray
import iotbx

sys.path.append(str(Path(Path.home(), "xray/src")))
import miller_ops


def get_f_model(
        pdb_file,
        uc_dim,
        sg_symbol,
        cif_file,
        res
):
    crystal_symmetry = cctbx.crystal.symmetry(
        unit_cell=uc_dim,
        space_group_symbol=sg_symbol
    )

    xray_structure = iotbx.pdb.input(str(pdb_file)).xray_structure_simple(
        crystal_symmetry=crystal_symmetry
    )

    f_obs = miller_ops.get_f_obs(str(cif_file))
    f_obs = miller_ops.filter_f_obs_resolution(
        f_obs=f_obs,
        d_min=res,
        d_max=None
    )

    f_model_manager = mmtbx.f_model.manager(
        xray_structure=xray_structure,
        f_obs=f_obs,
        target_name="ls"
    )

    f_model = f_model_manager.f_model()

    return f_model


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/remove_S_side_chains_3ca7_clean.pdb")

    cif_file = Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")
    model_cif_file = Path(Path.home(), "xray/data/reflections/3ca7/remove_S_side_chains_3ca7_clean.cif")

    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg_symbol = "C 1 2 1"

    f_model = get_f_model(
        pdb_file=pdb_file,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        cif_file=cif_file,
        res=None
    )
    print(pdb_file)
    print(model_cif_file)
    print(f_model.size())

    f = open(model_cif_file, "w")
    f_model.amplitudes().as_cif_simple(
        array_type="meas",
        out=f
    )
    f.close()


