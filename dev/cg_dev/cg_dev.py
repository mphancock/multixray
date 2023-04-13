import sys
from pathlib import Path
import shutil
import IMP
import IMP.core
import IMP.atom
sys.path.append(str(Path(Path.home(), "xray/src")))
from log_statistics import XYZTracker, dXYZTracker, fTracker, dfTracker, RMSDTracker, LogStatistics
from process_pdb import get_crystal_info
from xray_restraint import XtalRestraint
from charmm import charmm_restraints
from derivatives import get_dfs, set_dfs, get_df_mag_ratio


if __name__ == "__main__":
    job_dir = Path(Path.home(), "xray/dev/cg_dev")
    output_dir = Path(job_dir, "output")
    output_dir.mkdir(exist_ok=True)

    data_dir = Path(Path.home(), "xray/data")
    # input_pdb_file = Path(job_dir, "4n7f_heavy_perturb.pdb")
    input_pdb_file = Path(data_dir, "pdbs/4n7f_heavy.pdb")
    ref_pdb_file = Path(data_dir, "pdbs/4n7f_heavy.pdb")

    cif_file = Path(data_dir, "reflections/4n7f_heavy.cif")
    # uc_dim, sg_symbol = get_crystal_info(pdb_file)
    uc_dim, sg_symbol = (68.411, 68.411, 37.248, 90.00, 90.00,  90.00), "P 41 2 2"

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(input_pdb_file), m, s)

    for pid in m.get_particle_indexes():
        if IMP.atom.Atom.get_is_setup(m, pid):
            IMP.core.XYZR(m, pid).set_coordinates_are_optimized(True)

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    r_xtal = XtalRestraint(
        m=m,
        # ps=ps,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        scale=True,
        target="ls"
    )

    rs_xtal = IMP.RestraintSet(m, 1.0)
    rs_xtal.add_restraint(r_xtal)

    rs_charmm = IMP.RestraintSet(m, 1.0)
    rs_charmm.add_restraints(
        charmm_restraints(
            m,
            h,
            eps=False
        )
    )

    mag_ratio = get_df_mag_ratio(
        m=m,
        pids=pids,
        r1=rs_charmm,
        r2=r_xtal
    )
    r_xtal.set_weight(1*mag_ratio)
    print(r_xtal.get_weight())
    print(rs_xtal.get_weight())

    rs = list()
    # rs.extend(rs_xtal.get_restraints())
    rs.extend(rs_charmm.get_restraints())

    sf = IMP.core.RestraintsScoringFunction([rs_xtal, rs_charmm])
    # sf = IMP.core.RestraintsScoringFunction([rs_xtal])
    # sf = IMP.core.RestraintsScoringFunction([rs_charmm])

    cg = IMP.core.ConjugateGradients(m)

    cg.set_scoring_function(sf)
    cg.set_has_required_score_states(True)

    # Add logs.
    pids_track = IMP.atom.Selection(
        hierarchy=h,
        atom_type=IMP.atom.AtomType("CA")
    ).get_selected_particle_indexes()

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)
    pids_0 = list(IMP.atom.Selection(h_0).get_selected_particle_indexes())

    trackers = list()
    xyz_tracker = XYZTracker(
        name="AT1 x",
        pid=pids[0],
        at_id=0
    )
    trackers.append(xyz_tracker)

    dxyz_tracker = dXYZTracker(
        name="AT1 dx",
        pid=pids[0],
        at_id=0
    )
    trackers.append(dxyz_tracker)

    ff_tracker = fTracker(
        name="ff",
        r=rs_charmm
    )
    trackers.append(ff_tracker)

    dff_tracker = dfTracker(
        name="AT1 dffdx",
        r=rs_charmm,
        pid=pids[0],
        at_id=0
    )
    trackers.append(dff_tracker)

    xray_tracker = fTracker(
        name="xray",
        r=r_xtal
    )
    trackers.append(xray_tracker)

    dxray_tracker = dfTracker(
        name="AT1 dxraydx",
        r=r_xtal,
        pid=pids[0],
        at_id=0
    )
    trackers.append(dxray_tracker)

    rmsd_tracker = RMSDTracker(
        name="RMSD_tracker",
        m_0=m_0,
        pids_0=pids_0
    )

    log_ostate = LogStatistics(
        m=m,
        pids=pids,
        trackers=[dff_tracker, dxray_tracker, dxyz_tracker, rmsd_tracker]
    )
    cg.add_optimizer_state(log_ostate)

    # pdb_dir = Path(output_dir, "pdbs")
    # pdb_dir.mkdir(exist_ok=True)
    # out_file = Path(pdb_dir, "%1%.pdb")
    # pdb_ostate = IMP.atom.WritePDBOptimizerState(
    #     m,
    #     pids,
    #     str(out_file)
    # )
    # pdb_ostate.set_period(1)
    # cg.add_optimizer_state(pdb_ostate)

    # cg.set_log_level(IMP.PROGRESS)
    # cg.set_gradient_threshold(1e-6)
    cg.do_optimize(100)
    # log = log_ostate.get_log()
    # log_file = Path(output_dir, "log.csv")
    # log_df = pd.DataFrame(log)
    # log_df.to_csv(log_file)

    # this is needed to clean up memory properly for some reason
    # print(log_ostate.get_n_steps())
    print(log_ostate.get_number_of_updates())
    cg.remove_optimizer_state(log_ostate)

    del log_ostate
    del m