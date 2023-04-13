import sys
from pathlib import Path
import shutil

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import log_statistics
import xray_restraint
import charmm
import derivatives
import params


if __name__ == "__main__":
    job_dir = Path(Path.home(), "xray/dev/md_dev")
    output_dir = Path(job_dir, "output_0")
    if output_dir.exists():
        shutil.rmtree(output_dir)
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
    ps = [m.get_particle(pid) for pid in pids]
    r_xtal = xray_restraint.XtalRestraint(
        m=m,
        # ps=ps,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        scale=True,
        target="ml"
    )

    rs_xtal = IMP.RestraintSet(m, 1.0)
    rs_xtal.add_restraint(r_xtal)

    rset_charmm = IMP.RestraintSet(m, 1.0)
    rset_charmm.add_restraints(
        charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
    )

    # mag_ratio = get_df_mag_ratio(
    #     m=m,
    #     pids=pids,
    #     r1=rset_charmm,
    #     r2=r_xtal
    # )
    # r_xtal.set_weight(1*mag_ratio)
    # print(r_xtal.get_weight())
    # print(rs_xtal.get_weight())
    r_xtal.set_weight(100)

    rs = list()
    # rs.extend(rs_xtal.get_restraints())
    rs.extend(rset_charmm.get_restraints())

    sf = IMP.core.RestraintsScoringFunction([rs_xtal, rset_charmm])
    # sf = IMP.core.RestraintsScoringFunction([rs_xtal])
    # sf = IMP.core.RestraintsScoringFunction([rset_charmm])

    T, t_step = 600, 2
    md = IMP.atom.MolecularDynamics(m)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, T)
    md.add_optimizer_state(s_v)
    md.set_particles(ps)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    # Add all optimizer states.
    pids_track = IMP.atom.Selection(
        hierarchy=h,
        atom_type=IMP.atom.AtomType("CA")
    ).get_selected_particle_indexes()

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(ref_pdb_file), m_0, s_0)
    pids_0 = list(IMP.atom.Selection(h_0).get_selected_particle_indexes())

    trackers = list()
    xyz_tracker = log_statistics.XYZTracker(
        name="AT1 x",
        pid=pids[0],
        at_id=0
    )
    trackers.append(xyz_tracker)

    dxyz_tracker = log_statistics.dXYZTracker(
        name="AT1 dx",
        pid=pids[0],
        at_id=0
    )
    trackers.append(dxyz_tracker)

    ff_tracker = log_statistics.fTracker(
        name="ff",
        r=rset_charmm
    )
    trackers.append(ff_tracker)

    dff_tracker = log_statistics.dfTracker(
        name="AT1 dffdx",
        r=rset_charmm,
        pid=pids[0],
        at_id=0
    )
    trackers.append(dff_tracker)

    xray_tracker = log_statistics.fTracker(
        name="xray",
        r=r_xtal
    )
    trackers.append(xray_tracker)

    dxray_tracker = log_statistics.dfTracker(
        name="AT1 dxraydx",
        r=r_xtal,
        pid=pids[0],
        at_id=0
    )
    trackers.append(dxray_tracker)

    rmsd_tracker = log_statistics.RMSDTracker(
        name="RMSD",
        m_0=m_0,
        pids_0=pids_0
    )
    trackers.append(rmsd_tracker)

    df_mag_tracker = log_statistics.dfMagTracker(
        name="df mag",
        r1=rset_charmm,
        r2=r_xtal
    )
    trackers.append(df_mag_tracker)

    # log_ostate = log_statistics.LogStatistics(
    #     m=m,
    #     pids=pids,
    #     trackers=[xray_tracker]
    # )
    # md.add_optimizer_state(log_ostate)

    # weight_ostate = derivatives.RestraintWeightOptimizerState(
    #     m=m,
    #     pids=pids,
    #     r1=rset_charmm,
    #     r2=r_xtal
    # )
    # weight_ostate.set_period(1)
    # md.add_optimizer_state(weight_ostate)

    pdb_dir = Path(output_dir, "pdbs")
    pdb_dir.mkdir(exist_ok=True)
    out_file = Path(pdb_dir, "%1%.pdb")
    pdb_ostate = IMP.atom.WritePDBOptimizerState(
        m,
        pids,
        str(out_file)
    )
    pdb_ostate.set_period(1)
    md.add_optimizer_state(pdb_ostate)

    steps = 5
    md.setup(ps)
    md.set_temperature(T)
    [IMP.atom.LinearVelocity.setup_particle(m, pid) for pid in pids]
    md.assign_velocities(T)
    md.set_maximum_time_step(t_step)
    md.simulate(steps * t_step)

    # This is needed to clean up memory properly for some reason.
    # md.remove_optimizer_state(pdb_ostate)
    # md.remove_optimizer_state(weight_ostate)

    # del weight_ostate
    # del log_ostate

    # Clean up memory
    for tracker in trackers:
        del tracker
    # md.remove_optimizer_state(log_ostate)
    md.remove_optimizer_state(pdb_ostate)
    del pdb_ostate
    del md
    del m_0
    del m