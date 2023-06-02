from pathlib import Path
import sys
import argparse
import shutil

import IMP
import IMP.atom
import IMP.core

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import xray_restraint
import trackers
import com_restraint
import log_statistics
import align_imp
import pdb_optimizer_state
sys.path.append(str(Path(Path.home(), "xray/sample_bench/src")))
import molecular_dynamics


if __name__ == "__main__":
    output_dir = Path(Path.home(), "xray/dev/04_multi_state_md/output_0")
    pdb_dir = Path(output_dir, "pdbs")
    if pdb_dir.exists():
        shutil.rmtree(pdb_dir)
    pdb_dir.mkdir(exist_ok=True)
    local_pdb_dir = pdb_dir

    ref_pdb_file = Path(Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean_2_state.pdb"))
    pdb_file = ref_pdb_file

    uc_dim = (58.305, 36.154, 25.362, 90.00, 103.09, 90.00)
    sg_symbol = "C 1 2 1"
    cif_file = Path(Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif"))

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, s)

    pids = list()
    for h in hs:
        pids.extend(IMP.atom.Selection(h).get_selected_particle_indexes())

    for pid in pids:
        IMP.atom.Atom(m, pid).set_occupancy(1/len(hs))
    ps = [m.get_particle(pid) for pid in pids]

    n_states = len(hs)

    rset_charmm = IMP.RestraintSet(m, 1.0)
    for h in hs:
        charmm_rs = charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
        rset_charmm.add_restraints(charmm_rs)
    rset_charmm.set_weight(1/n_states)

    r_xtal = xray_restraint.XtalRestraint(
        m=m,
        pids=pids,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        d_min=0,
        d_max=None,
        scale=True,
        target="ml",
        w_xray=30000,
        dynamic_w=0
    )

    rs_xtal = IMP.RestraintSet(m, 1.0)
    rs_xtal.add_restraint(r_xtal)

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h_0 = IMP.atom.read_pdb(str(Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")), m_0, s_0)

    pids_ca_0 = list(IMP.atom.Selection(h_0, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    com_0 = IMP.atom.CenterOfMass.setup_particle(
        IMP.Particle(m_0),
        pids_ca_0
    )
    print(com_0.get_coordinates())

    rs = [rset_charmm, r_xtal]
    for h in hs:
        pids_ca = list(IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())

        r_com = com_restraint.CenterOfMassRestraint(
            m=m,
            pids=pids_ca,
            k=10,
            xyz_0=com_0.get_coordinates()
        )

        rs.append(r_com)

    sf = IMP.core.RestraintsScoringFunction(rs)

    all_trackers = list()
    pids_0 = list(IMP.atom.Selection(h_0).get_selected_particle_indexes())

    step_tracker = trackers.StepTracker(
        name="step",
        m=m
    )
    all_trackers.append(step_tracker)

    time_tracker = trackers.TimeTracker(
        name="time",
        m=m
    )
    all_trackers.append(time_tracker)

    ff_tracker = trackers.fTracker(
        name="ff",
        r=rset_charmm
    )
    all_trackers.append(ff_tracker)

    xray_tracker = trackers.fTracker(
        name="xray",
        r=r_xtal
    )
    all_trackers.append(xray_tracker)

    rmsd_tracker = trackers.RMSDTracker(
        name="rmsd",
        hs=hs,
        hs_0=[h_0],
        align=False
    )
    all_trackers.append(rmsd_tracker)

    rmsd = rmsd_tracker.evaluate()

    md = IMP.atom.MolecularDynamics(m)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, 300)
    md.add_optimizer_state(s_v)
    md.set_particles(ps)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    # Add all optimizer states.
    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=Path(Path.home(), "xray/dev/04_multi_state_md/output_0/log.csv"),
        log_freq=250
    )
    md.add_optimizer_state(log_ostate)

    IMP.atom.write_multimodel_pdb(hs, str(Path(local_pdb_dir, "-1.pdb")))
    # for i in range(n_states):
    #     h = hs[i]
    #     IMP.atom.write_pdb(h, str(Path(local_pdb_dir, "-1.pdb")))

    pdb_ostate = pdb_optimizer_state.WriteMultiStatePDBOptimizerState(
        m=m,
        hs=hs,
        pdb_dir=local_pdb_dir
    )
    pdb_ostate.set_period(10)
    md.add_optimizer_state(pdb_ostate)

    # for i in range(n_states):
    #     out_file = Path(local_pdb_dir, "%1%_{}.pdb".format(i))
    #     h = hs[i]
    #     state_pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    #     pdb_ostate = IMP.atom.WritePDBOptimizerState(
    #         m,
    #         state_pids,
    #         str(out_file)
    #     )
    #     pdb_ostate.set_period(10)
    #     md.add_optimizer_state(pdb_ostate)

    md.setup(ps)
    md.set_temperature(300)
    for pid in pids:
        IMP.atom.LinearVelocity.setup_particle(m, pid)

    md.assign_velocities(300)
    md.set_maximum_time_step(2)

    sf.evaluate(derivatives=True)
    log_ostate.update()

    md.simulate(1000)