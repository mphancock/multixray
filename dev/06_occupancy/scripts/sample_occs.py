from pathlib import Path
import sys
import argparse
import shutil

import IMP
import IMP.atom
import IMP.core
import IMP.isd
import IMP.algebra

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import xray_restraint
import trackers
import com_restraint
import log_statistics
import align_imp
import pdb_writer
import com_optimizer_state
import copy_pdbs
import update_weights_optimizer_state
import molecular_dynamics
import params


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_refine_2.pdb")
    ref_pdb_file = pdb_file

    m = IMP.Model()
    s = IMP.atom.AllPDBSelector()
    n_states = 2
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, s)

    occs = list()
    occs_0 = [.5,.5]

    for i in range(len(hs)):
        h = hs[i]
        pids = IMP.atom.Selection(h).get_selected_particle_indexes()
        for pid in pids:
            IMP.atom.Atom(m, pid).set_occupancy(occs_0[i])

    print("weights")
    w_p = IMP.Particle(m, "weights")
    w_pid = IMP.isd.Weight.setup_particle(w_p, IMP.algebra.VectorKD(occs_0))
    w = IMP.isd.Weight(m, w_pid)

    for i in range(n_states):
        print(w.get_weight(i))

    ## SCORING
    rs = list()
    rset_charmm = IMP.RestraintSet(m, 1.0)
    for h in hs:
        charmm_rs = charmm.charmm_restraints(
            m,
            h,
            eps=False
        )
        rset_charmm.add_restraints(charmm_rs)
    rset_charmm.set_weight(1/n_states)
    rs.append(rset_charmm)

    r_xray = xray_restraint.XtalRestraint(
        m=m,
        n_state=n_states,
        pids=pids,
        f_obs_file=Path(Path.home(), "xray/data/reflections/3ca7/3ca7_refine_2.cif"),
        d_min=0,
        d_max=None,
        scale=True,
        target="ml",
        w_xray=1.0,
        dynamic_w=1
    )

    rs_xray = IMP.RestraintSet(m, 1.0)
    rs_xray.add_restraint(r_xray)
    rs.append(rs_xray)

    m_0 = IMP.Model()
    s_0 = IMP.atom.AllPDBSelector()
    h0_s = IMP.atom.read_multimodel_pdb(str(ref_pdb_file), m_0, s_0)

    pids_ca_0 = list(IMP.atom.Selection(h0_s, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes())
    com_0 = IMP.atom.CenterOfMass.setup_particle(
        IMP.Particle(m_0),
        pids_ca_0
    )
    print(com_0.get_coordinates())

    all_trackers = list()
    pids_0 = list(IMP.atom.Selection(h0_s).get_selected_particle_indexes())

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
        r=r_xray
    )
    all_trackers.append(xray_tracker)

    r_work_tracker = trackers.RFactorTracker(
        name="r_work",
        r_xray=r_xray,
        stat="r_work"
    )
    all_trackers.append(r_work_tracker)

    r_free_tracker = trackers.RFactorTracker(
        name="r_free",
        r_xray=r_xray,
        stat="r_free"
    )
    all_trackers.append(r_free_tracker)

    rmsd_tracker = trackers.RMSDTracker(
        name="rmsd",
        hs=hs,
        hs_0=h0_s,
        ca_only=True
    )
    all_trackers.append(rmsd_tracker)

    rmsd_all_tracker = trackers.RMSDTracker(
        name="rmsd_all",
        hs=hs,
        hs_0=h0_s,
        ca_only=False
    )
    all_trackers.append(rmsd_all_tracker)

    for i in range(n_states):
        state_pids = IMP.atom.Selection(hs[i]).get_selected_particle_indexes()
        occ_tracker = trackers.OccTracker(
            name="occ_{}".format(i),
            m=m,
            at=IMP.atom.Atom(m, state_pids[0])
        )
        all_trackers.append(occ_tracker)

    out_dir = Path("/wynton/group/sali/mhancock/xray/dev/06_occupancy/data/output_0")
    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=Path(out_dir, "log.csv"),
        log_freq=100
    )
    for r in rs:
        r.evaluate(calc_derivs=True)
    log_ostate.update()

    o_states = list()
    o_states.append(log_ostate)

    w_os = update_weights_optimizer_state.UpdateWeightsOptimizerState(
        m=m,
        hs=hs,
        w=w,
        r_xtal=r_xray
    )
    w_os.set_period(1)
    o_states.append(w_os)

    molecular_dynamics.molecular_dynamics(
        output_dir=out_dir,
        hs=hs,
        rs=rs,
        T=300,
        t_step=2,
        steps=-1,
        sa_sched=None,
        o_states=o_states,
        md_ps=[]
    )