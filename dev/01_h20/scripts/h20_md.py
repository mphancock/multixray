from pathlib import Path
import sys
import shutil

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import pdb_optimizer_state
import trackers
import log_statistics
import xray_restraint


if __name__ == "__main__":
    ## READ IN WATER.
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhk/7mhk_clean_h20.pdb")
    s = IMP.atom.AllPDBSelector()
    m = IMP.Model()

    h = IMP.atom.read_pdb(str(pdb_file), m, s)

    h20s = IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("HET: O  ")).get_selected_particle_indexes()
    print(len(h20s))

    for pid in h20s:
        IMP.atom.CHARMMAtom.setup_particle(m, pid, "O")
        print(IMP.atom.CHARMMAtom(m, pid).get_charmm_type())

    print(len(IMP.atom.Selection(h).get_selected_particle_indexes()))
    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmm.charmm_restraints(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)
    print(len(IMP.atom.Selection(h).get_selected_particle_indexes()))

    # Configure the IMP model based on the CHARMM parameterization.
    # ff = IMP.atom.get_heavy_atom_CHARMM_parameters()
    # topology = ff.create_topology(h)

    # # topology.apply_default_patches()
    # topology.setup_hierarchy(h)

    # print(len(IMP.atom.Selection(h).get_selected_particle_indexes()))
    # IMP.atom.remove_charmm_untyped_atoms(h)
    # print(len(IMP.atom.Selection(h).get_selected_particle_indexes()))
    # # topology.add_missing_atoms(h)
    # topology.add_coordinates(h)
    # bonds = topology.add_bonds(h)
    # angles = ff.create_angles(bonds)
    # dihedrals = ff.create_dihedrals(bonds)
    # impropers = topology.add_impropers(h)
    # charges = topology.add_charges(h)

    # for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
    #     at = IMP.atom.Atom(m, pid)
    #     print(at.get_atom_type())

    # RUN MOLECULAR DYNAMICS.
    output_dir = Path(Path.home(), "xray/dev/01_h20/data/output_0")
    shutil.rmtree(output_dir, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    T = 300
    # pids = h20s
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    ps = [m.get_particle(pid) for pid in pids]

    r_xtal = xray_restraint.XtalRestraint(
        m=m,
        pids=pids,
        uc_dim=(58.305,36.154,25.362,90.00,103.09,90.00),
        sg_symbol="C 1 2 1",
        f_obs_file=str(Path(Path.home(), "xray/data/reflections/3ca7/3ca7.cif")),
        d_min=0,
        d_max=None,
        scale=True,
        target="ml",
        w_xray=1,
        dynamic_w=1
    )

    rs_xtal = IMP.RestraintSet(m, 1.0)
    rs_xtal.add_restraint(r_xtal)

    # rs = [rset_charmm]
    rs = [rset_charmm, r_xtal]
    sf = IMP.core.RestraintsScoringFunction(rs)

    md = IMP.atom.MolecularDynamics(m)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, T)
    md.add_optimizer_state(s_v)
    md.set_particles(ps)
    md.set_scoring_function(sf)
    md.set_has_required_score_states(True)

    all_trackers = list()

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
    for i in range(len(charmm_rs)):
        charmm_r = charmm_rs[i]
        ff_tracker = trackers.fTracker(
            name=charmm_r.get_name(),
            r=charmm_r
        )
        all_trackers.append(ff_tracker)

    xray_tracker = trackers.fTracker(
        name="xray",
        r=r_xtal
    )
    all_trackers.append(xray_tracker)

    r_work_tracker = trackers.RFactorTracker(
        name="r_work",
        r_xray=r_xtal,
        stat="r_work"
    )
    all_trackers.append(r_work_tracker)

    r_free_tracker = trackers.RFactorTracker(
        name="r_free",
        r_xray=r_xtal,
        stat="r_free"
    )
    all_trackers.append(r_free_tracker)

    o_states = list()
    pdb_dir = Path(output_dir, "pdbs")
    pdb_dir.mkdir(parents=False, exist_ok=True)
    pdb_ostate = pdb_optimizer_state.WriteMultiStatePDBOptimizerState(
        m=m,
        hs=[h],
        pdb_dir=pdb_dir
    )
    pdb_ostate.update()
    pdb_ostate.set_period(10)
    o_states.append(pdb_ostate)

    log_file=Path(output_dir, "log.txt")
    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=log_file,
        log_freq=100
    )
    o_states.append(log_ostate)

    for r in rs:
        r.evaluate(calc_derivs=True)
    log_ostate.update()

    for o_state in o_states:
        md.add_optimizer_state(o_state)

    md.setup(ps)
    md.set_temperature(T)
    for pid in pids:
        IMP.atom.LinearVelocity.setup_particle(m, pid)

    md.assign_velocities(T)
    md.set_maximum_time_step(2)
    md.simulate(1e99)