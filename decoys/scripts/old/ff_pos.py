import IMP, IMP.core, IMP.atom
import math, time, sys, pickle
from pathlib import Path
sys.path.append(str(Path(Path.home(), "xray/src")))
from log_statistics import LogStatistics
from log_statistics import LogRMSD
from xray_restraint import XtalRestraint
from charmm_restraint import charmm_restraints
from position_restraint import position_restraints
import shutil


def md_sim(
        pdb_file,
        m, h,
        uc_dim,
        sg_symbol,
        cif_file,
        output_dir,
        k
):
    # CHARMM restraint needs to be first because it adds an atom.
    rs_charmm = IMP.RestraintSet(m, 1.0)
    [rs_charmm.add_restraint(r_charmm) for r_charmm in charmm_restraints(m,h)]

    rs_position = IMP.RestraintSet(m, 1.0)
    [rs_position.add_restraint(r_pos) for r_pos in position_restraints(m,h,k)]

    atoms = IMP.atom.Selection(h)
    pids = atoms.get_selected_particle_indexes()
    [IMP.atom.LinearVelocity.setup_particle(m, pid) for pid in pids]
    ps = [m.get_particle(pid) for pid in pids]

    rs_xray = IMP.RestraintSet(m, 1.0)
    r_xray = XtalRestraint(
        m=m,
        ps=ps,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file)
    )
    rs_xray.add_restraint(r_xray)

    sf = IMP.core.RestraintsScoringFunction([rs_charmm, rs_position])

    # # Slightly relax the structure.
    # # log_ostate = LogStatistics(m, rs, print=False)
    # # cg = IMP.core.ConjugateGradients(m)
    # # cg.add_optimizer_state(log_ostate)
    # # cg.set_scoring_function(sf)
    # # cg.optimize(50)

    # Sample with molecular dynamics.
    md = IMP.atom.MolecularDynamics(m)
    md.set_particles(ps)
    md.set_scoring_function(sf)
    md.setup(ps)
    md.set_has_required_score_states(True)
    s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, 100)
    md.add_optimizer_state(s_v)

    # Add logs.
    log_ostate = LogStatistics(
        m,
        [r for r in rs_charmm.get_restraints()],
        print=True,
        freq=50
    )
    md.add_optimizer_state(log_ostate)

    log_rmsd_ostate = LogRMSD(
        m=m,
        h=h,
        ref_pdb_file=pdb_file,
        print=True,
        freq=50
    )
    md.add_optimizer_state(log_rmsd_ostate)

    # Warmup system.
    md.set_maximum_time_step(1)
    md.set_temperature(100)
    md.assign_velocities(100)
    md.simulate(100)

    # Add pdb writer.
    save_file = Path(output_dir, "pdbs_{}/md_%1%.pdb".format(k))
    pdb_ostate = IMP.atom.WritePDBOptimizerState(
        m,
        pids,
        str(save_file)
    )
    pdb_ostate.set_period(1000)
    md.add_optimizer_state(pdb_ostate)

    # Run simulation.
    steps = 100000
    t0 = time.time()
    s_v.set_temperature(298)
    md.set_maximum_time_step(2)
    md.set_temperature(298)
    md.assign_velocities(298)
    md.simulate(steps * 2)
    print(time.time()-t0)

    log = log_ostate.get_log()
    log_file = Path(output_dir, "md_{}.log".format(k))
    with open(log_file, 'wb') as file:
        pickle.dump(log, file)


if __name__ == "__main__":
    job_name = "4n7f_ff_pos"
    xray_dir = Path(Path.home(), "xray")
    pdb_file = Path(xray_dir, "refinement/4n7f_clean/pdbs/refine_6.pdb")
    cif_file = Path(xray_dir, "data/reflections/4n7f_clean_1.cif")
    output_dir = Path(xray_dir, "md_xray/{}".format(job_name))

    # if output_dir.exists():
    #     shutil.rmtree(path=output_dir)

    # output_dir.mkdir(
    #     exist_ok=True
    # )
    k = int(sys.argv[1])
    Path(output_dir, "pdbs_{}".format(k)).mkdir()

    print("job_name:   {}".format(job_name))
    print("pdb_file:   {}".format(pdb_file))
    print("cif_file:   {}".format(cif_file))
    print("output_dir: {}".format(output_dir))
    print("k:          {}".format(k))

    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.ATOMPDBSelector())

    md_sim(
        pdb_file=pdb_file,
        m=m,
        h=h,
        uc_dim=(30,30,30,90,90),
        sg_symbol="P 1",
        cif_file=cif_file,
        output_dir=output_dir,
        k=k
    )