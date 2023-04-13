import IMP, IMP.core, IMP.atom
import math, time, sys, pickle
from pathlib import Path
import shutil
sys.path.append(str(Path(Path.home(), "xray/src")))
from log_statistics import LogStatistics
from process_pdb import get_crystal_info
from xray_restraint import XtalRestraint
from charmm_restraint import charmm_restraints
import pandas as pd
import os


def single_refinement(
        output_dir,
        m,
        h,
        uc_dim,
        sg_symbol,
        cif_file,
        xray_sf,
        w_xray,
        algo,
        steps,
        ref_pdb_file,
        params_file
):
    param_dict = locals()
    param_f = open(params_file, "a")
    for key in param_dict.keys():
        print("{:<15}{}\n".format(key, param_dict[key]))
        param_f.write("{:<15}{}\n".format(key, param_dict[key]))
    param_f.write("\n\n")
    param_f.close()

    # w_xray = math.pow(10, 4)
    # T = 300

    ps = list()
    pids = list(IMP.atom.Selection(h).get_selected_particle_indexes())

    r_xtal = XtalRestraint(
        m=m,
        ps=ps,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        f_obs_file=str(cif_file),
        scale=True,
        target=xray_sf
    )

    rs_xtal = IMP.RestraintSet(m, w_xray)
    rs_xtal.add_restraint(r_xtal)

    rs_charmm = IMP.RestraintSet(m, 1.0)
    rs_charmm.add_restraints(
        charmm_restraints(
            m,
            h,
            eps=False
        )
    )

    rs = list()
    rs.extend(rs_xtal.get_restraints())
    rs.extend(rs_charmm.get_restraints())

    sf = IMP.core.RestraintsScoringFunction([rs_xtal, rs_charmm])
    # sf = IMP.core.RestraintsScoringFunction([rs_xtal])

    if algo == "cg":
        sampler = IMP.core.ConjugateGradients(m)
    elif algo == "md":
        T, t_step = 300, 2
        sampler = IMP.atom.MolecularDynamics(m)
        s_v = IMP.atom.VelocityScalingOptimizerState(m, ps, T)
        sampler.add_optimizer_state(s_v)
        sampler.set_particles(ps)

    sampler.set_scoring_function(sf)
    sampler.set_has_required_score_states(True)

    # Add logs.
    log_ostate = LogStatistics(
        m=m,
        h=h,
        ref_pdb_file=ref_pdb_file,
        restraints=rs,
        freq=1
    )
    sampler.add_optimizer_state(log_ostate)

    out_file = Path(output_dir, "pdbs/%1%.pdb")
    pdb_ostate = IMP.atom.WritePDBOptimizerState(m, pids, str(out_file))
    pdb_ostate.set_period(1)
    sampler.add_optimizer_state(pdb_ostate)

    if algo == "cg":
        sampler.optimize(steps)
    elif algo == "md":
        sampler.setup(ps)
        sampler.set_temperature(T)
        sampler.assign_velocities(T)
        sampler.set_maximum_time_step(t_step)
        sampler.simulate(steps * t_step)

    log = log_ostate.get_log()
    log_file = Path(output_dir, "log.csv")
    log_df = pd.DataFrame(log)
    log_df.to_csv(log_file)

    return


if __name__ == "__main__":
    job_dir = Path(Path.home(), "xray/sample_benchmark/trivial_refine_weights")

    runs = dict()

    for i in range(6):
        runs[i] = (i, "ls", "cg")
    for i in range(6,12):
        runs[i] = (i%6, "ml", "cg")

    # w_xray, xray_sf, algo = 1, "ml", "cg"
    # output_dir = Path(job_dir, "output_0")

    w_xray, xray_sf, algo = runs[int(sys.argv[1])]
    steps = 100
    # if algo == "cg":
    #     steps = 1000
    # elif algo == "md":
    #     steps = 1

    print(w_xray, xray_sf, algo)
    output_dir = Path(job_dir, "output_{}".format(int(sys.argv[1])))

    input_pdb_file = Path(Path.home(), "xray/data/pdbs/2h2z_clean_c121.pdb")
    ref_pdb_file = Path(Path.home(), "xray/data/pdbs/2h2z_clean_c121.pdb")

    cif_file = Path(Path.home(), "xray/data/reflections/2h2z_clean_c121_synth_exp.cif")
    # uc_dim, sg_symbol = get_crystal_info(pdb_file)
    uc_dim, sg_symbol = (108.406, 81.750, 53.553, 90.00, 104.70, 90.00), "C 1 2 1"

    m = IMP.Model()
    sel = IMP.atom.AllPDBSelector()
    h = IMP.atom.read_pdb(str(input_pdb_file), m, sel)

    for pid in m.get_particle_indexes():
        if IMP.atom.Atom.get_is_setup(m, pid):
            IMP.core.XYZR(m, pid).set_coordinates_are_optimized(True)

    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    Path(output_dir, "pdbs").mkdir(parents=True)
    params_file = Path(output_dir, "params.txt")

    single_refinement(
        output_dir=output_dir,
        m=m,
        h=h,
        uc_dim=uc_dim,
        sg_symbol=sg_symbol,
        cif_file=cif_file,
        xray_sf=xray_sf,
        w_xray=w_xray,
        algo=algo,
        steps=steps,
        ref_pdb_file=ref_pdb_file,
        params_file=params_file
    )
