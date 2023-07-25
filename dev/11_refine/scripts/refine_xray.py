import sys
from pathlib import Path
import shutil

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import params
import trackers
import log_statistics
import pdb_writer
import xray_restraint


if __name__ == "__main__":
    pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/51_synth_sa_5/9485465/output_31/pdbs/1195.pdb")
    pdb_file_ref = Path("/wynton/home/sali/mhancock/xray/data/pdbs/7mhf/7mhf_refine.pdb")
    n_steps = 100
    m = IMP.Model()
    hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    h = hs[0]
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    print(len(pids))
    ps = [m.get_particle(pid) for pid in pids]

    out_dir = Path(Path.home(), "xray/dev/11_refine/data/1_state_xray/output_0")
    shutil.rmtree(out_dir, ignore_errors=True)
    out_dir.mkdir(exist_ok=True, parents=True)
    pdb_dir = Path(out_dir, "pdbs")
    pdb_dir.mkdir()

    params_dict = dict()
    params_dict["pdb_file"] = pdb_file
    params_dict["pdb_id"] = 0
    params_dict["pdb_file_ref"] = pdb_file_ref
    params_dict["n_steps"] = n_steps
    params.write_params(
        param_dict=params_dict,
        param_file=str(Path(out_dir, "params.txt"))
    )

    rs = list()
    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmm.charmm_restraints(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)
    rset_charmm.set_weight(1)
    rs.append(rset_charmm)

    cif_file = Path(Path.home(), "xray/data/reflections/7mhf/7mhf_refine.cif")
    d_min = 0
    dyn_w_xray = 0
    w_xray = 100000

    # cif file here is a string.
    r_xray = xray_restraint.XtalRestraint(
        m=m,
        n_state=1,
        pids=pids,
        f_obs_file=cif_file,
        d_min=d_min,
        d_max=None,
        scale=True,
        target="ml",
        w_xray=w_xray,
        dynamic_w=dyn_w_xray
    )

    rs_xray = IMP.RestraintSet(m, 1.0)
    rs_xray.add_restraint(r_xray)
    rs.append(rs_xray)

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

    m_0 = IMP.Model()
    h_0 = IMP.atom.read_pdb(str(pdb_file_ref), m_0, IMP.atom.AllPDBSelector())
    rmsd_tracker = trackers.RMSDTracker(
        name="rmsd",
        hs=[h],
        hs_0=[h_0],
        align=False
    )
    all_trackers.append(rmsd_tracker)

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

    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=str(Path(out_dir, "log.csv")),
        log_freq=1
    )

    pdb_ostate = pdb_writer.WriteMultiStatePDBOptimizerState(
            m=m,
            hs=[h],
            pdb_dir=pdb_dir
    )
    pdb_ostate.set_period(1)
    pdb_ostate.do_update(None)

    sf = IMP.core.RestraintsScoringFunction(rs)
    cg = IMP.core.ConjugateGradients(m)
    cg.set_scoring_function(sf)
    cg.add_optimizer_state(log_ostate)
    cg.add_optimizer_state(pdb_ostate)
    cg.optimize(n_steps)

