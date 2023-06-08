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
import pdb_optimizer_state


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/7mhf/7mhf_no_H20_alt_H_ion.pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    ps = [m.get_particle(pid) for pid in pids]

    out_dir = Path(Path.home(), "xray/dev/11_refine/data/output_0")
    shutil.rmtree(out_dir)
    out_dir.mkdir(exist_ok=True)
    pdb_dir = Path(out_dir, "pdbs")
    pdb_dir.mkdir()

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
    h_0 = IMP.atom.read_pdb(str(pdb_file), m_0, IMP.atom.AllPDBSelector())
    rmsd_tracker = trackers.RMSDTracker(
        name="rmsd",
        hs=[h],
        hs_0=[h_0],
        align=False
    )
    all_trackers.append(rmsd_tracker)

    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=all_trackers,
        log_file=str(Path(out_dir, "log.csv")),
        log_freq=1
    )

    pdb_ostate = pdb_optimizer_state.WriteMultiStatePDBOptimizerState(
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
    cg.optimize(100)

