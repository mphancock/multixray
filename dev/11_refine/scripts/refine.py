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


if __name__ == "__main__":
    pdb_files = list(Path(Path.home(), "xray/dev/17_synthetic_native/data/pdbs/1_state").glob("*.pdb"))

    for pdb_file in pdb_files:
        ref_pdb_file = Path(Path.home(), "xray/dev/17_synthetic_native/data/pdbs/1_state_ref/{}.pdb".format(pdb_file.stem))
        # pdb_file = Path("/wynton/group/sali/mhancock/xray/decoys/data/3ca7/53_100/rand_1000_2x/67.pdb")
        pdb_file_ref = Path("/wynton/home/sali/mhancock/xray/data/pdbs/3ca7/3ca7_refine.pdb")
        n_steps = 100
        m = IMP.Model()
        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

        # job_dir = Path(Path.home(), "xray/dev/11_refine/data/3ca7_4_state")
        # job_dir.mkdir(exist_ok=True)

        for i in range(len(hs)):
            h = hs[i]
            pids = IMP.atom.Selection(h).get_selected_particle_indexes()
            print(len(pids))
            ps = [m.get_particle(pid) for pid in pids]

            # out_dir = Path(job_dir, "output_{}".format(i))
            # shutil.rmtree(out_dir, ignore_errors=True)
            # out_dir.mkdir(exist_ok=True)
            # pdb_dir = Path(out_dir, "pdbs")
            # pdb_dir.mkdir()

            params_dict = dict()
            params_dict["pdb_file"] = pdb_file
            params_dict["pdb_id"] = i
            params_dict["pdb_file_ref"] = pdb_file_ref
            params_dict["n_steps"] = n_steps

            # params.write_params(
            #     param_dict=params_dict,
            #     param_file=str(Path(out_dir, "params.txt"))
            # )

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

            o_states = list()

            # m_0 = IMP.Model()
            # h_0 = IMP.atom.read_pdb(str(pdb_file_ref), m_0, IMP.atom.AllPDBSelector())
            # rmsd_tracker = trackers.RMSDTracker(
            #     name="rmsd",
            #     hs=[h],
            #     hs_0=[h_0],
            #     ca_only=True
            # )
            # all_trackers.append(rmsd_tracker)

            # log_ostate = log_statistics.LogStatistics(
            #     m=m,
            #     all_trackers=all_trackers,
            #     log_file=str(Path(out_dir, "log.csv")),
            #     log_freq=1
            # )
            # o_states.append(log_ostate)
            # pdb_ostate = pdb_writer.WriteMultiStatePDBOptimizerState(
            #         m=m,
            #         hs=[h],
            #         pdb_dir=pdb_dir
            # )
            # pdb_ostate.set_period(1)
            # pdb_ostate.do_update(None)
            # o_states.append(pdb_ostate)

            sf = IMP.core.RestraintsScoringFunction(rs)
            cg = IMP.core.ConjugateGradients(m)
            cg.set_scoring_function(sf)

            for o_state in o_states:
                cg.add_optimizer_state(o_state)

            cg.optimize(n_steps)

        IMP.atom.write_multimodel_pdb(hs, str(ref_pdb_file))

