from pathlib import Path
import sys
import random

import IMP
import IMP.core
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import charmm
import params
import trackers
import log_statistics
import pdb_writer
import molecular_dynamics


class PlaneRestraint(IMP.Restraint):
    # take the list of particles and the key to use

    def __init__(self, m, k, ps):
        IMP.Restraint.__init__(self, m, "MyRestraint %1%")
        self.ps = ps
        self.k = k

    def unprotected_evaluate(self, da):
        score = 0
        for i in range(1, len(self.ps)):
            d = IMP.core.XYZ(self.ps[i])

            score = score+self.k*(d.get_z())**2
            d.add_to_derivative(2, 2*self.k*d.get_z(), da)

        print(score)
        return score

    def do_get_inputs(self):
        return self.ps


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/dev/21_degeneracy/data/AA.pdb")

    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    pids = IMP.atom.Selection(h).get_selected_particle_indexes()

    rs = list()
    rset_charmm = IMP.RestraintSet(m, 1.0)
    charmm_rs = charmm.charmm_restraints(
        m,
        h,
        eps=False
    )
    rset_charmm.add_restraints(charmm_rs)
    rs.append(rset_charmm)

    out_dir = Path(Path.home(), "xray/dev/21_degeneracy/data/output_0")
    pdb_dir = Path(out_dir, "pdbs")
    pdb_dir.mkdir(exist_ok=True, parents=True)
    pdb_tracker = pdb_writer.PDBWriterTracker(
        name="pdb",
        hs=[h],
        pdb_dir=pdb_dir,
        freq=10,
        log_pdb_dir=None
    )

    o_states = list()
    log_ostate = log_statistics.LogStatistics(
        m=m,
        all_trackers=[pdb_tracker],
        log_file=Path(out_dir, "log.csv"),
        log_freq=100
    )
    for r in rs:
        r.evaluate(calc_derivs=True)
    log_ostate.update()
    o_states.append(log_ostate)

    # md_sel = IMP.atom.Selection(h) - IMP.atom.Selection(h, residue_index=72, atom_type=IMP.atom.AtomType("N"))
    md_sel = IMP.atom.Selection(h, atom_types=[IMP.atom.AtomType("CA"), IMP.atom.AtomType("CB")])
    md_ps = [m.get_particle(pid) for pid in md_sel.get_selected_particle_indexes()]

    molecular_dynamics.molecular_dynamics(
        output_dir=Path(Path.home(), "xray/dev/21_degeneracy/data/output_0"),
        hs=[h],
        rs=rs,
        T=3000,
        t_step=2,
        steps=1000,
        sa_sched=None,
        o_states=o_states,
        md_ps=md_ps
    )



