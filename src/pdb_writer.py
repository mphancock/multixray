from pathlib import Path
import numpy as np

import IMP
import IMP.atom

import trackers

class WriteMultiStatePDBOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            hs,
            pdb_dir
    ):
        IMP.OptimizerState.__init__(self, m, "PDBWriter%1%")
        self.hs = hs
        self.pdb_dir = pdb_dir
        self.cur_pdb_id = 0
        self.cur_pdb_file = Path(self.pdb_dir, "{}.pdb".format(self.cur_pdb_id))

    def do_update(self, call):
        IMP.atom.write_multimodel_pdb(self.hs, str(self.cur_pdb_file))

        self.cur_pdb_id = self.cur_pdb_id+1
        self.cur_pdb_file = Path(self.pdb_dir, "{}.pdb".format(self.cur_pdb_id))


class WriteBestMultiStatePDBOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            hs,
            pdb_dir,
            tracker,
            N,
            n_skip
    ):
        IMP.OptimizerState.__init__(self, m, "PDBWriter%1%")
        self.hs = hs
        self.pdb_dir = pdb_dir
        self.tracker = tracker
        self.cache = [np.infty]*N
        self.n_skip = n_skip

        self.step = 0

    def do_update(self, call):
        score = self.tracker.evaluate()
        if score < np.max(self.cache) and self.step > self.n_skip:
            print("Updating cache")
            max_id = np.argmax(self.cache)
            self.cache[max_id] = score
            print(self.cache)

            pdb_file = Path(self.pdb_dir, "{}.pdb".format(max_id))
            IMP.atom.write_multimodel_pdb(self.hs, str(pdb_file))

        self.step = self.step+1


class PDBWriterTracker(trackers.Tracker):
    def __init__(
            self,
            name,
            hs,
            pdb_dir,
            freq=10
    ):
        trackers.Tracker.__init__(
            self,
            name=name,
            m=hs[0].get_model(),
            n=1
        )
        self.hs = hs
        self.pdb_dir = pdb_dir
        self.step = 0
        self.freq = freq

        self.cur_pdb_id = 0

    def evaluate(
            self
    ):
        if self.step % self.freq == 0 and self.writing:
            cur_pdb_file = Path(self.pdb_dir, "{}.pdb".format(self.cur_pdb_id))
            IMP.atom.write_multimodel_pdb(self.hs, str(cur_pdb_file))
            self.cur_pdb_id = self.cur_pdb_id+1
            pdb_name = cur_pdb_file.name
        else:
            pdb_name = np.nan

        self.step = self.step+1

        return pdb_name