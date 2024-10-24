import shutil
from pathlib import Path
import time
import numpy as np

import IMP
import IMP.atom

import trackers
import align_imp


"""
This tracker writes a PDB file every freq steps.

**********
Parameters

    name: the name of the tracker.

    hs: the IMP hierarchies to write to PDB.

    pdb_dir: the directory to write the PDB files to.

    freq: the frequency to write PDB files.

    log_pdb_dir: the directory to record the PDB files in the log. This is necessary because the PDB files are written locally to a temp directory but then moved to the permanent pdb directory.

**********
Returns

    return_val: the path to the PDB file written.

"""
class PDBWriterTracker(trackers.Tracker):
    def __init__(
            self,
            name,
            msmc_m,
            pdb_dir,
            log_pdb_dir=None
    ):
        trackers.Tracker.__init__(
            self,
            name=name,
            m=msmc_m.get_m(),
            n=1
        )
        self.msmc_m = msmc_m
        self.pdb_dir = pdb_dir
        self.step = 0
        self.log_pdb_dir = log_pdb_dir
        self.cur_pdb_id = 0

    def do_evaluate(self):
        cur_pdb_file = Path(self.pdb_dir, "{}.pdb".format(self.cur_pdb_id))

        # print(cur_pdb_file)
        # IMP.atom.write_multimodel_pdb(self.hs, str(cur_pdb_file))
        self.msmc_m.write_pdb_file(cur_pdb_file)
        self.cur_pdb_id = self.cur_pdb_id+1

        # final_pdb_file is the pdb file after copying.
        if self.log_pdb_dir:
            final_pdb_file = Path(self.log_pdb_dir, cur_pdb_file.name)
        else:
            final_pdb_file = cur_pdb_file

        return final_pdb_file

    def evaluate(
            self
    ):
        if self.step % self.get_period() == 0 and self.get_on():
            final_pdb_file = self.do_evaluate()
        else:
            final_pdb_file = np.nan

        self.step = self.step+1

        return [final_pdb_file]


