from pathlib import Path

import IMP
import IMP.atom


class WriteMultiStatePDBOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            hs,
            pdb_dir
    ):
        IMP.OptimizerState.__init__(self, m, "LogStatistics%1%")
        self.hs = hs
        self.pdb_dir = pdb_dir
        self.cur_pdb_id = 0
        self.cur_pdb_file = Path(self.pdb_dir, "{}.pdb".format(self.cur_pdb_id))

    def do_update(self, call):
        IMP.atom.write_multimodel_pdb(self.hs, str(self.cur_pdb_file))

        self.cur_pdb_id = self.cur_pdb_id+1
        self.cur_pdb_file = Path(self.pdb_dir, "{}.pdb".format(self.cur_pdb_id))

