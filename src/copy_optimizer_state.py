import shutil
from pathlib import Path
import time

import IMP
import IMP.core
import IMP.atom


class CopyOptimizerState(IMP.OptimizerState):
    def __init__(
            self,
            m,
            source_dir,
            dest_dir
    ):
        IMP.OptimizerState.__init__(self, m, "CopyOptimizerState%1%")
        self.m = m
        self.source_dir = source_dir
        self.dest_dir = dest_dir

    def do_update(self, call):
        for file in self.source_dir.glob("*"):
            dest_file = Path(self.dest_dir, file.name)
            shutil.move(file, dest_file)