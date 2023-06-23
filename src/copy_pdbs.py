import shutil
from pathlib import Path
import time

import IMP
import IMP.core
import IMP.atom

import trackers


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


class PDBCopyTracker(trackers.Tracker):
    def __init__(
            self,
            name,
            m,
            source_dir,
            dest_dir,
            freq=100
    ):
        trackers.Tracker.__init__(
            self,
            name=name,
            m=m,
            n=1
        )
        self.source_dir = source_dir
        self.dest_dir = dest_dir
        self.step = 0
        self.freq = freq

    def evaluate(
            self
    ):
        if self.step % self.freq == 0:
            for file in self.source_dir.glob("*.pdb"):
                dest_file = Path(self.dest_dir, file.name)
                shutil.move(file, dest_file)
            copied = 1
        else:
            copied = 0

        self.step = self.step+1

        return copied
