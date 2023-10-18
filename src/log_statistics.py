import time
import pandas as pd
import numpy as np

import IMP
import IMP.core
import IMP.atom


class LogStatistics(IMP.OptimizerState):
    def __init__(
            self,
            m,
            all_trackers,
            log_file,
            log_freq
    ):
        IMP.OptimizerState.__init__(self, m, "LogStatistics%1%")
        self.m = m
        self.t_0 = time.time()
        self.all_trackers = all_trackers
        self.log_file = log_file
        self.log_freq = log_freq

        columns = list()
        for tracker in self.all_trackers:
            if tracker.get_n() > 1:
                for i in range(tracker.get_n()):
                    columns.append(tracker.get_name()+"_"+str(i))
            else:
                columns.append(tracker.get_name())

        self.log_df = pd.DataFrame(columns=columns)
        self.print_header()


    def print_header(self):
        log_line = ""
        for entry in self.log_df.columns:
            log_line = log_line + entry + ","

        print(log_line)

    def print_last_entry(self):
        log_line = ""
        for column in self.log_df.columns:
            val = self.log_df.iloc[-1][column]
            if type(val) == float or type(val) == np.float64:
                last_entry = str(round(val, 3))
            else:
                last_entry = str(val)
            log_line = log_line + last_entry  + ","

        print(log_line)

    def get_log(
            self
    ):
        return self.log_df

    def get_tracker(
            self,
            name
    ):
        for tracker in self.all_trackers:
            if tracker.get_name() == name:
                return tracker

        raise RuntimeError("Tracker {} not found".format(name))

    def get_trackers(
            self
    ):
        return self.all_trackers

    def do_update(self, call):
        entries = list()
        for tracker in self.all_trackers:
            result = tracker.evaluate()
            if tracker.get_n() > 1:
                for i in range(tracker.get_n()):
                    entries.append(result[i])
            else:
                entries.append(result)

        self.log_df.loc[len(self.log_df)] = entries
        self.print_last_entry()

        if self.get_number_of_updates()%self.log_freq == 0:
            log_df = self.get_log()
            log_df.to_csv(self.log_file)