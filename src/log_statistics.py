import time
import pandas as pd

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

        log = dict()
        self.log = log
        for tracker in self.all_trackers:
            self.add_tracker(
                tracker=tracker
            )

        self.print_header()

    def add_tracker(
            self,
            tracker
    ):
        if tracker.get_n() > 1:
            for i in range(tracker.get_n()):
                self.log[(tracker.get_name()+str(i))] = list()
        else:
            self.log[(tracker.get_name())] = list()

    def print_header(self):
        log_line = ""
        for entry in self.log.keys():
            log_line = log_line + entry + ","

        print(log_line)

    def print_last_entry(self):
        log_line = ""
        for entry in self.log.keys():
            log_line = log_line + str(round(self.log[entry][-1], 3)) + ","
        print(log_line)

    def get_log(
            self
    ):
        log_df = pd.DataFrame(self.log)
        return log_df

    def do_update(self, call):
        for tracker in self.all_trackers:
            result = tracker.evaluate()

            if tracker.get_n() > 1:
                for i in range(tracker.get_n()):
                    entry = tracker.get_name() + str(i)
                    self.log[entry].append(result[i])
            else:
                self.log[tracker.get_name()].append(result)

        self.print_last_entry()

        if self.get_number_of_updates()%self.log_freq == 0:
            t0 = time.time()
            log_df = self.get_log()
            log_df.to_csv(self.log_file)
            # print("Log write took: {}s".format(time.time()-t0))