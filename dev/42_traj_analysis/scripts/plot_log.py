import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
import matplotlib.colors as mcolors

if __name__ == "__main__":
    log_dfs = list()
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/263_sb_sa")
    exp_num = exp_dir.stem.split("_")[0]
    # for job_dir in Path(exp_dir, "logs").glob("*"):

    for start,end,xray_name in [(0,5, "native_0"),(6,11, "native_0"),(12,17, "native_0"),(18,23, "native_0"),(24,29, "native_0"),(30,35, "native_0")]:
        log_dfs = list()
        job_ids = list(range(start, end+1))

        for job_id in job_ids:
        # for job_id in job_ids:
            job_dir = Path(exp_dir, str(job_id))
            job_log_dfs = list()
            # job_dir = Path(exp_dir, "logs/{}".format(job_id))

            print(job_dir)

            steps = list()

            for output_dir in job_dir.glob("output*"):
                log_file = Path(output_dir, "log.csv")

                if not log_file.exists():
                    continue

                log_df = pd.read_csv(log_file, index_col=0)
                job_log_dfs.append(log_df)
                steps.append(log_df["step"].iloc[-1])

            log_dfs.append(job_log_dfs)

            print(steps)


        print(len(log_dfs))

        colors = ["tab:blue", "tab:orange", "tab:red"]

        fields = ["ff", "r_free_{}".format(xray_name), "r_work_{}".format(xray_name), "dcharmm_mag", "dxray_{}_mag".format(xray_name), "temp", "vel_mag", "wxray", "rmsd_0", "com_delta_mag", "w_0_{}".format(xray_name)]

        fig, axs = plt.subplots(len(fields), len(log_dfs), figsize=(len(log_dfs)*10, 5*len(fields)))
        start, end, offset = 0, 100000, 1

        for job_id in range(len(log_dfs)):
            if len(log_dfs[job_id]) == 0:
                continue

            for i in range(len(fields)):
                ax = axs[i][job_id]
                job_log_dfs = log_dfs[job_id]

                lns = list()

                # fields = fields[i]
                field = fields[i]

                all_ys = list()
                for log_df in job_log_dfs:
                    if field not in log_df.columns:
                        continue

                    ln = ax.plot(log_df["step"][start:end:offset], log_df[field][start:end:offset], c=colors[0], label=field, alpha=0.5)

                    all_ys.extend(log_df[field][start:end:offset])

                all_ys = [x for x in all_ys if str(x) != 'nan']


                print(len(all_ys))

                if len(all_ys) == 0:
                    continue

                print(field, min(all_ys), max(all_ys))
                # ax.set_ylim(, np.percentile(all_ys, 95)*2)

                lns.extend(ln)

                labs = [l.get_label() for l in lns]
                ax.legend(lns, labs)

        plot_dir = Path("/wynton/home/sali/mhancock/xray/dev/42_traj_analysis/data/{}/plots".format(exp_num))
        plot_dir.mkdir(exist_ok=True, parents=True)

        plt.savefig(Path(plot_dir, "{}.png".format(job_ids[0])))