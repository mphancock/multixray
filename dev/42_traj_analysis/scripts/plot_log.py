import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_theme()
import matplotlib.colors as mcolors

if __name__ == "__main__":
    log_dfs = list()
    exp_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/247_res")
    exp_num = exp_dir.stem.split("_")[0]
    # for job_dir in Path(exp_dir, "logs").glob("*"):

    for start in [12, 36, 60]:
        log_dfs = list()
        job_ids = list(range(start,start+12))
        xray_name = "3k0n"

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


        all_fields = [["ff"], ["r_free_{}".format(xray_name)], ["r_work_{}".format(xray_name)], ["dcharmm_mag"], ["dxray_{}_mag".format(xray_name)], ["temp"], ["vel_mag"], ["wxray"], ["rmsd_0"], ["com_delta_mag"]]

        colors = ["tab:blue", "tab:orange", "tab:red"]

        fig, axs = plt.subplots(len(all_fields),len(log_dfs), figsize=(len(log_dfs)*10, 5*len(all_fields)))
        start, end, offset = 0, 100000, 1

        for i in range(len(all_fields)):
            for job_id in range(len(log_dfs)):
                if len(log_dfs[job_id]) == 0:
                    continue

                fields = all_fields[i]

                ax = axs[i][job_id]
                job_log_dfs = log_dfs[job_id]

                lns = list()

                for j in range(len(fields)):
                    field = fields[j]

                    for log_df in job_log_dfs:
                        if field not in log_df.columns:
                            continue

                        ln = ax.plot(log_df["step"][start:end:offset], log_df[fields[j]][start:end:offset], c=colors[j], label=fields[j], alpha=0.5)

                    lns.extend(ln)

                labs = [l.get_label() for l in lns]
                ax.legend(lns, labs)

        plot_dir = Path("/wynton/home/sali/mhancock/xray/dev/42_traj_analysis/data/{}/plots".format(exp_num))
        plot_dir.mkdir(exist_ok=True, parents=True)

        plt.savefig(Path(plot_dir, "{}.png".format(job_ids[0])))