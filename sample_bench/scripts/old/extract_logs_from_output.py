from pathlib import Path
import pandas as pd


if __name__ == "__main__":
    output_file = Path(Path.home(), "Documents/xray/dev/plot_time/md_300.o395440.1")

    with open(output_file) as f:
        lines = f.readlines()

    i = 0
    for line in lines:
        print(line)

        if "step," in line:
            break

        i = i+1

    print(i)

    log_df = pd.read_csv(output_file, skiprows=list(range(i)))
    log_df = log_df.drop(columns=["Unnamed: 5"])
    # print(log_df.head())
    # print(log_df.tail())

    log_df.to_csv(Path(Path.home(), "Documents/xray/dev/plot_time/log_395440.1.csv"))