import os
import datetime
from pathlib import Path


def file_created_before_date(file_path, year, month, day):
    # Get the creation time of the file
    creation_time_epoch = os.path.getctime(file_path)

    # Convert the creation time to a datetime object
    creation_time = datetime.datetime.fromtimestamp(creation_time_epoch)

    # Define the specific date to compare against
    specific_date = datetime.datetime(year, month, day)

    # Compare the creation time with the specific date
    return creation_time < specific_date


def get_creation_time(file_path):
    creation_time_epoch = os.path.getctime(file_path)

    # Convert the creation time to a datetime object
    creation_time = datetime.datetime.fromtimestamp(creation_time_epoch)

    return creation_time


if __name__ == "__main__":
    # Example usage:
    data_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/tmp")
    log_files = list(data_dir.glob("*"))

    print(len(log_files))
    for log_file in log_files:
        # print(log_file, get_creation_time(log_file))
        delete_file = file_created_before_date(log_file, 2024, 4, 23)
        print("{} {}".format(log_file, delete_file))
        if delete_file:
            log_file.unlink()

