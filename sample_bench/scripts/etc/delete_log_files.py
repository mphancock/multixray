import os
import datetime
from pathlib import Path
import multiprocessing
import random


def get_creation_time(file_path):
    creation_time_epoch = os.path.getctime(file_path)

    # Convert the creation time to a datetime object
    creation_time = datetime.datetime.fromtimestamp(creation_time_epoch)

    return creation_time


def delete_file(log_file):
    """
    Determine if the file should be deleted based on its creation time.
    Returns a tuple of the file path and a boolean indicating if it should be deleted.
    """
    year, month, day = 2025, 2, 2

    creation_time = get_creation_time(log_file)
    delete_before_date = datetime.datetime(year, month, day)
    ## YYYY, MM, DD

    if creation_time < delete_before_date:
        ## 1/1000 chance to print the file path and creation time
        if random.random() < .001:
            print(log_file, creation_time)

        log_file.unlink()
        return 1  # Increment count for deleted file
    else:
        return 0  # No deletion

if __name__ == "__main__":
    # Example usage:
    data_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/tmp")
    log_files = list(data_dir.glob("*"))
    # log_files = log_files[:10]
    # print(log_files)

    print(len(log_files))

    n_deleted = 0
    pool_obj = multiprocessing.Pool(multiprocessing.cpu_count())
    # Determine which files to delete in parallel
    result = pool_obj.imap(delete_file, log_files)
    for deleted in result:
        n_deleted += deleted

    print(n_deleted)
