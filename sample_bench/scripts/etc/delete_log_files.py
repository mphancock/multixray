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


from multiprocessing import Pool
from pathlib import Path

def should_delete_file(log_file):
    """
    Determine if the file should be deleted based on its creation time.
    Returns a tuple of the file path and a boolean indicating if it should be deleted.
    """
    creation_time = get_creation_time(log_file)
    delete_file = file_created_before_date(log_file, 2024, 9, 22)
    print(f"{log_file} {creation_time} {delete_file}")
    return log_file, delete_file

def delete_files(file_info):
    """
    Deletes the file if it should be deleted.
    """
    log_file, delete_file = file_info
    if delete_file:
        log_file.unlink()
        return 1  # Increment count for deleted file
    return 0  # No deletion

if __name__ == "__main__":
    # Example usage:
    data_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/tmp")
    log_files = list(data_dir.glob("*"))

    print(len(log_files))

    with Pool() as pool:
        # Determine which files to delete in parallel
        files_to_delete = pool.map(should_delete_file, log_files)

        # Delete files in parallel
        deleted_counts = pool.map(delete_files, files_to_delete)

    n_deleted = sum(deleted_counts)
    print(n_deleted)

# if __name__ == "__main__":
#     # Example usage:
#     data_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/tmp")
#     log_files = list(data_dir.glob("*"))

#     print(len(log_files))
#     n_deleted = 0
#     for log_file in log_files:
#         creation_time = get_creation_time(log_file)
#         delete_file = file_created_before_date(log_file, 2024, 9, 3)
#         print("{} {} {}".format(log_file, creation_time, delete_file))
#         if delete_file:
#             log_file.unlink()
#             n_deleted += 1
#     print(n_deleted)

