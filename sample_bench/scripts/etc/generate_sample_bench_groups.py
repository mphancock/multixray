import random
import multiprocessing
import pickle
from pathlib import Path


def get_m_rand_groups_of_size_n(
    params_dict
):
    elements = params_dict["elements"]
    m = params_dict["m"]
    n = params_dict["n"]
    groups_dir = params_dict["groups_dir"]

    groups = list()

    for i in range(m):
        elements_copy = elements.copy()
        group = list()
        for j in range(n):
            rand_element = random.choice(elements_copy)
            elements_copy.remove(rand_element)
            group.append(rand_element)

        groups.append(group)

        groups_file = Path(groups_dir, "{}.pl".format(n))
        with open(groups_file, "wb") as f:
            pickle.dump(groups, f)

    return n


if __name__ == "__main__":
    groups_dir = Path(Path.home(), "xray/sample_bench/data/groups/m_n_1000")

    pool_params = list()
    # for i in range(1,1000):
    for i in [1000]:
        params_dict = dict()
        params_dict["elements"] = list(range(1000))
        params_dict["m"] = 1000
        params_dict["n"] = i
        params_dict["groups_dir"] = groups_dir

        pool_params.append(params_dict)

        # break

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        get_m_rand_groups_of_size_n,
        pool_params
    )

    for pool_result in pool_results:
        n = pool_result
        print(n)

