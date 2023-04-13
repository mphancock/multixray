from pathlib import Path
import multiprocessing

import IMP
import IMP.atom


def pool_test(
        params_dict
):
    h = params_dict["h"]
    print(len(IMP.atom.Selection(h).get_selected_particle_indexes()))


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())

    pool_params = list()
    pool_params.append({"h": h})

    print("CPUs: {}".format(multiprocessing.cpu_count()))
    pool_obj = multiprocessing.Pool(
        multiprocessing.cpu_count()
    )

    pool_results = pool_obj.imap(
        pool_test,
        pool_params
    )

    for pool_result in pool_results:
        continue



