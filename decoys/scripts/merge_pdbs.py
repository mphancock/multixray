import IMP
import IMP.atom
from pathlib import Path
import numpy as np


def write_merge_pdf_file_pool(
        param_dict
):
    write_merge_pdb_file(
        merge_pdb_file=param_dict["merge_pdb_file"],
        pdb_files=param_dict["pdb_files"],
        occs=param_dict["occs"],
        n=param_dict["n"],
        order=param_dict["order"]
    )

    return param_dict["id"]


def write_merge_pdb_file(
        merge_pdb_file,
        pdb_files,
        occs,
        n,
        order
):
    # print(pdb_files)
    # pdb_dir = pdb_files[0].parents[0]
    # print(pdb_dir, len(pdb_files))

    hs_all = list()
    ms_all = list()
    # Need to iterate through by index in order to ensure that MD frames are ordered.

    if order:
        pdb_ids = [int(pdb_file.stem) for pdb_file in pdb_files]
        pdb_min = np.min(pdb_ids)
        pdb_max = np.max(pdb_ids)

        if n < 0 or n > pdb_max:
            max_it = pdb_max
        else:
            max_it = n

    for i in range(len(pdb_files)):
        pdb_file = pdb_files[i]
        print(pdb_file)
    # for i in range(pdb_min, max_it+1):
    #     pdb_file = Path(pdb_dir, "{}.pdb".format(i))
        m = IMP.Model()
        # h = IMP.atom.read_pdb(str(pdb_file), m)

        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m)

        for h in hs:
            pids = IMP.atom.Selection(h).get_selected_particle_indexes()
            for pid in pids:
                at = IMP.atom.Atom(m, pid)
                at.set_occupancy(occs[i])
            hs_all.append(h)

        ms_all.append(m)

    # Cannot return hs, because the hierarchy objects will get deallocated.
    IMP.atom.write_multimodel_pdb(hs_all, str(merge_pdb_file))


def write_merge_pdb_file_from_multistate_output(
        merge_pdb_file,
        pdb_files,
        occs,
        n,
        order
):
    return 0


if __name__ == "__main__":
    for n_state in [2,4,8,16,32]:
        pdb_files = [Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean.pdb")]*n_state
        out_file = Path(Path.home(), "xray/data/pdbs/3ca7/3ca7_clean_x{}.pdb".format(n_state))
        write_merge_pdb_file(
            merge_pdb_file=out_file,
            pdb_files=pdb_files,
            occs=[1]*n_state,
            n=-1,
            order=False
        )
