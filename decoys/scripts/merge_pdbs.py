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
    hs_all = list()
    ms_all = list()

    pdb_files_order = list()

    if order:
        pdb_ids = [int(pdb_file.stem) for pdb_file in pdb_files]
        pdb_min = np.min(pdb_ids)
        pdb_max = np.max(pdb_ids)

        if n < 0 or n > pdb_max:
            max_it = pdb_max
        else:
            max_it = n

        pdb_dir = pdb_files[0].parents[0]
        for i in range(pdb_min, max_it):
            pdb_files_order.append(Path(pdb_dir, "{}.pdb".format(i)))

    else:
        pdb_files_order = pdb_files

    for i in range(len(pdb_files_order)):
        pdb_file = pdb_files_order[i]
        print(pdb_file)
        m = IMP.Model()

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


if __name__ == "__main__":

    pdb_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/39_32/2504728/output_0/pdbs")
    pdb_files = list(pdb_dir.glob("*.pdb"))
    out_file = Path(Path.home(), "xray/tmp/0.pdb")
    write_merge_pdb_file(
        merge_pdb_file=out_file,
        pdb_files=pdb_files,
        occs=[1]*len(pdb_files),
        n=-1,
        order=True
    )
