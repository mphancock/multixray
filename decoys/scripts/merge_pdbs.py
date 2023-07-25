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
        order=param_dict["order"],
        state=param_dict["state"]
    )

    return param_dict["id"]


"""
This function is used to merge multiple pdb files together in a multi-model pdb file. This is useful in creating MD trajectory movies from sampling or merging together multiple 1-state structures into a multi-state structure.

Params:
    merge_pdb_file: Path to the output multi-model pdb file.

    pdb_files: List of paths to the pdb files to merge together.

    occs: List of occupancy values for each pdb file. The length of this list must be equal to the length of pdb_files.

    n: Number of pdb files to merge together. If n < 0, then all pdb files will be merged together. If n > len(pdb_files), then all pdb files will be merged together.

    order: If True, then the pdb files will be merged together in order of their file names. If False, then the pdb files will be merged together in the order they are given in the pdb_files list.

    state: The state of the model to read from the pdb files.

Returns:
    None

"""
def write_merge_pdb_file(
        merge_pdb_file,
        pdb_files,
        occs,
        n,
        order,
        state
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

        hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        h = hs[state]

        pids = IMP.atom.Selection(h).get_selected_particle_indexes()

        if occs:
            for pid in pids:
                at = IMP.atom.Atom(m, pid)
                at.set_occupancy(occs[i])

        hs_all.append(h)
        ms_all.append(m)

    # Cannot return hs, because the hierarchy objects will get deallocated.
    IMP.atom.write_multimodel_pdb(hs_all, str(merge_pdb_file))


if __name__ == "__main__":
    pdb_dir = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/3ca7/63_native_2x/39/output_0/pdbs")
    pdb_files = list(pdb_dir.glob("*.pdb"))
    print(pdb_files)

    out_file = Path(Path.home(), "xray/tmp/merge_0.pdb")
    write_merge_pdb_file(
        merge_pdb_file=out_file,
        pdb_files=pdb_files,
        occs=None ,
        n=500,
        order=True,
        state=0
    )
