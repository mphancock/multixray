import IMP
import IMP.atom


def get_n_state_from_pdb_file(pdb_file):
    pdb_return = pool_read_pdb(pdb_file)

    if isinstance(pdb_return, RuntimeError):
        raise pdb_return
    else:
        m, hs = pdb_return

    return len(hs)

def pool_read_pdb(
    pdb_file
):
    print(pdb_file)

    if pdb_file.exists():
        m = IMP.Model()

        try:
            hs = IMP.atom.read_multimodel_pdb(str(pdb_file), m, IMP.atom.AllPDBSelector())
        except IMP.ValueException:
            return RuntimeError("{} cannot be read".format(pdb_file))

        return m, hs
    else:
        return RuntimeError("{} file does not exist".format(pdb_file))


if __name__ == "__main__":
    from pathlib import Path

    pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/168_N4_ref/14833.pdb")
    print(pool_read_pdb(pdb_file))