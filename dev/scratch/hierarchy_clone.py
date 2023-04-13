from pathlib import Path
import IMP
import IMP.core
import IMP.atom


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/dev/scratch/one_res.pdb")
    m = IMP.Model()
    h = IMP.atom.read_pdb(str(pdb_file), m)
    pids = IMP.atom.Selection(h).get_selected_particle_indexes()
    xyz = IMP.core.XYZ(m, pids[0])
    print(xyz.get_coordinate(0))
    print(len(m.get_particle_indexes()))

    h_clone = IMP.atom.create_clone(h)
    pids_clone = IMP.atom.Selection(h_clone).get_selected_particle_indexes()
    print(len(m.get_particle_indexes()))

    xyz_clone = IMP.core.XYZ(m, pids_clone[0])
    xyz_clone.set_coordinate(0,19.3)
    print(xyz_clone.get_coordinate(0))
    print(xyz.get_coordinate(0))

    print(h.get_model())
    print(h_clone.get_model())
    print(h.get_model() == h_clone.get_model())
