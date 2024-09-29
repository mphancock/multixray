import sys
from pathlib import Path

import IMP
import IMP.atom

import mmtbx.f_model
import mmtbx.model
import cctbx.crystal
import cctbx.xray
from iotbx import pdb
from iotbx.pdb.hierarchy import root
from scitbx.array_family import flex

sys.path.append(str(Path(Path.home(), "xray/src")))
import utility
from u_aniso import get_u_anisos_from_file



## we need to read in only non altloc atoms to get u aniso becasue that is what IMP will be doing
def read_non_altloc_structure(
    pdb_file
):
    pdb_inp = pdb.input(file_name=str(pdb_file))
    hierarchy = pdb_inp.construct_hierarchy()

    asc = hierarchy.atom_selection_cache()

    # Create a selection for atoms where altloc is not 'B'
    sel = asc.selection("not altloc 'B'")

    # Alternatively, if you want to select atoms where altloc is blank or any character except 'B'
    # sel = asc.selection("altloc '?' and not altloc 'B'")

    # Apply the selection to get a new hierarchy with only the selected atoms
    hierarchy_not_altloc_B = hierarchy.select(sel)

    atoms = hierarchy_not_altloc_B.atoms()

    crystal_symmetry = pdb_inp.crystal_symmetry()

    # xray_structure = hierarchy_not_altloc_B.extract_xray_structure(crystal_symmetry=crystal_symmetry)

    return hierarchy_not_altloc_B, crystal_symmetry


class MultiStateMultiConditionModel:
    def __init__(
            self,
            pdb_file,
            w_mat
    ):
        self.m = IMP.Model()
        self.set_w_mat(w_mat)

        try:
            pdb_file_n_state = utility.get_n_state_from_pdb_file(pdb_file)
        except RuntimeError as e:
            raise e

        self.hs = list()

        # pdb_selectors = [IMP.atom.ChainPDBSelector(["A"]), IMP.atom.NonWaterNonHydrogenPDBSelector(), IMP.atom.NonAlternativePDBSelector(), IMP.atom.ATOMPDBSelector()]
        pdb_selectors = [IMP.atom.NonAlternativePDBSelector()]
        sel = pdb_selectors[0]
        for i in range(1, len(pdb_selectors)):
            sel = IMP.atom.AndPDBSelector(sel, pdb_selectors[i])

        ## create a temp list to store xray structures
        self.xray_hs = list()

        if pdb_file_n_state == 1 and self.n_state >= 1:
            for i in range(self.n_state):
                self.hs.append(IMP.atom.read_pdb(str(pdb_file), self.m, sel))
                hierarchy, crystal_symmetry = read_non_altloc_structure(pdb_file)

                self.xray_hs.append(hierarchy)
                self.crystal_symmetry = crystal_symmetry
        # elif pdb_file_n_state == self.n_state:
        #     self.hs.extend(IMP.atom.read_multimodel_pdb(str(pdb_file), self.m, sel))
        else:
            raise RuntimeError("Number of states in pdb file ({}) does not match the number of states in the model ({}).".format(pdb_file_n_state, self.n_state))

        ## merge the xray hierarchies as well
        self.merge_xray_h = root()

        for i, xray_hierarchy in enumerate(self.xray_hs):
            model = xray_hierarchy.only_model()
            model.id = str(i + 1)
            self.merge_xray_h.append_model(model.detached_copy())

        self.all_pids = list()
        self.all_atoms = list()
        self.state_atoms = dict()

        for i in range(self.n_state):
            pids = IMP.atom.Selection(self.hs[i]).get_selected_particle_indexes()
            self.all_pids.extend(pids)

            atoms = self.xray_hs[i].atoms()
            self.all_atoms.extend(atoms)

            self.state_atoms[i] = atoms

        print("SETTING UP MODEL", pdb_file)
        print("NO STATES: ", len(self.hs))
        pids = IMP.atom.Selection(self.hs[0]).get_selected_particle_indexes()
        print("NO OF ATOMS PER STATE: ", len(pids))
        # i = 0
        # for pid in pids:
        #     print(i, IMP.atom.Atom(self.m, pid).get_name())
        #     i += 1

        self.pids_dict = dict()
        self.prot_pids_dict = dict()
        self.water_pids_dict = dict()
        self.backbone_pids_dict = dict()
        self.side_pids_dict = dict()
        self.ca_pids_dict = dict()
        water_at_type = IMP.atom.AtomType("HET: O  ")

        for i in range(self.n_state):
            h = self.hs[i]
            self.pids_dict[i] = IMP.atom.Selection(h).get_selected_particle_indexes()
            self.prot_pids_dict[i] = (IMP.atom.Selection(h) - IMP.atom.Selection(h, atom_type=water_at_type)).get_selected_particle_indexes()
            self.water_pids_dict[i] = IMP.atom.Selection(h, atom_type=water_at_type).get_selected_particle_indexes()
            self.backbone_pids_dict[i] = IMP.atom.Selection(h, atom_types=[IMP.atom.AT_CA, IMP.atom.AT_C, IMP.atom.AT_O, IMP.atom.AT_N]).get_selected_particle_indexes()
            self.side_pids_dict[i] = list(set(self.prot_pids_dict[i]) - set(self.backbone_pids_dict[i]))
            self.ca_pids_dict[i] = IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()

        ## imp setup
        # Setup waters
        for i in range(len(self.hs)):
            for pid in self.water_pids_dict[i]:
                IMP.atom.CHARMMAtom.setup_particle(self.m, pid, "O")

        # Setup coordinates for md
        for h in self.hs:
            for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
                d = IMP.core.XYZR(self.m, pid)
                d.set_coordinates_are_optimized(True)

                IMP.atom.LinearVelocity.setup_particle(self.m, pid)

        ## cctbx setup
        self.perform_checks()

        self.create_lookup_table()

        self.multi_xray_structure = self.merge_xray_h.extract_xray_structure(crystal_symmetry=self.crystal_symmetry)

        self.multi_xray_structure.scatterers().flags_set_grads(
            state=False
        )
        self.multi_xray_structure.scatterers().flags_set_grad_site(
            iselection=self.multi_xray_structure.all_selection().iselection()
        )
        self.multi_xray_structure.scatterers().flags_set_grad_occupancy(
            iselection=self.multi_xray_structure.all_selection().iselection()
        )

    def get_m(self):
        return self.m

    def get_hs(self):
        return self.hs

    def get_n_state(self):
        return self.n_state

    def get_w_mat(self):
        return self.w_mat

    def get_u_anisos(self):
        return self.u_anisos

    def get_all_pids(self):
        return self.all_pids

    def get_pids_in_state(self, i):
        return self.pids_dict[i]

    def get_ca_pids(self, i):
        return self.ca_pids_dict[i]

    def get_pids_in_res_range(self, start, end):
        pids = list()
        for h in self.hs:
            pids.extend(IMP.atom.Selection(h, residue_indexes=range(start, end)).get_selected_particle_indexes())

        return pids

    def get_all_ca_pids(self):
        ca_pids = list()
        for i in range(self.n_state):
            ca_pids.extend(self.ca_pids_dict[i])

        return ca_pids

    def get_all_side_pids(self):
        side_pids = list()
        for i in range(self.n_state):
            side_pids.extend(self.side_pids_dict[i])

        return side_pids

    def get_xray_hierarchy(self, i):
        return self.xray_hs[i]

    def get_all_atoms(self):
        return self.all_atoms

    def get_atoms_in_state(self, i):
        return self.state_atoms[i]

    def get_cryustal_symmetry(self):
        return self.crystal_symmetry

    def get_atom(self, pid):
        return self.pid_to_atom[pid]

    def get_pid(self, atom):
        return self.atom_to_pid[atom]

    def get_com(self):
        ## ConterOfMass does not dynamically update so need to reinstatiant
        self.com = IMP.atom.CenterOfMass.setup_particle(
            IMP.Particle(self.m),
            self.get_all_ca_pids()
        )
        return self.com

    def get_occs_for_condition_i(self, i):
        return self.w_mat[:,i]

    def set_w_mat(self, w_mat):
        self.w_mat = w_mat
        self.normalize_w_mat()

        self.n_state = self.w_mat.shape[0]
        self.n_cond = self.w_mat.shape[1]

    def set_occs_for_condition_i(self, occs, i):
        self.w_mat[:,i] = occs
        self.normalize_w_mat()

    def normalize_w_mat(self):
        column_sums = self.w_mat.sum(axis=0)
        self.w_mat = self.w_mat / column_sums

    ## create a lookup table between IMP pid and cctbx atom
    def create_lookup_table(self):
        self.pid_to_atom = dict()
        self.atom_to_pid = dict()
        for i in range(self.n_state):
            h = self.hs[i]
            h_xray = self.xray_hs[i]

            for pid in self.get_pids_in_state(i):
                at = IMP.atom.Atom(self.m, pid)
                atom_type = at.get_atom_type()
                atom_type = str(atom_type).strip("\"")

                if atom_type == "HET: O  ":
                    atom_type = "O"

                res = IMP.atom.Residue(at.get_parent())
                res_id = res.get_index()

                asc = h_xray.atom_selection_cache()
                sel_str = "resseq {} and name {}".format(res_id, atom_type)
                sel = asc.selection(sel_str)
                sel_atoms = h_xray.select(sel).atoms()

                if len(sel_atoms) != 1:
                    raise RuntimeError("{} atoms found for {}".format(len(sel_atoms), sel_str))

                self.pid_to_atom[pid] = sel_atoms[0]
                self.atom_to_pid[sel_atoms[0]] = pid

    def update_xray_hs(self):
        for state in range(self.n_state):
            for pid in self.get_pids_in_state(state):
                atom = self.get_atom(pid)
                d = IMP.core.XYZR(self.m, pid)
                coords = (d.get_x(), d.get_y(), d.get_z())
                atom.set_xyz(coords)

        xyzs = list()
        for state in range(self.n_state):
            xyzs.extend(self.get_xray_hierarchy(state).atoms().extract_xyz())
        xyzs = flex.vec3_double(xyzs)
        self.merge_xray_h.atoms().set_xyz(xyzs)

    ## update the individual xray structures then get the multi xray structure which updates automatically
    def get_multi_xray_structure(self, cond):
        if cond == 0:
            ## this is the relatively expensive step
            self.update_xray_hs()

        ## update the coordinates of the multi xray structure from the merged xray hierarchy
        self.multi_xray_structure.set_sites_cart(self.merge_xray_h.atoms().extract_xyz())

        ## update the occs of the multi xray structure from the w_mat directly
        occs = list()
        for state in range(self.n_state):
            occs.extend([self.w_mat[state, cond]]*len(self.get_pids_in_state(state)))
        occs = flex.double(occs)

        self.multi_xray_structure.set_occupancies(occs)

        ## build the multi xray structure
        return self.multi_xray_structure

    def perform_checks(self):
        ## check for each state the number of pids and the number of atoms are equal
        for i in range(self.n_state):
            n_pids = len(self.get_pids_in_state(i))
            n_atoms = len(self.get_atoms_in_state(i))

            if n_pids != n_atoms:
                raise RuntimeError("Number of atoms in state {} xray structure ({}) does not match the number of atoms in the hierarchy ({})".format(i, n_pids, n_atoms))

    def write_pdb_file(self, pdb_file):
        # Write the PDB content to a file
        # new_sites_cart = self.multi_xray_structure.sites_cart()
        # Get the atoms from the hierarchy
        # atoms = self.multi_xray_hierarchy.atoms()
        # Update the atomic positions in the hierarchy
        # atoms.set_xyz(new_sites_cart)

        self.update_xray_hs()

        pdb_content = self.merge_xray_h.as_pdb_string(crystal_symmetry=self.get_cryustal_symmetry())

        with open(pdb_file, 'w') as f:
            f.write(pdb_content)