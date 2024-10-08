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
from utility import get_residue_indexes



## we need to read in only non altloc atoms to get u aniso becasue that is what IMP will be doing
def read_non_altloc_structure(
    pdb_file
):
    pdb_inp = pdb.input(file_name=str(pdb_file))
    hierarchy = pdb_inp.construct_hierarchy()
    asc = hierarchy.atom_selection_cache()
    sel = asc.selection("not altloc 'B'")
    hierarchy_not_altloc_B = hierarchy.select(sel)

    crystal_symmetry = pdb_inp.crystal_symmetry()

    return hierarchy_not_altloc_B, crystal_symmetry


## maintains invariance between the imp hierarchy list (hs) and a multi state xray hierarchy (multi_xray_h)
class MultiStateMultiConditionModel:
    def __init__(
            self,
            pdb_files,
            w_mat,
            crystal_symmetries
    ):
        self.m = IMP.Model()
        self.set_w_mat(w_mat)

        ## for testing purposes
        if not crystal_symmetries:
            self.crystal_symmetries = [cctbx.crystal.symmetry(
                unit_cell=(100, 100, 100, 90, 90, 90),
                space_group_symbol="P 1"
            )]


        ## pdb files needs to be a list of pdb files or contain a single multi state pdb file
        if not (len(pdb_files) == 1 or len(pdb_files) == self.n_state):
            raise RuntimeError("Invalid number of pdb files ({}) for the number of states ({})".format(len(pdb_files), self.n_state))

        try:
            pdb_file_n_state = utility.get_n_state_from_pdb_file(pdb_files[0])
        except RuntimeError as e:
            raise e

        self.hs = list()

        pdb_selectors = [IMP.atom.NonAlternativePDBSelector()]
        sel = pdb_selectors[0]
        for i in range(1, len(pdb_selectors)):
            sel = IMP.atom.AndPDBSelector(sel, pdb_selectors[i])

        ## if it is a single state pdb file then need to deal with altlocs
        if pdb_file_n_state == 1 and self.n_state >= 1:
            ## create a tmp list to store the individual xray hierarchies
            xray_hs = list()
            self.multi_xray_h = root()

            for state in range(self.n_state):
                pdb_file = pdb_files[state]
                self.hs.append(IMP.atom.read_pdb(str(pdb_file), self.m, sel))
                xray_h, crystal_symmetry = read_non_altloc_structure(pdb_file)
                # self.crystal_symmetry = crystal_symmetry

                ## merge the xray hierarchies as well into the multi state xray hierarchy
                model = xray_h.only_model()
                model.id = str(state+1)
                self.multi_xray_h.append_model(model.detached_copy())
        ## there's an assumption that the multistate pdb file doesn't contain altlocs
        elif pdb_file_n_state == self.n_state:
            pdb_file = pdb_files[0]
            self.hs.extend(IMP.atom.read_multimodel_pdb(str(pdb_file), self.m, sel))
            pdb_inp = pdb.input(file_name=str(pdb_file))
            self.multi_xray_h = pdb_inp.construct_hierarchy()
            # self.crystal_symmetry = pdb_inp.crystal_symmetry()
        else:
            raise RuntimeError("Number of states in pdb file ({}) does not match the number of states in the model ({}).".format(pdb_file_n_state, self.n_state))

        ## the models can only be used to modify atoms but share references with the multi hierarchy
        ## not hierarchies and can't be used to write pdbs or make atom selections
        if len(self.hs) != len(self.multi_xray_h.models()):
            raise RuntimeError("Number of IMP hierarchies ({}) does not match the number of xray hierarchies ({})".format(len(self.hs), len(self.multi_xray_h.models())))

        print("SETTING UP MODEL", pdb_files)
        print("NO STATES: ", len(self.hs))
        print("NO CONDITIONS: ", self.n_cond)
        pids = IMP.atom.Selection(self.hs[0]).get_selected_particle_indexes()
        print("NO OF ATOMS PER STATE: ", len(pids))

        self.water_at_type = IMP.atom.AtomType("HET: O  ")
        # for pid in pids:
        #     if IMP.atom.Atom(self.m, pid).get_atom_type() == self.water_at_type:
        #         IMP.atom.CHARMMAtom.setup_particle(self.m, pid, "O")
        #     # elif IMP.atom.Atom(self.m, pid).get_atom_type() == IMP.atom.AtomType("HET:ZN"):
        #         # IMP.atom.Atom(self.m, pid).set_element(IMP.atom.AtomType("ZN"))
        #         # IMP.atom.CHARMMAtom.setup_particle(self.m, pid, "ZN")

        ## if ligand then create rigid body
        self.ligands = list()
        for h in self.hs:
            ress = IMP.atom.get_by_type(h, IMP.atom.RESIDUE_TYPE)
            for res in ress:
                ## if the res 3 letter code not in the list of standard amino acids and more than 1 atom then it is a ligand
                if res.get_name() not in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HOH"] and len(IMP.atom.Selection(res).get_selected_particle_indexes()) > 1:
                    ## setup rigid body
                    res_pids = IMP.atom.Selection(res).get_selected_particle_indexes()
                    print("LIGAND: ", res.get_name(), len(res_pids))

                    atoms = IMP.core.get_leaves(res)
                    rb = IMP.core.RigidBody.setup_particle(res, atoms)
                    rb.set_coordinates_are_optimized(True)

                    # prb = IMP.core.RigidBody.setup_particle(IMP.Particle(self.m), IMP.algebra.ReferenceFrame3D())

                    # for pid in res_pids:
                    #     d = IMP.core.XYZR(self.m, pid)
                    #     prb.add_member(d)

                    rb.set_coordinates_are_optimized(True)
                    self.ligands.append(rb)
                    # prb.add_to_rotational_derivatives(IMP.algebra.Vector4D(1,0,0,0), IMP.DerivativeAccumulator(1))
                    # print(prb.get_rotational_derivatives())
                # else:
                #     for pid in IMP.atom.Selection(res).get_selected_particle_indexes():
                #         d = IMP.core.XYZR(self.m, pid)
                #         d.set_coordinates_are_optimized(True)


        # Setup coordinates for md
        for h in self.hs:
            for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
                d = IMP.core.XYZR(self.m, pid)
                d.set_coordinates_are_optimized(True)

                IMP.atom.LinearVelocity.setup_particle(self.m, pid)

        ## cctbx setup
        self.perform_checks()
        self.create_lookup_table()

        self.multi_xray_structures = list()
        self.crystal_symmetries = crystal_symmetries
        for i in range(self.n_cond):
            multi_xray_structure = self.multi_xray_h.extract_xray_structure(crystal_symmetry=crystal_symmetries[i])

            multi_xray_structure.scatterers().flags_set_grads(
                state=False
            )
            multi_xray_structure.scatterers().flags_set_grad_site(
                iselection=multi_xray_structure.all_selection().iselection()
            )
            multi_xray_structure.scatterers().flags_set_grad_occupancy(
                iselection=multi_xray_structure.all_selection().iselection()
            )

            self.multi_xray_structures.append(multi_xray_structure)

    def get_m(self):
        return self.m

    def get_hs(self):
        return self.hs

    def get_n_state(self):
        return self.n_state

    def get_w_mat(self):
        return self.w_mat

    def get_pids_in_state(self, state):
        return IMP.atom.Selection(self.hs[state]).get_selected_particle_indexes()

    def get_ca_pids_in_state(self, state):
        return IMP.atom.Selection(self.hs[state], atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()

    def get_water_pids_in_state(self, state):
        return IMP.atom.Selection(self.hs[state], atom_type=self.water_at_type).get_selected_particle_indexes()

    def get_pids(self):
        pids = list()
        for state in range(self.n_state):
            pids.extend(self.get_pids_in_state(state))

        return pids

    def get_ca_pids(self):
        pids = list()
        for state in range(self.n_state):
            pids.extend(self.get_ca_pids_in_state(state))

        return pids

    def get_water_pids(self):
        pids = list()
        for state in range(self.n_state):
            pids.extend(self.get_water_pids_in_state(state))

        return pids

    def get_atoms(self):
        atoms = list()
        for state in range(self.n_state):
            atoms.extend(self.get_atoms_in_state(state))

        return atoms

    def get_atoms_in_state(self, state):
        return self.multi_xray_h.models()[state].atoms()

    def get_crystal_symmetry(self):
        return self.crystal_symmetry

    def get_atom(self, pid):
        return self.pid_to_atom[pid]

    def get_pid(self, atom):
        return self.atom_to_pid[atom]

    def get_com(self):
        ## ConterOfMass does not dynamically update so need to reinstatiant
        self.com = IMP.atom.CenterOfMass.setup_particle(
            IMP.Particle(self.m),
            self.get_ca_pids()
        )
        return self.com

    def get_occs_for_condition_i(self, i):
        return self.w_mat[:,i]

    def get_ligands(self):
        return self.ligands

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
        for state in range(self.n_state):
            h = self.hs[state]
            # h_xray = self.xray_hs[i]

            for pid in self.get_pids_in_state(state):
                at = IMP.atom.Atom(self.m, pid)
                atom_type = at.get_atom_type()
                atom_type = str(atom_type).strip("\"")

                if atom_type == "HET: O  ":
                    atom_type = "O"

                if "HET: " in atom_type:
                    atom_type = atom_type.split("HET: ")[1]
                elif "HET:" in atom_type:
                    atom_type = atom_type.split("HET:")[1]

                # elif atom_type == "HET: S":
                #     atom_type = "S"

                res = IMP.atom.Residue(at.get_parent())
                res_id = res.get_index()

                asc = self.multi_xray_h.atom_selection_cache()
                sel_str = "model {} and resseq {} and name {}".format(state+1, res_id, atom_type)
                sel = asc.selection(sel_str)
                sel_atoms = self.multi_xray_h.select(sel).atoms()

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

    ## update the individual xray structures then get the multi xray structure which updates automatically
    def get_multi_xray_structure(self, cond):
        ## this is the relatively expensive step and theoretically only needs to be done for cond 0
        self.update_xray_hs()

        ## update the coordinates of the multi xray structure from the merged xray hierarchy
        self.multi_xray_structures[cond].set_sites_cart(self.multi_xray_h.atoms().extract_xyz())

        ## update the occs of the multi xray structure from the w_mat directly
        occs = list()
        for state in range(self.n_state):
            occs.extend([self.w_mat[state, cond]]*len(self.get_pids_in_state(state)))
        occs = flex.double(occs)

        self.multi_xray_structures[cond].set_occupancies(occs)

        ## return the multi xray structure
        return self.multi_xray_structures[cond]

    def perform_checks(self):
        ## check for each state the number of pids and the number of atoms are equal
        for i in range(self.n_state):
            n_pids = len(self.get_pids_in_state(i))
            n_atoms = len(self.get_atoms_in_state(i))

            if n_pids != n_atoms:
                raise RuntimeError("Number of atoms in state {} xray structure ({}) does not match the number of atoms in the hierarchy ({})".format(i, n_pids, n_atoms))

    def write_state_pdb_file(self, state, pdb_file):
        ## update the positions before writing
        self.update_xray_hs()

        # need to create a temporary root hierarchy
        tmp_hierarchy = pdb.hierarchy.root()
        tmp_hierarchy.append_model(self.multi_xray_h.models()[state].detached_copy())

        pdb_content = tmp_hierarchy.as_pdb_string(crystal_symmetry=self.crystal_symmetries[0])

        with open(pdb_file, 'w') as f:
            f.write(pdb_content)

    def write_pdb_file(self, pdb_file):
        ## update the positions before writing
        self.update_xray_hs()

        pdb_content = self.multi_xray_h.as_pdb_string(crystal_symmetry=self.crystal_symmetries[0])

        with open(pdb_file, 'w') as f:
            f.write(pdb_content)