import sys
from pathlib import Path

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
import utility


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
        pdb_selectors = [IMP.atom.ChainPDBSelector(["A"]), IMP.atom.NonWaterPDBSelector(), IMP.atom.NonAlternativePDBSelector(), IMP.atom.ATOMPDBSelector()]
        sel = pdb_selectors[0]
        for i in range(1, len(pdb_selectors)):
            sel = IMP.atom.AndPDBSelector(sel, pdb_selectors[i])

        if pdb_file_n_state == 1 and self.n_state > 1:
            for i in range(self.n_state):
                self.hs.append(IMP.atom.read_pdb(str(pdb_file), self.m, sel))
        elif pdb_file_n_state == self.n_state:
            self.hs.extend(IMP.atom.read_multimodel_pdb(str(pdb_file), self.m, sel))
        else:
            raise RuntimeError("Number of states in pdb file ({}) does not match the number of states in the model ({}).".format(pdb_file_n_state, self.n_state))

        # Set b factors
        avg_b_factor = 0
        for h in self.hs:
            for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
                avg_b_factor += IMP.atom.Atom(self.m, pid).get_temperature_factor()

        avg_b_factor /= len(self.hs) * len(IMP.atom.Selection(self.hs[0]).get_selected_particle_indexes())

        print("AVG B FACTOR: ", avg_b_factor)

        for h in self.hs:
            for pid in IMP.atom.Selection(h).get_selected_particle_indexes():
                # IMP.atom.Atom(self.m, pid).set_temperature_factor(15)
                IMP.atom.Atom(self.m, pid).set_temperature_factor(avg_b_factor)

        self.pids_dict = dict()
        self.prot_pids_dict = dict()
        self.water_pids_dict = dict()
        self.main_pids_dict = dict()
        self.side_pids_dict = dict()
        self.ca_pids_dict = dict()
        water_at_type = IMP.atom.AtomType("HET: O  ")

        # print(self.n_state, len(self.hs))
        for i in range(self.n_state):
            h = self.hs[i]
            self.pids_dict[i] = IMP.atom.Selection(h).get_selected_particle_indexes()
            self.prot_pids_dict[i] = (IMP.atom.Selection(h) - IMP.atom.Selection(h, atom_type=water_at_type)).get_selected_particle_indexes()
            self.water_pids_dict[i] = IMP.atom.Selection(h, atom_type=water_at_type).get_selected_particle_indexes()
            self.main_pids_dict[i] = IMP.atom.Selection(h, atom_types=[IMP.atom.AT_CA, IMP.atom.AT_C, IMP.atom.AT_O, IMP.atom.AT_N]).get_selected_particle_indexes()
            self.side_pids_dict[i] = list(set(self.prot_pids_dict[i]) - set(self.main_pids_dict[i]))
            self.ca_pids_dict[i] = IMP.atom.Selection(h, atom_type=IMP.atom.AtomType("CA")).get_selected_particle_indexes()

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

    def get_m(self):
        return self.m

    def get_hs(self):
        return self.hs

    def get_w_mat(self):
        return self.w_mat

    def get_all_pids(self):
        pids = list()
        for i in range(self.n_state):
            pids.extend(self.pids_dict[i])

        return pids

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
