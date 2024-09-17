import math
from pathlib import Path
import pickle
import random

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import pandas as pd
import numpy as np

import IMP
import IMP.atom

import sys
sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import refine_hs_to_max_ff


if __name__ == "__main__":
    pdb_1 = Path(Path.home(), "xray/data/pdbs/7mhf/7mhf_refine.pdb")
    pdb_2 = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/7mhf/217_no_ref_com/5/output_12/pdbs/500.pdb")

    m_1, m_2 = IMP.Model(), IMP.Model()
    h_1 = IMP.atom.read_pdb(str(pdb_1), m_1, IMP.atom.AllPDBSelector())
    h_2 = IMP.atom.read_pdb(str(pdb_2), m_2, IMP.atom.AllPDBSelector())

    pids_1 = IMP.atom.Selection(h_1, atom_type=IMP.atom.AT_CA).get_selected_particle_indexes()
    pids_2 = IMP.atom.Selection(h_2, atom_type=IMP.atom.AT_CA).get_selected_particle_indexes()

    com_1 = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(m_1), pids_1)
    com_2 = IMP.atom.CenterOfMass.setup_particle(IMP.Particle(m_2), pids_2)

    print(com_1.get_coordinates(), com_2.get_coordinates())
    print(np.linalg.norm(com_1.get_coordinates()-com_2.get_coordinates()))