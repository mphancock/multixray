from pathlib import Path
import sys
import pandas as pd
import numpy as np
import multiprocessing
import argparse
import time
import shutil
import os

import IMP
import IMP.atom

sys.path.append(str(Path(Path.home(), "xray/src")))
from refine import read_pdb_and_refine, read_pdb_and_refine_to_max_ff, refine_posterior
from params import read_job_csv, build_weights_matrix
from align_imp import align_one_to_two
sys.path.append(str(Path(Path.home(), "xray/score_bench/src")))
from score import pool_score
from files import pdb_to_df, write_pdb_from_df, update_occs, update_model_based_on_altconf, update_alt_loc_by_model, insert_single_atom, renumber_hetero_residues, duplicate_heteroatoms_for_all_altlocs
from miller_ops import get_crystal_symmetry
from etc import get_zn_coords_and_occ_after_align


if __name__ == "__main__":
    pdb_file = Path(Path.home(), "xray/dev/26_phenix_refine/data/ref_out.pdb")
    tmp_pdb_file = Path(Path.home(), "xray/tmp/ref_out.pdb")
    decoy_w_mat = np.array([[0.1835615230208595, 0.3552770352176667],[0.2934976053478322,0.0332077987252963],[0.1920311647301441,0.1573457315428919],[0.3309097069011641,0.454169434514145]])

    ## refine against all cif files and save as PDB_CIF.pdb
    cif_files = [Path("/wynton/home/sali/mhancock/xray/dev/38_standard_flags/data/7mhj.cif"), Path("/wynton/home/sali/mhancock/xray/dev/38_standard_flags/data/7mhk.cif")]
    cif_names = [cif_file.stem for cif_file in cif_files]

    cond = 0

    cif_file = cif_files[cond]
    cif_name = cif_names[cond]

    ## get the pdb csv
    pdb_df = pdb_to_df(pdb_file)

    ## update confs and occs
    pdb_df = update_occs(pdb_df, decoy_w_mat[:,cond])
    pdb_df = update_alt_loc_by_model(pdb_df)

    ## get zn coordinates and occs and add to pdb csv
    zn_coords, zn_occ = get_zn_coords_and_occ_after_align(
        pdb_file_1=Path(Path.home(), "xray/data/pdbs/7mhf/{}.pdb".format(cif_name)),
        pdb_file_2=pdb_file
    )
    pdb_df = insert_single_atom(
        df=pdb_df,
        atom_data={"model": 1, "record": "HETATM","atom_serial": 3000, "atom_name": "ZN", "alt_loc": "", "residue_name": "ZN", "chain_id": "A","residue_seq": 401, "insertion": "", "x": zn_coords[0], "y": zn_coords[1], "z": zn_coords[2], "occupancy": zn_occ, "temp_factor": 53.43, "element": "ZN", "charge": ""},
        index=None
    )

    ## save as altconf
    write_pdb_from_df(
        df=pdb_df,
        out_pdb_file=tmp_pdb_file,
        single_model=True
    )

    ## convert the orig multistate pdb file to altconf tmp file (jth column of the decoy_w_mat)
    cif_file = cif_files[cond]
    crystal_symmetry = get_crystal_symmetry(cif_file)
    refined_pdb_file = Path(Path.home(), "xray/tmp/{}_{}.pdb".format(pdb_file.stem, cif_names[cond]))
    print(pdb_file, tmp_pdb_file, refined_pdb_file)
    sg = crystal_symmetry.space_group_info().group().info()
    uc = crystal_symmetry.unit_cell()

    refine_command = "phenix.refine {} {} strategy=individual_sites+individual_adp ordered_solvent=true ordered_solvent.mode=every_macro_cycle  refinement.input.xray_data.labels=_refln.F_meas_au,_refln.F_meas_sigma_au refinement.pdb_interpretation.clash_guard.nonbonded_distance_threshold=None crystal_symmetry.unit_cell={},{},{},{},{},{} crystal_symmetry.space_group='{}' write_eff_file=false write_geo_file=false write_def_file=false write_maps=false write_map_coefficients=false write_model_cif_file=false".format(cif_file, tmp_pdb_file, uc.parameters()[0], uc.parameters()[1], uc.parameters()[2], uc.parameters()[3], uc.parameters()[4], uc.parameters()[5], sg)

    print(refine_command)
    os.system(refine_command)
    phenix_out_pdb_file = Path(Path.home(), "xray/tmp/{}_refine_001.pdb".format(pdb_file.stem))

    # ## convert all models back to multistate
    df = pdb_to_df(phenix_out_pdb_file)
    df = duplicate_heteroatoms_for_all_altlocs(df)
    df = update_model_based_on_altconf(df)
    df = renumber_hetero_residues(df)
    write_pdb_from_df(
        df=df,
        out_pdb_file=refined_pdb_file,
        single_model=False
    )

    # ## clean up temporary files from phenix
    # os.system("rm *.pdb")
    # os.system("rm *.mtz")
    # os.system("rm *.log")

