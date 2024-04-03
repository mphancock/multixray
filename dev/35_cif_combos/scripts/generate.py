import pandas as pd
from pathlib import Path
from itertools import chain, combinations

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))


if __name__ == "__main__":
    cif_files = list()
    ref_files = list()
    for cif_name in ["7mhf", "7mhg", "7mhh", "7mhi", "7mhj", "7mhk"]:
        cif_files.append(str(Path(Path.home(), "xray/data/cifs/7mhf/{}.cif".format(cif_name))))
        ref_files.append(str(Path(Path.home(), "xray/data/pdbs/7mhf/{}.pdb".format(cif_name))))

    cif_combos = list(powerset(cif_files))
    ref_combos = list(powerset(ref_files))
    cif_combo_df = pd.DataFrame(dtype="string", columns=["cifs", "refs", "ref_occs"])
    for i in range(1, len(cif_combos)):
        cif_combo = cif_combos[i]
        ref_combo = ref_combos[i]
        cif_str = cif_combo[0]
        ref_str = ref_combo[0]
        ref_occ_str = "1"
        for j in range(1, len(cif_combo)):
            cif_str = cif_str + "," + cif_combo[j]
            ref_str = ref_str + "," + ref_combo[j]
            ref_occ_str = ref_occ_str + ";1"

        cif_combo_df.loc[i-1, "cifs"] = cif_str
        cif_combo_df.loc[i-1, "refs"] = ref_str
        cif_combo_df.loc[i-1, "ref_occs"] = ref_occ_str

    cif_combo_df.to_csv(Path(Path.home(), "xray/dev/35_cif_combos/data/7mhf.csv"))

