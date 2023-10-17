from pathlib import Path


if __name__ == "__main__":
    cxc_file = Path(Path.home(), "Documents/xray/dev/20_images/view_sn.cxc")

    f = open(cxc_file, "w")

    for pdb_dir in ["1_state_ref", "2_state_ref", "4_state_ref"]:
        for pdb_id in range(40):
            f.write("open /Users/matthew/Documents/xray/dev/17_synthetic_native/data/pdbs/{}/{}.pdb\n".format(pdb_dir, pdb_id))
            f.write("open /Users/matthew/Documents/xray/dev/20_images/view.cxc\n")
            f.write("lighting soft\n")
            f.write("graphics silhouettes true\n")
            f.write("hide cartoons\n")
            f.write("show atoms\n")
            f.write("save /Users/matthew/Documents/xray_ms/sn/{}_{}.png\n".format(pdb_dir, pdb_id))
            f.write("close #1\n")
