from pathlib import Path


if __name__ == "__main__":
    cxc_file = Path(Path.home(), "Documents/xray/dev/20_images/view_sn.cxc")
    f = open(cxc_file, "w")

    pdb_dir = Path(Path.home(), "Documents/xray/dev/29_synthetic_native_3/data/pdbs/7mhf_30")
    for pdb_id in range(10):
        f.write("open {}/{}.pdb\n".format(pdb_dir, pdb_id))
        f.write("open /Users/matthew/Documents/xray/dev/20_images/view_mpro.cxc\n")
        f.write("save /Users/matthew/Documents/xray_ms/figure_3/7mhf_30/{}.png\n".format(pdb_id))
        f.write("close #1\n")
