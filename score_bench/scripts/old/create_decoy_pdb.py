from pathlib import Path


def extract_file(
        file
):
    decoy_f = open(file, "r")
    lines = decoy_f.read().splitlines()
    lines = [line.rstrip('\n') for line in lines]

    items = list()
    for line in lines:
        print(line)
        item_line = line.split(",")[0:-1]
        items.append([item for item in item_line])

    return items


if __name__ == "__main__":
    decoy_recipe = "4n7f_N_1000_x1"

    bench_dir = Path(Path.home(), "xray/score_rmsd_benchmark")
    data_dir = Path(bench_dir, "data")
    decoy_dir = Path(data_dir, "decoys/decoy_recipes")

    decoy_recipe_file = Path(data_dir, "decoy_recipes/{}.txt".format(decoy_recipe))

    decoy_list = extract_file(
        decoy_recipe_file
    )

    print(decoy_list)