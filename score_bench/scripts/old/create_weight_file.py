from pathlib import Path
import random


def create_weight_file(
        weights_file,
        N_decoys,
        N_struct,
        w_type
):
    decoy_f = open(weights_file, "w")
    for i in range(N_decoys):
        weights = list()
        for j in range(N_struct):
            if w_type == "uni":
                w = 1 / N_struct
            elif w_type == "rand":
                w = random.random()
            else:
                raise RuntimeError("Invalid w_type")

            weights.append(w)

        # Normalize to sum to 1.
        w_sum = 0
        for w in weights:
            w_sum = w_sum + w

        for i in range(len(weights)):
            weights[i] = weights[i] / w_sum

        decoy_line = ""
        for w in weights:
            decoy_line = decoy_line + "{},".format(w)
        decoy_line = decoy_line + "\n"
        decoy_f.write(decoy_line)

    decoy_f.close()


if __name__ == "__main__":
    weights_dir = Path(Path.home(), "xray/benchmark_score_rmsd/data/weights")
    N_decoys = 1000
    N_struct = 5
    w_type = "rand"

    weights_file = Path(weights_dir, "N_{}_N_struct_{}_{}.txt".format(N_decoys, N_struct, w_type))

    create_weight_file(
        weights_file=weights_file,
        N_decoys=N_decoys,
        N_struct=N_struct,
        w_type=w_type
    )


