

"""
    IMP sampling produces a multi model pdb file but phenix.refine requires a single model pdb file with altconfs. Need to fix the occupancies.
"""
def multi_to_altconf(
    in_pdb_file,
    occs,
    out_pdb_file
):
    with open(in_pdb_file, 'r') as file:
        lines = file.readlines()

    cnt = 0
    modified_lines = []
    in_section = False  # Flag to indicate if we are between 'MODEL' and 'TER'
    for line in lines:
        # Check if the line starts with 'MODEL'
        if line.startswith('MODEL'):
            in_section = True
            modified_lines.append(line)  # Keep the 'MODEL' line as is
            continue  # Proceed to the next line

        # Check if the line starts with 'TER'
        elif line.startswith('ENDMDL'):
            in_section = False
            cnt += 1
            modified_lines.append(line)  # Keep the 'TER' line as is
            continue  # Proceed to the next line

        # Modify lines only if we are within the 'MODEL' and 'TER' section
        ## Also modify the occupancies (56-59)
        if in_section and len(line) >= 17:
            # Replace the 17th character (index 16) with the replacement character
            line = line[:16] + chr(ord("A")+cnt) + line[17:]
            occ = f"{occs[cnt]:.2f}"
            line = line[:55] + occ + line[60:]

        # Append the modified or unmodified line to the list
        modified_lines.append(line)

    # Find all indices of lines containing "MODEL" and "ENDMDL"
    model_indices = [i for i, line in enumerate(modified_lines) if 'MODEL' in line]
    endmdl_indices = [i for i, line in enumerate(modified_lines) if 'ENDMDL' in line]

    # Determine which indices to keep
    keep_indices = set()
    if model_indices:
        keep_indices.add(model_indices[0])  # First "MODEL"
    if endmdl_indices:
        if len(endmdl_indices) > 1:
            keep_indices.add(endmdl_indices[-1])  # Last "ENDMDL"

    # Write the processed lines to the output file
    with open(out_pdb_file, 'w') as f:
        for i, line in enumerate(modified_lines):
            if ('MODEL' in line or 'ENDMDL' in line):
                if i in keep_indices:
                    f.write(line)
                else:
                    continue  # Skip the line
            else:
                f.write(line)


def altconf_to_multi(
    in_pdb_file,
    out_pdb_file,
    n_state
):
    with open(in_pdb_file, 'r') as file:
        lines = file.readlines()

    # Create the conformation dictionary with dynamic keys
    conformation_dict = {chr(65 + i): [] for i in range(n_state)}

    # Iterate over each line and sort based on the conformation identifier (character at column 17)
    for line in lines:
        if line.startswith('ATOM'):
            conformation_id = line[16]  # Character at column 17 (0-indexed position 16)

            mod_line = line[:16] + ' ' + line[17:]  # Replace the conformation identifier with a space

            mod_line = mod_line.strip("\n")
            conformation_dict[conformation_id].append(mod_line)
        ## if water then add to all conformations
        elif line.startswith("HETATM"):
            for conformation_id in conformation_dict.keys():
                conformation_dict[conformation_id].append(line)

    ## Fix the atom counts (8 - 11)
    for conf_id in conformation_dict.keys():
        atom_count = 1
        for i in range(len(conformation_dict[conf_id])):
            line = conformation_dict[conf_id][i]

            line = line[:7] + str(atom_count).rjust(4) + line[11:]
            conformation_dict[conf_id][i] = line
            atom_count += 1

    ## Fix the water residue numbers (23-25)
    for conf_id in conformation_dict.keys():
        res_id = None
        for i in range(len(conformation_dict[conf_id])):
            line = conformation_dict[conf_id][i]
            if line.startswith("ATOM"):
                res_id = int(line[23:26].strip())

            if line.startswith("HETATM"):
                res_id += 1 # Increment the residue number
                line = line[:23] + str(res_id).rjust(3) + line[26:]
                conformation_dict[conf_id][i] = line


    # Combine the sorted lines based on conformation order (A, B, C, D)
    sorted_entries = []
    for i in range(n_state):
        key = chr(ord('A') + i)
        sorted_entries.append(f"MODEL    {i+1}")
        sorted_entries.extend(conformation_dict[key])
        sorted_entries.append("ENDMDL")

    # Join the sorted lines and return as a single string
    new_str = "\n".join(sorted_entries)

    with open(out_pdb_file, 'w') as file:
        file.writelines(new_str)


if __name__ == "__main__":
    from pathlib import Path

    pdb_file = Path("/wynton/group/sali/mhancock/xray/sample_bench/out/280_exp_all_2/9/output_0/pdbs/499.pdb")
    multi_to_altconf(pdb_file, [0.35, 0.65], Path(Path.home(), "xray/tmp/out.pdb"))

    # altconf_to_multi(Path(Path.home(), "xray/tmp/499_refine_001.pdb"), Path(Path.home(), "xray/tmp/out_multi.pdb"), 2)
