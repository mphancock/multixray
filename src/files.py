import pandas as pd


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
            line = line[:56] + occ + line[60:]

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


"""
Reads a PDB file (which may contain multiple models) and returns a pandas DataFrame.

The DataFrame will include the following columns:
    - model: Model number (if available; otherwise, defaults to None)
    - record: Record type (ATOM or HETATM)
    - atom_serial: Atom serial number
    - atom_name: Atom name
    - alt_loc: Alternate location indicator
    - residue_name: Residue name
    - chain_id: Chain identifier
    - residue_seq: Residue sequence number
    - insertion: Insertion code
    - x, y, z: Atom coordinates
    - occupancy: Occupancy
    - temp_factor: Temperature factor
    - element: Element symbol (if available)
    - charge: Charge (if available)

Parameters:
    filename (str): Path to the input PDB file.

Returns:
    pd.DataFrame: DataFrame containing the parsed PDB data.
"""
def pdb_to_df(pdb_file):
    data = []
    current_model = None  # Will hold the model number (if present)

    with open(pdb_file, 'r') as file:
        for line in file:
            # Check for the start of a new model.
            if line.startswith("MODEL"):
                try:
                    # Model number is typically in columns 11-14 (1-indexed)
                    current_model = int(line[10:14].strip())
                except ValueError:
                    current_model = None

            # Process ATOM and HETATM records.
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    record       = line[0:6].strip()         # Columns 1-6
                    atom_serial  = int(line[6:11].strip())     # Columns 7-11
                    atom_name    = line[12:16].strip()         # Columns 13-16
                    alt_loc      = line[16].strip()            # Column 17
                    residue_name = line[17:20].strip()         # Columns 18-20
                    chain_id     = line[21].strip()            # Column 22
                    residue_seq  = int(line[22:26].strip())     # Columns 23-26
                    insertion    = line[26].strip()            # Column 27
                    x            = float(line[30:38].strip())    # Columns 31-38
                    y            = float(line[38:46].strip())    # Columns 39-46
                    z            = float(line[46:54].strip())    # Columns 47-54
                    occupancy    = float(line[54:60].strip())    # Columns 55-60
                    temp_factor  = float(line[60:66].strip())    # Columns 61-66
                except ValueError:
                    # Skip the line if there's a parsing error.
                    continue

                # Element symbol is typically in columns 77-78.
                element = line[76:78].strip() if len(line) >= 78 else ""
                # Charge is typically in columns 79-80.
                charge  = line[78:80].strip() if len(line) >= 80 else ""

                # Append the parsed data as a dictionary, including the current model.
                data.append({
                    "model": current_model,
                    "record": record,
                    "atom_serial": atom_serial,
                    "atom_name": atom_name,
                    "alt_loc": alt_loc,
                    "residue_name": residue_name,
                    "chain_id": chain_id,
                    "residue_seq": residue_seq,
                    "insertion": insertion,
                    "x": x,
                    "y": y,
                    "z": z,
                    "occupancy": occupancy,
                    "temp_factor": temp_factor,
                    "element": element,
                    "charge": charge
                })

            # When an ENDMDL record is encountered, reset the model (optional).
            elif line.startswith("ENDMDL"):
                current_model = None

    return pd.DataFrame(data)


"""
Write a PDB file from a DataFrame that includes a 'model' column.

The DataFrame is expected to contain the following columns:
    'model', 'record', 'atom_serial', 'atom_name', 'alt_loc',
    'residue_name', 'chain_id', 'residue_seq', 'insertion',
    'x', 'y', 'z', 'occupancy', 'temp_factor', 'element', 'charge'

The function groups the atoms by model, and for each model writes:
    MODEL     <model number>
    ... atom records ...
    ENDMDL

Parameters:
    df (pd.DataFrame): DataFrame containing PDB data.
    filename (str): Path to the output PDB file.
    single_model (bool): If True, only write a single model (default: False).
"""
def write_pdb_from_df_with_models(
    df,
    filename,
    single_model=False
):
    with open(filename, "w") as f:
        if single_model:
            f.write("MODEL     {:4d}\n".format(1))

        # Group by model (sorted by the model value)
        for model_number, group in df.groupby("model", sort=True):
            try:
                # Try to convert model_number to an integer.
                model_int = int(model_number)
            except (ValueError, TypeError):
                # If conversion fails (or model is missing), default to 1.
                model_int = 1

            # Write the MODEL header.
            # Format: Columns 1-6: "MODEL " (left-justified), columns 11-14: model number.
            # Here we use a simple formatting string.
            if not single_model:
                f.write("MODEL     {:4d}\n".format(model_int))

            # Write each atom/HETATM record.
            for idx, row in group.iterrows():
                # Format a PDB ATOM/HETATM record using fixed-width fields.
                # The typical PDB columns for an ATOM record are:
                #   1-6   : Record name (ATOM/HETATM, left-justified)
                #   7-11  : Atom serial number (right-justified)
                #   13-16 : Atom name (right-justified)
                #   17    : Alternate location indicator
                #   18-20 : Residue name (right-justified)
                #   22    : Chain identifier
                #   23-26 : Residue sequence number (right-justified)
                #   27    : Insertion code
                #   31-38: X coordinate (8.3f)
                #   39-46: Y coordinate (8.3f)
                #   47-54: Z coordinate (8.3f)
                #   55-60: Occupancy (6.2f)
                #   61-66: Temperature factor (6.2f)
                #   77-78: Element symbol (right-justified)
                #   79-80: Charge (right-justified)
                line = (
                    "{:<6}"    # record name (e.g., ATOM)
                    "{:>5d}"   # atom serial (integer)
                    " "        # spacer (column 12)
                    "{:>4}"    # atom name
                    "{:1}"     # alt_loc
                    "{:>3}"    # residue name
                    " "        # spacer (column 21)
                    "{:1}"     # chain id
                    "{:>4d}"   # residue sequence number (integer)
                    "{:1}"     # insertion code
                    "   "      # three spaces (columns 28-30)
                    "{:>8.3f}" # x coordinate
                    "{:>8.3f}" # y coordinate
                    "{:>8.3f}" # z coordinate
                    "{:>6.2f}" # occupancy
                    "{:>6.2f}" # temp factor
                    "          "  # ten spaces (columns 67-76)
                    "{:>2}"    # element
                    "{:>2}"    # charge
                    "\n"
                ).format(
                    row["record"],
                    int(row["atom_serial"]),
                    row["atom_name"],
                    row["alt_loc"],
                    row["residue_name"],
                    row["chain_id"],
                    int(row["residue_seq"]),
                    row["insertion"],
                    float(row["x"]),
                    float(row["y"]),
                    float(row["z"]),
                    float(row["occupancy"]),
                    float(row["temp_factor"]),
                    row["element"],
                    row["charge"]
                )
                f.write(line)
            # End of model
            if not single_model:
                f.write("ENDMDL\n")

        if single_model:
            f.write("ENDMDL\n")


if __name__ == "__main__":
    from pathlib import Path

    pdb_file = Path(Path.home(), "Documents/xray/tmp/tmp.pdb")
    df = pdb_to_df(pdb_file)
    print(df.head())
    print(df.tail())
    print(len(df))

    write_pdb_from_df_with_models(df, Path(Path.home(), "Documents/xray/tmp/tmp_out.pdb"), single_model=True)
