import pandas as pd


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
    single_model (bool): If True, only write a single model (default: False). If single model, then dont write altconf.
"""
def write_pdb_from_df(
    df,
    out_pdb_file,
    single_model=False
):
    with open(out_pdb_file, "w") as f:
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
                if single_model:
                    alt_loc = row["alt_loc"]
                else:
                    alt_loc = " "

                line = (
                    "{:<6}"    # record name (e.g., ATOM)
                    "{:>5d}"   # atom serial (integer)
                    " "        # spacer (column 12)
                    "{:<4}"    # atom name
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
                    alt_loc,
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


"""
IMP sampling produces a multi model pdb file but phenix.refine requires a single model pdb file with altconfs. Need to fix the occupancies.

Update the 'alt_loc' column in the PDB DataFrame so that each row's alt_loc is
set to a letter corresponding to its model number:
    Model 1 -> 'A'
    Model 2 -> 'B'
    Model 3 -> 'C'
    etc.

Parameters:
    df (pd.DataFrame): DataFrame containing at least the 'model' and 'alt_loc' columns.

Returns:
    pd.DataFrame: The DataFrame with updated 'alt_loc' values.
"""
def update_alt_loc_by_model(df):
    def model_to_letter(model):
        try:
            # Convert the model number to an integer.
            model_int = int(model)
            # Convert model 1 to 'A', model 2 to 'B', etc.
            return chr(64 + model_int)
        except (ValueError, TypeError):
            # If model cannot be interpreted as an integer, return an empty string.
            return ""

    # Update the 'alt_loc' column using the helper function.
    df['alt_loc'] = df['model'].apply(model_to_letter)
    return df


"""
Update the 'occupancy' column in the PDB DataFrame based on an occupancy value provided
for each model.

Parameters:
    df (pd.DataFrame): DataFrame containing at least a 'model' and an 'occupancy' column.
    occupancy_array (list or array): Occupancy values for each model. For example,
        occupancy_array[0] is the occupancy for model 1,
        occupancy_array[1] for model 2, etc.

Returns:
    pd.DataFrame: The DataFrame with updated occupancy values.
"""
def update_occs(
    df,
    occs
):
    def get_occ_for_model(row):
        try:
            # Convert the model value to an integer.
            model_int = int(row['model'])
            # Use the model number (1-indexed) to fetch the occupancy from the occupancy_array.
            if 0 <= model_int - 1 < len(occs):
                return occs[model_int - 1]
            else:
                # If model number is out-of-range, keep the original occupancy.
                return row['occupancy']
        except Exception:
            # In case of any conversion error, keep the original occupancy.
            return row['occupancy']

    # Apply the occupancy update for each row.
    df['occupancy'] = df.apply(get_occ_for_model, axis=1)
    return df


"""
Update (or create) the 'model' column in a DataFrame based on the 'alt_loc' column,
converting letters to numbers such that:
    'A' becomes 1, 'B' becomes 2, etc.

If a row's 'alt_loc' value is empty or missing, the 'model' column is set to None.

Parameters:
    df (pd.DataFrame): DataFrame containing at least the 'alt_loc' column.

Returns:
    pd.DataFrame: The DataFrame with the updated 'model' column.
"""
def update_model_based_on_altconf(df):
    def alt_to_model(alt_char):
        if isinstance(alt_char, str) and alt_char.strip():
            # Convert the character to uppercase and then to a number:
            # ord('A') = 65 so subtracting 64 gives 1 for 'A', 2 for 'B', etc.
            return ord(alt_char.upper()) - 64
        else:
            # Return None if the alternate location indicator is empty or not a string.
            return None

    # Apply the helper function to the 'alt_loc' column.
    df['model'] = df['alt_loc'].apply(alt_to_model)
    return df


"""
Insert a single atom (row) into the PDB DataFrame.

Parameters:
    df (pd.DataFrame): The existing DataFrame containing PDB data.
        Expected columns include:
        'model', 'record', 'atom_serial', 'atom_name', 'alt_loc',
        'residue_name', 'chain_id', 'residue_seq', 'insertion',
        'x', 'y', 'z', 'occupancy', 'temp_factor', 'element', 'charge'
    atom_data (dict): A dictionary containing the new atom's data.
        Example:
            {
                "model": 1,
                "record": "ATOM",
                "atom_serial": 3,
                "atom_name": "C",
                "alt_loc": "A",
                "residue_name": "MET",
                "chain_id": "A",
                "residue_seq": 1,
                "insertion": "",
                "x": 40.000,
                "y": 14.500,
                "z": 3.200,
                "occupancy": 1.00,
                "temp_factor": 50.00,
                "element": "C",
                "charge": ""
            }
    index (int, optional): The index at which to insert the new row.
        If not provided or if index is greater than or equal to the number of rows,
        the atom is appended at the end.

Returns:
    pd.DataFrame: A new DataFrame with the new atom row inserted.
"""
def insert_single_atom(
    df,
    atom_data,
    index=None
):
    # Check that atom_data is a dictionary
    if not isinstance(atom_data, dict):
        raise ValueError("atom_data must be a dictionary representing a single atom.")

    # Create a one-row DataFrame from the dictionary.
    new_row_df = pd.DataFrame([atom_data])

    if index is None or index >= len(df):
        # Append the new row at the end.
        new_df = pd.concat([df, new_row_df], ignore_index=True)
    else:
        # Insert the new row at the specified index.
        df_top = df.iloc[:index].copy()
        df_bottom = df.iloc[index:].copy()
        new_df = pd.concat([df_top, new_row_df, df_bottom], ignore_index=True)

    return new_df


"""
Reads a CSV file containing PDB data and duplicates every heteroatom row
(i.e. where 'record' == 'HETATM') for each unique alternate location (alt_loc)
found in the CSV. The new duplicate rows will have their 'alt_loc' field replaced
with one of the unique alternate location values.

Parameters:
    csv_file (str): Path to the CSV file containing the PDB data.

Returns:
    pd.DataFrame: A new DataFrame with heteroatom rows duplicated for each unique alt_loc.
"""
def duplicate_heteroatoms_for_all_altlocs(df):
    # Extract all unique alt_loc values from the dataframe,
    # ignoring empty strings or NaN values.
    unique_altlocs = [loc for loc in df['alt_loc'].unique()
                      if pd.notnull(loc) and loc != ""]

    # If no valid alternate location values are found, return the original dataframe.
    if not unique_altlocs:
        return df

    # Filter heteroatom rows (where record is "HETATM")
    hetero_df = df[df['record'] == "HETATM"]
    # Keep all non-heteroatom rows unchanged.
    non_hetero_df = df[df['record'] != "HETATM"]

    # List to store the duplicated heteroatom rows.
    duplicated_rows = []
    for _, row in hetero_df.iterrows():
        for alt in unique_altlocs:
            new_row = row.copy()
            new_row['alt_loc'] = alt  # Replace the alt_loc with the current value.
            duplicated_rows.append(new_row)

    # Create a DataFrame from the duplicated rows.
    hetero_dup_df = pd.DataFrame(duplicated_rows)

    # Combine the non-heteroatom rows with the duplicated heteroatom rows.
    result_df = pd.concat([non_hetero_df, hetero_dup_df], ignore_index=True)

    return result_df


"""
Renumber the heteroatom residues in a PDB DataFrame based on the highest residue number
found in standard (ATOM) entries.

The procedure is:
    1. Find the maximum residue number among rows where record == "ATOM".
    2. Compute the starting hetero residue number as:
        new_start = ((max_atom // 100) + 1) * 100 + 1
        For example, if max_atom is 306, then new_start = 401.
    3. Group the heteroatoms (rows with record "HETATM") by a combination of
        ('chain_id', 'residue_seq', 'insertion', 'residue_name') and assign new sequential residue numbers.
    4. Update the heteroatom rows in the DataFrame with the new residue numbers.

Parameters:
    df (pd.DataFrame): A DataFrame containing PDB data with at least the following columns:
        - 'record' (e.g., "ATOM" or "HETATM")
        - 'residue_seq'
        - 'chain_id'
        - 'insertion'
        - 'residue_name'

Returns:
    pd.DataFrame: The DataFrame with heteroatom residues renumbered.
"""
def renumber_hetero_residues(df):
    # Step 1: Find the maximum residue number among standard atoms.
    atom_mask = df['record'] == "ATOM"
    if atom_mask.any():
        max_atom = df.loc[atom_mask, 'residue_seq'].max()
    else:
        raise ValueError("No standard ATOM records found in the DataFrame.")

    # Step 2: Compute the new starting number for heteroatom residues.
    new_start = ((max_atom // 100) + 1) * 100 + 1

    # Step 3: Process heteroatom rows.
    hetero_mask = df['record'] == "HETATM"
    hetero_df = df.loc[hetero_mask].copy()

    # Define the grouping columns that uniquely identify a residue.
    group_cols = ['chain_id', 'residue_seq', 'insertion', 'residue_name']

    # Create a mapping from the original residue identifier to a new residue number.
    mapping = {}
    new_residue_number = new_start

    # Group heteroatom rows by the specified columns.
    for group_key, group_data in hetero_df.groupby(group_cols, sort=False):
        # Assign a new residue number for this group.
        mapping[group_key] = new_residue_number
        new_residue_number += 1

    # Step 4: Update heteroatom rows in the DataFrame.
    # Define a helper function to look up the new residue number for a given row.
    def get_new_residue(row):
        key = (row['chain_id'], row['residue_seq'], row['insertion'], row['residue_name'])
        return mapping.get(key, row['residue_seq'])

    # Apply the helper function to heteroatom rows.
    df.loc[hetero_mask, 'residue_seq'] = hetero_df.apply(get_new_residue, axis=1)

    return df


if __name__ == "__main__":
    from pathlib import Path

    df = pdb_to_df(Path(Path.home(), "Documents/xray/tmp/440_refine_001.pdb"))
    print(len(df))
    df = duplicate_heteroatoms_for_all_altlocs(df)
    df = update_model_based_on_altconf(df)
    df = renumber_hetero_residues(df)
    write_pdb_from_df(df, Path(Path.home(), "Documents/xray/tmp/tmp.pdb"), single_model=False)
