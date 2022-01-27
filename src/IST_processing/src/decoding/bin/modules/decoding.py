import pandas as pd
import sys
import csv
import re


def decodeSequentialMaxIntensity(path_to_max_intensity_csv: str, path_to_codebook_csv: str, tile_nr: int) -> pd.DataFrame:
    """decodeSequentialMaxIntensity.

    Parameters
    ----------
    path_to_max_intensity_csv : str
        path_to_max_intensity_csv
    path_to_codebook_csv : str
        path_to_codebook_csv
    tile_nr : int
        tile_nr

    Returns
    -------
    pd.DataFrame
        Dataframe that countains each decoded spot.

    """
    # Create empty dict
    codebook_dict = {}
    # First we parse the codebook into a usable dict
    with open(path_to_codebook_csv, 'r') as csvfile:
        # Check if the csv file has a header
        sniffer = csv.Sniffer()
        has_header = sniffer.has_header(csvfile.read(2048))
        csvfile.seek(0)
        if has_header:
            reader = csv.DictReader(csvfile)
            for row in reader:
                codebook_dict[str(row['Barcode']).strip()] = row['Gene'].strip()
        else :
            print("No header was found in the codebook, this might result in unexpected behaviour.")
            for line in csvfile:
                split_line = line.split(',')
                codebook_dict[str(split_line[1])]=split_line[0].strip()
    # Calculate number of rounds based on length of the barcode
    n_rounds = len(list(codebook_dict.values())[0])

    df = pd.read_csv(path_to_max_intensity_csv)
    result_df = pd.DataFrame(columns=['X', 'Y', 'Barcode', 'Gene'])
    # Group df by unique coordinate pair.
    coordinate_grouped_df = df.groupby(['X', 'Y'])
    barcode_list=[]
    #Iterate over each grouped coordinate pair with: name = tuple of X-Ycoordinates, group is a dataframe of only that group's rows.
    for name, group in coordinate_grouped_df:
        barcode = ""
        min_intensity_ratio=1
        # Sort based on rounds, from top to bottom the rows are sequential.
        for row in group.sort_values('Round').itertuples():
            barcode += str(row.Channel)
            if row.Intensity_ratio < min_intensity_ratio:
                min_intensity_ratio = row.Intensity_ratio
        try:
            gene_name = codebook_dict[barcode]
        except  KeyError:
            gene_name=""
        temp_dict = {'Tile': tile_nr, 'X':name[0], 'Y':name[1], 'Barcode':barcode, 'Gene':gene_name, "Intensity_ratio": min_intensity_ratio}
        barcode_list.append(temp_dict)

    result_df = pd.DataFrame(barcode_list)
    return result_df


