import sys
import csv
import re
from modules.nearestNeighbourPixelBasedDecoding import decodePixelBased

# Input parsing:

# needs to know up front how large the image is going to be
x_dim = int(sys.argv[1])
y_dim = int(sys.argv[2])

min_area = int(sys.argv[3])
# Extract tile nr
tile_nr =  sys.argv[4]
tile_nr_int = int(re.findall(r"\d+", tile_nr)[0])

codebook = sys.argv[5]
bit_len =  int(sys.argv[6])
threshold = float(sys.argv[7])

# Prefix to be able to sort the images in the correct order
image_prefix= sys.argv[8]
image_path_list = [sys.argv[i] for i in range(9, len(sys.argv))]

# Decode pixelbase
decoded_df = decodePixelBased(x_dim,y_dim, codebook, bit_len, image_path_list, image_prefix,threshold, min_area=min_area)

# Add an extra rown with tile number to the dataframe
decoded_df['Tile'] = [tile_nr_int for i in range(0,len(decoded_df))]
decoded_df.to_csv(f"decoded_{tile_nr}.csv", index=False)
