import os
import sys
import pandas
from modules.transformTileCoordinateSystem import transformTileCoordinateSystem

csv_file = sys.argv[1]
prefix = os.path.splitext(csv_file)[0]
tile_grid_shape = (int(sys.argv[2]),int(sys.argv[3]))
tile_size_x = int(sys.argv[4])
tile_size_y = int(sys.argv[5])
csv_df = transformTileCoordinateSystem(csv_file, tile_grid_shape, tile_size_x, tile_size_y)
try:
    prune = sys.argv[6] == "prune"
except:
    prune = False

if prune:
    csv_df = csv_df.loc[:,["Cell_Label", "global_X", "global_Y"]]
    csv_df.rename(columns = {'global_X':'X', 'global_Y':'Y'}, inplace = True)
    csv_df.to_csv("coords.csv", index=False)
else:
    csv_df.to_csv(f"{prefix}_transformed.csv", index=False)




