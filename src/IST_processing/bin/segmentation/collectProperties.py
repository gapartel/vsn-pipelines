import sys
import pandas as pd
from skimage import io
import numpy as np
from modules.thresholding import collectProperties

properties = [sys.argv[i] for i in range(1,len(sys.argv))]
concat_df = collectProperties(properties)
concat_df.to_csv(f"concat_segmented_properties.csv", index=False) 

outs_df =concat_df.loc[:,['Cell_Label', 'X', 'Y', 'Tile']]
outs_df.to_csv("local_coords.csv", index = False)

