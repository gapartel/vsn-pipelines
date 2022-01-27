import sys
import re
import os
import pandas as pd
from skimage import io
import numpy as np
from modules.stardistSegmentation import segment, getProperties
from stardist.models import StarDist2D

image_path = sys.argv[1]
img = io.imread(image_path)
prefix = os.path.splitext(image_path)[0]
model_path = sys.argv[2]
model_versatile = StarDist2D(None, name='2D_versatile_fluo', basedir=model_path)
model_versatile.load_weights('weights_best.h5')
try:
    # Check if there's a tile number in the image
    tile_nr = re.findall(r"\d+", re.findall(r"tile\d+", prefix)[0])[0]
except:
    tile_nr=""
labeled_image= segment(img, model_versatile)
attributes_df = getProperties(labeled_image, img, tile_nr)
# labeled_image is not yet viewable, it's just an int64 image with values 0-#objects, it has no pixel meaning.
# For that you need to call skimage.label2rgb
if 'tile_nr' in locals():
    attributes_df['Tile'] = tile_nr
io.imsave(f"{prefix}_labeled.tif", labeled_image, check_contrast=False)
attributes_df.to_csv(f"{prefix}_properties.csv", index=False)

