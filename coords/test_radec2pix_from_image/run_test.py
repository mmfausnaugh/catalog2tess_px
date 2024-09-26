# Convert RA, DEC to pixel coordinates and generate a phot.reg file

import sys
sys.path.insert(0,'..')
import coords
import numpy as np
import pandas as pd

# Read the RA, DEC from the hipparcos tsv file
df = pd.read_csv('hipparcos.tsv', sep='\t')
ra = df['_RAJ2000']
dec = df['_DEJ2000']

# Convert RA, DEC to image pixel coordinates
px, py = coords.radec2pix_from_image(ra,dec,'ref.new')

# Write a phot.reg file with the pixel coordinates
f = open("phot.reg", "w")
for i in range(len(px)):
    reg = 'circle ' + str(px[i]) + ' ' + str(py[i]) + ' 0.5\n'
    reg = reg + 'circle ' + str(px[i]) + ' ' + str(py[i]) + ' 5\n'
    f.write(reg)
f.close()
