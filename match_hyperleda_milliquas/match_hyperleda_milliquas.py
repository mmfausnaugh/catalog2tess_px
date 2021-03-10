import numpy as np
import os
import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from catalogs.MilliquasCsv import MilliquasCsv
from catalogs.HyperLedaCsv import HyperLedaCsv
from coords.match_coords import *


#match primary mission to start

sectors = np.r_[np.r_[1:8],19,20,21]
cams = [1,2,3,4]


dstem = os.path.abspath(__file__)
dstem = os.path.dirname(
    os.path.dirname(dstem))
print(dstem)

for sector in sectors:
    for cam in cams:
        hyperleda_cat = HyperLedaCsv( 
            os.path.join(dstem, 'HyperLEDA/s{:02d}/hyperleda_s{:02d}_cam{:1d}.txt'.format(sector, sector, cam) )
        )
        milliquas_cat = MilliquasCsv(
            os.path.join(dstem,'Milliquas/s{:02d}/milliquas_s{:02d}_cam{:1d}.txt'.format(sector, sector, cam)  )
        )
        
        for ccd in range(4):
            ccd_m1 = (hyperleda_cat.ccd == ccd + 1) & (hyperleda_cat.imag <= 20)
            ccd_m2 = milliquas_cat.ccd == ccd + 1
            
            
            hyperleda_idx, milliquas_idx, dr = chunk_match(hyperleda_cat.ccdcol[ccd_m1],
                                                           hyperleda_cat.ccdrow[ccd_m1],
                                                           milliquas_cat.ccdcol[ccd_m2],
                                                           milliquas_cat.ccdrow[ccd_m2],
                                                           threshold = 0.1) #matching threshold in pixels


            outdir = 's{:02d}'.format(sector)
            if not os.path.isdir(outdir):
                os.makedirs(outdir)

            with open( os.path.join( outdir ,  
                                     'match_hyperleda_milliquas_cam{}_ccd{}.txt'.format(cam, ccd + 1)
                                 ),'w') as fout:
                fout.write('{:>35s} {:>35s} '
                           '{:>12s} {:>12s} '
                           '{:>10s} {:>10s} {:>10s} '
                           '{:>7s} {:>7s} '
                           '{:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n'.format('hyperleda_id',
                                                                         'milliquas_id',
                                                                         'hyp_class',
                                                                         'mil_class',
                                                                         'hyp_BV',
                                                                         'hyp_z',
                                                                         'mil_z',
                                                                         'hyp_mag',
                                                                         'mil_mag',
                                                                         'hyp_col','hyp_row',
                                                                         'mil_col','mil_row',
                                                                         'delta_px'))
                for ii in range(len(hyperleda_idx)):            
                    fout.write('{:>35s} {:>35s} '
                               '{:>12s} {:>12s} '
                               '{:>10.3f} {:>10.6f} {:>10.6f} '
                               '{:>7.3f} {:>7.3f} '
                               '{:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.4f}\n'.format(                   
                                   hyperleda_cat.objname[ccd_m1][hyperleda_idx[ii]],
                                   milliquas_cat.name[ccd_m2][milliquas_idx[ii]],

                                   hyperleda_cat.agnclass[ccd_m1][hyperleda_idx[ii]],
                                   milliquas_cat.qso_class[ccd_m2][milliquas_idx[ii]],

                                   hyperleda_cat.corrected_BV[ccd_m1][hyperleda_idx[ii]],
                                   hyperleda_cat.velocity[ccd_m1][hyperleda_idx[ii]]/2.99792e10,
                                   milliquas_cat.redshift[ccd_m2][milliquas_idx[ii]],
                           
                                   hyperleda_cat.imag[ccd_m1][hyperleda_idx[ii]],
                                   milliquas_cat.rmag[ccd_m2][milliquas_idx[ii]],
                                   
                                   hyperleda_cat.ccdcol[ccd_m1][hyperleda_idx[ii]],
                                   hyperleda_cat.ccdrow[ccd_m1][hyperleda_idx[ii]],
                                   
                                   milliquas_cat.ccdcol[ccd_m2][milliquas_idx[ii]],
                                   milliquas_cat.ccdrow[ccd_m2][milliquas_idx[ii]],
                                   dr[ii],
                               ))
        
