import numpy as np
import os
import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from catalogs.VSX import VSX
from catalogs.GaiaCsv import GaiaCsv
from coords.match_coords import *

import os

#match primary mission to start

sectors = np.r_[1:27]
cams = [1,2,3,4]


dstem = os.path.abspath(__file__)
dstem = os.path.dirname(
    os.path.dirname(dstem))
print(dstem)

for sector in sectors:
    for cam in cams:
        gaia_cat = GaiaCsv( 
            os.path.join(dstem, 'Gaia/s{:02d}/gaia_dr2_s{:02d}_cam{:1d}.txt'.format(sector, sector, cam) )
        )
        vsx_cat = VSX(
            os.path.join(dstem,'VSX/s{:02d}/vsx_s{:02d}_cam{:1d}.txt'.format(sector, sector, cam)  )
        )
        
        for ccd in range(4):
            ccd_m1 = (gaia_cat.ccd == ccd + 1) & (gaia_cat.phot_rp_mean_mag <= 17)
            ccd_m2 = vsx_cat.ccd == ccd + 1
            
            
            gaia_idx, vsx_idx, dr = chunk_match(gaia_cat.ccdcol[ccd_m1],
                                                gaia_cat.ccdrow[ccd_m1],
                                                vsx_cat.ccdcol[ccd_m2],
                                                vsx_cat.ccdrow[ccd_m2],
                                                threshold = 0.1)


            outdir = 's{:02d}'.format(sector)
            if not os.path.isdir(outdir):
                os.makedirs(outdir)

            with open( os.path.join( outdir ,  
                                     'match_gaia_VSX_cam{}_ccd{}.txt'.format(cam, ccd + 1)
                                 ),'w') as fout:
                fout.write('{:>25s} {:>30s} {:>25s} '
                           '{:>5s} {:>5s} {:>5s} '
                           '{:>8s} {:>8s} {:>8s} {:>8s} {:>8s}\n'.format('gaia_id',
                                                                         'VSX_id',
                                                                         'VSX_class',
                                                                         'RP',
                                                                         'mag',
                                                                         'filt',
                                                                         'gaia_col','gaia_row',
                                                                         'VSX_col','VSX_row',
                                                                         'delta_px'))
                for ii in range(len(gaia_idx)):            
                    fout.write('{:>25d} {:>30s} {:>25s} '
                               '{:>5.2f} {:>5.2f} {:>5s} '
                               '{:>8.2f} {:>8.2f} {:>8.2f} {:>8.2f} {:>8.4f}\n'.format(                   
                                   gaia_cat.source_id[ccd_m1][gaia_idx[ii]],
                                   vsx_cat.name[ccd_m2][vsx_idx[ii]],
                                   vsx_cat.var_type[ccd_m2][vsx_idx[ii]],
                           
                                   gaia_cat.phot_rp_mean_mag[ccd_m1][gaia_idx[ii]],
                                   vsx_cat.mean_mag[ccd_m2][vsx_idx[ii]],
                                   vsx_cat.bright_passband[ccd_m2][vsx_idx[ii]],
                                   
                                   gaia_cat.ccdcol[ccd_m1][gaia_idx[ii]],
                                   gaia_cat.ccdrow[ccd_m1][gaia_idx[ii]],
                                   
                                   vsx_cat.ccdcol[ccd_m2][vsx_idx[ii]],
                                   vsx_cat.ccdrow[ccd_m2][vsx_idx[ii]],
                                   dr[ii],
                               ))
        
