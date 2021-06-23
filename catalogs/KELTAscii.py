import scipy as sp
from astropy.time import Time
import sys
import os
import pandas as pd
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog



class KELTAscii(Catalog):
    keys = ['KELT_ID', 'MASS_ID', 'TIC_ID', 'In_CTL', 'TESS_priority',
            'FP_type_name', 'FP_type', 'RA_hours', 'ra', 'RA_hms',
            'Dec_dms', 'dec', 'Galactic_long', 'Galactic_lat',
            'Ecliptic_long', 'Ecliptic_lat', 'Vmag', 'Tc', 'Tc_err',
            'Period_days', 'Period_err', 'Duration_hrs', 'KELT_depth_mmag',
            'EB_K_kms', 'EB_depth_mmag', 'NEB_RA', 'NEB_Dec', 'NEB_dist_text',
            'NEB_dist_arcsec', 'NEB_dist_is_approx_flag', 'NEB_direction',
            'NEB_depth_text', 'NEB_depth_percent', 'NEB_depth_is_approx_flag',
            'NEB_depth_is_lower_limit', 'NEB_obs_epoch', 'NEB_obs_filter',
            
            'camcol','camrow','ccd','ccdcol','ccdrow']
      

    #special call for KELT because the catalog files are not broken up by sector/cam
    def __init__(self, ifile, sector, cam,ignore_image_buffer=False):
        d = pd.read_csv(ifile, dtype = {'0': str, '1': float}, skiprows = 1, usecols = range(37), sep='\t',
                        engine='python', na_values = "N/A", header=None)
        print(d)
        d = pd.DataFrame(d).to_numpy()
        
        ra = d[:,8]
        dec = d[:,11]
        
        out = self.get_pix('s'+str(sector), int(cam - 1.0) ,ra,dec, d[:,3],
                           ignore_image_buffer=ignore_image_buffer)
        camcol = out[0]
        camrow = out[1]
        ccd    = out[2]
        ccdcol = out[3]
        ccdrow = out[4]
        idx = out[8]
        
        # print sector, cam, idx
        # print d
        d = d[idx]
        # print d, ccd
        ra = ra[idx]
        dec = dec[idx]
        
        
        super(KELTAscii,self).__init__(self.keys,
                                         [d[:,0], d[:,1], d[:,2], d[:,3], d[:,4],
                                          d[:,5], d[:,6], d[:,7], ra, d[:,9], d[:,10],
                                          dec, d[:,12], d[:,13], d[:,14], d[:,15],
                                          d[:,16], d[:,17], d[:,18], d[:,19], d[:,20],
                                          d[:,21], d[:,22], d[:,23], d[:,24], d[:,25],
                                          d[:,26], d[:,27], d[:,28], d[:,29], d[:,30],
                                          d[:,31], d[:,32], d[:,33], d[:,34], d[:,35],
                                          d[:,36],
                                          
                                          camcol, camrow, ccd, ccdcol, ccdrow])
        
        #all classes must have obj_name attribute
        self.obj_name = self.KELT_ID
