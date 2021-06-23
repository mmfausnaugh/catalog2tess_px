import numpy as np
from astropy.time import Time
import sys
import os
import pandas as pd
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog

class FermiCsv(Catalog): 
    keys = ['source_id', 'ra', 'dec', 'flux', 
            'flux_error', 'detection_significance', 'spect_type', 
            'gammaray_altname', 'assoc_name', 'analysis_flags', 'semimajor_axis',
            'semiminor_axis',  'ability_ind', 'fractional_var',
            'fractional_var_error', 'source_type', 'source_type_alt', 'assoc_name_alt',
            'camcol', 'camrow', 'ccd', 'ccdcol', 'ccdrow']
      

    def __init__(self, ifile, sector, cam,ignore_image_buffer=False):
        d = pd.read_csv(ifile, dtype = {'0': str, '1': float}, skiprows = 4, usecols = range(1, 19), sep='|',
                        engine='python', na_values = "N/A", header=0)
        d = pd.DataFrame(d).to_numpy()
        ra = [[float(e) for e in r.split(' ')] for r  in d[:,1]]
        dec = [[float(e) for e in r.split(' ')] for r in d[:,2]]
        coords = np.array([ self.sexigesimal_to_decimal(z[0],z[1]) for z in zip(ra,dec)])
        ra  = coords[:,0]
        dec = coords[:,1]
        
        out = self.get_pix('s'+str(sector), int(cam - 1.0), ra,dec, d[:,3], ignore_image_buffer=ignore_image_buffer)
        
        # any errors to correct? 
        camcol = out[0]
        camrow = out[1]
        ccd    = out[2]
        ccdcol = out[3]
        ccdrow = out[4]
        idx = out[8]
        
        d = d[idx]
        ra = ra[idx]
        dec = dec[idx]
        
        super(FermiCsv,self).__init__(self.keys,                                     
                                         [d[:,0], d[:,1], d[:,2], d[:,3], d[:,4],
                                          d[:,5], d[:,6], d[:,7], d[:,8], d[:,9], d[:,10],
                                          d[:,11], d[:,12], d[:,13], d[:,14], d[:,15],
                                          d[:,16], d[:,17], 
                                          camcol, camrow, ccd, ccdcol, ccdrow])

        #self.source_id = self.source_id.astype(int)

        self.obj_name = self.source_id
