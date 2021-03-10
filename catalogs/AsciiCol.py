import scipy as sp
from astropy.time import Time
import sys
import os
import pandas as pd
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog
        
class AsciiCol(Catalog):
    keys = ['obj_name','ra','dec','mag',
            'camcol','camrow','ccd','ccdcol','ccdrow']

    def __init__(self, ifile, sector , cam, sexagesimal=False, ignore_image_buffer=False ):
        d = sp.genfromtxt(ifile,dtype=str)
        if d.ndim == 1:
            d = sp.array([d])

        if sexagesimal:
            ra = [[float(e) for e in r.split(':')] for r  in d[:,1]]
            dec = [[float(e) for e in r.split(':')] for r in d[:,2]]
            coords = sp.array([ self.sexigesimal_to_decimal(z[0],z[1]) for z in zip(ra,dec)])
            ra  = coords[:,0]
            dec = coords[:,1]

        else:
            ra  = d[:,1]
            dec = d[:,2]


            
        out = self.get_pix('s'+str(int(sector)), int(cam - 1.0) ,ra,dec, d[:,3],
                           ignore_image_buffer=ignore_image_buffer)
        camcol = out[0]
        camrow = out[1]
        ccd    = out[2]
        ccdcol = out[3]
        ccdrow = out[4]
        idx = out[8]
        
#        print sector, cam, idx
#        print d
        d = d[idx]
#        print d, ccd
        ra = ra[idx]
        dec = dec[idx]

        super(AsciiCol,self).__init__(self.keys,
                                      [d[:,0], ra, dec, d[:,3],
                                       camcol, camrow, ccd, ccdcol, ccdrow])

