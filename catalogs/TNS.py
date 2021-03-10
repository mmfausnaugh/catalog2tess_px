import scipy as sp
from astropy.time import Time
import sys
import os
import pandas as pd
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog


class TNS(Catalog):

    keys = ['prefix','name','group','internal_name',
            'ra_sexagesimal','dec_sexagesimal','ra','dec',
            't_disc','tjd', 'obj_type', 'mag', 'fband', 'z', 'host','host_z',
            'camcol','camrow','ccd','ccdcol','ccdrow']
    
    def __init__(self,ifile,ignore_image_buffer=False):
        try:
            d = sp.genfromtxt(ifile,dtype=str)
        except:
            print('error in {}'.format(ifile))
            raise

        ra = [[float(e) for e in r.split(':')] for r  in d[:,4]]
        dec = [[float(e) for e in r.split(':')] for r in d[:,5]]
        coords = sp.array([ self.sexigesimal_to_decimal(z[0],z[1]) for z in zip(ra,dec)])
        ra  = coords[:,0]
        dec = coords[:,1]

        t1 = d[:,6]
        t2 = d[:,7]
        t = sp.array([ z[0]+'T'+z[1] for z in zip(t1,t2)])
        try:
            t = Time(t, format='isot',scale='utc')
        except:
            print(ifile)
            for tuse in t:
                print(tuse)
                Time(tuse,format='isot',scale='utc')
        tjd = t.jd - 2457000.0


        fstem = os.path.basename(ifile)
        sector = 's' + fstem.split('_')[0][6:]
        cam    = int(fstem.split('_')[1][-1]) - 1

        out = self.get_pix(sector,cam,ra,dec, d[:,9], ignore_image_buffer=ignore_image_buffer)
        camcol = out[0]
        camrow = out[1]
        ccd    = out[2]
        ccdcol = out[3]
        ccdrow = out[4]
#        ra_trim= out[5]
#        dec_trim=out[6]
#        mag_trim=out[7]
        idx    = out[8]

        d = d[idx]

        ra = ra[idx]
        dec = dec[idx]
        t = t[idx]
        tjd = tjd[idx]

        
        super(TNS,self).__init__(self.keys,
                                 [d[:,0], d[:,1], d[:,2], d[:,3], d[:,4], d[:,5],
                                  ra, dec, t, tjd, d[:,8], d[:,9], 
                                  d[:,10], d[:,11],d[:,12], d[:,13],
                                  camcol, camrow, ccd, ccdcol, ccdrow] )


        self.obj_name = self.name
        self.mag = self.mag.astype(float)
        self.z   = self.z.astype(float)
        self.host_z   = self.host_z.astype(float)
