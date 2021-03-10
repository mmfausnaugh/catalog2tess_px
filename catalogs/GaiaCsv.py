import scipy as sp
from astropy.time import Time
import sys
import os
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog

class GaiaCsv(Catalog):
    keys = ['source_id', 'ra', 'ra_error', 'dec', 'dec_error', 'parallax', 
            'parallax_error', 'astrometric_gof_al', 'phot_g_mean_mag', 
            'phot_bp_mean_mag', 'phot_rp_mean_mag', 'bp_rpl', 'radial_velocity',
            'radial_velocity_error',  'ecl_lon', 'ecl_lat',
            'teff_val', 'a_g_val', 'e_bp_min_rp_val', 'radius_val', 'lum_val',            
            'camcol','camrow','ccd','ccdcol','ccdrow']
      

    def __init__(self, ifile, ignore_image_buffer=False):
        id_use = sp.genfromtxt(ifile,dtype=sp.uint64,usecols=(0))

        d = sp.genfromtxt(ifile,usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))

        ra = d[:,0]
        dec = d[:,2]

        fstem = os.path.basename(ifile)
        sector = 's' + str(int(fstem.split('_')[2][1:]))
        cam = int(fstem.split('_')[3][-5]) - 1
        
        out = self.get_pix(sector, cam, ra,dec, d[:,9], ignore_image_buffer=ignore_image_buffer)
        
        # any errors to correct? 
        camcol = out[0]
        camrow = out[1]
        ccd    = out[2]
        ccdcol = out[3]
        ccdrow = out[4]
        idx = out[8]
        
        # print sector, cam, idx
        # print d
        id_use = id_use[idx]
        d = d[idx]
        # print d, ccd
        ra = ra[idx]
        dec = dec[idx]
        
        super(GaiaCsv,self).__init__(self.keys,                                     
                                         [id_use,
                                          d[:,0], d[:,1], d[:,2], d[:,3], d[:,4],
                                          d[:,5], d[:,6], d[:,7], d[:,8], d[:,9], d[:,10],
                                          d[:,11], d[:,12], d[:,13], d[:,14], d[:,15],
                                          d[:,16], d[:,17], d[:,18], d[:,19],
                                          camcol, camrow, ccd, ccdcol, ccdrow])

        #self.source_id = self.source_id.astype(int)

        self.obj_name = self.source_id
