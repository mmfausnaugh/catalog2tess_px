import numpy as np
import scipy as sp
from astropy.time import Time
import sys
import os
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog


class MilliquasCsv(Catalog):
    keys = ['name', 'ra', 'dec', 'qso_class', 'rmag', 'bmag', 
            'optical_flag', 'red_psf_flag', 'blue_psf_flag', 'redshift',
            'ref_name', 'ref_redshift', 'qso_prob', 'xray_name','radio_name', 
            'radio_lobe_id1','radio_lobe_id2',
            'camcol','camrow','ccd','ccdcol','ccdrow']
      

    def __init__(self, ifile, ignore_image_buffer=False):
        #d = pd.read_csv(ifile, dtype = {'0': str, '1': float}, sep=',',
        #                engine='python', na_values = sp.isnan, header=None)


        ra,dec,name, \
            qso_class,\
            red_mag, blue_mag,\
            optical_flag,\
            red_psf_flag, blue_psf_flag,\
            redshift, ref_name, ref_redshift,\
            qso_probability, \
            xray_name, radio_name, \
            radio_lobe_id1, radio_lobe_id2  = np.genfromtxt(ifile,
                                                         unpack=1,
                                                         delimiter=',',
                                                         dtype=str)

        ra,dec,red_mag, blue_mag,\
            redshift, qso_probability = np.genfromtxt(ifile,
                                                      unpack=1,
                                                      usecols=(0,1,4,5,9,12),
                                                      delimiter=',')
        
        red_mag[np.isnan(red_mag)] = 999
        blue_mag[np.isnan(blue_mag)] = 999


        fstem = os.path.basename(ifile)
        sector = 's' + str(int(fstem.split('_')[1][1:])) #str and int to strip leading zeros
        cam = int(fstem.split('_')[2][-5]) - 1
        out = self.get_pix(sector,cam,ra,dec, red_mag,
                           ignore_image_buffer=ignore_image_buffer)
        camcol = out[0]
        camrow = out[1]
        ccd    = out[2]
        ccdcol = out[3]
        ccdrow = out[4]
        idx = out[8]
        
        
        super(MilliquasCsv,self).__init__(self.keys,
                                         [name[idx], ra[idx], dec[idx],
                                          qso_class[idx], red_mag[idx], blue_mag[idx],
                                          optical_flag[idx], red_psf_flag[idx], blue_psf_flag[idx],
                                          redshift[idx],
                                          ref_name[idx], ref_redshift[idx],
                                          qso_probability[idx],
                                          xray_name[idx],radio_name[idx], 
                                          radio_lobe_id1[idx], radio_lobe_id2[idx],
                                          camcol, camrow, ccd, ccdcol, ccdrow])
        #all classes must have obj_name attriburte
        self.obj_name = self.name
