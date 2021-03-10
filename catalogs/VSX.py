import numpy as np
from astropy.time import Time
import sys
import os
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog

class VSX(Catalog):
    keys = ['name', 
            'ra', 'dec', 
            'var_flag','var_type',
            'mean_mag',
            'bright_mag','bright_passband',
            'faint_mag','faint_passband',
            'epoch_bright_mag',
            'period',
            'camcol','camrow','ccd','ccdcol','ccdrow']
      

    def __init__(self, ifile, ignore_image_buffer=False):
        id_use, var_type,bright_passband, amp_flag, \
            faint_passband = np.genfromtxt(ifile,dtype=str,
                                           usecols=(1,5,7,8,10),
                                           unpack=1,
                                           delimiter=',')

        var_flag, ra, dec, bright_mag, \
            faint_mag, epoch_bright_mag, period = np.genfromtxt(ifile,
                                                                unpack=1,
                                                                usecols=(2,3,4,6,9,11,12),
                                                                delimiter=',')

        idx_amplitudes = [ii for ii in range(len(amp_flag)) if len(amp_flag[ii]) == 1]
        faint_mag[idx_amplitudes] = bright_mag[idx_amplitudes] + faint_mag[idx_amplitudes]

        mean_mag = np.mean([bright_mag,faint_mag],axis=0)

        fstem = os.path.basename(ifile)
        sector = 's' + str(int(fstem.split('_')[1][1:]))
        cam = int(fstem.split('_')[2][-5]) - 1
        
        out = self.get_pix(sector, cam, ra,dec, mean_mag, ignore_image_buffer=ignore_image_buffer)
        
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
        
        super(VSX,self).__init__(self.keys,                                     
                                 [id_use,
                                  ra[idx], dec[idx],
                                  var_flag[idx], var_type[idx],
                                  mean_mag[idx],
                                  bright_mag[idx],bright_passband[idx],
                                  faint_mag[idx],faint_passband[idx],
                                  epoch_bright_mag[idx],period[idx],
                                  camcol, camrow, ccd, ccdcol, ccdrow])
        
        #self.source_id = self.source_id.astype(int)

        self.obj_name = self.name
