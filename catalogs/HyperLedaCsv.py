import scipy as sp
from astropy.time import Time
import sys
import os
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog


class HyperLedaCsv(Catalog):
    keys = ["objname", "pgc", "objtype", "ra", "dec",
            "l", "b", "type", "bar", "ring", "multiple",
            "agnclass", "umag", "e_umag",
            "bmag", "e_bmag", "vmag", "e_vmag",
            "imag", "e_imag", "kmag", "e_kmag",
            "mfir", "velocity", "e_velocity", "A_galaxy",
            "A_internal", "incl",
            "corrected_bmag", "corrected_imag", "corrected_umag", "corrected_BV",
            "modz", "e_modz", "moddist", "e_moddist",
            "modavg", "e_modavg", "mabs",
            
            'camcol','camrow','ccd','ccdcol','ccdrow']
      

    def __init__(self, ifile, ignore_image_buffer=False):
        id_use, objtype, morph_class,bar,ring,multi,agn = sp.genfromtxt(ifile, usecols=(0,2,12,13,14,15,19),unpack=1,dtype=str,delimiter=',')
        d = sp.genfromtxt(ifile,delimiter=',')
        
        ra = d[:,5] * (180.0/12.0) # converting to decimals
        dec = d[:,6]

        fstem = os.path.basename(ifile)
        sector = 's' + str(int(fstem.split('_')[1][1:]))
        cam = int(fstem.split('_')[2][-5]) - 1

        out = self.get_pix(sector, cam, ra,dec, d[:,33], ignore_image_buffer=False)

        camcol = out[0]
        camrow = out[1]
        ccd    = out[2]
        ccdcol = out[3]
        ccdrow = out[4]
        idx = out[8]
        
        # print sector, cam, idx
        # print d
        d = d[idx]
        id_use = id_use[idx]
        objtype = objtype[idx]
        morph_class = morph_class[idx]
        bar = bar[idx]
        ring = ring[idx]
        multi = multi[idx]
        agn = agn[idx]
        # print d, ccd
        ra = ra[idx]
        dec = dec[idx]
        
        
        super(HyperLedaCsv,self).__init__(self.keys,
                                         [id_use, d[:,1], objtype, d[:,5], d[:,6],
                                          d[:,7], d[:,8], 
                                          morph_class, bar, ring, multi, agn,
                                          d[:,27], d[:,28],
                                          d[:,29], d[:,30], d[:,31], d[:,32],
                                          d[:,33], d[:,34], d[:,35], d[:,36],
                                          d[:,39], d[:,52], d[:,53], d[:,54],
                                          d[:,55], d[:,56], d[:,59], d[:,60],
                                          d[:,61], d[:,62], d[:,74], d[:,75],
                                          d[:,76], d[:,77], d[:,78], d[:,79], d[:,80], 

                                          camcol, camrow, ccd, ccdcol, ccdrow])
        #classes from all catalogs must have obj_name
        self.obj_name = self.objname
