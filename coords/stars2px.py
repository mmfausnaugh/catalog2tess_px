#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  3 22:32:51 2018

@author: cjburke
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
from astropy.coordinates import SkyCoord
from astropy import units as u
import sys
import time
import datetime
import json
try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve
try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib  


class Levine_FPG():
    """Al Levine Focal Plane Geometry Methods
        Translated from starspx6.c
        INPUT:
            sc_ra_dec_roll = numpy array of the SpaceCraft boresite (sc Z-axis)
            ra, dec, and roll [deg]
            The roll angle is in RA, Dec space clockwise relative to the celestial
            pole.  roll angle = 0 [deg] implies space craft X-axis points N celestial (increasing dec)
            roll angle = 90 [deg] implies sc X-axis points towards increasing/decreasing (?) RA
        *** In practice there is a separate fpg file for each of the four cameras ***
        rmat1[3,3] = is the rotation matrix from ra&dec to spacecraft boresite coords
        rmat4[NCAM,3,3] - is the rotation matrix from ra&dec to NCAM coords
    """
    parm_dict_list = [{}, {}, {}, {}]
    NCAM = 4 # Number of Cameras
    NCCD = 4 # Number of CCDs per Camera
    
    def __init__(self, sc_ra_dec_roll=None, fpg_file_list=None):
        self.eulcam = np.zeros((self.NCAM,3), dtype=np.double)
        self.optcon = np.zeros((self.NCAM,6), dtype=np.double)
        self.ccdxy0 = np.zeros((self.NCAM, self.NCCD, 2), dtype=np.double)
        self.pixsz = np.zeros((self.NCAM, self.NCCD, 2), dtype=np.double)
        self.ccdang = np.zeros((self.NCAM, self.NCCD), dtype=np.double)
        self.ccdtilt = np.zeros((self.NCAM, self.NCCD, 2), dtype=np.double)
        self.asymang = np.zeros((self.NCAM,), dtype=np.double)
        self.asymfac = np.zeros((self.NCAM,), dtype=np.double)
        self.rmat1 = np.zeros((3,3), dtype=np.double)
        self.rmat4 = np.zeros((self.NCAM,3,3), dtype=np.double)
        self.havePointing = False
        # Read in the fpg parameter files        
        self.read_all_levine_fpg_files(fpg_file_list)
        # Generate rotation matrices if ra dec and roll values given
        if not sc_ra_dec_roll is None:
            # go from sky to spacecraft
            self.sky_to_sc_mat(sc_ra_dec_roll)
            # Go from spacecraft to each camera's coords
            for icam in range(self.NCAM):
                cureul = self.eulcam[icam,:]
                rmat2 = self.sc_to_cam_mat(cureul)
                self.rmat4[icam] = np.matmul(rmat2, self.rmat1)
            self.havePointing = True
        
    def read_all_levine_fpg_files(self, fpg_file_list=None):                
        default_fpg_file_list = ['fpg_pars.txt-', \
                                 'fpg_pars.txt-', \
                                 'fpg_pars.txt-', \
                                 'fpg_pars.txt-']
        # For each camera read in the separate fpg parameter file
        for icam in range(self.NCAM):
            if fpg_file_list == None:
                fpg_file = default_fpg_file_list[icam]
            else:
                fpg_file = fpg_file_list[icam]
            self.read_levine_fpg_file(icam, fpg_file)
        # We now have parameters for all 4 cameras in the parm_dict_list
        # parse the dictionary values into the working numpy arrays
        for icam in range(self.NCAM):
            pd = self.parm_dict_list[icam]
            self.eulcam[icam][0] = pd['ang1_cam1']
            self.eulcam[icam][1] = pd['ang2_cam1']
            self.eulcam[icam][2] = pd['ang3_cam1']
            self.optcon[icam][0] = pd['fl_cam1']
            self.optcon[icam][1] = pd['opt_coef1_cam1']	
            self.optcon[icam][2] = pd['opt_coef2_cam1']	
            self.optcon[icam][3] = pd['opt_coef3_cam1']	
            self.optcon[icam][4] = pd['opt_coef4_cam1']	
            self.optcon[icam][5] = pd['opt_coef5_cam1']
            self.asymang[icam] = pd['asymang_cam1']
            self.asymfac[icam] = pd['asymfac_cam1']
            for iccd in range(self.NCCD):
                self.ccdxy0[icam][iccd][0] = pd['x0_ccd{0:1d}_cam1'.format(iccd+1)]
                self.ccdxy0[icam][iccd][1] = pd['y0_ccd{0:1d}_cam1'.format(iccd+1)]
                self.pixsz[icam][iccd][0] = pd['pix_x_ccd{0:1d}_cam1'.format(iccd+1)]
                self.pixsz[icam][iccd][1] = pd['pix_y_ccd{0:1d}_cam1'.format(iccd+1)]
                self.ccdang[icam][iccd] = pd['ang_ccd{0:1d}_cam1'.format(iccd+1)]
                self.ccdtilt[icam][iccd][0] = pd['tilt_x_ccd{0:1d}_cam1'.format(iccd+1)]
                self.ccdtilt[icam][iccd][1] = pd['tilt_y_ccd{0:1d}_cam1'.format(iccd+1)]
        
            
    def read_levine_fpg_file(self, icam, fpg_file):
        gotParm = False
        parm_dict = {}
        if os.path.isfile(fpg_file):        
            try:    
                fpin = open(fpg_file, 'r')
                # Read in parameters 
                dtypeseq = ['U20','i4','f16']
                dataBlock = np.genfromtxt(fpin, dtype=dtypeseq)
                parm_keys = dataBlock['f0']
                parm_fitted_flags = dataBlock['f1']
                parm_values = dataBlock['f2']
                # Now build dictionary of the parameters
                for i in range(len(parm_keys)):
                    parm_dict[parm_keys[i]] = parm_values[i]
                self.parm_dict_list[icam] = parm_dict
                gotParm = True
                print('Successful Focal Plane Geometry Read From {0}'.format(fpg_file))
            except:
                print('Could not open {0}!  Using Hard-coded Focal Plane Geometry from Levine_FPG read_levine_fpg_file()'.format(fpg_file))
        # If anything goes wrong with reading in parameters revert to hard coded version
        # or file was never given and default_fpg_file does not exist
        if not gotParm:
            #print('Using Hard-coded Focal Plane Geometry from Levine_FPG read_levine_fpg_file')
            # *** For now this hard code is just a filler need to actually fill in values for all cameras separately
            if icam == 0:
                parm_dict = {'ang1_cam1': 0.116992, \
                             'ang2_cam1': -36.029617, \
                             'ang3_cam1': 90.033899, \
                             'fl_cam1': 145.846659, \
                             'opt_coef1_cam1': 1.0000014, \
                             'opt_coef2_cam1': 0.28174612, \
                             'opt_coef3_cam1': -0.59667259, \
                             'opt_coef4_cam1': 9.17151267, \
                             'opt_coef5_cam1': -4.36928235, \
                             'asymang_cam1': 0.0, \
                             'asymfac_cam1': 1.0, \
                             'x0_ccd1_cam1': 31.525684, \
                             'y0_ccd1_cam1': 31.56166, \
                             'pix_x_ccd1_cam1': 0.015, \
                             'pix_y_ccd1_cam1': 0.015, \
                             'ang_ccd1_cam1': 179.979265, \
                             'tilt_x_ccd1_cam1': 0.0, \
                             'tilt_y_ccd1_cam1': 0.0, \
                             'x0_ccd2_cam1': -0.948971, \
                             'y0_ccd2_cam1': 31.545675, \
                             'pix_x_ccd2_cam1': 0.015, \
                             'pix_y_ccd2_cam1': 0.015, \
                             'ang_ccd2_cam1': 180.0, \
                             'tilt_x_ccd2_cam1': 0.0, \
                             'tilt_y_ccd2_cam1': 0.0, \
                             'x0_ccd3_cam1': -31.695077, \
                             'y0_ccd3_cam1': -31.423998, \
                             'pix_x_ccd3_cam1': 0.015, \
                             'pix_y_ccd3_cam1': 0.015, \
                             'ang_ccd3_cam1': -0.024352, \
                             'tilt_x_ccd3_cam1': 0.0, \
                             'tilt_y_ccd3_cam1': 0.0, \
                             'x0_ccd4_cam1': 0.785759, \
                             'y0_ccd4_cam1': -31.442781, \
                             'pix_x_ccd4_cam1': 0.015, \
                             'pix_y_ccd4_cam1': 0.015, \
                             'ang_ccd4_cam1': 0.002997, \
                             'tilt_x_ccd4_cam1': 0.0, \
                             'tilt_y_ccd4_cam1': 0.0}
            if icam == 1:
                parm_dict =  {'ang1_cam1': -0.141259, \
                              'ang2_cam1': -12.030006, \
                              'ang3_cam1': 90.007001, \
                              'fl_cam1': 145.881826, \
                              'opt_coef1_cam1': 1.0000014, \
                              'opt_coef2_cam1': 0.28174612, \
                              'opt_coef3_cam1': -0.59667259, \
                              'opt_coef4_cam1': 9.17151267, \
                              'opt_coef5_cam1': -4.36928235, \
                              'asymang_cam1': 0.0, \
                              'asymfac_cam1': 1.0, \
                              'x0_ccd1_cam1': 31.621308, \
                              'y0_ccd1_cam1': 31.478634, \
                              'pix_x_ccd1_cam1': 0.015, \
                              'pix_y_ccd1_cam1': 0.015, \
                              'ang_ccd1_cam1': 180.008334, \
                              'tilt_x_ccd1_cam1': 0.0, \
                              'tilt_y_ccd1_cam1': 0.0, \
                              'x0_ccd2_cam1': -0.85576, \
                              'y0_ccd2_cam1': 31.500341, \
                              'pix_x_ccd2_cam1': 0.015, \
                              'pix_y_ccd2_cam1': 0.015, \
                              'ang_ccd2_cam1': 180.0, \
                              'tilt_x_ccd2_cam1': 0.0, \
                              'tilt_y_ccd2_cam1': 0.0, \
                              'x0_ccd3_cam1': -31.571863, \
                              'y0_ccd3_cam1': -31.536262, \
                              'pix_x_ccd3_cam1': 0.015, \
                              'pix_y_ccd3_cam1': 0.015, \
                              'ang_ccd3_cam1': -0.006384,\
                              'tilt_x_ccd3_cam1': 0.0, \
                              'tilt_y_ccd3_cam1': 0.0, \
                              'x0_ccd4_cam1': 0.889169, \
                              'y0_ccd4_cam1': -31.543301, \
                              'pix_x_ccd4_cam1': 0.015, \
                              'pix_y_ccd4_cam1': 0.015, \
                              'ang_ccd4_cam1': -0.015664, \
                              'tilt_x_ccd4_cam1': 0.0, \
                              'tilt_y_ccd4_cam1': 0.0}
            if icam == 2:
                parm_dict = {'ang1_cam1': 0.109909, \
                             'ang2_cam1': 12.000815, \
                             'ang3_cam1': -89.93835, \
                             'fl_cam1': 145.896267, \
                             'opt_coef1_cam1': 1.0000014, \
                             'opt_coef2_cam1': 0.28174612, \
                             'opt_coef3_cam1': -0.59667259, \
                             'opt_coef4_cam1': 9.17151267, \
                             'opt_coef5_cam1': -4.36928235, \
                             'asymang_cam1': 0.0, \
                             'asymfac_cam1': 1.0, \
                             'x0_ccd1_cam1': 31.593984, \
                             'y0_ccd1_cam1': 31.392635, \
                             'pix_x_ccd1_cam1': 0.015, \
                             'pix_y_ccd1_cam1': 0.015, \
                             'ang_ccd1_cam1': 179.994458,\
                             'tilt_x_ccd1_cam1': 0.0, \
                             'tilt_y_ccd1_cam1': 0.0, \
                             'x0_ccd2_cam1': -0.848621, \
                             'y0_ccd2_cam1': 31.404205, \
                             'pix_x_ccd2_cam1': 0.015, \
                             'pix_y_ccd2_cam1': 0.015, \
                             'ang_ccd2_cam1': 180.0, \
                             'tilt_x_ccd2_cam1': 0.0, \
                             'tilt_y_ccd2_cam1': 0.0, \
                             'x0_ccd3_cam1': -31.562291, \
                             'y0_ccd3_cam1': -31.624249, \
                             'pix_x_ccd3_cam1': 0.015, \
                             'pix_y_ccd3_cam1': 0.015, \
                             'ang_ccd3_cam1': 0.002481, \
                             'tilt_x_ccd3_cam1': 0.0, \
                             'tilt_y_ccd3_cam1': 0.0, \
                             'x0_ccd4_cam1': 0.876557, \
                             'y0_ccd4_cam1': -31.586095, \
                             'pix_x_ccd4_cam1': 0.015, \
                             'pix_y_ccd4_cam1': 0.015, \
                             'ang_ccd4_cam1': -0.003476, \
                             'tilt_x_ccd4_cam1': 0.0, \
                             'tilt_y_ccd4_cam1': 0.0}
            if icam == 3:
                parm_dict = {'ang1_cam1': 0.025385, \
                             'ang2_cam1': 35.963801, \
                             'ang3_cam1': -89.978603, \
                             'fl_cam1': 145.928875, \
                             'opt_coef1_cam1': 1.0000014, \
                             'opt_coef2_cam1': 0.28174612, \
                             'opt_coef3_cam1': -0.59667259, \
                             'opt_coef4_cam1': 9.17151267, \
                             'opt_coef5_cam1': -4.36928235, \
                             'asymang_cam1': 0.0, \
                             'asymfac_cam1': 1.0, \
                             'x0_ccd1_cam1': 31.56822, \
                             'y0_ccd1_cam1': 31.285512, \
                             'pix_x_ccd1_cam1': 0.015, \
                             'pix_y_ccd1_cam1': 0.015, \
                             'ang_ccd1_cam1': 179.970123, \
                             'tilt_x_ccd1_cam1': 0.0, \
                             'tilt_y_ccd1_cam1': 0.0, \
                             'x0_ccd2_cam1': -0.892286, \
                             'y0_ccd2_cam1': 31.332011, \
                             'pix_x_ccd2_cam1': 0.015, \
                             'pix_y_ccd2_cam1': 0.015, \
                             'ang_ccd2_cam1': 180.0, \
                             'tilt_x_ccd2_cam1': 0.0, \
                             'tilt_y_ccd2_cam1': 0.0, \
                             'x0_ccd3_cam1': -31.632113, \
                             'y0_ccd3_cam1': -31.742848, \
                             'pix_x_ccd3_cam1': 0.015, \
                             'pix_y_ccd3_cam1': 0.015, \
                             'ang_ccd3_cam1': -0.025901, \
                             'tilt_x_ccd3_cam1': 0.0, \
                             'tilt_y_ccd3_cam1': 0.0, \
                             'x0_ccd4_cam1': 0.81764, \
                             'y0_ccd4_cam1': -31.754969, \
                             'pix_x_ccd4_cam1': 0.015, \
                             'pix_y_ccd4_cam1': 0.015, \
                             'ang_ccd4_cam1': -0.024328, \
                             'tilt_x_ccd4_cam1': 0.0, \
                             'tilt_y_ccd4_cam1': 0.0}
            self.parm_dict_list[icam] = parm_dict

    def sky_to_sc_mat(self, sc_ra_dec_roll):
        """Calculate the rotation matrix that will convert a vector in ra&dec
            into the spacecraft boresite frame
        """
        deg2rad = np.pi / 180.0
        # Define the 3 euler angles of rotation
        xeul = np.zeros((3,), dtype=np.double)
        xeul[0] = deg2rad * sc_ra_dec_roll[0]
        xeul[1] = np.pi/2.0 - deg2rad*sc_ra_dec_roll[1]
        xeul[2] = deg2rad * sc_ra_dec_roll[2] + np.pi
        # Generate the rotation matrix from the 3 euler angles
        self.rmat1 = self.eulerm323(xeul)

    def sc_to_cam_mat(self, eul):
        """Calculate the rotation matrix that will convert a vector in spacecraft
            into the a camera's coords
        """
        deg2rad = np.pi / 180.0
        # Generate the rotation matrix from the 3 euler angles
        xeul = deg2rad * eul
        return self.eulerm323(xeul)
        
    def eulerm323(self, eul):
        mat1 = self.rotm1(2, eul[0])
        mat2 = self.rotm1(1, eul[1])
        mata = np.matmul(mat2, mat1)
        mat1 = self.rotm1(2, eul[2])
        rmat = np.matmul(mat1, mata)
        return rmat
        
    def rotm1(self, ax, ang):
        mat = np.zeros((3,3), dtype=np.double)
        n1 = ax
        n2 = np.mod((n1+1), 3)
        n3 = np.mod((n2+1), 3)
        sinang = np.sin(ang)
        cosang = np.cos(ang) 
        mat[n1][n1] = 1.0
        mat[n2][n2] = cosang
        mat[n3][n3] = cosang
        mat[n2][n3] = sinang
        mat[n3][n2] = -sinang
        return mat

    def sphereToCart(self, ras, decs):
        """ Convert 3d spherical coordinates to cartesian
        """
        deg2rad = np.pi / 180.0
        rarads = deg2rad * ras
        decrads = deg2rad * decs
        sinras = np.sin(rarads)
        cosras = np.cos(rarads)
        sindecs = np.sin(decrads)
        cosdecs = np.cos(decrads)
        vec0s = cosras * cosdecs
        vec1s = sinras * cosdecs
        vec2s = sindecs
        return vec0s, vec1s, vec2s

    def cartToSphere(self, vec):
        ra = 0.0
        dec = 0.0
        norm = np.sqrt(np.sum(vec*vec))
        if (norm > 0.0):
            dec = np.arcsin(vec[2] / norm)
            if (not vec[0] == 0.0) or (not vec[1] == 0.0):
                ra = np.arctan2(vec[1], vec[0])
                ra = np.mod(ra, 2.0*np.pi)
        return ra, dec

    def star_in_fov(self, lng, lat):
        deg2rad = np.pi / 180.0
        inView = False
        if lat > 73.0:
            vec0, vec1, vec2 = self.sphereToCart(lng, lat)
            vec = np.array([vec0, vec1, vec2], dtype=np.double)
            norm = np.sqrt(np.sum(vec*vec))
            if norm > 0.0:
                vec = vec / norm
                xlen = np.abs(np.arctan(vec[0]/vec[2]))
                ylen = np.abs(np.arctan(vec[1]/vec[2]))
                if (xlen <= (12.5 * deg2rad)) and (ylen <= (12.5 * deg2rad)):
                    inView = True
        return inView

    def optics_fp(self, icam, lng_deg, lat_deg):
        deg2rad = np.pi / 180.0
        thetar = np.pi / 2.0 - (lat_deg * deg2rad)
        tanth = np.tan(thetar)
        cphi = np.cos(deg2rad*lng_deg)
        sphi = np.sin(deg2rad*lng_deg)
        rfp0 = self.optcon[icam][0]*tanth
        noptcon = len(self.optcon[icam])
        ii = np.arange(1, noptcon)
        rfp = np.sum(self.optcon[icam][1:] * np.power(tanth, 2.0*(ii-1)))
        xytmp = np.zeros((2,), dtype=np.double)
        xytmp[0] = -cphi*rfp0*rfp
        xytmp[1] = -sphi*rfp0*rfp
        return self.make_az_asym(icam, xytmp)
        
    def make_az_asym(self, icam, xy):
        xyp = self.xyrotate(self.asymang[icam], xy)
        xypa = np.zeros_like(xyp)
        xypa[0] = self.asymfac[icam] * xyp[0]
        xypa[1] = xyp[1]
        xyout = self.xyrotate(-self.asymang[icam], xypa)
        return xyout
        
    def xyrotate(self, angle_deg, xin):
        deg2rad = np.pi / 180.0
        ca = np.cos(deg2rad * angle_deg)
        sa = np.sin(deg2rad * angle_deg)
        xyout = np.zeros_like(xin)
        xyout[0] = ca*xin[0] + sa*xin[1]
        xyout[1] = -sa*xin[0] + ca*xin[1]
        return xyout

    def mm_to_pix(self, icam, xy):
        """Convert focal plane to pixel location also need to add in the
            auxillary pixels added into FFIs
        """
        CCDWD_T=2048
        CCDHT_T=2058
        ROWA=44
        ROWB=44
        COLDK_T=20
        xya = np.copy(xy)
        xyb = np.zeros_like(xya)
        ccdpx = np.zeros_like(xya)
        fitpx = np.zeros_like(xya)
        if xya[0] >= 0.0:
            if xya[1] >= 0.0:
                iccd = 0
                xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
                xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
                xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
                ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
                ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
                fitpx[0] = (CCDWD_T - ccdpx[0]) + CCDWD_T + 2*ROWA + ROWB - 1.0
                fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0
            else:
                iccd = 3
                xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
                xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
                xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
                ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
                ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
                fitpx[0] = ccdpx[0] + CCDWD_T + 2*ROWA + ROWB
                fitpx[1] = ccdpx[1]
        else:
            if xya[1] >= 0.0:
                iccd = 1
                xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
                xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
                xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
                ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
                ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
                fitpx[0] = (CCDWD_T - ccdpx[0]) + ROWA - 1.0
                fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0
            else:
                iccd = 2
                xyb[0] = xya[0] - self.ccdxy0[icam][iccd][0]
                xyb[1] = xya[1] - self.ccdxy0[icam][iccd][1]
                xyccd = self.xyrotate(self.ccdang[icam][iccd], xyb)
                ccdpx[0] = (xyccd[0] / self.pixsz[icam][iccd][0]) - 0.5
                ccdpx[1] = (xyccd[1] / self.pixsz[icam][iccd][1]) - 0.5
                fitpx[0] = ccdpx[0] + ROWA
                fitpx[1] = ccdpx[1]
                
        return iccd, ccdpx, fitpx
        
    def radec2pix(self, ras, decs):
        """ After the rotation matrices are defined to the actual
            ra and dec to pixel coords mapping
        """
        nStar = len(ras)
        inCamera = np.zeros((nStar,), dtype=np.int)
        ccdNum = np.zeros((nStar,), dtype=np.int)
        fitsxpos = np.zeros((nStar,), dtype=np.double)
        fitsypos = np.zeros((nStar,), dtype=np.double)
        ccdxpos = np.zeros((nStar,), dtype=np.double)
        ccdypos = np.zeros((nStar,), dtype=np.double)
        deg2rad = np.pi / 180.0
        if self.havePointing == True:
            # Convert ra and dec spherical coords to cartesian
            vec0s, vec1s, vec2s = self.sphereToCart(ras, decs)
            for i in range(nStar):
                curVec = np.array([vec0s[i], vec1s[i], vec2s[i]], dtype=np.double)
                # Find the new vector in all cameras
                for j in range(self.NCAM):
                    # Do the rotation from ra dec coords to camera coords
                    camVec = np.matmul(self.rmat4[j], curVec)
                    # Get the longitude and latitude of camera coords position
                    lng, lat = self.cartToSphere(camVec)
                    lng = lng / deg2rad
                    lat = lat / deg2rad
                    if self.star_in_fov(lng, lat):
                        # Get the xy focal plane position in mm
                        xyfp = self.optics_fp(j, lng, lat)
                        # Convert mm to pixels
                        iccd, ccdpx, fitpx = self.mm_to_pix(j, xyfp)
                        inCamera[i] = j+1 # Als code is base 0 convert to base 1
                        ccdNum[i] = iccd+1 # ""
                        fitsxpos[i] = fitpx[0]
                        fitsypos[i] = fitpx[1]
                        ccdxpos[i] = ccdpx[0]
                        ccdypos[i] = ccdpx[1]
        else:
            print('Spacecraft Pointing Not specified!')
        
        return inCamera, ccdNum, fitsxpos, fitsypos, ccdxpos, ccdypos
        
class TESS_Spacecraft_Pointing_Data:
    #Hard coded spacecraft pointings by Sector
    sectors = np.arange(1,14, dtype=np.int)
    ras = np.array([352.6844,16.5571,36.3138,55.0070,73.5382, \
                    92.0096,110.2559,128.1156,145.9071,\
                    165.0475,189.1247,229.5885,298.6671], dtype=np.float)
    decs = np.array([ -64.8531,-54.0160,-44.2590,-36.6420, -31.9349, \
                     -30.5839,-32.6344,-37.7370,-45.3044, \
                     -54.8165,-65.5369,-75.1256,-76.3281], dtype=np.float)
    rolls = np.array([-137.8468,-139.5665,-146.9616,-157.1698,-168.9483, \
                      178.6367, 166.4476,155.3091,145.9163,\
                      139.1724,138.0761,153.9773,-161.0622], dtype=np.float)
    # Which hemisphere is pointing; +1==South ; -1==North
    hemis = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1], dtype=np.int)
    # Actual observed pointings versus TBD/predicted pointings
    obsPoint = np.array([1,1,1,0,0,0,0,0,0,0,0,0,0], dtype=np.int)
    
    def __init__(self, trySector=None, fpgParmFileList=None):
        # Convert S/C boresite pointings to ecliptic coords for each camera
        # If trySector is set only keep the single requested sector
        if not trySector is None:
            idx = np.where(self.sectors == trySector)[0]
            self.sectors = self.sectors[idx]
            self.ras = self.ras[idx]
            self.decs = self.decs[idx]
            self.rolls = self.rolls[idx]
            self.hemis = self.hemis[idx]
            self.obsPoint = self.obsPoint[idx]
        nPoints = len(self.sectors)
        self.camEclipLat = np.zeros((4, nPoints), dtype=np.float)
        self.camEclipLon = np.zeros((4, nPoints), dtype=np.float)
        self.camRa = np.zeros((4, nPoints), dtype=np.float)
        self.camDec = np.zeros((4, nPoints), dtype=np.float)
        # Convert S/C boresite ra and dec to ecliptic coords
        for iPnt in range(nPoints):
            scCoords = SkyCoord(self.ras[iPnt], self.decs[iPnt], unit='deg')
            scEclipCoords = scCoords.transform_to(frame='barycentrictrueecliptic')
            scEclipLat = scEclipCoords.lat.deg
            scEclipLong = scEclipCoords.lon.deg
            # Camera ecliptic latitude offset sign depends on which hemisphere
            camLatOffsets =  np.array([36.0, 12.0, -12.0, -36.0])*self.hemis[iPnt]
            for iCam in range(4):
                self.camEclipLat[iCam,iPnt] = scEclipLat + camLatOffsets[iCam]
                self.camEclipLon[iCam,iPnt] = scEclipLong
                # Just for testing convert camera ecliptic coords to ra and dec
                # compare to published values
                #cc = SkyCoord(self.camEclipLon[iCam,iPnt], \
                #              np.max([-89.99999,self.camEclipLat[iCam,iPnt]]), \
                #              frame='barycentrictrueecliptic', unit='deg')
                #cc = cc.transform_to(frame='icrs')
                #print('{:d} {:d} {:f} {:f} {:f} {:f}'.format(self.sectors[iPnt],iCam+1,\
                #         self.camEclipLon[iCam,iPnt], self.camEclipLat[iCam,iPnt], \
                #         cc.ra.deg, cc.dec.deg))
        # For every pointing make a Levine pointing class object
        self.fpgObjs = []
        fpg_file_list=None
        if not fpgParmFileList is None:
            fpg_file_list=fpgParmFileList
        for iPnt in range(nPoints):
            sc_ra_dec_roll =  np.array([self.ras[iPnt], self.decs[iPnt], self.rolls[iPnt]])
            self.fpgObjs.append(Levine_FPG(sc_ra_dec_roll, fpg_file_list=fpg_file_list))
        
                

class target_info:
    def __init__(self):
        self.ticid = 0
        self.ra = 0.0
        self.dec = 0.0
        self.eclipLong = 0.0
        self.eclipLat = 0.0
        self.sectors = np.array([], dtype=np.int)
        self.onSiliconFlag = np.array([], dtype=np.int)
        self.possibleOnSiliconFlag = np.array([], dtype=np.int)
        self.cameras = np.array([], dtype=np.int)
        self.xpxs = np.array([], dtype=np.float)
        self.ypxs = np.array([], dtype=np.float)

def make_target_objects(tic, ra, dec):
    starList = []
    for i, curTic in enumerate(tic):
        curRa = ra[i]
        curDec = dec[i]
        # instantiate target object
        curTarg = target_info()
        curTarg.ra = curRa
        curTarg.dec = curDec
        curTarg.ticid = tic[i]
        # Convert ra and dec coords to ecliptic
        planCoords = SkyCoord(curRa, curDec, unit='deg')
        planEclipCoords = planCoords.transform_to(frame='barycentrictrueecliptic')
        curTarg.eclipLat = planEclipCoords.lat.deg
        curTarg.eclipLong = planEclipCoords.lon.deg
        starList.append(curTarg)
    return starList

def doRoughPosition(targinfo, scinfo):
    # Return the combinations of sector and detectors that can possibly observe
    #  target
    # go through each position in the spacecraft info class
    eLat = targinfo.eclipLat
    eLon = targinfo.eclipLong
    FOVDeg = 12.5
    nPoints = len(scinfo.sectors)
    for iPnt in range(nPoints):
        # First check target is in correct ecliptic hemisphere of pointing
        targEclipCoords = SkyCoord(eLon, eLat, \
                                   frame='barycentrictrueecliptic', unit='deg')
        if not np.sign(eLat) == scinfo.hemis[iPnt]:
            for iCam in range(4):
                camLon = scinfo.camEclipLon[iCam,iPnt]
                camLat = scinfo.camEclipLat[iCam,iPnt]
                camCenter = SkyCoord(camLon, np.max([-89.99999,camLat]), unit='deg', frame='barycentrictrueecliptic')
                posAngs = camCenter.position_angle(targEclipCoords)
                seps = camCenter.separation(targEclipCoords)
                xseps = (seps * np.cos(posAngs.rad)).deg
                yseps = (seps * np.sin(posAngs.rad)).deg
                # Check for potentially on silicon
                if np.abs(xseps)<FOVDeg and np.abs(yseps) < FOVDeg:
                    # Append potential pointing camera combo to targets list
                    targinfo.sectors = np.append(targinfo.sectors, scinfo.sectors[iPnt])
                    targinfo.onSiliconFlag = np.append(targinfo.onSiliconFlag, 0)
                    targinfo.possibleOnSiliconFlag = np.append(targinfo.possibleOnSiliconFlag, 1)
                    targinfo.cameras = np.append(targinfo.cameras, iCam+1)
                    targinfo.xpxs = np.append(targinfo.xpxs, 0.0)
                    targinfo.ypxs = np.append(targinfo.ypxs, 0.0)
    return targinfo

## [Mast Query]
def mastQuery(request):

    server='mast.stsci.edu'

    # Grab Python Version 
    version = ".".join(map(str, sys.version_info[:3]))

    # Create Http Header Variables
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept": "text/plain",
               "User-agent":"python-requests/"+version}

    # Encoding the request as a json string
    requestString = json.dumps(request)
    requestString = urlencode(requestString)
    
    # opening the https connection
    conn = httplib.HTTPSConnection(server)

    # Making the query
    conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)

    # Getting the response
    resp = conn.getresponse()
    head = resp.getheaders()
    content = resp.read().decode('utf-8')

    # Close the https connection
    conn.close()

    return head,content
## [Mast Query]

def fileOutputHeader(fp, fpgParmFileList=None):
    # output a header to the file
    fp.write('# stars2px.py - Convert Target RA and Dec to TESS spacecraft pixel coordinates\n')
    fp.write('#   Original starspx.c by Alan Levine (MIT Kavli Institute)\n')
    fp.write('#   Python translation by Christopher Burke (MIT Kavli Institute)\n')
    fp.write('# Output columns Pipe Delimited; 16 header lines\n')
    fp.write('# File Creation: {:}\n'.format(datetime.datetime.now()))
    if fpgParmFileList is None:
        fp.write('# FPG Model Default\n')
    else:
        fp.write('# FPG Model {:s} {:s} {:s} {:s}\n'.format(fpgParmFileList[0], \
                 fpgParmFileList[1],fpgParmFileList[2],fpgParmFileList[3]))
    fp.write('# 1 [int] Input TIC ID\n')
    fp.write('# 2 [degree] Input or MAST TIC query target RA\n')
    fp.write('# 3 [degree] Input or MAST TIC query target Dec\n')
    fp.write('# 4 [degree] Ecliptic Longitude\n')
    fp.write('# 5 [degree] Ecliptic Latitude\n')
    fp.write('# 6 [int] - Observing Sector number for target\n')
    fp.write('# 7 [int] - Camera number for target\n')
    fp.write('# 8 [int] - Detector number for target\n')
    fp.write('# 9 [float] - Row pixel location for target\n')
    fp.write('# 10 [float] - Column pixel location for target\n')
             
    
if __name__ == '__main__':
    # Parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--ticId", type=int, \
                        help="TIC Id [int] for MAST coordinate query.  MUST Be Online For this option to Work!")
    parser.add_argument("-c", "--coord", type=float, nargs=2, \
                        help="RA and Dec of target [deg]")
    parser.add_argument("-f", "--inputFile", type=argparse.FileType('r'), \
                        help="Filename for input Target TIC [int]; RA[deg]; Dec[dec]; in white space delimited text file Column 1, 2, and 3 respectively")
    parser.add_argument("-o", "--outputFile", type=argparse.FileType('w'), \
                        help="Optional filename for output.  Default is output to stdout ")
    parser.add_argument("-s", "--sector", type=int, choices=range(1,14),\
                        help="Search a single sector Number [int]")
    parser.add_argument("-x", "--combinedFits", action='store_true', \
                        help="Output detector pixel coordinates for the 'Big' multi-detector combined fits file format")
    parser.add_argument("-fpg", "--fpgParameterFiles", nargs=4,\
                        help="Instead of default focal plane geometry parameters, list the 4 filenames for the fpg files to use.  Expects files in Al's format in camera numerical order")
    args = parser.parse_args()
 
    # At least one Mode -t -c -f must have been specified
    if (args.ticId is None) and (args.coord is None) and (args.inputFile is None):
        print('You must specify one and only one mode -t, -c, -f')
        print('`python stars2px.py -h\' for help')
        sys.exit(1)
        
    # Do single coords first
    if not (args.coord is None):
        nTarg = 1
        starTics = np.array([0], dtype=np.int32)
        starRas = np.array([args.coord[0]], dtype=np.float)
        starDecs = np.array([args.coord[1]], dtype=np.float)
    else:
        if not (args.inputFile is None): # Check for input file list next
            # Read in star positions in input
            # Now go through stars
            starFile = args.inputFile
            dataBlock = np.genfromtxt(starFile, dtype=['i4','f8','f8'])
            starTics = np.atleast_1d(dataBlock['f0'])
            starRas = np.atleast_1d(dataBlock['f1'])
            starDecs = np.atleast_1d(dataBlock['f2'])
        else:
            # Must have requested MAST query with TIC ID
            # Make a list of TICs using strings
            starTics = np.array([args.ticId], dtype=np.int32)
            ticStringList = ['{0:d}'.format(x) for x in starTics]    
            # Setup mast query
            request = {'service':'Mast.Catalogs.Filtered.Tic', \
               'params':{'columns':'*', 'filters':[{ \
                        'paramName':'ID', 'values':ticStringList}]}, \
                'format':'json', 'removenullcolumns':True}
            headers, outString = mastQuery(request)
            outObject = json.loads(outString)
            starRas = np.array([x['ra'] for x in outObject['data']])
            starDecs = np.array([x['dec'] for x in outObject['data']])
    
    trySector = None
    if not (args.sector is None):
        trySector = args.sector
    fpgParmFileList = None
    if not (args.fpgParameterFiles is None):
        fpgParmFileList = [x for x in args.fpgParameterFiles]

        
    # Instantiate Spacecraft position info
    scinfo = TESS_Spacecraft_Pointing_Data(trySector=trySector, fpgParmFileList=fpgParmFileList)
    # Open output file if requested
#    if not (args.outputFile is None):
#        fout = open(args.outputFile, 'w')

    # Add header to outputfile
    if not (args.outputFile is None):
        fileOutputHeader(args.outputFile, fpgParmFileList=fpgParmFileList)
    else:
        # add single line header to stdout
        print('# TIC     |   RA      |   Dec     | EclipticLong | EclipticLat | Sector | Camera | Ccd | RowPix | ColPix')
    # Now make list of the star objects
    starList = make_target_objects(starTics, starRas, starDecs)
    #print('Finished converting coords to ecliptic')
    # Make rough determination as to which pointing camera combos are worth
    # Checking in detail and then do detailed checking
    findAny=False
    for i, curTarg in enumerate(starList):
        curTarg = doRoughPosition(curTarg, scinfo)
        #print('Rough Position Done: {:d} {:d}'.format(i, curTarg.ticid))
        # Look to see if target was in any sectors
        if len(curTarg.sectors)>0:
            uniqSectors = np.unique(curTarg.sectors)
            starRas = np.array([curTarg.ra])
            starDecs =  np.array([curTarg.dec])
            for curSec in uniqSectors:
                idxSec = np.where(scinfo.sectors == curSec)[0][0]
                starInCam, starCcdNum, starFitsXs, starFitsYs, starCcdXs, starCcdYs = scinfo.fpgObjs[idxSec].radec2pix(\
                           starRas, starDecs)
                for jj, cam in enumerate(starInCam):
                    xUse = starCcdXs[jj]
                    yUse = starCcdYs[jj]
                    maxCoord = 2049
                    if args.combinedFits:
                        xUse = starFitsXs[jj]
                        yUse = starFitsYs[jj]
                        maxCoord = 4097
                    if xUse>0 and yUse>0 and xUse<maxCoord and yUse<maxCoord:
                        findAny=True
                        strout = '{:09d} | {:10.6f} | {:10.6f} | {:10.6f} | {:10.6f} | {:2d} | {:1d} | {:1d} | {:8.3f} | {:8.3f}'.format(\
                           curTarg.ticid, curTarg.ra, curTarg.dec, curTarg.eclipLong,\
                           curTarg.eclipLat, curSec, starInCam[jj], starCcdNum[jj], xUse, yUse)
                        if not (args.outputFile is None):
                            args.outputFile.write('{:s}\n'.format(strout))
                        else:
                            print(strout)

    if not findAny:
        print('No Target/s were found to be on detectors')
        