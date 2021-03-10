import scipy as sp
import matplotlib.pyplot as plt
from astropy.io import fits
from copy import deepcopy
import os

import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from camera_pointings import cam_pointings


d = sp.genfromtxt('Gaia_DR2_variability_cat.txt',delimiter=',',names=True,
                  dtype=(sp.uint64, 
                         sp.float64, sp.float64, sp.float64, sp.float64,
                         sp.float64, sp.float64, sp.float64, sp.float64,
                         sp.float64, sp.float64, sp.float64, sp.float64,
                         sp.float64,str,
                         sp.float64, sp.float64, sp.float64, sp.float64,
                         sp.float64, sp.float64, sp.float64))


ra = d['ra']
dec = d['dec']
#put in cartesian coords, get unit vector.  radius = 1.0
x = sp.cos(dec*sp.pi/180.0)*sp.cos(ra*sp.pi/180.0)
y = sp.cos(dec*sp.pi/180.0)*sp.sin(ra*sp.pi/180.0)
z = sp.sin(dec*sp.pi/180.0)
norm = sp.sqrt(x**2 + y**2 + z**2)
x /= norm
y /= norm
z /= norm

cams  =[cam_pointings.cam1, 
        cam_pointings.cam2, 
        cam_pointings.cam3, 
        cam_pointings.cam4]

for jj,cam in enumerate(cams):
    for ii in range(len(cam)):
        ra1,dec1 = cam[ii]
        print(ra1,dec1, ii, jj)
        #put in cartesian coords, get unit vector.  radius = 1.0
        x1 = sp.cos(dec1*sp.pi/180.0)*sp.cos(ra1*sp.pi/180.0)
        y1 = sp.cos(dec1*sp.pi/180.0)*sp.sin(ra1*sp.pi/180.0)
        z1 = sp.sin(dec1*sp.pi/180.0)
        norm = sp.sqrt(x1**2 + y1**2 + z1**2)
        x1 /= norm
        y1 /= norm
        z1 /= norm



        dot_product = x1*x + y1*y + z1*z
        theta1 = sp.arccos(dot_product)
        m1 = (theta1 >= 0.0) & (theta1 < 18*sp.pi/180.0)

        dout = d[m1]

        if os.path.isdir('s{:02d}'.format(ii+1)) == False:
            os.makedirs('s{:02d}'.format(ii+1))
        sp.savetxt('s{:02d}/gaia_dr2_s{:02d}_cam{}.txt'.format(ii+1, ii+1, jj + 1),dout,
                   fmt='%d %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %f %f %f %f')


