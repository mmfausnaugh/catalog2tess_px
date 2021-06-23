#!/usr/bin/env python

import numpy as np
import os
import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..') 
import argparse
from catalogs.GaiaCsv import GaiaCsv


def get_inputs(args):
    parser = argparse.ArgumentParser(
        description="Specify Sector, Camera/CCD, and period, and make a phot.data file for DIA photometry.")
    parser.add_argument('--cam', type=int, nargs="*",  help="Camera Number (1--4)")
    parser.add_argument('--ccd', type=int, nargs="*",   help="CCD Number (1--4)")
    parser.add_argument('--sector', type=int, nargs="*", help="Sector Number (for start/stop time")
    parser.add_argument('--mag0', type=float, default=-20.0, help="Magnitude limit: collect all sources fainter than this value (default = -20)")
    parser.add_argument('--mag1', type=float, default=20.0, help="Magnitude limit: collect all sources brighter than this value (default = 20)")
    parser.add_argument('--outdir',default='./', help="Output directory stem: phot.data appears in subdirs sectorP/camN_ccdM")
    parser.add_argument('--doall', action='store_true', help="If set, do all cams and CCDs in the sectors")

    return parser.parse_args()


def main():

    args = get_inputs(sys.argv[1:])
    if args.doall:
        cam_use = [1,2,3,4]
        ccd_use = [1,2,3,4]
    else:
        cam_use = args.cam
        ccd_use = args.ccd
    for sector in args.sector:        
        for cam in cam_use:
            for ccd in ccd_use:
                #find the catalog file

                path1 = os.path.abspath(__file__)
                path2 = os.path.split(path1)
                path3 = os.path.split(path2[0])
                pathfinal = path3[0] + "/Gaia/"
                
                catfile = 's{s_use:02d}/gaia_dr2_s{s_use:02d}_cam{cam:d}.txt'.format(
                    s_use = sector,
                    cam   = cam )
                catfile = os.path.join(pathfinal, catfile)                

                #instantiate the catalog object
                cat = GaiaCsv(catfile)
                
                #select the masks
                mag_mask  = (cat.phot_rp_mean_mag < args.mag1) & ( cat.phot_rp_mean_mag > args.mag0)
                ccd_mask  = cat.ccd == ccd
                row_mask = (cat.ccdrow > 10 ) & (cat.ccdrow < 2038)
                col_mask = (cat.ccdcol > 54 ) & (cat.ccdcol < 2082)
                mask = mag_mask & ccd_mask & row_mask & col_mask

                #make the file
                outdir = os.path.join(args.outdir,'sector{:02d}'.format(sector),'cam'+str(cam)+'_ccd'+str(ccd))
                if not os.path.isdir(outdir):
                    os.makedirs(outdir)

                name_use = np.array(["{:019d}".format(np.uint64(i)) for i in cat.source_id])
                cat.make_phot_data(cat.ccdcol[mask],
                                   cat.ccdrow[mask],
                                   name_use[mask],
                                   outdir )      


if __name__ == '__main__':
    main()

