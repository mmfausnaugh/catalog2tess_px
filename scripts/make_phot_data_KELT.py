#!/usr/bin/env python

import scipy as sp
import os
import sys
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..') 
import argparse
from catalogs.KELTAscii import KELTAscii


def get_inputs(args):
    parser = argparse.ArgumentParser(
        description="Specify Sector, Camera/CCD, and period, and make a phot.data file for DIA photometry.")
    parser.add_argument('--cam', type=int, nargs="*",  help="Camera Number (1--4)")
    parser.add_argument('--ccd', type=int, nargs="*",   help="CCD Number (1--4)")
    parser.add_argument('--sector', type=int, nargs="*", help="Sector Number (for start/stop time")
    parser.add_argument('--mag', type=float, default=20.0, help="Magnitude limit: collect all sources brighter than this value (default = 20)")
    parser.add_argument('--outdir',default='./', help="Output directory stem: phot.data appears in subdirs sectorP/camN_ccdM")
    parser.add_argument('--period',type=float,default=13.5, help="Cut objects with periods greater than this value, in days")
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
                try:
                    #find the catalog file
                    
                    path1 = os.path.abspath(__file__)
                    path2 = os.path.split(path1)
                    path3 = os.path.split(path2[0])
                    pathfinal = path3[0] + "/KELT/kelt_catalog.txt"
                    catfile = os.path.join(pathfinal)

                    #instantiate the catalog object
                    cat = KELTAscii(catfile, sector, cam)
                
                
                    #select the masks
                    mag_mask  = cat.Vmag < args.mag
                    ccd_mask  = cat.ccd == ccd
                    period_mask = cat.Period_days <= args.period

                    mask = mag_mask & ccd_mask & period_mask

                    #make the file
                    outdir = os.path.join(args.outdir,'sector'+str(sector),'cam'+str(cam)+'_ccd'+str(ccd))
                    if not os.path.isdir(outdir):
                        os.makedirs(outdir)
                    cat.make_phot_data(cat.ccdcol[mask],
                                       cat.ccdrow[mask],
                                       cat.KELT_ID[mask],
                                       outdir , 'lc_KELT')      
                except:
                    continue


if __name__ == '__main__':
    main()

