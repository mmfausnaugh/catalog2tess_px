import scipy as sp
from astropy.time import Time
import sys
import os
import pandas as pd
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords

#in TJD
sector_times = {'s1':[1324.70843962, 1353.49455073],
                's2':[1353.49455073, 1381.66746739],
                's3':[1381.66746739, 1410.83413406],
                's4':[1410.83413406, 1437.91746739],
                's5':[1437.91746739, 1464.40288407],
                's6':[1464.40288407, 1491.62533921],
                's7':[1491.62533921, 1516.1],
                's8':[1516.1, 1542.04495359],
                's9':[1542.04495359, 1568.52827666],
                's10':[1568.52827666, 1595.78798727],
                's11':[1595.78798727, 1624.5323],
                's12':[1624.5323, 1652.91],
                's13':[1652.91, 1682.50],
                's14':[1682.50, 	1710.4],
                's15':[1710.4, 1737.90000],
                's16':[1737.9 ,  1764.0],
                's17':[1764.0,1789.71],
                's18':[1789.71, 1815.4174674],
                's19':[1815.4174674, 1841.30],
                's20':[1841.3, 1868.425595974713],
                's21':[1868.425595974713, 1898.29246742],
                's22':[1898.29246742, 1926.80],
                's23':[1926.80, 1955.75991],
                's24':[1955.75991, 1982.40],
                's25':[1982.40, 2010.20],
                's26':[2010.20,2036.5],
                's27':[2036.5, 2060.9],
                's28':[2060.9, 2087.35],
                's29':[2087.35, 2114.70],
                's30':[2114.70, 2143.4],
                's31':[2143.4, 2170.1],
                's32':[2170.1, 2201.93],
                's33':[2201.93,2227.77],
                's34':[2227.77,2254.26476],
                's35':[2254.26476, 2280.1],
                's36':[2280.1,  2306.1],
                's37':[2306.1,  2332.75790],
                's38':[2332.75790, 2360.741],
                's39':[2360.741, 2389.728],
                's40':[2389.728, 2419.05359],
                's41':[2419.05359, 2446.76956],
                's42':[2446.76956,2473.55],
                's43':[2473.55, 2499.0],
                's44':[2499.0, 2524.63965],
                's45':[2524.63965, 2550.92436],
                's46':[2550.92436, 2578.90074],
                's47':[2578.90074, 2607.13960],
                's48':[2607.13960, 2636.1 ],
                's49':[2636.1, 2664.55 ],
                's50':[2664.55, 2691.55 ],
                's51':[2691.55, 2717.8],
                's52':[2717.8, 2743.3],
                's53':[2743.3,2769.2],
                's54':[2769.2,2796.4],
                's55':[2796.4,2825.4],
                's56':[2825.4, 2853.3],
                's57':[2853.3, 2882.28],
                's58':[2882.28, 2910.25],
                's59':[2910.25, 2936.9],
                's60':[2936.9, 2962.80],
                's61':[2962.80, 2988.3],
                's62':[2988.3, 3014.29],
                's63':[3014.29,3040.92],
                's64':[3040.92, 3068.7],
                's65':[3068.7, 3096.60],
                's66':[3096.60, 3126.60],
                's67':[3126.60, 3154.40],
                's68':[3154.40,3182.15],
                's69':[3182.15,3208.33],
                's70':[3208.33,3233.99],
                's71':[3233.99,3260.06],
                's72':[3260.06, 3285.78],
                's73':[3285.78,3312.80],
                's74':[3312.80,3339.72],
                's75':[3339.72, 3367.65],
                's76':[3367.65, 3394.79],
                's77':[3394.79, 3424.0],
                's78':[3424.0, 3452.51],
                's79':[3452.51, 3479.68],
                's80':[3479.68,3506.35],
                's81':[3506.35, 3533.19],
                's82':[3533.19, 3558.89],
                's83':[3558.89,3584.59],
                's84':[3584.59,5000],

                }


class Catalog(object):
    """Base class for individual catalog, but knows how to do easy things
    like convert to/from sexigesimal and call 'coord' module to get
    pixel coordinates.


    """
    def __init__(self, keys, arrays):
        for z in zip(keys,arrays):
            setattr(self, z[0], z[1])

    
    @staticmethod
    def sexigesimal_to_decimal(ra,dec, return_string=True):
        RA = ((ra[2]/60. + ra[1])/60. + ra[0])*180/12.
        if dec[0] < 0:
            DEC = -(dec[2]/60. + dec[1])/60. + dec[0]
        else:
            DEC = (dec[2]/60. + dec[1])/60. + dec[0]
        return RA,DEC
        
    @staticmethod
    def decimal_to_sexigesimal(ra,dec, return_string=True):
        ra1 = (ra*12./180).astype(int)
        ra2 = ((ra*12./180 - ra1)*60).astype(int)
        ra3 = ((ra*12./180 - ra1)*60 - ra2)*60

        dec1 = dec.astype(int)
        dec2 = ((dec - dec1)*60).astype(int)
        dec3 = ((dec - dec1)*60 - dec2)*60

        #        RA  = sp.array(['%i:%02i:%02.4f'%(ra1[i],ra2[i],ra3[i]) for i in range(ra.size)])
        #        DEC = sp.array(['%i:%02i:%02.4f'%(dec1[i],dec2[i],dec3[i]) for i in range(dec.size)])
        #going to enforce simple data types for now
        RA = '%i:%02i:%02.4f'%(ra1,ra2,ra3)
        DEC = '%i:%02i:%02.4f'%(dec1,dec2,dec3)
        return RA,DEC

    def get_pix(self, sector, cam, ra, dec, mag,ignore_image_buffer=False):
        s_use = SectorCoords(sector)
        return s_use.radec2pix(cam,ra,dec,mag, ignore_image_buffer = ignore_image_buffer)

    def make_phot_data(self, ccdcol, ccdrow, obj_name, outdir):
        fmt = '{:.3f} {:.3f} {:d} {:d} lc/lc_{:s} 1\n'
        with open( os.path.join(outdir,'phot.data'),'w') as fout:
            for i in range(len(obj_name)):
                fout.write(fmt.format(
                    ccdcol[i], ccdrow[i], 
                    ccdcol[i].astype(int), ccdrow[i].astype(int), 
                    obj_name[i]
                ))

    def make_coa_input(self, ccdcol, ccdrow, obj_name, ra, dec, mag, outdir):
        #needs to add it to the catalogStruct, targetlist data and pixeltalbe data
        #catalogStruct:         catId: [380002x1 double]
	#	              raDegrees: [380002x1 double]
	#	             decDegrees: [380002x1 double]
	#	                   pmRa: [380002x1 double]
	#	                  pmDec: [380002x1 double]
	#	                tessMag: [380002x1 double]
	#	                 radius: [380002x1 double]
	#	          effectiveTemp: [380002x1 double]
	#	    log10SurfaceGravity: [380002x1 double]
	#	       log10Metallicity: [380002x1 double]
        #targetlist:  catId
        #             index (in matlab array)
        #             labels (from ctl)
        #             pixel table indices (col vector)
        #             sectorNumber
        #            3 dummy structs which will be copied from pipeline inputsStruct
        #pixeltable:  index, ccdRow, ccdCol
        #
        fmt = '{} {} {} {} {} {} {}\n'
        catID = 9999999999

        with open( os.path.join(outdir,'coa_input.txt'),'w') as fout:
#            fout.write('%catID label RA DEC Cols Rows\n')
            for ii in range(len(obj_name)):
                #obj_name will go into the label
                c1 = sp.r_[ccdcol[ii] - 5 : ccdcol[ii] + 6].astype(int)
                r1 = sp.r_[ccdrow[ii] - 5 : ccdrow[ii] + 6].astype(int)
                C,R = sp.meshgrid(c1,r1)
                ccdcol_array = C.reshape(C.size,1)
                ccdrow_array = R.reshape(R.size,1)
                ccdcol_array = ','.join([str(_[0]) for _ in ccdcol_array])
                ccdrow_array = ','.join([str(_[0]) for _ in ccdrow_array])

                fout.write(fmt.format(
                    catID,
                    obj_name[ii],
                    ra[ii],
                    dec[ii],
                    mag[ii],
                    ccdcol_array, ccdrow_array, 
                ))

                catID -= 1


    def filter_on_sector_time(self,sector, early_offset=0.0, late_offset=0.0):
         tstart = sector_times[sector][0] - early_offset
         tend   = sector_times[sector][1] + late_offset
         
         mask = (self.tjd >= tstart) & ( self.tjd <= tend)
         return mask

