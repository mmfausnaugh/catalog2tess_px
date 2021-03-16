import scipy as sp
from subprocess import call, check_call, Popen, PIPE
import os

def aber_stars_and_starspx(ra,dec,mag,                                                  
                           SC_coords,                                      
                           SC_vel,                                                  
                           cam):

    SC_ra, SC_dec, SC_roll = SC_coords
    v1,v2,v3 =  SC_vel
    
    path_use = os.path.dirname(os.path.realpath(__file__))
    aber_stars = os.path.join(path_use,'aber_stars1')
    starspx    =  os.path.join(path_use,'starspx8')
    
    geom_path = os.path.join(path_use,'../coords/')
    
    data_string = '\n'.join([' '.join(r.astype(str)) for r in sp.c_[ra,dec,mag] ])
    p1 = Popen([aber_stars, "1",str(v1), str(v2), str(v3)],
               stdin = PIPE, stdout = PIPE, bufsize=1,universal_newlines=True)
    aberrated_data_string = p1.communicate(input=data_string )[0]
    starpx_input_string = ''
    data_string_list = data_string.split('\n')
    ab_data_string_list = aberrated_data_string.split('\n')    
    for z in zip(data_string_list, ab_data_string_list):
        if len(z[0]) > 0 and len(z[1]) > 0:
            z0 = z[0].split()
            z1 = z[1].split()
            starpx_input_string += '\t'.join([z0[0], z0[1],z1[0], z1[1], z1[2] ]) + '\n'

    check_call(["ln", "-s",
                os.path.join(geom_path,"geometry_model_cam{:d}.txt".format(cam+1)),
                "fpg_pars.txt"])
    p1 = Popen([starspx, "1", str(SC_ra), str(SC_dec), str(SC_roll), "ZERO_BASE"],
               stdin=PIPE, stdout=PIPE, stderr=PIPE, bufsize=1, universal_newlines=True)
    starpx_output_string = p1.communicate(input=starpx_input_string.strip() )[0]
    check_call(["rm", "fpg_pars.txt"])
#    if len(starpx_output_string) == 0:
        #this condition handles, e.g., one source (that we overwrite
        #coordinates in the GI office) when we tryi it on the wrong CCD.
#        raise RuntimeError("Empty return of starspx8")
    idxrowcol =  sp.array([ s1.split() for s1 in starpx_output_string.split('\n') if len(s1) > 0])
    #return col,row index (to mask), and CCD
    if len(idxrowcol) == 0:
        idxrowcol = sp.empty((0,9))
    return (idxrowcol[:,7].astype(float),
            idxrowcol[:,8].astype(float),
            idxrowcol[:,0].astype(int),
            idxrowcol[:,6].astype(int))



class SectorCoords(object):
    """This close knows about TESS pointings and velocities, and can
    convert sky coordinates to pixel coordinates

    """

    #ra,dec,roll
    sector_coords = {'s1':[352.684371, -64.853067,-137.846820],
                     's2':[16.557116,-54.015977,-139.566543],
                     's3':[36.3138,-44.2590,-146.9616],
                     's4':[55.0070,-36.6420,-157.1698],
                     's5':[73.538201, -31.934855, -168.948284],
                     's6':[92.009629, -30.583860, 178.636723],
                     's7':[110.255949, -32.634366, 166.447620],
                     's8':[128.115605, -37.736986, 155.309071],
                     's9':[145.907085, -45.304357, 145.916291],
                     's10':[165.047459, -54.816523, 139.172410],
                     's11':[189.124736, -65.536906, 138.076125],
                     's12':[229.588462, -75.125547, 153.977330],
                     's13':[298.675214, -76.327466, -161.057110],
                     's14':[276.716943, 62.475619, 32.232893],
                     's15':[280.398495, 64.067133, 55.427671],
                     's16':[282.4425, 66.1416, 79.4629],
                     's17':[351.238098, 57.845603, 41.968602],
                     's18':[16.110329, 67.957534, 40.545339],
                     's19':[60.202647, 76.234293, 19.646262],
                     's20':[129.386732, 75.252007, -25.431070],
                     's21':[171.795057, 65.192354, -42.050483],
                     's22':[197.100781, 53.743420, -40.300770],
                     's23':[217.287871, 43.807416, -32.575431],
                     's24':[261.451571, 63.118073, -42.737575],
                     's25':[265.609754, 61.938287, -20.470660],
                     's26':[270.138139, 61.563699, 0.603849],
                     's27':[326.852535, -72.426456, -145.493925],
                     's28':[357.294358, -63.005592, -137.478445],
                     's29':[18.918982, -52.829616, -140.203034],
                     's30':[38.356373, -43.317792, -147.955943], 
                     's31':[57.635686, -35.783459, -158.766643],
                     's32':[77.189072, -31.395666, -171.373732],
                     's33':[96.599650, -30.784763, 175.536874],
                     's34':[115.295067, -33.779012, 163.191581],
                     's35':[133.2035, -39.6871, 152.4006],
                     's36':[150.9497, -47.7512, 143.7306],
                     's37':[170.2540, -57.3725, 138.1685],
                     's38':[195.7176, -67.8307, 139.3519],
                     's39':[242.1981, -76.3969, 161.5986],                     
                 }


    #x,y,z velocity
    sector_velocities = {'s1':[22.15, 16.87, 7.03],
                         's2':[10.84,24.75,10.46],
                         's3':[-3.21,26.99,11.42],
                         's4':[-16.72,23.06,9.69],
                         's5':[-26.61, 13.87, 5.71],
                         's6':[-30.74, 1.60, 0.42],
                         's7':[-28.41, -10.93, -4.98],
                         's8':[-20.33, -21.05, -9.34],
                         's9':[-8.39, -26.74, -11.79],
                         's10':[5.14, -26.96, -11.89],
                         's11':[7.68, -21.43, -9.50],
                         's12':[26.28, -11.04, -4.98],
                         's13':[28.66, 1.83, 0.62],
                         's14':[24.34, 14.09, 5.97],
                         's15':[14.59, 23.02, 9.86],
                         's16':[1.81, 26.98, 11.57],
                         's17':[-11.32, 25.47, 10.91],
                         's18':[-22.38, 18.91, 8.06],
                         's19':[-29.32, 8.38, 3.52],
                         's20':[-30.45, -4.45, -2.00],
                         's21':[-24.64, -16.92, -7.36],
                         's22':[-12.81, -25.39, -11.01],
                         's23':[1.70, -27.45, -11.91],
                         's24':[15.07, -23.14, -10.06],
                         's25':[24.52, -14.17, -6.17],
                         's26':[28.6230, -2.7972, -1.2106],
                         's27':[27.0905, 8.8213, 3.8577],
                         's28':[20.3233, 18.8078, 8.2145],
                         's29':[9.1700, 25.4191, 11.0993],
                         's30':[-4.6989, 26.9173, 11.7474],
                         's31':[-18.2747, 22.1664, 9.6756],
                         's32':[-27.8651, 11.8030, 5.1899],
                         's33':[-30.8131, -1.4522, -0.5262],
                         's34':[-26.7453, -14.0494, -5.9518],
                         's35':[-17.2433, -23.1860, -9.8915],
                         's36':[-4.6973, -27.3443, -11.6954],
#                         's37':[]
    }


    def __init__(self,sector_num):
        try:
            self.sector_num = sector_num
            self.coords   = self.sector_coords[sector_num]
            self.velocity = self.sector_velocities[sector_num]
        except:
            print(sector_num)
            print('To initialize, give this class a string \'s1\','
                  ' \'s2\'... \'s12\' \'s13\' etc.')

    def radec2pix(self, camnum, ra, dec, mag, ignore_image_buffer = False):
        col1, row1, ccd1, ra1, dec1, mag1,idx = self.get_cam_pix(camnum,
                                                                 ra,
                                                                 dec,
                                                                 mag)
        col2,row2 = self.convert_cam2ccd_pix(col1, row1, ccd1)
        #print(self.sector_num, camnum, col1,row1, col2, row2)
        #on imaging region, and maske sure bkg annulus will fit
        #settled on bkground region inner/outer 4/8 pixels
        #seems to need some extra buffer for things near the bottom

        #for checking why I don't find something that I expect---how to make this robust?
        #print( self.sector_num, camnum,ccd1, col1,row1, col2,row2)

        if ignore_image_buffer:
            print('ignoring pixel buffer at edge of images...')
            m = (col2 < 2092 ) & (col2 >44 ) & ( row2 < 2048 ) & (row2 > 0 )
        else:
            m = (col2 < 2092 - 8.0) & (col2 > 44 + 8.0) & ( row2 < 2048 - 8.0) & (row2 > 0 + 8.0 + 2)

        if len(idx[m]) == 0:
            raise ValueError('No objects on this CCD in this Sector')
        
        return col1[m], row1[m], ccd1[m], col2[m], row2[m], ra1[m], dec1[m], mag1[m], idx[m]
        

    def get_cam_pix(self, camnum, ra, dec, mag):
        """Calls Al's code and converts ra/dec to camera coordinates.  Note
        that mag is just for book keeping, can be an array of ones if
        necessary.

        """
        
        
        col, row, idx, ccd  = aber_stars_and_starspx(ra,dec,mag,
                                                     self.coords,
                                                     self.velocity,
                                                     camnum)

        if len(idx) == 0:
            return sp.empty(0),sp.empty(0),sp.empty(0), sp.empty(0), sp.empty(0), sp.empty(0), sp.empty(0)
        else:
            ra,dec,mag  = ra[idx],dec[idx],mag[idx]
            return col, row, ccd, ra, dec, mag, idx

    def convert_cam2ccd_pix(self, col, row, ccd):
        colout, rowout = sp.zeros(len(col)), sp.zeros(len(row))
        
        #
        m = ccd == 1
        #To make CCD 0.0  --> cam, 4271
        #we need 2*2136 - 1.0, etc.
        colout[m] = 2*2136 - col[m] - 1.0
        rowout[m] = 2*2078 - row[m] - 1.0


        m = ccd == 2
        #To make CCD 0.0  --> cam, 4271
        #we need 2*2136 - 1.0, etc.
        colout[m] = 2136 - col[m] - 1.0
        rowout[m] = 2*2078 - row[m] - 1.0

        m = ccd == 3
        colout[m] = col[m]
        rowout[m] = row[m]

        m = ccd == 4
        colout[m] = col[m] - 2136
        rowout[m] = row[m]

        return colout, rowout

        
