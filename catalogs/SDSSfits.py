import scipy as sp
from astropy.time import Time
from astropy.io import fits
import sys
import os
sys.path.insert(0, os.path.abspath(   os.path.dirname(__file__)) + '/..')
from coords.coords import SectorCoords
from .catalog import Catalog


class SDSSfits(Catalog):
    #update by MMF 2020-08 for DR16

    keys = [ #col 1, 2, 3
        'sdss_name', 'ra', 'dec',
        #col 66--74
        'boss_target1','eboss_target0','eboss_target1',
        'eboss_target2','ancillary_target1','ancillary_target2',
        'nspec_sdss','nspec_boss','nspec   ',
        #col 95/96
        'psfmag','psfmag_err',
        #col 97/98
        'extinction','abs_i_mag',
        #col 100--103
        'galex_fuv','galex_fuv_ivar','galex_nuv','galex_nuv_ivar',
        #col 105--112; units of flux, not super helpful
        #'ukids_y','ukids_y_err',
        #'ukids_j','ukids_j_err',
        #'ukids_h','ukids_h_err',
        #'ukids_k','ukids_k_err',
        #col 115/116
        'w1_mag','w1_mag_err',
        #col 125/126
        'w2_mag','w2_mag_err',
        #col 134/135/136
        'first_flux','first_snr','first_separation'
        '2mass_j','2mass_j_err',
        
            # id, RA/DEC, unique bit ID, MJD, spectrograph -- columns 1-4, 6, 8
            'Z', 'Z_ERR', 'SOURCE_Z', 'Z_VI', 'Z_PIPE', 'Z_PIPE_ERR',
            'ZWARNING', 'Z_PCA', 'Z_PCA_ER', 'Z_MGII',
            # All the redshift estimates! -- columns 9-18
            'N_SPEC_SDSS', 'N_SPEC_BOSS', 'N_SPEC', 
            # Number of spectra, total, SDSS, BOSS -- columns 25-27
            'BI_CIV', 'ERR_BI_CIV', 
            # Balnicity -- columns 32-33
            'PSFFLUX_u', 'PSFFLUX_g', 'PSFFLUX_r', 'PSFFLUX_i', 'PSFFLUX_z',
            'IVAR_IVAR_u', 'IVAR_IVAR_g', 'IVAR_IVAR_r', 'IVAR_IVAR_i', 'IVAR_IVAR_z',
            'PSFMAG_u', 'PSFMAG_g', 'PSFMAG_r', 'PSFMAG_i', 'PSFMAG_z',
            'ERR_PSFMAG_u', 'ERR_PSFMAG_g', 'ERR_PSFMAG_r', 'ERR_PSFMAG_i', 'ERRPSFMAG_z',
            'MI', 'GAL_EXT_u', 'GAL_EXT_g', 'GAL_EXT_r', 'GAL_EXT_i', 'GAL_EXT_z',
            # SDSS mags, abs_mag (K-corrected), E(B-V) -- columns 39-44
            
            # Matching Catalog Photometry: skipping ROSAT, units are counts
            'FLUX_0.2_2.0keV', 'FLUX_0.2_2.0keV_ERR', 'FLUX_2.0_12.0keV', 
            'FLUX_2.0_12.0keV_ERR', 'FLUX_0.2_12.0keV', 'FLUX_0.2_12.0keV_ERR',
            'LUM_0.2_12.0keV', 'SDSS2XMM_SEP',
            # XMM [X-ray source] -- columnns 48-55
            
            'GALEX_MATCHED', 'FUV', 'FUV_IVAR', 'NUV', 'NUV_IVAR',
            # GALEX [UV survey] -- columns 56-60
            
            'JMAG', 'ERR_JMAG', 'JSNR', 'JRDFLAG', 'HMAG', 'ERR_HMAG', 'HSNR',
            'HRDFLAG', 'KMAG', 'ERR_KMAG', 'KSNR', 'KRDFLAG', 'SDSS2MASS_SEP',
            # 2MASS [IR survey] -- columns 61-73
            
            'W1MAG', 'ERR_W1MAG', 'W1SNR', 'W1CHI2', 'W2MAG', 'ERR_W2MAG', 
            'W2SNR', 'W2CHI2', 'W3MAG', 'ERR_W3MAG', 'W3SNR', 'W3CHI2',
            'W4MAG', 'ERR_W4MAG', 'W4SNR', 'W4CHI2', 'CC_FLAGS', 'PH_FLAG',
            'SDSS2WISE_SEP',
            # WISE [mid IR survey] -- columns 74-92
            
            'UKIDSS_MATCHED', 'YFLUX', 'YFLUX_ERR', 'JFLUX', 'JFLUX_ERR', 
            'HFLUX', 'HFLUX_ERR', 'KFLUX', 'KFLUX_ERR', 
            # UKIDS [IR survey] -- columns 93-101
            
            'FIRST_MATCHED', 'FIRST_FLUX', 'FIRST_SNR', 'SDSS2FIRST_SEP',
            # FIRST [radio] -- columns 102-105
            
            'camcol','camrow','ccd','ccdcol','ccdrow']
      

    def __init__(self, ifile, ignore_image_buffer=False):
        
        f = fits.open(ifile)
        data = f[1].data
        
        for i in data['Z_PCA']:
            if i == -1:
                i == sp.isnan
        for i in data['Z_PCA_ER']:
            if i == -1:
                i == sp.isnan
        for i in data['N_SPEC_SDSS']:
            if i == -1:
                i == sp.isnan
        for i in data['N_SPEC_BOSS']:
            if i == -1:
                i == sp.isnan
        for i in data['N_SPEC']:
            if i == -1:
                i == sp.isnan
        for i in data['BI_CIV']:
            if i == -1:
                i == sp.isnan
        for i in data['ERR_BI_CIV']:
            if i == -1:
                i == sp.isnan

        
        ra = data['RA']
        dec = data['DEC']
        mag = data['PSFMAG'][:,3]

        fstem = os.path.basename(ifile)
        sector = 's' + str(int(fstem.split('_')[2][1:])) #str and int to strip leading zeros
        cam = int(fstem.split('_')[3][-6]) - 1

        out = self.get_pix(sector,cam,ra,dec, mag, 
                           ignore_image_buffer=ignore_image_buffer)
        camcol = out[0]
        camrow = out[1]
        ccd    = out[2]
        ccdcol = out[3]
        ccdrow = out[4]
        idx = out[8]

        # print sector, cam, idx
        # print d
        data = data[idx]
        # print d, ccd
        ra = ra[idx]
        dec = dec[idx]
        
        
        super(SDSSfits,self).__init__(self.keys,
                                         [data['SDSS_NAME'], data['RA'], data['DEC'], data['THING_ID'], 
                                          data['MJD'], data['SPECTRO'], data['Z'], data['Z_ERR'], data['SOURCE_Z'],
                                          data['Z_VI'], data['Z_PIPE'], data['Z_PIPE_ERR'], data['ZWARNING'], 
                                          data['Z_PCA'], data['Z_PCA_ER'], data['Z_MGII'], data['N_SPEC_SDSS'],
                                          data['N_SPEC_BOSS'], data['N_SPEC'], data['BI_CIV'], data['ERR_BI_CIV'], 
                                          
                                          data['PSFFLUX'][:,0], data['PSFFLUX'][:,1], data['PSFFLUX'][:,2], data['PSFFLUX'][:,3], data['PSFFLUX'][:,4],
                                          data['IVAR_PSFFLUX'][:,0], data['IVAR_PSFFLUX'][:,1], data['IVAR_PSFFLUX'][:,2], data['IVAR_PSFFLUX'][:,3], data['IVAR_PSFFLUX'][:,4],
                                          data['PSFMAG'][:,0], data['PSFMAG'][:,1], data['PSFMAG'][:,2], data['PSFMAG'][:,3], data['PSFMAG'][:,4],
                                          data['ERR_PSFMAG'][:,0], data['ERR_PSFMAG'][:,1], data['ERR_PSFMAG'][:,2], data['ERR_PSFMAG'][:,3], data['ERR_PSFMAG'][:,4],
                                          
                                          data['MI'],
                                          
                                          data['GAL_EXT'][:,0], data['GAL_EXT'][:,1], data['GAL_EXT'][:,2], data['GAL_EXT'][:,3], data['GAL_EXT'][:,4], 
                                          data['FLUX_0.2_2.0keV'], data['FLUX_0.2_2.0keV_ERR'],
                                          data['FLUX_2.0_12.0keV'], data['FLUX_2.0_12.0keV_ERR'], data['FLUX_0.2_12.0keV'],
                                          data['FLUX_0.2_12.0keV_ERR'], data['LUM_0.2_12.0keV'], data['SDSS2XMM_SEP'],
                                          data['GALEX_MATCHED'], data['FUV'], data['FUV_IVAR'], data['NUV'], 
                                          data['NUV_IVAR'], data['JMAG'], data['ERR_JMAG'], data['JSNR'], 
                                          data['JRDFLAG'], data['HMAG'], data['ERR_HMAG'], data['HSNR'],
                                          data['HRDFLAG'], data['KMAG'], data['ERR_KMAG'], data['KSNR'], data['KRDFLAG'],
                                          data['SDSS2MASS_SEP'], data['W1MAG'], data['ERR_W1MAG'], data['W1SNR'],
                                          data['W1CHI2'], data['W2MAG'], data['ERR_W2MAG'], data['W2SNR'],
                                          data['W2CHI2'], data['W3MAG'], data['ERR_W3MAG'], data['W3SNR'], data['W3CHI2'],
                                          data['W4MAG'], data['ERR_W4MAG'], data['W4SNR'], data['W4CHI2'], data['CC_FLAGS'],
                                          data['PH_FLAGS'], data['SDSS2WISE_SEP'], data['UKIDSS_MATCHED'], data['YFLUX'],
                                          data['YFLUX_ERR'], data['JFLUX'], data['JFLUX_ERR'], data['HFLUX'],
                                          data['HFLUX_ERR'], data['KFLUX'], data['KFLUX_ERR'], data['FIRST_MATCHED'],
                                          data['FIRST_FLUX'], data['FIRST_SNR'], data['SDSS2FIRST_SEP'],
                                          
                                          camcol, camrow, ccd, ccdcol, ccdrow])
        #all catalogs must have obj_name attribute
        self.obj_name = self.SDSS_NAME
