#####################################################
# Python script that reads betaV images which are
# FITS images that has pixel values of the flux ratios
# between V = g-0.59(g-r)-0.01 (Jester+05) and Spitzer 
# IRAC 1 band and analyze
# written by Duho Kim (2/19/18)        
######################################################
# from pyraf import iraf
from astropy.io import fits
from astropy.io import ascii
import numpy as np
# import numpy.ma as ma
from astropy.wcs import WCS
# from astropy.table import Table
from astropy.table import vstack
from astropy.visualization import make_lupton_rgb
from astropy.cosmology import Planck15 as Cosmo
# import os
# import os.path
# import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import colors
# import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
# from astropy.coordinates import SkyCoord
import math
from shutil import copyfile
# from pandas import DataFrame, read_csv
import pandas as pd
from operator import itemgetter
import copy

pscale = 0.6			# ## pixel scale ["]
snr = 3					# ## signal to noise ratio that we constrain

# bv0s = [0.645,0.64,0.695,0.74,0.895,1.075,1.41]
# ### beta_V,0 for E0,S0(SFH3),Sa,Sb,Sbc(SFH4),Sc,Sd(SFH5) from Table 5 Kim+17
# BetaVzero values for SFH3,4,5 and Z=.0001,.0004,.004,.008,.02,.05
bv0s_grid = [0.5, 0.6, 0.7, 0.7, 0.75, 0.85, 1.0]      # central value for +-0.05 grid for bv0


# Hubble types corresponding Nair+10 T-Type
htypes = ['E0', 'S0', 'S0/a', 'Sa', 'Sab', 'Sb', 'Sbc', 'Sc', 'Scd', 'Sd', 'Sdm', 'Sm', 'Im', '?']
# Kauffmann+03 from Nair+10
agntypes = ['SF', 'transition/mixed AGN', 'Seyfert', 'LINER']

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

ned1 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")
# row | input object name | object name | RA | Dec | Gal. Ext. Burstein & Heiles A_B mag | Object Type 
# (1) |       (2)         |    (3)      | (4)| (5) |                (6)                  |    (7)
# Redshift | Redshift Uncertainty | Mag./Filter | Major Diam | Minor Diam | Morph. | Ref.| 
# (8)      |        (9)           |    (10)     |   (11)     |    (12)    |  (13)  | (14)|
ned_tot = vstack([ned1, ned2, ned3, ned4, ned5])				# Total NED output catalog


#  sha_cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')
sha_cat = ascii.read(cat_dir + 'sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv')
#  Sample NGC/IC numbers
#  load figure numbers for each galaxies
sha2fig = np.load('/Users/dhk/work/data/NGC_IC/pdf/sha2fig.npy')

rc3 = ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
#  read A_V values for different Hubble types from CALIFA paper
rosa_AV = np.genfromtxt("/Users/dhk/work/data/rosa/AV.txt")

#  read NGC/IC catalogue
steineke = r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'
df = pd.read_excel(steineke)
#  read T-type
ttype_xls = r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(ttype_xls)

####################### Jarrett+03 ########################################################
#  Modified logarithmic visualization method P' = sqrt( log(1+P/(n*sig))  ), where P is the pixel
# intensity value, sig is the image rms "noise",
#  and n is a threshold throttle (w/ satisfactory values between 5 and 10
def jarrett(pix_values, sigma):
    n = 7.5
    return np.sqrt(np.log10(1+pix_values/(n*sigma)))
############################################################################################

gal_type = []
betav_med = []

########################################## ELLIPTICAL SURFACE BRIGHTNESS MEASURE ##################################################
def isophote(fits_f, cx, cy, pa_f, boa_f, rad_max, xd, yd):
    xlb = cx - rad_max if cx - rad_max > 0 else 0       # bound within a FITS file
    xhb = cx + rad_max if cx + rad_max < xd else xd
    ylb = cy - rad_max if cy - rad_max > 0 else 0
    yhb = cy + rad_max if cy + rad_max < yd else yd
    yy, xx = np.ogrid[ylb:yhb, xlb:xhb]                 # opposite convention between FITS & Python array indexing

    iso_mean_f = np.zeros(rad_max, dtype=float)
    iso_med_f = np.zeros(rad_max, dtype=float)
    iso_std_f = np.zeros(rad_max, dtype=float)
    tot_mean_f = np.zeros(rad_max, dtype=float)
    tot_med_f = np.zeros(rad_max, dtype=float)
    iso_mean_f.fill(np.nan)
    iso_med_f.fill(np.nan)
    iso_std_f.fill(np.nan)
    tot_mean_f.fill(np.nan)
    tot_med_f.fill(np.nan)

    pa_rad = math.radians(-pa_f)             # convert PA convention to functional form with units
    rads_f = np.arange(1, rad_max+1, dtype=float)
    rads_eff_f = np.sqrt(rads_f * rads_f * boa_f)

    for r in rads_f:
        # General Equation of an Ellipse : (x*cos(A) + y*sin(A))^2/a^2 + (x*sin(A)-y*cos(A))^2/b^2 = 1
        if r == 1:
            ind = 1 >= ((yy - cy) * math.cos(pa_rad) + (xx - cx) * math.sin(pa_rad)) ** 2 / r ** 2 + \
                  ((yy - cy) * math.sin(pa_rad) - (xx - cx) * math.cos(pa_rad)) ** 2 / (r * boa_f) ** 2
        else:
            ind1 = 1 < ((yy - cy) * math.cos(pa_rad) + (xx - cx) * math.sin(pa_rad)) ** 2 / (r-1) ** 2 + \
                   ((yy - cy) * math.sin(pa_rad) - (xx - cx) * math.cos(pa_rad)) ** 2 / ((r-1) * boa_f) ** 2
            ind2 = 1 >= ((yy - cy) * math.cos(pa_rad) + (xx - cx) * math.sin(pa_rad)) ** 2 / r ** 2 + \
                   ((yy - cy) * math.sin(pa_rad) - (xx - cx) * math.cos(pa_rad)) ** 2 / (r * boa_f) ** 2
            ind = ind1 * ind2
        anul = fits_f[ind]

        if np.sum(np.isnan(anul)) < np.sum(~np.isnan(anul)):
            iso_mean_f[int(r-1)] = np.nanmean(anul)
            iso_med_f[int(r-1)] = np.nanmedian(anul)
            iso_std_f[int(r-1)] = np.nanstd(anul)
        # If half of the anulus is masked but still inside
        elif r < rad_max / 5:
            iso_mean_f[int(r - 1)] = np.nanmean(anul)
            iso_med_f[int(r - 1)] = np.nanmedian(anul)
            iso_std_f[int(r - 1)] = np.nanstd(anul)
        # If half of the anulus is masked and outside then exclude
        else:
            break

        circ_ellip = 2 * np.pi * np.sqrt((r ** 2 + (r * boa_f) ** 2) / 2)

        tot_mean_f[int(r-1)] = iso_mean_f[0] if r == 1 else tot_mean_f[int(r-2)] + iso_mean_f[int(r-1)] * circ_ellip
        tot_med_f[int(r-1)] = iso_med_f[0] if r == 1 else tot_med_f[int(r-2)] + iso_med_f[int(r-1)] * circ_ellip

    idx_valid = ~np.isnan(iso_med_f)

    return  rads_eff_f[idx_valid],      iso_mean_f[idx_valid], iso_med_f[idx_valid], \
            iso_std_f[idx_valid],   tot_mean_f[idx_valid], tot_med_f[idx_valid]
########################################################################################################################


for x in range(0, len(sha_cat)):
#for x in range(0,1):
    name = sha_cat['col1'][x]      # NGC/IC name

    ########### NGC match ###########
    if name[0] == 'n':
        galnum = name[3:].strip()
        ngc_match = df.loc[(df['N'] == 'N') & (df['NI'] == int(galnum))]
        if len(galnum) == 3:
            galnum = '0' + galnum
        elif len(galnum) == 2:
            galnum = '00' + galnum
        elif len(galnum) == 1:
            galnum = '000' + galnum
        rc3name = 'NGC' + galnum
        table_name = 'NGC ' + galnum
    elif name[0] == 'i':
        galnum = name[2:].strip()
        ngc_match = df.loc[(df['N'] == 'I') & (df['NI'] == int(galnum))]
        if len(galnum) == 3:
            galnum = '0' + galnum
        elif len(galnum) == 2:
            galnum = '00' + galnum
        elif len(galnum) == 1:
            galnum = '000' + galnum
        rc3name = 'IC' + galnum
        table_name = 'IC ' + galnum

    pa = ngc_match.iloc[0,21]
    if math.isnan(pa):
        pa = 0
    elif pa > 90:
        pa = pa - 180.0

    ################ T-type match ###########
    ttype_match = ttype.loc[ttype['name'] == table_name]
    T = ttype_match.iloc[0, 5]
    if T < -3:                  # E
        T_idx = 0
    if T < 0 and T >= -3:       # S0
        T_idx = 1
    if T < 2 and T >= 0:        # Sa
        T_idx = 2
    if T < 4 and T >= 2:        # Sb
        T_idx = 3
    if T < 5 and T >= 4:        # Sbc
        T_idx = 4
    if T < 7 and T >=5:         # Sc
        T_idx = 5
    if T < 9 and T >= 7:        # Sd
        T_idx = 6
    if T >= 9:
        # no BetaVzero is available for these types
        bv0 = -1

    ######### NED match ############
    ned_match = ned_tot[ned_tot['col2'] == name]
    if x == 90:
        GalExtFactor = 10.0 ** (0.072 / 2.5)				# NGC 4561, A_V=0.072	(NED)
    else:
        GalExtFactor = 10.0 ** (ned_match['col6'][0] / 2.5)		# read Galactic extinction from NED [Vint/Vobs]
    gal_type.append(ned_match['col13'][0])

    semi = ned_match['col11'][0]		# read semi-major axis from NED [']
    thumb = int(semi * 100)                 # a half size of figure size
    z = ned_match['col8'][0]		# read redshift from NED
    boa = ned_match['col12'][0] / ned_match['col11'][0]

    ########### WCS read  #############
    w = WCS(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')

    i, j = w.all_world2pix(sha_cat['col2'][x],sha_cat['col3'][x],1)
    i = int(i)
    j = int(j)

    ############ READ FITS ##############
    hdu	= fits.open(work_dir+'SDSS/g/'+name+'-gir.fits')
    gir = hdu[0].data
    hdu = fits.open(work_dir+'SDSS/r/'+name+'-rir.fits')
    rir = hdu[0].data
    hdu	= fits.open(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')
    i1 = hdu[0].data
    hdu = fits.open(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.std.fits')
    i1u = hdu[0].data

    hdu = fits.open(work_dir + 'SDSS/gir_seg/' + name + '-gir.seg.fits')
    gir_seg = hdu[0].data
    hdu = fits.open(work_dir + 'SDSS/rir_seg/' + name + '-rir.seg.fits')
    rir_seg = hdu[0].data

    ############# Calculate V-band ###########
    V = 10 ** (np.log10(gir) + 0.59 * (np.log10(rir / gir)) + 0.004) 	# Jester+05 All stars R-I < 1.15

    ############# READ SE SEGMENTATION MAP #############
    hdu = fits.open(work_dir + 'SHA/segmap/' + name + '.fits')
    i1_seg = hdu[0].data

    ########### SEGMENATION & select SNR > 3 pixels #################
    gir_cent = gir[j - 500:j + 500, i - 500:i + 500]
    gir_seg_cent = gir_seg[j - 500:j + 500, i - 500:i + 500]
    rir_cent = rir[j - 500:j + 500, i - 500:i + 500]
    rir_seg_cent = rir_seg[j - 500:j + 500, i - 500:i + 500]
    gerr = np.nanstd(gir_cent[np.where(gir_seg_cent == 0)])
    rerr = np.nanstd(rir_cent[np.where(rir_seg_cent == 0)])

    betav_seg = V / (i1 * 2.330451129) * GalExtFactor
    betav_seg[np.where(i1_seg != 1)] = float('nan')
    betav_seg[(i1 / i1u < snr) | (gir / gerr < snr) | (rir / rerr < snr)] = float('nan')

    #####################################################################
    if T < 9:
        for ii in range(-5, 5, 1):
            dust_map = 2.5 * np.log10((bv0s_grid[T_idx] + ii * 0.01) / betav_seg)

            dust_map[np.where(betav_seg > (bv0s_grid[T_idx] + ii * 0.01))] = 0.0

            rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(dust_map[j-thumb:j+thumb,
                                        i-thumb:i+thumb], i, j, pa, boa, thumb, i1.shape[1], i1.shape[0])
            np.savetxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_'+str(ii)+'.txt', (rads, iso_mean, iso_med, iso_std))

