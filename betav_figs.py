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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# from astropy.coordinates import SkyCoord
import math
from shutil import copyfile
# from pandas import DataFrame, read_csv
import pandas as pd
from operator import itemgetter
import copy

AV_fit = True
bv_fit = True
fig_flag = True
sdss_flag = False
save_fits_flag = True

cmap = 'gist_rainbow'  	# ## color map for betaV images
cmap2 = colors.ListedColormap(['white', 'cyan', 'gray'])
pscale = 0.6			# ## pixel scale ["]
snr = 3					# ## signal to noise ratio that we constrain

# bv0s = [0.645,0.64,0.695,0.74,0.895,1.075,1.41]
# ### beta_V,0 for E0,S0(SFH3),Sa,Sb,Sbc(SFH4),Sc,Sd(SFH5) from Table 5 Kim+17
# BetaVzero values for SFH3,4,5 and Z=.0001,.0004,.004,.008,.02,.05
#bv0s = [0.5, 0.6, 0.73, 0.73, 0.81, 0.85, 0.97]
bv0s = [0.59, 0.73, 0.81, 0.78, 0.92, 1.07, 1.22]    # new betaVzero values from individual A_V fitting
# BetaVzero values for SFH3,4,5 and Z=.0001,.0004,.004,.008,.02,.05 w/ no Z evolution
#bv0s = [[1.9408058, 1.8461447, 1.1066952, 0.80323478, 0.59881037, 0.38305162],
#        [1.9878865, 1.8940383, 1.1417368, 0.84291972, 0.63275795, 0.40320143],
#        [2.3858643, 2.2975972, 1.4978560, 1.1777806, 0.93497056, 0.59137094]]


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

petroR50_V = np.load('/Users/dhk/work/py/pyraf/petroR50_V.npy')
bv0_tot = np.load('/Users/dhk/Documents/publish/ngcic_rev/bvzero_indi.npy')

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


betav_mean = np.zeros(len(sha_cat))
betav_upper = np.zeros(len(sha_cat))
betav_lower = np.zeros(len(sha_cat))

betav_257 = np.zeros(len(sha_cat))
betav_257_std = np.zeros(len(sha_cat))
ttype_257 = np.zeros(len(sha_cat))
ttype_257_err = np.zeros(len(sha_cat))
betav_257_center = np.zeros(len(sha_cat))
betav_257_center_std = np.zeros(len(sha_cat))
hlr_pix = np.zeros(len(sha_cat))
hlr_petro = np.zeros(len(sha_cat))
mag_petro = np.zeros(len(sha_cat))

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

########################################## CIRCULAR SURFACE BRIGHTNESS MEASURE ##################################################
def isophote_circ(fits_f, cx, cy, rad_max, xd, yd):
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

    rads_f = np.arange(1, rad_max+1, dtype=float)

    for r in rads_f:
        # General Equation of an Ellipse : (x*cos(A) + y*sin(A))^2/a^2 + (x*sin(A)-y*cos(A))^2/b^2 = 1
        if r == 1:
            ind = r ** 2 >= (yy - cy) ** 2 + (xx - cx) ** 2
        else:
            ind1 = (r - 1) ** 2 < (yy - cy) ** 2 + (xx - cx) ** 2
            ind2 = r ** 2 >= (yy - cy) ** 2 + (xx - cx) ** 2
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

        circ_ellip = 2 * np.pi * r

        tot_mean_f[int(r-1)] = iso_mean_f[0] if r == 1 else tot_mean_f[int(r-2)] + iso_mean_f[int(r-1)] * circ_ellip
        tot_med_f[int(r-1)] = iso_med_f[0] if r == 1 else tot_med_f[int(r-2)] + iso_med_f[int(r-1)] * circ_ellip

    idx_valid = ~np.isnan(iso_med_f)

    return  rads_f[idx_valid],      iso_mean_f[idx_valid], iso_med_f[idx_valid], \
            iso_std_f[idx_valid],   tot_mean_f[idx_valid], tot_med_f[idx_valid]
########################################################################################################################

for x in range(0, len(sha_cat)):
#for x in range(0, 1):
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

    #############  RC3 catalogue ###########
    rc3_match = rc3[[j for j, s in enumerate(rc3['name']) if s.strip() == rc3name]]
    if len(rc3_match) != 1:
        print('rc3 match is not one')
        ttype_257_err[x] = 1
    elif rc3_match['eT'][0] == '*':
        ttype_257_err[x] = 1
    else:
        ttype_257_err[x] = rc3_match['eT'][0]

    ################ T-type match ###########
    ttype_match = ttype.loc[ttype['name'] == table_name]
    T = ttype_match.iloc[0, 5]
    ttype_257[x] = T

    if np.isnan(bv0_tot[x]):
        bv0 = -1
    else:
        bv0 = bv0_tot[x]

    if T < -3:
        AV_prof = rosa_AV[0][2:]
    if T < 0 and T >= -3:
        AV_prof = rosa_AV[2][2:]
    if T < 2 and T >= 0:
        AV_prof = rosa_AV[4][2:]
    if T < 4 and T >= 2:
        AV_prof = rosa_AV[6][2:]
    if T < 5 and T >= 4:
        AV_prof = rosa_AV[8][2:]
    if T < 7 and T >=5:
        AV_prof = rosa_AV[10][2:]
    if T < 9 and T >= 7:
        AV_prof = rosa_AV[12][2:]
    #if T >= 9:
        # no BetaVzero is available for these types

    ######### NED match ############
    ned_match = ned_tot[ned_tot['col2'] == name]
    if x == 90:
        GalExtFactor = 10.0 ** (0.072 / 2.5)				# NGC 4561, A_V=0.072	(NED)
    else:
        GalExtFactor = 10.0 ** (ned_match['col6'][0] / 2.5)		# read Galactic extinction from NED [Vint/Vobs]
    gal_type.append(ned_match['col13'][0])

    semi = ned_match['col11'][0]		# read semi-major axis from NED [']
    if sdss_flag:
        thumb = int(semi*100*0.6/0.396)			# a half size of figure size
    else:
        thumb = int(semi * 100)                 # a half size of figure size
    z = ned_match['col8'][0]		# read redshift from NED
    boa = ned_match['col12'][0] / ned_match['col11'][0]

    ########### WCS read  #############
    if sdss_flag:
        w = WCS(work_dir + 'SDSS/r/' + name + '-r.fits')
    else:
        w = WCS(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')

    i, j = w.all_world2pix(sha_cat['col2'][x],sha_cat['col3'][x],1)
    i = int(i)
    j = int(j)

    ############ READ FITS ##############
    hdu	= fits.open(work_dir+'SDSS/g/'+name+'-gir.fits')
    gir = hdu[0].data
    hdu = fits.open(work_dir+'SDSS/gir_seg/'+name+'-gir.seg.fits')
    gir_seg = hdu[0].data
    hdu = fits.open(work_dir+'SDSS/r/'+name+'-rir.fits')
    rir = hdu[0].data
    hdu = fits.open(work_dir + 'SDSS/r/' + name + '-r.fits')
    r_sdss = hdu[0].data
    hdu = fits.open(work_dir+'SDSS/rir_seg/'+name+'-rir.seg.fits')
    rir_seg = hdu[0].data
    hdu	= fits.open(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')
    i1 = hdu[0].data
    hdu = fits.open(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.std.fits')
    i1u = hdu[0].data

    ############# Calculate V-band ###########
    V = 10 ** (np.log10(gir) + 0.59 * (np.log10(rir / gir)) + 0.004) 	# Jester+05 All stars R-I < 1.15

    ########### Calculate betaV, betag ########
    betav_raw = V / (i1 * 2.330451129) * GalExtFactor		# correct for zero point difference between SDSS & Spitzer
    betag_raw = gir / (i1 * 2.330451129) * GalExtFactor	# and Galactic Extinction

    ############# READ SE SEGMENTATION MAP #############
    hdu = fits.open(work_dir + 'SHA/segmap/' + name + '.fits')
    i1_seg = hdu[0].data

    ########### SEGMENATION & select SNR > 3 pixels #################
    betav_seg = V / (i1 * 2.330451129) * GalExtFactor
    betav_seg[np.where(i1_seg != 1)] = float('nan')
    V_seg = copy.deepcopy(V)
    r_seg = copy.deepcopy(rir)
    g_seg = copy.deepcopy(gir)
    V_seg[np.where(i1_seg != 1)] = float('nan')
    r_seg[np.where(i1_seg != 1)] = float('nan')
    g_seg[np.where(i1_seg != 1)] = float('nan')

    if fig_flag:
        # only use central region due to stripe regions after registering
        gir_cent = gir[j-500:j+500, i-500:i+500]
        gir_seg_cent = gir_seg[j-500:j+500, i-500:i+500]
        rir_cent = rir[j-500:j+500, i-500:i+500]
        rir_seg_cent = rir_seg[j-500:j+500, i-500:i+500]
        gerr = np.nanstd(gir_cent[np.where(gir_seg_cent == 0)])
        rerr = np.nanstd(rir_cent[np.where(rir_seg_cent == 0)])

        # segmap has holes for masking foreground stars
        #ierr    = np.nanstd(i1[np.where(i1_seg_nomask == 0)])

        betav_seg[(i1 / i1u < snr) | (gir / gerr < snr) | (rir / rerr < snr)] = float('nan')

        betav_mean[x] = np.nanpercentile(betav_seg, 50)         # central value
        betav_upper[x] = np.nanpercentile(betav_seg, 84.13)     # upper value of 1 sigma range
        betav_lower[x] = np.nanpercentile(betav_seg, 15.87)     # lower value of 1 sigma range

        betav_257[x] = np.nanmean(betav_seg)
        betav_257_std[x] = np.nanstd(betav_seg)

        betav_257_center[x] = np.nanmean(betav_seg[int(j) - 1 : int(j) + 2, int(i) - 1 : int(i) + 2])
        betav_257_center_std[x] = np.nanstd(betav_seg[int(j) - 1 : int(j) + 2, int(i) - 1 : int(i) + 2])

        ################# PLOT FIGURES for each sample ###############

        fig = plt.figure(figsize=(12, 7.7))

        h11 = fig.add_axes([0.0, 0.5, 0.315, 0.5])
        h12 = fig.add_axes([0.315, 0.5, 0.315, 0.5])
        h13 = fig.add_axes([0.63, 0.5, 0.315, 0.5])

        h21 = fig.add_axes([0.0, 0.0, 0.315, 0.5])
        h22 = fig.add_axes([0.315, 0.0, 0.315, 0.5])

        minsize = min([min(i1.shape), min(rir.shape), min(gir.shape)])
        if thumb * 2 > minsize:
            thumb = int(minsize / 2)

        img11 = make_lupton_rgb(i1[j - thumb: j + thumb, i - thumb: i + thumb] * 2.33,
                                rir[j - thumb: j + thumb, i - thumb: i + thumb],
                                gir[j - thumb: j + thumb, i - thumb: i + thumb], Q=10, stretch=0.5)

        h11.imshow(img11, origin = 'lower')
        h11.get_xaxis().set_visible(False)
        h11.get_yaxis().set_visible(False)
        h11.text(0.02, 0.95, r'R : IRAC 3.6$\mu$m', color='white', fontsize=15, transform=h11.transAxes, fontweight='bold')
        h11.text(0.02, 0.90, 'G : SDSS r', color='white', fontsize=15, transform=h11.transAxes, fontweight='bold')
        h11.text(0.02, 0.85, 'B : SDSS g', color='white', fontsize=15, transform=h11.transAxes, fontweight='bold')
        h11.text(0.6, 0.95, table_name, color='white', fontsize=15, transform=h11.transAxes, fontweight='bold')
        h11.text(0.3, 0.05, ned_match['col13'][0], color='white', fontsize=15, transform=h11.transAxes,
             fontweight='bold')

        kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)         # read xxx kpc / arcmin
        arcmin_5kpc = 5.0 / kpc_arcmin			            # calculate xxx arcmin / 5 kpc
        frac_5kpc = arcmin_5kpc * 100.0 / (2 * thumb)	    # calculate fraction of a length of 5 kpc in fig size
        h11.plot([0.05, 0.05 + frac_5kpc.value], [0.02, 0.02], color='white', transform=h11.transAxes)
        h11.text(0.02, 0.05, '5kpc, '+'{:4.2f}'.format(arcmin_5kpc.value) + '\'', color='white', fontsize=12, transform=h11.transAxes)

        img1 = h12.imshow(betav_raw[j - thumb: j + thumb, i - thumb: i + thumb], cmap=cmap, clim=[0, 2], origin='lower')
        h12.get_xaxis().set_visible(False)
        h12.get_yaxis().set_visible(False)
        h12.text(0.02, 0.9, r'$\beta_{V}$ (V/3.6$\mu$m)', color='black', fontsize=18, transform=h12.transAxes, fontweight='bold')

        img2 = h13.imshow(betag_raw[j - thumb: j + thumb, i - thumb: i + thumb], cmap=cmap, clim=[0, 2], origin='lower')
        axins = inset_axes(h13, width="5%", height="97%", loc='lower left', bbox_to_anchor=(1.0, 0.0, 1.0, 1.0), bbox_transform=h13.transAxes)
        img2.axes.figure.colorbar(img2, cax=axins)
        h13.get_xaxis().set_visible(False)
        h13.get_yaxis().set_visible(False)
        h13.text(0.02, 0.9, r'$\beta_{g}$ (g/3.6$\mu$m)', color='black', fontsize=18, transform=h13.transAxes, fontweight='bold')

        i1_seg[np.where(i1_seg > 2)] = 2
        h21.imshow(i1_seg[j - thumb: j + thumb, i - thumb: i + thumb], cmap=cmap2, origin='lower')
        h21.get_xaxis().set_visible(False)
        h21.get_yaxis().set_visible(False)
        h21.text(0.02, 0.95, 'Segmentation', fontsize=15, transform = h21.transAxes, fontweight='bold')
        h21.text(0.02, 0.89, '& Mask out', fontsize=15, transform = h21.transAxes, fontweight='bold')

        h22.imshow(betav_seg[j-thumb:j+thumb, i - thumb: i + thumb], cmap=cmap, clim=[0, 2], origin='lower')
        h22.get_xaxis().set_visible(False)
        h22.get_yaxis().set_visible(False)
        h22.text(0.02, 0.95, 'S/N > 3', fontsize=15, transform=h22.transAxes, fontweight='bold')

        dust_map = 2.5*np.log10(bv0/betav_seg)
        dust_map[np.where(betav_seg > bv0)] = 0.0			# assign A_V = 0 for bv > bv0

        if save_fits_flag:
            betav_fits = betav_seg[j - thumb: j + thumb, i - thumb: i + thumb]
            av_fits = dust_map[j - thumb: j + thumb, i - thumb: i + thumb]
            hdu_betav = fits.PrimaryHDU(betav_fits)
            hdu_av = fits.PrimaryHDU(av_fits)
            hdu_betav.writeto('/Users/dhk/work/data/NGC_IC/betav_mean_seg/' + name + '_betaV.fits')
            hdu_av.writeto('/Users/dhk/work/data/NGC_IC/Av_seg/' + name + '_AV.fits')

########## surface profile measure #############
    if sdss_flag:
        xlb2 = i - thumb if i - thumb > 0 else 0
        xhb2 = i + thumb if i + thumb < r_sdss.shape[1] else r_sdss.shape[1]
        ylb2 = j - thumb if j - thumb > 0 else 0
        yhb2 = j + thumb if j + thumb < r_sdss.shape[0] else r_sdss.shape[0]
    else:
        xlb2 = i - thumb if i - thumb > 0 else 0
        xhb2 = i + thumb if i + thumb < V_seg.shape[1] else V_seg.shape[1]
        ylb2 = j - thumb if j - thumb > 0 else 0
        yhb2 = j + thumb if j + thumb < V_seg.shape[0] else V_seg.shape[0]

    if sdss_flag:
        fits_seg = r_sdss[ylb2: yhb2, xlb2: xhb2]
        rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(fits_seg, i, j, pa, boa, thumb, r_sdss.shape[1], r_sdss.shape[0])
    else:
        fits_seg = V_seg[ylb2: yhb2, xlb2: xhb2]
        if np.isnan(fits_seg[thumb, thumb]):
            # In case of masked center, use galfit model image
            hdu = fits.open('/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask/' + name + '.out.fits')
            galfit_model = hdu[2].data
            fitsect = hdu[2].header['FITSECT']  # read FITS sector information	[x1:x2,y1:y2]
            fitsect_div = fitsect.split(',')
            fitsect_x = fitsect_div[0].replace('[', '')
            fitsect_y = fitsect_div[1]
            fitsect_x_div = fitsect_x.split(':')
            fitsect_y_div = fitsect_y.split(':')
            x_off = int(fitsect_x_div[0])
            y_off = int(fitsect_y_div[0])
            rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote_circ(galfit_model[j - y_off - thumb:j - y_off + thumb,
                i - x_off - thumb:i - x_off + thumb], i - x_off, j - y_off, thumb, galfit_model.shape[1],galfit_model.shape[0])
        else:
            rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote_circ(fits_seg, i, j, thumb, V_seg.shape[1], V_seg.shape[0])


    if fig_flag:
        # h25 = fig.add_axes([0.843-0.3,0.09,0.1,0.1])
        # h25.plot(rads[~np.isnan(tot_med)],tot_med[~np.isnan(tot_med)])
        # h25.plot([rads[hlr_idx],rads[hlr_idx]],[0,np.nanmax(tot_med)])

        if bv0 != -1:
            h23 = fig.add_axes([0.63, 0.0, 0.315, 0.5])
            img3=h23.imshow(dust_map[j-thumb:j+thumb, i-thumb:i+thumb], cmap='jet', clim=[0,1], origin='lower')
            axins = inset_axes(h23, width="5%", height="97%", loc='lower left', bbox_to_anchor=(1.0, 0.0, 1.0, 1.0),
                               bbox_transform=h23.transAxes)
            img3.axes.figure.colorbar(img3, cax=axins)
            h23.text(0.02, 0.95, 'DUST MAP (Av)', fontsize=15, transform=h23.transAxes, fontweight='bold')
            #h23.text(0.02, 0.89, r'2.5$\times$log('+'{:5.3f}'.format(bv0)+r'/$\beta_{V}$)', fontsize=15, transform=h23.transAxes, fontweight='bold')
            h23.text(0.02, 0.89, r'$\beta_{V,0}$=' + '{:5.3f}'.format(bv0), fontsize=15, transform=h23.transAxes, fontweight='bold')
            h23.get_xaxis().set_visible(False)
            h23.get_yaxis().set_visible(False)

            ellipse1 = Ellipse((thumb,thumb), petroR50_V[x], petroR50_V[x] * boa, pa+90, fill=False, color='yellow')
            ellipse2 = Ellipse((thumb,thumb), petroR50_V[x] * 3, petroR50_V[x] * 3 * boa, pa+90, fill=False, color='yellow')
            h23.add_artist(ellipse1)
            h23.add_artist(ellipse2)

            ######## surface betaV & A_V profile measure ##########
            if bv_fit:
                rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(betav_seg[j-thumb:j+thumb,
                                                i-thumb:i+thumb], i, j, pa, boa, thumb, i1.shape[1], i1.shape[0])
                np.savetxt(work_dir+'ellipse/scratch_ellip/'+name+'_bv.txt', (rads,iso_mean,iso_med,iso_std))

            if AV_fit:
                rads, iso_mean_map, iso_med_map, iso_std_map, tot_mean, tot_med = isophote(dust_map[j-thumb:j+thumb,
                                                i-thumb:i+thumb], i, j, pa, boa, thumb, i1.shape[1], i1.shape[0])

                h24 = fig.add_axes([0.843, 0.39, 0.1, 0.1])
                h24.plot(rads[~np.isnan(iso_med_map)] / petroR50_V[x], iso_med_map[~np.isnan(iso_med_map)], color='red')
                h24.fill_between(rads[~np.isnan(iso_med_map)] / petroR50_V[x], iso_med_map[~np.isnan(iso_med_map)] -
                                 iso_std_map[~np.isnan(iso_med_map)], iso_med_map[~np.isnan(iso_med_map)] +
                                 iso_std_map[~np.isnan(iso_med_map)], color='red', alpha=0.3)
                h24.plot(np.arange(0., 3., 0.1), AV_prof[1:], color='black', alpha=0.5)
                h24.errorbar(1, AV_prof[11], yerr=AV_prof[0], color='black', alpha=0.5, capsize=2)
                plt.xticks([1, 3])
                h24.set(xlabel=r'$R/R_{50,r}^{P}$', ylabel='Av')
                h24.get_xaxis().set_label_coords(0.6,-0.2)
                h24.get_yaxis().set_label_coords(-0.1,0.5)

        else:
            h23 = fig.add_axes([0.63, 0.0, 0.315, 0.5])
            h23.get_yaxis().set_visible(False)
            h23.get_xaxis().set_visible(False)
            h23.text(0.02, 0.9, 'No '+r'$\beta_{V,0}$ is available for this T-type', size=15)

    if fig_flag:
        fig.savefig(work_dir + 'pdf/' + name + '.tex.pdf')
        plt.close(fig)

        # if int((sha2fig[x] - 1) % 2):
        #     copyfile(work_dir + 'pdf/' + name + '.tex.pdf', '/Users/dhk/Documents/publish/ngcic/figset2_petro/figset2_' +
        #              str(int((sha2fig[x] - 1) / 2)) + 'b.pdf')
        # else:
        #     copyfile(work_dir + 'pdf/' + name + '.tex.pdf', '/Users/dhk/Documents/publish/ngcic/figset2_petro/figset2_' +
        #              str(int((sha2fig[x] - 1) / 2)) + 'a.pdf')

        copyfile(work_dir + 'pdf/' + name + '.tex.pdf', '/Users/dhk/Documents/publish/ngcic_rev2/figset2_2nd/figset2_' +
                     str(int((sha2fig[x]))) + '.pdf')



#np.save('/Users/dhk/work/py/pyraf/betav_ttype.npy',[betav_257,betav_257_std,ttype_257,ttype_257_err,betav_257_center,betav_257_center_std])
#np.save('/Users/dhk/work/py/pyraf/betav_range.npy',[betav_mean,betav_upper,betav_lower])
#np.save('/Users/dhk/work/py/pyraf/hlr_V_circ.npy', hlr_pix)
#np.save('/Users/dhk/work/py/pyraf/hlr_V_petro_circ.npy', hlr_petro)
#np.save('/Users/dhk/work/py/pyraf/mag_V_petro_circ.npy', mag_petro)
#np.save('/Users/dhk/work/py/pyraf/av.npy', [av_257_Z0_05, av_257_Z0_05_std, av_257_Z0_02, av_257_Z0_02_std])