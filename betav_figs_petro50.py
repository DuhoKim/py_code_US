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

petroR50_V = np.zeros(len(sha_cat))
petroR50_r = np.zeros(len(sha_cat))

petroM50_V = np.zeros(len(sha_cat))
petroM50_r = np.zeros(len(sha_cat))


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
#for x in range(48,49):
    name = sha_cat['col1'][x]      # NGC/IC name

    ######### NED match ############
    ned_match = ned_tot[ned_tot['col2'] == name]
    semi = ned_match['col11'][0]  # read semi-major axis from NED [']
    thumb = int(semi * 100)  # a half size of figure size

    ########### WCS read  #############
    w_SDSS = WCS(work_dir + 'SDSS/r/' + name + '-r.fits')
    w_SHA = WCS(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')

    i_SDSS, j_SDSS = w_SDSS.all_world2pix(sha_cat['col2'][x],sha_cat['col3'][x],1)
    i, j = w_SHA.all_world2pix(sha_cat['col2'][x],sha_cat['col3'][x],1)

    i_SDSS = int(i_SDSS)
    j_SDSS = int(j_SDSS)
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

    ############# READ SE SEGMENTATION MAP #############
    hdu = fits.open(work_dir + 'SHA/segmap/' + name + '.fits')
    i1_seg = hdu[0].data

    ########### SEGMENATION & select SNR > 3 pixels #################
    betav_seg = V / (i1 * 2.330451129)
    betav_seg[np.where(i1_seg != 1)] = float('nan')

    V_seg = copy.deepcopy(V)
    r_seg = copy.deepcopy(rir)
    V_seg[np.where(i1_seg != 1)] = float('nan')
    r_seg[np.where(i1_seg != 1)] = float('nan')

    xlb2 = i - thumb if i - thumb > 0 else 0
    xhb2 = i + thumb if i + thumb < V_seg.shape[1] else V_seg.shape[1]
    ylb2 = j - thumb if j - thumb > 0 else 0
    yhb2 = j + thumb if j + thumb < V_seg.shape[0] else V_seg.shape[0]

    fits_seg_V = V_seg[ylb2: yhb2, xlb2: xhb2]
    fits_seg_r = r_seg[ylb2: yhb2, xlb2: xhb2]

    if np.isnan(fits_seg_V[thumb, thumb]):
        fits_seg_V = V[ylb2: yhb2, xlb2: xhb2]
        fits_seg_r = rir[ylb2: yhb2, xlb2: xhb2]

    #########       V        #############
    rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote_circ(fits_seg_V, i, j, thumb, V_seg.shape[1], V_seg.shape[0])
    rad_indices = np.where(np.logical_and.reduce((rads > 5, rads < np.max(rads) * 0.8)))
    for ii in range(0, len(rad_indices[0])):
        rad_index = rad_indices[0][ii]
        anul_flux = 0
        for jj in range(int(rads[rad_index]*0.8), int(rads[rad_index]*1.25)):
            anul_idx = min(enumerate(np.abs(rads-jj)), key=itemgetter(1))[0]           # find closest radius
            anul_flux = anul_flux + (2 * np.pi * jj)  * iso_med[anul_idx]

        if anul_flux < tot_med[rad_index] * 0.2:
            half_light = tot_med[rad_index] * 0.5
            half_idx = min(enumerate(np.abs(tot_med - half_light)), key=itemgetter(1))[0]  # find closest radius
            petroR50_V[x] = rads[half_idx]
            double_idx = min(enumerate(np.abs(rads - rads[rad_index] * 2.0)), key=itemgetter(1))[0]  # find closest radius
            petroM50_V[x] = tot_med[double_idx] if double_idx < len(tot_med) else np.nanmax(tot_med) # Petrosian mag = mag within 2 * rp
            break

    if petroR50_V[x] == 0:
        half_light = np.nanmax(tot_med) / 2
        half_idx = min(enumerate(np.abs(tot_med - half_light)), key=itemgetter(1))[0]  # find closest radius
        petroR50_V[x] = rads[half_idx]
        petroM50_V[x] = np.nanmax(tot_med)  # Petrosian mag = mag within 2 * rp


    #########       r        #############
    rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote_circ(fits_seg_r, i, j, thumb, V_seg.shape[1], V_seg.shape[0])
    rad_indices = np.where(np.logical_and.reduce((rads > 5, rads < np.max(rads) * 0.8)))

    for ii in range(0, len(rad_indices[0])):
        rad_index = rad_indices[0][ii]
        anul_flux = 0
        for jj in range(int(rads[rad_index] * 0.8), int(rads[rad_index] * 1.25)):
            anul_idx = min(enumerate(np.abs(rads - jj)), key=itemgetter(1))[0]  # find closest radius
            anul_flux = anul_flux + (2 * np.pi * jj) * iso_med[anul_idx]

        if anul_flux < tot_med[rad_index] * 0.2:
            half_light = tot_med[rad_index] * 0.5
            half_idx = min(enumerate(np.abs(tot_med - half_light)), key=itemgetter(1))[0]  # find closest radius
            petroR50_r[x] = rads[half_idx]
            double_idx = min(enumerate(np.abs(rads - rads[rad_index] * 2.0)), key=itemgetter(1))[0]  # find closest radius
            petroM50_r[x] = tot_med[double_idx] if double_idx < len(tot_med) else np.nanmax(
                tot_med)  # Petrosian mag = mag within 2 * rp
            break

    if petroR50_r[x] == 0:
        half_light = np.nanmax(tot_med) / 2
        half_idx = min(enumerate(np.abs(tot_med - half_light)), key=itemgetter(1))[0]  # find closest radius
        petroR50_r[x] = rads[half_idx]
        petroM50_r[x] = np.nanmax(tot_med)  # Petrosian mag = mag within 2 * rp


    #####################################################################

np.save('/Users/dhk/work/py/pyraf/petroR50_r.npy', petroR50_r)
np.save('/Users/dhk/work/py/pyraf/petroR50_V.npy', petroR50_V)
np.save('/Users/dhk/work/py/pyraf/petroM50_r.npy', petroM50_r)
np.save('/Users/dhk/work/py/pyraf/petroM50_V.npy', petroM50_V)