#####################################################
# Python script that reads betaV images which are
# FITS images that has pixel values of the flux ratios
# between V = g-0.59(g-r)-0.01 (Jester+05) and Spitzer 
# IRAC 1 band and IRAF PHOT result and mask out bg sources
# written by Duho Kim (2/2/18)        
######################################################
from pyraf import iraf
from astropy.io import fits
from astropy.io import ascii
import numpy as np
from astropy.wcs import WCS
import os
import sys

snr = 1

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

#sha_cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')
sha_cat=ascii.read(cat_dir+'sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match.csv')

#for x in range(0,len(eha_cat)):
for x in range(1,21):
        name    = sha_cat['col1'][x]      # NGC/IC name

	##### WCS read  #####
	w	= WCS(work_dir+'SHA/SHA_NGC_IC_LONG/'+name+'.fits')
	i,j	= w.all_world2pix(sha_cat['col2'][x],sha_cat['col3'][x],1)

	##### READ FITS #####
	hdu	= fits.open(work_dir+'SDSS/g/'+name+'-gir.fits')
	gir      = hdu[0].data
        hdu     = fits.open(work_dir+'SDSS/gir_seg/'+name+'-gir.seg.fits')
        gir_seg = hdu[0].data
        hdu     = fits.open(work_dir+'SDSS/r/'+name+'-rir.fits')
        rir     = hdu[0].data
        hdu     = fits.open(work_dir+'SDSS/rir_seg/'+name+'-rir.seg.fits')
        rir_seg = hdu[0].data
	hdu	= fits.open(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')
	i1	= hdu[0].data
	hdu	= fits.open(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.unc.fits')
	i1u	= hdu[0].data
	hdu	= fits.open(work_dir+'betav_mean/'+name+'.betaV.fits')
	betav	= hdu[0].data	
	hdu	= fits.open(work_dir+'SHA/check_images/'+name+'.seg.back_off.thresh1.fits')
	check	= hdu[0].data

	try:
		##### LEAVE ONLY NGC/IC galaxy VALUES #####
		betav[np.where(check!=check[int(j),int(i)])] = float('nan')
	except:
		print name
		continue
	
	###### Masking any pixel doesn't meet our S/N limit ##########
	gerr	= np.nanstd(gir[np.where(gir_seg==0)])
	rerr    = np.nanstd(rir[np.where(rir_seg==0)])
	ierr	= np.nanstd(i1[np.where(check==0)])
	
	betav[(i1/i1u < snr) | (gir/gerr < snr) | (rir/rerr < snr)] = float('nan')	# mask out also the region with S/N < 1

	###### Saving betav image #############################
	fits.writeto(work_dir+'betav_mean/'+name+'.betav.snr1.back_off.thresh1.fits',betav,overwrite=True)

	betav[(i1/i1u < snr+2) | (gir/gerr < snr+2) | (rir/rerr < snr+2)] = float('nan')
	fits.writeto(work_dir+'betav_mean/'+name+'.betav.snr3.back_off.thresh1.fits',betav,overwrite=True)





