#####################################################
# Python script of building PSF for given FITS images 
# of SDSS data
# written by Duho Kim (12/26/17)	
######################################################
import numpy as np
from pyraf import iraf
from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import os
import matplotlib.pyplot as plt
from scipy.stats import mode
import math
import shutil
import sys

data_dir = '/Users/dhk/work/data/NGC_IC/SDSS/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'
work_dir = data_dir+'g_bench4/'
psf_dir = data_dir+'g_psf4'

os.chdir(work_dir)
#os.system("mkiraf")

sha_cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')

#### BASIC PARAMETERS ##########
pixel_size = 0.396      # Pixel size of SDSS data
irac_fwhm = 1.2         # Initial guess of FWHM of SDSS data ["]
zmag = 22.5		# Zero point of magnitude scale

#num_psf_stars = 10 (SHA)
num_psf_stars = 100      # Number of PSF stars
psf_radius = 30         # Radius of PSF [pixel]
nan_check = 2		# NaN check around a PSF stars with a radius of nan_check * psf_radius
tol_star_match = 1.0    # Tolerance of coordinate matching of PSF stars between sex catalog and IRAF tables

##### SELECT PSF STARS #########################################
# Select objects from the catalog resulted by Source Extractor
# above which is having 'FLAGS' 0 and 'CLASS_STAR' larger than 
# half and ELLIPTICITY (1-b/a) is less than 0.1 and FWHM_IMAGE
# is in between +- 1.5 pixels centered in mode value of FWHM_IMAGE
# values and sort with magnitudes 'MAG_AUTO' to build PSF from
################################################################
#class_range = 0.7 (SHA)
class_range = 0.8       # Maximum value of CLASS_STAR selecting PSF stars form sex catalog
ellip_range = 0.1       # Maximum value of ELLIPTICITY selecting PSF stars from sex catalog 

mag_limit = 17.0

##### PARAMETER SETTING FOR PyRAF ##############################################
# referenced Manuel 'PSF photometry using DAOPHOT' written by H. Hwang & M. Lee
################################################################################
iraf.daophot()
iraf.daopars.function = 'auto'          # fit all 6 analytic functions and use a function having minimum norm scatter value
iraf.daopars.varorder = 0               # -1 results simple not-realistic model, 1 and 2 are noisy
iraf.daopars.nclean = 10                 # showed improvement with an increament until 10 but not much difference compared to 100
                                        # The number of additional iterations the PSF task performs to compute the PSF look-up tables.
                                        # If nclean is > 0, stars which contribute deviant residuals to the PSF look-up tables in the 
                                        # first iteration, will be down-weighted in succeeding iterations. (default=0)
iraf.daopars.matchrad = 3.0             # Object matching radius in scale units (default=3.0)
iraf.daopars.psfrad = psf_radius        # Opt to use 60x60 size PSF image size (DK)
iraf.daopars.fitrad = 10.0              # The fitting radius in scale units. Only pixels within the fitting radius of the center of a star will contribute
                                        # to the fits computed by the PEAK, NSTAR and ALLSTAR tasks. For most images the fitting radius should be approximately
                                        # equal to the FWHM of the PSF. Under severely crowded conditions a somewhat smaller value may be used in order to 
                                        # improve the fit. If the PSF is variable, the FWHM is very small, or sky fitting is enabled in PEAK and NSTAR on the
                                        # other hand, it may be necessary to increase the fitting radius to achieve a good fit.


iraf.daopars.recenter = 'yes'           # Recenter stars during fit?
iraf.daopars.fitsky = 'yes'             # Recompute group sky value during fit? 
iraf.daopars.groupsky = 'yes'           # Use group rather than individual sky values? (default=yes)
iraf.daopars.sannulus = 0.0             # The inner radius of the sky annulus used by ALLSTAR to recompute the sky values. (default=0.0)
iraf.daopars.wsannulus = 11.0           # Width of sky fitting annulus in scale units (default=11.0)
iraf.daopars.critsnratio = 1.0          # The ratio of the model intensity of the brighter star computed at a distance of one fitting radius from the center
                                        # of the fainter star,  (default=1.0)
#iraf.findpars.threshold = 4.0           # Threshold in sigma for feature detection (default=4.0)
iraf.findpars.threshold = 20.0           # Threshold in sigma for feature detection (default=4.0)
iraf.findpars.nsigma = 1.5              # Width of convolution kernel in sigma  (default=1.5)
iraf.findpars.sharplo = 0.2             # Lower bound on sharpness for feature detection (default=0.2)
iraf.findpars.sharphi = 1.0             # Upper bound on sharpness for feature detection (default=1.0)



for x in range(0,len(sha_cat)):
	fn = sha_cat['id'][x]+'-g.fits'		# Image file name
	shutil.copy(data_dir+'g/'+fn,work_dir)

	##### READ FITS #####
	hdu = fits.open(fn)
	fits_data = hdu[0].data

	##### RUN SEXTRACTOR & GET BACKGROUND, FWHM #####
	os.system("sex "+fn+" -PIXEL_SCALE "+str(pixel_size)+" -SEEING_FWHM "+str(irac_fwhm)+' -CHECKIMAGE_TYPE BACKGROUND \
		-CHECKIMAGE_NAME bkg_'+fn+' -VERBOSE_TYPE QUIET')

	bkg_hdu = fits.open('bkg_'+fn)		# Read background FITS image
	bkg = bkg_hdu[0].data
	bg = np.median(bkg)			# Adopt median of background image as background value (mode function takes too much time)
	sig = np.std(bkg)			# Adopt std of background image as sigma value

	cat = ascii.read('test.cat')		# Read output catalog of SExtractor
	fwhm = mode(cat['FWHM_IMAGE'])[0][0]	# Adopt mode of FWHM as input FWHM for PyRAF [pixel] (doesn't take too much time for this)
	fwhm_sig = np.std(cat['FWHM_IMAGE'])	# std of FWHM [pixel]

	stars = cat[np.where((cat['FLAGS'] == 0) & (cat['CLASS_STAR'] > class_range) & (cat['ELLIPTICITY'] < \
	    ellip_range) & (cat['FWHM_IMAGE'] < fwhm+fwhm_sig) & (cat['FWHM_IMAGE'] > fwhm-fwhm_sig))]

	stars.sort('MAG_AUTO')	# Sort output catalog in the order of magnitude(luminosity)

	##### PARAMETER SETTING FOR PyRAF ##############################################
	# referenced Manuel 'PSF photometry using DAOPHOT' written by H. Hwang & M. Lee
	################################################################################
	iraf.datapars.fwhmpsf = fwhm
	iraf.datapars.sigma = sig

	iraf.centerpars.cbox = max(5, 2*fwhm)
	iraf.fitskypars.skyvalu = bg
	iraf.fitskypars.annulus = 4*fwhm
	iraf.fitskypars.dannulus = 3*fwhm
	iraf.photpars.aperture = max(3, fwhm)
	iraf.photpars.zmag = zmag
	iraf.daopars.fitrad = fwhm              # Fitting radius in scale units. Only pixels within the fitting radius of the center of a star will contribute to
                                        # the fits computed by the PEAK, NSTAR and ALLSTAR tasks.       

	##### Get Coordiates and Magnitudes Stars #################
	# Run DAOFIND to get coordinates and PHOT to get magnitudes
	###########################################################
	iraf.daofind(fn,fn[:-5]+'.coo.1',verify='no',verbose='no')			# generate coo.1 file
	iraf.phot(fn,fn[:-5]+'.coo.1',fn[:-5]+'.mag.1',skyfile='bkg_'+fn,verify='no',verbose='no')	# generate mag.1 file

	##### MAKE PSF STAR LIST FILE 'pst.1' #######
	# read coordinates and magnitudes of the stars selected based 
	# on the catalog generated by Source Extractor and make 
	# 'pst.1' file to use for running PSF task
	#######################################################
	coo1 = ascii.read(fn[:-5]+'.coo.1')
	mag1 = ascii.read(fn[:-5]+'.mag.1')
	coo1.sort('ID')
	mag1.sort('ID')
	pst1 = Table(names=('ID','XCENTER','YCENTER','MAG','MSKY'),dtype=('i4','f8','f8','f8','f8'))
	for i in range(min(num_psf_stars,len(stars))):
		ind = np.where((abs(stars['X_IMAGE'][i]-coo1['XCENTER']) < tol_star_match) & \
			(abs(stars['Y_IMAGE'][i]-coo1['YCENTER']) < tol_star_match))[0]
		if len(ind)==0 or np.isnan(fits_data[int(stars['Y_IMAGE'][i] - psf_radius * nan_check) : int(stars['Y_IMAGE'][i]+psf_radius * nan_check ), \
        	        int(stars['X_IMAGE'][i] - psf_radius * nan_check) : int(stars['X_IMAGE'][i] + psf_radius * nan_check)]).any() or \
			(mag1['MAG'][ind] < mag_limit).any():
			continue
		pst1.add_row( ( mag1['ID'][ind[0]], mag1['XINIT'][ind[0]], mag1['YINIT'][ind[0]], mag1['MAG'][ind[0]], mag1['MSKY'][ind[0]]) )
	
	pst1.write(fn[:-5]+'.pst.0', format='ascii.commented_header', delimiter='\t', comment= '#N ID    XCENTER   YCENTER \
	    MAG         MSKY\n#U ##    pixels    pixels    magnitudes  counts\n#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-15.7g\n#')

	original = sys.stdout
	try:
		sys.stdout = open(fn[:-5]+'.psf1.out','w')
		iraf.psf(fn[:-5],photfile=fn[:-5]+'.mag.1',pstfile=fn[:-5]+'.pst.0',psfimage=fn[:-5]+'.psf.1',opstfile=fn[:-5]+'.pst.1', \
			groupfile=fn[:-5]+'.psg.1',verify='no',interactive='no',plotfile=fn[:-5]+'.psf1.plots') 
		iraf.allstar(fn[:-5],photfile=fn[:-5]+'.psg.1',psfimage=fn[:-5]+'.psf.1',allstarfile=fn[:-5]+'.als.1',rejfile=fn[:-5]+ \
		    '.arj.1',subimage=fn[:-5]+'.sub.1',verify='no',verbos='no')
		iraf.substar(fn[:-5],photfile=fn[:-5]+'.als.1',exfile=fn[:-5]+'.pst.1',psfimage=fn[:-5]+'.psf.1',subimage=fn[:-5]+ \
		    '.sub.11',verify='no',verbose='no')

        	sys.stdout = open(fn[:-5]+'.psf2.out','w')
		iraf.psf(fn[:-5]+'.sub.11',photfile=fn[:-5]+'.mag.1',pstfile=fn[:-5]+'.pst.1',psfimage=fn[:-5]+'.psf.2', \
		    opstfile=fn[:-5]+'.pst.2', groupfile=fn[:-5]+'.psg.2',verify='no',interactive='no')
		iraf.allstar(fn[:-5],photfile=fn[:-5]+'.psg.2',psfimage=fn[:-5]+'.psf.2',allstarfile=fn[:-5]+'.als.2',rejfile=fn[:-5]+ \
		    '.arj.2',subimage=fn[:-5]+'.sub.2',verify='no',verbos='no')
		iraf.substar(fn[:-5],photfile=fn[:-5]+'.als.2',exfile=fn[:-5]+'.pst.2',psfimage=fn[:-5]+'.psf.2',subimage=fn[:-5]+ \
		    '.sub.22',verify='no',verbose='no')

        	sys.stdout = open(fn[:-5]+'.psf3.out','w')
		iraf.psf(fn[:-5]+'.sub.22',photfile=fn[:-5]+'.mag.1',pstfile=fn[:-5]+'.pst.2',psfimage=fn[:-5]+'.psf.3',opstfile \
		    =fn[:-5]+'.pst.3', groupfile=fn[:-5]+'.psg.3',verify='no',interactive='no')
	except iraf.IrafError, e:
		sys.stdout = original
		print fn
		print e
		os.system('mv test.cat '+fn[:-5]+'.cat')
		continue

	sys.stdout = original
	iraf.seepsf(fn[:-5]+'.psf.1',fn[:-5]+'.psf1')
	iraf.seepsf(fn[:-5]+'.psf.2',fn[:-5]+'.psf2')
	iraf.seepsf(fn[:-5]+'.psf.3',fn[:-5]+'.psf3')

	shutil.copy(fn[:-5]+'.psf1.fits',psf_dir)
	shutil.copy(fn[:-5]+'.psf2.fits',psf_dir)
	shutil.copy(fn[:-5]+'.psf3.fits',psf_dir)

	os.system('mv test.cat '+fn[:-5]+'.cat')
	os.system('rm *.fits')
#	os.system('rm '+fn[:-5]+'*')










