#####################################################
# Python script that reads betaV images which are
# FITS images that has pixel values of the flux ratios
# between V = g-0.59(g-r)-0.01 (Jester+05) and Spitzer 
# IRAC 1 band and analyze
# written by Duho Kim (2/19/18)        
######################################################
from pyraf import iraf
from astropy.io import fits
from astropy.io import ascii
import numpy as np
from astropy.wcs import WCS
from astropy.table import Table, vstack
from astropy.visualization import make_lupton_rgb
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

cmap = 'gist_rainbow'	### color map for betaV images

thresh = ['.thresh05','.thresh10','.thresh15','.thresh20']
thresh1 = ['.thresh25','.thresh30','.thresh35','.thresh40']
back = ['','.backoff','.local','.back128','.back32','.nthresh2','.nthresh1','.nthresh1.back256','.nthresh1.backoff']

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])


#sha_cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')
sha_cat=ascii.read(cat_dir+'sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match.csv')

#for x in range(0,len(sha_cat)):

######## Jarrett+03 #############
# Modified logarithmic visualization method
# P' = sqrt( log(1+P/(n*sig))  )
# where P is the pixel intensity value,
# sig is the image rms "noise",
# and n is a threshold throttle (w/ satisfactory
# values between 5 and 10
#################################
def Jarrett(pix_values,sigma):
	n = 7.5
	return np.sqrt(np.log10(1+pix_values/(n*sigma)))

def colorbar(mappable):
	ax = mappable.axes
	figg = ax.figure
	divider = make_axes_locatable(ax)
	cax = divider.append_axes("right",size="5%",pad=0.05)
	return figg.colorbar(mappable, cax=cax)

def draw(name,fig,i,j,row,col,thumb):
	if (name == 'ngc4736' or name == 'ngc2403') and row == 1:
		hdu = fits.open(work_dir+'SHA/check_images/'+name+thresh1[col]+back[8]+'.fits')
	else:
		hdu = fits.open(work_dir+'SHA/check_images/'+name+thresh[col]+back[row]+'.fits')
	i1_seg = hdu[0].data
	if name == 'ngc7674' or name == 'ngc2768':
		i1_seg[np.where(i1_seg==0)]=-1
		i1_seg[np.where(i1_seg==i1_seg[int(j+10),int(i+10)])]=0
		i1_seg[np.where(i1_seg==i1_seg[int(j-10),int(i-10)])]=0
		i1_seg[np.where(i1_seg>0)]=1
	elif name == 'ngc4494':
                i1_seg[np.where(i1_seg==0)]=-1
                i1_seg[np.where(i1_seg==i1_seg[int(j+10),int(i-10)])]=0
                i1_seg[np.where(i1_seg==i1_seg[int(j-10),int(i+10)])]=0
                i1_seg[np.where(i1_seg>0)]=1
	else:
		i1_seg[np.where(i1_seg==0)]=-1
		i1_seg[np.where(i1_seg==i1_seg[int(j),int(i)])]=0
		i1_seg[np.where(i1_seg>0)]=1
		
	h=fig.add_axes([0.25*col,0.777-0.111*row,0.25,0.111])
        h.imshow(i1_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb])
	if (name == 'ngc4736' or name == 'ngc2403') and row == 1:
	        h.text(10,50,thresh1[col],color='white')
        	h.text(10,100,back[8],color='white')
	else:
		h.text(10,50,thresh[col],color='white')
                h.text(10,100,back[row],color='white')

	

Gal_ext=[]
gal_type=[]
betav_med=[]

#for x in range(0,len(sha_cat)):
for x in range(255,257):
        name    = sha_cat['col1'][x]      # NGC/IC name
	
	###### NED match ####
	ned_match = ned_tot[ned_tot['col2']==name]
	Gal_ext.append(ned_match['col6'][0])
	gal_type.append(ned_match['col13'][0])

	semi = ned_match['col11'][0]	# read semi-major axis from NED [']
#	thumb = int(semi*100)		# a half size of figure size
	thumb = int(semi*50)		# a half size of figure size

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

	fig = plt.figure(figsize=(12,18))
	
	h11 = fig.add_axes([0.00,0.888,0.25,0.111])
#	h12 = fig.add_axes([0.35,0.85,0.25,0.111])
#	h13 = fig.add_axes([0.7,0.85,0.25,0.111])

	img11 = make_lupton_rgb	(	i1[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]*2.33,	\
					rir[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],	\
					gir[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]	, Q=10, stretch=0.5)

	h11.imshow(img11)
	h11.text(10,50,'R : I1',color='white')
	h11.text(10,100,'G : r',color='white')
	h11.text(10,150,'B : g',color='white')

#	img1=h12.imshow(betav_raw[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],cmap=cmap,clim=[0,2])
#	h12.text(10,50,'V/I1',color='blue')
#	colorbar(img1)

#       img2=h13.imshow(betag[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],cmap=cmap,clim=[0,2])
#        h13.text(10,50,'g/I1',color='blue')
#        colorbar(img2)

	for col in range(0,4):
		for row in range(0,8):
			draw(name,fig,i,j,row,col,thumb)

	fig.savefig(work_dir+'pdf/'+name+'.thresh.back.pdf')


