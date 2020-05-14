#####################################################
# Python script of building master PSF from rotation 
# and magnitude matched PSFs 
# written by Duho Kim (12/06/17)        
######################################################
from astropy.io import fits
from pyraf import iraf

data_dir = '/Users/dhk/work/data/NGC_IC/SHA/'

psfs = ['ic195','ic208','ic486','ic505','ic5298','ic692','ngc1265','ngc1700','ngc2552','ngc2604','ngc2906', \
	'ngc3184','ngc3192','ngc3212','ngc3215','ngc3310','ngc3351','ngc3395','ngc3419','ngc3622','ngc3741', \
	'ngc3921','ngc3928','ngc4369','ngc4412','ngc4420','ngc4449','ngc4457','ngc4561','ngc4580','ngc4651', \
	'ngc5068','ngc5127','ngc5273','ngc5313','ngc5363','ngc5515','ngc5520','ngc5576','ngc5596','ngc5668', \
	'ngc5820','ngc5936','ngc596','ngc5992','ngc6166','ngc6340','ngc636','ngc7077','ngc7080','ngc7177', \
	'ngc7385','ngc7674','ngc7742','ngc7743','ngc777','ngc985']

hdu = fits.open(data_dir+'psf4_rotate/med.fits')	# read master PSF file
master = hdu[0].data					# 

for i in range(0,len(psfs)):
	tmp = fits.open(data_dir+'psf4_rotate/'+psfs[i]+'_rot.fits')
	psf = tmp[0].data
	div = psf/master
	fits.writeto(data_dir+'psf4_div/'+psfs[i]+'_div.fits',div)

