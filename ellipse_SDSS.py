#####################################################
# Python script of PyRAF ELLIPSE fitting for given FITS images 
# of SDSS data
# written by Duho Kim (05/17/18)	
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
from astropy.wcs import WCS
from astropy.table import Table, vstack
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd
from copy import deepcopy

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'
cur_dir = '/Users/dhk/work/py/pyraf/'

ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")
# row | input object name | object name | RA | Dec | Gal. Ext. Burstein & Heiles A_B mag | Object Type 
# (1) |       (2)         |    (3)      | (4)| (5) |                (6)                  |    (7)
# Redshift | Redshift Uncertainty | Mag./Filter | Major Diam | Minor Diam | Morph. | Ref.| 
# (8)      |        (9)           |    (10)     |   (11)     |    (12)    |  (13)  | (14)|
ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])                              # Total NED output catalog

file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'
df = pd.read_excel(file)

sha_cat=ascii.read(cat_dir+'sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv') # Sample NGC/IC numbers

iraf.stsdas()
iraf.analysis()
iraf.isophote()

#for x in range(0,len(sha_cat)):
for x in range(0,1):
	name    = sha_cat['col1'][x]      # NGC/IC name

	###### NGC match ####
	if name[0]=='n':
		galnum = name[3:].strip()
		ngc_match = df.loc[(df['N']=='N') & (df['NI']==int(galnum))]
	elif name[0]=='i':
		galnum = name[2:].strip()
		ngc_match = df.loc[(df['N']=='I') & (df['NI']==int(galnum))]


        ###### NED match ####
        ned_match = ned_tot[ned_tot['col2']==name]

        semi = ned_match['col11'][0]            # read semi-major axis from NED [']
        thumb = int(semi*100)                   # a half size of figure size
	boa = ned_match['col12'][0] / ned_match['col11'][0]
	ellip = 1-boa
	if ellip < 0.05:
		ellip = 0.05

	##### READ FITS #####
	hdu     = fits.open(work_dir+'SDSS/g/'+name+'-gir.fits')
        gir      = hdu[0].data
        hdu     = fits.open(work_dir+'SDSS/r/'+name+'-rir.fits')
        rir     = hdu[0].data

        ###### Calculate V-band ##########
        V = 10**(np.log10(gir)+0.59*(np.log10(rir/gir))+0.004)  # Jester+05 All stars R-I < 1.15

        ######### READ WCS ###########
        w       = WCS(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')
        i,j     = w.all_world2pix(sha_cat['col2'][x],sha_cat['col3'][x],1)

        ###### READ SE SEGMENTATION MAP ####
	hdu = fits.open(work_dir+'SHA/segmap/'+name+'.fits')
        i1_seg = hdu[0].data

	mask = deepcopy(i1_seg)
	mask[np.where(mask==0)]=99
	mask[np.where(mask==1)]=0
	mask[np.where(mask!=0)]=32

	if os.path.exists(work_dir+'SHA/segmap/'+name+'.dqf'):
		os.remove(work_dir+'SHA/segmap/'+name+'.dqf')

	
	hdu_dqf = fits.PrimaryHDU(mask)
	hdu_dqf.writeto(work_dir+'SHA/segmap/'+name+'.dqf')

	###### Erase files ##########
	if os.path.exists(cur_dir+'tmp_V.fits'):
		os.remove(cur_dir+'tmp_V.fits')
	if os.path.exists(cur_dir+'tmp_I1.tab'):
                os.remove(cur_dir+'tmp_I1.tab')
        if os.path.exists(cur_dir+'tmp_V.tab'):
                os.remove(cur_dir+'tmp_V.tab')

	###### SAVE FITS for ELLIPSE run ########
	if name != 'ngc7674':						# due to disability of finding center in ELLIPSE fit attempt
		V[np.where(i1_seg!=1)]=float('nan')
	hdu1=fits.PrimaryHDU(V)
	hdu1.writeto(cur_dir+'tmp_V.fits')

	###### Run PyRAF ELLIPSE #########
	pa = ngc_match.iloc[0,21]
	if math.isnan(pa):
		pa=0
		ellip=0.05
	elif pa > 90:
		pa=pa-180.0
#	iraf.ellipse.controlpar.minit = 10
#	iraf.ellipse.controlpar.maxit = 100
#	iraf.ellipse.controlpar.olthresh = 0.0


#	iraf.ellipse(input=cur_dir+'tmp_V.fits',output=cur_dir+'tmp_V.tab',inter='no',x0=i,y0=j,ellip0=ellip,pa0=pa,sma0=thumb/10.0, \
#		minsma=1,maxsma=thumb,linear='no',recenter='no',hcenter='yes',hellip='yes',hpa='yes',integrmode='median',refer=1.0, \
#		zerolevel=0.0,mag0=22.5,nclip=0,step=0.3)

        iraf.ellipse(input=cur_dir+'tmp_V.fits',output=cur_dir+'tmp_V.tab',inter='no',x0=i,y0=j,ellip0=ellip,pa0=pa,sma0=thumb/10.0, \
                minsma=1,maxsma=thumb,linear='no',recenter='yes',hcenter='no',hellip='no',hpa='no',integrmode='median',refer=1.0, \
                zerolevel=0.0,mag0=22.5,nclip=0,step=0.3)


	iraf.tdump(table=cur_dir+'tmp_V.tab',datafile='/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'.txt',columns='@colnames.lis')
	copyfile(cur_dir+'tmp_V.tab','/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'.tab')

	iraf.ellipse(input=work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits', output=cur_dir+'tmp_I1.tab',inter='no',inellip=cur_dir+'tmp_V.tab',mag0=21.5814)
	iraf.tdump(table=cur_dir+'tmp_I1.tab',datafile='/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'_I1.txt',columns='@colnames.lis')
	copyfile(cur_dir+'tmp_I1.tab','/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'_I1.tab')


#input=work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits', output=cur_dir+'tmp_I1.tab',inter='no',x0=i,y0=j,ellip0=ellip,pa0=pa,sma0=thumb/10.0, \
#		minsma=1,maxsma=thumb,linear='no',recenter='no',hcenter='yes',hellip='yes',hpa='yes',integrmode='median',refer=1.0, \
#		zerolevel=0.0,mag0=21.5814,nclip=0,step=0.3)
#dqf=work_dir+'SHA/segmap/'+name+'.dqf')

#       iraf.tdump(table=cur_dir+'tmp_I1.tab',datafile='/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'_I1.txt',columns='@colnames.lis')
#       copyfile(cur_dir+'tmp_I1.tab','/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'_I1.tab')

#	iraf.ellipse(input=cur_dir+'tmp_V.fits', output=cur_dir+'tmp_V.tab',inter='no',inellip=cur_dir+'tmp_I1.tab',mag0=22.5,dqf=work_dir+'SHA/segmap/'+name+'.dqf')

#	iraf.tdump(table=cur_dir+'tmp_V.tab',datafile='/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'.txt',columns='@colnames.lis')
#	copyfile(cur_dir+'tmp_V.tab','/Users/dhk/work/data/NGC_IC/ellipse/pa_ellip/'+name+'.tab')


