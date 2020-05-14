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
import numpy.ma as ma
from astropy.wcs import WCS
from astropy.table import Table, vstack
from astropy.visualization import make_lupton_rgb
from astropy.cosmology import Planck15 as cosmo
import os
import os.path
import sys
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import colors
import matplotlib.image as mpimg
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.coordinates import SkyCoord
import math
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd
from operator import itemgetter
import copy

AV_fit = False
bv_fit = False

cmap = 'gist_rainbow'	### color map for betaV images
cmap2 = colors.ListedColormap(['white','cyan','gray'])
pscale = 0.6		### pixel scale ["]
snr = 3			### signal to noise ratio that we constrain

#bv0s = [0.645,0.64,0.695,0.74,0.895,1.075,1.41]			### beta_V,0 for E0,S0(SFH3),Sa,Sb,Sbc(SFH4),Sc,Sd(SFH5) from Table 5 Kim+17
bv0s=[   [1.9409089,1.9044532,1.3527486,1.1198820,0.82968015,0.58503551],        # BetaVzero values for SFH3,4,5 and Z=.0001,.0004,.004,.008,.02,.05
        [1.9860456,1.9576220,1.4390702,1.2023316,0.91737698,0.65453906],
        [2.3880801,2.3478914,1.6838646,1.4124115,1.1048444,0.77272439]]


htypes = ['E0','S0','S0/a','Sa','Sab','Sb','Sbc','Sc','Scd','Sd','Sdm','Sm','Im','?'] # Hubble types corresponding Nair+10 T-Type
agntypes = ['SF','transition/mixed AGN','Seyfert','LINER']	# Kauffmann+03 from Nair+10

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

ned1=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5=ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")
# row | input object name | object name | RA | Dec | Gal. Ext. Burstein & Heiles A_B mag | Object Type 
# (1) |       (2)         |    (3)      | (4)| (5) |                (6)                  |    (7)
# Redshift | Redshift Uncertainty | Mag./Filter | Major Diam | Minor Diam | Morph. | Ref.| 
# (8)      |        (9)           |    (10)     |   (11)     |    (12)    |  (13)  | (14)|
ned_tot=vstack([ned1,ned2,ned3,ned4,ned5])				# Total NED output catalog


#sha_cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')
sha_cat=ascii.read(cat_dir+'sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv')	# Sample NGC/IC numbers
sha2fig=np.load('/Users/dhk/work/data/NGC_IC/pdf/sha2fig.npy')		# load figure numbers for each galaxies

rc3=ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
rosa_AV=np.genfromtxt("/Users/dhk/work/data/rosa/AV.txt")		# read A_V values for different Hubble types from CALIFA paper

file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'			# read NGC/IC catalogue 
df = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'				# read T-type
ttype = pd.read_excel(file)

####################### Jarrett+03 ########################################################
# Modified logarithmic visualization method P' = sqrt( log(1+P/(n*sig))  ), where P is the pixel intensity value, sig is the image rms "noise",
# and n is a threshold throttle (w/ satisfactory values between 5 and 10
def Jarrett(pix_values,sigma):
	n = 7.5
	return np.sqrt(np.log10(1+pix_values/(n*sigma)))
############################################################################################

gal_type=[]
betav_med=[]

betav_257=np.zeros(len(sha_cat))
betav_257_std=np.zeros(len(sha_cat))
ttype_257=np.zeros(len(sha_cat))
ttype_257_err=np.zeros(len(sha_cat))
betav_257_center=np.zeros(len(sha_cat))
betav_257_center_std=np.zeros(len(sha_cat))

######################### SURFACE BRIGHTNESS MEASURE #######################################
def isophote(fits,cx,cy,pa,boa,rad_min,rad_max,rad_num):
	fits[:] = float('nan')
	yy,xx           = np.ogrid[cy-rad_max:cy+rad_max,cx-rad_max:cx+rad_max]         # opposite convention between FITS & Python array indexing
	r               = np.logspace(np.log10(rad_min),np.log10(rad_max),rad_num+1)
	iso_mean        = np.zeros(rad_num)
	iso_med         = np.zeros(rad_num)
	iso_std         = np.zeros(rad_num)
	tot_mean        = np.zeros(rad_num)
	tot_med         = np.zeros(rad_num)
	pa_rad          = math.radians(-pa)                                           # convert PA convention to functional form with units
	for ii in range(1,rad_num+1):
		# General Equation of an Ellipse : (x*cos(A) + y*sin(A))^2/a^2 + (x*sin(A)-y*cos(A))^2/b^2 = 1
		ind1    = 1     <       ((yy-cy)*math.cos(pa_rad)+(xx-cx)*math.sin(pa_rad))**2/(r[ii-1]+1)**2 + ((yy-cy)*math.sin(pa_rad)-(xx-cx)*math.cos(pa_rad))**2/((r[ii-1]+1)*boa)**2
		ind2    = 1     >       ((yy-cy)*math.cos(pa_rad)+(xx-cx)*math.sin(pa_rad))**2/(r[ii]-1)**2 + ((yy-cy)*math.sin(pa_rad)-(xx-cx)*math.cos(pa_rad))**2/((r[ii]-1)*boa)**2
		ind     = ind1 * ind2
		anul    = fits[ind]
		iso_mean[ii-1]   = np.nanmean(anul)
		iso_med[ii-1]    = np.nanmedian(anul)
		iso_std[ii-1]    = np.nanstd(anul)
		fits[ind] = ii*10
		if ii==1:
			tot_mean[0] = iso_mean[ii-1]*np.sum(ind)
			tot_med[0] = iso_med[ii-1]*np.sum(ind)
		else:
			tot_mean[ii-1] = tot_mean[ii-2] + iso_mean[ii-1]*np.sum(ind)
			tot_med[ii-1] = tot_med[ii-2] + iso_med[ii-1]*np.sum(ind)
	return fits, r[1:], iso_mean, iso_med, iso_std, tot_mean, tot_med
#############################################################################################
	

#for x in range(0,len(sha_cat)):
for x in range(0,5):
        name    = sha_cat['col1'][x]      # NGC/IC name
	
	###### NGC match ####
	if name[0]=='n':
		galnum = name[3:].strip()
		ngc_match = df.loc[(df['N']=='N') & (df['NI']==int(galnum))]
		if len(galnum) == 3:
			galnum='0'+galnum
		elif len(galnum) ==2:
			galnum='00'+galnum
		elif len(galnum) ==1:
			galnum='000'+galnum
		rc3name = 'NGC'+galnum
		table_name = 'NGC '+galnum
	elif name[0]=='i':
		galnum = name[2:].strip()
		ngc_match = df.loc[(df['N']=='I') & (df['NI']==int(galnum))]
		if len(galnum) == 3:
			galnum='0'+galnum
		elif len(galnum)==2:
			galnum='00'+galnum
		elif len(galnum)==1:
			galnum='000'+galnum
		rc3name = 'IC'+galnum
		table_name = 'IC '+galnum
	
	pa = ngc_match.iloc[0,21]
        if math.isnan(pa):
                pa=0
        elif pa > 90:
                pa=pa-180.0

	####  RC3 catalogue #########
	rc3_match = rc3[[j for j,s in enumerate(rc3['name']) if s.strip() == rc3name]]
	if len(rc3_match) != 1:
		print('rc3 match is not one')
		ttype_257_err[x] = 1
	elif rc3_match['eT'][0] == '*':
		ttype_257_err[x] = 1
	else:
		ttype_257_err[x] = rc3_match['eT'][0]

	###### T-type match ####
	ttype_match = ttype.loc[ttype['name']==table_name]
	T = ttype_match.iloc[0,5]
	ttype_257[x] = T
	if T < -3:
		bv0 = bv0s[0][5]		# for Dust profile on Figure
		bv0Zs = bv0s[0]			# for Dust profiles 6 Zs
		AV_prof=rosa_AV[0][2:]
	if T < 0 and T >= -3:
		bv0 = bv0s[0][5]
		bv0Zs = bv0s[0]
		AV_prof=rosa_AV[2][2:]
	if T < 2 and T >=0:
		bv0 = bv0s[1][4]
		bv0Zs = bv0s[1]
		AV_prof=rosa_AV[4][2:]
	if T < 4 and T >=2:
		bv0 = bv0s[1][4]
		bv0Zs = bv0s[1]
		AV_prof=rosa_AV[6][2:]
	if T < 5 and T >=4:
		bv0 = bv0s[1][3]
		bv0Zs = bv0s[1]
		AV_prof=rosa_AV[8][2:]
	if T < 7 and T >=5:
		bv0 = bv0s[2][3]
		bv0Zs = bv0s[2]
		AV_prof=rosa_AV[10][2:]
	if T < 9 and T >=7:
		bv0 = bv0s[2][3]
		bv0Zs = bv0s[2]
		AV_prof=rosa_AV[12][2:]
	if T>=9:
		bv0 = -1		# no BetaVzero is available for these types
	
	###### NED match ####
	ned_match = ned_tot[ned_tot['col2']==name]
	if x == 90:
		GalExtFactor=10.0**(0.072/2.5)				# NGC 4561, A_V=0.072	(NED)
	else:
		GalExtFactor=10.0**(ned_match['col6'][0]/2.5)		# read Galactic extinction from NED [Vint/Vobs]
	gal_type.append(ned_match['col13'][0])

	semi = ned_match['col11'][0]		# read semi-major axis from NED [']
	thumb = int(semi*100)			# a half size of figure size
	z = ned_match['col8'][0]		# read redshift from NED	
	boa = ned_match['col12'][0] / ned_match['col11'][0]
	
	##### WCS read  #####
	w	= WCS(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.fits')
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
	hdu     = fits.open(work_dir+'SHA/SHA_NGC_IC_LONG_mean/'+name+'.std.fits')
        i1u     = hdu[0].data

	###### Calculate V-band ##########
        V = 10**(np.log10(gir)+0.59*(np.log10(rir/gir))+0.004) 	# Jester+05 All stars R-I < 1.15

	###### Calculate betaV, betag ####
	betav_raw	=	V/(i1*2.330451129)*GalExtFactor		# correct for zero point difference between SDSS & Spitzer
	betag_raw	=	gir/(i1*2.330451129)*GalExtFactor	# and Galactic Extinction

	###### READ SE SEGMENTATION MAP ####
	hdu = fits.open(work_dir+'SHA/segmap/'+name+'.fits')
	i1_seg = hdu[0].data

	####### SEGMENATION & select SNR > 3 pixels #################	
	betav_seg	= 	V/(i1*2.330451129)*GalExtFactor
	betav_seg[np.where(i1_seg!=1)] = float('nan')	
	V_seg = copy.deepcopy(V)	
	V_seg[np.where(i1_seg!=1)] = float('nan')

	gir_cent = gir[int(j)-500:int(j)+500,int(i)-500:int(i)+500]			# only use central region due to stripe regions after registering
	gir_seg_cent = gir_seg[int(j)-500:int(j)+500,int(i)-500:int(i)+500]		
        rir_cent = rir[int(j)-500:int(j)+500,int(i)-500:int(i)+500]
        rir_seg_cent = rir_seg[int(j)-500:int(j)+500,int(i)-500:int(i)+500]

        gerr    = np.nanstd(gir_cent[np.where(gir_seg_cent==0)])
        rerr    = np.nanstd(rir_cent[np.where(rir_seg_cent==0)])
	betav_seg[(i1/i1u < snr) | (gir/gerr < snr) | (rir/rerr < snr)] = float('nan')

	betav_257[x] = np.nanmean(betav_seg)
	betav_257_std[x] = np.nanstd(betav_seg)

	betav_257_center[x] = np.nanmean(betav_seg[int(j)-1:int(j)+2,int(i)-1:int(i)+2])
	betav_257_center_std[x] = np.nanstd(betav_seg[int(j)-1:int(j)+2,int(i)-1:int(i)+2])

	####### PLOT FIGURES for each sample ###############
	
	fig = plt.figure(figsize=(12,7.8))
	
	h11 = fig.add_axes([0.0,0.5,0.315,0.5])
	h12 = fig.add_axes([0.315,0.5,0.315,0.5])
	h13 = fig.add_axes([0.63,0.5,0.331,0.5])

	h21 = fig.add_axes([0.0,0.0,0.315,0.5])
	h22 = fig.add_axes([0.315,0.0,0.315,0.5])

	minsize = min([min(i1.shape),min(rir.shape),min(gir.shape)]) 
	if thumb*2 > minsize:
		thumb=int(minsize/2)

	img11 = make_lupton_rgb	(	i1[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]*2.33,	\
					rir[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],	\
					gir[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]	, Q=10, stretch=0.5)

	h11.imshow(img11,origin='lower')
	h11.get_xaxis().set_visible(False)
	h11.get_yaxis().set_visible(False)
	h11.text(0.02,0.95,r'R : IRAC 3.6$\mu$m',color='white',fontsize=15,transform=h11.transAxes,fontweight='bold')
	h11.text(0.02,0.90,'G : SDSS r',color='white',fontsize=15,transform=h11.transAxes,fontweight='bold')
	h11.text(0.02,0.85,'B : SDSS g',color='white',fontsize=15,transform=h11.transAxes,fontweight='bold')
	h11.text(0.6,0.95,table_name,color='white',fontsize=15,transform=h11.transAxes,fontweight='bold')
	h11.text(0.3,0.05,ned_match['col13'][0],color='white',fontsize=15,transform=h11.transAxes,fontweight='bold')
		
	kpc_arcmin = cosmo.kpc_proper_per_arcmin(z)     # read xxx kpc / arcmin
	arcmin_5kpc = 5.0/kpc_arcmin			# calculate xxx arcmin / 5 kpc
	frac_5kpc = arcmin_5kpc*100.0/(2*thumb)		# calculate fraction of a length of 5 kpc in fig size
	h11.plot([0.05,0.05+frac_5kpc.value],[0.02,0.02],color='white',transform=h11.transAxes)
	h11.text(0.02,0.05,'5kpc, '+'{:4.2f}'.format(arcmin_5kpc.value)+'\'',color='white',fontsize=12,transform=h11.transAxes)

	img1=h12.imshow(betav_raw[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],cmap=cmap,clim=[0,2],origin='lower')
        h12.get_xaxis().set_visible(False)
        h12.get_yaxis().set_visible(False)
	h12.text(0.02,0.9,r'$\beta_{V}$ (V/3.6$\mu$m)',color='black',fontsize=18,transform=h12.transAxes,fontweight='bold')

        img2=h13.imshow(betag_raw[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],cmap=cmap,clim=[0,2],origin='lower')
	img2.axes.figure.colorbar(img2,cax=make_axes_locatable(img2.axes).append_axes("right",size="5%",pad=0.0))
        h13.get_xaxis().set_visible(False)
        h13.get_yaxis().set_visible(False)
        h13.text(0.02,0.9,r'$\beta_{g}$ (g/3.6$\mu$m)',color='black',fontsize=18,transform=h13.transAxes,fontweight='bold')

	i1_seg[np.where(i1_seg > 2)]=2
	h21.imshow(i1_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],cmap=cmap2,origin='lower')
        h21.get_xaxis().set_visible(False)
        h21.get_yaxis().set_visible(False)
	h21.text(0.02,0.95,'Segmentation',fontsize=15,transform=h21.transAxes,fontweight='bold')
	h21.text(0.02,0.90,'& Mask out',fontsize=15,transform=h21.transAxes,fontweight='bold')
	
	h22.imshow(betav_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],cmap=cmap,clim=[0,2],origin='lower')
        h22.get_xaxis().set_visible(False)
        h22.get_yaxis().set_visible(False)
	h22.text(0.02,0.95,'S/N > 3',fontsize=15,transform=h22.transAxes,fontweight='bold')

	dust_map = 2.5*np.log10(bv0/betav_seg)
	dust_map[np.where(betav_seg > bv0)] = 0.0			# assign A_V = 0 for bv > bv0

	########## surface profile measure #############
	fits_anul, rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(V_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],int(i),int(j),pa,boa,2,thumb,11)

	half_light = max(tot_med)/2.0
	hlr_idx = min(enumerate(np.abs(tot_med-half_light)),key=itemgetter(1))[0] 		# find half-light radius
	hlr3_idx = min(enumerate(np.abs(rads[hlr_idx]*3.0-rads)),key=itemgetter(1))[0]		# find half-light radius*3.0

	h25 = fig.add_axes([0.843-0.3,0.09,0.1,0.1])
	h25.plot(rads,tot_med)
	h25.plot([rads[hlr_idx],rads[hlr_idx]],[0,max(tot_med)])

	if bv0 != -1:
		h23 = fig.add_axes([0.63,0.0,0.331,0.5])
		img3=h23.imshow(fits_anul,origin='lower')
#		img3.axes.figure.colorbar(img3,cax=make_axes_locatable(img3.axes).append_axes("right",size="5%",pad=0.0))
#	        h23.text(0.02,0.95,'DUST MAP (Av)',fontsize=15,transform=h23.transAxes,fontweight='bold')
#        	h23.text(0.02,0.90,r'2.5$\times$log('+'{:5.3f}'.format(bv0)+r'/$\beta_{V}$)',fontsize=15,transform=h23.transAxes,fontweight='bold')
		h23.get_xaxis().set_visible(False)
        	h23.get_yaxis().set_visible(False)

		ellipse1 = Ellipse((thumb,thumb),rads[hlr_idx],rads[hlr_idx]*boa,pa+90,fill=False)
		ellipse2 = Ellipse((thumb,thumb),rads[hlr3_idx],rads[hlr3_idx]*boa,pa+90,fill=False)
		h23.add_artist(ellipse1)
		h23.add_artist(ellipse2)	

		######## surface betaV & A_V profile measure ##########
		if bv_fit:
			rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(betav_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],int(i),int(j),pa,boa,2,thumb,11)
			np.savetxt(work_dir+'ellipse/scratch/'+name+'_bv.txt',(rads,iso_mean,iso_med,iso_std))
		
		if AV_fit:
                        rads, iso_mean_map, iso_med_map, iso_std_map, tot_mean, tot_med = isophote(dust_map[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb],int(i),int(j),pa,boa,2,thumb,11)
                        np.savetxt(work_dir+'ellipse/scratch/'+name+'_AV.txt',(rads,iso_mean,iso_med,iso_std))
			rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(2.5*np.log10(bv0Zs[0]/betav_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]),int(i),int(j),pa,boa,2,thumb,11)
			np.savetxt(work_dir+'ellipse/scratch/'+name+'_AV1.txt',(rads,iso_mean,iso_med,iso_std))
                        rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(2.5*np.log10(bv0Zs[1]/betav_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]),int(i),int(j),pa,boa,2,thumb,11)
                        np.savetxt(work_dir+'ellipse/scratch/'+name+'_AV2.txt',(rads,iso_mean,iso_med,iso_std))
                        rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(2.5*np.log10(bv0Zs[2]/betav_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]),int(i),int(j),pa,boa,2,thumb,11)
                        np.savetxt(work_dir+'ellipse/scratch/'+name+'_AV3.txt',(rads,iso_mean,iso_med,iso_std))
                        rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(2.5*np.log10(bv0Zs[3]/betav_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]),int(i),int(j),pa,boa,2,thumb,11)
                        np.savetxt(work_dir+'ellipse/scratch/'+name+'_AV4.txt',(rads,iso_mean,iso_med,iso_std))
                        rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(2.5*np.log10(bv0Zs[4]/betav_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]),int(i),int(j),pa,boa,2,thumb,11)
                        np.savetxt(work_dir+'ellipse/scratch/'+name+'_AV5.txt',(rads,iso_mean,iso_med,iso_std))
                        rads, iso_mean, iso_med, iso_std, tot_mean, tot_med = isophote(2.5*np.log10(bv0Zs[5]/betav_seg[int(j)-thumb:int(j)+thumb,int(i)-thumb:int(i)+thumb]),int(i),int(j),pa,boa,2,thumb,11)
                        np.savetxt(work_dir+'ellipse/scratch/'+name+'_AV6.txt',(rads,iso_mean,iso_med,iso_std))

			h24 = fig.add_axes([0.843,0.39,0.1,0.1])
			h24.plot(rads/rads[hlr_idx],iso_med_map,color='red')
			h24.fill_between(rads/rads[hlr_idx],iso_med_map-iso_std_map,iso_med_map+iso_std_map,color='red',alpha=0.3)
			h24.plot(np.arange(0.,3.,0.1),AV_prof[1:],color='black',alpha=0.5)
			h24.fill_between(np.arange(0.,3.,0.1),AV_prof[1:]-AV_prof[0],AV_prof[1:]+AV_prof[0],alpha=0.2,color='black')
			plt.xticks([1,3])
			h24.set(xlabel='HLR',ylabel='Av')
        else:
		h23 = fig.add_axes([0.63,0.0075,0.331,0.485])
		h23.get_yaxis().set_visible(False)
		h23.get_xaxis().set_visible(False)
                h23.text(0.02,0.9,'No '+r'$\beta_{V,0}$ is available for this T-type',size=15)


	fig.savefig(work_dir+'pdf/'+name+'.tex.pdf')
	plt.close(fig)

	if int((sha2fig[x]-1)%2):
		copyfile(work_dir+'pdf/'+name+'.tex.pdf','/Users/dhk/Documents/publish/ngcic/figset2_'+str(int((sha2fig[x]-1)/2))+'b_test.pdf')
	else:
		copyfile(work_dir+'pdf/'+name+'.tex.pdf','/Users/dhk/Documents/publish/ngcic/figset2_'+str(int((sha2fig[x]-1)/2))+'a_test.pdf')


#np.save('betav_ttype.npy',[betav_257,betav_257_std,ttype_257,ttype_257_err,betav_257_center,betav_257_center_std])
