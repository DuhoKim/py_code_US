#####################################################
# Python script that reads AV profiles of 257 galaxies
# and stack as a function of SFHs (SFH3,4,5 for each
# E S0, Sa-Sbc, Sc-Sd) and Z (.0001.0004.004.008.02.05)
# written by Duho Kim (06/28/18)        
######################################################
from astropy.io import ascii
import numpy as np
from astropy.table import Table, vstack
import os
import os.path
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pandas import DataFrame, read_csv
import pandas as pd
from operator import itemgetter
from scipy.interpolate import interp1d

htypes = ['E0','S0','S0/a','Sa','Sab','Sb','Sbc','Sc','Scd','Sd','Sdm','Sm','Im','?'] # Hubble types corresponding Nair+10 T-Type
agntypes = ['SF','transition/mixed AGN','Seyfert','LINER']	# Kauffmann+03 from Nair+10

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

text_for_Z=[r'Z=0.0001','Z=0.0004','Z=0.004','Z=0.008','Z=0.02','Z=0.05']
cols_for_Z=['black','violet','blue','green','darkorange','red']

sha_cat=ascii.read(cat_dir+'sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv')	# Sample NGC/IC numbers
sha2fig=np.load('/Users/dhk/work/data/NGC_IC/pdf/sha2fig.npy')		# load figure numbers for each galaxies

rc3=ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
bt=ascii.read("/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask/btot.txt")
rosa_AV=np.genfromtxt("/Users/dhk/work/data/rosa/AV.txt")		# read A_V values for different Hubble types from CALIFA paper

file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'			# read NGC/IC catalogue 
df = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'				# read T-type
ttype = pd.read_excel(file)

av_prof_sfh3 = np.zeros((84,6,50))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_sfh4 = np.zeros((90,6,50))        # For stacking Dust profiles for SFH4 w/ 6 Zs
av_prof_sfh41 = np.zeros((90,6,50))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T < 0.5
av_prof_sfh42 = np.zeros((90,6,50))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T > 0.5
av_prof_sfh5 = np.zeros((60,6,50))        # For stacking Dust profiles for SFH5 w/ 6 Zs
av_prof_sfh51 = np.zeros((60,6,50))        # For stacking Dust profiles for SFH5 w/ 6 Zs with B/T < 0.5
av_prof_sfh52 = np.zeros((60,6,50))        # For stacking Dust profiles for SFH5 w/ 6 Zs with B/T > 0.5

cnt_sfh3 = 0	# counter for SFH3
cnt_sfh4 = 0	# counter for SFH4
cnt_sfh41 = 0	# counter for SFH4 with B/T < 0.5
cnt_sfh42 = 0	# counter for SFH4 with B/T > 0.5
cnt_sfh5 = 0	# counter for SFH5
cnt_sfh51 = 0	# counter for SFH5 with B/T < 0.5
cnt_sfh52 = 0	# counter for SFH5 with B/T > 0.5

for x in range(0,len(sha_cat)):
#for x in range(0,1):
        name    = sha_cat['col1'][x]      # NGC/IC name
	if name[0]=='n':
		galnum = name[3:].strip()
		
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
		if len(galnum) == 3:
			galnum='0'+galnum
		elif len(galnum)==2:
			galnum='00'+galnum
		elif len(galnum)==1:
			galnum='000'+galnum
		rc3name = 'IC'+galnum
		table_name = 'IC '+galnum

	###### HLR Calc ########
        ell = ascii.read(work_dir+'ellipse/pa_ellip/'+name+'.txt')
        # col1:SMA, col2:INTENS, col3:INT_ERR, col4:RMS, col5:ELLIP, col6:ELLIP_ERR, col7:RSMA, col8:MAG, col9:MAG_LERR, col10:MAG_UERR
        # col11:TFLUX_E, col12:TMAG_E, col13:PA, col14:PA_ERR, col15:X0, col16:X0_ERR, col17:Y0, col18:Y0_ERR

        min_idx=0
        hlr_idx=0

        if ell['col2'][0] != 0:
                min_idx = min(enumerate(ell['col2']),key=itemgetter(1))[0]                                      # find the index where INTENS is minimum
                hlr_idx = min(enumerate(np.abs(ell['col11'][min_idx]/2.0-ell['col11'])),key=itemgetter(1))[0]   # find half-light radius
        else:
                half_light = max(ell['col11'])/2.0
                hlr_idx = min(enumerate(np.abs(ell['col11']-half_light)),key=itemgetter(1))[0]   # find half-light radius

        hlr3_idx = min(enumerate(np.abs(ell['col1'][hlr_idx]*3.0-ell['col1'])),key=itemgetter(1))[0]   # find half-light radius*3.0
        hlr_pix = ell['col1'][hlr_idx]	

	###### T-type match ####
	ttype_match = ttype.loc[ttype['name']==table_name]
	T = ttype_match.iloc[0,5]
		
	###### B/T match #######
	ind = np.where(bt['col1']==name)
	if ind[0] >= 0:
		btot = bt['col2'][ind[0]]
		chi = 0.5
	else:
		btot = -99
		chi = -99
	#btot = ttype_match.iloc[0,7]
	#chi = ttype_match.iloc[0,8]

	if T >=9:
		continue

	ell1=ascii.read(work_dir+'ellipse/pa_ellip/'+name+'_AV1.txt')
	ell2=ascii.read(work_dir+'ellipse/pa_ellip/'+name+'_AV2.txt')
	ell3=ascii.read(work_dir+'ellipse/pa_ellip/'+name+'_AV3.txt')
	ell4=ascii.read(work_dir+'ellipse/pa_ellip/'+name+'_AV4.txt')
	ell5=ascii.read(work_dir+'ellipse/pa_ellip/'+name+'_AV5.txt')
	ell6=ascii.read(work_dir+'ellipse/pa_ellip/'+name+'_AV6.txt')

	xnew = np.linspace(0,5,num=50)
	x = ell1['col1']/hlr_pix	# SMA	

	y1 = ell1['col2']		# INTENS
	y2 = ell2['col2']		# INTENS
	y3 = ell3['col2']		# INTENS
	y4 = ell4['col2']		# INTENS
	y5 = ell5['col2']		# INTENS
	y6 = ell6['col2']		# INTENS
	
	f1 = interp1d(x,y1,kind='cubic',fill_value='extrapolate')	
	f2 = interp1d(x,y2,kind='cubic',fill_value='extrapolate')	
	f3 = interp1d(x,y3,kind='cubic',fill_value='extrapolate')	
	f4 = interp1d(x,y4,kind='cubic',fill_value='extrapolate')	
	f5 = interp1d(x,y5,kind='cubic',fill_value='extrapolate')	
	f6 = interp1d(x,y6,kind='cubic',fill_value='extrapolate')	

	y1new = f1(xnew)
	y2new = f2(xnew)
	y3new = f3(xnew)
	y4new = f4(xnew)
	y5new = f5(xnew)
	y6new = f6(xnew)

	y1new[np.where(y1new < 0)] = 0.0
	y2new[np.where(y2new < 0)] = 0.0
	y3new[np.where(y3new < 0)] = 0.0
	y4new[np.where(y4new < 0)] = 0.0
	y5new[np.where(y5new < 0)] = 0.0
	y6new[np.where(y6new < 0)] = 0.0

	if T < 0:
		av_prof_sfh3[cnt_sfh3,0,:] = y1new
		av_prof_sfh3[cnt_sfh3,1,:] = y2new
		av_prof_sfh3[cnt_sfh3,2,:] = y3new
		av_prof_sfh3[cnt_sfh3,3,:] = y4new
		av_prof_sfh3[cnt_sfh3,4,:] = y5new
		av_prof_sfh3[cnt_sfh3,5,:] = y6new
		cnt_sfh3 = cnt_sfh3 + 1
	if T < 5 and T >= 0:
                av_prof_sfh4[cnt_sfh4,0,:] = y1new
                av_prof_sfh4[cnt_sfh4,1,:] = y2new
                av_prof_sfh4[cnt_sfh4,2,:] = y3new
                av_prof_sfh4[cnt_sfh4,3,:] = y4new
                av_prof_sfh4[cnt_sfh4,4,:] = y5new
                av_prof_sfh4[cnt_sfh4,5,:] = y6new
                cnt_sfh4 = cnt_sfh4 + 1
        if T < 5 and T >= 0 and btot > 0.0 and btot < 0.5 and chi < 1.0:
                av_prof_sfh41[cnt_sfh41,0,:] = y1new
                av_prof_sfh41[cnt_sfh41,1,:] = y2new
                av_prof_sfh41[cnt_sfh41,2,:] = y3new
                av_prof_sfh41[cnt_sfh41,3,:] = y4new
                av_prof_sfh41[cnt_sfh41,4,:] = y5new
                av_prof_sfh41[cnt_sfh41,5,:] = y6new
                cnt_sfh41 = cnt_sfh41 + 1
        if T < 5 and T >= 0 and btot > 0.5 and btot < 1.0 and chi < 1.0:
                av_prof_sfh42[cnt_sfh42,0,:] = y1new
                av_prof_sfh42[cnt_sfh42,1,:] = y2new
                av_prof_sfh42[cnt_sfh42,2,:] = y3new
                av_prof_sfh42[cnt_sfh42,3,:] = y4new
                av_prof_sfh42[cnt_sfh42,4,:] = y5new
                av_prof_sfh42[cnt_sfh42,5,:] = y6new
                cnt_sfh42 = cnt_sfh42 + 1
	if T < 9 and T >=5:
                av_prof_sfh5[cnt_sfh5,0,:] = y1new
                av_prof_sfh5[cnt_sfh5,1,:] = y2new
                av_prof_sfh5[cnt_sfh5,2,:] = y3new
                av_prof_sfh5[cnt_sfh5,3,:] = y4new
                av_prof_sfh5[cnt_sfh5,4,:] = y5new
                av_prof_sfh5[cnt_sfh5,5,:] = y6new
                cnt_sfh5 = cnt_sfh5 + 1
        if T < 9 and T >=5 and btot > 0.0 and btot < 0.5 and chi < 1.0:
                av_prof_sfh51[cnt_sfh51,0,:] = y1new
                av_prof_sfh51[cnt_sfh51,1,:] = y2new
                av_prof_sfh51[cnt_sfh51,2,:] = y3new
                av_prof_sfh51[cnt_sfh51,3,:] = y4new
                av_prof_sfh51[cnt_sfh51,4,:] = y5new
                av_prof_sfh51[cnt_sfh51,5,:] = y6new
                cnt_sfh51 = cnt_sfh51 + 1
        if T < 9 and T >=5 and btot > 0.5 and btot < 1.0 and chi < 1.0:
                av_prof_sfh52[cnt_sfh52,0,:] = y1new
                av_prof_sfh52[cnt_sfh52,1,:] = y2new
                av_prof_sfh52[cnt_sfh52,2,:] = y3new
                av_prof_sfh52[cnt_sfh52,3,:] = y4new
                av_prof_sfh52[cnt_sfh52,4,:] = y5new
                av_prof_sfh52[cnt_sfh52,5,:] = y6new
                cnt_sfh52 = cnt_sfh52 + 1


fig = plt.figure(figsize=(10,12))

h1 = fig.add_axes([0.1,0.65,0.275,0.275])
h2 = fig.add_axes([0.375,0.65,0.275,0.275])
h3 = fig.add_axes([0.65,0.65,0.275,0.275])
h21 = fig.add_axes([0.375,0.375,0.275,0.275])
h31 = fig.add_axes([0.65,0.375,0.275,0.275])
h22 = fig.add_axes([0.375,0.1,0.275,0.275])
h32 = fig.add_axes([0.65,0.1,0.275,0.275])

alcalline=0.3	# alpha values for data from CALIFA
alcalfill=0.1
h1.plot(np.arange(0.,3.,0.1),rosa_AV[0][3:],color='black',alpha=alcalline)
h1.fill_between(np.arange(0.,3.,0.1),rosa_AV[0][3:]-rosa_AV[0][2],rosa_AV[0][3:]+rosa_AV[0][2],alpha=alcalfill,color='black')
h1.plot(np.arange(0.,3.,0.1),rosa_AV[2][3:],color='black',alpha=alcalline,marker='.')
h1.fill_between(np.arange(0.,3.,0.1),rosa_AV[2][3:]-rosa_AV[2][2],rosa_AV[2][3:]+rosa_AV[2][2],alpha=alcalfill,color='black')

h2.plot(np.arange(0.,3.,0.1),rosa_AV[4][3:],color='black',alpha=alcalline)
h2.fill_between(np.arange(0.,3.,0.1),rosa_AV[4][3:]-rosa_AV[4][2],rosa_AV[4][3:]+rosa_AV[4][2],alpha=alcalfill,color='black')
h2.plot(np.arange(0.,3.,0.1),rosa_AV[6][3:],color='black',alpha=alcalline,marker='.')
h2.fill_between(np.arange(0.,3.,0.1),rosa_AV[6][3:]-rosa_AV[6][2],rosa_AV[6][3:]+rosa_AV[6][2],alpha=alcalfill,color='black')
h2.plot(np.arange(0.,3.,0.1),rosa_AV[8][3:],color='black',alpha=alcalline,marker='+')
h2.fill_between(np.arange(0.,3.,0.1),rosa_AV[8][3:]-rosa_AV[8][2],rosa_AV[8][3:]+rosa_AV[8][2],alpha=alcalfill,color='black')

h21.plot(np.arange(0.,3.,0.1),rosa_AV[4][3:],color='black',alpha=alcalline)
h21.fill_between(np.arange(0.,3.,0.1),rosa_AV[4][3:]-rosa_AV[4][2],rosa_AV[4][3:]+rosa_AV[4][2],alpha=alcalfill,color='black')
h21.plot(np.arange(0.,3.,0.1),rosa_AV[6][3:],color='black',alpha=alcalline,marker='.')
h21.fill_between(np.arange(0.,3.,0.1),rosa_AV[6][3:]-rosa_AV[6][2],rosa_AV[6][3:]+rosa_AV[6][2],alpha=alcalfill,color='black')
h21.plot(np.arange(0.,3.,0.1),rosa_AV[8][3:],color='black',alpha=alcalline,marker='+')
h21.fill_between(np.arange(0.,3.,0.1),rosa_AV[8][3:]-rosa_AV[8][2],rosa_AV[8][3:]+rosa_AV[8][2],alpha=alcalfill,color='black')

h22.plot(np.arange(0.,3.,0.1),rosa_AV[4][3:],color='black',alpha=alcalline)
h22.fill_between(np.arange(0.,3.,0.1),rosa_AV[4][3:]-rosa_AV[4][2],rosa_AV[4][3:]+rosa_AV[4][2],alpha=alcalfill,color='black')
h22.plot(np.arange(0.,3.,0.1),rosa_AV[6][3:],color='black',alpha=alcalline,marker='.')
h22.fill_between(np.arange(0.,3.,0.1),rosa_AV[6][3:]-rosa_AV[6][2],rosa_AV[6][3:]+rosa_AV[6][2],alpha=alcalfill,color='black')
h22.plot(np.arange(0.,3.,0.1),rosa_AV[8][3:],color='black',alpha=alcalline,marker='+')
h22.fill_between(np.arange(0.,3.,0.1),rosa_AV[8][3:]-rosa_AV[8][2],rosa_AV[8][3:]+rosa_AV[8][2],alpha=alcalfill,color='black')

h3.plot(np.arange(0.,3.,0.1),rosa_AV[10][3:],color='black',alpha=alcalline)
h3.fill_between(np.arange(0.,3.,0.1),rosa_AV[10][3:]-rosa_AV[10][2],rosa_AV[10][3:]+rosa_AV[10][2],alpha=alcalfill,color='black')
h3.plot(np.arange(0.,3.,0.1),rosa_AV[12][3:],color='black',alpha=alcalline,marker='.')
h3.fill_between(np.arange(0.,3.,0.1),rosa_AV[12][3:]-rosa_AV[12][2],rosa_AV[12][3:]+rosa_AV[12][2],alpha=alcalfill,color='black')

h31.plot(np.arange(0.,3.,0.1),rosa_AV[10][3:],color='black',alpha=alcalline)
h31.fill_between(np.arange(0.,3.,0.1),rosa_AV[10][3:]-rosa_AV[10][2],rosa_AV[10][3:]+rosa_AV[10][2],alpha=alcalfill,color='black')
h31.plot(np.arange(0.,3.,0.1),rosa_AV[12][3:],color='black',alpha=alcalline,marker='.')
h31.fill_between(np.arange(0.,3.,0.1),rosa_AV[12][3:]-rosa_AV[12][2],rosa_AV[12][3:]+rosa_AV[12][2],alpha=alcalfill,color='black')

h32.plot(np.arange(0.,3.,0.1),rosa_AV[10][3:],color='black',alpha=alcalline)
h32.fill_between(np.arange(0.,3.,0.1),rosa_AV[10][3:]-rosa_AV[10][2],rosa_AV[10][3:]+rosa_AV[10][2],alpha=alcalfill,color='black')
h32.plot(np.arange(0.,3.,0.1),rosa_AV[12][3:],color='black',alpha=alcalline,marker='.')
h32.fill_between(np.arange(0.,3.,0.1),rosa_AV[12][3:]-rosa_AV[12][2],rosa_AV[12][3:]+rosa_AV[12][2],alpha=alcalfill,color='black')


alsfhline=1.0	# alpha values for model Kim+17
alsfhfill=0.1
linest=[':','-.','--']
for i in range(3,6):	# for each metallicities
	h1.plot(np.arange(0.,5.,0.1),np.mean(av_prof_sfh3[:,i,:],axis=0),alpha=alsfhline,color=cols_for_Z[i],linestyle=linest[i-3])	
	h1.fill_between(np.arange(0.,5.,0.1),np.mean(av_prof_sfh3[:,i,:],axis=0)-np.std(av_prof_sfh3[:,i,:],axis=0), \
		np.mean(av_prof_sfh3[:,i,:],axis=0)+np.std(av_prof_sfh3[:,i,:],axis=0),alpha=alsfhfill,color=cols_for_Z[i])

	h2.plot(np.arange(0.,5.,0.1),np.mean(av_prof_sfh4[:cnt_sfh4,i,:],axis=0),alpha=alsfhline,color=cols_for_Z[i],linestyle=linest[i-3])	
	h2.fill_between(np.arange(0.,5.,0.1),np.mean(av_prof_sfh4[:cnt_sfh4,i,:],axis=0)-np.std(av_prof_sfh4[:cnt_sfh4,i,:],axis=0), \
		np.mean(av_prof_sfh4[:cnt_sfh4,i,:],axis=0)+np.std(av_prof_sfh4[:cnt_sfh4,i,:],axis=0),alpha=alsfhfill,color=cols_for_Z[i])

        h21.plot(np.arange(0.,5.,0.1),np.mean(av_prof_sfh41[:cnt_sfh41,i,:],axis=0),alpha=alsfhline,color=cols_for_Z[i],linestyle=linest[i-3])
        h21.fill_between(np.arange(0.,5.,0.1),np.mean(av_prof_sfh41[:cnt_sfh41,i,:],axis=0)-np.std(av_prof_sfh41[:cnt_sfh41,i,:],axis=0), \
                np.mean(av_prof_sfh41[:cnt_sfh41,i,:],axis=0)+np.std(av_prof_sfh41[:cnt_sfh41,i,:],axis=0),alpha=alsfhfill,color=cols_for_Z[i])

        h22.plot(np.arange(0.,5.,0.1),np.mean(av_prof_sfh42[:cnt_sfh42,i,:],axis=0),alpha=alsfhline,color=cols_for_Z[i],linestyle=linest[i-3])
        h22.fill_between(np.arange(0.,5.,0.1),np.mean(av_prof_sfh42[:cnt_sfh42,i,:],axis=0)-np.std(av_prof_sfh42[:cnt_sfh42,i,:],axis=0), \
                np.mean(av_prof_sfh42[:cnt_sfh42,i,:],axis=0)+np.std(av_prof_sfh42[:cnt_sfh42,i,:],axis=0),alpha=alsfhfill,color=cols_for_Z[i])

	h3.plot(np.arange(0.,5.,0.1),np.mean(av_prof_sfh5[:cnt_sfh5,i,:],axis=0),alpha=alsfhline,color=cols_for_Z[i],linestyle=linest[i-3])	
	h3.fill_between(np.arange(0.,5.,0.1),np.mean(av_prof_sfh5[:cnt_sfh5,i,:],axis=0)-np.std(av_prof_sfh5[:cnt_sfh5,i,:],axis=0), \
		np.mean(av_prof_sfh5[:cnt_sfh5,i,:],axis=0)+np.std(av_prof_sfh5[:cnt_sfh5,i,:],axis=0),alpha=alsfhfill,color=cols_for_Z[i])

        h31.plot(np.arange(0.,5.,0.1),np.mean(av_prof_sfh51[:cnt_sfh51,i,:],axis=0),alpha=alsfhline,color=cols_for_Z[i],linestyle=linest[i-3])
        h31.fill_between(np.arange(0.,5.,0.1),np.mean(av_prof_sfh51[:cnt_sfh51,i,:],axis=0)-np.std(av_prof_sfh51[:cnt_sfh51,i,:],axis=0), \
                np.mean(av_prof_sfh51[:cnt_sfh51,i,:],axis=0)+np.std(av_prof_sfh51[:cnt_sfh51,i,:],axis=0),alpha=alsfhfill,color=cols_for_Z[i])

        h32.plot(np.arange(0.,5.,0.1),np.mean(av_prof_sfh52[:cnt_sfh52,i,:],axis=0),alpha=alsfhline,color=cols_for_Z[i],linestyle=linest[i-3])
        h32.fill_between(np.arange(0.,5.,0.1),np.mean(av_prof_sfh52[:cnt_sfh52,i,:],axis=0)-np.std(av_prof_sfh52[:cnt_sfh52,i,:],axis=0), \
                np.mean(av_prof_sfh52[:cnt_sfh52,i,:],axis=0)+np.std(av_prof_sfh52[:cnt_sfh52,i,:],axis=0),alpha=alsfhfill,color=cols_for_Z[i])



custom_lines = [        Line2D([0],[0],color=cols_for_Z[0],alpha=alsfhline),
                        Line2D([0],[0],color=cols_for_Z[1],alpha=alsfhline),
                        Line2D([0],[0],color=cols_for_Z[2],alpha=alsfhline),
                        Line2D([0],[0],color=cols_for_Z[3],alpha=alsfhline,linestyle=':'),
                        Line2D([0],[0],color=cols_for_Z[4],alpha=alsfhline,linestyle='-.'),
                        Line2D([0],[0],color=cols_for_Z[5],alpha=alsfhline,linestyle='--')   ]

h1.set_ylim([0,1.5])
h2.set_ylim([0,1.5])
h3.set_ylim([0,1.5])
h21.set_ylim([0,1.5])
h22.set_ylim([0,1.5])
h31.set_ylim([0,1.5])
h32.set_ylim([0,1.5])
h1.set_xlim([0,3])
h2.set_xlim([0,3])
h3.set_xlim([0,3])
h21.set_xlim([0,3])
h22.set_xlim([0,3])
h31.set_xlim([0,3])
h32.set_xlim([0,3])
h1.text(0.6,1.35,'E, S0 (SFH3)',size=20)
h2.text(0.11,1.35,'Sa, Sb, Sbc (SFH4)',size=20)
h3.text(0.5,1.35,'Sc, Sd (SFH5)',size=20)
h21.text(1.5,1.35,'B/T < 0.5',size=20,color='blue')
h31.text(1.5,1.35,'B/T < 0.5',size=20,color='blue')
h22.text(1.5,1.35,'B/T > 0.5',size=20,color='red')
h32.text(1.5,1.35,'B/T > 0.5',size=20,color='red')
h1.text(1.9,0.9,'#=%d' % cnt_sfh3,size=15)
h2.text(1.9,0.9,'#=%d' % cnt_sfh4,size=15)
h3.text(1.9,0.9,'#=%d' % cnt_sfh5,size=15)
h21.text(1.9,1.1,'#=%d' % cnt_sfh41,size=15)
h22.text(1.9,1.1,'#=%d' % cnt_sfh42,size=15)
h31.text(1.9,1.1,'#=%d' % cnt_sfh51,size=15)
h32.text(1.9,1.1,'#=%d' % cnt_sfh52,size=15)
h1.set_ylabel(r'$A_{V}$',size=20)
h21.set_ylabel(r'$A_{V}$',size=20)
h22.set_ylabel(r'$A_{V}$',size=20)
h22.set_xlabel(r'$r/r_{50,V}$',size=20)
h1.set_xlabel(r'$r/r_{50,V}$',size=20)
h32.set_xlabel(r'$r/r_{50,V}$',size=20)
h1.set_xticks([0,1,2])
h22.set_xticks([0,1,2,3])
h32.set_xticks([1,2,3])
h1.tick_params(labelsize='large',direction='in',top=True,right=True)
h2.tick_params(direction='in',top=True,right=True,bottom=True,left=True)
h3.tick_params(direction='in',top=True,right=True,bottom=True,left=True)
h21.tick_params(labelsize='large',direction='in',top=True,right=True,bottom=True,left=True)
h22.tick_params(labelsize='large',direction='in',top=True,right=True,bottom=True,left=True)
h31.tick_params(labelsize='large',direction='in',top=True,right=True,bottom=True,left=True)
h32.tick_params(labelsize='large',direction='in',top=True,right=True,bottom=True,left=True)

h1.legend(custom_lines[3:],text_for_Z[3:],loc='upper left',bbox_to_anchor=(0.51,0.905),frameon=False)
h1.plot([0.55,0.85],[1.3,1.3],color='k',alpha=alcalline)
h1.plot([1.05,1.35],[1.3,1.3],color='k',marker='.',alpha=alcalline)
h2.legend(custom_lines[3:],text_for_Z[3:],loc='upper left',bbox_to_anchor=(0.57,0.905),frameon=False)
h2.plot([0.15,0.45],[1.3,1.3],color='k',alpha=alcalline)
h2.plot([0.7,1.0],[1.3,1.3],color='k',marker='.',alpha=alcalline)
h2.plot([1.35,1.65],[1.3,1.3],color='k',marker='+',alpha=alcalline)
h3.legend(custom_lines[3:],text_for_Z[3:],loc='upper left',bbox_to_anchor=(0.51,0.905),frameon=False)
h3.plot([0.52,0.82],[1.3,1.3],color='k',alpha=alcalline)
h3.plot([1.1,1.4],[1.3,1.3],color='k',marker='.',alpha=alcalline)

h2.set_xticklabels([])
h2.set_yticklabels([])
h3.set_xticklabels([])
h3.set_yticklabels([])
h21.set_xticklabels([])
h31.set_xticklabels([])
h31.set_yticklabels([])
h32.set_yticklabels([])


fig.savefig('/Users/dhk/Documents/publish/ngcic/figs/dust_prof_all.pdf')
plt.close(fig)

