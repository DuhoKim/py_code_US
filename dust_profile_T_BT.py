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

text_for_bv0s=['-0.05', '-0.04', '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03', '0.04']
cols_for_bv0s=['black', 'gray', 'silver', 'firebrick', 'sienna', 'sandybrown', 'tan', 'olivedrab', 'palegreen', 'navy']
cols_for_AV = 'maroon'
linest_for_CALIFA = ':'

sha_cat=ascii.read(cat_dir+'sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv')	# Sample NGC/IC numbers
sha2fig=np.load('/Users/dhk/work/data/NGC_IC/pdf/sha2fig.npy')		# load figure numbers for each galaxies

rc3=ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
bt=ascii.read("/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask/btot_3rd.txt")
rosa_AV=np.genfromtxt("/Users/dhk/work/data/rosa/AV.txt")		# read A_V values for different Hubble types from CALIFA paper

file=r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'			# read NGC/IC catalogue 
df = pd.read_excel(file)
file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'				# read T-type
ttype = pd.read_excel(file)

# hlr_file = np.load('/Users/dhk/Documents/publish/ngcic/hlr_V_comb.npz')
# hlr_name = hlr_file['name']
# hlr_pix = hlr_file['hlr']

hlr_pix = np.load('/Users/dhk/work/py/pyraf/hlr_V_petro_circ.npy')

av_prof_E1 = np.zeros((31,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_E1[:] = np.nan
av_prof_S01 = np.zeros((20,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_S01[:] = np.nan
av_prof_Sa1 = np.zeros((17,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_Sa1[:] = np.nan
av_prof_Sb1 = np.zeros((40,10,30))        # For stacking Dust profiles for SFH4 w/ 6 Zs
av_prof_Sb1[:] = np.nan

av_prof_E2 = np.zeros((31,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_E2[:] = np.nan
av_prof_S02 = np.zeros((18,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_S02[:] = np.nan
av_prof_Sa2 = np.zeros((14,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_Sa2[:] = np.nan
av_prof_Sb2 = np.zeros((40,10,30))        # For stacking Dust profiles for SFH4 w/ 6 Zs
av_prof_Sb2[:] = np.nan

cnt_E1 = 0	# counter for E type and B/T > 0.6
cnt_S01 = 0	# counter for S0 type and B/T > 0.5
cnt_Sa1 = 0	# counter for Sa type and B/T > 0.45
cnt_Sb1 = 0	# counter for Sb type and B/T > 0.35
cnt_E2 = 0	# counter for E type and B/T < 0.6
cnt_S02 = 0	# counter for S0 type and B/T < 0.5
cnt_Sa2 = 0	# counter for Sa type and B/T < 0.45
cnt_Sb2 = 0	# counter for Sb type and B/T < 0.35


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

    ###### T-type match ####
    ttype_match = ttype.loc[ttype['name']==table_name]
    T = ttype_match.iloc[0,5]

    ###### B/T match #######
    ind = np.where(bt['col1']==name)
    if ind[0] >= 0:
        btot = bt['col2'][ind[0]]
        chi = 0.5
    else:
        print('no match for btot')
    #btot = ttype_match.iloc[0,7]
    #chi = ttype_match.iloc[0,8]

    if T >=9:
        continue

    try:
        x1, m1, y1, s1 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_-5.txt')
    except IOError:
        continue

    x2, m2, y2, s2 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_-4.txt')
    x3, m3, y3, s3 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_-3.txt')
    x4, m4, y4, s4 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_-2.txt')
    x5, m5, y5, s5 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_-1.txt')
    x6, m6, y6, s6 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_0.txt')
    x7, m7, y7, s7 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_1.txt')
    x8, m8, y8, s8 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_2.txt')
    x9, m9, y9, s9 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_3.txt')
    x10, m10, y10, s10 = np.loadtxt(work_dir+'ellipse/scratch_grid2/'+name+'_AV_4.txt')

    xnew = np.linspace(0, 3, num=30)
    # x_hlr = x1 / hlr_pix[np.where(hlr_name == name)]  # SMA
    x1_hlr = x1 / hlr_pix[x]
    x2_hlr = x2 / hlr_pix[x]
    x3_hlr = x3 / hlr_pix[x]
    x4_hlr = x4 / hlr_pix[x]
    x5_hlr = x5 / hlr_pix[x]
    x6_hlr = x6 / hlr_pix[x]
    x7_hlr = x7 / hlr_pix[x]
    x8_hlr = x8 / hlr_pix[x]
    x9_hlr = x9 / hlr_pix[x]
    x10_hlr = x10 / hlr_pix[x]

    y1[np.where(y1 < 0)] = 0.0
    y2[np.where(y2 < 0)] = 0.0
    y3[np.where(y3 < 0)] = 0.0
    y4[np.where(y4 < 0)] = 0.0
    y5[np.where(y5 < 0)] = 0.0
    y6[np.where(y6 < 0)] = 0.0
    y7[np.where(y7 < 0)] = 0.0
    y8[np.where(y8 < 0)] = 0.0
    y9[np.where(y9 < 0)] = 0.0
    y10[np.where(y10 < 0)] = 0.0

    f1 = interp1d(x1_hlr, y1, fill_value=np.nan, bounds_error=False)
    f2 = interp1d(x2_hlr, y2, fill_value=np.nan, bounds_error=False)
    f3 = interp1d(x3_hlr, y3, fill_value=np.nan, bounds_error=False)
    f4 = interp1d(x4_hlr, y4, fill_value=np.nan, bounds_error=False)
    f5 = interp1d(x5_hlr, y5, fill_value=np.nan, bounds_error=False)
    f6 = interp1d(x6_hlr, y6, fill_value=np.nan, bounds_error=False)
    f7 = interp1d(x7_hlr, y7, fill_value=np.nan, bounds_error=False)
    f8 = interp1d(x8_hlr, y8, fill_value=np.nan, bounds_error=False)
    f9 = interp1d(x9_hlr, y9, fill_value=np.nan, bounds_error=False)
    f10 = interp1d(x10_hlr, y10, fill_value=np.nan, bounds_error=False)

    if T < -3 and btot > 0.7:
        av_prof_E1[cnt_E1, 0, :] = f1(xnew)
        av_prof_E1[cnt_E1, 1, :] = f2(xnew)
        av_prof_E1[cnt_E1, 2, :] = f3(xnew)
        av_prof_E1[cnt_E1, 3, :] = f4(xnew)
        av_prof_E1[cnt_E1, 4, :] = f5(xnew)
        av_prof_E1[cnt_E1, 5, :] = f6(xnew)
        av_prof_E1[cnt_E1, 6, :] = f7(xnew)
        av_prof_E1[cnt_E1, 7, :] = f8(xnew)
        av_prof_E1[cnt_E1, 8, :] = f9(xnew)
        av_prof_E1[cnt_E1, 9, :] = f10(xnew)
        cnt_E1 = cnt_E1 + 1
    if T < -3 and btot < 0.7:
        av_prof_E2[cnt_E2, 0, :] = f1(xnew)
        av_prof_E2[cnt_E2, 1, :] = f2(xnew)
        av_prof_E2[cnt_E2, 2, :] = f3(xnew)
        av_prof_E2[cnt_E2, 3, :] = f4(xnew)
        av_prof_E2[cnt_E2, 4, :] = f5(xnew)
        av_prof_E2[cnt_E2, 5, :] = f6(xnew)
        av_prof_E2[cnt_E2, 6, :] = f7(xnew)
        av_prof_E2[cnt_E2, 7, :] = f8(xnew)
        av_prof_E2[cnt_E2, 8, :] = f9(xnew)
        av_prof_E2[cnt_E2, 9, :] = f10(xnew)
        cnt_E2 = cnt_E2 + 1
    if T < 0 and T >= -3 and btot > 0.5:
        av_prof_S01[cnt_S01, 0, :] = f1(xnew)
        av_prof_S01[cnt_S01, 1, :] = f2(xnew)
        av_prof_S01[cnt_S01, 2, :] = f3(xnew)
        av_prof_S01[cnt_S01, 3, :] = f4(xnew)
        av_prof_S01[cnt_S01, 4, :] = f5(xnew)
        av_prof_S01[cnt_S01, 5, :] = f6(xnew)
        av_prof_S01[cnt_S01, 6, :] = f7(xnew)
        av_prof_S01[cnt_S01, 7, :] = f8(xnew)
        av_prof_S01[cnt_S01, 8, :] = f9(xnew)
        av_prof_S01[cnt_S01, 9, :] = f10(xnew)
        cnt_S01 = cnt_S01 + 1
    if T < 0 and T >= -3 and btot < 0.5:
        av_prof_S02[cnt_S02, 0, :] = f1(xnew)
        av_prof_S02[cnt_S02, 1, :] = f2(xnew)
        av_prof_S02[cnt_S02, 2, :] = f3(xnew)
        av_prof_S02[cnt_S02, 3, :] = f4(xnew)
        av_prof_S02[cnt_S02, 4, :] = f5(xnew)
        av_prof_S02[cnt_S02, 5, :] = f6(xnew)
        av_prof_S02[cnt_S02, 6, :] = f7(xnew)
        av_prof_S02[cnt_S02, 7, :] = f8(xnew)
        av_prof_S02[cnt_S02, 8, :] = f9(xnew)
        av_prof_S02[cnt_S02, 9, :] = f10(xnew)
        cnt_S02 = cnt_S02 + 1
    if T < 2 and T >= 0 and btot > 0.45:
        av_prof_Sa1[cnt_Sa1, 0, :] = f1(xnew)
        av_prof_Sa1[cnt_Sa1, 1, :] = f2(xnew)
        av_prof_Sa1[cnt_Sa1, 2, :] = f3(xnew)
        av_prof_Sa1[cnt_Sa1, 3, :] = f4(xnew)
        av_prof_Sa1[cnt_Sa1, 4, :] = f5(xnew)
        av_prof_Sa1[cnt_Sa1, 5, :] = f6(xnew)
        av_prof_Sa1[cnt_Sa1, 6, :] = f7(xnew)
        av_prof_Sa1[cnt_Sa1, 7, :] = f8(xnew)
        av_prof_Sa1[cnt_Sa1, 8, :] = f9(xnew)
        av_prof_Sa1[cnt_Sa1, 9, :] = f10(xnew)
        cnt_Sa1 = cnt_Sa1 + 1
    if T < 2 and T >= 0 and btot < 0.45:
        av_prof_Sa2[cnt_Sa2, 0, :] = f1(xnew)
        av_prof_Sa2[cnt_Sa2, 1, :] = f2(xnew)
        av_prof_Sa2[cnt_Sa2, 2, :] = f3(xnew)
        av_prof_Sa2[cnt_Sa2, 3, :] = f4(xnew)
        av_prof_Sa2[cnt_Sa2, 4, :] = f5(xnew)
        av_prof_Sa2[cnt_Sa2, 5, :] = f6(xnew)
        av_prof_Sa2[cnt_Sa2, 6, :] = f7(xnew)
        av_prof_Sa2[cnt_Sa2, 7, :] = f8(xnew)
        av_prof_Sa2[cnt_Sa2, 8, :] = f9(xnew)
        av_prof_Sa2[cnt_Sa2, 9, :] = f10(xnew)
        cnt_Sa2 = cnt_Sa2 + 1
    if T < 4 and T >= 2 and btot > 0.25:
        av_prof_Sb1[cnt_Sb1, 0, :] = f1(xnew)
        av_prof_Sb1[cnt_Sb1, 1, :] = f2(xnew)
        av_prof_Sb1[cnt_Sb1, 2, :] = f3(xnew)
        av_prof_Sb1[cnt_Sb1, 3, :] = f4(xnew)
        av_prof_Sb1[cnt_Sb1, 4, :] = f5(xnew)
        av_prof_Sb1[cnt_Sb1, 5, :] = f6(xnew)
        av_prof_Sb1[cnt_Sb1, 6, :] = f7(xnew)
        av_prof_Sb1[cnt_Sb1, 7, :] = f8(xnew)
        av_prof_Sb1[cnt_Sb1, 8, :] = f9(xnew)
        av_prof_Sb1[cnt_Sb1, 9, :] = f10(xnew)
        cnt_Sb1 = cnt_Sb1 + 1
    if T < 4 and T >= 2 and btot < 0.25:
        av_prof_Sb2[cnt_Sb2, 0, :] = f1(xnew)
        av_prof_Sb2[cnt_Sb2, 1, :] = f2(xnew)
        av_prof_Sb2[cnt_Sb2, 2, :] = f3(xnew)
        av_prof_Sb2[cnt_Sb2, 3, :] = f4(xnew)
        av_prof_Sb2[cnt_Sb2, 4, :] = f5(xnew)
        av_prof_Sb2[cnt_Sb2, 5, :] = f6(xnew)
        av_prof_Sb2[cnt_Sb2, 6, :] = f7(xnew)
        av_prof_Sb2[cnt_Sb2, 7, :] = f8(xnew)
        av_prof_Sb2[cnt_Sb2, 8, :] = f9(xnew)
        av_prof_Sb2[cnt_Sb2, 9, :] = f10(xnew)
        cnt_Sb2 = cnt_Sb2 + 1


fig = plt.figure(figsize=(10,7))

h11 = fig.add_axes([0.1,0.525,0.2,0.375])
h21 = fig.add_axes([0.3,0.525,0.2,0.375])
h31 = fig.add_axes([0.5,0.525,0.2,0.375])
h41 = fig.add_axes([0.7,0.525,0.2,0.375])
h12 = fig.add_axes([0.1,0.15,0.2,0.375])
h22 = fig.add_axes([0.3,0.15,0.2,0.375])
h32 = fig.add_axes([0.5,0.15,0.2,0.375])
h42 = fig.add_axes([0.7,0.15,0.2,0.375])

alcalline=0.3	# alpha values for data from CALIFA
alcalfill=0.1
h11.plot(np.arange(0.,3.,0.1),rosa_AV[0][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h11.fill_between(np.arange(0.,3.,0.1),rosa_AV[0][3:]-rosa_AV[0][2],rosa_AV[0][3:]+rosa_AV[0][2],alpha=alcalfill,color='black')
h11.errorbar(1, rosa_AV[0][12], yerr=rosa_AV[0][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h21.plot(np.arange(0.,3.,0.1),rosa_AV[2][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h21.fill_between(np.arange(0.,3.,0.1),rosa_AV[2][3:]-rosa_AV[2][2],rosa_AV[2][3:]+rosa_AV[2][2],alpha=alcalfill,color='black')
h21.errorbar(1, rosa_AV[2][12], yerr=rosa_AV[2][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h31.plot(np.arange(0.,3.,0.1),rosa_AV[4][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h31.fill_between(np.arange(0.,3.,0.1),rosa_AV[4][3:]-rosa_AV[4][2],rosa_AV[4][3:]+rosa_AV[4][2],alpha=alcalfill,color='black')
h31.errorbar(1, rosa_AV[4][12], yerr=rosa_AV[4][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h41.plot(np.arange(0.,3.,0.1),rosa_AV[6][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h41.fill_between(np.arange(0.,3.,0.1),rosa_AV[6][3:]-rosa_AV[6][2],rosa_AV[6][3:]+rosa_AV[6][2],alpha=alcalfill,color='black')
h41.errorbar(1, rosa_AV[6][12], yerr=rosa_AV[6][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h12.plot(np.arange(0.,3.,0.1),rosa_AV[0][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h12.fill_between(np.arange(0.,3.,0.1),rosa_AV[0][3:]-rosa_AV[0][2],rosa_AV[0][3:]+rosa_AV[0][2],alpha=alcalfill,color='black')
h12.errorbar(1, rosa_AV[0][12], yerr=rosa_AV[0][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h22.plot(np.arange(0.,3.,0.1),rosa_AV[2][3:],color='black',alpha=alcalline, linewidth=3, label='GD15 (CALIFA)', linestyle=linest_for_CALIFA)
#h22.fill_between(np.arange(0.,3.,0.1),rosa_AV[2][3:]-rosa_AV[2][2],rosa_AV[2][3:]+rosa_AV[2][2],alpha=alcalfill,color='black')
h22.errorbar(1, rosa_AV[2][12], yerr=rosa_AV[2][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h32.plot(np.arange(0.,3.,0.1),rosa_AV[4][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h32.fill_between(np.arange(0.,3.,0.1),rosa_AV[4][3:]-rosa_AV[4][2],rosa_AV[4][3:]+rosa_AV[4][2],alpha=alcalfill,color='black')
h32.errorbar(1, rosa_AV[4][12], yerr=rosa_AV[4][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h42.plot(np.arange(0.,3.,0.1),rosa_AV[6][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h42.fill_between(np.arange(0.,3.,0.1),rosa_AV[6][3:]-rosa_AV[6][2],rosa_AV[6][3:]+rosa_AV[6][2],alpha=alcalfill,color='black')
h42.errorbar(1, rosa_AV[6][12], yerr=rosa_AV[6][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

min_idx = np.zeros(8)
min_chi = np.zeros(8)
min_chi[:] = 1e6
for i in range(0, 10):
    prof_E1 = np.nanmean(av_prof_E1[:, i, :], axis=0)
    prof_S01 = np.nanmean(av_prof_S01[:, i, :], axis=0)
    prof_Sa1 = np.nanmean(av_prof_Sa1[:, i, :], axis=0)
    prof_Sb1 = np.nanmean(av_prof_Sb1[:, i, :], axis=0)
    prof_E2 = np.nanmean(av_prof_E2[:, i, :], axis=0)
    prof_S02 = np.nanmean(av_prof_S02[:, i, :], axis=0)
    prof_Sa2 = np.nanmean(av_prof_Sa2[:, i, :], axis=0)
    prof_Sb2 = np.nanmean(av_prof_Sb2[:, i, :], axis=0)

    prof_E1_err = np.nanstd(av_prof_E1[:, i, :], axis=0)
    prof_S01_err = np.nanstd(av_prof_S01[:, i, :], axis=0)
    prof_Sa1_err = np.nanstd(av_prof_Sa1[:, i, :], axis=0)
    prof_Sb1_err = np.nanstd(av_prof_Sb1[:, i, :], axis=0)
    prof_E2_err = np.nanstd(av_prof_E2[:, i, :], axis=0)
    prof_S02_err = np.nanstd(av_prof_S02[:, i, :], axis=0)
    prof_Sa2_err = np.nanstd(av_prof_Sa2[:, i, :], axis=0)
    prof_Sb2_err = np.nanstd(av_prof_Sb2[:, i, :], axis=0)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_E1[j]) and ~np.isnan(prof_E1_err[j]):
            chi += np.abs(prof_E1[j] - rosa_AV[0][3+j]) / prof_E1_err[j]
    if chi < min_chi[0]:
        min_chi[0] = chi
        min_idx[0] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_S01[j]) and ~np.isnan(prof_S01_err[j]):
            chi += np.abs(prof_S01[j] - rosa_AV[2][3 + j]) / prof_S01_err[j]
    if chi < min_chi[1]:
        min_chi[1] = chi
        min_idx[1] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sa1[j]) and ~np.isnan(prof_Sa1_err[j]):
            chi += np.abs(prof_Sa1[j] - rosa_AV[4][3 + j]) / prof_Sa1_err[j]
    if chi < min_chi[2]:
        min_chi[2] = chi
        min_idx[2] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sb1[j]) and ~np.isnan(prof_Sb1_err[j]):
            chi += np.abs(prof_Sb1[j] - rosa_AV[6][3 + j]) / prof_Sb1_err[j]
    if chi < min_chi[3]:
        min_chi[3] = chi
        min_idx[3] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_E2[j]) and ~np.isnan(prof_E2_err[j]):
            chi += np.abs(prof_E2[j] - rosa_AV[0][3 + j]) / prof_E2_err[j]
    if chi < min_chi[4]:
        min_chi[4] = chi
        min_idx[4] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_S02[j]) and ~np.isnan(prof_S02_err[j]):
            chi += np.abs(prof_S02[j] - rosa_AV[2][3 + j]) / prof_S02_err[j]
    if chi < min_chi[5]:
        min_chi[5] = chi
        min_idx[5] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sa2[j]) and ~np.isnan(prof_Sa2_err[j]):
            chi += np.abs(prof_Sa2[j] - rosa_AV[4][3 + j]) / prof_Sa2_err[j]
    if chi < min_chi[6]:
        min_chi[6] = chi
        min_idx[6] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sb2[j]) and ~np.isnan(prof_Sb2_err[j]):
            chi += np.abs(prof_Sb2[j] - rosa_AV[6][3 + j]) / prof_Sb2_err[j]
    if chi < min_chi[7]:
        min_chi[7] = chi
        min_idx[7] = int(i)


alsfhline=0.1	# alpha values for model Kim+17
alsfhfill=0.1
linest=[':','-.','--']
#for i in range(0,10):	# for each metallicities
#    h1.plot(np.arange(0., 3., 0.1),np.nanmean(av_prof_E[:, i, :],axis=0),alpha=alsfhline, label=text_for_bv0s[i], color='k')
#    h2.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_S0[:, i, :], axis=0), alpha=alsfhline, label=text_for_bv0s[i], color='k')
#    h3.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sa[:, i, :], axis=0), alpha=alsfhline, label=text_for_bv0s[i], color='k')
#    h4.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sb[:, i, :], axis=0), alpha=alsfhline, label=text_for_bv0s[i], color='k')
#    h5.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sbc[:, i, :], axis=0), alpha=alsfhline, label=text_for_bv0s[i], color='k')
#    h6.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sc[:, i, :], axis=0), alpha=alsfhline, label=text_for_bv0s[i], color='k')
#    h7.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sd[:, i, :], axis=0), alpha=alsfhline, label=text_for_bv0s[i], color='k')

h11.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_E1[:, int(min_idx[0]), :], axis=0), color=cols_for_AV, linewidth=3)
h21.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_S01[:, int(min_idx[1]), :], axis=0), color=cols_for_AV, linewidth=3)
h31.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sa1[:, int(min_idx[2]), :], axis=0), color=cols_for_AV, linewidth=3)
h41.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sb1[:, int(min_idx[3]), :], axis=0), color=cols_for_AV, linewidth=3)
h12.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_E2[:, int(min_idx[4]), :], axis=0), color=cols_for_AV, linewidth=3)
h22.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_S02[:, int(min_idx[5]), :], axis=0), color=cols_for_AV, linewidth=3)
h32.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sa2[:, int(min_idx[6]), :], axis=0), color=cols_for_AV, linewidth=3)
h42.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sb2[:, int(min_idx[7]), :], axis=0), color=cols_for_AV, linewidth=3, label=r'$\beta_{V}$ method')


h11.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_E1[:, int(min_idx[0]), :], axis=0) - np.nanstd(av_prof_E1[:, int(min_idx[0]), :], axis=0),
    np.nanmean(av_prof_E1[:, int(min_idx[0]), :], axis=0) + np.nanstd(av_prof_E1[:, int(min_idx[0]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h21.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_S01[:, int(min_idx[1]), :], axis=0) - np.nanstd(av_prof_S01[:, int(min_idx[1]), :], axis=0),
    np.nanmean(av_prof_S01[:, int(min_idx[1]), :], axis=0) + np.nanstd(av_prof_S01[:, int(min_idx[1]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h31.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sa1[:, int(min_idx[2]), :], axis=0) - np.nanstd(av_prof_Sa1[:, int(min_idx[2]), :], axis=0),
    np.nanmean(av_prof_Sa1[:, int(min_idx[2]), :], axis=0) + np.nanstd(av_prof_Sa1[:, int(min_idx[2]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h41.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sb1[:, int(min_idx[3]), :], axis=0) - np.nanstd(av_prof_Sb1[:, int(min_idx[3]), :], axis=0),
    np.nanmean(av_prof_Sb1[:, int(min_idx[3]), :], axis=0) + np.nanstd(av_prof_Sb1[:, int(min_idx[3]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h12.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_E2[:, int(min_idx[4]), :], axis=0) - np.nanstd(av_prof_E2[:, int(min_idx[4]), :], axis=0),
    np.nanmean(av_prof_E2[:, int(min_idx[4]), :], axis=0) + np.nanstd(av_prof_E2[:, int(min_idx[4]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h22.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_S02[:, int(min_idx[5]), :], axis=0) - np.nanstd(av_prof_S02[:, int(min_idx[5]), :], axis=0),
    np.nanmean(av_prof_S02[:, int(min_idx[5]), :], axis=0) + np.nanstd(av_prof_S02[:, int(min_idx[5]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h32.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sa2[:, int(min_idx[6]), :], axis=0) - np.nanstd(av_prof_Sa2[:, int(min_idx[6]), :], axis=0),
    np.nanmean(av_prof_Sa2[:, int(min_idx[6]), :], axis=0) + np.nanstd(av_prof_Sa2[:, int(min_idx[6]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h42.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sb2[:, int(min_idx[7]), :], axis=0) - np.nanstd(av_prof_Sb2[:, int(min_idx[7]), :], axis=0),
    np.nanmean(av_prof_Sb2[:, int(min_idx[7]), :], axis=0) + np.nanstd(av_prof_Sb2[:, int(min_idx[7]), :], axis=0), alpha=alcalfill, color=cols_for_AV)


h11.set_ylim([0,1])
h21.set_ylim([0,1])
h31.set_ylim([0,1])
h41.set_ylim([0,1])
h12.set_ylim([0,1])
h22.set_ylim([0,1])
h32.set_ylim([0,1])
h42.set_ylim([0,1])

h11.set_xlim([0,3])
h21.set_xlim([0,3])
h31.set_xlim([0,3])
h41.set_xlim([0,3])
h12.set_xlim([0,3])
h22.set_xlim([0,3])
h32.set_xlim([0,3])
h42.set_xlim([0,3])

h11.text(1.4,1.05,'E',size=30)
h21.text(1.3,1.05,'S0',size=30)
h31.text(1.3,1.05,'Sa',size=30)
h41.text(1.3,1.05,'Sb',size=30)

h41.text(3.1,0.7,'Large bulge',size=20, rotation=270, color='red')
h42.text(3.1,0.7,'Small bulge',size=20, rotation=270, color='blue')

h11.text(0.7,0.85,'B/T > 0.7',size=20, color=cols_for_AV)
h21.text(0.7,0.85,'B/T > 0.5',size=20, color=cols_for_AV)
h31.text(0.5,0.85,'B/T > 0.45',size=20, color=cols_for_AV)
h41.text(0.5,0.85,'B/T > 0.25',size=20, color=cols_for_AV)
h12.text(0.7,0.85,'B/T < 0.7',size=20, color=cols_for_AV)
h22.text(0.7,0.85,'B/T < 0.5',size=20, color=cols_for_AV)
h32.text(0.5,.85,'B/T < 0.45',size=20, color=cols_for_AV)
h42.text(0.5,0.85,'B/T < 0.25',size=20, color=cols_for_AV)

h11.text(1.5,0.6,'N=%d' % cnt_E1,size=20, color=cols_for_AV)
h21.text(1.5,0.6,'N=%d' % cnt_S01,size=20, color=cols_for_AV)
h31.text(1.5,0.6,'N=%d' % cnt_Sa1,size=20, color=cols_for_AV)
h41.text(1.5,0.6,'N=%d' % cnt_Sb1,size=20, color=cols_for_AV)
h12.text(1.5,0.6,'N=%d' % cnt_E2,size=20, color=cols_for_AV)
h22.text(1.5,0.6,'N=%d' % cnt_S02,size=20, color=cols_for_AV)
h32.text(1.5,0.6,'N=%d' % cnt_Sa2,size=20, color=cols_for_AV)
h42.text(1.5,0.6,'N=%d' % cnt_Sb2,size=20, color=cols_for_AV)

h11.text(0.6,0.725,r'$\beta_{V,0}$=0.51',size=20, color=cols_for_AV)
h21.text(0.6,0.725,r'$\beta_{V,0}$=0.58',size=20, color=cols_for_AV)
h31.text(0.6,0.725,r'$\beta_{V,0}$=0.68',size=20, color=cols_for_AV)
h41.text(0.6,0.725,r'$\beta_{V,0}$=0.69',size=20, color=cols_for_AV)
h12.text(0.8,0.725,r'$\beta_{V,0}$=0.5',size=20, color=cols_for_AV)
h22.text(0.6,0.725,r'$\beta_{V,0}$=0.64',size=20, color=cols_for_AV)
h32.text(0.6,0.725,r'$\beta_{V,0}$=0.76',size=20, color=cols_for_AV)
h42.text(0.6,0.725,r'$\beta_{V,0}$=0.79',size=20, color=cols_for_AV)


h11.set_ylabel(r'$\langle A_{V}\rangle$',size=30)
h11.yaxis.set_label_coords(-0.2,0)
h11.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
h11.set_yticklabels([0, '', 0.5, '', 1.0])
h22.set_xlabel(r'$R/R_{50,V}^{P}$',size=30)
h22.xaxis.set_label_coords(1.0,-0.1)
h12.set_yticks([0,0.2,0.4,0.6,0.8])
h12.set_xticks([0,1,2,3])
h22.set_xticks([1,2,3])
h32.set_xticks([1,2,3])
h42.set_xticks([1,2,3])
h12.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
h12.set_yticklabels([0, '', 0.5, '', ''])
h11.tick_params(labelsize='20',direction='in',top=True,right=True)
h21.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h31.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h41.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h12.tick_params(labelsize='20',direction='in',top=True,right=True)
h22.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h32.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h42.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)

h22.legend(bbox_to_anchor=(0.6, -0.05), prop={'size': 20}, frameon=False)
h42.legend(bbox_to_anchor=(1.0, -0.05), prop={'size': 20}, frameon=False)

h11.set_xticklabels([])
h21.set_xticklabels([])
h31.set_xticklabels([])
h41.set_xticklabels([])
h21.set_yticklabels([])
h31.set_yticklabels([])
h41.set_yticklabels([])
h22.set_yticklabels([])
h32.set_yticklabels([])
h42.set_yticklabels([])

fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/dust_prof_all_T_BT.pdf')
plt.close(fig)

