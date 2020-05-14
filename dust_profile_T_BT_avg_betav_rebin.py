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

hubble_bv0 = [0.55, 0.61, 0.70, 0.73, 0.56, 0.75, 0.76, 0.81]

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

bv_prof_E1 = np.zeros((19, 15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_E1[:] = np.nan
bv_prof_S01 = np.zeros((20, 15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_S01[:] = np.nan
bv_prof_Sa1 = np.zeros((17, 15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_Sa1[:] = np.nan
bv_prof_Sb1 = np.zeros((22, 15))        # For stacking Dust profiles for SFH4 w/ 6 Zs
bv_prof_Sb1[:] = np.nan

bv_prof_E2 = np.zeros((27, 15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_E2[:] = np.nan
bv_prof_S02 = np.zeros((18, 15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_S02[:] = np.nan
bv_prof_Sa2 = np.zeros((14, 15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_Sa2[:] = np.nan
bv_prof_Sb2 = np.zeros((21, 15))        # For stacking Dust profiles for SFH4 w/ 6 Zs
bv_prof_Sb2[:] = np.nan

nan_cnt_E1 = np.zeros(15)        # for counting NaN values in the stack
nan_cnt_E2 = np.zeros(15)        # for counting NaN values in the stack
nan_cnt_S01 = np.zeros(15)        # for counting NaN values in the stack
nan_cnt_S02 = np.zeros(15)        # for counting NaN values in the stack
nan_cnt_Sa1 = np.zeros(15)        # for counting NaN values in the stack
nan_cnt_Sa2 = np.zeros(15)        # for counting NaN values in the stack
nan_cnt_Sb1 = np.zeros(15)        # for counting NaN values in the stack
nan_cnt_Sb2 = np.zeros(15)        # for counting NaN values in the stack


cnt_E1 = 0	# counter for E type and B/T > 0.6
cnt_S01 = 0	# counter for S0 type and B/T > 0.5
cnt_Sa1 = 0	# counter for Sa type and B/T > 0.45
cnt_Sb1 = 0	# counter for Sb type and B/T > 0.35
cnt_E2 = 0	# counter for E type and B/T < 0.6
cnt_S02 = 0	# counter for S0 type and B/T < 0.5
cnt_Sa2 = 0	# counter for Sa type and B/T < 0.45
cnt_Sb2 = 0	# counter for Sb type and B/T < 0.35

xnew = np.linspace(0, 3, num=15)

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
        x1, m1, y1, s1 = np.loadtxt(work_dir+'ellipse/scratch_ellip/'+name+'_bv.txt')        # rads,iso_mean,iso_med,iso_std
    except IOError:
        continue

    # x_hlr = x1 / hlr_pix[np.where(hlr_name == name)]  # SMA
    x1_hlr = x1 / hlr_pix[x]



    f1 = interp1d(x1_hlr, y1, fill_value=np.nan, bounds_error=False)

    if T < -3 and btot > 0.7:
        bv_prof_E1[cnt_E1, :] = f1(xnew)
        nan_cnt_E1 += np.isnan(f1(xnew))
        cnt_E1 = cnt_E1 + 1
    if T < -3 and btot < 0.7:
        bv_prof_E2[cnt_E2, :] = f1(xnew)
        nan_cnt_E2 += np.isnan(f1(xnew))
        cnt_E2 = cnt_E2 + 1
    if T < 0 and T >= -3 and btot > 0.5:
        bv_prof_S01[cnt_S01, :] = f1(xnew)
        nan_cnt_S01 += np.isnan(f1(xnew))
        cnt_S01 = cnt_S01 + 1
    if T < 0 and T >= -3 and btot < 0.5:
        bv_prof_S02[cnt_S02, :] = f1(xnew)
        nan_cnt_S02 += np.isnan(f1(xnew))
        cnt_S02 = cnt_S02 + 1
    if T < 2 and T >= 0 and btot > 0.45:
        bv_prof_Sa1[cnt_Sa1, :] = f1(xnew)
        nan_cnt_Sa1 += np.isnan(f1(xnew))
        cnt_Sa1 = cnt_Sa1 + 1
    if T < 2 and T >= 0 and btot < 0.45:
        bv_prof_Sa2[cnt_Sa2, :] = f1(xnew)
        nan_cnt_Sa2 += np.isnan(f1(xnew))
        cnt_Sa2 = cnt_Sa2 + 1
    if T < 4 and T >= 2 and btot > 0.25:
        bv_prof_Sb1[cnt_Sb1, :] = f1(xnew)
        nan_cnt_Sb1 += np.isnan(f1(xnew))
        cnt_Sb1 = cnt_Sb1 + 1
    if T < 4 and T >= 2 and btot < 0.25:
        bv_prof_Sb2[cnt_Sb2, :] = f1(xnew)
        nan_cnt_Sb2 += np.isnan(f1(xnew))
        cnt_Sb2 = cnt_Sb2 + 1

flag_E1 = np.where(nan_cnt_E1 <= (cnt_E1 - 5))[0]
flag_S01 = np.where(nan_cnt_S01 <= (cnt_S01 - 5))[0]
flag_Sa1 = np.where(nan_cnt_Sa1 <= (cnt_Sa1 - 5))[0]
flag_Sb1 = np.where(nan_cnt_Sb1 <= (cnt_Sb1 - 5))[0]
flag_E2 = np.where(nan_cnt_E2 <= (cnt_E2 - 5))[0]
flag_S02 = np.where(nan_cnt_S02 <= (cnt_S02 - 5))[0]
flag_Sa2 = np.where(nan_cnt_Sa2 <= (cnt_Sa2 - 5))[0]
flag_Sb2 = np.where(nan_cnt_Sb2 <= (cnt_Sb2 - 5))[0]

prof_E1 = np.nanmedian(bv_prof_E1[:, :], axis=0)
prof_S01 = np.nanmedian(bv_prof_S01[:, :], axis=0)
prof_Sa1 = np.nanmedian(bv_prof_Sa1[:, :], axis=0)
prof_Sb1 = np.nanmedian(bv_prof_Sb1[:, :], axis=0)
prof_E2 = np.nanmedian(bv_prof_E2[:, :], axis=0)
prof_S02 = np.nanmedian(bv_prof_S02[:, :], axis=0)
prof_Sa2 = np.nanmedian(bv_prof_Sa2[:, :], axis=0)
prof_Sb2 = np.nanmedian(bv_prof_Sb2[:, :], axis=0)

prof_E1_low = np.nanpercentile(bv_prof_E1[:, :],  15.87, axis=0)
prof_S01_low = np.nanpercentile(bv_prof_S01[:, :], 15.87, axis=0)
prof_Sa1_low = np.nanpercentile(bv_prof_Sa1[:, :], 15.87, axis=0)
prof_Sb1_low = np.nanpercentile(bv_prof_Sb1[:, :], 15.87, axis=0)
prof_E2_low = np.nanpercentile(bv_prof_E2[:, :],  15.87, axis=0)
prof_S02_low = np.nanpercentile(bv_prof_S02[:, :], 15.87, axis=0)
prof_Sa2_low = np.nanpercentile(bv_prof_Sa2[:, :], 15.87, axis=0)
prof_Sb2_low = np.nanpercentile(bv_prof_Sb2[:, :], 15.87, axis=0)

prof_E1_hi = np.nanpercentile(bv_prof_E1[:, :], 84.13, axis=0)
prof_S01_hi = np.nanpercentile(bv_prof_S01[:, :], 84.13, axis=0)
prof_Sa1_hi = np.nanpercentile(bv_prof_Sa1[:, :], 84.13, axis=0)
prof_Sb1_hi = np.nanpercentile(bv_prof_Sb1[:, :], 84.13, axis=0)
prof_E2_hi = np.nanpercentile(bv_prof_E2[:, :], 84.13, axis=0)
prof_S02_hi = np.nanpercentile(bv_prof_S02[:, :], 84.13, axis=0)
prof_Sa2_hi = np.nanpercentile(bv_prof_Sa2[:, :], 84.13, axis=0)
prof_Sb2_hi = np.nanpercentile(bv_prof_Sb2[:, :], 84.13, axis=0)

min_bv0 = np.zeros(8)
min_chi = np.zeros(8)
min_chi[:] = 1e6

for i in np.arange(-0.1, 0.1, 0.01):
    av_prof_E1 = 2.5 * np.log10((hubble_bv0[0] + i) / prof_E1)
    av_prof_S01 = 2.5 * np.log10((hubble_bv0[1] + i) / prof_S01)
    av_prof_Sa1 = 2.5 * np.log10((hubble_bv0[2] + i) / prof_Sa1)
    av_prof_Sb1 = 2.5 * np.log10((hubble_bv0[3] + i) / prof_Sb1)
    av_prof_E2 = 2.5 * np.log10((hubble_bv0[4] + i) / prof_E2)
    av_prof_S02 = 2.5 * np.log10((hubble_bv0[5] + i) / prof_S02)
    av_prof_Sa2 = 2.5 * np.log10((hubble_bv0[6] + i) / prof_Sa2)
    av_prof_Sb2 = 2.5 * np.log10((hubble_bv0[7] + i) / prof_Sb2)

    av_prof_E1[np.where(av_prof_E1 < 0)] = 0.0
    av_prof_S01[np.where(av_prof_S01 < 0)] = 0.0
    av_prof_Sa1[np.where(av_prof_Sa1 < 0)] = 0.0
    av_prof_Sb1[np.where(av_prof_Sb1 < 0)] = 0.0
    av_prof_E2[np.where(av_prof_E2 < 0)] = 0.0
    av_prof_S02[np.where(av_prof_S02 < 0)] = 0.0
    av_prof_Sa2[np.where(av_prof_Sa2 < 0)] = 0.0
    av_prof_Sb2[np.where(av_prof_Sb2 < 0)] = 0.0

    chi = 0
    for j in flag_E1:
        sig = 1.25 * np.log10(prof_E1_low[j] / prof_E1_hi[j])  # half of magenta range in Figure 6
        chi += (av_prof_E1[j] - rosa_AV[0][3+j*2])**2 / sig**2
    if chi < min_chi[0]:
        min_chi[0] = chi
        min_bv0[0] = hubble_bv0[0] + i

    chi = 0
    for j in flag_S01:
        sig = 1.25 * np.log10(prof_S01_low[j] / prof_S01_hi[j])  # half of magenta range in Figure 6
        chi += (av_prof_S01[j] - rosa_AV[2][3 + j*2])**2 / sig**2
    if chi < min_chi[1]:
        min_chi[1] = chi
        min_bv0[1] = hubble_bv0[1] + i

    chi = 0
    for j in flag_Sa1:
        sig = 1.25 * np.log10(prof_Sa1_low[j] / prof_Sa1_hi[j])  # half of magenta range in Figure 6
        chi += (av_prof_Sa1[j] - rosa_AV[4][3 + j*2])**2 / sig**2
    if chi < min_chi[2]:
        min_chi[2] = chi
        min_bv0[2] = hubble_bv0[2] + i

    chi = 0
    for j in flag_Sb1:
        sig = 1.25 * np.log10(prof_Sb1_low[j] / prof_Sb1_hi[j])  # half of magenta range in Figure 6
        chi += (av_prof_Sb1[j] - rosa_AV[6][3 + j*2])**2 / sig**2
    if chi < min_chi[3]:
        min_chi[3] = chi
        min_bv0[3] = hubble_bv0[3] + i

    chi = 0
    for j in flag_E2:
        sig = 1.25 * np.log10(prof_E2_low[j] / prof_E2_hi[j])  # half of magenta range in Figure 6
        chi += (av_prof_E2[j] - rosa_AV[0][3+j*2])**2 / sig**2
    if chi < min_chi[4]:
        min_chi[4] = chi
        min_bv0[4] = hubble_bv0[4] + i

    chi = 0
    for j in flag_S02:
        sig = 1.25 * np.log10(prof_S02_low[j] / prof_S02_hi[j])  # half of magenta range in Figure 6
        chi += (av_prof_S02[j] - rosa_AV[2][3 + j*2])**2 / sig**2
    if chi < min_chi[5]:
        min_chi[5] = chi
        min_bv0[5] = hubble_bv0[5] + i

    chi = 0
    for j in flag_Sa2:
        sig = 1.25 * np.log10(prof_Sa2_low[j] / prof_Sa2_hi[j])  # half of magenta range in Figure 6
        chi += (av_prof_Sa2[j] - rosa_AV[4][3 + j*2])**2 / sig**2
    if chi < min_chi[6]:
        min_chi[6] = chi
        min_bv0[6] = hubble_bv0[6] + i

    chi = 0
    for j in flag_Sb2:
        sig = 1.25 * np.log10(prof_Sb2_low[j] / prof_Sb2_hi[j])  # half of magenta range in Figure 6
        chi += (av_prof_Sb2[j] - rosa_AV[6][3 + j*2])**2 / sig**2
    if chi < min_chi[7]:
        min_chi[7] = chi
        min_bv0[7] = hubble_bv0[7] + i



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

xgrid = np.arange(0.,3.,0.1)

h11.plot(xgrid,rosa_AV[0][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h11.errorbar(1, rosa_AV[0][12], yerr=rosa_AV[0][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h21.plot(xgrid,rosa_AV[2][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h21.errorbar(1, rosa_AV[2][12], yerr=rosa_AV[2][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h31.plot(xgrid,rosa_AV[4][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h31.errorbar(1, rosa_AV[4][12], yerr=rosa_AV[4][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h41.plot(xgrid,rosa_AV[6][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h41.errorbar(1, rosa_AV[6][12], yerr=rosa_AV[6][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h12.plot(xgrid,rosa_AV[0][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h12.errorbar(1, rosa_AV[0][12], yerr=rosa_AV[0][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h22.plot(xgrid,rosa_AV[2][3:],color='black',alpha=alcalline, linewidth=3, label='GD15 (CALIFA)', linestyle=linest_for_CALIFA)
h22.errorbar(1, rosa_AV[2][12], yerr=rosa_AV[2][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h32.plot(xgrid,rosa_AV[4][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h32.errorbar(1, rosa_AV[4][12], yerr=rosa_AV[4][2], color='black', alpha=alcalline, linewidth=2, capsize=3)
h42.plot(xgrid,rosa_AV[6][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h42.errorbar(1, rosa_AV[6][12], yerr=rosa_AV[6][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

alsfhline=0.1	# alpha values for model Kim+17
alsfhfill=0.1
linest=[':','-.','--']
xgrid = np.arange(0.,3.,0.2)

h11.plot(xgrid[flag_E1], 2.5 * np.log10(min_bv0[0]/prof_E1[flag_E1]), color=cols_for_AV, linewidth=3)
h21.plot(xgrid[flag_S01], 2.5 * np.log10(min_bv0[1]/prof_S01[flag_S01]), color=cols_for_AV, linewidth=3)
h31.plot(xgrid[flag_Sa1], 2.5 * np.log10(min_bv0[2]/prof_Sa1[flag_Sa1]), color=cols_for_AV, linewidth=3)
h41.plot(xgrid[flag_Sb1], 2.5 * np.log10(min_bv0[3]/prof_Sb1[flag_Sb1]), color=cols_for_AV, linewidth=3)
h12.plot(xgrid[flag_E2], 2.5 * np.log10(min_bv0[4]/prof_E2[flag_E2]), color=cols_for_AV, linewidth=3)
h22.plot(xgrid[flag_S02], 2.5 * np.log10(min_bv0[5]/prof_S02[flag_S02]), color=cols_for_AV, linewidth=3)
h32.plot(xgrid[flag_Sa2], 2.5 * np.log10(min_bv0[6]/prof_Sa2[flag_Sa2]), color=cols_for_AV, linewidth=3)
h42.plot(xgrid[flag_Sb2], 2.5 * np.log10(min_bv0[7]/prof_Sb2[flag_Sb2]), color=cols_for_AV, linewidth=3,
         label=r'$\langle\beta_{V}\rangle$-profile fitting')

h11.fill_between(xgrid[flag_E1], 2.5 * np.log10(min_bv0[0]/prof_E1_hi[flag_E1]), 2.5 * np.log10(min_bv0[0] / prof_E1_low[flag_E1]), alpha=alcalfill, color=cols_for_AV)
h21.fill_between(xgrid[flag_S01], 2.5 * np.log10(min_bv0[1]/prof_S01_hi[flag_S01]), 2.5 * np.log10(min_bv0[1] / prof_S01_low[flag_S01]), alpha=alcalfill, color=cols_for_AV)
h31.fill_between(xgrid[flag_Sa1], 2.5 * np.log10(min_bv0[2]/prof_Sa1_hi[flag_Sa1]), 2.5 * np.log10(min_bv0[2] / prof_Sa1_low[flag_Sa1]), alpha=alcalfill, color=cols_for_AV)
h41.fill_between(xgrid[flag_Sb1], 2.5 * np.log10(min_bv0[3]/prof_Sb1_hi[flag_Sb1]), 2.5 * np.log10(min_bv0[3] / prof_Sb1_low[flag_Sb1]), alpha=alcalfill, color=cols_for_AV)
h12.fill_between(xgrid[flag_E2], 2.5 * np.log10(min_bv0[4]/prof_E2_hi[flag_E2]), 2.5 * np.log10(min_bv0[4] / prof_E2_low[flag_E2]), alpha=alcalfill, color=cols_for_AV)
h22.fill_between(xgrid[flag_S02], 2.5 * np.log10(min_bv0[5]/prof_S02_hi[flag_S02]), 2.5 * np.log10(min_bv0[5] / prof_S02_low[flag_S02]), alpha=alcalfill, color=cols_for_AV)
h32.fill_between(xgrid[flag_Sa2], 2.5 * np.log10(min_bv0[6]/prof_Sa2_hi[flag_Sa2]), 2.5 * np.log10(min_bv0[6] / prof_Sa2_low[flag_Sa2]), alpha=alcalfill, color=cols_for_AV)
h42.fill_between(xgrid[flag_Sb2], 2.5 * np.log10(min_bv0[7]/prof_Sb2_hi[flag_Sb2]), 2.5 * np.log10(min_bv0[7] / prof_Sb2_low[flag_Sb2]), alpha=alcalfill, color=cols_for_AV)


h11.set_ylim([0, 1])
h21.set_ylim([0, 1])
h31.set_ylim([0, 1])
h41.set_ylim([0, 1])
h12.set_ylim([0, 1])
h22.set_ylim([0, 1])
h32.set_ylim([0, 1])
h42.set_ylim([0, 1])

h11.set_xlim([0, 3])
h21.set_xlim([0, 3])
h31.set_xlim([0, 3])
h41.set_xlim([0, 3])
h12.set_xlim([0, 3])
h22.set_xlim([0, 3])
h32.set_xlim([0, 3])
h42.set_xlim([0, 3])

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

h11.text(1.5,0.475,'N=%d' % cnt_E1,size=20, color=cols_for_AV)
h21.text(1.5,0.475,'N=%d' % cnt_S01,size=20, color=cols_for_AV)
h31.text(1.5,0.475,'N=%d' % cnt_Sa1,size=20, color=cols_for_AV)
h41.text(1.5,0.475,'N=%d' % cnt_Sb1,size=20, color=cols_for_AV)
h12.text(1.5,0.475,'N=%d' % cnt_E2,size=20, color=cols_for_AV)
h22.text(1.5,0.475,'N=%d' % cnt_S02,size=20, color=cols_for_AV)
h32.text(1.5,0.475,'N=%d' % cnt_Sa2,size=20, color=cols_for_AV)
h42.text(1.5,0.475,'N=%d' % cnt_Sb2,size=20, color=cols_for_AV)

h11.text(0.6,0.725,r'$\beta_{V,0}$=%.2f' % min_bv0[0], size=20, color=cols_for_AV)
h21.text(0.6,0.725,r'$\beta_{V,0}$=%.2f' % min_bv0[1], size=20, color=cols_for_AV)
h31.text(0.6,0.725,r'$\beta_{V,0}$=%.2f' % min_bv0[2], size=20, color=cols_for_AV)
h41.text(0.6,0.725,r'$\beta_{V,0}$=%.2f' % min_bv0[3], size=20, color=cols_for_AV)
h12.text(0.8,0.725,r'$\beta_{V,0}$=%.2f' % min_bv0[4], size=20, color=cols_for_AV)
h22.text(0.6,0.725,r'$\beta_{V,0}$=%.2f' % min_bv0[5], size=20, color=cols_for_AV)
h32.text(0.6,0.725,r'$\beta_{V,0}$=%.2f' % min_bv0[6], size=20, color=cols_for_AV)
h42.text(0.6,0.725,r'$\beta_{V,0}$=%.2f' % min_bv0[7], size=20, color=cols_for_AV)

h11.text(1.0,0.6,r'$\chi^2$=%.2f' % min_chi[0], size=20, color=cols_for_AV)
h21.text(1.0,0.6,r'$\chi^2$=%.2f' % min_chi[1], size=20, color=cols_for_AV)
h31.text(1.0,0.6,r'$\chi^2$=%.2f' % min_chi[2], size=20, color=cols_for_AV)
h41.text(1.0,0.6,r'$\chi^2$=%.2f' % min_chi[3], size=20, color=cols_for_AV)
h12.text(1.0,0.6,r'$\chi^2$=%.2f' % min_chi[4], size=20, color=cols_for_AV)
h22.text(1.0,0.6,r'$\chi^2$=%.2f' % min_chi[5], size=20, color=cols_for_AV)
h32.text(1.0,0.6,r'$\chi^2$=%.2f' % min_chi[6], size=20, color=cols_for_AV)
h42.text(1.0,0.6,r'$\chi^2$=%.2f' % min_chi[7], size=20, color=cols_for_AV)

h11.set_ylabel(r'$A_{V}$',size=30)
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
h42.legend(bbox_to_anchor=(1.2, -0.05), prop={'size': 20}, frameon=False)

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

fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/dust_prof_all_T_BT_avg_betav_rebin.pdf')
plt.close(fig)

