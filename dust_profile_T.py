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

av_prof_E = np.zeros((46,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_E[:] = np.nan
av_prof_S0 = np.zeros((38,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_S0[:] = np.nan
av_prof_Sa = np.zeros((31,10,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_Sa[:] = np.nan
av_prof_Sb = np.zeros((43,10,30))        # For stacking Dust profiles for SFH4 w/ 6 Zs
av_prof_Sb[:] = np.nan
av_prof_Sbc = np.zeros((16,10,30))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T < 0.5
av_prof_Sbc[:] = np.nan
av_prof_Sc = np.zeros((44,10,30))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T > 0.5
av_prof_Sc[:] = np.nan
av_prof_Sd = np.zeros((15,10,30))        # For stacking Dust profiles for SFH5 w/ 6 Zs
av_prof_Sd[:] = np.nan

cnt_E = 0	# counter for SFH3
cnt_S0 = 0	# counter for SFH3 with B/T < div_btot_sfh3
cnt_Sa = 0	# counter for SFH3 with B/T < div_btot_sfh3
cnt_Sb = 0	# counter for SFH4
cnt_Sbc = 0	# counter for SFH4 with B/T < div_btot_sfh4
cnt_Sc = 0	# counter for SFH4 with B/T > div_btot_sfh4
cnt_Sd = 0	# counter for SFH5

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

    x2, m2, y2, s2 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_-4.txt')
    x3, m3, y3, s3 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_-3.txt')
    x4, m4, y4, s4 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_-2.txt')
    x5, m5, y5, s5 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_-1.txt')
    x6, m6, y6, s6 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_0.txt')
    x7, m7, y7, s7 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_1.txt')
    x8, m8, y8, s8 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_2.txt')
    x9, m9, y9, s9 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_3.txt')
    x10, m10, y10, s10 = np.loadtxt(work_dir+'ellipse/scratch_grid/'+name+'_AV_4.txt')

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

    if T < -3:
        av_prof_E[cnt_E, 0, :] = f1(xnew)
        av_prof_E[cnt_E, 1, :] = f2(xnew)
        av_prof_E[cnt_E, 2, :] = f3(xnew)
        av_prof_E[cnt_E, 3, :] = f4(xnew)
        av_prof_E[cnt_E, 4, :] = f5(xnew)
        av_prof_E[cnt_E, 5, :] = f6(xnew)
        av_prof_E[cnt_E, 6, :] = f7(xnew)
        av_prof_E[cnt_E, 7, :] = f8(xnew)
        av_prof_E[cnt_E, 8, :] = f9(xnew)
        av_prof_E[cnt_E, 9, :] = f10(xnew)
        cnt_E = cnt_E + 1
    if T < 0 and T >= -3:
        av_prof_S0[cnt_S0, 0, :] = f1(xnew)
        av_prof_S0[cnt_S0, 1, :] = f2(xnew)
        av_prof_S0[cnt_S0, 2, :] = f3(xnew)
        av_prof_S0[cnt_S0, 3, :] = f4(xnew)
        av_prof_S0[cnt_S0, 4, :] = f5(xnew)
        av_prof_S0[cnt_S0, 5, :] = f6(xnew)
        av_prof_S0[cnt_S0, 6, :] = f7(xnew)
        av_prof_S0[cnt_S0, 7, :] = f8(xnew)
        av_prof_S0[cnt_S0, 8, :] = f9(xnew)
        av_prof_S0[cnt_S0, 9, :] = f10(xnew)
        cnt_S0 = cnt_S0 + 1
    if T < 2 and T >= 0:
        av_prof_Sa[cnt_Sa, 0, :] = f1(xnew)
        av_prof_Sa[cnt_Sa, 1, :] = f2(xnew)
        av_prof_Sa[cnt_Sa, 2, :] = f3(xnew)
        av_prof_Sa[cnt_Sa, 3, :] = f4(xnew)
        av_prof_Sa[cnt_Sa, 4, :] = f5(xnew)
        av_prof_Sa[cnt_Sa, 5, :] = f6(xnew)
        av_prof_Sa[cnt_Sa, 6, :] = f7(xnew)
        av_prof_Sa[cnt_Sa, 7, :] = f8(xnew)
        av_prof_Sa[cnt_Sa, 8, :] = f9(xnew)
        av_prof_Sa[cnt_Sa, 9, :] = f10(xnew)
        cnt_Sa = cnt_Sa + 1
    if T < 4 and T >= 2:
        av_prof_Sb[cnt_Sb, 0, :] = f1(xnew)
        av_prof_Sb[cnt_Sb, 1, :] = f2(xnew)
        av_prof_Sb[cnt_Sb, 2, :] = f3(xnew)
        av_prof_Sb[cnt_Sb, 3, :] = f4(xnew)
        av_prof_Sb[cnt_Sb, 4, :] = f5(xnew)
        av_prof_Sb[cnt_Sb, 5, :] = f6(xnew)
        av_prof_Sb[cnt_Sb, 6, :] = f7(xnew)
        av_prof_Sb[cnt_Sb, 7, :] = f8(xnew)
        av_prof_Sb[cnt_Sb, 8, :] = f9(xnew)
        av_prof_Sb[cnt_Sb, 9, :] = f10(xnew)
        cnt_Sb = cnt_Sb + 1
    if T < 5 and T >= 4:
        av_prof_Sbc[cnt_Sbc, 0, :] = f1(xnew)
        av_prof_Sbc[cnt_Sbc, 1, :] = f2(xnew)
        av_prof_Sbc[cnt_Sbc, 2, :] = f3(xnew)
        av_prof_Sbc[cnt_Sbc, 3, :] = f4(xnew)
        av_prof_Sbc[cnt_Sbc, 4, :] = f5(xnew)
        av_prof_Sbc[cnt_Sbc, 5, :] = f6(xnew)
        av_prof_Sbc[cnt_Sbc, 6, :] = f7(xnew)
        av_prof_Sbc[cnt_Sbc, 7, :] = f8(xnew)
        av_prof_Sbc[cnt_Sbc, 8, :] = f9(xnew)
        av_prof_Sbc[cnt_Sbc, 9, :] = f10(xnew)
        cnt_Sbc = cnt_Sbc + 1
    if T < 7 and T >= 5:
        av_prof_Sc[cnt_Sc, 0, :] = f1(xnew)
        av_prof_Sc[cnt_Sc, 1, :] = f2(xnew)
        av_prof_Sc[cnt_Sc, 2, :] = f3(xnew)
        av_prof_Sc[cnt_Sc, 3, :] = f4(xnew)
        av_prof_Sc[cnt_Sc, 4, :] = f5(xnew)
        av_prof_Sc[cnt_Sc, 5, :] = f6(xnew)
        av_prof_Sc[cnt_Sc, 6, :] = f7(xnew)
        av_prof_Sc[cnt_Sc, 7, :] = f8(xnew)
        av_prof_Sc[cnt_Sc, 8, :] = f9(xnew)
        av_prof_Sc[cnt_Sc, 9, :] = f10(xnew)
        cnt_Sc = cnt_Sc + 1
    if T < 9 and T >= 7:
        av_prof_Sd[cnt_Sd, 0, :] = f1(xnew)
        av_prof_Sd[cnt_Sd, 1, :] = f2(xnew)
        av_prof_Sd[cnt_Sd, 2, :] = f3(xnew)
        av_prof_Sd[cnt_Sd, 3, :] = f4(xnew)
        av_prof_Sd[cnt_Sd, 4, :] = f5(xnew)
        av_prof_Sd[cnt_Sd, 5, :] = f6(xnew)
        av_prof_Sd[cnt_Sd, 6, :] = f7(xnew)
        av_prof_Sd[cnt_Sd, 7, :] = f8(xnew)
        av_prof_Sd[cnt_Sd, 8, :] = f9(xnew)
        av_prof_Sd[cnt_Sd, 9, :] = f10(xnew)
        cnt_Sd = cnt_Sd + 1

fig = plt.figure(figsize=(10,12))

h1 = fig.add_axes([0.1,0.65,0.275,0.275])
h2 = fig.add_axes([0.375,0.65,0.275,0.275])
h3 = fig.add_axes([0.65,0.65,0.275,0.275])
h4 = fig.add_axes([0.1,0.375,0.275,0.275])
h5 = fig.add_axes([0.375,0.375,0.275,0.275])
h6 = fig.add_axes([0.65,0.375,0.275,0.275])
h7 = fig.add_axes([0.1,0.1,0.275,0.275])

alcalline=0.3	# alpha values for data from CALIFA
alcalfill=0.1
h1.plot(np.arange(0.,3.,0.1),rosa_AV[0][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h1.fill_between(np.arange(0.,3.,0.1),rosa_AV[0][3:]-rosa_AV[0][2],rosa_AV[0][3:]+rosa_AV[0][2],alpha=alcalfill,color='black')
h1.errorbar(1, rosa_AV[0][12], yerr=rosa_AV[0][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h2.plot(np.arange(0.,3.,0.1),rosa_AV[2][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h2.fill_between(np.arange(0.,3.,0.1),rosa_AV[2][3:]-rosa_AV[2][2],rosa_AV[2][3:]+rosa_AV[2][2],alpha=alcalfill,color='black')
h2.errorbar(1, rosa_AV[2][12], yerr=rosa_AV[2][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h3.plot(np.arange(0.,3.,0.1),rosa_AV[4][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h3.fill_between(np.arange(0.,3.,0.1),rosa_AV[4][3:]-rosa_AV[4][2],rosa_AV[4][3:]+rosa_AV[4][2],alpha=alcalfill,color='black')
h3.errorbar(1, rosa_AV[4][12], yerr=rosa_AV[4][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h4.plot(np.arange(0.,3.,0.1),rosa_AV[6][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h4.fill_between(np.arange(0.,3.,0.1),rosa_AV[6][3:]-rosa_AV[6][2],rosa_AV[6][3:]+rosa_AV[6][2],alpha=alcalfill,color='black')
h4.errorbar(1, rosa_AV[6][12], yerr=rosa_AV[6][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h5.plot(np.arange(0.,3.,0.1),rosa_AV[8][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h5.fill_between(np.arange(0.,3.,0.1),rosa_AV[8][3:]-rosa_AV[8][2],rosa_AV[8][3:]+rosa_AV[8][2],alpha=alcalfill,color='black')
h5.errorbar(1, rosa_AV[8][12], yerr=rosa_AV[8][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h6.plot(np.arange(0.,3.,0.1),rosa_AV[10][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
#h6.fill_between(np.arange(0.,3.,0.1),rosa_AV[10][3:]-rosa_AV[10][2],rosa_AV[10][3:]+rosa_AV[10][2],alpha=alcalfill,color='black')
h6.errorbar(1, rosa_AV[10][12], yerr=rosa_AV[10][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h7.plot(np.arange(0.,3.,0.1),rosa_AV[12][3:],color='black',alpha=alcalline, linewidth=3, label='GD15 (CALIFA)', linestyle=linest_for_CALIFA)
#h7.fill_between(np.arange(0.,3.,0.1),rosa_AV[12][3:]-rosa_AV[12][2],rosa_AV[12][3:]+rosa_AV[12][2],alpha=alcalfill,color='black')
h7.errorbar(1, rosa_AV[12][12], yerr=rosa_AV[12][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

# h1.plot([0.6, 0.85], [0.85, 0.85], color='k', alpha=alcalline)
# h2.plot([1.0, 1.3], [0.85, 0.85], color='k', alpha=alcalline, marker='.')
# h3.plot([0.2, 0.5], [0.85, 0.85], color='k', alpha=alcalline, marker='+')
# h4.plot([0.7, 1.1], [0.85, 0.85], color='k', alpha=alcalline, marker='*')
# h5.plot([1.3, 1.7], [0.85, 0.85], color='k', alpha=alcalline, marker='s')
# h6.plot([0.5, 0.85], [0.85, 0.85], color='k', alpha=alcalline, marker='d')
# h7.plot([1.1, 1.4], [0.85, 0.85], color='k', alpha=alcalline, marker='p')

min_idx = np.zeros(7)
min_chi = np.zeros(7)
min_chi[:] = 1e6
for i in range(0, 10):
    prof_E = np.nanmean(av_prof_E[:, i, :], axis=0)
    prof_S0 = np.nanmean(av_prof_S0[:, i, :], axis=0)
    prof_Sa = np.nanmean(av_prof_Sa[:, i, :], axis=0)
    prof_Sb = np.nanmean(av_prof_Sb[:, i, :], axis=0)
    prof_Sbc = np.nanmean(av_prof_Sbc[:, i, :], axis=0)
    prof_Sc = np.nanmean(av_prof_Sc[:, i, :], axis=0)
    prof_Sd = np.nanmean(av_prof_Sd[:, i, :], axis=0)

    prof_E_err = np.nanstd(av_prof_E[:, i, :], axis=0)
    prof_S0_err = np.nanstd(av_prof_S0[:, i, :], axis=0)
    prof_Sa_err = np.nanstd(av_prof_Sa[:, i, :], axis=0)
    prof_Sb_err = np.nanstd(av_prof_Sb[:, i, :], axis=0)
    prof_Sbc_err = np.nanstd(av_prof_Sbc[:, i, :], axis=0)
    prof_Sc_err = np.nanstd(av_prof_Sc[:, i, :], axis=0)
    prof_Sd_err = np.nanstd(av_prof_Sd[:, i, :], axis=0)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_E[j]) and ~np.isnan(prof_E_err[j]):
            chi += np.abs(prof_E[j] - rosa_AV[0][3+j]) / prof_E_err[j]
    if chi < min_chi[0]:
        min_chi[0] = chi
        min_idx[0] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_S0[j]) and ~np.isnan(prof_S0_err[j]):
            chi += np.abs(prof_S0[j] - rosa_AV[2][3 + j]) / prof_S0_err[j]
    if chi < min_chi[1]:
        min_chi[1] = chi
        min_idx[1] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sa[j]) and ~np.isnan(prof_Sa_err[j]):
            chi += np.abs(prof_Sa[j] - rosa_AV[4][3 + j]) / prof_Sa_err[j]
    if chi < min_chi[2]:
        min_chi[2] = chi
        min_idx[2] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sb[j]) and ~np.isnan(prof_Sb_err[j]):
            chi += np.abs(prof_Sb[j] - rosa_AV[6][3 + j]) / prof_Sb_err[j]
    if chi < min_chi[3]:
        min_chi[3] = chi
        min_idx[3] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sbc[j]) and ~np.isnan(prof_Sbc_err[j]):
            chi += np.abs(prof_Sbc[j] - rosa_AV[8][3 + j]) / prof_Sbc_err[j]
    if chi < min_chi[4]:
        min_chi[4] = chi
        min_idx[4] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sc[j]) and ~np.isnan(prof_Sc_err[j]):
            chi += np.abs(prof_Sc[j] - rosa_AV[10][3 + j]) / prof_Sc_err[j]
    if chi < min_chi[5]:
        min_chi[5] = chi
        min_idx[5] = int(i)

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sd[j]) and ~np.isnan(prof_Sd_err[j]):
            chi += np.abs(prof_Sd[j] - rosa_AV[12][3 + j]) / prof_Sd_err[j]
    if chi < min_chi[6]:
        min_chi[6] = chi
        min_idx[6] = int(i)

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

h1.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_E[:, int(min_idx[0]), :], axis=0), color=cols_for_AV, linewidth=3)
h2.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_S0[:, int(min_idx[1]), :], axis=0), color=cols_for_AV, linewidth=3)
h3.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sa[:, int(min_idx[2]), :], axis=0), color=cols_for_AV, linewidth=3)
h4.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sb[:, int(min_idx[3]), :], axis=0), color=cols_for_AV, linewidth=3)
h5.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sbc[:, int(min_idx[4]), :], axis=0), color=cols_for_AV, linewidth=3)
h6.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sc[:, int(min_idx[5]), :], axis=0), color=cols_for_AV, linewidth=3)
h7.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sd[:, int(min_idx[6]), :], axis=0), color=cols_for_AV, linewidth=3, label=r'$\beta_{V}$ method')

h1.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_E[:, int(min_idx[0]), :], axis=0) - np.nanstd(av_prof_E[:, int(min_idx[0]), :], axis=0),
    np.nanmean(av_prof_E[:, int(min_idx[0]), :], axis=0) + np.nanstd(av_prof_E[:, int(min_idx[0]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h2.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_S0[:, int(min_idx[1]), :], axis=0) - np.nanstd(av_prof_S0[:, int(min_idx[1]), :], axis=0),
    np.nanmean(av_prof_S0[:, int(min_idx[1]), :], axis=0) + np.nanstd(av_prof_S0[:, int(min_idx[1]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h3.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sa[:, int(min_idx[2]), :], axis=0) - np.nanstd(av_prof_Sa[:, int(min_idx[2]), :], axis=0),
    np.nanmean(av_prof_Sa[:, int(min_idx[2]), :], axis=0) + np.nanstd(av_prof_Sa[:, int(min_idx[2]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h4.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sb[:, int(min_idx[3]), :], axis=0) - np.nanstd(av_prof_Sb[:, int(min_idx[3]), :], axis=0),
    np.nanmean(av_prof_Sb[:, int(min_idx[3]), :], axis=0) + np.nanstd(av_prof_Sb[:, int(min_idx[3]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h5.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sbc[:, int(min_idx[4]), :], axis=0) - np.nanstd(av_prof_Sbc[:, int(min_idx[4]), :], axis=0),
    np.nanmean(av_prof_Sbc[:, int(min_idx[4]), :], axis=0) + np.nanstd(av_prof_Sbc[:, int(min_idx[4]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h6.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sc[:, int(min_idx[5]), :], axis=0) - np.nanstd(av_prof_Sc[:, int(min_idx[5]), :], axis=0),
    np.nanmean(av_prof_Sc[:, int(min_idx[5]), :], axis=0) + np.nanstd(av_prof_Sc[:, int(min_idx[5]), :], axis=0), alpha=alcalfill, color=cols_for_AV)
h7.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sd[:, int(min_idx[6]), :], axis=0) - np.nanstd(av_prof_Sd[:, int(min_idx[6]), :], axis=0),
    np.nanmean(av_prof_Sd[:, int(min_idx[6]), :], axis=0) + np.nanstd(av_prof_Sd[:, int(min_idx[6]), :], axis=0), alpha=alcalfill, color=cols_for_AV)


h1.set_ylim([0,1])
h2.set_ylim([0,1])
h3.set_ylim([0,1])
h4.set_ylim([0,1])
h5.set_ylim([0,1])
h6.set_ylim([0,1])
h7.set_ylim([0,1])

h1.set_xlim([0,3])
h2.set_xlim([0,3])
h3.set_xlim([0,3])
h4.set_xlim([0,3])
h5.set_xlim([0,3])
h6.set_xlim([0,3])
h7.set_xlim([0,3])

h1.text(0.3,0.85,'E',size=30)
h2.text(0.3,0.85,'S0',size=30)
h3.text(0.3,0.85,'Sa',size=30)
h4.text(0.3,0.85,'Sb',size=30)
h5.text(0.3,0.85,'Sbc',size=30)
h6.text(0.3,0.85,'Sc',size=30)
h7.text(0.3,0.85,'Sd',size=30)

h1.text(0.3,0.65,'N=%d' % cnt_E,size=20, color=cols_for_AV)
h2.text(0.3,0.65,'N=%d' % cnt_S0,size=20, color=cols_for_AV)
h3.text(0.3,0.65,'N=%d' % cnt_Sa,size=20, color=cols_for_AV)
h4.text(0.3,0.65,'N=%d' % cnt_Sb,size=20, color=cols_for_AV)
h5.text(0.3,0.65,'N=%d' % cnt_Sbc,size=20, color=cols_for_AV)
h6.text(0.3,0.65,'N=%d' % cnt_Sc,size=20, color=cols_for_AV)
h7.text(0.3,0.65,'N=%d' % cnt_Sd,size=20, color=cols_for_AV)

h1.text(0.3,0.75,r'$\beta_{V,0}$=0.5',size=20, color=cols_for_AV)
h2.text(0.3,0.75,r'$\beta_{V,0}$=0.6',size=20, color=cols_for_AV)
h3.text(0.3,0.75,r'$\beta_{V,0}$=0.73',size=20, color=cols_for_AV)
h4.text(0.3,0.75,r'$\beta_{V,0}$=0.73',size=20, color=cols_for_AV)
h5.text(0.3,0.75,r'$\beta_{V,0}$=0.81',size=20, color=cols_for_AV)
h6.text(0.3,0.75,r'$\beta_{V,0}$=0.85',size=20, color=cols_for_AV)
h7.text(0.3,0.75,r'$\beta_{V,0}$=0.97',size=20, color=cols_for_AV)

h4.set_ylabel(r'$\langle A_{V}\rangle$',size=30)
h5.set_xlabel(r'$R/R_{50,V}^{P}$',size=30)
h6.set_xlabel(r'$R/R_{50,V}^{P}$',size=30)
h7.set_xlabel(r'$R/R_{50,V}^{P}$',size=30)
h4.set_yticks([0,0.2,0.4,0.6,0.8])
h7.set_yticks([0,0.2,0.4,0.6,0.8])
h7.set_xticks([0,1,2,3])
h5.set_xticks([1,2,3])
h6.set_xticks([1,2,3])
h1.tick_params(labelsize='20',direction='in',top=True,right=True)
h2.tick_params(direction='in',top=True,right=True,bottom=True,left=True)
h3.tick_params(direction='in',top=True,right=True,bottom=True,left=True)
h4.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h5.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h6.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h7.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)

h7.legend(bbox_to_anchor=(1.05, 0.5), prop={'size': 30})

h1.set_xticklabels([])
h2.set_xticklabels([])
h2.set_yticklabels([])
h3.set_xticklabels([])
h3.set_yticklabels([])
h4.set_xticklabels([])
h5.set_yticklabels([])
h6.set_yticklabels([])

fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/dust_prof_all_T.pdf')
plt.close(fig)

