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
from matplotlib.lines import Line2D
from pandas import DataFrame, read_csv
import pandas as pd
from operator import itemgetter
from scipy.interpolate import interp1d
from astropy.stats import sigma_clipped_stats

import matplotlib.pyplot as plt

htypes = ['E0','S0','S0/a','Sa','Sab','Sb','Sbc','Sc','Scd','Sd','Sdm','Sm','Im','?'] # Hubble types corresponding Nair+10 T-Type
agntypes = ['SF','transition/mixed AGN','Seyfert','LINER']	# Kauffmann+03 from Nair+10

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

text_for_bv0s=['-0.05', '-0.04', '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03', '0.04']
cols_for_bv0s=['black', 'gray', 'silver', 'firebrick', 'sienna', 'sandybrown', 'tan', 'olivedrab', 'palegreen', 'navy']
cols_for_AV = 'maroon'
linest_for_CALIFA = ':'

#hubble_bv0 = [0.56, 0.66, 0.68, 0.69, 0.87, 1.00, 1.16]
hubble_bv0 = [0.56, 0.66, 0.74, 0.77, 0.94, 0.97, 1.20]
#hubble_bv0_err = [0.14, 0.2, 0.28, 0.25, 0.21, 0.34, 0.34]

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

bv_prof_E = np.zeros((46,15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_E[:] = np.nan
nan_cnt_E = np.zeros(15)        # for counting NaN values in the stack
bv_prof_S0 = np.zeros((38,15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_S0[:] = np.nan
nan_cnt_S0 = np.zeros(15)        # for counting NaN values in the stack
bv_prof_Sa = np.zeros((31,15))	# For stacking Dust profiles for SFH3 w/ 6 Zs
bv_prof_Sa[:] = np.nan
nan_cnt_Sa = np.zeros(15)        # for counting NaN values in the stack
bv_prof_Sb = np.zeros((43,15))        # For stacking Dust profiles for SFH4 w/ 6 Zs
bv_prof_Sb[:] = np.nan
nan_cnt_Sb = np.zeros(15)        # for counting NaN values in the stack
bv_prof_Sbc = np.zeros((16,15))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T < 0.5
bv_prof_Sbc[:] = np.nan
nan_cnt_Sbc = np.zeros(15)        # for counting NaN values in the stack
bv_prof_Sc = np.zeros((45,15))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T > 0.5
bv_prof_Sc[:] = np.nan
nan_cnt_Sc = np.zeros(15)        # for counting NaN values in the stack
bv_prof_Sd = np.zeros((15,15))        # For stacking Dust profiles for SFH5 w/ 6 Zs
bv_prof_Sd[:] = np.nan
nan_cnt_Sd = np.zeros(15)        # for counting NaN values in the stack

cnt_E = 0	# counter for SFH3
cnt_S0 = 0	# counter for SFH3 with B/T < div_btot_sfh3
cnt_Sa = 0	# counter for SFH3 with B/T < div_btot_sfh3
cnt_Sb = 0	# counter for SFH4
cnt_Sbc = 0	# counter for SFH4 with B/T < div_btot_sfh4
cnt_Sc = 0	# counter for SFH4 with B/T > div_btot_sfh4
cnt_Sd = 0	# counter for SFH5

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
    btot = ttype_match.iloc[0,7]
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

    if T < -3:
        bv_prof_E[cnt_E,:] = f1(xnew)
        nan_cnt_E += np.isnan(f1(xnew))
        cnt_E = cnt_E + 1
    if T < 0 and T >= -3:
        bv_prof_S0[cnt_S0, :] = f1(xnew)
        nan_cnt_S0 += np.isnan(f1(xnew))
        cnt_S0 = cnt_S0 + 1
    if T < 2 and T >= 0:
        bv_prof_Sa[cnt_Sa, :] = f1(xnew)
        nan_cnt_Sa += np.isnan(f1(xnew))
        cnt_Sa = cnt_Sa + 1
    if T < 4 and T >= 2:
        bv_prof_Sb[cnt_Sb, :] = f1(xnew)
        nan_cnt_Sb += np.isnan(f1(xnew))
        cnt_Sb = cnt_Sb + 1
    if T < 5 and T >= 4:
        bv_prof_Sbc[cnt_Sbc, :] = f1(xnew)
        nan_cnt_Sbc += np.isnan(f1(xnew))
        cnt_Sbc = cnt_Sbc + 1
        print(name)
    if T < 7 and T >= 5:
        bv_prof_Sc[cnt_Sc, :] = f1(xnew)
        nan_cnt_Sc += np.isnan(f1(xnew))
        cnt_Sc = cnt_Sc + 1
    if T < 9 and T >= 7:
        bv_prof_Sd[cnt_Sd, :] = f1(xnew)
        nan_cnt_Sd += np.isnan(f1(xnew))
        cnt_Sd = cnt_Sd + 1

flag_E = np.where(nan_cnt_E <= (cnt_E-5))[0]
flag_S0 = np.where(nan_cnt_S0 <= (cnt_S0-5))[0]
flag_Sa = np.where(nan_cnt_Sa <= (cnt_Sa-5))[0]
flag_Sb = np.where(nan_cnt_Sb <= (cnt_Sb-5))[0]
flag_Sbc = np.where(nan_cnt_Sbc <= (cnt_Sbc-5))[0]
flag_Sc = np.where(nan_cnt_Sc <= (cnt_Sc-5))[0]
flag_Sd = np.where(nan_cnt_Sd <= (cnt_Sd-5))[0]

sig_E = sigma_clipped_stats(bv_prof_E[:, :], sigma=2, axis=0)
sig_S0 = sigma_clipped_stats(bv_prof_S0[:, :], sigma=2, axis=0)
sig_Sa = sigma_clipped_stats(bv_prof_Sa[:, :], sigma=2, axis=0)
sig_Sb = sigma_clipped_stats(bv_prof_Sb[:, :], sigma=2, axis=0)
sig_Sbc = sigma_clipped_stats(bv_prof_Sbc[:, :], sigma=2, axis=0)
sig_Sc = sigma_clipped_stats(bv_prof_Sc[:, :], sigma=2, axis=0)
sig_Sd = sigma_clipped_stats(bv_prof_Sd[:, :], sigma=2, axis=0)

med_E = np.nanmedian(bv_prof_E[:, :], axis=0)
med_S0 = np.nanmedian(bv_prof_S0[:, :], axis=0)
med_Sa = np.nanmedian(bv_prof_Sa[:, :], axis=0)
med_Sb = np.nanmedian(bv_prof_Sb[:, :], axis=0)
med_Sbc = np.nanmedian(bv_prof_Sbc[:, :], axis=0)
med_Sc = np.nanmedian(bv_prof_Sc[:, :], axis=0)
med_Sd = np.nanmedian(bv_prof_Sd[:, :], axis=0)

prof_E = med_E
prof_S0 = med_S0
prof_Sa = med_Sa
prof_Sb = med_Sb
prof_Sbc = med_Sbc
prof_Sc = med_Sc
prof_Sd = med_Sd

prof_E_low = np.nanpercentile(bv_prof_E[:, :],  15.87, axis=0)
prof_S0_low = np.nanpercentile(bv_prof_S0[:, :], 15.87, axis=0)
prof_Sa_low = np.nanpercentile(bv_prof_Sa[:, :], 15.87, axis=0)
prof_Sb_low = np.nanpercentile(bv_prof_Sb[:, :], 15.87, axis=0)
prof_Sbc_low = np.nanpercentile(bv_prof_Sbc[:, :], 15.87, axis=0)
prof_Sc_low = np.nanpercentile(bv_prof_Sc[:, :], 15.87, axis=0)
prof_Sd_low = np.nanpercentile(bv_prof_Sd[:, :], 15.87, axis=0)

prof_E_hi = np.nanpercentile(bv_prof_E[:, :], 84.13, axis=0)
prof_S0_hi = np.nanpercentile(bv_prof_S0[:, :], 84.13, axis=0)
prof_Sa_hi = np.nanpercentile(bv_prof_Sa[:, :], 84.13, axis=0)
prof_Sb_hi = np.nanpercentile(bv_prof_Sb[:, :], 84.13, axis=0)
prof_Sbc_hi = np.nanpercentile(bv_prof_Sbc[:, :], 84.13, axis=0)
prof_Sc_hi = np.nanpercentile(bv_prof_Sc[:, :], 84.13, axis=0)
prof_Sd_hi = np.nanpercentile(bv_prof_Sd[:, :], 84.13, axis=0)

min_bv0 = np.zeros(7)
min_chi = np.zeros(7)
min_chi[:] = 1e6

for i in np.arange(-0.1, 0.1, 0.01):
    av_prof_E = 2.5 * np.log10((hubble_bv0[0] + i) / prof_E)
    av_prof_S0 = 2.5 * np.log10((hubble_bv0[1] + i) / prof_S0)
    av_prof_Sa = 2.5 * np.log10((hubble_bv0[2] + i) / prof_Sa)
    av_prof_Sb = 2.5 * np.log10((hubble_bv0[3] + i) / prof_Sb)
    av_prof_Sbc = 2.5 * np.log10((hubble_bv0[4] + i) / prof_Sbc)
    av_prof_Sc = 2.5 * np.log10((hubble_bv0[5] + i) / prof_Sc)
    av_prof_Sd = 2.5 * np.log10((hubble_bv0[6] + i) / prof_Sd)

    av_prof_E[np.where(av_prof_E < 0)] = 0.0
    av_prof_S0[np.where(av_prof_S0 < 0)] = 0.0
    av_prof_Sa[np.where(av_prof_Sa < 0)] = 0.0
    av_prof_Sb[np.where(av_prof_Sb < 0)] = 0.0
    av_prof_Sbc[np.where(av_prof_Sbc < 0)] = 0.0
    av_prof_Sc[np.where(av_prof_Sc < 0)] = 0.0
    av_prof_Sd[np.where(av_prof_Sd < 0)] = 0.0

    chi = 0
    for j in flag_E:
        sig = 1.25 * np.log10(prof_E_low[j] / prof_E_hi[j]) # half of magenta range in Figure 6
        chi += (av_prof_E[j] - rosa_AV[0][3+j*2])**2 / sig**2
    if chi < min_chi[0]:
        min_chi[0] = chi
        min_bv0[0] = hubble_bv0[0] + i

    chi = 0
    for j in flag_S0:
        sig = 1.25 * np.log10(prof_S0_low[j] / prof_S0_hi[j])
        chi += (av_prof_S0[j] - rosa_AV[2][3 + j*2])**2 / sig**2
    if chi < min_chi[1]:
        min_chi[1] = chi
        min_bv0[1] = hubble_bv0[1] + i

    chi = 0
    for j in flag_Sa:
        sig = 1.25 * np.log10(prof_Sa_low[j] / prof_Sa_hi[j])
        chi += (av_prof_Sa[j] - rosa_AV[4][3 + j*2])**2 / sig**2
    if chi < min_chi[2]:
        min_chi[2] = chi
        min_bv0[2] = hubble_bv0[2] + i

    chi = 0
    for j in flag_Sb:
        sig = 1.25 * np.log10(prof_Sb_low[j] / prof_Sb_hi[j])
        chi += (av_prof_Sb[j] - rosa_AV[6][3 + j*2])**2 / sig**2
    if chi < min_chi[3]:
        min_chi[3] = chi
        min_bv0[3] = hubble_bv0[3] + i

    chi = 0
    for j in flag_Sbc:
        sig = 1.25 * np.log10(prof_Sbc_low[j] / prof_Sbc_hi[j])
        chi += (av_prof_Sbc[j] - rosa_AV[8][3 + j*2])**2 / sig**2
    if chi < min_chi[4]:
        min_chi[4] = chi
        min_bv0[4] = hubble_bv0[4] + i

    chi = 0
    for j in flag_Sc:
        sig = 1.25 * np.log10(prof_Sc_low[j] / prof_Sc_hi[j])
        chi += (av_prof_Sc[j] - rosa_AV[10][3 + j*2])**2 / sig**2
    if chi < min_chi[5]:
        min_chi[5] = chi
        min_bv0[5] = hubble_bv0[5] + i

    chi = 0
    for j in flag_Sd:
        sig = 1.25 * np.log10(prof_Sd_low[j] / prof_Sd_hi[j])
        chi += (av_prof_Sd[j] - rosa_AV[12][3 + j*2])**2 / sig**2
    if chi < min_chi[6]:
        min_chi[6] = chi
        min_bv0[6] = hubble_bv0[6] + i



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

xgrid = np.arange(0.,3.,0.1)

h1.plot(xgrid,rosa_AV[0][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h1.errorbar(1, rosa_AV[0][12], yerr=rosa_AV[0][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h2.plot(xgrid,rosa_AV[2][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h2.errorbar(1, rosa_AV[2][12], yerr=rosa_AV[2][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h3.plot(xgrid,rosa_AV[4][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h3.errorbar(1, rosa_AV[4][12], yerr=rosa_AV[4][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h4.plot(xgrid,rosa_AV[6][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h4.errorbar(1, rosa_AV[6][12], yerr=rosa_AV[6][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h5.plot(xgrid,rosa_AV[8][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h5.errorbar(1, rosa_AV[8][12], yerr=rosa_AV[8][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h6.plot(xgrid,rosa_AV[10][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h6.errorbar(1, rosa_AV[10][12], yerr=rosa_AV[10][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h7.plot(xgrid,rosa_AV[12][3:],color='black',alpha=alcalline, linewidth=3, label='GD15 (CALIFA)', linestyle=linest_for_CALIFA)
h7.errorbar(1, rosa_AV[12][12], yerr=rosa_AV[12][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

xgrid = np.arange(0.,3.,0.2)

h1.plot(xgrid[flag_E], 2.5 * np.log10(min_bv0[0]/prof_E[flag_E]), color=cols_for_AV, linewidth=3)
h2.plot(xgrid[flag_S0], 2.5 * np.log10(min_bv0[1]/prof_S0[flag_S0]), color=cols_for_AV, linewidth=3)
h3.plot(xgrid[flag_Sa], 2.5 * np.log10(min_bv0[2]/prof_Sa[flag_Sa]), color=cols_for_AV, linewidth=3)
h4.plot(xgrid[flag_Sb], 2.5 * np.log10(min_bv0[3]/prof_Sb[flag_Sb]), color=cols_for_AV, linewidth=3)
h5.plot(xgrid[flag_Sbc], 2.5 * np.log10(min_bv0[4]/prof_Sbc[flag_Sbc]), color=cols_for_AV, linewidth=3)
h6.plot(xgrid[flag_Sc], 2.5 * np.log10(min_bv0[5]/prof_Sc[flag_Sc]), color=cols_for_AV, linewidth=3)
h7.plot(xgrid[flag_Sd], 2.5 * np.log10(min_bv0[6]/prof_Sd[flag_Sd]), color=cols_for_AV, linewidth=3,
        label=r'$\langle\beta_{V}\rangle$-profile fitting')

h1.fill_between(xgrid[flag_E], 2.5 * np.log10(min_bv0[0]/prof_E_hi[flag_E]), 2.5 * np.log10(min_bv0[0] / prof_E_low[flag_E]), alpha=alcalfill, color=cols_for_AV)
h2.fill_between(xgrid[flag_S0], 2.5 * np.log10(min_bv0[1]/prof_S0_hi[flag_S0]), 2.5 * np.log10(min_bv0[1] / prof_S0_low[flag_S0]), alpha=alcalfill, color=cols_for_AV)
h3.fill_between(xgrid[flag_Sa], 2.5 * np.log10(min_bv0[2]/prof_Sa_hi[flag_Sa]), 2.5 * np.log10(min_bv0[2] / prof_Sa_low[flag_Sa]), alpha=alcalfill, color=cols_for_AV)
h4.fill_between(xgrid[flag_Sb], 2.5 * np.log10(min_bv0[3]/prof_Sb_hi[flag_Sb]), 2.5 * np.log10(min_bv0[3] / prof_Sb_low[flag_Sb]), alpha=alcalfill, color=cols_for_AV)
h5.fill_between(xgrid[flag_Sbc], 2.5 * np.log10(min_bv0[4]/prof_Sbc_hi[flag_Sbc]), 2.5 * np.log10(min_bv0[4] / prof_Sbc_low[flag_Sbc]), alpha=alcalfill, color=cols_for_AV)
h6.fill_between(xgrid[flag_Sc], 2.5 * np.log10(min_bv0[5]/prof_Sc_hi[flag_Sc]), 2.5 * np.log10(min_bv0[5] / prof_Sc_low[flag_Sc]), alpha=alcalfill, color=cols_for_AV)
h7.fill_between(xgrid[flag_Sd], 2.5 * np.log10(min_bv0[6]/prof_Sd_hi[flag_Sd]), 2.5 * np.log10(min_bv0[6] / prof_Sd_low[flag_Sd]), alpha=alcalfill, color=cols_for_AV)

h1.set_ylim([0, 1])
h2.set_ylim([0, 1])
h3.set_ylim([0, 1])
h4.set_ylim([0, 1])
h5.set_ylim([0, 1])
h6.set_ylim([0, 1])
h7.set_ylim([0, 1])

h1.set_xlim([0, 3])
h2.set_xlim([0, 3])
h3.set_xlim([0, 3])
h4.set_xlim([0, 3])
h5.set_xlim([0, 3])
h6.set_xlim([0, 3])
h7.set_xlim([0, 3])

h1.text(0.3,0.85,'E',size=30)
h2.text(0.3,0.85,'S0',size=30)
h3.text(0.3,0.85,'Sa',size=30)
h4.text(0.3,0.85,'Sb',size=30)
h5.text(0.3,0.85,'Sbc',size=30)
h6.text(0.3,0.85,'Sc',size=30)
h7.text(0.3,0.85,'Sd',size=30)

h1.text(0.3,0.55,'N=%d' % cnt_E,size=20, color=cols_for_AV)
h2.text(0.3,0.55,'N=%d' % cnt_S0,size=20, color=cols_for_AV)
h3.text(0.3,0.55,'N=%d' % cnt_Sa,size=20, color=cols_for_AV)
h4.text(0.3,0.55,'N=%d' % cnt_Sb,size=20, color=cols_for_AV)
h5.text(0.3,0.55,'N=%d' % cnt_Sbc,size=20, color=cols_for_AV)
h6.text(0.3,0.55,'N=%d' % cnt_Sc,size=20, color=cols_for_AV)
h7.text(0.3,0.55,'N=%d' % cnt_Sd,size=20, color=cols_for_AV)

h1.text(0.3,0.75,r'$\beta_{V,0}$=%.2f' % min_bv0[0],size=20, color=cols_for_AV)
h2.text(0.3,0.75,r'$\beta_{V,0}$=%.2f' % min_bv0[1],size=20, color=cols_for_AV)
h3.text(0.3,0.75,r'$\beta_{V,0}$=%.2f' % min_bv0[2],size=20, color=cols_for_AV)
h4.text(0.3,0.75,r'$\beta_{V,0}$=%.2f' % min_bv0[3],size=20, color=cols_for_AV)
h5.text(0.3,0.75,r'$\beta_{V,0}$=%.2f' % min_bv0[4],size=20, color=cols_for_AV)
h6.text(0.3,0.75,r'$\beta_{V,0}$=%.2f' % min_bv0[5],size=20, color=cols_for_AV)
h7.text(0.3,0.75,r'$\beta_{V,0}$=%.2f' % min_bv0[6],size=20, color=cols_for_AV)

h1.text(0.3, 0.65, r'$\chi^2$=%.2f' % min_chi[0],size=20, color=cols_for_AV)
h2.text(0.3, 0.65, r'$\chi^2$=%.2f' % min_chi[1],size=20, color=cols_for_AV)
h3.text(0.3, 0.65, r'$\chi^2$=%.2f' % min_chi[2],size=20, color=cols_for_AV)
h4.text(0.3, 0.65, r'$\chi^2$=%.2f' % min_chi[3],size=20, color=cols_for_AV)
h5.text(0.3, 0.65, r'$\chi^2$=%.2f' % min_chi[4],size=20, color=cols_for_AV)
h6.text(0.3, 0.65, r'$\chi^2$=%.2f' % min_chi[5],size=20, color=cols_for_AV)
h7.text(0.3, 0.65, r'$\chi^2$=%.2f' % min_chi[6],size=20, color=cols_for_AV)

h4.set_ylabel(r'$A_{V}$',size=30)
h5.set_xlabel(r'$R/R_{50,V}^{P}$',size=30)
h6.set_xlabel(r'$R/R_{50,V}^{P}$',size=30)
h7.set_xlabel(r'$R/R_{50,V}^{P}$',size=30)
h4.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
h7.set_yticks([0, 0.2, 0.4, 0.6, 0.8])
h7.set_xticks([0, 1, 2, 3])
h5.set_xticks([1, 2, 3])
h6.set_xticks([1, 2, 3])
h1.tick_params(labelsize='20',direction='in',top=True,right=True)
h2.tick_params(direction='in',top=True,right=True,bottom=True,left=True)
h3.tick_params(direction='in',top=True,right=True,bottom=True,left=True)
h4.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h5.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h6.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)
h7.tick_params(labelsize='20',direction='in',top=True,right=True,bottom=True,left=True)

h7.legend(bbox_to_anchor=(1.05, 0.6), prop={'size': 25})

h1.set_xticklabels([])
h2.set_xticklabels([])
h2.set_yticklabels([])
h3.set_xticklabels([])
h3.set_yticklabels([])
h4.set_xticklabels([])
h5.set_yticklabels([])
h6.set_yticklabels([])

fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/dust_prof_T_avg_betav_median_rebin.pdf')
plt.close(fig)

