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

import matplotlib.pyplot as plt

htypes = ['E0','S0','S0/a','Sa','Sab','Sb','Sbc','Sc','Scd','Sd','Sdm','Sm','Im','?'] # Hubble types corresponding Nair+10 T-Type
agntypes = ['SF','transition/mixed AGN','Seyfert','LINER']	# Kauffmann+03 from Nair+10

work_dir = '/Users/dhk/work/data/NGC_IC/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

text_for_bv0s=['-0.05', '-0.04', '-0.03', '-0.02', '-0.01', '0', '0.01', '0.02', '0.03', '0.04']
cols_for_bv0s=['black', 'gray', 'silver', 'firebrick', 'sienna', 'sandybrown', 'tan', 'olivedrab', 'palegreen', 'navy']
cols_for_AV = 'maroon'
linest_for_CALIFA = ':'

hubble_bv0 = [0.56, 0.66, 0.68, 0.69, 0.87, 1.00, 1.16]
hubble_bv0_err = [0.14, 0.2, 0.28, 0.25, 0.21, 0.34, 0.34]

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

av_prof_E = np.zeros((46,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_E[:] = np.nan
av_prof_S0 = np.zeros((38,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_S0[:] = np.nan
av_prof_Sa = np.zeros((31,30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_Sa[:] = np.nan
av_prof_Sb = np.zeros((43,30))        # For stacking Dust profiles for SFH4 w/ 6 Zs
av_prof_Sb[:] = np.nan
av_prof_Sbc = np.zeros((16,30))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T < 0.5
av_prof_Sbc[:] = np.nan
av_prof_Sc = np.zeros((45,30))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T > 0.5
av_prof_Sc[:] = np.nan
av_prof_Sd = np.zeros((15,30))        # For stacking Dust profiles for SFH5 w/ 6 Zs
av_prof_Sd[:] = np.nan

cnt_E = 0	# counter for SFH3
cnt_S0 = 0	# counter for SFH3 with B/T < div_btot_sfh3
cnt_Sa = 0	# counter for SFH3 with B/T < div_btot_sfh3
cnt_Sb = 0	# counter for SFH4
cnt_Sbc = 0	# counter for SFH4 with B/T < div_btot_sfh4
cnt_Sc = 0	# counter for SFH4 with B/T > div_btot_sfh4
cnt_Sd = 0	# counter for SFH5

xnew = np.linspace(0, 3, num=30)

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


    if T < -3:
        bv0_zero = hubble_bv0[0]
        Av1 = 2.5 * np.log10(bv0_zero / y1)
        Av1[np.where(Av1 < 0)] = 0.0
        f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
        av_prof_E[cnt_E,:] = f1(xnew)
        cnt_E = cnt_E + 1
    if T < 0 and T >= -3:
        bv0_zero = hubble_bv0[1]
        Av1 = 2.5 * np.log10(bv0_zero / y1)
        Av1[np.where(Av1 < 0)] = 0.0
        f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
        av_prof_S0[cnt_S0, :] = f1(xnew)
        cnt_S0 = cnt_S0 + 1
    if T < 2 and T >= 0:
        bv0_zero = hubble_bv0[2]
        Av1 = 2.5 * np.log10(bv0_zero / y1)
        Av1[np.where(Av1 < 0)] = 0.0
        f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
        av_prof_Sa[cnt_Sa, :] = f1(xnew)
        cnt_Sa = cnt_Sa + 1
    if T < 4 and T >= 2:
        bv0_zero = hubble_bv0[3]
        Av1 = 2.5 * np.log10(bv0_zero / y1)
        Av1[np.where(Av1 < 0)] = 0.0
        f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
        av_prof_Sb[cnt_Sb, :] = f1(xnew)
        cnt_Sb = cnt_Sb + 1
    if T < 5 and T >= 4:
        bv0_zero = hubble_bv0[4]
        Av1 = 2.5 * np.log10(bv0_zero / y1)
        Av1[np.where(Av1 < 0)] = 0.0
        f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
        av_prof_Sbc[cnt_Sbc, :] = f1(xnew)
        cnt_Sbc = cnt_Sbc + 1
    if T < 7 and T >= 5:
        bv0_zero = hubble_bv0[5]
        Av1 = 2.5 * np.log10(bv0_zero / y1)
        Av1[np.where(Av1 < 0)] = 0.0
        f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
        av_prof_Sc[cnt_Sc, :] = f1(xnew)
        cnt_Sc = cnt_Sc + 1
    if T < 9 and T >= 7:
        bv0_zero = hubble_bv0[6]
        Av1 = 2.5 * np.log10(bv0_zero / y1)
        Av1[np.where(Av1 < 0)] = 0.0
        f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
        av_prof_Sd[cnt_Sd, :] = f1(xnew)
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
h1.errorbar(1, rosa_AV[0][12], yerr=rosa_AV[0][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h2.plot(np.arange(0.,3.,0.1),rosa_AV[2][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h2.errorbar(1, rosa_AV[2][12], yerr=rosa_AV[2][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h3.plot(np.arange(0.,3.,0.1),rosa_AV[4][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h3.errorbar(1, rosa_AV[4][12], yerr=rosa_AV[4][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h4.plot(np.arange(0.,3.,0.1),rosa_AV[6][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h4.errorbar(1, rosa_AV[6][12], yerr=rosa_AV[6][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h5.plot(np.arange(0.,3.,0.1),rosa_AV[8][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h5.errorbar(1, rosa_AV[8][12], yerr=rosa_AV[8][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h6.plot(np.arange(0.,3.,0.1),rosa_AV[10][3:],color='black',alpha=alcalline, linewidth=3, linestyle=linest_for_CALIFA)
h6.errorbar(1, rosa_AV[10][12], yerr=rosa_AV[10][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

h7.plot(np.arange(0.,3.,0.1),rosa_AV[12][3:],color='black',alpha=alcalline, linewidth=3, label='GD15 (CALIFA)', linestyle=linest_for_CALIFA)
h7.errorbar(1, rosa_AV[12][12], yerr=rosa_AV[12][2], color='black', alpha=alcalline, linewidth=2, capsize=3)

# h1.plot(np.arange(0., 3., 0.1), np.nanmedian(av_prof_E, axis=0),  linewidth=3)
# h2.plot(np.arange(0., 3., 0.1), np.nanmedian(av_prof_S0, axis=0), linewidth=3)
# h3.plot(np.arange(0., 3., 0.1), np.nanmedian(av_prof_Sa, axis=0), linewidth=3)
# h4.plot(np.arange(0., 3., 0.1), np.nanmedian(av_prof_Sb, axis=0), linewidth=3)
# h5.plot(np.arange(0., 3., 0.1), np.nanmedian(av_prof_Sbc, axis=0), linewidth=3)
# h6.plot(np.arange(0., 3., 0.1), np.nanmedian(av_prof_Sc, axis=0), linewidth=3)
# h7.plot(np.arange(0., 3., 0.1), np.nanmedian(av_prof_Sd, axis=0), linewidth=3, label=r'$\beta_{V}$ method (median)')


h1.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_E, axis=0), color=cols_for_AV, linewidth=3)
h2.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_S0, axis=0), color=cols_for_AV, linewidth=3)
h3.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sa, axis=0), color=cols_for_AV, linewidth=3)
h4.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sb, axis=0), color=cols_for_AV, linewidth=3)
h5.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sbc, axis=0), color=cols_for_AV, linewidth=3)
h6.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sc, axis=0), color=cols_for_AV, linewidth=3)
h7.plot(np.arange(0., 3., 0.1), np.nanmean(av_prof_Sd, axis=0), color=cols_for_AV, linewidth=3, label=r'$\beta_{V}$ method')


h1.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_E, axis=0) - np.nanstd(av_prof_E, axis=0),
    np.nanmean(av_prof_E, axis=0) + np.nanstd(av_prof_E, axis=0), alpha=alcalfill, color=cols_for_AV)
h2.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_S0, axis=0) - np.nanstd(av_prof_S0, axis=0),
    np.nanmean(av_prof_S0, axis=0) + np.nanstd(av_prof_S0, axis=0), alpha=alcalfill, color=cols_for_AV)
h3.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sa, axis=0) - np.nanstd(av_prof_Sa, axis=0),
    np.nanmean(av_prof_Sa, axis=0) + np.nanstd(av_prof_Sa, axis=0), alpha=alcalfill, color=cols_for_AV)
h4.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sb, axis=0) - np.nanstd(av_prof_Sb, axis=0),
    np.nanmean(av_prof_Sb, axis=0) + np.nanstd(av_prof_Sb, axis=0), alpha=alcalfill, color=cols_for_AV)
h5.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sbc, axis=0) - np.nanstd(av_prof_Sbc, axis=0),
    np.nanmean(av_prof_Sbc, axis=0) + np.nanstd(av_prof_Sbc, axis=0), alpha=alcalfill, color=cols_for_AV)
h6.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sc, axis=0) - np.nanstd(av_prof_Sc, axis=0),
    np.nanmean(av_prof_Sc, axis=0) + np.nanstd(av_prof_Sc, axis=0), alpha=alcalfill, color=cols_for_AV)
h7.fill_between(np.arange(0, 3, 0.1), np.nanmean(av_prof_Sd, axis=0) - np.nanstd(av_prof_Sd, axis=0),
    np.nanmean(av_prof_Sd, axis=0) + np.nanstd(av_prof_Sd, axis=0), alpha=alcalfill, color=cols_for_AV)


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

h1.text(0.3,0.75,r'$\beta_{V,0}$=0.56',size=20, color=cols_for_AV)
h2.text(0.3,0.75,r'$\beta_{V,0}$=0.66',size=20, color=cols_for_AV)
h3.text(0.3,0.75,r'$\beta_{V,0}$=0.68',size=20, color=cols_for_AV)
h4.text(0.3,0.75,r'$\beta_{V,0}$=0.69',size=20, color=cols_for_AV)
h5.text(0.3,0.75,r'$\beta_{V,0}$=0.87',size=20, color=cols_for_AV)
h6.text(0.3,0.75,r'$\beta_{V,0}$=1.00',size=20, color=cols_for_AV)
h7.text(0.3,0.75,r'$\beta_{V,0}$=1.16',size=20, color=cols_for_AV)

h4.set_ylabel(r'$A_{V}$',size=30)
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

fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/dust_prof_T_indivi.pdf')
plt.close(fig)

