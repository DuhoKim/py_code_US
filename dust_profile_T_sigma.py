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

hubble_bv0 = [0.5, 0.6, 0.73, 0.73, 0.81, 0.85, 0.97]
hubble_bv0_bh = [0.513, 0.58, 0.68, 0.69]
hubble_bv0_bl = [0.507, 0.64, 0.76, 0.79]

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

bv0_step = np.linspace(-0.4, 0.4, 161)

av_prof_E = np.zeros((46, len(bv0_step), 30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_E[:] = np.nan
av_prof_S0 = np.zeros((38, len(bv0_step), 30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_S0[:] = np.nan
av_prof_Sa = np.zeros((31, len(bv0_step), 30))	# For stacking Dust profiles for SFH3 w/ 6 Zs
av_prof_Sa[:] = np.nan
av_prof_Sb = np.zeros((43, len(bv0_step), 30))        # For stacking Dust profiles for SFH4 w/ 6 Zs
av_prof_Sb[:] = np.nan
av_prof_Sbc = np.zeros((16, len(bv0_step), 30))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T < 0.5
av_prof_Sbc[:] = np.nan
av_prof_Sc = np.zeros((45, len(bv0_step), 30))        # For stacking Dust profiles for SFH4 w/ 6 Zs with B/T > 0.5
av_prof_Sc[:] = np.nan
av_prof_Sd = np.zeros((15, len(bv0_step), 30))        # For stacking Dust profiles for SFH5 w/ 6 Zs
av_prof_Sd[:] = np.nan

cnt_E = 0	# counter for SFH3
cnt_S0 = 0	# counter for SFH3 with B/T < div_btot_sfh3
cnt_Sa = 0	# counter for SFH3 with B/T < div_btot_sfh3
cnt_Sb = 0	# counter for SFH4
cnt_Sbc = 0	# counter for SFH4 with B/T < div_btot_sfh4
cnt_Sc = 0	# counter for SFH4 with B/T > div_btot_sfh4
cnt_Sd = 0	# counter for SFH5

unc_fac = 1.3

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

    if T < -3:
        for i in range(0, len(bv0_step)):
            Av1 = 2.5 * np.log10( (hubble_bv0[0] + bv0_step[i]) / y1 )
            Av1[np.where(Av1 < 0)] = 0.0
            f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
            av_prof_E[cnt_E, i, :] = f1(xnew)
        cnt_E = cnt_E + 1
    if T < 0 and T >= -3:
        for i in range(0, len(bv0_step)):
            Av1 = 2.5 * np.log10( (hubble_bv0[1] + bv0_step[i]) / y1)
            Av1[np.where(Av1 < 0)] = 0.0
            f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
            av_prof_S0[cnt_S0, i, :] = f1(xnew)
        cnt_S0 = cnt_S0 + 1
    if T < 2 and T >= 0:
        for i in range(0, len(bv0_step)):
            Av1 = 2.5 * np.log10( (hubble_bv0[2] + bv0_step[i]) / y1)
            Av1[np.where(Av1 < 0)] = 0.0
            f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
            av_prof_Sa[cnt_Sa, i, :] = f1(xnew)
        cnt_Sa = cnt_Sa + 1
    if T < 4 and T >= 2:
        for i in range(0, len(bv0_step)):
            Av1 = 2.5 * np.log10( (hubble_bv0[3] + bv0_step[i]) / y1)
            Av1[np.where(Av1 < 0)] = 0.0
            f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
            av_prof_Sb[cnt_Sb, i, :] = f1(xnew)
        cnt_Sb = cnt_Sb + 1
    if T < 5 and T >= 4:
        for i in range(0, len(bv0_step)):
            Av1 = 2.5 * np.log10( (hubble_bv0[4] + bv0_step[i]) / y1)
            Av1[np.where(Av1 < 0)] = 0.0
            f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
            av_prof_Sbc[cnt_Sbc, i, :] = f1(xnew)
        cnt_Sbc = cnt_Sbc + 1
    if T < 7 and T >= 5:
        for i in range(0, len(bv0_step)):
            Av1 = 2.5 * np.log10( (hubble_bv0[5] + bv0_step[i]) / y1)
            Av1[np.where(Av1 < 0)] = 0.0
            f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
            av_prof_Sc[cnt_Sc, i, :] = f1(xnew)
        cnt_Sc = cnt_Sc + 1
    if T < 9 and T >= 7:
        for i in range(0, len(bv0_step)):
            Av1 = 2.5 * np.log10( (hubble_bv0[6] + bv0_step[i]) / y1)
            Av1[np.where(Av1 < 0)] = 0.0
            f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
            av_prof_Sd[cnt_Sd, i, :] = f1(xnew)
        cnt_Sd = cnt_Sd + 1

fig = plt.figure(figsize=(10,12))

h1 = fig.add_axes([0.1,0.65,0.275,0.275])
h2 = fig.add_axes([0.375,0.65,0.275,0.275])
h3 = fig.add_axes([0.65,0.65,0.275,0.275])
h4 = fig.add_axes([0.1,0.375,0.275,0.275])
h5 = fig.add_axes([0.375,0.375,0.275,0.275])
h6 = fig.add_axes([0.65,0.375,0.275,0.275])
h7 = fig.add_axes([0.1,0.1,0.275,0.275])

all_chi = np.zeros((7, len(bv0_step)))
min_chi = np.zeros(7)
min_chi[:] = 1e6
for i in range(0, len(bv0_step)):
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
            chi += (prof_E[j] - rosa_AV[0][3+j])**2 / (prof_E_err[j]**2 + rosa_AV[0][2]**2)
    all_chi[0, i] = chi

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_S0[j]) and ~np.isnan(prof_S0_err[j]):
            chi += (prof_S0[j] - rosa_AV[2][3 + j])**2 / (prof_S0_err[j]**2 + rosa_AV[2][2]**2)
    all_chi[1, i] = chi

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sa[j]) and ~np.isnan(prof_Sa_err[j]):
            chi += (prof_Sa[j] - rosa_AV[4][3 + j])**2 / (prof_Sa_err[j]**2 + rosa_AV[4][2]**2)
    all_chi[2, i] = chi

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sb[j]) and ~np.isnan(prof_Sb_err[j]):
            chi += (prof_Sb[j] - rosa_AV[6][3 + j])**2 / (prof_Sb_err[j]**2 + rosa_AV[6][2]**2)
    all_chi[3, i] = chi

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sbc[j]) and ~np.isnan(prof_Sbc_err[j]):
            chi += (prof_Sbc[j] - rosa_AV[8][3 + j])**2 / (prof_Sbc_err[j]**2 + rosa_AV[8][2]**2)
    all_chi[4, i] = chi

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sc[j]) and ~np.isnan(prof_Sc_err[j]):
            chi += (prof_Sc[j] - rosa_AV[10][3 + j])**2 / (prof_Sc_err[j]**2 + rosa_AV[10][2]**2)
    all_chi[5, i] = chi

    chi = 0
    for j in range(0, 30):
        if ~np.isnan(prof_Sd[j]) and ~np.isnan(prof_Sd_err[j]):
            chi += (prof_Sd[j] - rosa_AV[12][3 + j])**2 / (prof_Sd_err[j]**2 + rosa_AV[12][2]**2)
    all_chi[6, i] = chi

min_chi = np.zeros(7)
min_idx = range(7)
low_idx = range(7)
high_idx = range(7)
for i in range(0, 7):
    min_chi[i] = min(all_chi[i, :])
    min_idx[i] = (np.abs(all_chi[i, :] - min_chi[i])).argmin()
    low_idx[i] = (np.abs(all_chi[i, : min_idx[i]] - min_chi[i] * unc_fac)).argmin()
    high_idx[i] = (np.abs(all_chi[i, min_idx[i] :] - min_chi[i] * unc_fac)).argmin()

alsfhline=0.1	# alpha values for model Kim+17
alsfhfill=0.1
linest=[':','-.','--']

hub_types = ['E', 'S0', 'Sa', 'Sb', 'Sbc', 'Sc', 'Sd']
hub_cnt = [cnt_E, cnt_S0, cnt_Sa, cnt_Sb, cnt_Sbc, cnt_Sc, cnt_Sd]
for i, ax in enumerate(fig.axes):
    ax.plot(bv0_step + hubble_bv0[i], all_chi[i, :]/29, color=cols_for_AV, linewidth=3)
    ax.set_xlim(0.4, 1.5)
    ax.set_xticks(np.arange(0.4, 1.5, 0.1), minor=True)
    ax.set_xticks([0.5,  1.0, 1.5])
    ax.set_ylim(0, 1)
    #ax.set_xscale('log')
    ax.text(0.1, 0.85, hub_types[i], size=30, transform=ax.transAxes)
    ax.text(0.1, 0.65, 'N=%d' % hub_cnt[i], size=20, color=cols_for_AV, transform=ax.transAxes)
    ax.text(0.1, 0.75, r'$\beta_{V,0}=%4.2f^{+%4.2f}_{-%4.2f}$' % (bv0_step[min_idx[i]] + hubble_bv0[i],
                                              bv0_step[min_idx[i]+high_idx[i]] - bv0_step[min_idx[i]],
                                              bv0_step[min_idx[i]] - bv0_step[low_idx[i]]),
            size=20, color=cols_for_AV, transform=ax.transAxes)
    print("%i, {%4.2f}, {%4.2f}, {%4.2f}" % (i, bv0_step[min_idx[i]] + hubble_bv0[i],
                                              bv0_step[min_idx[i]+high_idx[i]] - bv0_step[min_idx[i]],
                                              bv0_step[min_idx[i]] - bv0_step[low_idx[i]]))
    ax.plot([hubble_bv0[i] + bv0_step[low_idx[i]], hubble_bv0[i] + bv0_step[min_idx[i]+high_idx[i]]],
            [min_chi[i] * unc_fac /29, min_chi[i] * unc_fac /29], 'o', markersize=6, color='blue')
    ax.plot(hubble_bv0[i] + bv0_step[min_idx[i]], min_chi[i]/29, 'o', markersize=6, color='red')
    ax.tick_params(labelsize='20', direction='in', top=True, right=True, bottom=True, left=True)
    ax.tick_params(which='minor', direction='in', top=True, bottom=True)

h4.set_ylabel(r'$\chi^{2}/(r_{n}-1)$',size=30)
h5.set_xlabel(r'$\beta_{V,0}$',size=30)
h6.set_xlabel(r'$\beta_{V,0}$',size=30)
h7.set_xlabel(r'$\beta_{V,0}}$',size=30)

#h7.legend(bbox_to_anchor=(1.05, 0.5), prop={'size': 30})

h1.set_xticklabels([])
h2.set_xticklabels([])
h2.set_yticklabels([])
h3.set_xticklabels([])
h3.set_yticklabels([])
h4.set_xticklabels([])
h4.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8])
h5.set_xticklabels([0.5, 1.0, ''])
h5.set_yticklabels([])
#h6.set_xticklabels([0.5,'', 0.7,'','', 1.0,''])
h6.set_yticklabels([])
#h7.set_xticklabels([0.5,'', 0.7,'','', 1.0,''])
h7.set_yticklabels([0.0, 0.2, 0.4, 0.6, 0.8])

fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/betaVzero_chi.pdf')
plt.close(fig)

