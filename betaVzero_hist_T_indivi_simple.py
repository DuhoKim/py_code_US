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

hubble_bv0 = [0.56, 0.66, 0.74, 0.77, 0.94, 0.97, 1.20]
hubble_bv0_bh = [0.55, 0.61, 0.70, 0.73]
hubble_bv0_bl = [0.56, 0.75, 0.76, 0.81]
hubble_btot = [0.7, 0.5, 0.45, 0.25]

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

bv0_step = np.linspace(-0.5, 1.0, 161)

bv0s = np.zeros((7, 50))	# bv0 individual
bv0_tot = np.zeros(257)	# bv0 individual
chis = np.zeros((7, 50))	# chi individual
bv0s[:] = np.nan
bv0_tot[:] = np.nan
chis[:] = np.nan
cnts = np.zeros(7)

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

    if T < -3:
        bv0_zero = hubble_bv0[0]
        califa = rosa_AV[0][3:]
        ttype_no = 0
    if T < 0 and T >= -3:
        bv0_zero = hubble_bv0[1]
        califa = rosa_AV[2][3:]
        ttype_no = 1
    if T < 2 and T >= 0:
        bv0_zero = hubble_bv0[2]
        califa = rosa_AV[4][3:]
        ttype_no = 2
    if T < 4 and T >= 2:
        bv0_zero = hubble_bv0[3]
        califa = rosa_AV[6][3:]
        ttype_no = 3
    if T < 5 and T >= 4:
        bv0_zero = hubble_bv0[4]
        califa = rosa_AV[8][3:]
        ttype_no = 4
    if T < 7 and T >= 5:
        bv0_zero = hubble_bv0[5]
        califa = rosa_AV[10][3:]
        ttype_no = 5
    if T < 9 and T >= 7:
        bv0_zero = hubble_bv0[6]
        califa = rosa_AV[12][3:]
        ttype_no = 6

    cnts[ttype_no] = cnts[ttype_no] + 1

    min_chi = 1e6
    for i in range(0, len(bv0_step)):
        Av1 = 2.5 * np.log10( (bv0_zero + bv0_step[i]) / y1)
        Av1_low = 2.5 * np.log10((bv0_zero + bv0_step[i]) / (y1 + s1))
        Av1_hi = 2.5 * np.log10((bv0_zero + bv0_step[i]) / (y1 - s1))
        Av1[np.where(Av1 < 0)] = 0.0
        #Av1_low[np.where(Av1_low < 0)] = 0.0
        #Av1_hi[np.where(Av1_hi < 0)] = 0.0
        f1 = interp1d(x1_hlr, Av1, fill_value=np.nan, bounds_error=False)
        f1_low = interp1d(x1_hlr, Av1_low, fill_value=np.nan, bounds_error=False)
        f1_hi = interp1d(x1_hlr, Av1_hi, fill_value=np.nan, bounds_error=False)
        av_prof = f1(xnew)
        av_low = f1_low(xnew)
        av_hi = f1_hi(xnew)
        av_err = (av_hi - av_low) / 2
        chi = 0
        rn = 0
        for j in range(0, 15):
            if ~np.isnan(av_prof[j]) and ~np.isnan(av_err[j]):
                chi += (av_prof[j] - califa[j*2])**2 / (av_err[j]**2)
                rn = rn + 1
        if chi < min_chi:
            min_chi = chi
            bv0s[ttype_no, int(cnts[ttype_no]) - 1] = bv0_zero + bv0_step[i]
            chis[ttype_no, int(cnts[ttype_no]) - 1] = chi / (rn - 1)
            bv0_tot[x] = bv0_zero + bv0_step[i]

fig = plt.figure(figsize=(10,12))

h1 = fig.add_axes([0.1,0.65,0.275,0.275])
h2 = fig.add_axes([0.375,0.65,0.275,0.275])
h3 = fig.add_axes([0.65,0.65,0.275,0.275])
h4 = fig.add_axes([0.1,0.375,0.275,0.275])
h5 = fig.add_axes([0.375,0.375,0.275,0.275])
h6 = fig.add_axes([0.65,0.375,0.275,0.275])
h7 = fig.add_axes([0.1,0.1,0.275,0.275])

hub_types = ['E', 'S0', 'Sa', 'Sb', 'Sbc', 'Sc', 'Sd']
bins = np.arange(0.3, 2.0, 0.1)
for i, ax in enumerate(fig.axes):
    bv0s_notnan = bv0s[i, :int(cnts[i])]
    bv0s_notnan = bv0s_notnan[~np.isnan(bv0s_notnan)]
    ax.hist(bv0s_notnan, bins=bins, label=r'Individual A$_V$-profile fitting')

    ax.set_xlim(0.2, 2.0)
    ax.set_ylim(0, 25)
    ax.set_xticks([0.5, 1, 1.5, 2])
    ax.set_xticks(np.arange(0.3, 2.0, 0.1), minor=True)
    ax.set_yticks(np.arange(0, 25, 1), minor=True)
    if i == 4:
        ax.text(0.7, 0.85, hub_types[i], size=30, transform=ax.transAxes)
    else:
        ax.text(0.78, 0.85, hub_types[i], size=30, transform=ax.transAxes)
    ax.plot([np.percentile(bv0s_notnan, 16.5), np.percentile(bv0s_notnan, 16.5)], [0, 25], linestyle=':', linewidth=3.0,
            color='seagreen', label=r'1$\sigma$ range of individual fit')
    ax.plot([np.percentile(bv0s_notnan, 50), np.percentile(bv0s_notnan, 50)], [0, 25], linestyle='-.', linewidth=3.0,
            color='darkgreen', label='median of individual fit')
    ax.plot([np.percentile(bv0s_notnan, 83.5), np.percentile(bv0s_notnan, 83.5)], [0, 25], linestyle=':', linewidth=3.0,
            color='seagreen')

    ax.plot([hubble_bv0[i], hubble_bv0[i]], [0, 25], color='maroon', linestyle='--', linewidth=3.0,
            label=r'$\langle\beta_V\rangle$-profile fitting')
    if i == 6:
        ax.legend(bbox_to_anchor=(1.05, 0.75), prop={'size': 25}, frameon=False)
    print("%i, {%4.2f}, {%4.2f}, {%4.2f}" % (i, np.percentile(bv0s_notnan, 16.5),
                                             np.percentile(bv0s_notnan, 50),
                                             np.percentile(bv0s_notnan, 83.5)))

    ax.tick_params(labelsize='20', direction='in', top=True, right=True, bottom=True, left=True, length=5, width=1)
    ax.tick_params(which='minor', direction='in', top=True, bottom=True, right=True)

h4.set_ylabel('N',size=30)
h5.set_xlabel(r'$\beta_{V,0}$',size=30)
h6.set_xlabel(r'$\beta_{V,0}$',size=30)
h7.set_xlabel(r'$\beta_{V,0}}$',size=30)

#h7.legend(bbox_to_anchor=(1.05, 0.5), prop={'size': 20})

h1.set_xticklabels([])
h2.set_xticklabels([])
h2.set_yticklabels([])
h3.set_xticklabels([])
h3.set_yticklabels([])
h4.set_xticklabels([])
h4.set_yticklabels([0, 5, 10, 15, 20])
#h5.set_xticklabels(['', '', '', 0.5, '', '', '', '', 1, '', '', '', '', 1.5], minor=True)
h5.set_yticklabels([])
#h6.set_xticklabels([0.2, '', '', '', '', '', '', '', 1, '', '', '', '', 1.5, '', '', '', '', 2.0], minor=True)
h6.set_yticklabels([])
#h7.set_xticklabels([0.2, '', '', '', '', '', '', '', 1, '', '', '', '', 1.5, '', '', '', '', 2.0], minor=True)
h7.set_yticklabels([0, 5, 10, 15, 20])

fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/betaVzero_hist_rebin.pdf')
plt.close(fig)

np.save('/Users/dhk/Documents/publish/ngcic_rev/bvzero_indi_rebin.npy', bv0_tot)

