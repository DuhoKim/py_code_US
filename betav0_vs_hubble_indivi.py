#####################################################
# Python script that reads betaV values and T-types
# of 257 NGC/IC galaxies and plot
# written by Duho Kim (06/21/18)        
######################################################
import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Ellipse
from matplotlib.lines import Line2D
from matplotlib import colors
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd
from operator import itemgetter
import os.path

blank_flag = True

file=r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)

bv0=[	[1.9409089,1.9044532,1.3527486,1.1198820,0.82968015,0.58503551],	# BetaVzero values for SFH3,4,5 and Z=.0001,.0004,.004,.008,.02,.05
	[1.9860456,1.9576220,1.4390702,1.2023316,0.91737698,0.65453906],
	[2.3880801,2.3478914,1.6838646,1.4124115,1.1048444,0.77272439]]

# BetaVzero values for SFH3,4,5 and Z=.0001,.0004,.004,.008,.02,.05 w/ no Z evolution
bv0_noevol = [[1.9408058, 1.8461447, 1.1066952, 0.80323478, 0.59881037, 0.38305162],
        [1.9878865, 1.8940383, 1.1417368, 0.84291972, 0.63275795, 0.40320143],
        [2.3858643, 2.2975972, 1.4978560, 1.1777806, 0.93497056, 0.59137094]]

Zs = [	['0.13', '0.34', '0.85'],		# mass-weighted Z for SFH3,4,5 w/ Z0 = 0.05
		  ['0.11', '0.28', '0.69'],
		  ['0.15', '0.37', '0.93']]

Zs_noevol = ['0.4', '', '2.5' ]		# Z

xoff = [[-0.1, -0.1, -0.1],
		[0.25, 0.25, 0.25],
		[-0.35, -0.35, -0.35]]

xoff_noevol = [	[0, 0.2, 0],
				[0.4, 0.5, 0.4],
				[0, 0.2, 0]]

yoff = [[0.01, 0.01, -0.05],
		[-0.05, 0.01, 0.01],
		[-0.05, 0.01, 0.01]]

yoff_noevol = [	[-0.05, 0.01, 0.01],
				[0.01, -0.05, 0.01],
				[0.01, 0.01, 0.01]]

hubble_bv0_med=[0.55, 0.64, 0.72, 0.68, 0.89, 0.96, 1.14]
hubble_bv0_med_low=[0.48, 0.52, 0.45, 0.49, 0.69, 0.70, 0.98]
hubble_bv0_med_high=[0.63, 0.77, 0.83, 0.87, 1.04, 1.33, 1.43]

hubble_bv0 = [0.56, 0.66, 0.68, 0.69, 0.87, 1.00, 1.16]
hubble_bv0_err = [0.14, 0.20, 0.28, 0.25, 0.21, 0.34, 0.34]
hubble_bv0_bh = [0.54, 0.66, 0.65, 0.66]
hubble_bv0_bh_err = [0.04, 0.11, 0.2, 0.16]
hubble_bv0_bl = [0.57, 0.65, 0.72, 0.78]
hubble_bv0_bl_err = [0.15, 0.28, 0.35, 0.40]
hubble_bv0_med_comb = [0.56, 0.66, 0.74, 0.77, 0.94, 0.97, 1.20]
hubble_bv0_bh_orig = [0.51, 0.58, 0.68, 0.69]
hubble_bv0_bl_orig = [0.5, 0.64, 0.76, 0.79]

hubble_label = ('E', 'S0', 'Sa' ,'Sb', 'Sbc', 'Sc', 'Sd')
hubble_coords = range(7)

sfh_ttype_bgn=[0, 2, 5]
sfh_ttype_end=[1, 4, 6]
cols_for_Z=['black','gray','violet','blue','green','darkorange']
linest_for_SFH=['-', '--', ':']
#cols_for_Z_noevol=['black','green','red']
text_for_Z=[r'Z=0.0001','Z=0.0004','Z=0.004','Z=0.008','Z=0.02 (Z$_{\odot}$)','Z=0.05']

agn=['Not active','Starburst','LIRG','LINER','Sy 2','Sy 1']
cm=['black','blue','red','darkorange','darkgreen','darkturquoise']

bt=np.load('/Users/dhk/work/py/pyraf/betav_ttype.npy')
btot=ascii.read("/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask/btot_3rd.txt")
sha2fig=np.load('/Users/dhk/work/data/NGC_IC/pdf/sha2fig.npy')
sha_cat=ascii.read('/Users/dhk/work/cat/NGC_IC/sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv') # Sample NGC/IC numbers

########################################  BetaV0 vs. hubble ########################################

fig = plt.figure(figsize=(5,5))

a = plt.subplot(111)

#plt.fill_between(hubble_coords, np.array([0.65, 0.64, 0.65, 0.72, 0.77, 1.06, 1.39]) - np.array([0.07, 0.07, 0.03, 0.09, 0.12, 0.31, 0.34]),
#				 np.array([0.65, 0.64, 0.65, 0.72, 0.77, 1.06, 1.39]) + np.array([0.03, 0.03, 0.14, 0.06, 0.26, 0.3, 0.43]), alpha=0.1, color='k', label='')


for i in range(3):		# for SFH3,4,5
	for j in range(3,6):	# for Z=.0001,.0004,.004,.008,.02,.05
		plt.plot([sfh_ttype_bgn[i], sfh_ttype_end[i]],[bv0[i][j],bv0[i][j]],'-',color=cols_for_Z[j],alpha=1.0, linestyle=linest_for_SFH[i])
		plt.plot([sfh_ttype_bgn[i], sfh_ttype_end[i]], [bv0_noevol[i][j], bv0_noevol[i][j]], '-', alpha=0.3, linestyle=':', color='k')
		plt.text(sfh_ttype_bgn[i] + xoff[i][j-3], bv0[i][j] + yoff[i][j-3], r'$\langle$Z$\rangle$$_{\mathrm{M}}$='+Zs[i][j-3]+r'$\,$Z$\odot$', size=10, color=cols_for_Z[j], alpha=1.0)
		plt.text(sfh_ttype_bgn[i] + xoff_noevol[i][j-3] , bv0_noevol[i][j] + yoff_noevol[i][j-3], r'Z='+Zs_noevol[j-3]+r'$\,$Z$\odot$', size=10, alpha=0.3, color='k')


#plt.scatter(hubble_coords[:4], hubble_bv0_bl, color='blue', label='Small bulge')
#plt.scatter(hubble_coords[:4], hubble_bv0_bh, color='red', label='Large bulge')

if not blank_flag:
	plt.fill_between(hubble_coords, hubble_bv0_med_low, hubble_bv0_med_high, alpha=0.3, color='seagreen')
	plt.plot(hubble_coords, hubble_bv0_med, color='darkgreen', linestyle='-.', alpha=0.6, linewidth=5.0)
	plt.plot(hubble_coords, hubble_bv0_med_comb, color='maroon', linestyle='--', alpha=0.6, linewidth=5.0)

# plt.errorbar(hubble_coords, hubble_bv0, yerr=hubble_bv0_err, fmt='o', color='black', label=r'All (Individual $A_{V,r}$ fitting) ', )
# plt.scatter(hubble_coords, hubble_bv0_av, color='k', label=r'All ($\langle A_{V,r} \rangle$ fitting)', marker='s')
# plt.errorbar([-0.2, 0.8, 1.8, 2.8], hubble_bv0_bh, yerr=hubble_bv0_bh_err, fmt='o', color='red', label=r'Large bulge ($\langle A_{V,r} \rangle$ in square)')
# plt.errorbar([-0.2, 0.8, 1.8, 2.8], hubble_bv0_bh_orig, yerr=[0, 0, 0, 0], fmt='s', color='red')
# plt.errorbar([0.2, 1.2, 2.2, 3.2], hubble_bv0_bl, yerr=hubble_bv0_bl_err, fmt='o', color='blue', label=r'Small bulge ($\langle A_{V,r} \rangle$ in square)')
# plt.errorbar([0.2, 1.2, 2.2, 3.2], hubble_bv0_bl_orig, yerr=[0, 0, 0, 0], fmt='s', color='blue')


# plt.legend(fontsize='10', frameon='False')

#plt.ylim(0,3.5)
plt.xlabel('Hubble type',size=15)
plt.ylabel(r'$\beta_{V,0,\mathrm{z\simeq0}}$',size=15)
plt.tick_params(direction='in',labelsize='large',top=True,right=True)
#plt.xticks([-5,0,5,10])
plt.setp(a, xticks=range(7),xticklabels=hubble_label)
fig.savefig('/Users/dhk/Documents/publish/ngcic_rev/betav0_hubble_blank.pdf')
