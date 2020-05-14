##########################################################
# Python script that compare Petrosian Half-light radius
# from SDSS and our measurement, and save combined radius
# written by Duho Kim (08/23/18)
##########################################################
from astropy.io import ascii
from astropy.table import Table, vstack
from astropy import units as u
from astropy.coordinates import SkyCoord
import numpy as np
import os
from shutil import copyfile
from pandas import DataFrame, read_csv
import pandas as pd
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from operator import itemgetter

hist_flag = False
alpha_type = 0.7

work_dir = '/Users/dhk/work/data/NGC_IC/'

ned1 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot = vstack([ned1, ned2, ned3, ned4, ned5])

s257 = ascii.read("/Users/dhk/work/cat/NGC_IC/sha_quarry_batch_257_without_prefix.txt")

rc3 = ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")
file = r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype = pd.read_excel(file)
file = r'/Users/dhk/work/cat/NGC_IC/SDSS/sdss_query_result_petroMag_R50.xls'
sdss = pd.read_excel(file)
file = r'/Users/dhk/work/cat/NGC_IC/Steinicke/NI2018.xls'
df = pd.read_excel(file)

hlr_r_pix = np.load('/Users/dhk/work/py/pyraf/petroR50_r.npy')
pmag_r = np.load('/Users/dhk/work/py/pyraf/petroM50_r.npy')

pmag_r = 22.5 - 2.5 * np.log10(pmag_r)


z = np.zeros(len(s257))
V = np.zeros(len(s257))
T = np.zeros(len(s257))

reff_g_scratch = np.zeros(len(s257))
pmag_g_scratch = np.zeros(len(s257))
pmag_r_scratch = np.zeros(len(s257))
reff_r_scratch = np.zeros(len(s257))
reff_V_scratch = np.zeros(len(s257))
reff_g_sdss = np.zeros(len(s257))
reff_r_sdss = np.zeros(len(s257))
reff_V_sdss = np.zeros(len(s257))
pmag_g_sdss = np.zeros(len(s257))
pmag_r_sdss = np.zeros(len(s257))

hubble = np.zeros(9)	# E S0 Sa Sb Sbc Sc Sd Irr Pec
hubble_label = ('E', 'S0', 'Sa' ,'Sb', 'Sbc', 'Sc', 'Sd', 'Irr', 'Pec')

i = 0
for x in range(0, len(ned_tot)):

    if ned_tot['col2'][x] in s257['id']:

        name = ned_tot['col2'][x]
        #reff_name.append(name)

        if name[0] == 'n':
            galnum = name[3:].strip()
            ngc_match = df.loc[(df['N'] == 'N') & (df['NI'] == int(galnum))]

        elif name[0] == 'i':
            galnum = name[2:].strip()
            ngc_match = df.loc[(df['N'] == 'I') & (df['NI'] == int(galnum))]

        z[i] = ned_tot['col8'][x]
        T[i] = ttype.iloc[i, 5]

        ########## read Reff from SDSS ##############
        sdss_match = sdss.loc[sdss['name'] == name]

        if len(sdss_match):
            g = sdss_match.iloc[0, 9]
            r = sdss_match.iloc[0, 10]
            sdss_V = g - 0.59 * (g - r) - 0.01

            pmag_g_sdss[i] = g
            pmag_r_sdss[i] = r

            if np.abs(ngc_match.iloc[0, 16] - sdss_V) < 1:
                reff_V_sdss[i] = (sdss_match.iloc[0, 13] * 2 + sdss_match.iloc[0, 14]) / 3.0
                reff_g_sdss[i] = sdss_match.iloc[0, 13]
                reff_r_sdss[i] = sdss_match.iloc[0, 14]

        ########## SHA match #############
        sha_match_idx = np.where(s257['id'] == name)
        reff_r_scratch[i] = hlr_r_pix[sha_match_idx[0]] * 0.6       # pixel to arcsec
        pmag_r_scratch[i] = pmag_r[sha_match_idx[0]]  # petrosian magnitude we measured

        i = i + 1

idx_E = np.where(np.logical_and.reduce((reff_r_sdss > 0.01, T < 0)))
idx_S = np.where(np.logical_and.reduce((reff_r_sdss > 0.01, T >= 0, T < 5)))
idx_L = np.where(np.logical_and.reduce((reff_r_sdss > 0.01, T >= 5)))

plt.figure(figsize=(5, 5))
plt.scatter(reff_r_scratch[idx_E], reff_r_sdss[idx_E], s=5, color='red', alpha=alpha_type, label='E, S0')
plt.scatter(reff_r_scratch[idx_S], reff_r_sdss[idx_S], s=5, color='green', alpha=alpha_type, label='Sa, Sb, Sbc')
plt.scatter(reff_r_scratch[idx_L], reff_r_sdss[idx_L], s=5, color='blue', alpha=alpha_type, label='Sc, Sd, Irr, Pec')
plt.plot([0, 50], [0, 50], ':', alpha=0.5, color='black')
plt.xlim(0, 50)
plt.ylim(0, 50)
plt.ylabel(r'$R_{50,r}^{P}$ from SDSS Database ["]', fontsize=15)
plt.xlabel(r'$R_{50,r}^{P}$ (This study) ["]', fontsize=15)
plt.tick_params(direction='in',top=True,right=True,labelsize='large')
#plt.gca().set_aspect('equal')
#plt.axis('scaled')
legend = plt.legend(fontsize=15)
#plt.setp(legend.)
plt.tick_params(direction='in')
plt.savefig('/Users/dhk/Documents/publish/ngcic/reff_r_petro_reff2.pdf')

# plt.figure(figsize=(5, 5))
# plt.scatter(reff_g_scratch[idx_E], reff_g_sdss[idx_E], s=3.0, color='red', alpha=0.5, label='E, S0')
# plt.scatter(reff_g_scratch[idx_S], reff_g_sdss[idx_S], s=3.0, color='green', alpha=0.5, label='Sa, Sb, Sbc')
# plt.scatter(reff_g_scratch[idx_L], reff_g_sdss[idx_L], s=3.0, color='blue', alpha=0.5, label='Sc, Sd, Irr, Pec')
# plt.plot([0, 50], [0, 50], ':', alpha=0.5, color='black')
# plt.xlim(0, 50)
# plt.ylim(0, 50)
# plt.ylabel(r'$R_{50,g}^{P}$ from SDSS Database ["]')
# plt.xlabel(r'$R_{50,g}^{P}$ (This study) ["]')
# plt.tick_params(direction='in',top=True,right=True,labelsize='large')
# #plt.gca().set_aspect('equal')
# #plt.axis('scaled')
# plt.legend()
# plt.tick_params(direction='in')
# plt.savefig('/Users/dhk/Documents/publish/ngcic/reff_g_petro.pdf')

# plt.figure(figsize=(5, 5))
# plt.scatter(pmag_g_scratch[idx_E], pmag_g_sdss[idx_E], s=3.0, color='red', alpha=0.5, label='E, S0')
# plt.scatter(pmag_g_scratch[idx_S], pmag_g_sdss[idx_S], s=3.0, color='green', alpha=0.5, label='Sa, Sb, Sbc')
# plt.scatter(pmag_g_scratch[idx_L], pmag_g_sdss[idx_L], s=3.0, color='blue', alpha=0.5, label='Sc, Sd, Irr, Pec')
# plt.plot([10, 15], [10, 15], ':', alpha=0.5, color='black')
# plt.xlim(15, 10)
# plt.ylim(15, 10)
# plt.ylabel(r'$mag_{g}^{Petro}$ from SDSS Database [AB mag]')
# plt.xlabel(r'$mag_{g}^{Petro}$ (This study) [AB mag]')
# #plt.gca().set_aspect('equal')
# #plt.axis('scaled')
# plt.tick_params(direction='in',top=True,right=True,labelsize='large')
# plt.legend()
# plt.tick_params(direction='in')
# plt.savefig('/Users/dhk/Documents/publish/ngcic/mag_g_petro.pdf')

plt.figure(figsize=(5, 5))
plt.scatter(pmag_r_scratch[idx_E], pmag_r_sdss[idx_E], s=5, color='red', alpha=alpha_type, label='E, S0')
plt.scatter(pmag_r_scratch[idx_S], pmag_r_sdss[idx_S], s=5, color='green', alpha=alpha_type, label='Sa, Sb, Sbc')
plt.scatter(pmag_r_scratch[idx_L], pmag_r_sdss[idx_L], s=5, color='blue', alpha=alpha_type, label='Sc, Sd, Irr, Pec')
plt.plot([10, 15], [10, 15], ':', alpha=0.5, color='black')
plt.xlim(15, 10)
plt.ylim(15, 10)
plt.ylabel(r'$m_{r}^{P}$ from SDSS Database', fontsize=15)
plt.xlabel(r'$m_{r}^{P}$ (This study)', fontsize=15)
#plt.gca().set_aspect('equal')
#plt.axis('scaled')
plt.tick_params(direction='in',top=True,right=True,labelsize='large')
plt.legend(fontsize=15)
plt.tick_params(direction='in')
plt.savefig('/Users/dhk/Documents/publish/ngcic/mag_r_petro2.pdf')

# plt.figure(figsize=(5, 5))
# plt.scatter(pmag_r2_scratch[idx_E], pmag_r_sdss[idx_E], s=3.0, color='red', alpha=0.5, label='E, S0')
# plt.scatter(pmag_r2_scratch[idx_S], pmag_r_sdss[idx_S], s=3.0, color='green', alpha=0.5, label='Sa, Sb, Sbc')
# plt.scatter(pmag_r2_scratch[idx_L], pmag_r_sdss[idx_L], s=3.0, color='blue', alpha=0.5, label='Sc, Sd, Irr, Pec')
# plt.plot([10, 15], [10, 15], ':', alpha=0.5, color='black')
# plt.xlim(15, 10)
# plt.ylim(15, 10)
# plt.ylabel(r'$mag_{r}^{Petro}$ from SDSS Database [AB mag]')
# plt.xlabel(r'$mag_{r}^{Petro}$ (This study) [AB mag]')
# #plt.gca().set_aspect('equal')
# #plt.axis('scaled')
# plt.tick_params(direction='in',top=True,right=True,labelsize='large')
# plt.legend()
# plt.tick_params(direction='in')
# plt.savefig('/Users/dhk/Documents/publish/ngcic/mag_r_petro2.pdf')

if hist_flag:

    bins = np.linspace(0,250,10)

    plt.figure(figsize=(5,5))
    plt.hist(reff_r_scratch[np.where(reff_sdss!=0)],bins,alpha=0.5,label='SDSS')
    plt.hist(reff_r_scratch[np.where(reff_sdss==0)],bins,alpha=0.5,label='No SDSS')
    plt.legend(loc='upper right')
    plt.savefig('/Users/dhk/Documents/publish/ngcic/reff_V_petro_hist2.pdf')

#np.savez('/Users/dhk/Documents/publish/ngcic/hlr_V_comb.npz', name = reff_name, hlr = reff_comb / 0.6)
