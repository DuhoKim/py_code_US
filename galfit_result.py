from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table, vstack
import numpy as np
import os
import shutil
import glob
from pandas import DataFrame, read_csv
import pandas as pd
from astropy.cosmology import Planck15 as Cosmo
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.backends.backend_pdf as pdf
from matplotlib.colors import LogNorm
import math

ned1 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_1.txt")
ned2 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_2.txt")
ned3 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_3.txt")
ned4 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_4.txt")
ned5 = ascii.read("/Users/dhk/work/cat/NGC_IC/ned_result_5.txt")

ned_tot = vstack([ned1, ned2, ned3, ned4, ned5])

file = r'/Users/dhk/work/cat/NGC_IC/ttype.xls'
ttype_cat = pd.read_excel(file)

rc3 = ascii.read("/Users/dhk/work/cat/NGC_IC/RC3/myrc3.dat")

out_dir = '/Users/dhk/Documents/publish/ngcic_rev/'

work_dir = '/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask/'

fn = sorted(glob.glob(work_dir + 'galfit.*'), key=os.path.getmtime, reverse=True)

btot = np.zeros(len(fn))  # bulge-to-total light ratios
btot[:] = np.nan  # Init w/ nan

btot_disk = np.zeros(len(fn))  # bulge-to-total light ratios
btot_disk[:] = np.nan  # Init w/ nan
btot_bulge = np.zeros(len(fn))  # bulge-to-total light ratios
btot_bulge[:] = np.nan  # Init w/ nan

mtot_disk = np.zeros(len(fn))  # bulge-to-total light ratios
mtot_disk[:] = np.nan  # Init w/ nan
mtot_bulge = np.zeros(len(fn))  # bulge-to-total light ratios
mtot_bulge[:] = np.nan  # Init w/ nan
mtot_disk_err = np.zeros(len(fn))  # bulge-to-total light ratios
mtot_disk_err[:] = np.nan  # Init w/ nan
mtot_bulge_err = np.zeros(len(fn))  # bulge-to-total light ratios
mtot_bulge_err[:] = np.nan  # Init w/ nan

mtot = np.zeros(len(fn))  # Total magnitude
mtot[:] = np.nan  # Init w/ nan
mtot_err = np.zeros(len(fn))  # Total magnitude
mtot_err[:] = np.nan  # Init w/ nan

ttype = np.zeros(len(fn))  # RC3 T-type
ttype[:] = np.nan  # Init w/ nan
ttype_err = np.zeros(len(fn))  # RC3 T-type error
ttype_err[:] = np.nan  # Init w/ nan

chi = np.zeros(len(fn))  # Reduced Chi squre
chi[:] = np.nan  # Init w/ nan
oned = np.zeros(len(fn))
oned[:] = np.nan
names = ["" for i in range(len(fn))]

# subimg = pdf.PdfPages(work_dir+'subs_3rd.pdf')
btot_fn = open(work_dir + 'btot_3rd.txt', 'w')
fit_log = open(work_dir + 'fit.log', 'r')
fit_log_lines = fit_log.readlines()

skip = False  # flag for check duplication

for i in range(0, len(fn)):
    # for i in range(36,37):
    res = open(fn[i])
    comp = 0
    disk = np.nan
    bulge = np.nan

    for line in res:
        line_spl = line.split(' ')
        if skip:
            break
        if '0) expdisk' in line:
            comp = 1
            continue
        if '0) devauc' in line:
            comp = 2
            continue

        for j in range(0, len(line_spl)):
            if 'feedme' in line_spl[j]:
                name_spl = line_spl[j].split('.')
                name = name_spl[0]

                same_name_idx = np.isin(names, name)
                if sum(same_name_idx):
                    skip = True
                    break
                else:
                    names[i] = name

                if name[0] == 'n':
                    galnum = name[3:].strip()
                    if len(galnum) == 3:
                        galnum = '0' + galnum
                    elif len(galnum) == 2:
                        galnum = '00' + galnum
                    elif len(galnum) == 1:
                        galnum = '000' + galnum
                    table_name = 'NGC ' + galnum
                elif name[0] == 'i':
                    galnum = name[2:].strip()
                    if len(galnum) == 3:
                        galnum = '0' + galnum
                    elif len(galnum) == 2:
                        galnum = '00' + galnum
                    elif len(galnum) == 1:
                        galnum = '000' + galnum
                    table_name = 'IC ' + galnum

                match = np.where(ttype_cat.iloc[:, 0] == table_name)
                ttype[i] = ttype_cat.iloc[match[0], 5]
                oned[i] = ttype_cat.iloc[match[0], 7]

            if 'Chi^2/nu' in line_spl[j]:
                chi[i] = float(line_spl[4][:-1])

            if 'Integrated' in line_spl[j] and comp == 1:
                disk = float(line_spl[2])
                btot_disk[i] = disk
                continue
            if 'Integrated' in line_spl[j] and comp == 2:
                bulge = float(line_spl[2])
                btot_bulge[i] = bulge
                continue

    if skip:
        skip = False
        continue

    for j in range(len(fit_log_lines) - 1, 0, -1):
        if name in fit_log_lines[j]:
            for k in range(0, 10):
                if 'expdisk' in fit_log_lines[j + k]:
                    try:
                        mtot_disk[i] = float(fit_log_lines[j + k].split(')')[1].split()[0])
                        mtot_disk_err[i] = float(fit_log_lines[j + k + 1].split(')')[1].split()[0])
                    except ValueError:
                        mtot_disk[i] = float(fit_log_lines[j + k].split(')')[1].split()[0][1:-1])
                        mtot_disk_err[i] = float(fit_log_lines[j + k + 1].split(')')[1].split()[0][1:-1])
                if 'devauc' in fit_log_lines[j + k]:
                    try:
                        mtot_bulge[i] = float(fit_log_lines[j + k].split(')')[1].split()[0])
                        mtot_bulge_err[i] = float(fit_log_lines[j + k + 1].split(')')[1].split()[0])
                    except ValueError:
                        mtot_bulge[i] = float(fit_log_lines[j + k].split(')')[1].split()[0][1:-1])
                        mtot_bulge_err[i] = float(fit_log_lines[j + k + 1].split(')')[1].split()[0][1:-1])
            break

    if name[0] == 'n':
        galnum = name[3:].strip()
        if len(galnum) == 3:
            galnum = '0' + galnum
        elif len(galnum) == 2:
            galnum = '00' + galnum
        elif len(galnum) == 1:
            galnum = '000' + galnum
        rc3_name = 'NGC' + galnum
    elif name[0] == 'i':
        galnum = name[2:].strip()
        if len(galnum) == 3:
            galnum = '0' + galnum
        elif len(galnum) == 2:
            galnum = '00' + galnum
        elif len(galnum) == 1:
            galnum = '000' + galnum
        rc3_name = 'IC' + galnum

    rc3_match = rc3[[j for j, s in enumerate(rc3['name']) if s.strip() == rc3_name]]
    if len(rc3_match) != 1:
        print('rc3 match is not one')
        ttype_err[i] = 0
    elif rc3_match['T'][0] == '*':
        ttype_err[i] = 0
    else:
        ttype_err[i] = rc3_match['T'][0]

    ned_idx = np.where(ned_tot['col2'] == name)

    if name == 'ngc6789':
        lum_dist = 3.6e5
        lum_dist_low = 3.55e5
        lum_dist_up = 3.65e5
    else:
        z = ned_tot['col8'][ned_idx[0]]
        z_unc = ned_tot['col9'][ned_idx[0]]
        lum_dist = Cosmo.luminosity_distance(z).value * 1e5  # in 10pc
        lum_dist_up = Cosmo.luminosity_distance(z + z_unc).value * 1e5  # in 10pc
        lum_dist_low = Cosmo.luminosity_distance(z - z_unc).value * 1e5  # in 10pc

    if np.isnan(bulge) and ~np.isnan(disk):
        btot[i] = 0
        mtot_apparent = disk
        mtot_bulge_err[i] = 99
    else:
        if disk > bulge:
            btot[i] = (10 ** ((disk - bulge) / 2.5) / (10 ** ((disk - bulge) / 2.5) + 1.0))
        if disk == bulge:
            btot[i] = 0.5
        if disk < bulge:
            btot[i] = (1.0 / (10 ** ((bulge - disk) / 2.5) + 1.0))
        mtot_apparent = 21.581 - 2.5 * np.log10(10 ** ((21.581 - disk) / 2.5) + 10 ** ((21.581 - bulge) / 2.5))

    mtot[i] = mtot_apparent - 5 * np.log10(lum_dist)
    mtot_zup = mtot_apparent - 5 * np.log10(lum_dist_up)
    mtot_zdown = mtot_apparent - 5 * np.log10(lum_dist_low)

    mtot_err_z = mtot_zdown - mtot_zup
    mtot_err_phot = np.min([mtot_disk_err[i], mtot_bulge_err[i]])
    mtot_err[i] = np.sqrt(mtot_err_z ** 2 + mtot_err_phot ** 2)

    btot_fn.writelines(
        name + ' ' + str(btot[i]) + ' ' + str(chi[i]) + ' ' + str(mtot[i]) + ' ' + str(mtot_err[i]) + ' ' +
        str(mtot_apparent) + ' \n')

    ########## ADD SUBTRACT PLOT in PDF ###########
    hdu = fits.open(work_dir + name + '.out.fits')
    img_orig = hdu[1].data
    img_model = hdu[2].data
    img_sub = hdu[3].data

    fig = plt.figure(figsize=(11, 3.85))

    h1 = fig.add_axes([0, 0, 0.33, 1])
    h2 = fig.add_axes([0.33, 0, 0.33, 1])
    h3 = fig.add_axes([0.66, 0, 0.33, 1])

    img1 = h1.imshow(img_orig, norm=LogNorm(vmin=0.001, vmax=5), origin='lower')
    img2 = h2.imshow(img_model, norm=LogNorm(vmin=0.001, vmax=5), origin='lower')
    img3 = h3.imshow(img_sub, norm=LogNorm(vmin=0.001, vmax=5), origin='lower')

    kpc_arcmin = Cosmo.kpc_proper_per_arcmin(z)
    arcmin_5kpc = 5.0 / kpc_arcmin
    frac_5kpc = arcmin_5kpc * 100.0 / img_orig.shape[0]

    props = dict(boxstyle='round', facecolor='wheat')
    h1.text(0.03, 0.92, r'IRAC 3.6$\mu$m', color='g', size=14, transform=h1.transAxes, bbox=props)
    h2.text(0.03, 0.92, 'Model with B/T=' + str('{:5.2}'.format(btot[i])), color='g', size=14, transform=h2.transAxes, bbox=props)
    h3.text(0.03, 0.92, r'Residual ($\chi^{2}$=' + str(chi[i])+')', color='g', size=14, transform=h3.transAxes, bbox=props)
    #h1.plot([0.05, 0.05 + frac_5kpc.value], [0.02, 0.02], color='darkslategray', transform=h1.transAxes, linewidth=5)
    #h1.text(0.05, 0.05, '5kpc, '+'{:4.2f}'.format(float(arcmin_5kpc.value)) + '\'', color='g', fontsize=24,
    #        transform=h1.transAxes, bbox=props)

    h1.get_xaxis().set_visible(False)
    h1.get_yaxis().set_visible(False)
    h2.get_xaxis().set_visible(False)
    h2.get_yaxis().set_visible(False)
    h3.get_xaxis().set_visible(False)
    h3.get_yaxis().set_visible(False)

    fig.savefig(work_dir + 'pdf/' + name + '.pdf')
    plt.close(fig)

btot_fn.close()
# subimg.close()

idx = np.where(ttype < 9)

fig = plt.figure(figsize=(8, 4))

h1 = fig.add_axes([0.5, 0.15, 0.45, 0.8])

cax = h1.scatter(btot[idx], ttype[idx], c=chi[idx], s=20, alpha=0.5)
fig.colorbar(cax, ticks=[0, 0.5, 1.0, 1.5], label=r'$\chi^{2}/N_{DOF}$')
# h1.set_ylabel('T-type')
h1.set_xlim([-0.1, 1.1])
h1.set_ylim([-6, 10])
h1.set_yticks([-5, -3, 0, 3, 6, 9])
h1.tick_params(direction='in', top=True, right=True, labelsize='large')
h1.set_xlabel('B/T (GALFIT)')
rect1 = patches.Rectangle((0.7, -5), 0.3, 2, linewidth=1, edgecolor='red', facecolor='none')
rect2 = patches.Rectangle((0.25, -5), 0.45, 2, linewidth=1, edgecolor='blue', facecolor='none')
rect3 = patches.Rectangle((0.5, -3), 0.3, 3, linewidth=1, edgecolor='red', facecolor='none')
rect4 = patches.Rectangle((0.13, -3), 0.37, 3, linewidth=1, edgecolor='blue', facecolor='none')
rect5 = patches.Rectangle((0.45, 0), 0.45, 2, linewidth=1, edgecolor='red', facecolor='none')
rect6 = patches.Rectangle((0.0, 0), 0.45, 2, linewidth=1, edgecolor='blue', facecolor='none')
rect7 = patches.Rectangle((0.25, 2), 0.65, 2, linewidth=1, edgecolor='red', facecolor='none')
rect8 = patches.Rectangle((0.0, 2), 0.25, 2, linewidth=1, edgecolor='blue', facecolor='none')

h1.add_patch(rect1)
h1.add_patch(rect2)
h1.add_patch(rect3)
h1.add_patch(rect4)
h1.add_patch(rect5)
h1.add_patch(rect6)
h1.add_patch(rect7)
h1.add_patch(rect8)

h1.text(-0.05, 8.5, '(b)', size=15)

# h1.errorbar(1.0,7.5,xerr=np.mean((mtot_d[idx]-mtot_u[idx])/2),yerr=np.mean(ttype_err[idx]), fmt='o', color='k')

h2 = fig.add_axes([0.1, 0.15, 0.35, 0.8])
h2.scatter(mtot[idx], ttype[idx], s=20, alpha=0.5)
h2.set_ylabel('T-type')
h2.set_xlim(-15, -27)
h2.set_ylim([-6, 10])
h2.set_yticks([-5, -3, 0, 3, 6, 9])
h2.set_xticks([-25, -20, -15])
h2.set_xlabel(r'"$M_{3.6\,\mu m}$"')
h2.tick_params(direction='in', top=True, right=True, labelsize='large')
h2.errorbar(-26, 8, xerr=np.mean(mtot_err[idx]), yerr=np.mean(ttype_err[idx]), color='k')

h2.text(-15.5, 8.5, '(a)', size=15)

# axs[1].hist(chi[idx])
# axs[1].set_xlabel(r'$\chi^{2}/N_{DOF}$')
# axs[1].set_ylabel('N')
# axs[1].set_xlim([0,3])
# axs[1].tick_params(direction='in',top=True,right=True,labelsize='large')
fig.savefig(out_dir + 'ttype_M_L_B2T.pdf')
plt.close(fig)
