from astropy.io import ascii
import os
import shutil
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table

data_dir = '/Users/dhk/work/data/NGC_IC/SHA/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')
pbcd=ascii.read(data_dir+'pbcdByPosition.tbl',format='ipac')
bcd=ascii.read(data_dir+'BCD_batch_download/bcdByPosition.tbl',format='ipac')

fig = plt.figure(figsize=(8,8))

h=fig.add_axes([0.1,0.6,0.4,0.4])
h2=fig.add_axes([0.6,0.6,0.4,0.4])
h3=fig.add_axes([0.1,0.1,0.4,0.4])
h4=fig.add_axes([0.6,0.1,0.4,0.4])

t = Table(names=('ID','NPBCD','VPBCD','NBCD','VBCD'),dtype=('S7','i4','f8','i4','f8'))

for x in range(0,len(cat)):	
	# select PBCDs for the galaxy in the catalog
        pbcd_id = pbcd[pbcd['Search_Tgt'] == cat['id'][x]]
	pbcd_sel = pbcd_id[[i for i,s in enumerate(pbcd_id['Instrument/Mode']) if 'PC' not in s]]
	
	pbcd_pa=[]
        for y in range(0,len(pbcd_sel)):
                try:
                        pbcd_hdu = fits.open(data_dir+'PBCD_batch_download/'+pbcd_sel['externalname'][y])
			pbcd_pa.append(float(pbcd_hdu[0].header['PA']))
                except IOError:
                        print 'No PBCD file of ', cat['id'][x]
        
	# select BCDs for the galaxy in the catalog
	bcd_id = bcd[bcd['Search_Tgt'] == cat['id'][x]]
	bcd_sel = bcd_id[[i for i,s in enumerate(bcd_id['Instrument/Mode']) if 'PC' not in s]]

	bcd_pa=[]
	for y in range(0,len(bcd_sel)):
		try:
			bcd_hdu = fits.open(data_dir+'BCD_batch_download/'+bcd_sel['externalname'][y])
			bcd_pa.append(float(bcd_hdu[0].header['PA']))
		except IOError:
			print 'No BCD file of ', cat['id'][x]
	
	
	t.add_row( (cat['id'][x], len(pbcd_sel), np.std(pbcd_pa), len(bcd_sel), np.std(bcd_pa)) )

h.hist((t['NPBCD'][np.where(t['NPBCD']==1)],t['NPBCD'][np.where(t['NPBCD']>1)]),bins=np.arange(0,20,1), \
	histtype='bar',label=('# of PBCD = 1 : '+str(sum(t['NPBCD']==1)),'# of PBCD > 1 : '+str(sum(t['NPBCD']>1))))
h.set_xlabel('# of Level 2 PBCDs for each mosaic')
h.set_ylabel('Frequency')
h.legend()
h.text(3,180,'PBCD : post-BCD (Level 2)',fontsize=10)
h.text(3,170,'BCD : Basic Calibrated Data (Level 1)',fontsize=10)


h2.hist((t['VPBCD'][np.where(t['NPBCD']==1)],t['VPBCD'][np.where(t['NPBCD']>1)]),bins=np.arange(0,120,5),histtype='bar')
h2.set_xlabel('Stdev of PAs in PBCDs')
h2.set_ylabel('Frequency')

h3.hist((t['NBCD'][np.where(t['NPBCD']==1)],t['NPBCD'][np.where(t['NPBCD']>1)]),bins=np.arange(0,200,10),histtype='bar')
h3.set_xlabel('# of Level 1 BCDs for each mosaic')
h3.set_ylabel('Frequency')


h4.hist((t['VBCD'][np.where(t['NPBCD']==1)],t['VPBCD'][np.where(t['NPBCD']>1)]),bins=np.arange(0,120,5),histtype='bar')
h4.set_xlabel('Stdev of PAs in BCDs')
h4.set_ylabel('Frequency')


fig.savefig("/Users/dhk/work/data/NGC_IC/SHA/psf_pa_var.pdf")



