from astropy.io import ascii
import os
import shutil
import numpy as np

data_dir = '/Users/dhk/work/data/NGC_IC/SHA/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

def copy_rename(old_path_and_fn, old_fn, new_fn):
	shutil.copy(data_dir+'PBCD_batch_download/'+old_path_and_fn, data_dir+'PBCD_NGC_IC')
	os.rename(data_dir+'PBCD_NGC_IC/'+old_fn, data_dir+'PBCD_NGC_IC/'+new_fn)

cat=ascii.read(cat_dir+'sha_quarry_batch_257_without_prefix.txt')
sha=ascii.read(data_dir+'pbcdByPosition.tbl',format='ipac')

for x in range(2,len(cat)):
	# select mosaic for the galaxy in the catalog
	match = sha[(sha['Search_Tgt'] == cat['id'][x])]

	# copy and rename only for long mosaic
#	long_mosaic = [row for row in match['Instrument/Mode'] if 'short' not in row]
	ind = [i for i, s in enumerate(match['Instrument/Mode']) if 'PC' not in s]
	pbcd_path_and_fn = match['externalname'][ind]
	pbcd_fn = pbcd_path_and_fn[0].split('/')[3]
	copy_rename(pbcd_path_and_fn[0], pbcd_fn, 'PBCD_'+cat['id'][x]+'.fits')

