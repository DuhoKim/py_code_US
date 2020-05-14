from astropy.io import ascii
from astropy.coordinates import SkyCoord
import os
import shutil
import numpy as np

data_dir = '/Users/dhk/work/data/NGC_IC/SHA/'
cat_dir = '/Users/dhk/work/cat/NGC_IC/'

def copy_rename(old_fn, new_fn):
	shutil.copy(data_dir+'SHA_batch_download/'+old_fn, data_dir+'SHA_NGC_IC_LONG')
	os.rename(data_dir+'SHA_NGC_IC_LONG/'+old_fn, data_dir+'SHA_NGC_IC_LONG/'+new_fn)

cat=ascii.read(cat_dir+'sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match.csv')
sha=ascii.read(data_dir+'smByPosition.tbl',format='ipac')

for x in range(0,len(cat)):
	# select long mosaic for the galaxy in the catalog
	id_match = sha[(sha['Search_Tgt'] == cat['col1'][x])]
	long_bool = []
	for y in range(0,len(id_match)):
		long_bool.append('short' not in id_match['fname'][y])
	long_mosaic = id_match[long_bool]

	# lowest coordinate difference
	coord_ned = SkyCoord(cat['col2'][x],cat['col3'][x],unit='deg')
	coord_sel = SkyCoord(long_mosaic['ra'],long_mosaic['dec'],unit='deg')
	offs = coord_ned.separation(coord_sel)
	match = long_mosaic[offs == min(offs)]

	# copy and rename mosaic
	fn = match['fname'][0].split('/')[6]
	copy_rename(fn, cat['col1'][x]+'.fits')

	# copy and rename uncertainty mosaic
	fn_spl = fn.split('.')
	fn_spl[4] = fn_spl[4]+'_unc'
	unc_fn = '.'.join(fn_spl)
	copy_rename(unc_fn, cat['col1'][x]+'.unc.fits')	
