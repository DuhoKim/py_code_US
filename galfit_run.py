from astropy.io import ascii
import os

os.chdir('/Users/dhk/work/data/NGC_IC/SHA/galfit/9th_mask')
sha_cat=ascii.read('/Users/dhk/work/cat/NGC_IC/sha_ned_match_off_lt_1_SDSS_100psfs_b_over_a_gt_05_from_Tom_2_erase_Tamura_match_segmentation.csv') # Sample NGC/IC numbers

for x in range(0,len(sha_cat)):
	name = sha_cat['col1'][x]
	os.system('galfit '+name+'.feedme')
