from astropy.io import fits
from astropy.io import ascii
import numpy as np
import os

cat=ascii.read("/Users/dhk/work/py/sha_ned_match_w_mag_off_lt_1_b_over_a_lt_half_409.txt")

counter=0
with open('/Users/dhk/work/cat/NGC_IC/14840.rtf') as fin:
	for line in fin:
		currentline = line.split(",")
		for word in currentline:
			comp1 = word.strip()
			for comp2 in cat['id']:
				if str(comp1).lower() == str(comp2).lower():
					print str(comp1).lower()
					print str(comp2).lower()
					counter=counter+1
print counter
