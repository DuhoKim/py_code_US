from astropy.io import fits
from astropy.io import ascii
import numpy as np
import os
import shutil

#tot=ascii.read("/Users/dhk/work/data/NGC_IC/SDSS/id_scrfn.txt")
tot=ascii.read("/Users/dhk/work/data/NGC_IC/SDSS/id_scrfn_missing3.txt")

for x in range(0,len(tot)):
#for x in range(0,1):
	fn='/Users/dhk/work/data/NGC_IC/SDSS/scripts/'+tot['scr'][x]
	with open(fn+'.0.sh') as old, open(fn+'.1.sh','w') as new:
		for line in old:
			index = line.find('https://')
			if index != -1:
				line = line[:index] + ' --no-check-certificate ' + line[index:]
			new.write(line)
	os.chmod(fn+'.1.sh',0o755)
	os.makedirs('/Users/dhk/work/data/NGC_IC/SDSS/run')
	os.chdir('/Users/dhk/work/data/NGC_IC/SDSS/run')
        os.system(fn+'.1.sh')
	for filename in os.listdir("."):
		if filename.startswith(tot['scr'][x]):
			os.rename(filename, '../'+filename[20]+'/'+tot['id'][x]+filename[19:])
	shutil.rmtree('/Users/dhk/work/data/NGC_IC/SDSS/run')
			

	
	
