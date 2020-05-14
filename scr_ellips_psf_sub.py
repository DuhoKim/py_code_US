from astropy.io import fits
import numpy as np

with open('/Users/dhk/work/data/NGC_IC/pilot5/sdss/J134744.00+381816.0.sh') as old, open('/Users/dhk/work/data/NGC_IC/pilot5/sdss/new.sh','w') as new:
	for line in old:
		index = line.find('https://')
		if index != -1:
			line = line[:index] + ' --no-check-certificate ' + line[index:]
		new.write(line)

