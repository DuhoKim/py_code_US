from astropy.io import fits
import numpy as np
	
hdulist=fits.open('/Users/dhk/work/cat/NGC_IC/VII_118.fits')
tb=hdulist[1].data
for x in range(0,len(tb)/1000+1):
	f=open("sha_quarry_batch_%d.txt" % (x),"w")
	f.write("COORD_SYSTEM: Equatorial\n")
	f.write("EQUINOX: J2000\n")
	f.write("NAME-RESOLVER: NED\n")
	for y in range(x*1000,(x+1)*1000):
		if y == len(tb) :
			break
		if tb[y][1]==' Gx':
			if tb[y][0][0]=='I':
				f.write('ic'+tb[y][0][1:].strip()+'\n')
			else:
				f.write('ngc'+tb[y][0].strip()+'\n')
			
	f.close()
			
