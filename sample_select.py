from astropy.io import ascii
import numpy as np

i1=ascii.read("/Users/dhk/work/cat/NGC_IC/irac1_gals.txt")
irac1_gals=[]
cnt=0
for x in range(0,len(i1)):
	data=ascii.read(
"/Users/dhk/work/cat/NGC_IC/SHA_mosaic_quarry_result_%d.tbl" % (x),format='ipac')
	tmp=""
	for y in range(0,len(data)):
		if data[y][10] == 'IRAC1':
#			if y == 0 and data[y][0].upper() in rc3['name']:
			if tmp == "" :
				tmp=data[y][0]
				cnt=cnt+1
				irac1_gals.append(tmp)
			elif data[y][0] != tmp:
				tmp=data[y][0]
				cnt=cnt+1
				irac1_gals.append(tmp)
print(cnt)	
thefile = open('irac1_gals.txt','w')
for item in irac1_gals:
	thefile.write("%s\n" % item)
	
