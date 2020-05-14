import random as rd
import numpy as np

b = range(9994)
c = np.zeros(1000)
for x in range(0,1000):
	b1 = rd.sample(b,275)
	b2 = rd.sample(b,491)
	b3 = [val for val in b1 if val in b2]
	c[x] = len(b3)

print mean(c)
