from matplotlib.pylab import *
from numpy import *

n = 100
rhoN = 5.0

rho = linspace(1e-10, rhoN, n)
u = []

f = open("u1.txt", "r")

for line in f:
	u.append(float(line))

u = array(u)

plot(rho, u)
show()