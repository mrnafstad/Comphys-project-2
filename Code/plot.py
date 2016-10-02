from matplotlib.pylab import *
from numpy import *

n = 100
rhoN = 4.0

rho = linspace(1e-10, rhoN, n)
u = []


f = open("nonint.txt", "r")

for line in f:
	u.append(float(line))

u = array(u)
v = []

f = open("interacting1.txt", "r")

for line in f:
	v.append(float(line))

v = array(v)

milli = []

f = open("intmilli.txt", "r")

for line in f:
	milli.append(float(line))

milli= array(milli)

half = []

f = open("inthalf.txt", "r")

for line in f:
	half.append(float(line))

half = array(half)



five = []

f = open("intfive.txt", "r")

for line in f:
	five.append(float(line))

five = array(five)

#just a random comment
plot(rho, u, "r")
hold("On")
plot(rho, v)
plot(rho, milli)
plot(rho, half)
plot(rho, five)
xlabel(r"$\rho$")
ylabel(r"$|u(\rho)|^2$")
legend(["Non interactive", r"interactive $\omega_r = 1.0$" r"interactive $\omega_r = 0.001$", r"interactive $\omega_r = 0.5$", r"interactive $\omega_r = 5$"])
show()