
import numpy as np
import matplotlib
import math
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

def findDensity(l):
    den = []
    for i in l:
        val = i/(0.01**2)
        den.append(val)
    return den

npart = [50., 100., 200., 400., 800., 1000., 1200., 1500.]

boxL = findDensity(npart)

Kn = [0.78155, 0.30759, 0.265316, 0.17978, 0.138879, 0.090627, 0.088954, 0.0750727]

def fit(x, a, b, c):
    return a*np.exp(-b*x) + c 

popt, _ = curve_fit(fit, boxL, Kn)

a, b, c = popt
print a
print b
print c

xvals = np.linspace(int(boxL[0]), int(boxL[-1]), 200)

print xvals
yvals = []
for x in xvals:
    y = -b*np.exp(x) + 0.8
    print y
    yvals.append(y)

plt.plot(boxL, Kn, color='black')
plt.plot(xvals, yvals, color='red')
plt.xlabel('Number Density)')
plt.ylabel('Kn')
plt.yscale('log')
plt.xscale('log')
plt.show()
