from astropy.io import ascii
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt
import scipy
from scipy.optimize import curve_fit
import scipy.stats
from astropy.time import Time

def single_harm(x,a,b,c):
    return a * np.sin(b*x)+c

def double_harm(x, aa,bb,cc,dd,ee):
    return aa*np.sin(bb*x) + cc*np.sin(dd*x) + ee

input = ascii.read('stars.py')
starnames = input['starnames']
nstars = len(starnames)

for i in range(0, nstars):
    data = ascii.read(starnames[i])
    ISOT = data['YYYY-MM-DDTHH:MM:SS']
    mag = data['mag']
    mag_err = data['mag_err']
    zpt = data['zpt']
    zpt_err = data['zpt_err']
    t = Time(ISOT, format = 'isot')
    mjd = t.mjd

    popt, _ = curve_fit(single_harm, mjd, mag)
    popt2, _ = curve_fit(double_harm, mjd,mag)
    a,b,c = popt
    aa,bb,cc,dd,ee = popt2
    print(popt)
    print(popt2)

    plt.errorbar(mjd,mag, yerr=mag_err, fmt="o")
    full_mjd = np.arange(np.min(mjd), np.max(mjd), 0.1)

    single = single_harm(full_mjd, a,b,c)
    
    plt.plot(full_mjd, single, label = 'Single Harmonic')

    double = double_harm(full_mjd, aa,bb,cc,dd,ee)
    plt.plot(full_mjd,double, label = 'Two Harmonics')
    plt.show()
