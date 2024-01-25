from astropy.io import ascii
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt
import scipy
from scipy.optimize import curve_fit
import scipy.stats
from astropy.time import Time
from astropy.table import Table

def abs_mag(m, d):
    abs_mag = m - 5*np.log10(d) + 5
    return abs_mag

def linear(log_period, param_a, param_b):
    return param_a + param_b*log_period



'''
input = ascii.read('abs_mag_vs_para_dist_single_harm')
data = input
'''
#starnames_list = np.loadtxt("abs_mag_vs_para_dist_single_harm.txt", skiprows = 1, usecols=0)
#av_mag_list = np.loadtxt("abs_mag_vs_para_dist_single_harm.txt", skiprows = 1, usecols=1)
#dist_list = np.loadtxt("abs_mag_vs_para_dist_single_harm.txt", skiprows = 1, usecols=2)
#periods_list = np.loadtxt("abs_mag_vs_para_dist_single_harm.txt", skiprows = 1, usecols=3)

av_mag_list2 = np.loadtxt("optimised_params_two_harms.txt", skiprows = 1, usecols=3)
periods_list2 = np.loadtxt("optimised_params_two_harms.txt", skiprows = 1, usecols=2)
dist_list = np.loadtxt("abs_mag_vs_para_dist_two_harm.txt", skiprows = 1, usecols=2)

dist_err_list = np.loadtxt("abs_mag_vs_para_dist_two_harm.txt", skiprows = 1, usecols=4)
av_mag_err_list2 = np.loadtxt("optimised_params_two_harms.txt", skiprows = 1, usecols=9)
periods_err_list2 = np.loadtxt("optimised_params_two_harms.txt", skiprows = 1, usecols=8)

nstars = len(av_mag_list2)
#abs_mag_list = []
abs_mag_list2 = []
for i in range(0, nstars):
    #star = starnames_list[i]
    #av_mag = av_mag_list[i]
    av_mag2 = av_mag_list2[i]
    
    dist = dist_list[i]
    
    #M = abs_mag(av_mag, dist)

    #abs_mag_list += [M]

    M2 = abs_mag(av_mag2, dist)
    abs_mag_list2 += [M2]

    
#log_periods = np.log10(periods_list)
log_periods2 = np.log10(periods_list2)


plt.figure()
#plt.scatter(log_periods, abs_mag_list)
plt.scatter(log_periods2, abs_mag_list2, color = 'orange')
plt.gca().invert_yaxis()
plt.show()
