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

av_mag_list2 = np.loadtxt("optimised_params_two_harms_no_co_aur.txt", skiprows = 1, usecols=3)
periods_list2 = np.loadtxt("optimised_params_two_harms_no_co_aur.txt", skiprows = 1, usecols=2)
dist_list = np.loadtxt("abs_mag_vs_para_dist_two_harm_no_co_aur.txt", skiprows = 1, usecols=2)

dist_err_list = np.loadtxt("abs_mag_vs_para_dist_two_harm_no_co_aur.txt", skiprows = 1, usecols=4)
av_mag_err_list2 = np.loadtxt("optimised_params_two_harms_no_co_aur.txt", skiprows = 1, usecols=9)
periods_err_list2 = np.loadtxt("optimised_params_two_harms_no_co_aur.txt", skiprows = 1, usecols=8)

nstars = len(av_mag_list2)
#abs_mag_list = []
abs_mag_list2 = []
abs_mag_err_list2 = []


for i in range(0, nstars):
    #star = starnames_list[i]
    #av_mag = av_mag_list[i]
    av_mag2 = av_mag_list2[i]
    av_mag_err2 = av_mag_err_list2[i]
    
    dist = dist_list[i]
    dist_err = dist_err_list[i]
    #M = abs_mag(av_mag, dist)

    #abs_mag_list += [M]

    M2 = abs_mag(av_mag2, dist)
    abs_mag_list2 += [M2]

    # for the errors:
    av_mag_err_on_M2 = abs(abs_mag(av_mag2 + av_mag_err2, dist) - M2)
    dist_err_on_M2 = abs(abs_mag(av_mag2, dist + dist_err) - M2)

    M2_err = np.sqrt(av_mag_err_on_M2 ** 2 + dist_err_on_M2 ** 2)
    abs_mag_err_list2 += [M2_err]
    
#log_periods = np.log10(periods_list)
log_periods2 = np.log10(periods_list2)
err_log_per2 = periods_err_list2 / (np.log(10)*periods_list2)

linear_guess = np.array([1.0,4.0])
deg_freedom = log_periods2.size - linear_guess.size

from scipy.optimize import curve_fit

popt, pcov = curve_fit(linear, log_periods2, abs_mag_list2, sigma = abs_mag_err_list2, absolute_sigma = True, p0 = linear_guess)

def chi_squared(model_params,model,x_data,y_data,y_error):
        return np.sum(((y_data - model(x_data, *model_params))/y_error)**2)

chisq_min = chi_squared(popt, linear, log_periods2, abs_mag_list2,abs_mag_err_list2)
print(err_log_per2)
print(f"DoF = {deg_freedom}")
print(f"Minimised Chi-Squared = {chisq_min}")
print(f"a = {popt[0]}\n b = {popt[1]}")

plt.figure()
#plt.scatter(log_periods, abs_mag_list)
plt.scatter(log_periods2, abs_mag_list2, color = 'orange')
plt.errorbar(log_periods2, abs_mag_list2, xerr=err_log_per2, yerr=abs_mag_err_list2, fmt='o')
plt.plot(log_periods2, linear(log_periods2, *popt))
plt.gca().invert_yaxis()
plt.savefig("Abs Mag vs Log Period (NO co_aur).png", bbox_inches = 'tight', dpi = 1000)
plt.show()
