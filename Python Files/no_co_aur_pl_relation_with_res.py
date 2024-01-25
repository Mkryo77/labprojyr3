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
from matplotlib.ticker import AutoMinorLocator

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


chisq_reduced = chisq_min/deg_freedom
print('reduced chi^2 = {}'.format(chisq_reduced))

import scipy.stats

P_value = scipy.stats.chi2.sf(chisq_min, deg_freedom)
print('P(chi^2_min, DoF) = {}'.format(P_value))

print('pcov =')
print(pcov)

errs_cov = np.sqrt(np.diag(pcov))
print('errs_cov = {}'.format(errs_cov))

a_err, b_err = errs_cov

# Equivalently
# a_err = errs_cov[0]
# b_err = errs_cov[1]

print('Parameter a = {} +/- {}'.format(popt[0], a_err))           
print('Parameter b = {} +/- {}'.format(popt[1], b_err))  

'''
plt.figure()
#plt.scatter(log_periods, abs_mag_list)
plt.scatter(log_periods2, abs_mag_list2, color = 'orange')
plt.errorbar(log_periods2, abs_mag_list2, xerr=err_log_per2, yerr=abs_mag_err_list2, fmt='o')
plt.plot(log_periods2, linear(log_periods2, *popt))
plt.gca().invert_yaxis()
'''

fig, axs = plt.subplots(2,1,sharex=True, gridspec_kw={'height_ratios': [3, 1]})
fig.subplots_adjust(hspace=0)
axs[0].errorbar(log_periods2,abs_mag_list2,xerr=err_log_per2, yerr=abs_mag_err_list2,fmt='o', label = 'Data', color = 'royalblue')
axs[0].plot(log_periods2, linear(log_periods2, *popt),color = 'orange',label='Linear Model')
axs[0].legend()

box0 = axs[0].get_position()
axs[0].set_position([box0.x0, box0.y0, box0.width * 0.8, box0.height])

#axs[0].legend(bbox_to_anchor=(1, 0.5))
axs[0].tick_params(top = False, bottom = True, left=True, right=True)
axs[0].tick_params(axis ='x', direction = 'in', top = False, bottom = True)
axs[0].tick_params(axis ='y', direction = 'inout', left = 'on', right = 'on')
axs[0].xaxis.set_minor_locator(AutoMinorLocator())
axs[0].yaxis.set_minor_locator(AutoMinorLocator())
axs[0].invert_yaxis()
axs[0].set_ylabel("Absolute Magnitude")
res = (abs_mag_list2 - linear(log_periods2, popt[0],popt[1]))/abs_mag_err_list2
y_res_lim = [-np.max(abs(res))-1,np.max(abs(res))+1]

full_logp = np.arange(0.2, 1.6, 0.1)

axs[1].scatter(log_periods2, res, s = 8, color = 'royalblue')
zeros = np.zeros(len(full_logp))
axs[1].plot(full_logp, zeros, color = 'grey', linestyle = '--')
stdev = np.std(res)

print(f"Mean normalised res = {np.mean(res)}")
print(stdev)

axs[1].fill_between(full_logp, stdev, -stdev, color = 'grey', alpha = 0.25, label = '1-sigma Interval')
axs[1].set_ylabel("Normalised" '\n' "Residual" '\n' "(A)", fontsize=10)
axs[1].set_xlim([0.2, 1.5])
axs[1].set_ylim(y_res_lim)
box1 = axs[1].get_position()
axs[1].set_position([box1.x0, box1.y0, box1.width * 0.8, box1.height])

axs[1].set_xlabel("log(Period [Days])")
axs[1].legend(loc = 'upper center', bbox_to_anchor = (0.3,0.375))
axs[1].tick_params(top = True, bottom = True, left=True, right=True)
#axs[1].set_xticklabels([59930,59950,59970,59990,60010])
axs[1].tick_params(axis = 'x')
axs[1].tick_params(axis ='x', direction = 'inout', top = True, bottom = True)
axs[1].tick_params(axis ='y', direction = 'inout', left = True, right = False)
#axs[1].tick_params(axis='y', direction = 'in', left = False, right = True)
axs[1].xaxis.set_minor_locator(AutoMinorLocator())

#norm res histograms
fig.add_axes((box1.width ,box1.y0,box1.width * 0.8 * 0.25,box1.height))
plt.gca().yaxis.set_ticks_position("left")
plt.tick_params(labeltop = False, labelbottom = True, labelleft=False, labelright=False)
plt.tick_params(top = False, bottom = True, left=True, right=True)
plt.hist(res, bins = 7, orientation = 'horizontal', edgecolor = 'black', linewidth = 0.8, color = 'royalblue')
plt.ylim(y_res_lim)
plt.xlabel("Occurence")
plt.xlim([0,8])
plt.xticks([0,2,4,6,8])
mean_gauss = 0
stdev_gauss = 1

x_gauss = np.arange(-2.5,2.5,0.1)
y_gauss = scipy.stats.norm(mean_gauss, stdev_gauss)
# plot opposite for horizontal plot
plt.plot(12*y_gauss.pdf(x_gauss), x_gauss, color = 'black', linestyle = 'dashed')
plt.gca().xaxis.set_minor_locator(AutoMinorLocator())

plt.savefig("Abs Mag vs Log Period (NO co_aur).png", bbox_inches = 'tight', dpi = 1000)
plt.show()


print(abs_mag_list2)
print(abs_mag_err_list2)