from astropy.io import ascii
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt
import scipy
from scipy.optimize import curve_fit
import scipy.stats
from astropy.time import Time

input = ascii.read('stars.py')
starnames = input['starnames']
nstars = len(starnames)

def single_harmonic(mjd, *param_vals):
    a = param_vals[0]
    b = param_vals[1]
    return av_mag + a * np.sin(b * mjd)

def chi_square(model_params,model,mjd,mag,mag_err):
    return np.sum(((mag-model(mjd,*model_params))/mag_err)**2)

def double_harmonic(mjd, *param_vals2):
    a = param_vals2[0]
    b = param_vals2[1]
    c = param_vals2[2]
    d = param_vals2[3]
    h1 = a * np.sin(b * mjd)
    h2 = c * np.sin(d * mjd)
    return av_mag + h1 + h2

for i in range(0, nstars):
    data = ascii.read(starnames[i])
    ISOT = data['YYYY-MM-DDTHH:MM:SS']
    mag = data['mag']
    mag_err = data['mag_err']
    zpt = data['zpt']
    zpt_err = data['zpt_err']
    t = Time(ISOT, format = 'isot')
    mjd = t.mjd
    #plt.plot(mjd, mag)
    
    plt.scatter(mjd, mag)
    plt.errorbar(mjd,mag, yerr=mag_err, fmt="o")
    plt.gca().invert_yaxis()
    plt.xlabel('Time (Modified Julian Time)')
    plt.ylabel('Apparent Magnitude (Mags)')
    plt.title(f'{starnames[i]}')
    plt.show()
    plt.savefig(f'{starnames[i]}.png', bbox_inches = 'tight', dpi = 500)
    av_mag = np.mean(mag)
    initial_values = np.array([0.0,1.0])

    deg_freedom = mjd.size - initial_values.size
    print(f"{starnames[i]} Single Harmonic1 DoF = {deg_freedom}")
    popt, pcov = curve_fit(single_harmonic, mjd, mag, sigma = mag_err, absolute_sigma = True, p0 = initial_values)

    print(f"{starnames[i]} Single Harmonic1 Optimised parameters = {popt}")

    a_solution, b_solution = popt

    chisq_min = chi_square(popt, single_harmonic, mjd, mag, mag_err)
    print(f"{starnames[i]} Single Harmonic1 Minimised chi square = {chisq_min}")

    chisq_reduced = chisq_min/deg_freedom
    print(f"{starnames[i]} Single Harmonic1 Reduced chi squared = {chisq_reduced}")

    p_value = scipy.stats.chi2.sf(chisq_min, deg_freedom)
    print(f"{starnames[i]} Single Harmonic1 p value = {p_value}")

    errs_cov = np.sqrt(np.diag(pcov))
    a_err, b_err = errs_cov
    print(f"{starnames[i]} Single Harmonic1: a = {a_solution} +/- {a_err}")
    print(f"{starnames[i]} Single Harmonic1: b = {b_solution} +/- {b_err}")   
  
    params, params_covariance = curve_fit(single_harmonic, mjd, mag, p0 = initial_values)
    print(f"{starnames[i]} Single Harmonic2 {params}")



    initial_values2 = np.array([1.0,1.0,1.0,1.0])

    deg_freedom2 = mjd.size - initial_values2.size
    print(f"{starnames[i]} Double Harmonic1 DoF = {deg_freedom2}")
    popt2, pcov2 = curve_fit(double_harmonic, mjd, mag, sigma = mag_err, absolute_sigma = True, p0 = initial_values2)

    print(f"{starnames[i]} double Harmonic1 Optimised parameters = {popt2}")

    a_solution2, b_solution2, c_solution2, d_solution2 = popt2

    chisq_min2 = chi_square(popt2, single_harmonic, mjd, mag, mag_err)
    print(f"{starnames[i]} Double Harmonic1 Minimised chi square = {chisq_min2}")

    chisq_reduced2 = chisq_min2/deg_freedom2
    print(f"{starnames[i]} Double Harmonic1 Reduced chi squared = {chisq_reduced2}")

    p_value2 = scipy.stats.chi2.sf(chisq_min2, deg_freedom2)
    print(f"{starnames[i]} Double Harmonic1 p value = {p_value2}")

    errs_cov2 = np.sqrt(np.diag(pcov2))
    a_err2, b_err2, c_err2, d_err2 = errs_cov2
    print(f"{starnames[i]} Doulbe Harmonic1: a = {a_solution2} +/- {a_err2}")
    print(f"{starnames[i]} Double Harmonic1: b = {b_solution2} +/- {b_err2}")
    print(f"{starnames[i]} Double Harmonic1: c = {c_solution2} +/- {c_err2}")
    print(f"{starnames[i]} Double Harmonic1: d = {d_solution2} +/- {d_err2}") 
  
    params2, params_covariance2 = curve_fit(double_harmonic, mjd, mag, p0 = initial_values2)
    print(f"{starnames[i]} Double Harmonic2 {params2}")

    full_times = np.arange(np.min(mjd), np.max(mjd), 0.1)
    
    single_1 = single_harmonic(full_times, *popt)
    plt.plot(full_times, single_1, label = 'Single Harmonic 1')
    plt.errorbar(mjd,mag, yerr=mag_err, fmt="o", label = 'Data')

    double_1 = double_harmonic(full_times, *popt2)
    plt.plot(full_times, double_1, label = 'Double Harmonic 1')
    
    plt.plot(full_times, single_harmonic(full_times, params[0], params[1]), label = 'Single Harmonic 2')
    plt.plot(full_times, double_harmonic(full_times, params2[0], params2[1], params2[2], params2[3]), label = 'Double Harmonic 2')
    
    plt.gca().invert_yaxis()
    plt.legend(loc = 'best')
    plt.title(f'{starnames[i]}')
    plt.show()
    plt.savefig(f'{starnames[i]}_With_Harmonics.png', bbox_inches = 'tight', dpi = 500)
