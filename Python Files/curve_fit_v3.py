from astropy.io import ascii
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt
import scipy
from scipy.optimize import curve_fit
import scipy.stats
from astropy.time import Time

def single_harm(t,A0,T0,m0,phi):
    return A0 * np.sin((2*np.pi/T0)*t + phi) + m0

def two_harm(t,A1,T1,m0,phi1,A2,T2,phi2):
    return A1 *np.sin((2*np.pi/T1)*t + phi1) + A2 *np.sin((2*np.pi/T2)*t + phi2) + m0

input = ascii.read('stars.py')
starnames = input['starnames']
nstars = len(starnames)
initial_values = ascii.read('initial_values')

for i in range(0, nstars):
    data = ascii.read(starnames[i])
    ISOT = data['YYYY-MM-DDTHH:MM:SS']
    mag = data['mag']
    mag_err = data['mag_err']
    zpt = data['zpt']
    zpt_err = data['zpt_err']
    t = Time(ISOT, format = 'isot')
    mjd = t.mjd

    av_mag = np.mean(mag)
    m0_initial = av_mag
    ''' need below to work to loop through list'''
    A0_initial = initial_values['A0'][i]  #not sure what to guess
    T0_initial = initial_values['T0'][i]
    phi_initial = initial_values['phi'][i]
    '''
    print(A0_initial)
    print(T0_initial)
    print(phi_initial)
    '''

    # bounds
    T0low    = T0_initial / 1.2
    T0high   = T0_initial * 1.2
    m0low    = m0_initial - 1.  # bound m0 to +/-1 mag
    m0high   = m0_initial + 1.
    A0low    = 0.
    A0high   = 100
    philow  = 0.0
    phihigh = 1.2 * T0_initial

    
    
    initial_guess_sing = [A0_initial, T0_initial, m0_initial, phi_initial]
    bounds_sing = ([A0low, T0low, m0low, philow], [A0high, T0high, m0high, phihigh])
    #print(bounds)
    #print(initial_values)
    popt, pcov = curve_fit(single_harm, mjd, mag, sigma = mag_err, absolute_sigma = True, bounds = bounds_sing, p0 = initial_guess_sing)
    
    print(f'{starnames[i]} A0,T0,m0,phi0 = {popt}')
    plt.figure()
    plt.errorbar(mjd,mag,yerr=mag_err,fmt='o', label = 'Data')
    
    full_mjd = np.arange(np.min(mjd), np.max(mjd), 0.1)
    single = single_harm(full_mjd, popt[0],popt[1],popt[2],popt[3])
    
    plt.plot(full_mjd, single, label = 'Single Harmonics')
    plt.legend()
    plt.title(f'{starnames[i]}')
    plt.savefig(f'{starnames[i]}_with_curve.png',bbox_inches = 'tight', dpi=800)
    plt.show()
