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
#def single_harm(t,A0,T0,m0,phi):
   # return A0 * np.sin((2*np.pi/T0)*t + phi) + m0

def single_harm(t,A0,T0,m0,phi):
    return A0 * np.sin((2*np.pi/T0)*t + phi) + m0

def two_harm(t,A1,T0,m0,phi1,A2,phi2):
    return A1 *np.sin((2*np.pi/T0)*t + phi1) + A2 *np.sin((6*np.pi/T0)*t + phi2) + m0

def three_harm(t,A3,T0,m0,phi3,A4,phi4,A5,phi5):
     return A3*np.sin((2*np.pi/T0)*t + phi3) + A4*np.sin((6*np.pi/T0)*t + phi4) + A5*np.sin((4*np.pi/T0)*t + phi5) + m0


input = ascii.read('stars.py')
starnames = input['starnames']
nstars = len(starnames)
initial_values = ascii.read('initial_values')

lines = ['star A1 T0 m0 phi1 A2 phi2 A1err T0err m0err phi1err A2err phi2err \n']
with open('optimised_params_two_harms.txt', 'w') as f:
    f.writelines(lines)

for i in range(0, nstars):
    data = ascii.read(starnames[i])
    ISOT = data['YYYY-MM-DDTHH:MM:SS']
    mag = data['mag']
    mag_err = data['mag_err']
    zpt = data['zpt']
    zpt_err = data['zpt_err']
    t = Time(ISOT, format = 'isot')
    mjd = t.mjd
    print(f'{starnames[i]}')
   # if starnames[i] == 'sy_aur':
            
    # reset mag error to be sum of errors in quadrature
    mag_err = np.sqrt(mag_err**2 + zpt_err**2)
    
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
    A1_initial = A0_initial
    phi1_initial = phi_initial
    A2_initial = 0.5
    phi2_initial = 0.0

    A3_initial = A0_initial
    phi3_initial = phi_initial
    A4_initial = 0.5
    phi4_initial = 0.0
    A5_initial = 0.25
    phi5_initial = 0.0

    # bounds single harmonic
    T0low    = T0_initial / 1.5
    T0high   = T0_initial * 1.5
    m0low    = m0_initial - 1.  # bound m0 to +/-1 mag
    m0high   = m0_initial + 1.
    A0low    = 0.
    A0high   = 100
    philow  = 0.0
    phihigh = 2*np.pi #1.2 * T0_initial

    # two_harmonics
    A1low    = 0.
    A1high   = 100
    phi1low  = 0.0
    phi1high = 2*np.pi
    A2low    = 0.
    A2high   = 100
    phi2low  = 0.0
    phi2high = 2*np.pi
    
    #three harmonics
    A3low    = 0.
    A3high   = 100
    phi3low  = 0.0
    phi3high = 2*np.pi
    A4low    = 0.
    A4high   = 100
    phi4low  = 0.0
    phi4high = 2*np.pi
    A5low    = 0.
    A5high   = 100
    phi5low  = 0.0
    phi5high = 2*np.pi


    initial_guess_sing = [A0_initial, T0_initial, m0_initial, phi_initial]
    #bounds_sing = ([A0low, T0low, m0low, philow], [A0high, T0high, m0high, phihigh])
    #print(bounds)
    #print(initial_values)
    #popt, pcov = curve_fit(single_harm, mjd, mag, sigma = mag_err, absolute_sigma = True, bounds = bounds_sing, p0 = initial_guess_sing)
    popt, pcov = curve_fit(single_harm, mjd, mag, sigma = mag_err, absolute_sigma = True, p0 = initial_guess_sing)

    errs_cov = np.sqrt(np.diag(pcov))

    initial_guess_two = [A1_initial, T0_initial, m0_initial, phi1_initial, A2_initial, phi2_initial]
    bounds_two = ([A1low, T0low, m0low, phi1low, A2low, phi2low], [A1high, T0high, m0high, phi1high, A2high, phi2high])

    popt_two, pcov_two = curve_fit(two_harm, mjd, mag, sigma = mag_err, absolute_sigma = True,bounds=bounds_two, p0 = initial_guess_two)
    #popt_two, pcov_two = curve_fit(two_harm, mjd, mag, sigma = mag_err, absolute_sigma = True, bounds = bounds_two, p0 = initial_guess_two)

    errs_cov_two = np.sqrt(np.diag(pcov_two))
    
    initial_guess_three = [A3_initial, T0_initial, m0_initial, phi3_initial, A4_initial, phi4_initial, A5_initial, phi5_initial]
    bounds_three = ([A3low, T0low, m0low, phi3low, A4low, phi4low, A5low, phi5low], [A3high, T0high, m0high, phi3high, A4high, phi4high, A5high, phi5high])

    popt_three, pcov_three = curve_fit(three_harm, mjd, mag, sigma = mag_err, absolute_sigma = True, bounds = bounds_three, p0 = initial_guess_three)

    errs_cov_three = np.sqrt(np.diag(pcov_three))

    #print(f'{starnames[i]} Single Harmonic A0,T0,m0,phi0 = {popt}')
    #print(f'{starnames[i]} Two Harmonics A1,T0,m0,phi1,A2,phi2 = {popt_two}')

    #lines = [f'{starnames[i]} {popt_two[0]} {popt_two[1]} {popt_two[2]} {popt_two[3]} {popt_two[4]} {popt_two[5]} {errs_cov_two[0]} {errs_cov_two[1]} {errs_cov_two[2]} {errs_cov_two[3]} {errs_cov_two[4]} {errs_cov_two[5]}\n\n']
    #with open('optimised_params_two_harms.txt', 'a') as f:
    #    f.writelines(lines)
    '''
    fig, axs = plt.subplots(2,1,sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    fig.subplots_adjust(hspace=0)
    axs[0].set_title(f'{starnames[i]}')
    axs[0].errorbar(mjd,mag,yerr=mag_err,fmt='o', label = 'Data', color = 'royalblue')
    '''    
    
    full_mjd = np.arange(np.min(mjd)-5, np.max(mjd)+5, 0.1)
    
    
    single = single_harm(full_mjd, popt[0],popt[1],popt[2],popt[3])
    second_harmonic = two_harm(full_mjd, popt_two[0],popt_two[1],popt_two[2],popt_two[3],popt_two[4],popt_two[5])
    third_harmonic = three_harm(full_mjd, popt_three[0],popt_three[1],popt_three[2],popt_three[3],popt_three[4],popt_three[5],popt_three[6],popt_three[7])
    '''
    #plt.plot(full_mjd, single, label = 'Single Harmonics')
    axs[0].plot(full_mjd, single, label = "Single Harmonic", color = 'green')
    axs[0].plot(full_mjd, second_harmonic, label = 'Two Harmonics', color = 'orange')
    axs[0].set_xlim([np.min(full_mjd), np.max(full_mjd)])
    axs[0].set_ylabel("Apparent Magnitude")

    box0 = axs[0].get_position()
    axs[0].set_position([box0.x0, box0.y0, box0.width * 0.8, box0.height])

    axs[0].legend(bbox_to_anchor=(1, 0.5))
    axs[0].tick_params(top = False, bottom = True, left=True, right=True)
    axs[0].tick_params(axis ='x', direction = 'in', top = False, bottom = True)
    axs[0].tick_params(axis ='y', direction = 'inout', left = 'on', right = 'on')

    axs[0].xaxis.set_minor_locator(AutoMinorLocator())
    axs[0].yaxis.set_minor_locator(AutoMinorLocator())
    '''
    '''
    #axs[0].tick_params(axis='y', direction = 'out', left = 'on', right = 'off')
    #norm res plot
    res = (mag - two_harm(mjd, popt_two[0],popt_two[1],popt_two[2],popt_two[3],popt_two[4],popt_two[5]))/mag_err

    y_res_lim = [-np.max(abs(res))-1,np.max(abs(res))+1]

    axs[1].scatter(mjd, res, s = 8, color = 'royalblue')
    zeros = np.zeros(len(full_mjd))
    axs[1].plot(full_mjd, zeros, color = 'grey', linestyle = '--')

    stdev = np.std(res)
    axs[1].fill_between(full_mjd, stdev, -stdev, color = 'grey', alpha = 0.25, label = '1-sigma Interval')
    axs[1].set_ylabel("Normalised" '\n' "Residual" '\n' "(A)", fontsize=10)
    axs[1].set_xlim([np.min(full_mjd), np.max(full_mjd)])
    axs[1].set_ylim(y_res_lim)
    box1 = axs[1].get_position()
    axs[1].set_position([box1.x0, box1.y0, box1.width * 0.8, box1.height])
    axs[1].set_xlabel("Time (MJD)")
    axs[1].legend(loc = 'upper center', bbox_to_anchor = (0.5,1.05))
    axs[1].tick_params(top = True, bottom = True, left=True, right=True)
    #axs[1].set_xticklabels([59930,59950,59970,59990,60010])
    axs[1].tick_params(axis = 'x', labelrotation = 45)
    axs[1].tick_params(axis ='x', direction = 'inout', top = True, bottom = True)
    axs[1].tick_params(axis ='y', direction = 'inout', left = True, right = False)
    #axs[1].tick_params(axis='y', direction = 'in', left = False, right = True)
    axs[1].xaxis.set_minor_locator(AutoMinorLocator())

    #norm res histograms
    fig.add_axes((box1.width ,box1.y0,box1.width * 0.8 * 0.25,box1.height))
    plt.gca().yaxis.set_ticks_position("left")
    plt.tick_params(labeltop = False, labelbottom = True, labelleft=False, labelright=False)
    plt.tick_params(top = False, bottom = True, left=True, right=True)
    plt.hist(res, bins = 4, orientation = 'horizontal', edgecolor = 'black', linewidth = 0.8, color = 'royalblue')
    plt.ylim(y_res_lim)
    plt.xlabel("Occurence")
    plt.xlim([0,6])
    plt.xticks([0,2,4,6])
    mean_gauss = 0
    stdev_gauss = 1
    
    x_gauss = np.arange(-2.5,2.5,0.1)
    y_gauss = scipy.stats.norm(mean_gauss, stdev_gauss)
    # plot opposite for horizontal plot
    plt.plot(12*y_gauss.pdf(x_gauss), x_gauss, color = 'black', linestyle = 'dashed')
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
    '''

    import scipy.stats

    def chi_squared(model_params,model,x_data,y_data,y_error):
            return np.sum(((y_data - model(x_data, *model_params))/y_error)**2)

    chisq_min_two = chi_squared(popt_two, two_harm, mjd,mag,mag_err)
    deg_freedom_two = len(mjd) - len(initial_guess_two)
    print(f"DoF = {deg_freedom_two}")
    print(f"Minimised Chi-Squared = {chisq_min_two}")
    #print(f"a = {popt_two[0]}\n b = {popt_two[1]}")


    chisq_reduced_two = chisq_min_two/deg_freedom_two
    print('reduced chi^2 = {}'.format(chisq_reduced_two))


    P_value_two = scipy.stats.chi2.sf(chisq_min_two, deg_freedom_two)
    print('P(chi^2_min, DoF) = {}'.format(P_value_two))

    print('pcov =')
    print(pcov_two)

    errs_cov_two = np.sqrt(np.diag(pcov_two))
    print('errs_cov = {}'.format(errs_cov_two))


    #single harm
    chisq_min_single = chi_squared(popt, single_harm  , mjd,mag,mag_err)
    deg_freedom_single = len(mjd) - len(initial_guess_sing)
    print(f"DoF = {deg_freedom_single}")
    print(f"Minimised Chi-Squared = {chisq_min_single}")
    #print(f"a = {popt_two[0]}\n b = {popt_two[1]}")


    chisq_reduced_single = chisq_min_single/deg_freedom_single
    print('reduced chi^2 = {}'.format(chisq_reduced_single))


    P_value_single = scipy.stats.chi2.sf(chisq_min_single, deg_freedom_single)
    print('P(chi^2_min, DoF) = {}'.format(P_value_single))

    print('pcov =')
    print(pcov)

    errs_cov = np.sqrt(np.diag(pcov))
    print('errs_cov = {}'.format(errs_cov))


    # three harms


    chisq_min_three = chi_squared(popt_three, three_harm, mjd,mag,mag_err)
    deg_freedom_three = len(mjd) - len(initial_guess_three)
    print(f"DoF = {deg_freedom_three}")
    print(f"Minimised Chi-Squared = {chisq_min_three}")
    #print(f"a = {popt_two[0]}\n b = {popt_two[1]}")


    chisq_reduced_three = chisq_min_three/deg_freedom_three
    print('reduced chi^2 = {}'.format(chisq_reduced_three))


    P_value_three = scipy.stats.chi2.sf(chisq_min_three, deg_freedom_three)
    print('P(chi^2_min, DoF) = {}'.format(P_value_three))

    print('pcov =')
    print(pcov_three)

    errs_cov_three = np.sqrt(np.diag(pcov_three))
    print('errs_cov = {}'.format(errs_cov_three))

    harm_number = [1,2,3]
    chisq_min_list = [chisq_min_single, chisq_min_two, chisq_min_three]
    red_chisq_list = [chisq_reduced_single, chisq_reduced_two, chisq_reduced_three]

    plt.figure()
    plt.scatter(harm_number,chisq_min_list)
    plt.title(f"{starnames[i]}")
    plt.xticks([1,2,3])
    plt.xlabel("Number of Harmonics")
    plt.ylabel("Minimised Chi-Square")



    plt.savefig(f'{starnames[i]}_chisq_min_vs_harm.png',bbox_inches = 'tight', dpi=1000)


    plt.figure()
    plt.scatter(harm_number, red_chisq_list)
    plt.title(f"{starnames[i]}")
    plt.xticks([1,2,3])
    plt.xlabel("Number of Harmonics")
    plt.ylabel("Reduced Chi-Square")

    plt.savefig(f'{starnames[i]}_red_chisq_vs_harm.png',bbox_inches = 'tight', dpi=1000)
    
    plt.figure()
    plt.scatter(harm_number, chisq_min_list, color = 'royalblue', label = "Minimised $\u03C7^2$", marker = 'x', s = 200)
    plt.scatter(harm_number, red_chisq_list, color = 'orange', label = "Reduced $\u03C7^2$", marker = '+', s=200)
    plt.xticks([1,2,3])
    plt.xlabel("Number of Harmonics")
    plt.ylabel("$\u03C7^2$ Value")
    plt.legend()
    plt.savefig(f'{starnames[i]}_both_chisq_vs_harm.png',bbox_inches = 'tight', dpi=1000)
    
    plt.show()
