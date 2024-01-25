def func(x, P0, m0, A1, phi1, A2, phi2):
    f0 = A1 * np.sin(2*pi/P0*x+phi1)
    f1 = A2 * np.sin(5*pi/P0*x+phi2)
    f = m0 + f0 + f1
    return f

import matplotlib
from astropy.io import ascii
import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.optimize import curve_fit
from astropy.time import Time

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'font.weight': 'normal'})
matplotlib.rcParams.update({'axes.linewidth': 3})
matplotlib.rcParams.update({'axes.labelweight': 'normal'})
matplotlib.rcParams.update({'xtick.direction': 'in'})
matplotlib.rcParams.update({'ytick.direction': 'in'})
matplotlib.rcParams.update({'xtick.major.width': '3'})
matplotlib.rcParams.update({'ytick.major.width': '3'})
matplotlib.rcParams.update({'xtick.minor.width': '3'})
matplotlib.rcParams.update({'ytick.minor.width': '3'})
matplotlib.rcParams.update({'xtick.major.size': '10'})
matplotlib.rcParams.update({'ytick.major.size': '10'})
matplotlib.rcParams.update({'xtick.minor.size': '10'})
matplotlib.rcParams.update({'ytick.minor.size': '10'})
matplotlib.rcParams.update({'axes.titlesize': 18})
matplotlib.rcParams.update({'axes.labelsize': 18})
matplotlib.rcParams.update({'xtick.labelsize': 18})
matplotlib.rcParams.update({'ytick.labelsize': 18})

## data file:
input = ascii.read('dl_cas')
P0 = 8.00067

### constants ###
pi = math.pi
#err_mag_lim = 0.0015  # magnitude error lower limit

### read data ###
data= input
ISOT = data['YYYY-MM-DDTHH:MM:SS']
mag = data['mag']
mag_err = data['mag_err']
zpt = data['zpt']
zpt_err = data['zpt_err']
t = Time(ISOT, format = 'isot')
mjd = t.mjd

time = mjd
err= mag_err
'''
# set y-error bars to minimum of err_mag_lim
for i in range(0,len(err)):
    if err[i] <= err_mag_lim:
        err[i] = err_mag_lim
'''

# sort data
i = np.argsort(time)
time = time[i]
mag = mag[i]
err = err[i]

# subtract JD of earliest observation
time0 = np.min(time)
time = time-time0

# plot data
plt.figure()
plt.scatter(time,mag)
plt.errorbar(time,mag,err,fmt='o')
plt.ylim(np.max(mag)+0.1,np.min(mag)-0.1)
plt.xlim(-5,np.max(time)+5)


#initial guess
# get average magnitude (starting guess for m0)
m0 = np.mean(mag)
p0 =      [P0, m0, 0.5, 0,  0.1, 0]

# bounds.  bound P0 to within 20% of initial guess
P0low    = P0 / 1.2
P0high   = P0 * 1.2
m0low    = m0 - 1.  # bound m0 to +/-1 mag
m0high   = m0 + 1.
A1low    = 0.
A1high   = 100
A2low    = 0.
A2high   = 100
phi1low  = 0.
phi1high = 1.2 * P0
phi2low  = 0.
phi2high = 1.2 * P0

bounds = ([P0low, m0low, A1low, phi1low, A2low, phi2low], [P0high, m0high, A1high, phi1high, A2high, phi2high])
popt, pcov = curve_fit(func, time, mag, sigma=err,bounds=bounds)
print(*popt)

# build finely sampled model
t = np.arange(0,500,0.1)
model = func(t,*popt)
plt.plot(t,model,color='red')
plt.xlabel('MJD- %5.3f' % (time0))
plt.ylabel('V mag')

# calculate chi^2
bestP0 = popt[0]
Tchi2 = np.sum((mag - func(time,*popt))**2/err**2)
Rchi2 = 1/(len(time)-len(popt)) * Tchi2
print('Best period: ',bestP0)
print('chi2_nu: ',Rchi2)

plt.show()


## build and plot folded LC for best-fit P0
#attempt to build and fit a light curve
# divide by period and subtract phase
time = time/bestP0
for i in range(0,len(time)):
    time[i] = time[i] - int(time[i])

t = t/bestP0
# just need one phase
gd = np.where(t <=1)
t = t[gd]
model = model[gd]
for i in range(0,len(t)):
    t[i] = t[i] - int(t[i])

plt.scatter(time,mag,color='blue')
plt.errorbar(time,mag,err,fmt='o')
plt.plot(t,model)
plt.xlabel('Phase')
plt.ylabel('V mag')
plt.ylim(np.max(mag)+0.1,np.min(mag)-0.1)
plt.xlim(-0.1,1.1)
plt.show()

# QED.
