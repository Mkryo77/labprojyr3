from astropy.io import ascii
import numpy
import sys
import matplotlib
from matplotlib import pyplot as plt
from astropy.time import Time


input = ascii.read('stars')
starnames = input['starnames']
nstars = len(starnames)

for i in range(0,nstars):
    data = ascii.read(starnames[i])
    date_time = data['dateTtime']
    mag = data['mag']
    err = data['mag_err']
    t = Time('date_time', format='isot')
    mjd = t.mjd
    print(mjd,mag,err)


    plt.scatter(mjd, mag)
#plt.errorbar(date,mag,yerr=err)
    plt.gca().invert_yaxis()
    plt.show()
    
sys.exit()
