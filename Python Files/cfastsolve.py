import numpy as np
#from pdmpy import pdm
import matplotlib
from matplotlib import pyplot as plt
from datetime import datetime
#import pwkit
#from pwkit import pdm
import sys
from astropy.io import ascii


def __to_format(jd: float, fmt: str) -> float:
    """
    Converts a Julian Day object into a specific format.  For
    example, Modified Julian Day.
    Parameters
    ----------
    jd: float
    fmt: str

    Returns
    -------
    jd: float
    """
    if fmt.lower() == 'jd':
        return jd
    elif fmt.lower() == 'mjd':
        return jd - 2400000.5
    elif fmt.lower() == 'rjd':
        return jd - 2400000
    else:
        raise ValueError('Invalid Format')


def __from_format(jd: float, fmt: str) -> (int, float):
    """
    Converts a Julian Day format into the "standard" Julian
    day format.
    Parameters
    ----------
    jd
    fmt

    Returns
    -------
    (jd, fractional): (int, float)
         A tuple representing a Julian day.  The first number is the
         Julian Day Number, and the second is the fractional component of the
         day.  A fractional component of 0.5 represents noon.  Therefore
         the standard julian day would be (jd + fractional + 0.5)
    """
    if fmt.lower() == 'jd':
        # If jd has a fractional component of 0, then we are 12 hours into
        # the day
        return math.floor(jd + 0.5), jd + 0.5 - math.floor(jd + 0.5)
    elif fmt.lower() == 'mjd':
        return __from_format(jd + 2400000.5, 'jd')
    elif fmt.lower() == 'rjd':
        return __from_format(jd + 2400000, 'jd')
    else:
        raise ValueError('Invalid Format')


def to_jd(dt: datetime, fmt: str = 'jd') -> float:
    """
    Converts a given datetime object to Julian date.
    Algorithm is copied from https://en.wikipedia.org/wiki/Julian_day
    All variable names are consistent with the notation on the wiki page.

    Parameters
    ----------
    fmt
    dt: datetime
        Datetime object to convert to MJD

    Returns
    -------
    jd: float
    """
    a = math.floor((14-dt.month)/12)
    y = dt.year + 4800 - a
    m = dt.month + 12*a - 3

    jdn = dt.day + math.floor((153*m + 2)/5) + 365*y + math.floor(y/4) - math.floor(y/100) + math.floor(y/400) - 32045

    jd = jdn + (dt.hour - 12) / 24 + dt.minute / 1440 + dt.second / 86400 + dt.microsecond / 86400000000

    return __to_format(jd, fmt)


def from_jd(jd: float, fmt: str = 'jd') -> datetime:
    """
    Converts a Julian Date to a datetime object.
    Algorithm is from Fliegel and van Flandern (1968)

    Parameters
    ----------
    jd: float
        Julian Date as type specified in the string fmt

    fmt: str

    Returns
    -------
    dt: datetime

    """
    jd, jdf = __from_format(jd, fmt)

    l = jd+68569
    n = 4*l//146097
    l = l-(146097*n+3)//4
    i = 4000*(l+1)//1461001
    l = l-1461*i//4+31
    j = 80*l//2447
    k = l-2447*j//80
    l = j//11
    j = j+2-12*l
    i = 100*(n-49)+i+l

    year = int(i)
    month = int(j)
    day = int(k)

    # in microseconds
    frac_component = int(jdf * (1e6*24*3600))

    hours = int(frac_component // (1e6*3600))
    frac_component -= hours * 1e6*3600

    minutes = int(frac_component // (1e6*60))
    frac_component -= minutes * 1e6*60

    seconds = int(frac_component // 1e6)
    frac_component -= seconds*1e6

    frac_component = int(frac_component)

    dt = datetime(year=year, month=month, day=day,
                  hour=hours, minute=minutes, second=seconds, microsecond=frac_component)
    return dt

def counter(i0,i1):
    print(str(int(np.round(i0/i1*100)))+'%', end='\r')



if len(sys.argv) != 4:
    print(' ')
    print('please input an object name, suggested period and number of bins')
    print(' ' )
    print('e.g. python test.py ct_cas 3.81 100')
    print(' ')
    sys.exit()

elif len(sys.argv) == 4:
    P0 = int(sys.argv[1])  # days
    f_period_scan = 0.005 # percentage.  0.1=10%
    nP = int(sys.argv[2]) # number of bins to scan
    nbins = 100 # number of bins to output?




data = ascii.read('ASASSN-V_' starnames[i] '.csv')
mjd = data['hjd']
mag = data['mag']
err = data['mag_err']

t = mjd.copy()

tmin = np.min(t)
t = t - tmin

#phase_min = 0.312
#phase_max = 0.317
#nphase = 100
nbins = 100

# set up period array to scan across
Pmin = P0 - f_period_scan * P0
Pmax = P0 + f_period_scan * P0
dP = (Pmax - Pmin)/nP
phase = Pmin + dP * np.arange(0,nP)
PD = np.zeros(int(nP))

print(P0)

#dp = (phase_max-phase_min)/nP
#phase = phase_min + dp*np.arange(nP)
#PD = np.zeros(nP)

for i in range(0,nP):
    p = t/phase[i]

    integers = p.astype(int)
    p = p-integers

    plt.clf()
    plt.subplot(1,2,1)
    plt.scatter(p,mag,s=1)
    #plt.ylim(
    plt.ylim(np.max(mag)+0.2,np.min(mag)-0.2)
#    plt.show(block=False)
#    plt.pause(0.2)
    
    phase_bin = np.arange(nbins)/nbins
    phase_scatter = np.zeros(nbins)
    for j in range(0,nbins-1):
        gd = np.where((p>=phase_bin[j]) & (p<=phase_bin[j+1]))
        if len(gd[0])>=2:
            this_mag = mag[gd]
            this_err = err[gd]

            #sig = np.std(1./this_err*this_mag)/np.sum(1./this_err)
            av = np.mean(this_mag)
            sig = np.std(this_mag)
            #off = av-this_err
            #sig = np.std(off)
            phase_scatter[j] = sig
            #plt.plot([phase_bin[j],phase_bin[j]],[7.4,8],c='k',alpha=0.1)
            #plt.plot([phase_bin[j+1],phase_bin[j+1]],[7.4,8],c='k',alpha=0.1)

            #plt.errorbar(phase_bin[j],av,sig,color='red')

    gd = np.where(phase_scatter)        
    PD[i] = np.mean(phase_scatter[gd[0]])

    print(i,PD[i])
    
    plt.subplot(1,2,2)
    plt.scatter(phase[0:i],PD[0:i])
    plt.show(block=False)
    plt.pause(0.02)

# find minimum phase dispersion:
gd = np.where(PD==np.min(PD))
bestP = phase[gd]
#print('best phase: ',bestP)
plt.plot([bestP,bestP],[np.min(PD),np.max(PD)])

plt.scatter(phase,PD)
plt.show(block=False)
plt.pause(0.01)

p = t/bestP
integers = p.astype(int)
p = p-integers
PD = np.zeros(nbins)
PDerr = np.zeros(nbins)
for j in range(0,nbins-1):
    gd = np.where((p>=phase_bin[j]) & (p<=phase_bin[j+1]))
    if len(gd[0])>=2:
        this_mag = mag[gd]
        this_err = err[gd]

        #sig = np.std(1./this_err*this_mag)/np.sum(1./this_err)
        av = np.mean(this_mag)
        sig = np.std(this_mag)
    
        PD[j] = av
        PDerr[j] = sig
PD[nbins-1] = PD[0]
PDerr[nbins-1] = PDerr[0]
plt.figure()
plt.scatter(phase_bin,PD)
plt.errorbar(phase_bin,PD,PDerr,ls='none')
plt.ylim(np.max(PD)+0.2,np.min(PD)-0.2)
plt.show()

a = open('model.dat','w')
a.write('# phase mag err')
plt.ylim(np.max(mag)+0.2,np.min(mag)-0.2)
for i in range(0, len(PD)):
    string = str(phase_bin[i])+' '+str(PD[i])+' '+str(PDerr[i])+'\n'
    a.write(string)
a.close()

