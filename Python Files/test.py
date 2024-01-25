from astropy.io import ascii
import numpy
import sys
import matplotlib
from matplotlib import pyplot as plt



input = ascii.read('stars')
starnames = input['starnames']
nstars = len(starnames)

for i in range(0,nstars):
    data = ascii.read(starnames[i])
    date = data['col1']
    time = data['col2']
    mag = data['col3']
    err = data['col4']
    print(date,time,mag,err)

plt. plot(mag, date)

    
    sys.exit()
