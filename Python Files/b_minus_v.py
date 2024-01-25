

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import math
from tqdm import tqdm
import sys
import os


if len(sys.argv) == 2:

    starname = sys.argv[1]
    
    if os.path.isfile('B/' + starname) == False:

        print(' ')
        print('File does not exist: B/' + starname + '.csv')
        print(' ')
        sys.exit()
else:

    print(' ')
    print(' Missing Starname ')
    sys.exit()


b_data = ascii.read('B/' + starname)
b_mag = b_data['mag']

v_data = ascii.read('V/' + starname)
v_mag = v_data['mag']

b_minus_v = b_mag - v_mag

with open('b_minus_v_results.txt', 'a') as f:
    f.write(starname + ' ' +  str(b_minus_v) + '\n')