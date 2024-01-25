from astropy.io import ascii
import numpy as np
from astropy.table import Table

data = ascii.read('ASASSN-V_ap_cas.csv')

hjd = data['hjd']
mag = data['mag']
mag_err = data['mag_err']

newdata = Table()
newdata['hjd'] = hjd
newdata['mag'] = mag
newdata['mag_err'] = mag_err

ascii.write(newdata, 'results.diff', overwrite = True)
