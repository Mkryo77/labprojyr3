from astropy.io import ascii
import numpy as np
import sys
import matplotlib
from matplotlib import pyplot as plt

#input = ascii.read('aperture_size_vs_StN.txt')
data = ascii.read('aperture_size_vs_StN.txt')

aperture_size = data['radius']
mean_counts = data['mean_count']
error_in_counts = data['err_count']
sky_bkg = data['sky_bkg']
sum_ap = data['sum_ap']

'''
ap_sn_file = 'aperture_size_vs_StN.txt'
ap_sn_file_r = np.genfromtxt(ap_sn_file,dtype='float',delimiter=' ')
'''


sig_to_noise = mean_counts/error_in_counts

sig_to_noise_2 = mean_counts/sky_bkg

plt.figure(1)
plt.plot(aperture_size, sig_to_noise)
plt.xlabel('Aperture Semi-major Axis (arcsec)')
plt.ylabel('Signal-to-Noise Ratio')
plt.xticks([0,2,4,6,8,10,12,14,16,18,20])
#plt.xlim([0,10])
plt.savefig("Signal_to_Noise_vs_Aperture_Size.png", bbox_inches = 'tight', dpi = 1000)

'''
plt.figure(2)
plt.plot(aperture_size, sig_to_noise_2)
plt.xlabel('Aperture Semi-major Axis (arcsec)')
plt.ylabel('Signal to Noise Ratio')
#plt.xlim([0,10])
plt.savefig("Signal_to_Noise_vs_Aperture_Size2.png", bbox_inches = 'tight', dpi = 400)
'''

plt.figure(3)
plt.errorbar(aperture_size, sum_ap, yerr=error_in_counts, ecolor = 'black')
plt.plot(aperture_size, sum_ap)
plt.xlabel('Aperture Semi-major Axis (arcsec)')
plt.ylabel('Total Counts in Aperture')
plt.xticks([0,2,4,6,8,10,12,14,16,18,20])
plt.savefig("Total_Counts_vs_Aperture_Size.png", bbox_inches = 'tight', dpi = 1000)

plt.show()
sys.exit()
