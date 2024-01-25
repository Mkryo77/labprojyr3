import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import math
from tqdm import tqdm


# import asas-sn data

data = ascii.read('asas-sn_data/ASASSN-V_ct_cas.csv')
hjd = data['hjd']
mag = data['mag']

bin_size = 0.05    # sets bin size as a fraction of the period
folding_interval = 0.001    # folds in n day steps
fold_fraction = 0.1    # sets a maximum period to try folding as a fraction of the data size




# define a data folding function

def folding(hjds, period):
    hjdfold = np.empty(len(hjds))
    for n in range(len(hjds)):
        number = math.floor(hjds[n]/period)    # calculate the number of periods offset
        hjdfold[n] = hjds[n] - number * period    # subtract the number of periods * n
    return hjdfold




# fold data over a wide range of values

max_val = (max(hjd) - min(hjd)) * fold_fraction    # finds the max hjd to use fold as
periods = np.arange(folding_interval, 20 , folding_interval)     # creates a list of periods to try folding to 
master_variances = np.empty(len(periods))    # creates an empty array to store the variance for each period
folded_data = np.empty((len(periods),  2, len(hjd)))     # creates an empty 3D array to store all folded data in




for i in tqdm(range(0, len(periods)), desc = 'Folding...'):
    folded_hjds = folding(hjd, periods[i])    # fold data for all potential periods
    hjds_and_mags = np.vstack((folded_hjds, mag))
    folded_data[i] = hjds_and_mags

    # find variance within bins for all potential periods

    
    actual_bin_size = bin_size * periods[i]
    no_bins = math.floor(1 / bin_size)
    variances = []   # resets variance value
    bin_list = []
    for k in range(0, no_bins+1):    # creates a list of bin boundaries
        bin_list = np.append(bin_list, k*actual_bin_size)
        
    for k in range(0, len(bin_list)-1):    # iterates over list of bin boundaries

        lower_bound = bin_list[k]
        upper_bound = bin_list[k+1]
   
        points_in_bin = np.array([])
            
        for n in range(0, len(hjds_and_mags.T)):
            if hjds_and_mags[0][n] >= lower_bound and hjds_and_mags[0][n] < upper_bound:    # checks for points within the time bin
              points_in_bin = np.append(points_in_bin, hjds_and_mags[1][n])    # appends corresponding magnitudes to bin list
          
        if len(points_in_bin) >= 2:    # if there are fewer than two elements in the bin the variance is 0
            variances = np.append(variances, np.var(points_in_bin))
    
            
    total_variance = sum(variances)    # calculates the total variance for this folding size
        
    master_variances[i] = total_variance    # adds to list of variances for each folding size

min_variance = min(master_variances)    # finds minimum variance

mvs = master_variances.tolist()

min_index = mvs.index(min_variance)    # finds the corresponding period

optimal_period = periods[min_index]

print('The Optimal Period is ' + str(optimal_period) + ' days')

a, variances_around_mean, b = np.split(master_variances, [min_index-50, min_index+50])

a, varxs, b = np.split(periods, [min_index-50, min_index+50])


#plot folded data

plt.scatter(folded_data[min_index][0], folded_data[min_index][1])
plt.gca().invert_yaxis()
plt.title('Data Folded over a Period of ' + str(optimal_period) + ' days')

#plot variances

plt.scatter(varxs, variances_around_mean)
plt.xlabel("Period")
plt.ylabel("Variance")
plt.title("Variances of Several Period Values Centred Around the Chosen Period Value")
plt.savefig('variances_ctcas')
plt.show()


plt.scatter(periods, master_variances)
plt.xlabel("Period")
plt.ylabel("Variance")
plt.title("Variances of all period values")
plt.show()
