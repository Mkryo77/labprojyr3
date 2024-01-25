import matplotlib.pyplot as plt
import numpy as np

observed = np.array([6.847,12.14,3.811,13.64,9.81,8,6.24,26.1,23,18.19,16.4,4.973,10.15,2.042,2.660,13.85,11.6,15.69,8.64,8.003])
observed_err = np.array([0.004,0.05,0.009,0.04,0.03,0.01,0.03,0.3,0.3,0.04,0.1,0.004,0.01,0.006,0.003,0.04,0.1,0.03,0.01,0.006])
known = np.array([6.847,12.13888,3.81061,13.637747,9.812027,8.000669,6.23705,26.091,22.9519,18.193212,16.414812,4.972516,10.144698,2.0416,2.65956,13.84765,11.623537,15.69073,8.63793,8.002432])

o_c = observed - known

starname = ['ap_cas','ry_cas','ct_cas','sz_cas','dd_cas','dl_cas','fw_cas','ot_per','bm_per','yz_aur','rw_cam','as_per','sy_aur','gp_per','ew_aur','cy_aur','rx_aur','er_aur','hq_per','bk_aur']

plt.figure()
plt.errorbar(starname, o_c,yerr = observed_err, color = 'royalblue',fmt='o')
plt.tick_params(axis = 'x', labelrotation = 90)
plt.xlabel('Star Name')
plt.ylabel('O-C (Days)')
plt.savefig("O_C_vs_starname.png", bbox_inches = 'tight',dpi=1000)
plt.show()

plt.figure()
plt.errorbar(known, o_c, yerr = observed_err, color = 'royalblue',fmt='o')
plt.tick_params(axis = 'x', labelrotation = 90)
plt.xlabel('Known Period (Days)')
plt.ylabel('O-C (Days)')
plt.savefig("O_C_vs_knwon.png", bbox_inches = 'tight',dpi=1000)
plt.show()

print(np.mean(o_c))
print(o_c-observed_err)