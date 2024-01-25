import matplotlib.pyplot as plt

starname = ['ct_cas','ry_cas','dl_cas','bm_per','hq_per','sy_aur','bk_aur','cy_aur','er_aur','ew_aur','fw_cas','rw_cam','sz_cas','yz_aur','rx_aur','ap_cas','as_per','dd_cas','gp_per','ot_per']
red_chi = [436.9399283798078,104.50435163895418,321.0040287078952,4.030237096191452,19.335377374201023,177.9783262691824,129.1477112011485,105.05261187813915,129.91587271780392,23.2706395677778,41.2862573822432,131.88301540144823,13.991780722691841,82.89029904830196,45.7777801483816,42.82394027739102,111.02229550196719,64.3020481090004,6.039368431114419,1.935723963038507]

plt.scatter(starname,red_chi, color = 'royalblue', label = 'Data')
plt.tick_params(axis = 'x', labelrotation = 90)
plt.xlabel('Star Name')
plt.ylabel('Reduced $\u03C7^2$')


plt.savefig("Red_chi_each_star.png", bbox_inches = 'tight', dpi = 1000)
plt.show()



average = [sum(red_chi) / 20]

print(average)
'''
av = average * 20

plt.plot(starname,av, color = 'black', linestyle = '--', label = 'Average')
'''