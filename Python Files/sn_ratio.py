import numpy as np
import matplotlib.pyplot as plt

zpt_errs = [0.010, 0.008, 0.007, 0.003]
stacked = [1,2,4,8]

full_list = np.arange(1,9,0.1)

root_n = zpt_errs[0] / np.sqrt(full_list)

plt.plot(full_list, root_n, linestyle = '--', color = 'black', label = r"Predicted $\dfrac{1}{\sqrt{N}}$")
plt.scatter(stacked, zpt_errs, color = 'royalblue', label = "Measured Data")
plt.xlabel("Number of Images Stacked")
plt.ylabel("Zero-Point Error")
plt.legend()

plt.savefig('StN_Proof.png', bbox_inches = 'tight', dpi=1000)
plt.show()