import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.usetex'] = True

fig, ax = plt.subplots(figsize=(6,3), dpi=300)
for sigma in [0.1, 0.03, 0.01, 0.003, 0.001]:
  suffix = f"_sigma{sigma}" if sigma != 0.01 else ""
  data = np.loadtxt(f"tmp_cbc{suffix}/gnuplotData/data/l2_absolute.dat")
  ax.plot(data[:,0], data[:,1], label=f"cbc, sigma={sigma}")

ax.set_title("Comparing scaling factor for K1 in characteristics-based BC")
ax.set_yscale('log')
ax.set_ylim((1e-3,1.1))
ax.set_ylabel(r"$\frac{L_p}{L_{p0}}$", rotation=0, ha='right', fontsize=15)
ax.set_xlabel(r"$t_\mathrm{LU}$", fontsize=11, va='top')
ax.legend()
plt.tight_layout()
plt.show()
fig.savefig("l2ComparisonPlot.png")
