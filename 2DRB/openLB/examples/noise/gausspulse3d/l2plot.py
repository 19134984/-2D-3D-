import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['text.usetex'] = True

fig, ax = plt.subplots(figsize=(6,3), dpi=300)
outdirSuffixes = ["eternal_scale3", "periodic", "local", "interpolated", "sponge_bd20x1", "cbcSponge", "cbc"]
labels = ["ideal case (larger domain)", "periodic BC", "local BC", "interpolated BC",
          "sponge layer", "cbc outlet plus sponge", "cbc BC"]
for outdirSuffix, label in zip(outdirSuffixes, labels):
  filename = "tmp_"+outdirSuffix+"/gnuplotData/data/l2_absolute.dat"
  data = np.loadtxt(filename)
  ax.plot(data[:,0], data[:,1], label=label)
  
filename = "../refinedGausspulse3d/tmp_coarse/gnuplotData/data/l2_absolute_fine.dat"
data = np.loadtxt(filename)
ax.plot(data[:,0], data[:,1], label="coarse")

ax.set_title("Comparing far-field boundary conditions")
ax.set_yscale('log')
ax.set_ylim((1e-3,1.1))
ax.set_ylabel(r"$\frac{L_p}{L_{p0}}$", rotation=0, ha='right', fontsize=15)
ax.set_xlabel(r"$t_\mathrm{LU}$", fontsize=11, va='top')
ax.legend()
plt.tight_layout()
plt.show()
fig.savefig("l2ComparisonPlot.png")
