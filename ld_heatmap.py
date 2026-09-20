import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
m = np.loadtxt("results/ld_matrix.ld")
pos = [int(l.split()[3]) for l in open("results/ld_sub.bim")]
fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(m, cmap="Reds", vmin=0, vmax=1)
fig.colorbar(im, ax=ax, label="r$^2$")
if 39913696 in pos:
    i = pos.index(39913696)
    ax.axvline(i, color="blue", lw=0.6)
    ax.axhline(i, color="blue", lw=0.6)
ax.set_title("LD (r$^2$), PJL, 17q21 window, MAF>=0.10")
ax.set_xlabel("SNP number (ordered by position)")
fig.savefig("results/ld_heatmap.png", dpi=200, bbox_inches="tight")
print("saved", m.shape)
