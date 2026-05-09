import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep

# )
from utils.helper import (clean_null_values, draw_fit, draw_fitted_xs,
                          get_background_label_list, get_canvas,
                          get_color_list, get_histograms_ratio, get_x_label,
                          save_figure)

# import numpy as np

plt.style.use(hep.style.CMS)
fig, ax = get_canvas()

hep.cms.label(
        data=True,
        label="",
        com="7.0",
        lumi=None,
        ax=ax,
    # Add lumi manually at the standard top-right position
ax.annotate(
        r"50 pb$^{-1}$               ",
        xy=(1, 1),
        xycoords="axes fraction",
        ha="right",
        va="bottom",
        fontsize=27,
    )
scan = np.loadtxt("mtop_likelihood_scan.txt", comments="#")
scan_signal = np.loadtxt("mtop_likelihood_scan_signal.txt", comments="#")

mtop = scan[:, 0]
delta = scan[:, 2]
delta_signal = scan_signal[:, 2]

delta = delta - np.min(delta)
delta_signal = delta_signal - np.min(delta_signal)

best_idx = np.argmin(delta)
best_idx_signal = np.argmin(delta_signal)
mt_best = mtop[best_idx]
mt_best_signal = mtop[best_idx_signal]


ax.plot(mtop, delta, color="blue", label="Data Stat. only")
ax.plot(mtop, delta_signal, color="red", label="Signal Stat. only")

ax.axhline(1.0, color="gray", linestyle="--", label=r"$1\sigma$")
ax.axhline(4.0, color="gray", linestyle=":", label=r"$2\sigma$")

ax.set_xlabel(r"$\mu$ [GeV]")
ax.set_ylabel(r"$-2\Delta\ln L$")

y_min, y_max = plt.ylim()

ax.set_ylim(0, 10)
plt.legend(bbox_to_anchor=(0.55, 0.7))

output_directory = "../plots/likelihood/"

output_name = "min_likelihood"


save_figure(fig, output_directory, output_name)

print("Best-fit mtop =", mt_best)
