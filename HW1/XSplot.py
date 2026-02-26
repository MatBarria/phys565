import matplotlib.pyplot as plt
import numpy as np
import uproot as ur
import os
import mplhep as hep

plt.style.use(hep.style.CMS)

def get_canvas(draw_ratio=False):
    if draw_ratio:
        fig, axs = plt.subplots(2, 1, height_ratios=[10, 2])
        fig.subplots_adjust(hspace=0.1)
        return fig, axs

    fig, axs = plt.subplots(1, 1)
    return fig, axs


def save_figure(fig, outputDirectory, name):

    os.makedirs(outputDirectory, exist_ok=True)
    fig.savefig(outputDirectory + name + ".pdf", bbox_inches="tight")
    fig.savefig(outputDirectory + name + ".png", bbox_inches="tight", dpi=300)
    fig.savefig(outputDirectory + name + ".pdf")
    fig.savefig(outputDirectory + name + ".png", dpi=300)
    print(outputDirectory + name + " Has been created")



colors = [
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:olive",
    "tab:cyan",
]
# Constants
alpha = 1/137     # fine-structure constant
s_list = [8,12,16,20,24,28,32]
s_labels= ["8 [GeV]" ,"12 [GeV]" ,"16 [GeV]","20 [GeV]","24 [GeV]","28 [GeV]","32 [GeV]"]
# Angular range (avoid theta = 0 to prevent divergence)
theta = np.linspace(1e-3, np.pi - 1e-3, 1000)
ctheta = np.linspace(-.9, .9 , 1000)

to_mbarn = 0.3894
# Differential cross section
fig, ax = get_canvas()

i =0
for s in s_list:
    dsigma_dOmega = (alpha**2 / (4*s**2)) * ((3 +np.cos(theta)**2) / (1 - np.cos(theta)))**2
    dsigma_dOmega = (alpha**2 / (4*s**2)) * ((3 +np.cos(theta)**2) / (1 - np.cos(theta)))**2
    # dsigma_dOmega = (alpha**2 / (4*s**2)) * (1 + np.cos(theta)**2) 
    ax.plot(ctheta, 1e6*to_mbarn*dsigma_dOmega, color = colors[i], label = s_labels[i])
    i = i +1
ax.set_ylabel(r"d$\sigma$/d(cos($\theta$)[nb] ", loc="center")
ax.set_xlabel(r"cos($\theta$) ", loc="center")
# ax.set_yscale("log")
# ax.set_ylim([0.01,10])
ax.legend()


# save_figure(fig, "./", "XStheta-log")
save_figure(fig, "./", "XStheta")
