import uproot as ur
import ROOT as root
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt

from utils.helper import (
    get_canvas,
    save_figure,
    get_background_label_list,
    get_color_list,
)
from utils.constants import (
    ALL_BRANCHES,
    BACKGROUND_SOURCES,
    SIGNAL_SOURCES,
    SOURCES_LABEL,
    N_BINS,
    X_RANGE,
)

plt.style.use(hep.style.CMS)


def draw_sig_and_bg(variable):

    tuple_path = "../data/proccess_tuples/"

    print("PLOTTING: ", variable)
    background_tree = ur.open(tuple_path + "background.root:tree_output")
    background_branches = background_tree.arrays([variable, "weight", "triggerIsoMu24", "Nlep_valid"], library="np")  # type: ignore

    signal_tree = ur.open(tuple_path + "signal.root:tree_output")
    signal_branches = signal_tree.arrays([variable, "weight", "triggerIsoMu24", "Nlep_valid"], library="np")  # type: ignore

    bool_list_bkg = np.ones_like(background_branches[variable], dtype=bool)
    bool_list_bkg = (bool_list_bkg) & (background_branches["triggerIsoMu24"] == 1)
    bool_list_bkg = (bool_list_bkg) & (background_branches["Nlep_valid"] > 0)

    bool_list_sig = np.ones_like(signal_branches[variable], dtype=bool)
    bool_list_sig = (bool_list_sig) & (signal_branches["triggerIsoMu24"] == 1)
    bool_list_sig = (bool_list_sig) & (signal_branches["Nlep_valid"] > 0)

    bkg_histogram, bins = np.histogram(
        background_branches[variable][bool_list_bkg],
        bins=N_BINS[variable],
        range=X_RANGE[variable],
        weights=background_branches["weight"][bool_list_bkg],
    )
    signal_histogram, _ = np.histogram(
        signal_branches[variable][bool_list_sig],
        bins=N_BINS[variable],
        range=X_RANGE[variable],
        weights=signal_branches["weight"][bool_list_sig],
    )
    fig, ax = get_canvas()

    bkg_histogram_norm = bkg_histogram / np.sum(bkg_histogram)
    signal_histogram_norm = signal_histogram / np.sum(signal_histogram)

    hep.histplot(
        signal_histogram_norm,
        bins,
        yerr=False,
        histtype="fill",
        label="Signal Evnts: " + str(np.sum(signal_histogram)),
        ax=ax,
        color="blue",
        alpha=0.4,
    )

    hep.histplot(
        signal_histogram_norm,
        bins,
        yerr=False,
        ax=ax,
        color="blue",
        linewidth=2,
    )

    hep.histplot(
        bkg_histogram_norm,
        bins,
        yerr=False,
        histtype="fill",
        label="Background Evnts: " + str(np.sum(bkg_histogram)),
        ax=ax,
        color="red",
        alpha=0.4,
    )

    hep.histplot(
        bkg_histogram_norm,
        bins,
        yerr=False,
        ax=ax,
        color="red",
        linewidth=2,
    )
    hep.cms.label(
        data="True",
        label="",
        year=2011,
        com="7",
        lumi="50",
        ax=ax,
    )

    # plot_log_variables = []

    ax.set_ylabel(r"Events / Total events")  # type: ignore
    ax.set_xlabel(variable)  # type: ignore
    y_max = max(
        np.max(bkg_histogram_norm),
        np.max(signal_histogram_norm),
    )
    if np.isnan(y_max):
        y_max = 1
    ax.set_ylim(0.0, 1.3 * y_max)  # type: ignore
    ax.set_xlim(bins[0], bins[-1])  # type: ignore
    ax.legend(frameon=False, loc="upper right", ncols=1)  # type: ignore

    output_directory = "../plots/sig_vs_bkg/"

    save_figure(fig, output_directory, variable + "_sig_vs_bkg")

    plt.close()


for var in ALL_BRANCHES:
    draw_sig_and_bg(var)
