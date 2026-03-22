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

    tuple_path = "../data/tuples/"

    print("PLOTTING: ", variable)
    background_tree = ur.open(tuple_path + "background.root:events")
    background_branches = background_tree.arrays([variable, "EventWeight"], library="np")  # type: ignore

    signal_tree = ur.open(tuple_path + "signal.root:events")
    signal_branches = signal_tree.arrays([variable, "EventWeight", "triggerIsoMu24"], library="np")  # type: ignore
    if isinstance(signal_branches[variable][0], (list, np.ndarray)):
        print("array of arrays")
        # print("before: ", signal_branches[variable])

        signal_branches_bool = [
            True if len(arr)> 0  else False for arr in signal_branches[variable]
        ]
        signal_branches[variable] = [
            arr[0] for arr in signal_branches[variable] if len(arr) > 0
        ]
        background_branches_bool = [
            True if len(arr)> 0  else False for arr in background_branches[variable]
        ]
        background_branches[variable] = [
            arr[0] for arr in background_branches[variable] if len(arr) > 0
        ]
        signal_branches[variable] = np.array(signal_branches[variable])
        background_branches[variable] = np.array(background_branches[variable])
        signal_branches["triggerIsoMu24"] = signal_branches["triggerIsoMu24"][signal_branches_bool]
        signal_branches["EventWeight"] = signal_branches["EventWeight"][signal_branches_bool]
        background_branches["EventWeight"] = background_branches["EventWeight"][background_branches_bool]
        # print("after: ", signal_branches[variable])
       
    else:
        print("array of floats")

    if variable != "triggerIsoMu24":
        signal_branches["EventWeight"] = signal_branches["EventWeight"][
            signal_branches["triggerIsoMu24"] == 1
        ]
        signal_branches[variable] = signal_branches[variable][
            signal_branches["triggerIsoMu24"] == 1
        ]

    if variable[0] != "N":
        signal_branches["EventWeight"] = signal_branches["EventWeight"][
            signal_branches[variable] != 0
        ]
        signal_branches[variable] = signal_branches[variable][
            signal_branches[variable] != 0
        ]
        background_branches["EventWeight"] = background_branches["EventWeight"][
            background_branches[variable] != 0
        ]
        background_branches[variable] = background_branches[variable][
            background_branches[variable] != 0
        ]

    bkg_histogram, bins = np.histogram(
        background_branches[variable],
        bins=N_BINS[variable],
        range=X_RANGE[variable],
        weights=background_branches["EventWeight"],
    )
    signal_histogram, _ = np.histogram(
        signal_branches[variable],
        bins=N_BINS[variable],
        range=X_RANGE[variable],
        weights=signal_branches["EventWeight"],
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
