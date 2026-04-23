import math
import sys

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot as ur
from utils.constants import N_BINS, X_RANGE
from utils.helper import get_canvas, get_histograms_ratio, save_figure

plt.style.use(hep.style.CMS)
# "Bjet1_Btag" = "MET_pt"
# "Bjet1_Btag" = "diMuon_mass"
path = "../data/proccess_tuples/"
N_max_iterations = 1
variables = [
    "bjet1_Btag",
    "bjet2_Btag",
    "bjet3_Btag",
    "bjet4_Btag",
    "N_valid_jets",
    "N_valid_b_jets",
    "weight",
    "NMuon_valid",
]


cut_min = 1.0
cut_max = 2


def build_btag_matrix(arrays):
    return np.column_stack(
        [
            arrays["bjet1_Btag"],
            arrays["bjet2_Btag"],
            arrays["bjet3_Btag"],
            arrays["bjet4_Btag"],
        ]
    )


def count_bjets_passing_cut(btag_matrix, cut):
    return np.sum(btag_matrix > cut, axis=1)


def find_best_cut(cuts_list, cut_min=cut_min, cut_max=cut_max):

    with ur.open(path + "background.root:tree_output") as file:
        bkg_branches = file.arrays(variables, library="np")
    with ur.open(path + "signal.root:tree_output") as file:
        signal_branches = file.arrays(variables, library="np")
    signal_btags = build_btag_matrix(signal_branches)
    bkg_btags = build_btag_matrix(bkg_branches)

    signal = []
    bkg_sqrt = []
    significance = []
    bins = []
    cuts = np.linspace(cut_min, cut_max, 15)
    for cut in cuts:

        print("----------------------------------------------------------------")
        print("cut:", cut)
        print("----------------------------------------------------------------")
    signal_btags = build_btag_matrix(signal_branches)
    bkg_btags = build_btag_matrix(bkg_branches)

        signal_n_valid_bjets_cut = count_bjets_passing_cut(signal_btags, cut)
        bkg_n_valid_bjets_cut = count_bjets_passing_cut(bkg_btags, cut)

        signal_n_valid_jets = signal_branches["N_valid_jets"] + (
            4 - signal_n_valid_bjets_cut
        )
        bkg_n_valid_jets = bkg_branches["N_valid_jets"] + (4 - bkg_n_valid_bjets_cut)

        bkg_bool_list = (
            (bkg_n_valid_jets >= 2)
            & (bkg_n_valid_bjets_cut >= 2)
            & (bkg_branches["NMuon_valid"] == 1)
        )
        # signal_bool_list = (signal_branches["Bjet1_Btag"] < cut_min) & (
        signal_bool_list = (
            (signal_n_valid_jets >= 2)
            & (signal_n_valid_bjets_cut >= 2)
            & (signal_branches["NMuon_valid"] == 1)
        )
        signal_events = np.sum(signal_branches["weight"][signal_bool_list])
        bkg_events = np.sum(bkg_branches["weight"][bkg_bool_list])

        print("----------------------------------------------------------------")
        print("signal = ", signal_events)
        print("bkg = ", bkg_events)
        print("S/sqrt(bkg) = ", signal_events / math.sqrt(bkg_events))
        print("cut= ", cut)
        print("----------------------------------------------------------------")
        signal.append(signal_events)
        bkg_sqrt.append(math.sqrt(bkg_events))
        if bkg_events <= 0:
            continue
        # bins.append(cut_min)
        bins.append(cut)
        significance.append(signal_events / math.sqrt(bkg_events))

    fig, ax = get_canvas()

    ax.errorbar(
        bins,
        significance,
        # ey,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=5,
        # label=labelList[j],
    )

    hep.cms.label(data="True", ax=ax, com="13.6")

    max_height = max(significance)
    best_cut = bins[significance.index(max(significance))]
    print("MAX significance: ", max(significance))
    ax.axvline(
        x=best_cut,
        # ymin=0.0,
        # ymax=max_height,
        color="red",
        linestyle="--",
        alpha=0.5,
    )

    ax.set_ylim(0.0, 1.3 * max_height)
    ax.set_xlim(bins[0], bins[-1])
    ax.set_ylabel(r"S/$\sqrt{B}$", loc="center")
    # ax.legend(frameon=False, loc="upper right")
    ax.set_xlabel("Threshold")
    output_directory = "../plots/cuts-scans/"

    save_name = "cuts_Btag"
    save_figure(fig, output_directory, save_name)

    cuts_list.append(round(best_cut, 3))


def draw_bdt_categories(cut_min=cut_min, cut_max=cut_max):

    plt.style.use(hep.style.CMS)
    cuts_list = []
    find_best_cut(cuts_list)
    cuts_list.append(cut_min)
    cuts_list.reverse()
    cuts_list.append(cut_max)

    print("Categories: ", cuts_list)
    

    plot_variable= "bjet4_Btag"
    with ur.open(path + "background.root:tree_output") as file:
        bkg_branches = file.arrays(variables, library="np")
        bkg_hist, bkg_bins = np.histogram(
            bkg_branches[plot_variable],
            bins=N_BINS[plot_variable],
            range=X_RANGE[plot_variable],
            weights=bkg_branches["weight"],
        )
    with ur.open(path + "signal.root:tree_output") as file:
        signal_branches = file.arrays(variables, library="np")

        signal_hist, signal_bins = np.histogram(
            signal_branches[plot_variable],
            bins=N_BINS[plot_variable],
            range=X_RANGE[plot_variable],
            weights=signal_branches["weight"],
        )

    fig, axs = get_canvas(True)

    hep.histplot(
        bkg_hist / np.sum(bkg_hist),
        bkg_bins,
        label="Background",
        ax=axs[0],
        stack=True,
        linewidth=2,
        color="blue",
    )

    hep.histplot(
        signal_hist / np.sum(signal_hist),
        signal_bins,
        label="Signal",
        color="red",
        linewidth=2,
        ax=axs[0],
    )

    hep.cms.label(
        data="True",
        label="",
        # year=era,
        com="13.6",
        # lumi=luminosity[era],
        ax=axs[0],
    )

    axs[0].set_ylim(
        0.0,
        1.3
        * max(
            np.max(signal_hist / np.sum(signal_hist)),
            np.max(bkg_hist / np.sum(bkg_hist)),
        ),
    )
    axs[0].set_xlim(signal_bins[0], signal_bins[-1])
    axs[0].set_ylabel(r"Events/ Total events", loc="center")
    axs[0].legend(frameon=False, loc="upper right")

    signal = []
    bkg_sqrt = []

    for category in range(len(cuts_list) - 1):
        if category != 0:
            axs[0].axvline(
                x=cuts_list[category],
                color="grey",
                linestyle="--",
                alpha=0.5,
            )
        bkg_bool_list = bkg_branches[plot_variable] < cuts_list[category]
        signal_bool_list = signal_branches[plot_variable] < cuts_list[category]

        signal_events = np.sum(signal_branches["weight"][signal_bool_list])
        bkg_events = np.sum(bkg_branches["weight"][bkg_bool_list])
        signal.append(signal_events)
        bkg_sqrt.append(math.sqrt(bkg_events))

    ratio_hist, ratio_error = get_histograms_ratio(np.array(signal), np.array(bkg_sqrt))

    hep.histplot(
        ratio_hist,
        cuts_list,
        yerr=ratio_error,
        xerr=True,
        histtype="errorbar",
        color="black",
        ax=axs[1],
    )

    print(ratio_hist)

    axs[1].set_ylabel(r"S/$\sqrt{B}$", loc="center")
    axs[1].set_ylim(0.0, 2.000)
    axs[1].set_xlim(signal_bins[0], signal_bins[-1])
    axs[1].set_xlabel("Bjet1_Btag")

    output_directory = "../plots/cuts/"
    save_name = "cut_" + plot_variable
    save_figure(fig, output_directory, save_name)


draw_bdt_categories()
