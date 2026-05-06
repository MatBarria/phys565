import math
import sys

import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot as ur
from utils.constants import N_BINS, X_RANGE
from utils.helper import get_canvas, get_histograms_ratio, get_x_label, save_figure

plt.style.use(hep.style.CMS)
path = "../data/proccess_tuples/no_cuts/"
# path = "../data/proccess_tuples/"
N_max_iterations = 1

N_min_b_jets = 1
cut_min = 0
cut_max = 7
n_steps = 140

variables = list(
    set(
        [
            "bjet1_Btag",
            "bjet2_Btag",
            "bjet3_Btag",
            "weight",
            "ME"
            "NMuon_valid",
            "N_valid_jets_tot",
            "N_valid_b_jets",
            "triggerIsoMu24",
        ]
    )
)


def find_best_cut(cuts_list, cut_min=cut_min, cut_max=cut_max):

    with ur.open(path + "background.root:tree_output") as file:
        bkg_branches = file.arrays(variables, library="np")
    with ur.open(path + "signal.root:tree_output") as file:
        signal_branches = file.arrays(variables, library="np")

    print(bkg_branches["NMuon_valid"])
    print(signal_branches["NMuon_valid"])
    # print(bkg_branches["bjet1_Btag"])
    # print(signal_branches["bjet1_Btag"])
    print(bkg_branches["weight"])
    print(signal_branches["weight"])

    signal = []
    bkg_sqrt = []
    significance = []
    bins = []
    cuts = np.linspace(cut_min, cut_max, n_steps)
    # cuts = np.linspace(cut_max, cut_min, 400)
    for cut in cuts:

        # print("----------------------------------------------------------------")
        # print("cut:", cut)
        # print("----------------------------------------------------------------")
        # bkg_bool_list = (bkg_branches["bjet1_Btag"] < cut_min) & (
        # if cut < 4.5:
        # continue
        #
        #
        #

        bkg_n_btags = (
            (bkg_branches["bjet1_Btag"] > cut).astype(int)
            + (bkg_branches["bjet2_Btag"] > cut).astype(int)
            + (bkg_branches["bjet3_Btag"] > cut).astype(int)
        )

        signal_n_btags = (
            (signal_branches["bjet1_Btag"] > cut).astype(int)
            + (signal_branches["bjet2_Btag"] > cut).astype(int)
            + (signal_branches["bjet3_Btag"] > cut).astype(int)
        )

        bkg_bool_list = (
            (bkg_branches["NMuon_valid"] > 0)
            & (bkg_branches["N_valid_jets_tot"] >= 4)
            & (bkg_n_btags >= N_min_b_jets)
        )
        # signal_bool_list = (signal_branches["bjet1_Btag"] < cut_min) & (
        signal_bool_list = (
            (signal_branches["NMuon_valid"] > 0)
            & (signal_branches["triggerIsoMu24"] == 1)
            & (signal_branches["N_valid_jets_tot"] >= 4)
            & (signal_n_btags >= N_min_b_jets)
        )

        signal_events = np.sum(signal_branches["weight"][signal_bool_list])
        bkg_events = np.sum(bkg_branches["weight"][bkg_bool_list])
        # print(signal_branches["weight"][signal_bool_list])
        # print(bkg_branches["weight"][bkg_bool_list])
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
        significance.append(signal_events / math.sqrt(bkg_events + signal_events))

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

    # hep.cms.label(data="True", ax=ax, com="13.6")

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

    save_name = "cuts_" + "bjet1_Btag"
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

    with ur.open(path + "background.root:tree_output") as file:
        bkg_branches = file.arrays(variables, library="np")

        bkg_bool_list = (
            # (bkg_branches["bjet1_Btag"] < cu)
            (bkg_branches["NMuon_valid"] > 0)
            & (bkg_branches["N_valid_jets_tot"] >= 4)
        )

        bkg_hist, bkg_bins = np.histogram(
            bkg_branches["bjet1_Btag"][bkg_bool_list],
            # bkg_branches["bjet1_Btag"],
            bins=N_BINS["bjet1_Btag"],
            range=X_RANGE["bjet1_Btag"],
            weights=bkg_branches["weight"][bkg_bool_list],
        )
    with ur.open(path + "signal.root:tree_output") as file:
        signal_branches = file.arrays(variables, library="np")
        # signal_bool_list = (signal_branches["bjet1_Btag"] < cut_min) & (

        signal_bool_list = (
            # (signal_branches["bjet1_Btag"] < cut)
            (signal_branches["NMuon_valid"] > 0)
            & (signal_branches["N_valid_jets_tot"] >= 4)
            & (signal_branches["triggerIsoMu24"] == 1)
        )

        signal_hist, signal_bins = np.histogram(
            signal_branches["bjet1_Btag"][signal_bool_list],
            # signal_branches["bjet1_Btag"],
            bins=N_BINS["bjet1_Btag"],
            range=X_RANGE["bjet1_Btag"],
            weights=signal_branches["weight"][signal_bool_list],
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
        data=False,
        label="",
        # year=era,
        com="7.0",
        lumi=50,
        ax=axs[0],
    )

    axs[0].axvline(
        x=cuts_list[1],
        color="grey",
        linestyle="--",
        # alpha=0.8,
        label="Threshold",
    )

    axs[0].legend(frameon=False, loc="upper right")
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

    signal = []
    bkg_sqrt = []

    bkg_bool_list = (
        # (bkg_branches["bjet1_Btag"] < cu)
        (
            (bkg_branches["bjet1_Btag"] > cuts_list[1])
            | (bkg_branches["bjet2_Btag"] > cuts_list[1])
            | (bkg_branches["bjet3_Btag"] > cuts_list[1])
        )
        & (bkg_branches["NMuon_valid"] > 0)
        & (bkg_branches["N_valid_jets_tot"] >= 4)
        & (bkg_branches["N_valid_b_jets"] >= N_min_b_jets)
    )

    signal_bool_list = (
        (
            (signal_branches["bjet1_Btag"] > cuts_list[1])
            | (signal_branches["bjet2_Btag"] > cuts_list[1])
            | (signal_branches["bjet3_Btag"] > cuts_list[1])
        )
        & (signal_branches["NMuon_valid"] > 0)
        & (signal_branches["N_valid_jets_tot"] >= 4)
        & (signal_branches["N_valid_b_jets"] >= N_min_b_jets)
        & (signal_branches["triggerIsoMu24"] == 1)
    )

    # bkg_bool_list = bkg_branches["bjet1_Btag"] > cuts_list[1]s_list[category]
    # signal_bool_list = signal_branches["bjet1_Btag"] > cuts_list[1]s_list[category]

    signal_events = np.sum(signal_branches["weight"][signal_bool_list])
    bkg_events = np.sum(bkg_branches["weight"][bkg_bool_list])
    signal.append(signal_events)
    bkg_sqrt.append(math.sqrt(bkg_events))

    bkg_bool_list = (
        (
            (bkg_branches["bjet1_Btag"] > cuts_list[1])
            | (bkg_branches["bjet2_Btag"] > cuts_list[1])
            | (bkg_branches["bjet3_Btag"] > cuts_list[1])
        )
        & (bkg_branches["NMuon_valid"] > 0)
        & (bkg_branches["N_valid_jets_tot"] >= 4)
        & (bkg_branches["N_valid_b_jets"] >= N_min_b_jets)
    )

    signal_bool_list = (
        (
            (signal_branches["bjet1_Btag"] > cuts_list[1])
            | (signal_branches["bjet2_Btag"] > cuts_list[1])
            | (signal_branches["bjet3_Btag"] > cuts_list[1])
        )
        & (signal_branches["NMuon_valid"] > 0)
        & (signal_branches["N_valid_jets_tot"] >= 4)
        & (signal_branches["N_valid_b_jets"] >= N_min_b_jets)
        & (signal_branches["triggerIsoMu24"] == 1)
    )

    signal_events = np.sum(signal_branches["weight"][signal_bool_list])
    bkg_events = np.sum(bkg_branches["weight"][bkg_bool_list])
    signal.append(signal_events)
    bkg_sqrt.append(math.sqrt(bkg_events))

    ratio_hist, ratio_error = get_histograms_ratio(
        np.array(signal), np.sqrt(np.array(signal) + np.array(bkg_sqrt))
    )

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

    axs[1].set_ylabel(r"S/$\sqrt{S + B}$", loc="center")
    axs[1].set_ylim(0.0, 1.2 * np.max(ratio_hist))
    axs[1].set_xlim(signal_bins[0], signal_bins[-1])
    axs[1].set_xlabel(get_x_label("bjet1_Btag"))

    output_directory = "../plots/cuts/"
    save_name = "cut" + "btagging"
    save_figure(fig, output_directory, save_name)

def draw_bdt_categoriesv2(cut_min=cut_min, cut_max=cut_max):

    plt.style.use(hep.style.CMS)

    cuts_list = []
    find_best_cut(cuts_list)
    cuts_list.append(cut_min)
    cuts_list.reverse()
    cuts_list.append(cut_max)

    print("Categories: ", cuts_list)

    with ur.open(path + "background.root:tree_output") as file:
        bkg_branches = file.arrays(variables, library="np")

    with ur.open(path + "signal.root:tree_output") as file:
        signal_branches = file.arrays(variables, library="np")

    # ------------------------------------------------------------
    # Event selections
    # ------------------------------------------------------------
    bkg_event_sel = (
        (bkg_branches["NMuon_valid"] > 0)
        & (bkg_branches["N_valid_jets_tot"] >= 4)
        # & (bkg_branches["N_valid_b_jets"] >= N_min_b_jets)
    )

    signal_event_sel = (
        (signal_branches["NMuon_valid"] > 0)
        & (signal_branches["N_valid_jets_tot"] >= 4)
        # & (signal_branches["N_valid_b_jets"] >= N_min_b_jets)
        & (signal_branches["triggerIsoMu24"] == 1)
    )

    # ------------------------------------------------------------
    # Concatenate the b-tag scores of the three jets
    # One single distribution, not three curves
    # ------------------------------------------------------------
    bkg_btag_scores = np.concatenate(
        [
            bkg_branches["bjet1_Btag"][bkg_event_sel],
            bkg_branches["bjet2_Btag"][bkg_event_sel],
            bkg_branches["bjet3_Btag"][bkg_event_sel],
        ]
    )

    signal_btag_scores = np.concatenate(
        [
            signal_branches["bjet1_Btag"][signal_event_sel],
            signal_branches["bjet2_Btag"][signal_event_sel],
            signal_branches["bjet3_Btag"][signal_event_sel],
        ]
    )

    # Repeat the event weights once per jet
    bkg_weights = np.concatenate(
        [
            bkg_branches["weight"][bkg_event_sel],
            bkg_branches["weight"][bkg_event_sel],
            bkg_branches["weight"][bkg_event_sel],
        ]
    )

    signal_weights = np.concatenate(
        [
            signal_branches["weight"][signal_event_sel],
            signal_branches["weight"][signal_event_sel],
            signal_branches["weight"][signal_event_sel],
        ]
    )

    bkg_hist, bkg_bins = np.histogram(
        bkg_btag_scores,
        bins=N_BINS["bjet1_Btag"],
        range=X_RANGE["bjet1_Btag"],
        weights=bkg_weights,
    )

    signal_hist, signal_bins = np.histogram(
        signal_btag_scores,
        bins=N_BINS["bjet1_Btag"],
        range=X_RANGE["bjet1_Btag"],
        weights=signal_weights,
    )

    fig, axs = get_canvas(True)

    # Avoid division by zero
    bkg_hist_norm = bkg_hist / np.sum(bkg_hist) if np.sum(bkg_hist) > 0 else bkg_hist
    signal_hist_norm = signal_hist / np.sum(signal_hist) if np.sum(signal_hist) > 0 else signal_hist

    hep.histplot(
        bkg_hist_norm,
        bkg_bins,
        label="Background",
        ax=axs[0],
        linewidth=2,
        color="blue",
    )

    hep.histplot(
        signal_hist_norm,
        signal_bins,
        label="Signal",
        color="red",
        linewidth=2,
        ax=axs[0],
    )

    hep.cms.label(
        data=False,
        label="",
        com="7.0",
        lumi=50,
        ax=axs[0],
    )

    axs[0].axvline(
        x=cuts_list[1],
        color="grey",
        linestyle="--",
        label="Threshold",
    )

    axs[0].legend(frameon=False, loc="upper right")

    axs[0].set_ylim(
        0.0,
        1.3 * max(np.max(signal_hist_norm), np.max(bkg_hist_norm)),
    )

    axs[0].set_xlim(signal_bins[0], signal_bins[-1])
    axs[0].set_ylabel(r"Jets / Total jets", loc="center")

    # ------------------------------------------------------------
    # Compute significance after the chosen cut
    # Event-level cut: at least one of the three jets passes threshold
    # ------------------------------------------------------------
    cut = cuts_list[1]

    bkg_cut_sel = (
        (
            (bkg_branches["bjet1_Btag"] > cut)
            | (bkg_branches["bjet2_Btag"] > cut)
            | (bkg_branches["bjet3_Btag"] > cut)
        )
        & bkg_event_sel
    )

    signal_cut_sel = (
        (
            (signal_branches["bjet1_Btag"] > cut)
            | (signal_branches["bjet2_Btag"] > cut)
            | (signal_branches["bjet3_Btag"] > cut)
        )
        & signal_event_sel
    )

    S = np.sum(signal_branches["weight"][signal_cut_sel])
    B = np.sum(bkg_branches["weight"][bkg_cut_sel])

    Z = S / np.sqrt(S + B) if (S + B) > 0 else 0.0

    print("Cut:", cut)
    print("Signal:", S)
    print("Background:", B)
    print("S/sqrt(S+B):", Z)

    # For the lower panel, draw one point/bin showing the significance
    ratio_hist = np.array([Z])
    ratio_bins = np.array([cut_min, cut_max])

    hep.histplot(
        ratio_hist,
        ratio_bins,
        histtype="step",
        color="black",
        ax=axs[1],
    )

    axs[1].set_ylabel(r"S/$\sqrt{S + B}$", loc="center")
    axs[1].set_ylim(0.0, 1.2 * Z if Z > 0 else 1.0)
    axs[1].set_xlim(signal_bins[0], signal_bins[-1])
    axs[1].set_xlabel("b-tag discriminator")

    output_directory = "../plots/cuts/"
    save_name = "cut_btagging_all_three_jets"
    save_figure(fig, output_directory, save_name)




draw_bdt_categoriesv2()
