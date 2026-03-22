import uproot as ur
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt
import awkward as ak

from utils.helper import (
    get_canvas,
    save_figure,
    get_histograms_ratio,
    # clean_null_values,
    get_background_label_list,
    get_color_list,
)
from utils.constants import (
    SOURCES_LABEL,
    N_BINS,
    ALL_BRANCHES,
    X_RANGE,
    BACKGROUND_SOURCES,
    SIGNAL_SOURCES,
)


def get_histograms_from_tuple(
    sources,
    variables,
    # is_background,
):

    # if not use_puweight:
    # variables.append("pileup_weight")

    histograms_list = []
    bins_list = []

    # variable_bin = variables[0]

    tuple_path = "../data/tuples/"

    for source in sources:
        file_name = source + ".root:events"

        with ur.open(tuple_path + file_name) as file:
            branches = file.arrays(variables, library="np")

        # If branch is "array per event" (numpy object array), take the first element per event
        if isinstance(branches[variables[0]][0], (list, np.ndarray)):
            # print("before: ", signal_branches[variable])

            branches_bool= [
                True if len(arr)> 0  else False for arr in branches[variables[0]]
            ]
            branches[variables[0]] = [
                arr[0] for arr in branches[variables[0]] if len(arr) > 0
            ]
            branches["triggerIsoMu24"] = branches["triggerIsoMu24"][branches_bool]
            branches["EventWeight"] = branches["EventWeight"][branches_bool]
            branches[variables[0]] = np.array(branches[variables[0]])

        if variables[0] != "triggerIsoMu24":
            branches[variables[0]] =branches[variables[0]][branches["triggerIsoMu24"]==1] 
            if variables[0] != "EventWeight":
                branches["EventWeight"] =branches["EventWeight"][branches["triggerIsoMu24"]==1] 

        histogram, bins = np.histogram(
            branches[variables[0]],
            bins=N_BINS[variables[0]],
            range=X_RANGE[variables[0]],
            weights=(branches["EventWeight"]),
        )
        histograms_list.append(histogram)
        bins_list.append(bins)
    return histograms_list, bins_list


def draw_data_and_simul_and_ratio(
    variable,
    background_sources,
    signal_sources,
):
    plt.style.use(hep.style.CMS)

    print("*" * len("****** PLOTTING " + variable + " *****"))
    print("****** PLOTTING " + variable + " *****")
    print("*" * len("****** PLOTTING " + variable + " *****"))

    variables = [variable, "EventWeight", "triggerIsoMu24"]

    data_histogram, data_bins = get_histograms_from_tuple(
        ["data"],
        variables,
    )

    bkg_histograms_list, bkg_bins_list = get_histograms_from_tuple(
        background_sources,
        variables,
    )

    signal_histograms_list, signal_bins_list = get_histograms_from_tuple(
        signal_sources,
        variables,
    )

    fig, axs = get_canvas(True)

    hep.histplot(
        bkg_histograms_list,
        bins=bkg_bins_list[0],
        yerr=True,
        histtype="fill",
        label=get_background_label_list(background_sources),
        ax=axs[0],
        stack=True,
        color=get_color_list(len(background_sources)),
    )

    hep.histplot(
        data_histogram[0],
        data_bins[0],
        yerr=True,
        histtype="errorbar",
        label="Data",
        color="black",
        ax=axs[0],
    )

    signal_scale_factor = 1
    for source, histogram in zip(signal_sources, signal_histograms_list):
        hep.histplot(
            histogram * signal_scale_factor,
            signal_bins_list[0],
            yerr=False,
            # yerr=True,
            # label=source + " (x" + str(signal_scale_factor) + ")",
            label=SOURCES_LABEL[source],
            color="red",
            ax=axs[0],
        )

    label = ""
    hep.cms.label(
        data="True",
        label=label,
        year=2011,
        com="7",
        lumi="50",
        ax=axs[0],
    )

    axs[0].set_ylabel(r"Events")
    axs[0].set_ylim(0.1, 1000 * np.max(data_histogram))
    axs[0].set_xlim(data_bins[0][0], data_bins[0][-1])
    axs[0].set_yscale("log")
    axs[0].legend(frameon=False, loc="upper right", ncols=2)
    axs[0].tick_params(axis="x", which="both", bottom=True, top=True, labelbottom=False)

    plt.axhline(y=1, color="grey", linestyle="--", alpha=0.5)

    tot_bg_numpy_hist = np.array([])
    for i, bg_hist in enumerate(bkg_histograms_list):
        if i == 0:
            tot_bg_numpy_hist = bg_hist
        else:
            tot_bg_numpy_hist = tot_bg_numpy_hist + bg_hist

    ratio_hist, ratio_error = get_histograms_ratio(data_histogram[0], tot_bg_numpy_hist)

    hep.histplot(
        ratio_hist,
        data_bins[0],
        yerr=ratio_error,
        histtype="errorbar",
        label="data",
        color="black",
        ax=axs[1],
    )

    axs[1].set_ylabel("Data/MC", loc="center")
    axs[1].set_ylim(0.5, 1.5)
    axs[1].set_xlim(data_bins[0][0], data_bins[0][-1])
    axs[1].set_xlabel(variable)

    output_directory = "../plots/ratio/"
    output_name = variable + "_MCData_ratio"

    # output_directory = get_output_directory(variable, output_directory, variables_type)

    save_figure(fig, output_directory, output_name)

    plt.close()


for variable in ALL_BRANCHES:
    draw_data_and_simul_and_ratio(
        variable,
        BACKGROUND_SOURCES,
        SIGNAL_SOURCES,
    )
