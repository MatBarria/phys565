import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot as ur
from utils.constants import (
    ALL_BRANCHES,
    BACKGROUND_SOURCES,
    N_BINS,
    SIGNAL_SOURCES,
    SOURCES_LABEL,
    USEFUL_BRANCHES,
    X_RANGE,
)
from utils.helper import (
    clean_null_values,
    get_background_label_list,
    get_canvas,
    get_color_list,
    get_histograms_ratio,
    save_figure,
)

final_bkg_sources = []


def get_histograms_from_tuple(
    sources,
    variables,
    # is_background,
    use_dimuon_mass_cut=False,
):

    # if not use_puweight:
    # variables.append("pileup_weight")

    histograms_list = []
    bins_list = []

    # variable_bin = variables[0]

    tuple_path = "../data/proccess_tuples/"

    for source in sources:
        file_name = source + "_tuples.root:tree_output"

        with ur.open(tuple_path + file_name) as file:
            branches = file.arrays(variables, library="np")

        bool_list = np.ones_like(branches[variables[0]], dtype=bool)

        if variables[0] != "triggerIsoMu24":
            bool_list = (bool_list) & (branches["triggerIsoMu24"] == 1)
        # if variables[0] != "NMuon_valid":
        bool_list = (bool_list) & (branches["NMuon_valid"] == 1)
        # bool_list = (bool_list) & (branches["N_valid_jets"] >= 2)
        # bool_list = (bool_list) & (branches["N_valid_b_jets"] >= 2)
        if use_dimuon_mass_cut:
            bool_list = (bool_list) & (branches["diMuon_mass"] < 20)

        for var in variables:
            # if var in ["top_hadronic_mass_1", "top_hadronic_mass_2"]:
                # continue
            branches[var] = branches[var][bool_list]
        if variables[0] != "diMuon_mass":
            clean_null_values(branches, variables)

        histogram, bins = np.histogram(
            branches[variables[0]],
            bins=N_BINS[variables[0]],
            range=X_RANGE[variables[0]],
            weights=(branches["weight"]),
        )
        # if np.sum(histogram) ==0:
        # sources.remove(source)
        # continue;
        histogram[histogram == 0] = 0.00000001
        histograms_list.append(histogram)
        bins_list.append(bins)
        # print("variable: ", variables[0])
        # print("source : ", source)
        # print("histo : ", histogram)
    return histograms_list, bins_list


def draw_data_and_simul_and_ratio(
    variable,
    background_sources,
    signal_sources,
    use_dimuon_mass_cut=False,
):
    plt.style.use(hep.style.CMS)

    print("*" * len("****** PLOTTING " + variable + " *****"))
    print("****** PLOTTING " + variable + " *****")
    print("*" * len("****** PLOTTING " + variable + " *****"))

    variables = [
        variable,
        "weight",
        "triggerIsoMu24",
        "diMuon_mass",
        "NMuon_valid",
        "Nlep_valid",
        "N_valid_jets",
        "N_valid_b_jets",
    ]
    variables=list(set(variables))
    variables.insert(0, variables.pop(variables.index(variable)))

    # print(variables)
    # if variable == "top_hadronic_mass":
        # variables.append("top_hadronic_mass_1")
        # variables.append("top_hadronic_mass_2")


    print(variables)
    
    data_histogram, data_bins = get_histograms_from_tuple(
        ["data"],
        variables,
        use_dimuon_mass_cut,
    )

    bkg_histograms_list, bkg_bins_list = get_histograms_from_tuple(
        background_sources,
        variables,
        use_dimuon_mass_cut,
    )

    signal_histograms_list, signal_bins_list = get_histograms_from_tuple(
        signal_sources,
        variables,
        use_dimuon_mass_cut,
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
        color=get_color_list(len(bkg_histograms_list)),
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
    axs[0].set_ylim(
        0.1,
        1000 * np.max(np.concatenate([data_histogram[0], signal_histograms_list[0]])),
    )
    axs[0].set_xlim(data_bins[0][0], data_bins[0][-1])
    # axs[0].set_yscale("log")
    axs[0].set_ylim(0.,1.4 * np.max(data_histogram))
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

    output_directory = "../plots/ratio-proccessed/"
    output_name = variable + "_MCData_ratio"

    # output_directory = get_output_directory(variable, output_directory, variables_type)

    save_figure(fig, output_directory, output_name)

    plt.close()


for variable in USEFUL_BRANCHES:
    draw_data_and_simul_and_ratio(
        variable,
        BACKGROUND_SOURCES,
        SIGNAL_SOURCES,
    )
