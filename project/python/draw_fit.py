import matplotlib
import pandas as pd

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
    draw_fit,
    draw_fitted_xs,
    get_background_label_list,
    get_canvas,
    get_color_list,
    get_histograms_ratio,
    get_x_label,
    save_figure,
)

final_bkg_sources = []


def get_histograms_from_tuple(
    sources,
    variables,
    # is_background,
    use_dimuon_mass_cut=False,
    use_permutation_weight=False,
):
    # if not use_puweight:
    # variables.append("pileup_weight")

    histograms_list = []
    bins_list = []
    errors_list = []

    # variable_bin = variables[0]

    # tuple_path = "../data/proccess_tuples/no_cuts/"
    # if use_permutation_weight:
    tuple_path = "../data/proccess_tuples/"

    for source in sources:
        file_name = source + "_tuples_JES1p00.root:tree_output"

        with ur.open(tuple_path + file_name) as file:
            branches = file.arrays(variables, library="np")

        # if use_permutation_weight:
        # branches["weight"] = (
        # branches["weight"]
        # * branches["permutation_weight"]
        # / branches["permutation_weight_sum"]
        # )
        weight_var = "weight"
        if use_permutation_weight:
            weight_var = "total_weight"

        bool_list = np.ones_like(branches[variables[0]], dtype=bool)

        if variables[0] != "triggerIsoMu24":
            bool_list = (bool_list) & (branches["triggerIsoMu24"] == 1)
        if variables[0] != "NMuon_valid":
            bool_list = (bool_list) & (branches["NMuon_valid"] > 0)
        bool_list = (bool_list) & (branches["N_valid_b_jets"] >= 1)
        # if variables[0] != "N_valid_jets_tot" and variables[0] != "N_valid_b_jets":
        bool_list = (bool_list) & (branches["N_valid_jets_tot"] >= 4)

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
            weights=(branches[weight_var]),
        )

        hist_variance, _ = np.histogram(
            branches[variables[0]],
            bins=N_BINS[variables[0]],
            range=X_RANGE[variables[0]],
            weights=(branches[weight_var] ** 2),
        )

        # if np.sum(histogram) ==0:
        # sources.remove(source)
        # continue;
        histogram[histogram == 0] = 0.00000001
        hist_variance[hist_variance == 0] = 0.00000001
        hist_error = np.sqrt(hist_variance)
        histograms_list.append(histogram)
        errors_list.append(hist_error)
        bins_list.append(bins)
        # print("variable: ", variables[0])
        # print("source : ", source)
        # print("histo : ", histogram)
    return histograms_list, errors_list, bins_list


def draw_distribution_and_fit(
    variable,
    sample,
    use_dimuon_mass_cut=False,
    use_permutation_weight=False,
    y_label=" ",
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
        "N_valid_jets_tot",
        "N_valid_b_jets",
        "permutation_weight",
        # "total_weight",
    ]
    if use_permutation_weight:
        variables.append("total_weight"),
        variables.append("total_weight_norm"),
    variables = list(set(variables))
    variables.insert(0, variables.pop(variables.index(variable)))

    # print(variables)
    fig, ax = get_canvas()
    axs = [ax]

    data_histogram, data_histogram_error, data_bins = get_histograms_from_tuple(
        [sample],
        variables,
        use_dimuon_mass_cut,
        use_permutation_weight,
    )
    legend_label = r"$t\bar{t}$ MC sample"
    if sample == "data":
        legend_label = r"data"


    hep.histplot(
        data_histogram[0],
        data_bins[0],
        # yerr=True,
        yerr=data_histogram_error[0],
        histtype="errorbar",
        label=legend_label,
        color="black",
        ax=axs[0],
    )

    is_data = False
    if sample == "data":
        is_data = True
    
    hep.cms.label(
        data=is_data,
        label="",
        com="7.0",
        lumi=None,  # suppress automatic lumi
        ax=axs[0],
    )
    # Add lumi manually at the standard top-right position
    axs[0].annotate(
        r"50 pb$^{-1}$               ",
        xy=(1, 1),
        xycoords="axes fraction",
        ha="right",
        va="bottom",
        fontsize=27,
    )

    #### DATA ####
    crystal_ball_params = {
        "data": [
            167.584,
            9.94351,  # mean, mean_err
            24.5185,
            5.61156,  # sigma, sigma_err
            -0.76271,
            0.988846,  # alpha, alpha_err
            2.59427,
            6.83795,  # n, n_err
            99.7227,
            0.0,  # sumW, sumW_err
            0,
            0.0,  # fit_status
            3,
            0.0,  # cov_quality
            515.61,
            0.0,  # minNll
        ],
        "ttbar": [
            165.358,
            1.51148,  # mean, mean_err
            21.037,
            1.05112,  # sigma, sigma_err
            -0.949964,
            0.177869,  # alpha, alpha_err
            1.97993,
            0.732002,  # n, n_err
            110.242,
            0.0,  # sumW, sumW_err
            0,
            0.0,  # fit_status
            3,
            0.0,  # cov_quality
            552.904,
            0.0,  # minNll
        ],
    }

    draw_fit(
        axs=axs[0],
        parameters=crystal_ball_params[sample],
        mass_min=X_RANGE[variable][0],
        mass_max=X_RANGE[variable][1],
        n_bins=N_BINS[variable],
        label="Crystal Ball: Best fit parameters",
    )

    axs[0].set_ylabel(r"Events")
    if y_label != " ":
        axs[0].set_ylabel(y_label)
    axs[0].set_ylim(
        0.1,
        1000 * np.max(data_histogram[0]),
    )
    axs[0].set_xlim(data_bins[0][0], data_bins[0][-1])
    # axs[0].set_yscale("log")
    axs[0].set_ylim(0.0, 1.4 * np.max(data_histogram))
    axs[0].legend(frameon=False, loc="upper right")
    axs[0].tick_params(axis="x", which="both", bottom=True, top=True, labelbottom=False)
    axs[0].set_xlabel(get_x_label(variable))

    # mass_reco = 172.46875
    # if "W" in variable:
    # mass_reco = 80
    # axs[0].axvline(
    # x=mass_reco,
    # # ymin=0.0,
    # # ymax=max_height,
    # color="red",
    # linestyle="--",
    # alpha=0.5,
    # )

    output_directory = "../plots/fit/"

    output_name = sample + "_fit"
    if use_permutation_weight:
        output_name += "_per_weight"

    # output_directory = get_output_directory(variable, output_directory, variables_type)

    save_figure(fig, output_directory, output_name)

    plt.close()


draw_distribution_and_fit(
    "top_hadronic_mass",
    "data",
    use_permutation_weight=True,
    y_label=r"Sum of permutations weights/ 15 GeV",
)
draw_distribution_and_fit(
    "top_hadronic_mass",
    "ttbar",
    use_permutation_weight=True,
    y_label=r"Sum of permutations weights/ 15 GeV",
)
