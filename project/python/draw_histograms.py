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


def draw_data_and_simul_and_ratio(
    variable,
    background_sources,
    signal_sources,
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

    data_histogram, _, data_bins = get_histograms_from_tuple(
        ["data"],
        variables,
        use_dimuon_mass_cut,
        use_permutation_weight,
    )

    bkg_histograms_list, bkg_histogram_error_list, bkg_bins_list = (
        get_histograms_from_tuple(
            background_sources,
            variables,
            use_dimuon_mass_cut,
            use_permutation_weight,
        )
    )

    signal_histograms_list, signal_histogram_error_list, signal_bins_list = (
        get_histograms_from_tuple(
            signal_sources,
            variables,
            use_dimuon_mass_cut,
            use_permutation_weight,
        )
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

    # Total MC statistical uncertainty band
    total_mc = np.sum(np.array(bkg_histograms_list), axis=0)
    total_mc_var = np.sum(np.array(bkg_histogram_error_list) ** 2, axis=0)
    total_mc_err = np.sqrt(total_mc_var)

    # If you want to include signal in the MC uncertainty band:

    bin_edges = bkg_bins_list[0]

    axs[0].stairs(
        total_mc + total_mc_err,
        bin_edges,
        baseline=total_mc - total_mc_err,
        fill=True,
        hatch="////",
        alpha=0.3,
        label="MC bkg stat. unc.",
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

    total_mc_sig = np.sum(np.array(signal_histograms_list), axis=0)
    total_mc_var_sig = np.sum(np.array(signal_histogram_error_list) ** 2, axis=0)
    total_mc_err_sig = np.sqrt(total_mc_var)

    total_mc_with_signal = total_mc + total_mc_sig
    total_mc_with_signal_err = np.sqrt(total_mc_err**2 + signal_histogram_error_list[0] ** 2)

    lower = total_mc_sig - total_mc_err_sig
    upper = total_mc_sig + total_mc_err_sig

    axs[0].fill_between(
        signal_bins_list[0],
        np.r_[lower, lower[-1]],
        np.r_[upper, upper[-1]],
        step="post",
        alpha=0.25,
        label="MC signal stat. unc.",
        color="red",
    )

    label = ""
    hep.cms.label(
        data=True,
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

    # if not use_permutation_weight and variable in [
    # "N_valid_jets_tot",
    # "N_valid_b_jets",
    # ]:

    # draw_fitted_xs(axs[0], variable, 506.941, 32.12, 1888.65, 49.1959)

    if variable == "top_hadronic_mass_2":
        crystal_ball_params = [
            163.753,
            1.42987,  # mean, mean_err
            23.3245,
            1.06441,  # sigma, sigma_err
            -0.860257,
            0.156361,  # alpha, alpha_err
            2.75,
            1.20396,  # n, n_err
            229.696,
            0.0,  # sumW, sumW_err
            0,
            0.0,  # fit_status
            3,
            0.0,  # cov_quality
            1152.91,
            0.0,  # minNll
        ]

        draw_fit(
            axs=axs[0],
            parameters=crystal_ball_params,
            mass_min=X_RANGE[variable][0],
            mass_max=X_RANGE[variable][1],
            n_bins=N_BINS[variable],
            label="Crystal Ball",
        )

    axs[0].set_ylabel(r"Events")
    if y_label != " ":
        axs[0].set_ylabel(y_label)
    axs[0].set_ylim(
        0.1,
        1000 * np.max(np.concatenate([data_histogram[0], signal_histograms_list[0]])),
    )
    axs[0].set_xlim(data_bins[0][0], data_bins[0][-1])
    # axs[0].set_yscale("log")
    axs[0].set_ylim(0.0, 1.4 * np.max(data_histogram))
    axs[0].legend(frameon=False, loc="upper right", ncols=2, fontsize=15)
    axs[0].tick_params(axis="x", which="both", bottom=True, top=True, labelbottom=False)

    plt.axhline(y=1, color="grey", linestyle="--", alpha=0.5)

    tot_bg_numpy_hist = np.array([])
    for i, bg_hist in enumerate(bkg_histograms_list):
        if i == 0:
            tot_bg_numpy_hist = bg_hist
        else:
            tot_bg_numpy_hist = tot_bg_numpy_hist + bg_hist

    tot_bg_numpy_hist = tot_bg_numpy_hist + signal_histograms_list[0]

    ratio_hist, ratio_error = get_histograms_ratio(
        data_histogram[0], tot_bg_numpy_hist, mc_err=total_mc_with_signal_err
    )

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
    axs[1].set_ylim(0.3, 2.5)
    axs[1].set_xlim(data_bins[0][0], data_bins[0][-1])
    axs[1].set_xlabel(get_x_label(variable))

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

    output_directory = "../plots/ratio-proccessed/"

    output_name = variable + "_MCData_ratio"
    if use_permutation_weight:
        output_name += "_per_weight"

    # output_directory = get_output_directory(variable, output_directory, variables_type)

    save_figure(fig, output_directory, output_name)

    plt.close()


# for variable in USEFUL_BRANCHES:
# draw_data_and_simul_and_ratio(
# variable,
# BACKGROUND_SOURCES,
# SIGNAL_SOURCES,
# )
# draw_data_and_simul_and_ratio(
# "N_valid_jets_tot",
# BACKGROUND_SOURCES,
# SIGNAL_SOURCES,
# )
# draw_data_and_simul_and_ratio(
# "N_valid_b_jets",
# BACKGROUND_SOURCES,
# SIGNAL_SOURCES,
# )
draw_data_and_simul_and_ratio(
    "top_hadronic_mass",
    BACKGROUND_SOURCES,
    SIGNAL_SOURCES,
    use_permutation_weight=True,
    y_label=r"Sum of permutations weights/ 15 GeV",
)
# draw_data_and_simul_and_ratio(
# "top_hadronic_mass",
# BACKGROUND_SOURCES,
# SIGNAL_SOURCES,
# # use_permutation_weight=True,
# y_label=r"Permutations / 5 GeV",
# )
draw_data_and_simul_and_ratio(
    "top_hadronic_mass_reco",
    BACKGROUND_SOURCES,
    SIGNAL_SOURCES,
    # use_permutation_weight=True,
    y_label=r"Permutations / 15 GeV",
)
draw_data_and_simul_and_ratio(
    "W_hadronic_mass_reco",
    BACKGROUND_SOURCES,
    SIGNAL_SOURCES,
    # use_permutation_weight=true,
    y_label=r"Permutations / 10 gev",
)
draw_data_and_simul_and_ratio(
    "W_hadronic_mass_reco",
    BACKGROUND_SOURCES,
    SIGNAL_SOURCES,
    use_permutation_weight=True,
    y_label=r"Sum of permutations weights/ 10 gev",
)
# draw_data_and_simul_and_ratio(
# "W_leptonic_mass_reco",
# BACKGROUND_SOURCES,
# SIGNAL_SOURCES,
# use_permutation_weight=True,
# )
# draw_data_and_simul_and_ratio(
# "NMuon_valid",
# BACKGROUND_SOURCES,
# SIGNAL_SOURCES,
# use_permutation_weight=True,
# )
# draw_data_and_simul_and_ratio(
# "chi2",
# BACKGROUND_SOURCES,
# SIGNAL_SOURCES,
# use_permutation_weight=False,
# )
