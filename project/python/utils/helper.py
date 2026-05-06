import os

import matplotlib.pyplot as plt
import numpy as np
import uproot as ur

from .constants import ALL_BRANCHES_ORIGINAL, COLORS, SOURCES_LABEL, X_LABEL


def get_canvas(draw_ratio=False):
    """
    Creates a matplotlib figure and axes for plotting.

    Parameters:
    ----------
    draw_ratio : bool, optional
        If True, creates a figure with two subplots in a vertical layout
        for ratio plots (default is False).

    Returns:
    -------
    fig : matplotlib.figure.Figure
        The created figure for plotting.
    axs : matplotlib.axes.Axes or array of Axes
        The axes for plotting. If draw_ratio is False, a single Axes object
        is returned; otherwise, an array of two Axes objects.
    """
    if draw_ratio:
        fig, axs = plt.subplots(2, 1, height_ratios=[10, 2])
        fig.subplots_adjust(hspace=0.1)
        return fig, axs

    fig, axs = plt.subplots(1, 1)
    return fig, axs


def save_figure(fig, outputDirectory, name):
    """
    Saves a matplotlib figure in PDF and PNG.

    Parameters:
    ----------
    fig : matplotlib.figure.Figure
        The figure to be saved.
    outputDirectory : str
        Directory where the figure will be saved.
    name : str
        The base filename for the saved figure.

    Returns:
    -------
    None
    """

    os.makedirs(outputDirectory + "pdf/", exist_ok=True)
    os.makedirs(outputDirectory + "png/", exist_ok=True)
    fig.savefig(outputDirectory + "pdf/" + name + ".pdf", bbox_inches="tight")
    fig.savefig(outputDirectory + "png/" + name + ".png", bbox_inches="tight", dpi=300)
    # fig.savefig(outputDirectory + name + ".pdf")
    # fig.savefig(outputDirectory + "png/"name + ".png", dpi=300)
    print(outputDirectory + name + " Has been created")


def get_histograms(variable, sources, era):
    """
    Loads histograms from ROOT files using the uproot library.

    Parameters:
    ----------
    variable : str
        The name of the variable or histogram to retrieve.
    sources : list of str
        List of proccess sources (e.g., "DY50to120", "TTto2L2Q", etc) to be loaded.
    era : str
        A string representing the data-taking period ("2022", "2022EE", etc.).

    Returns:
    -------
    histograms : List of ROOT histograms.
    """
    histograms = []
    for source in sources:
        with ur.open("../root_io/" + source + "_" + era + "_histograms.root") as file:
            histograms.append(file[variable])

    return histograms


def get_histograms_ratio(numerator_histogram, denominator_histogram):
    """
    Calculates the ratio of two histograms with error propagation.

    Parameters:
    ----------
    numerator_histogram : numpy.ndarray
        The histogram values for the numerator.
    denominator_histogram : numpy.ndarray
        The histogram values for the denominator.

    Returns:
    -------
    ratio : numpy.ndarray
        Numpy array with the ratio of the two histograms.
    error : numpy.ndarray
        Numpy array with the erros of ratio of the two histograms.

    Notes:
    ------
    If any elements in the denominator are zero, they are excluded from the calculation.
    """
    ratio = np.divide(
        numerator_histogram, denominator_histogram, where=(denominator_histogram != 0)
    )
    error = np.divide(
        numerator_histogram * np.sqrt(denominator_histogram)
        + denominator_histogram * np.sqrt(numerator_histogram),
        np.power(denominator_histogram, 2),
        where=(denominator_histogram != 0),
    )
    if len(error[error < 0.0]):
        msg = (
            "Unexpected negative ratio-error values found."
            "Setting them to zero.\n > Re-run. If the value changes, "
            "it might have to do with the minimum subnormal number bug!"
        )
        print(msg)
        print("Error; length")
        print(error[error < 0.0], len(error[error < 0.0]))
        print("Ratio values")
        print(ratio[np.where(error < 0.0)])
        error[error < 0.0] = 0.0

    return ratio, error


def get_output_directory(variable, base_directory, variables_dic):
    for var_type, var_list in variables_dic.items():
        if variable in var_list:
            return base_directory + var_type + "/"
    return base_directory


def clean_null_values(branches, variables):

    null_value = -1
    if "pt" in variables[0] or "Pt" in variables[0]:
        null_value = 1
    if variables[0] in ALL_BRANCHES_ORIGINAL:
        null_value = 0

    bool_list = branches[variables[0]] != null_value
    for var in variables:
        branches[var] = branches[var][bool_list]


def get_background_label_list(background_sources):
    labels_list = []
    for source in background_sources:
        labels_list.append(SOURCES_LABEL[source])
    return labels_list


def get_color_list(number_of_histograms):

    color_list = []
    for i in range(number_of_histograms):
        color_list.append(COLORS[i])

    return color_list


def draw_fitted_xs(ax, variable, n_sig, n_sig_sigma, n_bkg, n_bkg_sigma):

    sig_fractions_2d = (
        (0.060687332, 0),
        (0.17761989, 0.076765424),
        (0.19453496, 0.13924528),
        (0.12662257, 0.10835013),
        (0.056224331, 0.059950072),
    )

    bkg_fractions_2d = (
        (0.61658524, 0),
        (0.25721943, 0.013623876),
        (0.075947519, 0.009301414),
        (0.018317789, 0.0021972499),
        (0.0047723009, 0.0020351877),
    )

    sig_fraction_1d = []
    bkg_fraction_1d = []

    if "tot" in variable:
        for i in range(len(sig_fractions_2d)):
            sig_tot = 0
            bkg_tot = 0
            for j in range(len(sig_fractions_2d[0])):
                sig_tot = sig_tot + sig_fractions_2d[i][j]
                bkg_tot = bkg_tot + bkg_fractions_2d[i][j]

            sig_fraction_1d.append(sig_tot)
            bkg_fraction_1d.append(bkg_tot)

        bin_edges = np.array([1, 2, 3, 4, 5, 6], dtype=float)
    else:
        for i in range(len(sig_fractions_2d[0])):
            sig_tot = 0
            bkg_tot = 0
            for j in range(len(sig_fractions_2d)):
                sig_tot = sig_tot + sig_fractions_2d[j][i]
                bkg_tot = bkg_tot + bkg_fractions_2d[j][i]

            sig_fraction_1d.append(sig_tot)
            bkg_fraction_1d.append(bkg_tot)

        bin_edges = np.array([1, 2, 3], dtype=float)

    print(sig_fraction_1d)
    print(bkg_fraction_1d)
    # sig_fraction =np.array(sig_fraction_1d)/np.sum(np.array(sig_fraction_1d))
    # bkg_fraction =np.array(bkg_fraction_1d)/np.sum(np.array(bkg_fraction_1d))
    sig_fit = n_sig * np.array(sig_fraction_1d)
    bkg_fit = n_bkg * np.array(bkg_fraction_1d)
    total_fit = sig_fit + bkg_fit
    if "tot" in variable:
        bkg_print = bkg_fit[3] + bkg_fit[3]
        sig_print = sig_fit[3] + sig_fit[3]
        print("bkg N fit = ", bkg_print)
        print("sig N fit = ", sig_print)

    # Simple uncertainty propagation from fitted yields only.
    # This does not include MC template statistical/systematic uncertainty.
    sig_unc = n_sig_sigma * np.array(sig_fraction_1d)
    bkg_unc = n_bkg_sigma * np.array(bkg_fraction_1d)
    total_unc = np.sqrt(sig_unc**2 + bkg_unc**2)

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_widths = np.diff(bin_edges)

    # Total fit as a step histogram
    total_step = np.r_[total_fit, total_fit[-1]]

    print(total_fit)
    # print()

    ax.step(
        bin_edges,
        total_step,
        where="post",
        label="Total fit",
        linewidth=2,
    )

    # Uncertainty band as step-like fill
    lower = total_fit - total_unc
    upper = total_fit + total_unc

    ax.fill_between(
        bin_edges,
        np.r_[lower, lower[-1]],
        np.r_[upper, upper[-1]],
        step="post",
        alpha=0.25,
        label="Fit uncertainty",
    )


def get_x_label(variable):
    print(X_LABEL.keys())
    if variable in X_LABEL.keys():
        print("in keys")
        return X_LABEL[variable]

    return variable


def crystal_ball_right_tail(x, mean, sigma, alpha, n):
    x = np.asarray(x, dtype=float)

    if sigma <= 0:
        raise ValueError(f"sigma must be positive, got sigma={sigma}")
    if n <= 0:
        raise ValueError(f"n must be positive, got n={n}")
    if alpha == 0:
        raise ValueError("alpha cannot be zero")

    a = abs(alpha)
    t = (x - mean) / sigma

    A = (n / a) ** n * np.exp(-0.5 * a**2)
    B = n / a - a

    y = np.zeros_like(x, dtype=float)

    # For alpha < 0 in RooCBShape, the tail is on the right.
    core = t < a
    tail = ~core

    y[core] = np.exp(-0.5 * t[core] ** 2)
    y[tail] = A * (B + t[tail]) ** (-n)

    return y


def draw_fit(
    axs,
    parameters,
    mass_min=None,
    mass_max=None,
    n_bins=40,
    bin_edges=None,
    label="Crystal Ball fit",
):
    mean, mean_err = parameters[0], parameters[1]
    sigma, sigma_err = parameters[2], parameters[3]
    alpha, alpha_err = parameters[4], parameters[5]
    n, n_err = parameters[6], parameters[7]
    sumW, sumW_err = parameters[8], parameters[9]

    # Prefer using the same bin edges as the histogram
    if bin_edges is not None:
        mass_min = float(bin_edges[0])
        mass_max = float(bin_edges[-1])
        n_bins = len(bin_edges) - 1

    if mass_min is None:
        mass_min = mean - 5.0 * sigma
    if mass_max is None:
        mass_max = mean + 8.0 * sigma

    mass_min = float(mass_min)
    mass_max = float(mass_max)

    if mass_max <= mass_min:
        raise ValueError(
            f"Bad fit range: mass_min={mass_min}, mass_max={mass_max}. "
            "Pass the same mass_min/mass_max or bin_edges used for the histogram."
        )

    x = np.linspace(mass_min, mass_max, 2000)

    y = crystal_ball_right_tail(
        x,
        mean=mean,
        sigma=sigma,
        alpha=alpha,
        n=n,
    )

    y = np.nan_to_num(y, nan=0.0, posinf=0.0, neginf=0.0)

    area = np.trapz(y, x)

    if area <= 0 or not np.isfinite(area):
        print("DEBUG Crystal Ball")
        print("mass_min, mass_max =", mass_min, mass_max)
        print("mean, sigma, alpha, n =", mean, sigma, alpha, n)
        print("y min/max =", np.min(y), np.max(y))
        print("area =", area)
        raise ValueError(f"Bad Crystal Ball normalization area = {area}")

    bin_width = (mass_max - mass_min) / n_bins
    y_counts = y / area * sumW * bin_width

    axs.plot(
        x,
        y_counts,
        linewidth=2,
        label=(
            rf"{label}: "
            # rf"$m={mean:.2f}\pm{mean_err:.2f}$ GeV, "
            # rf"$\sigma={sigma:.2f}\pm{sigma_err:.2f}$ GeV"
        ),
    )

    return x, y_counts
