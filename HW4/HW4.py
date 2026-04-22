from math import sqrt
import uproot as ur
import ROOT as root
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.integrate import quad
from scipy.optimize import curve_fit

plt.style.use(hep.style.CMS)
tuple_path = "../project/data/tuples/"
tuples = ["Z-5", "Z", "Z+5"]

file_name = {"Z": "ee2mm_Z.root", "Z-5": "ee2mm_Z-5.root", "Z+5": "ee2mm_Z+5.root"}
labels = {
    "Z": r"$\sqrt{s}= M_z$",
    "Z-5": r"$\sqrt{s}= M_z-5[GeV]$",
    "Z+5": r"$\sqrt{s}= M_z+5[GeV]$",
}
tree_name = ":Ntuple"

variables = ["partOUT_px", "partOUT_py", "partOUT_pz", "partOUT_E", "partOUT_id"]


alpha = 1 / 137.0
mass_Z = 91.1876
gamma_Z = 2.4952
sqrt_s_values = {"Z": mass_Z, "Z-5": mass_Z - 5, "Z+5": mass_Z + 5}

# initial guess
sin2_thetaW = 0.22305
cos2_thetaW = 1.0 - sin2_thetaW
c_a_sm = -0.5
c_v_sm = -0.5 + 2.0 * sin2_thetaW


def get_a_gamma(cos):
    return 1 + cos**2


# def get_c_w(sin2):
# return 1.0 / (4.0 * sin2 * (1 - sin2))


def get_c_w(c_v=c_v_sm):
    return 4.0 / ((2 * c_v + 1) * (3 - 2 * c_v))


def get_a_z(cos, c_a=c_a_sm, c_v=c_v_sm):
    const = get_c_w(c_v)
    return const**2 * (
        (c_v**2 + c_a**2) ** 2 * (1 + cos**2) + 8 * c_a**2 * c_v**2 * cos
    )


def get_c_z(s):
    return s**2 / ((s - mass_Z**2) ** 2 + mass_Z**2 * gamma_Z**2)


def get_a_z_term(cos, s, c_a=c_a_sm, c_v=c_v_sm):
    const_z = get_c_z(s)
    return get_a_z(cos, c_a, c_v) * const_z


def get_a_zgamma(cos, c_a=c_a_sm, c_v=c_v_sm):
    const = get_c_w(c_v)
    return const * (2 * c_v**2 * (1 + cos**2) + 4 * c_a**2 * cos)


def get_c_zgamma(s):
    return (s * (s - mass_Z**2)) / ((s - mass_Z**2) ** 2 + mass_Z**2 * gamma_Z**2)


def get_a_zgamma_term(cos, s, c_a=c_a_sm, c_v=c_v_sm):
    const_zgamma = get_c_zgamma(s)
    return get_a_zgamma(cos, c_a, c_v) * const_zgamma


def get_sigma_dcos(cos, s, c_a=c_a_sm, c_v=c_v_sm):

    if c_v <= -0.5 or c_v >= 1.5:
        return np.ones_like(cos) * 1e10

    return (alpha**2 / (4 * s)) * (
        get_a_gamma(cos)
        + get_a_z_term(cos, s, c_a, c_v)
        + get_a_zgamma_term(cos, s, c_a, c_v)
    )


def get_afb(sqrt_s, c_a, c_v):
    c_w = get_c_w(c_v)
    c_z = get_c_z(sqrt_s**2)
    c_zgamma = get_c_zgamma(sqrt_s**2)

    numerator = c_a**2 * (c_w * 4 * c_z * c_v**2 + 2.0 * c_zgamma)
    denominator = (
        1.0 / c_w + (c_z * c_w) * (c_v**2 + c_a**2) ** 2 + 2.0 * c_zgamma * c_v**2
    )

    return (3.0 / 4.0) * numerator / denominator


def save_figure(fig, outputDirectory, name):

    os.makedirs(outputDirectory, exist_ok=True)
    fig.savefig(outputDirectory + name + ".pdf", bbox_inches="tight")
    fig.savefig(outputDirectory + name + ".png", bbox_inches="tight", dpi=300)
    # fig.savefig(outputDirectory + name + ".pdf")
    # fig.savefig(outputDirectory + name + ".png", dpi=300)
    print(outputDirectory + name + " Has been created")


def get_cos_from_sim(tuple_name, sqrt_s, skip_cut=False):
    with ur.open(tuple_path + file_name[tuple_name] + tree_name) as file:  # type: ignore
        branches = file.arrays(variables, library="np")

    cos_theta = []
    skipped_events = 0

    for idx in range(len(branches["partOUT_E"])):

        if idx % 10000 == 0:
            print("are we moving?", idx)
        if not skip_cut:
            vector_1 = root.Math.PxPyPzEVector(
                branches["partOUT_px"][idx][0],
                branches["partOUT_py"][idx][0],
                branches["partOUT_pz"][idx][0],
                branches["partOUT_E"][idx][0],
            )
            vector_2 = root.Math.PxPyPzEVector(
                branches["partOUT_px"][idx][1],
                branches["partOUT_py"][idx][1],
                branches["partOUT_pz"][idx][1],
                branches["partOUT_E"][idx][1],
            )

            # Di-Muon four momentum vector
            vector_sum = vector_1 + vector_2
            if abs(vector_sum.M() - sqrt_s) > 0.5:
                skipped_events = skipped_events + 1
                continue

        idx_mu = 0
        if branches["partOUT_id"][idx][0] != 13:
            idx_mu = 1

        px = branches["partOUT_px"][idx][idx_mu]
        py = branches["partOUT_py"][idx][idx_mu]
        pz = branches["partOUT_pz"][idx][idx_mu]

        p = (px**2 + py**2 + pz**2) ** 0.5
        if p == 0:
            continue

        cos_theta.append(pz / p)
    print("skipped_events: ", skipped_events)

    return np.array(cos_theta)


def draw_parameters():
    fig, axs = plt.subplots(1, 3, sharey="row")
    # fig, axs = plt.subplots(1, 3)
    width = 16  # 7.056870070568701
    height = 4
    fig.set_size_inches(width, height)
    fig.subplots_adjust(
        left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02
    )
    for idx, tuple in enumerate(tuples):

        s = sqrt_s_values[tuple] ** 2
        cos_theta = np.linspace(-1, 1, 200)
        sigma = get_sigma_dcos(cos_theta, s)
        a_gamma = get_a_gamma(cos_theta)
        a_z = get_a_z_term(cos_theta, s)
        a_zgamma = get_a_zgamma_term(cos_theta, s)
        axs[idx].plot(
            cos_theta, (alpha**2 / (4 * s)) * a_gamma, "blue", label=r"$A_\gamma$"
        )
        axs[idx].plot(cos_theta, (alpha**2 / (4 * s)) * a_z, "g", label=r"$A_Z$")
        axs[idx].plot(
            cos_theta, (alpha**2 / (4 * s)) * a_zgamma, "r", label=r"$A_{z\gamma}$"
        )
        axs[idx].plot(cos_theta, sigma, "black", label=r"$d\sigma/d\Omega$")
        axs[idx].set_xlabel(r"$\cos(\theta)$", fontsize=14)

        axs[idx].set_title(labels[tuple], fontsize=16)

    axs[0].legend(frameon=False, loc="best", fontsize=11)
    axs[0].set_ylabel(r"Distribution", loc="center", fontsize=15)
    axs[0].set_ylim((1.3*np.min((alpha**2 / (4 * s)) * a_zgamma),1.2*np.max(sigma)))

    save_figure(fig, "./", "dSigma-zoom")
    plt.close()


def draw_asymmetry():
    fig, axs = plt.subplots(1, 3, sharey="row")
    width = 16  # 7.056870070568701
    height = 4
    fig.set_size_inches(width, height)
    fig.subplots_adjust(
        left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02
    )

    s = np.linspace(mass_Z - 90, mass_Z + 90, 200)
    c_a = np.linspace(c_a_sm - 0.4, c_a_sm + 0.4, 200)
    c_v = np.linspace(c_v_sm - 0.4, c_v_sm + 0.4, 200)
    asymmetry_s = get_afb(s, c_a_sm, c_v_sm)
    asymmetry_c_a_s = get_afb(mass_Z, c_a, c_v_sm)
    asymmetry_c_a_sp5 = get_afb(mass_Z + 5, c_a, c_v_sm)
    asymmetry_c_a_sm5 = get_afb(mass_Z - 5, c_a, c_v_sm)
    asymmetry_c_v_s = get_afb(mass_Z, c_a_sm, c_v)
    asymmetry_c_v_sp5 = get_afb(mass_Z + 5, c_a_sm, c_v)
    asymmetry_c_v_sm5 = get_afb(mass_Z - 5, c_a_sm, c_v)

    cos_z = get_cos_from_sim("Z", mass_Z)
    cos_zp5 = get_cos_from_sim("Z+5", mass_Z + 5)
    cos_zm5 = get_cos_from_sim("Z-5", mass_Z - 5)
    forward = np.sum(cos_z > 0)
    backward = np.sum(cos_z < 0)
    afb_measured_z = (forward - backward) / (forward + backward)
    forward = np.sum(cos_zp5 > 0)
    backward = np.sum(cos_zp5 < 0)
    afb_measured_zp5 = (forward - backward) / (forward + backward)
    forward = np.sum(cos_zm5 > 0)
    backward = np.sum(cos_zm5 < 0)
    afb_measured_zm5 = (forward - backward) / (forward + backward)

    print(afb_measured_z)
    axs[0].plot(s, asymmetry_s, "b-", linewidth=1.7, label="SM pred")
    axs[0].plot(
        mass_Z,
        afb_measured_z,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=6,
        label="Simulated data",
    )
    axs[0].plot(
        mass_Z + 5,
        afb_measured_zp5,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=6,
    )
    axs[0].plot(
        mass_Z - 5,
        afb_measured_zm5,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=6,
    )
    axs[1].plot(
        c_a, asymmetry_c_a_sm5, "r-", linewidth=1.7, label=r"$\sqrt{s}=M_z-5GeV$"
    )
    axs[1].plot(c_a, asymmetry_c_a_s, "g-", linewidth=1.7, label=r"$\sqrt{s}=M_z$")
    axs[1].plot(
        c_a, asymmetry_c_a_sp5, "b-", linewidth=1.7, label=r"$\sqrt{s}=M_z+5GeV$"
    )
    axs[1].plot(
        c_a_sm,
        afb_measured_z,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        label="Simulated data",
        markersize=6,
    )
    axs[1].plot(
        c_a_sm,
        afb_measured_zp5,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=6,
    )
    axs[1].plot(
        c_a_sm,
        afb_measured_zm5,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=6,
    )
    axs[2].plot(c_v, asymmetry_c_v_sm5, "red", label=r"$\sqrt{s}=M_z-5GeV$")
    axs[2].plot(c_v, asymmetry_c_v_s, "green", label=r"$\sqrt{s}=M_z$")
    axs[2].plot(c_v, asymmetry_c_v_sp5, "blue", label=r"$\sqrt{s}=M_z+5GeV$")
    axs[2].plot(
        c_v_sm,
        afb_measured_z,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=6,
        label="Simulated",
    )
    axs[2].plot(
        c_v_sm,
        afb_measured_zp5,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=6,
    )
    axs[2].plot(
        c_v_sm,
        afb_measured_zm5,
        marker="o",
        linestyle="",
        markerfacecolor="black",
        color="black",
        markersize=6,
    )
    # axs[idx].plot(cos_theta, sigma, "black", label=r"$\d\sigma/d\Omega$")
    axs[0].set_xlabel(r"$\sqrt{s}$", fontsize=14)
    axs[1].set_xlabel(r"$c_a$", fontsize=14)
    axs[2].set_xlabel(r"$c_v$", fontsize=14)

    axs[0].set_title(r"SM $c_A$ and $c_v$", fontsize=16)
    axs[1].set_title(r"$\sqrt{s}=M_z$ SM $c_v$", fontsize=16)
    axs[2].set_title(r"$\sqrt{s}=M_z$ SM $c_a$", fontsize=16)

    axs[0].legend(frameon=False, loc="upper left", fontsize=11)
    axs[1].legend(frameon=False, loc="upper left", fontsize=11)
    axs[2].legend(frameon=False, loc="lower left", fontsize=11)
    axs[0].set_ylabel(r"$A_{FB}$", loc="center", fontsize=15)

    save_figure(fig, "./", "asymmetrys")
    plt.close()


def draw_cos_from_sim():

    fig, axs = plt.subplots(1, 3, sharey="row")
    width = 16  # 7.056870070568701
    height = 4
    fig.set_size_inches(width, height)
    fig.subplots_adjust(
        left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02
    )
    for idx, tuple in enumerate(tuples):

        data = get_cos_from_sim(tuple, sqrt_s_values[tuple])
        events, bins = np.histogram(data, bins=40, range=(-1, 1), density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        # 2. Plot Simulated Data as points with errors
        errors = np.sqrt(len(data) * events * (bins[1] - bins[0])) / (
            len(data) * (bins[1] - bins[0])
        )
        axs[idx].errorbar(
            bin_centers,
            events,
            yerr=errors,
            fmt="ok",
            markersize=4,
            label="Simulated Data",
        )

        cos_theta = np.linspace(-1, 1, 200)
        sigma = get_sigma_dcos(cos_theta, sqrt_s_values[tuple] ** 2)
        sigma_norm = sigma / np.trapz(sigma, cos_theta)
        axs[idx].plot(cos_theta, sigma_norm, "b-", linewidth=1.7, label=r"SM pred")
        axs[idx].set_xlabel(r"$\cos(\theta)$", fontsize=14)
        axs[idx].legend(frameon=False, loc="best", fontsize=11)
        axs[idx].set_title(labels[tuple], fontsize=16)
    
    axs[0].set_ylabel(r"Normalized $d\theta/d\cos(\theta)$", fontsize=14)
    axs[0].set_xlim(bins[0], bins[-1])

    save_figure(fig, "./", "cos_dist")
    plt.close()


def get_normalized_dsigma(cos_data, cv, ca, s_val):

    theo = get_sigma_dcos(cos_data, s_val, cv, ca)
    return theo / np.trapz(theo, cos_data)


def get_fit_model(x_combined, cv, ca):

    n = len(x_combined) // 3
    x1, x2, x3 = x_combined[:n], x_combined[n : 2 * n], x_combined[2 * n :]

    fit1 = get_normalized_dsigma(x1, cv, ca, (mass_Z - 5) ** 2)
    fit2 = get_normalized_dsigma(x2, cv, ca, (mass_Z) ** 2)
    fit3 = get_normalized_dsigma(x3, cv, ca, (mass_Z + 5) ** 2)

    return np.concatenate([fit1, fit2, fit3])


def do_fit():

    data_zm5 = get_cos_from_sim("Z-5", mass_Z - 5)
    data_z = get_cos_from_sim("Z", mass_Z)
    data_zp5 = get_cos_from_sim("Z+5", mass_Z + 5)
    histogram_zm5, bins = np.histogram(data_zm5, bins=50, range=(-1, 1), density=True)
    histogram_z, _ = np.histogram(data_z, bins=50, range=(-1, 1), density=True)
    histogram_zp5, _ = np.histogram(data_zp5, bins=50, range=(-1, 1), density=True)

    bin_centers = (bins[:-1] + bins[1:]) / 2

    x_global = np.concatenate([bin_centers, bin_centers, bin_centers])
    y_global = np.concatenate([histogram_zm5, histogram_z, histogram_zp5])

    # Fit
    popt, pcov = curve_fit(
        get_fit_model, x_global, y_global, p0=[c_v_sm, c_a_sm]
    )

    cv_fit, ca_fit = popt
    cv_err, ca_err = np.sqrt(np.diag(pcov))

    sin2_measured = (cv_fit + 0.5) / 2.0
    sin2_err = cv_err / 2.0

    print("-" * 30)
    print(f"FIT RESULTS:")
    print(f"cV = {cv_fit:.5f} +/- {cv_err:.5f}")
    print(f"cA = {ca_fit:.5f} +/- {ca_err:.5f}")
    print(f"MEASURED sin^2(theta_W) = {sin2_measured:.5f} +/- {sin2_err:.5f}")
    print("-" * 30)

    return ca_fit, cv_fit


def draw_fits(cv_fit, ca_fit):
    fig, axs = plt.subplots(1, 3, sharey="row")
    width = 16  # 7.056870070568701
    height = 4
    fig.set_size_inches(width, height)
    fig.subplots_adjust(
        left=None, bottom=None, right=None, top=None, wspace=0.02, hspace=0.02
    )

    cos_plot = np.linspace(-1, 1, 100)

    for idx, tuple in enumerate(tuples):
        sqrt_s = sqrt_s_values[tuple]
        data = get_cos_from_sim(tuple, sqrt_s)
        events, bins = np.histogram(data, bins=40, range=(-1, 1), density=True)
        bin_centers = (bins[:-1] + bins[1:]) / 2

        ## Scale the error because distribution is normalized
        errors = np.sqrt(len(data) * events * (bins[1] - bins[0])) / (
            len(data) * (bins[1] - bins[0])
        )
        axs[idx].errorbar(
            bin_centers,
            events,
            yerr=errors,
            fmt="ok",
            markersize=4,
            label="Simulated Data",
        )

        fit_curve = get_normalized_dsigma(cos_plot, cv_fit, ca_fit, sqrt_s**2)
        axs[idx].plot(cos_plot, fit_curve, "r-", linewidth=1.7, label="Best Fit")

        axs[idx].set_title(labels[tuple], fontsize=16)
        axs[idx].set_xlabel(r"$\cos(\theta)$", fontsize=14)
    
    axs[0].set_ylabel(r"Normalized $d\theta/d\cos(\theta)$", fontsize=14)
    axs[0].legend(frameon=False, fontsize=11)

    save_figure(fig, "./", "Fits")
    plt.close()


# draw_parameters()
# draw_asymmetry()
# draw_cos_from_sim()
c_a_fit, c_v_fit = do_fit()
# draw_fits(c_a_fit, c_v_fit)
