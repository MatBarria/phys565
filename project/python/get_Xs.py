import math

import numpy as np
import uproot as ur


def compute_cross_section(
    n_obs,
    n_bkg,
    lumi_pb,
    trigger_efficiency,
    acceptance_efficiency,
    n_trigger_total,
    n_acceptance_total,
    dn_obs_stat=None,
    dn_bkg_stat=0.0,
    dlumi=0.0,
    dn_bkg_syst=0.0,
    cross_section_btag_syst_error=0.0,
):

    if dn_obs_stat is None:
        dn_obs_stat = math.sqrt(max(n_obs, 0.0))

    n_signal = n_obs - n_bkg

    cross_section_pb = n_signal / (lumi_pb * trigger_efficiency * acceptance_efficiency)

    trigger_efficiency_error = math.sqrt(
        trigger_efficiency * (1.0 - trigger_efficiency) / n_trigger_total
    )

    acceptance_efficiency_error = math.sqrt(
        acceptance_efficiency * (1.0 - acceptance_efficiency) / n_acceptance_total
    )

    dn_signal_stat = math.sqrt(dn_obs_stat**2 + dn_bkg_stat**2)

    # Statistical uncertainty on cross section
    cross_section_stat_error = dn_signal_stat / (
        lumi_pb * trigger_efficiency * acceptance_efficiency
    )

    cross_section_bkg_syst_error = dn_bkg_syst / (
        lumi_pb * trigger_efficiency * acceptance_efficiency
    )

    relative_trigger_error = trigger_efficiency_error / trigger_efficiency
    relative_acceptance_error = acceptance_efficiency_error / acceptance_efficiency
    relative_lumi_error = dlumi / lumi_pb

    cross_section_trigger_error = cross_section_pb * relative_trigger_error
    cross_section_acceptance_error = cross_section_pb * relative_acceptance_error
    cross_section_lumi_error = cross_section_pb * relative_lumi_error

    cross_section_syst_error = math.sqrt(
        cross_section_bkg_syst_error**2
        + cross_section_btag_syst_error**2
        + cross_section_trigger_error**2
        + cross_section_acceptance_error**2
        + cross_section_lumi_error**2
    )
    cross_section_total_error = math.sqrt(
        cross_section_stat_error**2 + cross_section_syst_error**2
    )

    return {
        "n_obs": n_obs,
        "n_bkg": n_bkg,
        "n_signal": n_signal,
        "cross_section_pb": cross_section_pb,
        "trigger_efficiency": trigger_efficiency,
        "trigger_efficiency_error": trigger_efficiency_error,
        "acceptance_efficiency": acceptance_efficiency,
        "acceptance_efficiency_error": acceptance_efficiency_error,
        "dn_obs_stat": dn_obs_stat,
        "dn_bkg_stat": dn_bkg_stat,
        "dn_signal_stat": dn_signal_stat,
        "cross_section_stat_error": cross_section_stat_error,
        "dn_bkg_syst": dn_bkg_syst,
        "cross_section_bkg_syst_error": cross_section_bkg_syst_error,
        "cross_section_btag_syst_error": cross_section_btag_syst_error,
        "dlumi": dlumi,
        "cross_section_trigger_error": cross_section_trigger_error,
        "cross_section_acceptance_error": cross_section_acceptance_error,
        "cross_section_lumi_error": cross_section_lumi_error,
        "cross_section_syst_error": cross_section_syst_error,
        "cross_section_total_error": cross_section_total_error,
    }



def print_latex_uncertainty_table(results):
    sigma = results["cross_section_pb"]

    rows = [
        ("Statistical", results["cross_section_stat_error"]),
        ("Background normalization", results["cross_section_bkg_syst_error"]),
        ("b-tagging", results["cross_section_btag_syst_error"]),
        ("Trigger efficiency", results["cross_section_trigger_error"]),
        ("Acceptance", results["cross_section_acceptance_error"]),
        ("Luminosity", results["cross_section_lumi_error"]),
        ("Total systematic", results["cross_section_syst_error"]),
        ("Total", results["cross_section_total_error"]),
    ]

    print("\n==========================================")
    print("LaTeX uncertainty table")
    print("==========================================")
    print(r"\begin{table}[htbp]")
    print(r"\centering")
    print(r"\begin{tabular}{lcc}")
    print(r"\hline")
    print(r"Source & Uncertainty [pb] & Relative uncertainty [\%] \\")
    print(r"\hline")

    for name, err in rows:
        rel = 100.0 * err / sigma if sigma != 0 else 0.0
        print(f"{name} & {err:.2f} & {rel:.1f} \\\\")

    print(r"\hline")
    print(r"\end{tabular}")
    print(r"\caption{Uncertainty contributions to the measured $t\bar{t}$ production cross section.}")
    print(r"\label{tab:cross_section_uncertainties}")
    print(r"\end{table}")
    print("==========================================")

if __name__ == "__main__":

    variables = [
        "weight",
        "NMuon_valid",
        "N_valid_b_jets",
        "N_valid_jets_tot",
        "MET_pt",
        "triggerIsoMu24",
        "mu1_Pt",
    ]
    tuple_path = "../data/proccess_tuples/no_cuts/"
    # tuple_path = "../data/proccess_tuples/no_cuts_btagsys/"
    # tuple_path = "../data/proccess_tuples/"
    background_tree = ur.open(tuple_path + "background.root:tree_output")
    background_branches = background_tree.arrays(variables, library="np")  # type: ignore

    signal_tree = ur.open(tuple_path + "signal.root:tree_output")
    signal_branches = signal_tree.arrays(variables, library="np")  # type: ignore

    data_tree = ur.open(tuple_path + "data_tuples_JES1p00.root:tree_output")
    data_branches = data_tree.arrays(variables, library="np")  # type: ignore

    bool_list_bkg = np.ones_like(background_branches[variables[0]], dtype=bool)
    bool_list_bkg = (bool_list_bkg) & (background_branches["triggerIsoMu24"] == 1)
    bool_list_bkg = (bool_list_bkg) & (background_branches["NMuon_valid"] > 0)
    bool_list_bkg = (bool_list_bkg) & (background_branches["N_valid_jets_tot"] >= 4)
    bool_list_bkg = (bool_list_bkg) & (background_branches["N_valid_b_jets"] >= 1)
    bool_list_bkg = (bool_list_bkg) & (background_branches["MET_pt"] >= 7)

    bool_list_sig = np.ones_like(signal_branches[variables[0]], dtype=bool)
    bool_list_sig = (bool_list_sig) & (signal_branches["NMuon_valid"] > 0)
    bool_list_sig_trigger = (bool_list_sig) & (signal_branches["triggerIsoMu24"] == 1)
    bool_list_sig_Muon = bool_list_sig
    bool_list_sig = (bool_list_sig) & (signal_branches["N_valid_jets_tot"] >= 4)
    bool_list_sig = (bool_list_sig) & (signal_branches["N_valid_b_jets"] >= 1)
    bool_list_sig = (bool_list_sig) & (signal_branches["MET_pt"] >= 7)

    bool_list_data = np.ones_like(data_branches[variables[0]], dtype=bool)
    bool_list_data = (bool_list_data) & (data_branches["triggerIsoMu24"] == 1)
    bool_list_data = (bool_list_data) & (data_branches["NMuon_valid"] > 0)
    bool_list_data = (bool_list_data) & (data_branches["N_valid_jets_tot"] >= 4)
    bool_list_data = (bool_list_data) & (data_branches["N_valid_b_jets"] >= 1)
    bool_list_data = (bool_list_data) & (data_branches["MET_pt"] >= 7)

    background_branches["weight"] = background_branches["weight"] * 0.9
    # signal_branches["weight"] = signal_branches["weight"] * 0.9

    n_bkg = np.sum(background_branches["weight"][bool_list_bkg])
    n_obs = len(data_branches["mu1_Pt"][bool_list_data])
    dn_bkg_stat = math.sqrt(np.sum(background_branches["weight"][bool_list_bkg] ** 2))

    lumi_pb = 50.0
    # 525.9 +/- 26.5843

    # dn_obs_stat = 16.7797
    # n_bkg = 0
    # n_obs = 133.048

    # dn_obs_stat = 17.6393 * ((59.0248 + 98.7158 + 37.486 + 34.20564) / 506.941)
    # n_bkg = 0
    # n_obs = (59.0248 + 98.7158 + 37.486 + 34.20564)
    # dn_bkg_stat = 0

    n_trigger_total = np.sum(signal_branches["weight"][bool_list_sig_Muon])
    n_trigger = np.sum(signal_branches["weight"][bool_list_sig_trigger])
    n_acceptance_total = np.sum(signal_branches["weight"])
    n_acceptance = np.sum(signal_branches["weight"][bool_list_sig])
    print(n_trigger_total)
    print(n_acceptance_total)
    # Denominators used to calculate efficiencies
    # n_trigger_total = 7928.8203
    # n_acceptance_total = 1001.7932

    # trigger_efficiency = 0.12634833  # 1001.7932/7928.8203
    # acceptance_efficiency = 0.056316022  # 56.417007/1001.7932
    trigger_efficiency = n_trigger / n_trigger_total  # 1001.7932/7928.8203
    acceptance_efficiency = n_acceptance / n_acceptance_total  # 56.417007/1001.7932
    print(trigger_efficiency)
    print(acceptance_efficiency)

    dn_bkg_syst = 0.20 * n_bkg
    dlumi = 50 * 0.03
    cross_section_btag_syst_error = 160.614012 - 158.977347
    results = compute_cross_section(
        n_obs=n_obs,
        n_bkg=n_bkg,
        lumi_pb=lumi_pb,
        trigger_efficiency=trigger_efficiency,
        acceptance_efficiency=acceptance_efficiency,
        n_trigger_total=n_trigger_total,
        n_acceptance_total=n_acceptance_total,
        # dn_obs_stat=dn_obs_stat,  # uses sqrt(n_obs)
        dn_bkg_stat=dn_bkg_stat,
        dlumi=dlumi,
        dn_bkg_syst=dn_bkg_syst,
        cross_section_btag_syst_error=cross_section_btag_syst_error,
    )

    print("==========================================")
    print("Cross section result")
    print("==========================================")
    print(f"N_obs                         = {results['n_obs']:.6f}")
    print(f"N_bkg                         = {results['n_bkg']:.6f}")
    print(f"N_signal                      = {results['n_signal']:.6f}")
    print(f"cross section                 = {results['cross_section_pb']:.6f} pb")
    print("------------------------------------------")
    print(
        f"trigger efficiency            = "
        f"{results['trigger_efficiency']:.6f} ± {results['trigger_efficiency_error']:.6f}"
    )
    print(
        f"acceptance efficiency         = "
        f"{results['acceptance_efficiency']:.6f} ± {results['acceptance_efficiency_error']:.6f}"
    )
    print("------------------------------------------")
    print(
        f"statistical uncertainty       = {results['cross_section_stat_error']:.6f} pb"
    )
    print(
        f"systematic uncertainty        = {results['cross_section_syst_error']:.6f} pb"
    )
    print(
        f"total uncertainty             = {results['cross_section_total_error']:.6f} pb"
    )
    print("------------------------------------------")
    print(
        f"background syst contribution  = {results['cross_section_bkg_syst_error']:.6f} pb"
    )
    print(
        f"trigger contribution          = {results['cross_section_trigger_error']:.6f} pb"
    )
    print(
        f"acceptance contribution       = {results['cross_section_acceptance_error']:.6f} pb"
    )
    print(
        f"luminosity contribution       = {results['cross_section_lumi_error']:.6f} pb"
    )
    print(
        f"b-tag syst contribution       = {results['cross_section_btag_syst_error']:.6f} pb"
    )
    print("==========================================")


    print_latex_uncertainty_table(results)
