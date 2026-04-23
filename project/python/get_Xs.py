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
):
    
    if dn_obs_stat is None:
        dn_obs_stat = math.sqrt(max(n_obs, 0.0))

    # ----------------------------------------
    # Signal yield
    # ----------------------------------------
    n_signal = n_obs - n_bkg
    if n_signal <= 0:
        raise ValueError(f"n_signal = n_obs - n_bkg = {n_signal:.6f} <= 0")

    # ----------------------------------------
    # Cross section
    # ----------------------------------------
    cross_section_pb = n_signal / (
        lumi_pb * trigger_efficiency * acceptance_efficiency
    )

    # ----------------------------------------
    # Efficiency errors from binomial formula
    # ----------------------------------------
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
    relative_lumi_error = dlumi / lumi_pb if lumi_pb > 0 else 0.0

    cross_section_trigger_error = cross_section_pb * relative_trigger_error
    cross_section_acceptance_error = cross_section_pb * relative_acceptance_error
    cross_section_lumi_error = cross_section_pb * relative_lumi_error

    cross_section_syst_error = math.sqrt(
        cross_section_bkg_syst_error**2
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
        "dlumi": dlumi,
        "cross_section_trigger_error": cross_section_trigger_error,
        "cross_section_acceptance_error": cross_section_acceptance_error,
        "cross_section_lumi_error": cross_section_lumi_error,
        "cross_section_syst_error": cross_section_syst_error,
        "cross_section_total_error": cross_section_total_error,
    }


if __name__ == "__main__":

    tuple_path = "../data/proccess_tuples/"
    background_tree = ur.open(tuple_path + "background.root:tree_output")
    background_branches = background_tree.arrays(["weight"], library="np")  # type: ignore
    background_tree = ur.open(tuple_path + "background.root:tree_output")

    data_tree = ur.open(tuple_path + "data_tuples.root:tree_output")
    data_branches = data_tree.arrays(["mu1_Pt"], library="np")  # type: ignore
    n_bkg = np.sum(background_branches["weight"])
    n_obs = len(data_branches["mu1_Pt"])
    dn_bkg_stat = math.sqrt(np.sum(background_branches["weight"]))
    # n_obs = 120
    # n_bkg = 35.0
    lumi_pb = 50.0


    # Denominators used to calculate efficiencies
    n_trigger_total = 7928.8203
    n_acceptance_total = 1001.7932
    
    trigger_efficiency = 0.12634833 # 1001.7932/7928.8203
    acceptance_efficiency =  0.056316022 #56.417007/1001.7932
    
    # dn_bkg_syst = 3.5
    dlumi = 50*0.025

    results = compute_cross_section(
        n_obs=n_obs,
        n_bkg=n_bkg,
        lumi_pb=lumi_pb,
        trigger_efficiency=trigger_efficiency,
        acceptance_efficiency=acceptance_efficiency,
        n_trigger_total=n_trigger_total,
        n_acceptance_total=n_acceptance_total,
        dn_obs_stat=None,   # uses sqrt(n_obs)
        # dn_bkg_stat=dn_bkg_stat,
        dlumi=dlumi,
        # dn_bkg_syst=dn_bkg_syst,
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
    print(f"statistical uncertainty       = {results['cross_section_stat_error']:.6f} pb")
    print(f"systematic uncertainty        = {results['cross_section_syst_error']:.6f} pb")
    print(f"total uncertainty             = {results['cross_section_total_error']:.6f} pb")
    print("------------------------------------------")
    print(f"background syst contribution  = {results['cross_section_bkg_syst_error']:.6f} pb")
    print(f"trigger contribution          = {results['cross_section_trigger_error']:.6f} pb")
    print(f"acceptance contribution       = {results['cross_section_acceptance_error']:.6f} pb")
    print(f"luminosity contribution       = {results['cross_section_lumi_error']:.6f} pb")
    print("==========================================")

    


