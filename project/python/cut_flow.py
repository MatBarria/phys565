import numpy as np
import uproot as ur

path = "../data/proccess_tuples/no_cuts/"

def make_cutflow_table():
    cutflow_variables_mc = [
        "NMuon_valid",
        "N_valid_jets_tot",
        "MET_pt",
        "N_valid_b_jets",
        "weight",
    ]

    cutflow_variables_data = [
        "NMuon_valid",
        "N_valid_jets_tot",
        "MET_pt",
        "N_valid_b_jets",
    ]

    with ur.open(path + "signal.root:tree_output") as file:
        sig = file.arrays(cutflow_variables_mc, library="np")

    with ur.open(path + "background.root:tree_output") as file:
        bkg = file.arrays(cutflow_variables_mc, library="np")

    with ur.open(path + "data_tuples_JES1p00.root:tree_output") as file:
        data = file.arrays(cutflow_variables_data, library="np")

    cut_names = [
        r"All events",
        r"$N_{\mu} > 0$",
        r"$N_{\mathrm{jets}} \geq 4$",
        r"$E_{\mathrm{T}}^{\mathrm{miss}} > 7~\mathrm{GeV}$",
        r"$N_{b\mathrm{-jets}} \geq 1$",
    ]

    def get_masks(branches):
        n_events = len(branches["NMuon_valid"])

        masks = []

        all_events = np.ones(n_events, dtype=bool)
        masks.append(all_events)

        cut1 = all_events & (branches["NMuon_valid"] > 0)
        masks.append(cut1)

        cut2 = cut1 & (branches["N_valid_jets_tot"] >= 4)
        masks.append(cut2)

        cut3 = cut2 & (branches["MET_pt"] > 7)
        masks.append(cut3)

        cut4 = cut3 & (branches["N_valid_b_jets"] >= 1)
        masks.append(cut4)

        return masks

    sig_masks = get_masks(sig)
    bkg_masks = get_masks(bkg)
    data_masks = get_masks(data)

    print("")
    print(r"\begin{table}[htbp]")
    print(r"\centering")
    print(r"\begin{tabular}{lccc}")
    print(r"\hline")
    print(r"Selection & Signal MC & Background MC & Data \\")
    print(r"\hline")

    for cut_name, sig_mask, bkg_mask, data_mask in zip(
        cut_names, sig_masks, bkg_masks, data_masks
    ):
        signal_yield = np.sum(sig["weight"][sig_mask])
        background_yield = np.sum(bkg["weight"][bkg_mask])
        data_yield = np.sum(data_mask)

        print(
            f"{cut_name} & "
            f"{signal_yield:.2f} & "
            f"{background_yield:.2f} & "
            f"{data_yield:.0f} \\\\"
        )

    print(r"\hline")
    print(r"\end{tabular}")
    print(
        r"\caption{Cut-flow table for the event selection. Signal and background yields are taken from weighted Monte Carlo simulation, while data yields are unweighted event counts.}"
    )
    print(r"\label{tab:cutflow}")
    print(r"\end{table}")
    print("")


make_cutflow_table()
