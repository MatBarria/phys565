import math


MT_GEN = 172.46875  

mu_mc_nominal = 165.358
dmu_mc_nominal = 1.51148

mu_data_nominal = 167.584
dmu_data_nominal = 9.94351


systematics = {
    "JES": {
        "mu_mc_up": 172.00,
        "mu_data_up": 173.00,
        "mu_mc_down": 169.50,
        "mu_data_down": 171.20,
    },
    "b-tag": {
        "mu_mc_up": 165.574,
        "mu_data_up": 168.233,
        "mu_mc_down": 165.574,
        "mu_data_down": 168.233,
    },
    "MET": {
        "mu_mc_up": 165.804,
        "mu_data_up": 166.99,
        "mu_mc_down": 165.804,
        "mu_data_down": 166.99,
    },
}



def extract_mass(mu_data, mu_mc, mt_gen=MT_GEN):
    return mt_gen + (mu_data - mu_mc)


def symmetric_systematic(nominal, up, down):
    shift_up = abs(up - nominal)
    shift_down = abs(down - nominal)
    return max(shift_up, shift_down), shift_up, shift_down



mt_nominal = extract_mass(mu_data_nominal, mu_mc_nominal)

mt_stat = math.sqrt(dmu_data_nominal**2 + dmu_mc_nominal**2)


syst_results = {}

for name, vals in systematics.items():

    mt_up = extract_mass(vals["mu_data_up"], vals["mu_mc_up"])
    mt_down = extract_mass(vals["mu_data_down"], vals["mu_mc_down"])

    syst, shift_up, shift_down = symmetric_systematic(
        mt_nominal,
        mt_up,
        mt_down,
    )

    syst_results[name] = {
        "mt_up": mt_up,
        "mt_down": mt_down,
        "shift_up": shift_up,
        "shift_down": shift_down,
        "syst": syst,
    }


# Total systematic uncertainty
mt_syst = math.sqrt(sum(v["syst"] ** 2 for v in syst_results.values()))

# Total uncertainty
mt_total = math.sqrt(mt_stat**2 + mt_syst**2)

print("                 TOP MASS MEASUREMENT")
print("============================================================")

print(f"MC fitted mean      = {mu_mc_nominal:.3f} ± {dmu_mc_nominal:.3f} GeV")
print(f"Data fitted mean    = {mu_data_nominal:.3f} ± {dmu_data_nominal:.3f} GeV")
print(f"MC generated mass   = {MT_GEN:.3f} GeV")

print("\nCalibrated result:")
print(f"m_t = {mt_nominal:.3f} ± {mt_stat:.3f} (stat.) GeV")

print("\nSystematic uncertainties:")
for name, vals in syst_results.items():
    print(
        f"{name:8s}: "
        f"+ shift = {vals['shift_up']:.3f} GeV, "
        f"- shift = {vals['shift_down']:.3f} GeV, "
        f"assigned = {vals['syst']:.3f} GeV"
    )

print("\nFinal result:")
print(f"m_t = {mt_nominal:.3f} ± {mt_stat:.3f} (stat.) " f"± {mt_syst:.3f} (syst.) GeV")

print(f"m_t = {mt_nominal:.3f} ± {mt_total:.3f} (total) GeV")



print("\n============================================================")
print("                 LATEX TABLE")
print("============================================================\n")

print(r"\begin{table}[htbp]")
print(r"\centering")
print(r"\caption{Summary of the uncertainties in the top-quark mass measurement.}")
print(r"\begin{tabular}{lc}")
print(r"\hline")
print(r"Source & Uncertainty [GeV] \\")
print(r"\hline")
print(f"Statistical & {mt_stat:.2f} \\\\")

for name, vals in syst_results.items():
    print(f"{name} & {vals['syst']:.2f} \\\\")

print(r"\hline")
print(f"Total systematic & {mt_syst:.2f} \\\\")
print(f"Total & {mt_total:.2f} \\\\")
print(r"\hline")
print(r"\end{tabular}")
print(r"\label{tab:top_mass_uncertainties}")
print(r"\end{table}")



print("\n============================================================")
print("                 LATEX FINAL RESULT")
print("============================================================\n")

print(
    r"\["
    + f"\n"
    + f"m_t = {mt_nominal:.2f} \\pm {mt_stat:.2f}\\,\\mathrm{{(stat.)}} "
    + f"\\pm {mt_syst:.2f}\\,\\mathrm{{(syst.)}}~\\mathrm{{GeV}}.\n"
    + r"\]"
)
