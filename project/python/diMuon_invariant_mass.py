import uproot as ur
import ROOT as root
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt

from utils.helper import (
    get_canvas,
    save_figure,
    get_background_label_list,
    get_color_list,
)
from utils.constants import (
    BACKGROUND_SOURCES,
    SIGNAL_SOURCES,
    SOURCES_LABEL,
)

plt.style.use(hep.style.CMS)
variables = [
    "NMuon",
    "Muon_Px",
    "Muon_Py",
    "Muon_Pz",
    "Muon_E",
    "Muon_Charge",
    "EventWeight",
]


tuple_path = "../data/tuples/"

fig, ax = get_canvas()
sources = BACKGROUND_SOURCES + SIGNAL_SOURCES + ["data"]
colors = get_color_list(len(sources))[::-1]
bkg_histograms = []
signal_histograms = []
data_histogram = 0
bins = [25,125]

for source_idx, source in enumerate(sources):
    file_name = source + ".root:events"

    print("***Working on " + source + "***")
    with ur.open(tuple_path + file_name) as file:
        branches = file.arrays(variables, library="np")

    invariant_mass = []
    weights = []

    for idx, n_muon in enumerate(branches["NMuon"]):

        if n_muon < 2:
            continue
        pair_idx = -1
        if branches["Muon_Charge"][idx][0] * branches["Muon_Charge"][idx][1] == -1:
            pair_idx = 1
        elif (
            n_muon > 2
            and branches["Muon_Charge"][idx][0] * branches["Muon_Charge"][idx][2] == -1
        ):
            pair_idx = 2
        else:
            continue
        vector_1 = root.Math.PxPyPzEVector(
            branches["Muon_Px"][idx][0],
            branches["Muon_Py"][idx][0],
            branches["Muon_Pz"][idx][0],
            branches["Muon_E"][idx][0],
        )
        vector_2 = root.Math.PxPyPzEVector(
            branches["Muon_Px"][idx][pair_idx],
            branches["Muon_Py"][idx][pair_idx],
            branches["Muon_Pz"][idx][pair_idx],
            branches["Muon_E"][idx][pair_idx],
        )

        vector_sum = vector_1 + vector_2
        invariant_mass.append(vector_sum.M())
        weights.append(branches["EventWeight"][idx])

    if len(invariant_mass) == 0:
        print("   No diMuon events here")
        # source_idx = source_idx -1
        BACKGROUND_SOURCES.remove(source)
        continue
    invariant_mass = np.array(invariant_mass)
    weights = np.array(weights)

    histogram, bins = np.histogram(
        invariant_mass,
        bins=100,
        range=(25, 125),
        weights=weights,
    )
    # histogram = histogram / np.sum(histogram)
    if source in BACKGROUND_SOURCES: 
        bkg_histograms.append(histogram)
    elif source in SIGNAL_SOURCES:
        signal_histograms.append(histogram)
    else:
        data_histogram = histogram

hep.histplot(
    bkg_histograms,
    bins=bins,
    yerr=True,
    histtype="fill",
    label=get_background_label_list(BACKGROUND_SOURCES),
    ax=ax,
    stack=True,
    color=get_color_list(len(BACKGROUND_SOURCES)),
)

hep.histplot(
    data_histogram,
    bins,
    yerr=True,
    histtype="errorbar",
    label=SOURCES_LABEL[source],
    color="black",
    ax=ax,
)

for source, histogram in zip(SIGNAL_SOURCES, signal_histograms):
    hep.histplot(
        histogram,
        bins,
        yerr=False,
        label=SOURCES_LABEL[source],
        color="red",
        ax=ax,
    )

hep.cms.label(
    data="True",
    label="",
    year=2011,
    com="7",
    lumi="50",
    ax=ax,
)
ax.set_ylabel(r"Events")
ax.set_ylim(0.1, 1000 * np.max(data_histogram))
ax.set_xlim(bins[0], bins[-1])
ax.set_yscale("log")
ax.legend(frameon=False, loc="upper right", ncols=2)
ax.tick_params(axis="x", which="both", bottom=True, top=True, labelbottom=False)


output_directory = "../plots/diMuon_invariant_mass/"
output_name = "diMuon_invariant_mass"

# output_directory = get_output_directory(variable, output_directory, variables_type)

save_figure(fig, output_directory, output_name)

plt.close()

