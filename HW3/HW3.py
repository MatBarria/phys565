import uproot as ur
import ROOT as root
import mplhep as hep
import numpy as np
import matplotlib.pyplot as plt
import os

plt.style.use(hep.style.CMS)
tuple_path = "../project/data/tuples/"
file_name = "ee2mm.root"
tree_name = ":Ntuple"

variables =  ["partOUT_px", "partOUT_py","partOUT_pz","partOUT_E",]



def save_figure(fig, outputDirectory, name):

    os.makedirs(outputDirectory, exist_ok=True)
    fig.savefig(outputDirectory + name + ".pdf", bbox_inches="tight")
    fig.savefig(outputDirectory + name + ".png", bbox_inches="tight", dpi=300)
    fig.savefig(outputDirectory + name + ".pdf")
    fig.savefig(outputDirectory + name + ".png", dpi=300)
    print(outputDirectory + name + " Has been created")


with ur.open(tuple_path + file_name+ tree_name) as file: # type ignore
    branches = file.arrays(variables, library="np")

invariant_mass = []
between_angle = []
azimuthal_angle_1 = []
azimuthal_angle_2 = []
weights = []

# Loop over all the events
for idx in range(len(branches["partOUT_E"])):

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
    between_angle.append(root.Math.VectorUtil.Angle(vector_1.Vect(), vector_2.Vect()))
    azimuthal_angle_1.append(vector_1.Phi())
    azimuthal_angle_2.append(vector_2.Phi())
    invariant_mass.append(vector_sum.M())
    

invariant_mass = np.array(invariant_mass)
between_angle = np.array(between_angle)
azimuthal_angle_1 = np.array(azimuthal_angle_1)
azimuthal_angle_2 = np.array(azimuthal_angle_2)

# create the 4 different histograms 
histogram_mass, bins_mass = np.histogram(
    invariant_mass,
    bins=100,
    range=(85, 95),
)


histogram_between_angle, bins_between_angle = np.histogram(
    between_angle,
    bins=100,
    range=(-3.5, 3.5),
)


histogram_azimuthal_1, bins_azimuthal_1 = np.histogram(
    azimuthal_angle_1,
    bins=100,
    range=(-3.5, 3.5),
)


histogram_azimuthal_2, bins_azimuthal_2 = np.histogram(
    azimuthal_angle_2,
    bins=100,
    range=(-3.5, 3.5),
)

# Drow di-muon mass
fig, ax = plt.subplots(1, 1)
hep.histplot(
    histogram_mass,
    bins_mass,
    yerr=False,
    histtype="fill",
    label="Di-Muon mass",
    ax=ax,
)
ax.set_ylabel(r"Events")
ax.set_xlabel(r"Di-Muon_mass")
ax.set_xlim(bins_mass[0], bins_mass[-1])
ax.legend(frameon=False, loc="upper right", ncols=2)

save_figure(fig, "./", "Di-Muon_mass")
plt.close()


# Drow angle between the muons
fig, ax = plt.subplots(1, 1)
hep.histplot(
    histogram_between_angle,
    bins_between_angle,
    yerr=False,
    histtype="fill",
    label= r"Angle $\mu_1$ and $\mu_2$ ",
    # color="Blue",
    ax=ax,
)
ax.set_ylabel(r"Events")
ax.set_xlabel(r"Angle[rad]")
ax.set_xlim(bins_between_angle[0], bins_between_angle[-1])
ax.legend(frameon=False, loc="upper right", ncols=2)

save_figure(fig, "./", "Between_angles")
plt.close()

# Draw azimuthal angle for both muons
fig, ax = plt.subplots(1, 1)
hep.histplot(
    histogram_azimuthal_1,
    bins_azimuthal_1,
    yerr=False,
    histtype="fill",
    label= r"Azimuthal Angle $\mu_1$ ",
    color="blue",
    alpha=0.4,
    ax=ax,
)
hep.histplot(
    histogram_azimuthal_2,
    bins_azimuthal_2,
    yerr=False,
    histtype="fill",
    label= r"Azimuthal Angle $\mu_2$ ",
    color="red",
    alpha=0.4,
    ax=ax,
)
ax.set_ylabel(r"Events")
ax.set_xlabel(r"Azimuthal Angle[rad]")
ax.set_ylim(0.,  1.3*np.max(histogram_azimuthal_2))
ax.set_xlim(bins_azimuthal_1[0], bins_azimuthal_1[-1])
ax.legend(frameon=False, loc="upper right", ncols=2)

save_figure(fig, "./", "azimuthal_angle")
plt.close()

