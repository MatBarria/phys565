#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <Math/Vector4D.h>
#include <ROOT/RVec.hxx>
#include <array>
#include <Math/Vector3D.h>

constexpr int N_MEAS = 17;
constexpr int N_PAR = 15;

constexpr double BOTTOM_MASS = 4.183;
constexpr double TOP_MASS = 173.0;
constexpr double W_MASS = 80.36;
constexpr double NU_MASS = 0.0;
constexpr double E_MASS = 0.000511;
constexpr double MU_MASS = 0.10566;
constexpr double TAU_MASS = 1.77693;
constexpr double Q_MASS = 0.0;

struct FitResult {
    bool valid = false;
    double minValue = 1e99;
    std::array<double, N_PAR> pfit{};
};
enum MeasParIndex {
    X_MU_PX = 0,
    X_MU_PY = 1,
    X_MU_PZ = 2,
    
    X_J1_E = 3,
    X_J1_ETA = 4,
    X_J1_PHI = 5,

    X_J2_E = 6,
    X_J2_ETA = 7,
    X_J2_PHI = 8,

    X_J3_E = 9,
    X_J3_ETA = 10,
    X_J3_PHI = 11,

    X_J4_E = 12,
    X_J4_ETA = 13,
    X_J4_PHI = 14,

    X_MET_X = 15,
    X_MET_Y = 16,

    // NU_PZ = 17
};

enum FitParIndex {
    MU_PX = 0,
    MU_PY = 1,
    MU_PZ = 2,

    J1_E = 3,
    J1_ETA = 4,
    J1_PHI = 5,

    J2_ETA = 6,
    J2_PHI = 7,

    J3_E = 8,
    J3_ETA = 9,
    J3_PHI = 10,

    //J4_E = 11,
    J4_ETA = 11,
    J4_PHI = 12,

    MET_X = 13,
    MET_Y = 14,
    // J2_E = 6,
    // J2_ETA = 7,
    // J2_PHI = 8,

    // J3_E = 9,
    // J3_ETA = 10,
    // J3_PHI = 11,

    // J4_E = 12,
    // J4_ETA = 13,
    // J4_PHI = 14,

    // MET_X = 15,
    // MET_Y = 16,

    // NU_PZ = 17
};

ROOT::Math::PxPyPzEVector muonFromPxPyPz(double px, double py, double pz);
ROOT::Math::PxPyPzEVector jetFromEetaPhi(double E, double eta, double phi,
                                         double mass);
ROOT::Math::PxPyPzEVector neutrinoFromMet(double metx, double mety,
                                          double nupz);
double CalculateNeutrinoPz(double mu_px, double mu_py, double mu_pz,
                           double nu_px, double nu_py);

FitResult minimizeEvent(const std::array<double, N_PAR> &p_start,
                        const std::array<double, N_MEAS> &x_meas,
                        const std::array<double, N_MEAS> &sigma);

double deltaPhi(double phi1, double phi2);
bool isPhiIndex(int i);

double chi2Diagonal(const double *p, const std::array<double, N_MEAS> &x_meas,
                    const std::array<double, N_MEAS> &sigma);

double totalFitFunction(const double *par,
                        const std::array<double, N_MEAS> &x_meas,
                        const std::array<double, N_MEAS> &sigma);

double jetSigmaE(double E, double eta);
double jetSigmaPhi();
double jetSigmaEta();
double metSigma();
double muSigma(double p);
double solveJ2E(double j1E, double j1eta, double j1phi, double j2eta,
                double j2phi);

ROOT::Math::XYZVector unitVectorFromEtaPhi(double eta, double phi);
double solveJ4EFromTopEquality(const ROOT::Math::PxPyPzEVector &Wlep,
                               double targetTopMass, double j4eta, double j4phi,
                               double measuredJ4E = -1.0, double bMass = 5.0);
#endif
