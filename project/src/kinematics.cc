#include "../lib/kinematics.h"

#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Vector3D.h>
#include <cmath>
#include <cstdio>
FitResult minimizeEvent(const std::array<double, N_PAR> &p_start,
                        const std::array<double, N_MEAS> &x_meas,
                        const std::array<double, N_MEAS> &sigma,
                        const float nu_pz) {
    FitResult result;

    ROOT::Math::Minimizer *min =
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    if (!min) {
        std::cerr << "ERROR: could not create Minuit2 minimizer\n";
        return result;
    }

    min->SetMaxFunctionCalls(100000);
    min->SetMaxIterations(10000);
    min->SetTolerance(1e-3);
    min->SetPrintLevel(0);

    ROOT::Math::Functor f(
        [&](const double *par) { return totalFitFunction(par, x_meas, sigma); },
        N_PAR);

    min->SetFunction(f);

    // Set starting values and step sizes
    min->SetVariable(MU_PX, "mu_px", p_start[MU_PX],
                     0.01 * std::abs(p_start[MU_PX]) + 1e-3);
    min->SetVariable(MU_PY, "mu_py", p_start[MU_PY],
                     0.01 * std::abs(p_start[MU_PY]) + 1e-3);
    min->SetVariable(MU_PZ, "mu_pz", p_start[MU_PZ],
                     0.01 * std::abs(p_start[MU_PZ]) + 1e-3);

    min->SetVariable(J1_E, "j1_E", p_start[J1_E], 0.05 * p_start[J1_E]);
    min->SetVariable(J1_ETA, "j1_eta", p_start[J1_ETA], 0.01);
    min->SetVariable(J1_PHI, "j1_phi", p_start[J1_PHI], 0.01);

    // min->SetVariable(J2_E, "j2_E", p_start[J2_E], 0.05 * p_start[J2_E]);
    min->SetVariable(J2_ETA, "j2_eta", p_start[J2_ETA], 0.01);
    min->SetVariable(J2_PHI, "j2_phi", p_start[J2_PHI], 0.01);

    min->SetVariable(J3_E, "j3_E", p_start[J3_E], 0.05 * p_start[J3_E]);
    min->SetVariable(J3_ETA, "j3_eta", p_start[J3_ETA], 0.01);
    min->SetVariable(J3_PHI, "j3_phi", p_start[J3_PHI], 0.01);

    //min->SetVariable(J4_E, "j4_E", p_start[J4_E], 0.05 * p_start[J4_E]);
    min->SetVariable(J4_ETA, "j4_eta", p_start[J4_ETA], 0.01);
    min->SetVariable(J4_PHI, "j4_phi", p_start[J4_PHI], 0.01);

    min->SetVariable(MET_X, "metx", p_start[MET_X], 5.0);
    min->SetVariable(MET_Y, "mety", p_start[MET_Y], 5.0);

    // min->SetVariable(NU_PZ, "nu_pz", p_start[NU_PZ], 10.0);

    bool ok = min->Minimize();

    if (!ok) {
        delete min;
        return result;
    }

    const double *xs = min->X();

    result.valid = true;
    result.minValue = min->MinValue();

    for (int i = 0; i < N_PAR; ++i) {
        result.pfit[i] = xs[i];
    }

    delete min;
    return result;
}

ROOT::Math::PxPyPzEVector muonFromPxPyPz(double px, double py, double pz) {
    const double mMu = 0.105658; // GeV

    double E = std::sqrt(px * px + py * py + pz * pz + mMu * mMu);

    return ROOT::Math::PxPyPzEVector(px, py, pz, E);
};

double deltaPhi(double phi1, double phi2) {
    double dphi = phi1 - phi2;

    while (dphi > M_PI) {
        dphi -= 2.0 * M_PI;
    }

    while (dphi <= -M_PI) {
        dphi += 2.0 * M_PI;
    }

    return dphi;
};

bool isPhiIndex(int i) {

    return i == J1_PHI || i == J2_PHI || i == J3_PHI || i == J4_PHI;
};

bool isPhiMeasIndex(int i) {

    return i == X_J1_PHI || i == X_J2_PHI || i == X_J3_PHI || i == X_J4_PHI;
};
// double chi2Diagonal(const double *p, const std::array<double, N_MEAS>
// &x_meas, const std::array<double, N_MEAS> &sigma) {
// double chi2 = 0.0;

// int j = 0;
// for (int i = 0; i < N_PAR; ++i) {
// if (sigma[i] <= 0.0) {
// continue;
//}
// double diff;

// if (i == X_J2_E) {
// j++;
//}

// if (isPhiIndex(i)) {
// diff = deltaPhi(p[i], x_meas[j]);
//} else {
// diff = p[i] - x_meas[j];
//}

// double pull = diff / sigma[j];
// j++;
// chi2 += pull * pull;
//}

// return chi2;
//};
double chi2Diagonal(const double *p, const std::array<double, N_MEAS> &x_meas,
                    const std::array<double, N_MEAS> &sigma) {
    double chi2 = 0.0;

    auto addPull = [&](double fitValue, int imeas) {
        if (sigma[imeas] <= 0.0)
            return;

        double diff;

        if (isPhiMeasIndex(imeas)) {
            diff = deltaPhi(fitValue, x_meas[imeas]);
        } else {
            diff = fitValue - x_meas[imeas];
        }

        double pull = diff / sigma[imeas];
        chi2 += pull * pull;
    };

    // Normal parameters
    addPull(p[MU_PX], X_MU_PX);
    addPull(p[MU_PY], X_MU_PY);
    addPull(p[MU_PZ], X_MU_PZ);

    addPull(p[J1_E], X_J1_E);
    addPull(p[J1_ETA], X_J1_ETA);
    addPull(p[J1_PHI], X_J1_PHI);

    // J2_E is not in p. It is solved.
    double j2E = solveJ2E(p[J1_E], p[J1_ETA], p[J1_PHI], p[J2_ETA], p[J2_PHI]);

    auto j1 = jetFromEetaPhi(p[J1_E], p[J1_ETA], p[J1_PHI], 0.);
    auto j2 = jetFromEetaPhi(j2E, p[J2_ETA], p[J2_PHI], 0.);
    auto j3 = jetFromEetaPhi(p[J3_E], p[J3_ETA], p[J3_PHI], BOTTOM_MASS);

    auto mu = muonFromPxPyPz(p[MU_PX], p[MU_PY], p[MU_PZ]);
    double nu_pz =
        CalculateNeutrinoPz(p[MU_PX], p[MU_PY], p[MU_PZ], p[MET_X], p[MET_Y]);
    auto nu = neutrinoFromMet(p[MET_X], p[MET_Y], nu_pz);

    auto W_had = j1 + j2;
    auto W_lep = mu + nu;
    auto top_had = W_had + j3;

    double j4E = solveJ4EFromTopEquality(W_lep, top_had.M(), p[J4_ETA],
                                         p[J4_PHI], x_meas[X_J4_E], BOTTOM_MASS);

    // auto j4 = jetFromEetaPhi(p[J4_E], p[J4_ETA], p[J4_PHI], W_MASS);
    //if (j2E < 0.0)
        //return 1e99;

    addPull(j2E, X_J2_E);
    addPull(p[J2_ETA], X_J2_ETA);
    addPull(p[J2_PHI], X_J2_PHI);

    addPull(p[J3_E], X_J3_E);
    addPull(p[J3_ETA], X_J3_ETA);
    addPull(p[J3_PHI], X_J3_PHI);

    addPull(j4E, X_J4_E);
    //addPull(p[J4_E], X_J4_E);
    addPull(p[J4_ETA], X_J4_ETA);
    addPull(p[J4_PHI], X_J4_PHI);

    addPull(p[MET_X], X_MET_X);
    addPull(p[MET_Y], X_MET_Y);

    return chi2;
}
// FitResult minimizeEvent(const std::array<double, N_PAR> &p_start,
// const std::array<double, N_MEAS> &x_meas,
// const std::array<double, N_MEAS> &sigma);

ROOT::Math::PxPyPzEVector jetFromEetaPhi(double E, double eta, double phi,
                                         double mass = 0.0) {
    double p2 = E * E - mass * mass;

    if (p2 < 0.0) {
        p2 = 0.0;
    }

    double p = std::sqrt(p2);
    double pt = p / std::cosh(eta);

    double px = pt * std::cos(phi);
    double py = pt * std::sin(phi);
    double pz = pt * std::sinh(eta);

    return ROOT::Math::PxPyPzEVector(px, py, pz, E);
}

ROOT::Math::PxPyPzEVector neutrinoFromMet(double metx, double mety,
                                          double nupz) {
    double Enu = std::sqrt(metx * metx + mety * mety + nupz * nupz);

    return ROOT::Math::PxPyPzEVector(metx, mety, nupz, Enu);
};

double totalFitFunction(const double *par, const std::array<double, 17> &x_meas,
                        const std::array<double, 17> &sigma) {
    // 1. detector chi2
    double chi2 = chi2Diagonal(par, x_meas, sigma);

    // 2. build fitted particles
    // auto mu = muonFromPxPyPz(par[MU_PX], par[MU_PY], par[MU_PZ]);
    // auto nu = neutrinoFromMet(par[MET_X], par[MET_Y], nu_pz);

    // double jet2_e =
    // solveJ2E(par[J1_E], par[J1_ETA], par[J1_PHI], par[J2_ETA], par[J1_PHI]);
    // auto j1 = jetFromEetaPhi(par[J1_E], par[J1_ETA], par[J1_PHI], 0.);
    // auto j2 = jetFromEetaPhi(jet2_e, par[J2_ETA], par[J2_PHI], 0.);
    // auto j3 = jetFromEetaPhi(par[J3_E], par[J3_ETA], par[J3_PHI], W_MASS);
    // double jet4_e =
    // sol(par[J1_E], par[J1_ETA], par[J1_PHI], par[J2_ETA], par[J1_PHI]);
    // auto j4 = jetFromEetaPhi(par[J4_E], par[J4_ETA], par[J4_PHI], W_MASS);

    // double nu_pz = CalculateNeutrinoPz(par[MU_PX], par[MU_PY], par[MU_PZ],
    // par[MET_X], par[MET_Y]);
    // auto nu = neutrinoFromMet(par[MET_X], par[MET_Y], nu_pz);

    //// 3. build W and top candidates
    // auto W_had = j1 + j2;
    // auto W_lep = mu + nu;

    // auto top_had = W_had + j3;
    // auto top_lep = W_lep + j4;

    // 4. constraints
    // double C1 = W_had.M() - W_MASS;
    // double C2 = W_lep.M() - W_MASS;
    // double C3 = top_had.M() - top_lep.M();

    // std::cout << "c2: " << C2 << std::endl;
    //  5. constraint penalty
    // double sigmaW = 1;   // start loose
    // double sigmaTop = 1; // start loose
    // chi2 += std::pow((jet2_e - x_meas[X_J2_E]), 2) / sigma[X_J2_E];
    // double constraintChi2 = std::pow(C3 / sigmaTop, 2);
    // std::pow(C1 / sigmaW, 2) +
    //  std::pow(C2 / sigmaW, 2) +

    return chi2;
    // return chi2 + constraintChi2;
};

double jetSigmaE(double E, double eta) {
    double a = 0.0;
    double b = 0.0;

    double abseta = std::abs(eta);

    if (abseta < 0.8) {
        a = 0.036;
        b = 1.145;
    } else if (abseta < 1.4) {
        a = 0.082;
        b = 1.264;
    } else if (abseta < 2.0) {
        a = 0.046;
        b = 1.305;
    } else {
        a = 0.046;
        b = 1.305;
    }

    double rel = std::sqrt(a * a + (b / std::sqrt(E)) * (b / std::sqrt(E)));
    return rel * 2 * E;
};

double jetSigmaPhi() {
    // if (std::abs(eta) < 0.8)
    // return 0.04;
    return 0.1;
};

double jetSigmaEta() {
    // if (std::abs(eta) < 0.8)
    // return 0.04;
    return 0.1;
};
double metSigma() {
    return 12.0; // GeV
};

double muSigma(const double p) {
    // return 0.1 * p; // GeV
    return 10; // GeV
};

double CalculateNeutrinoPz(double mu_px, double mu_py, double mu_pz,
                           double nu_px, double nu_py) {

    double mu_E = std::sqrt(mu_px * mu_px + mu_py * mu_py + mu_pz * mu_pz);

    double Lambda = (W_MASS * W_MASS - MU_MASS * MU_MASS) / 2.0 +
                    (mu_px * nu_px) + (mu_py * nu_py);

    // 2. Calculate the quadratic coefficients A, B, and C
    double A = (mu_E * mu_E) - (mu_pz * mu_pz);
    double B = -2.0 * Lambda * mu_pz;

    double pt_nu_sq = (nu_px * nu_px) + (nu_py * nu_py);
    double C = (mu_E * mu_E * pt_nu_sq) - (Lambda * Lambda);

    // 3. Calculate the Discriminant
    double D = (B * B) - (4.0 * A * C);

    double pz_nu = 0.0;

    // Scenario A: Discriminant is Positive (Real roots)
    if (D > 0) {
        double pz1 = (-B + std::sqrt(D)) / (2.0 * A);
        double pz2 = (-B - std::sqrt(D)) / (2.0 * A);

        // Pick the solution with the smallest absolute value
        if (std::abs(pz1) < std::abs(pz2)) {
            pz_nu = pz1;
        } else {
            pz_nu = pz2;
        }
    }
    // Scenario B: Discriminant is Negative (Complex roots,
    // mismeasurement)
    else {
        // Set D = 0 and just take the real part
        pz_nu = -B / (2.0 * A);
    }

    return pz_nu;
}

double solveJ2E(double j1E, double j1eta, double j1phi, double j2eta,
                double j2phi) {
    double cosTheta12 =
        std::tanh(j1eta) * std::tanh(j2eta) +
        std::cos(j1phi - j2phi) / (std::cosh(j1eta) * std::cosh(j2eta));
    double denom = 2.0 * j1E * (1.0 - cosTheta12);
    if (denom < 1e-6)
        return -1.0;
    return (W_MASS * W_MASS) / denom;
}

ROOT::Math::XYZVector unitVectorFromEtaPhi(double eta, double phi) {
    double ux = std::cos(phi) / std::cosh(eta);
    double uy = std::sin(phi) / std::cosh(eta);
    double uz = std::tanh(eta);

    return ROOT::Math::XYZVector(ux, uy, uz);
}

double solveJ4EFromTopEquality(const ROOT::Math::PxPyPzEVector &Wlep,
                               double targetTopMass, double j4eta, double j4phi,
                               double measuredJ4E, double bMass) {
    // Direction of J4
    auto n = unitVectorFromEtaPhi(j4eta, j4phi);

    double Ew = Wlep.E();
    double wx = Wlep.Px();
    double wy = Wlep.Py();
    double wz = Wlep.Pz();

    double Wmass2 = Wlep.M2();
    double mb2 = bMass * bMass;
    double Mt2 = targetTopMass * targetTopMass;

    // A = p_W dot n
    double A = wx * n.X() + wy * n.Y() + wz * n.Z();

    // Constraint:
    // (W + b)^2 = Mt^2
    //
    // W^2 + mb^2 + 2(Ew Eb - pW dot pb) = Mt^2
    //
    // with pb = sqrt(Eb^2 - mb^2) * n
    //
    // Ew Eb - A sqrt(Eb^2 - mb^2) = K

    double K = 0.5 * (Mt2 - Wmass2 - mb2);

    // Quadratic:
    // (Ew Eb - K)^2 = A^2 (Eb^2 - mb^2)
    //
    // (Ew^2 - A^2) Eb^2 - 2 Ew K Eb + (K^2 + A^2 mb^2) = 0

    double qa = Ew * Ew - A * A;
    double qb = -2.0 * Ew * K;
    double qc = K * K + A * A * mb2;

    if (std::abs(qa) < 1e-12)
        return -1.0;

    double disc = qb * qb - 4.0 * qa * qc;

    if (disc < 0.0)
        return -1.0;

    double sqrtDisc = std::sqrt(disc);

    double Esol1 = (-qb + sqrtDisc) / (2.0 * qa);
    double Esol2 = (-qb - sqrtDisc) / (2.0 * qa);

    bool ok1 = Esol1 >= bMass;
    bool ok2 = Esol2 >= bMass;

    if (!ok1 && !ok2)
        return -1.0;

    if (ok1 && !ok2)
        return Esol1;

    if (!ok1 && ok2)
        return Esol2;

    // If both are physical, choose the one closest to measured J4 energy.
    // If measuredJ4E was not provided, choose the smaller energy.
    if (measuredJ4E > 0.0) {
        return (std::abs(Esol1 - measuredJ4E) < std::abs(Esol2 - measuredJ4E))
                   ? Esol1
                   : Esol2;
    }

    return (Esol1 < Esol2) ? Esol1 : Esol2;
}
