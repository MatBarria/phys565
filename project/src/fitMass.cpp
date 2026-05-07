#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include <fstream>
// #include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooNovosibirsk.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooVoigtian.h"
// #include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <iostream>

using namespace RooFit;
const char *TREE_NAME = "tree_output";         // TTree name
const char *BRANCH_MASS = "top_hadronic_mass"; // fitted top mass branch
const char *BRANCH_WEIGHT = "total_weight";    // event weight branch
const double MASS_MIN = 100.0;                 // fit range low  [GeV]
const double MASS_MAX = 350.0;                 // fit range high [GeV]
const double MT_GEN = 172.5;                   // generated top mass [GeV]
const double MT_FIT_LO = 160.0;                // Mt search range low
const double MT_FIT_HI = 185.0;                // Mt search range high

void fitMass() {
    // ---------------------------------------------------------
    // 1. Define the observable
    // ---------------------------------------------------------
    // RooRealVar m_fit("m_fit", "Reconstructed Top Mass [GeV]", MASS_MIN,
    // MASS_MAX);

    //// Common peak position
    // RooRealVar mean("mean", "Common Gaussian Mean", 175.0, 140.0, 220.0);

    //// Narrow/core Gaussian
    // RooRealVar sigma1("sigma1", "Core Gaussian Width", 15.0, 3.0, 50.0);
    // RooGaussian gauss1("gauss1", "Core Gaussian", m_fit, mean, sigma1);

    //// Wide Gaussian for tails / wrong permutations
    // RooRealVar sigma2("sigma2", "Wide Gaussian Width", 45.0, 10.0, 150.0);
    // RooGaussian gauss2("gauss2", "Wide Gaussian", m_fit, mean, sigma2);

    //// Fraction of the narrow Gaussian
    // RooRealVar f_core("f_core", "Core Gaussian Fraction", 0.7, 0.0, 1.0);

    //// Combined model
    // RooAddPdf model("model", "Double Gaussian", RooArgList(gauss1, gauss2),
    // RooArgList(f_core));

    // ---------------------------------------------------------
    // 2. Breit-Wigner convolved with Gaussian = Voigtian
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // 1. Define observable
    // ---------------------------------------------------------
    // RooRealVar m_fit("m_fit", "Reconstructed Top Mass [GeV]", MASS_MIN,
    // MASS_MAX);

    //// Peak position
    // RooRealVar mean("mean", "Voigtian Mean", 175.0, 140.0, 220.0);

    //// Breit-Wigner width
    //// Start small/moderate. The real top width is ~1.4 GeV,
    //// but your reconstructed width will mostly come from resolution.
    // RooRealVar width("width", "Breit-Wigner Width", 5.0, 0.1, 50.0);

    //// Gaussian detector/reconstruction resolution
    // RooRealVar sigma("sigma", "Gaussian Resolution", 15.0, 1.0, 80.0);

    //// Voigtian model
    // RooVoigtian model("model", "Voigtian: Breit-Wigner #otimes Gaussian",
    // m_fit, mean, width, sigma);

    //// ---------------------------------------------------------
    RooRealVar m_fit("m_fit", "Reconstructed Top Mass [GeV]", MASS_MIN,
                     MASS_MAX);

    //---------------------------------------------------------
    // 2. Crystal Ball with right tail
    //---------------------------------------------------------

    RooRealVar mean("mean", "Crystal Ball Mean", 170.0, 130.0, 220.0);
    RooRealVar sigma("sigma", "Crystal Ball Sigma", 25.0, 5.0, 80.0);

    // alpha < 0 gives a right-side tail
    RooRealVar alpha("alpha", "Right-tail alpha", -1.5, -10.0, -0.05);

    // tail power
    RooRealVar n("n", "Tail exponent", 3.0, 0.5, 50.0);

    RooCBShape model("model", "Crystal Ball with right tail", m_fit, mean,
                     sigma, alpha, n);

    // ---------------------------------------------------------
    // 1. Observable
    // ---------------------------------------------------------
    // RooRealVar m_fit("m_fit", "Reconstructed Top Mass [GeV]", MASS_MIN,
    // MASS_MAX);

    // ---------------------------------------------------------
    // 2. Landau convolved with Gaussian
    // ---------------------------------------------------------

    // Landau most probable value.
    // This is not exactly the same as the mean.
    // RooRealVar mpv("mpv", "Landau MPV", 165.0, 120.0, 230.0);

    //// Landau width controls the asymmetric right tail.
    // RooRealVar landau_width("landau_width", "Landau Width", 20.0, 2.0,
    // 100.0);

    // RooLandau landau("landau", "Landau", m_fit, mpv, landau_width);

    //// Gaussian resolution.
    //// For a convolution kernel centered at zero, define a second observable
    //// variable.
    // RooRealVar gauss_mean("gauss_mean", "Gaussian Mean", 0.0);
    // gauss_mean.setConstant(true);

    // RooRealVar gauss_sigma("gauss_sigma", "Gaussian Sigma", 10.0, 1.0, 80.0);

    // RooGaussian gauss("gauss", "Gaussian Resolution", m_fit, gauss_mean,
    // gauss_sigma);

    //// Important for FFT convolution: use enough bins in the observable cache.
    // m_fit.setBins(10000, "cache");

    // RooFFTConvPdf model("model", "Landau #otimes Gaussian", m_fit, landau,
    // gauss);

    // TFile *file =
    // TFile::Open("../data/proccess_tuples/ttbar_tuples_JES1p00.root");

    TFile *file =
        TFile::Open("../data/proccess_tuples/data_tuples_JES1p00.root");

    // TFile *file = TFile::Open("../data/proccess_tuples/data_tuples.root");
    //  TFile *file =
    //  TFile::Open("../data/proccess_tuples/data_tuples_JES1p03.root"); TFile
    //  *file = TFile::Open("../data/proccess_tuples/signal.root");
    TTree *tree = (TTree *)file->Get("tree_output");
    // RooDataSet *data =
    // new RooDataSet("data", "dataset from tree", tree, RooArgSet(m_fit));

    float mass = 0., weight = 1.;
    tree->SetBranchAddress(BRANCH_MASS, &mass);
    tree->SetBranchAddress(BRANCH_WEIGHT, &weight);

    // RooDataSet needs a weight variable declared in the arg set
    // RooRealVar w("w", "event weight", 1., -1e15, 1e10);
    RooRealVar w("w", "event weight", 1., 0, 1e15);

    // RooDataSet *data = new RooDataSet("my_data_set", "my_data_set",
    // RooArgSet(x, w), RooFit::WeightVar(w));
    RooDataSet *data =
        new RooDataSet("my_data_set", "my_data_set", RooArgSet(m_fit, w),
                       RooFit::WeightVar(w));

    Long64_t nEntries = tree->GetEntries();
    Long64_t nPass = 0;

    double sumW = 0.;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (mass < MASS_MIN || mass > MASS_MAX)
            continue;
        m_fit.setVal(mass);
        w.setVal(weight);
        data->add(RooArgSet(m_fit, w), weight);
        sumW += weight;
        ++nPass;
    }

    std::cout << "[loadDataset] " << "my_data_set" << ": " << nPass
              << " events, sum(w) = " << sumW << "\n";

    // ---------------------------------------------------------
    // 6. Perform the Unbinned Maximum Likelihood Fit
    // ---------------------------------------------------------
    std::cout << "Starting unbinned fit..." << std::endl;

    RooFitResult *fitResult = model.fitTo(
        *data, Save(true),
        SumW2Error(true), // correcto para datos pesados
        // Minos(true), // errores asimétricos via perfil de likelihood
        PrintLevel(1));

    TCanvas *canvas = new TCanvas("canvas", "Mass Fit", 800, 600);

    RooPlot *frame =
        m_fit.frame(Title("Unbinned ML Fit: Double gaussian"),
                    RooFit::Range(MASS_MIN, MASS_MAX), RooFit::Bins(25));

    // bars as sqrt(sum(w^2)) instead of sqrt(N).
    data->plotOn(frame, RooFit::Name("Data"),
                 RooFit::DataError(RooAbsData::SumW2));

    // Plot the full model (Signal + Background)
    model.plotOn(frame, RooFit::Name("FullModel"), RooFit::LineColor(kBlue));

    // Draw everything on the canvas
    frame->Draw();

    // std::cout <<
    // "\n========================================================\n"; std::cout
    // << "                      FIT RESULTS                       \n";
    // std::cout <<
    // "========================================================\n";

    // std::cout << "Common peak mean       : " << mean.getVal() << " +/- "
    //<< mean.getError() << " GeV" << std::endl;

    // std::cout << "Core width sigma1      : " << sigma1.getVal() << " +/- "
    //<< sigma1.getError() << " GeV" << std::endl;

    // std::cout << "Wide width sigma2      : " << sigma2.getVal() << " +/- "
    //<< sigma2.getError() << " GeV" << std::endl;

    // std::cout << "Core fraction          : " << f_core.getVal() << " +/- "
    //<< f_core.getError() << std::endl;

    // std::cout <<
    // "========================================================\n";
    // fitResult->Print("v");

    // TString outPrefix = "double_gaussian_fit";

    //// ---------------------------------------------------------
    //// 1. Save fit parameters
    //// ---------------------------------------------------------
    //// ---------------------------------------------------------
    //// Save double Gaussian fit information for Python plotting
    //// ---------------------------------------------------------

    // std::ofstream fitInfo("../data/double_gaussian_fit_results.csv");

    // fitInfo << "parameter,value,error\n";

    // fitInfo << "mass_peak," << mean.getVal() << "," << mean.getError() <<
    // "\n";

    // fitInfo << "sigma1," << sigma1.getVal() << "," << sigma1.getError() <<
    // "\n";

    // fitInfo << "sigma2," << sigma2.getVal() << "," << sigma2.getError() <<
    // "\n";

    // fitInfo << "f_core," << f_core.getVal() << "," << f_core.getError() <<
    // "\n";

    // fitInfo << "fit_status," << fitResult->status() << ",0\n";

    // fitInfo << "cov_quality," << fitResult->covQual() << ",0\n";

    // fitInfo << "minNll," << fitResult->minNll() << ",0\n";

    // fitInfo.close();

    //// ---------------------------------------------------------
    //// Save fitted curve only
    //// ---------------------------------------------------------

    // std::ofstream curveFile("../data/double_gaussian_fit_curve.csv");

    // curveFile << "mass,total,gauss1,gauss2\n";

    // const int nPoints = 1000;
    // const int nBinsForPlot = 40; // use same binning as your Python histogram
    // const double binWidth = (MASS_MAX - MASS_MIN) / double(nBinsForPlot);

    //// This should match the total weighted entries used in the fit.
    //// If your dataset weights are already correct, sumW is okay.
    // const double norm = sumW;

    // for (int i = 0; i < nPoints; ++i) {
    // double xval =
    // MASS_MIN + i * (MASS_MAX - MASS_MIN) / double(nPoints - 1);

    // m_fit.setVal(xval);

    // double totalPdf = model.getVal(RooArgSet(m_fit));
    // double g1Pdf = gauss1.getVal(RooArgSet(m_fit));
    // double g2Pdf = gauss2.getVal(RooArgSet(m_fit));

    // double totalY = totalPdf * norm * binWidth;
    // double g1Y = f_core.getVal() * g1Pdf * norm * binWidth;
    // double g2Y = (1.0 - f_core.getVal()) * g2Pdf * norm * binWidth;

    // curveFile << std::setprecision(10) << xval << "," << totalY << ","
    //<< g1Y << "," << g2Y << "\n";
    //}

    // curveFile.close();

    std::cout << "\n========================================================\n";
    std::cout << "                      FIT RESULTS                       \n";
    std::cout << "========================================================\n";

    // std::cout << "Voigtian mean          : " << mean.getVal() << " +/- "
    //<< mean.getError() << " GeV" << std::endl;

    // std::cout << "Breit-Wigner width     : " << width.getVal() << " +/- "
    //<< width.getError() << " GeV" << std::endl;

    // std::cout << "Gaussian sigma         : " << sigma.getVal() << " +/- "
    //<< sigma.getError() << " GeV" << std::endl;

    // std::cout <<
    // "========================================================\n";
    // fitResult->Print("v");

    std::cout << "Crystal Ball mean      : " << mean.getVal() << " +/- "
              << mean.getError() << " GeV\n";

    std::cout << "Crystal Ball sigma     : " << sigma.getVal() << " +/- "
              << sigma.getError() << " GeV\n";

    std::cout << "Alpha right tail       : " << alpha.getVal() << " +/- "
              << alpha.getError() << "\n";

    std::cout << "Tail exponent n        : " << n.getVal() << " +/- "
              << n.getError() << "\n";

    // std::cout <<
    // "\n========================================================\n"; std::cout
    // << "                      FIT RESULTS                       \n";
    // std::cout <<
    // "========================================================\n";

    // std::cout << "Landau MPV            : " << mpv.getVal() << " +/- "
    //<< mpv.getError() << " GeV\n";

    // std::cout << "Landau width          : " << landau_width.getVal() << " +/-
    // "
    //<< landau_width.getError() << " GeV\n";

    // std::cout << "Gaussian sigma        : " << gauss_sigma.getVal() << " +/-
    // "
    //<< gauss_sigma.getError() << " GeV\n";

    // std::cout <<
    // "========================================================\n";

    fitResult->Print("v");
    std::cout << "\n# Python parameters for Crystal Ball fit\n";

    std::cout << "crystal_ball_params = [\n";

    std::cout << "    " << mean.getVal() << ", " << mean.getError()
              << ",  # mean, mean_err\n";

    std::cout << "    " << sigma.getVal() << ", " << sigma.getError()
              << ",  # sigma, sigma_err\n";

    std::cout << "    " << alpha.getVal() << ", " << alpha.getError()
              << ",  # alpha, alpha_err\n";

    std::cout << "    " << n.getVal() << ", " << n.getError()
              << ",  # n, n_err\n";

    std::cout << "    " << sumW << ", 0.0,  # sumW, sumW_err\n";

    std::cout << "    " << fitResult->status() << ", 0.0,  # fit_status\n";

    std::cout << "    " << fitResult->covQual() << ", 0.0,  # cov_quality\n";

    std::cout << "    " << fitResult->minNll() << ", 0.0,  # minNll\n";

    std::cout << "]\n";
}
