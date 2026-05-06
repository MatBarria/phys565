
#include "RooAddPdf.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooFFTConvPdf.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooLandau.h"
#include "RooNovosibirsk.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooVoigtian.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
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

void fit_voi() {
    // ---------------------------------------------------------
    // 1. Define the observable
    // ---------------------------------------------------------
    RooRealVar m_fit("m_fit", "Reconstructed Top Mass [GeV]", MASS_MIN,
                     MASS_MAX);

    // ---------------------------------------------------------
    // 2. Breit-Wigner convolved with Gaussian = Voigtian
    // ---------------------------------------------------------

    // Peak position: this is your reconstructed top-mass estimator
    RooRealVar voigt_mean("voigt_mean", "Voigtian Mean", 172.5, 140.0, 220.0);

    // Breit-Wigner width.
    // For the real top, Gamma_t ~ 1.4 GeV, but your reconstructed distribution
    // is much wider, so for a phenomenological fit we can let it float.
    RooRealVar bw_width("bw_width", "Breit-Wigner Width", 10.0, 0.5, 80.0);

    // Gaussian detector/reconstruction resolution
    RooRealVar gauss_sigma("gauss_sigma", "Gaussian Resolution", 20.0, 3.0,
                           80.0);

    // Voigtian = Breit-Wigner ⊗ Gaussian
    RooVoigtian voigt("voigt", "Breit-Wigner convolved with Gaussian", m_fit,
                      voigt_mean, bw_width, gauss_sigma);

    RooAbsPdf &model = voigt;

    // TFile *file =
    // TFile::Open("../data/proccess_tuples/ttbar_tuples_JES1p00.root");
    
     TFile *file =
     TFile::Open("../data/proccess_tuples/data_tuples_JES1p00.root");

    // TFile *file = TFile::Open("../data/proccess_tuples/data_tuples.root");
    // TFile *file =
    // TFile::Open("../data/proccess_tuples/data_tuples_JES1p03.root");
    // TFile
    //   *file = TFile::Open("../data/proccess_tuples/signal.root");
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
        m_fit.frame(Title("Unbinned ML Fit: voy"), RooFit::Range(100, MASS_MAX),
                    RooFit::Bins(25));

    // bars as sqrt(sum(w^2)) instead of sqrt(N).
    data->plotOn(frame, RooFit::Name("Data"),
                 RooFit::DataError(RooAbsData::SumW2));

    // Plot the full model (Signal + Background)
    model.plotOn(frame, RooFit::Name("FullModel"), RooFit::LineColor(kBlue));

    // Draw everything on the canvas
    frame->Draw();

    std::cout << "\n========================================================\n";
    std::cout << "                      FIT RESULTS                       \n";
    std::cout << "========================================================\n";

    std::cout << "Voigtian peak mean      : " << voigt_mean.getVal() << " +/- "
              << voigt_mean.getError() << " GeV" << std::endl;

    std::cout << "Breit-Wigner width      : " << bw_width.getVal() << " +/- "
              << bw_width.getError() << " GeV" << std::endl;

    std::cout << "Gaussian resolution     : " << gauss_sigma.getVal() << " +/- "
              << gauss_sigma.getError() << " GeV" << std::endl;

    std::cout << "========================================================\n";
    fitResult->Print("v");
}
