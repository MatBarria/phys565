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
const char *TREE_NAME = "tree_output";     
const char *BRANCH_MASS = "top_hadronic_mass";
const char *BRANCH_WEIGHT = "total_weight";  
const double MASS_MIN = 80.0;               
const double MASS_MAX = 380.0;             
const double MT_GEN = 172.5;              
const double MT_FIT_LO = 160.0;          
const double MT_FIT_HI = 185.0;         

void fitMass() {

    RooRealVar m_fit("m_fit", "Reconstructed Top Mass [GeV]", MASS_MIN,
                     MASS_MAX);


    RooRealVar mean("mean", "Crystal Ball Mean", 170.0, 130.0, 220.0);
    RooRealVar sigma("sigma", "Crystal Ball Sigma", 25.0, 5.0, 80.0);

    // alpha < 0 gives a right-side tail
    RooRealVar alpha("alpha", "Right-tail alpha", -1.5, -10.0, -0.05);

    // tail power
    RooRealVar n("n", "Tail exponent", 3.0, 0.5, 50.0);

    RooCBShape model("model", "Crystal Ball with right tail", m_fit, mean,
                     sigma, alpha, n);

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

    std::cout << "\n========================================================\n";
    std::cout << "                      FIT RESULTS                       \n";
    std::cout << "========================================================\n";
    std::cout << "Crystal Ball mean      : " << mean.getVal() << " +/- "
              << mean.getError() << " GeV\n";

    std::cout << "Crystal Ball sigma     : " << sigma.getVal() << " +/- "
              << sigma.getError() << " GeV\n";

    std::cout << "Alpha right tail       : " << alpha.getVal() << " +/- "
              << alpha.getError() << "\n";

    std::cout << "Tail exponent n        : " << n.getVal() << " +/- "
              << n.getError() << "\n";


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

    std::cout << "\n# Copy this into Python\n";
    std::cout << "crystal_ball_params = [" << mean.getVal() << ", "
              << sigma.getVal() << ", " << alpha.getVal() << ", " << n.getVal()
              << ", " << sumW << "]\n";

    std::cout << "crystal_ball_errors = [" << mean.getError() << ", "
              << sigma.getError() << ", " << alpha.getError() << ", "
              << n.getError() << "]\n";


    std::ofstream scanFile("mtop_likelihood_scan.txt");
    scanFile << "# mtop  nll  minus2DeltaNLL\n";

    // First find the global minimum from the nominal fit
    double nll_min = fitResult->minNll();

    const int nScan = 80;
    double mt_min = 150.0;
    double mt_max = 190.0;

    for (int i = 0; i <= nScan; ++i) {
        double mt_test = mt_min + (mt_max - mt_min) * i / nScan;

        double mean_test = mt_test;

        mean.setVal(mean_test);
        mean.setConstant(true);

        sigma.setConstant(false);
        alpha.setConstant(false);
        n.setConstant(false);

        RooFitResult *scanResult =
            model.fitTo(*data, Save(true), SumW2Error(true), PrintLevel(-1),
                        Verbose(false), Warnings(false));

        double nll = scanResult->minNll();
        double minus2DeltaNLL = 2.0 * (nll - nll_min);

        scanFile << mt_test << " " << nll << " " << minus2DeltaNLL << "\n";

        delete scanResult;
    }

    scanFile.close();

    // Restore mean
    mean.setConstant(false);
}
