
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TTree.h"

#include <cmath>
#include <iostream>

using namespace RooFit;

void get_Xs_2D() {
    // ---------------------------------------------------------
    // 1. Define observable
    // bins: 1, 2, 3, 4, >=5 jets
    // ---------------------------------------------------------
    RooRealVar nJets("nJets", "Jet multiplicity", 0.5, 5.5);
    nJets.setBins(5);

    // ---------------------------------------------------------
    // 2. Create histograms
    // ---------------------------------------------------------
    TH1D *h_mc_sig = new TH1D("h_mc_sig", "Signal template", 5, 0.5, 5.5);
    TH1D *h_mc_bkg = new TH1D("h_mc_bkg", "Background template", 5, 0.5, 5.5);
    TH1D *h_data = new TH1D("h_data", "Data", 5, 0.5, 5.5);

    h_mc_sig->Sumw2();
    h_mc_bkg->Sumw2();
    h_data->Sumw2();

    // ---------------------------------------------------------
    // 3. Fill MC templates using your fractions
    // ---------------------------------------------------------
    double sig_fractions[5] = {0.06348451, 0.2794353, 0.32464346, 0.20753953,
                               0.12489715};
    double bkg_fractions[5] = {6.2237692e-01, 2.8062779e-01, 6.7370966e-02,
                               2.7822454e-02, 0.00180189614};
    // Arbitrary normalization, because RooHistPdf will normalize the shapes.
    // Only the shape matters here.
    double template_norm = 100.0;

    for (Int_t i = 1; i <= 5; ++i) {
        h_mc_sig->SetBinContent(i, sig_fractions[i - 1] * template_norm);
        h_mc_bkg->SetBinContent(i, bkg_fractions[i - 1] * template_norm);

        h_mc_sig->SetBinError(i,
                              std::sqrt(sig_fractions[i - 1] * template_norm));
        h_mc_bkg->SetBinError(i,
                              std::sqrt(bkg_fractions[i - 1] * template_norm));
    }

    // ---------------------------------------------------------
    // 4. Read data from ROOT file
    // ---------------------------------------------------------
    TFile *file_data =
        TFile::Open("../data/proccess_tuples/no_cuts/data_tuples.root");

    if (!file_data || file_data->IsZombie()) {
        std::cerr << "ERROR: Could not open ./data/data.root" << std::endl;
        return;
    }

    // Change "tree_output" if your tree has a different name.
    TTree *tree = nullptr;
    file_data->GetObject("tree_output", tree);

    if (!tree) {
        std::cerr << "ERROR: Could not find tree named tree_output"
                  << std::endl;
        std::cerr << "Check the tree name inside ./data/data.root" << std::endl;
        file_data->Close();
        return;
    }

    UInt_t N_valid_jets_tot = 0;
    UInt_t NMuon_valid = 0;
    UInt_t N_valid_b_jets = 0;

    if (!tree->GetBranch("N_valid_jets_tot")) {
        std::cerr << "ERROR: Branch N_valid_jets_tot not found in tree_output"
                  << std::endl;
        file_data->Close();
        return;
    }

    tree->SetBranchAddress("N_valid_jets_tot", &N_valid_jets_tot);
    tree->SetBranchAddress("N_valid_b_jets", &N_valid_b_jets);
    tree->SetBranchAddress("NMuon_valid", &NMuon_valid);

    Long64_t nEntries = tree->GetEntries();

    std::cout << "Entries: " << nEntries << std::endl;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (N_valid_b_jets < 1 || NMuon_valid < 1)
            continue;

        Int_t jet_bin = N_valid_jets_tot;

        // Bin 5 means 5 or more jets
        if (jet_bin >= 5) {
            jet_bin = 5;
        }

        h_data->Fill(jet_bin);
    }

    std::cout << "Loaded data events: " << h_data->Integral() << std::endl;

    // ---------------------------------------------------------
    // 5. Convert histograms to RooFit objects
    // ---------------------------------------------------------
    RooDataHist rd_sig("rd_sig", "Signal template", RooArgList(nJets),
                       Import(*h_mc_sig));

    RooDataHist rd_bkg("rd_bkg", "Background template", RooArgList(nJets),
                       Import(*h_mc_bkg));

    RooDataHist data("data", "Observed data", RooArgList(nJets),
                     Import(*h_data));

    RooHistPdf pdf_sig("pdf_sig", "Signal PDF", RooArgSet(nJets), rd_sig);

    RooHistPdf pdf_bkg("pdf_bkg", "Background PDF", RooArgSet(nJets), rd_bkg);

    // ---------------------------------------------------------
    // 6. Build extended model
    // Data = n_sig * sig_shape + n_bkg * bkg_shape
    // ---------------------------------------------------------
    double nData = h_data->Integral();

    RooRealVar n_sig("n_sig", "Fitted signal yield", 0.5 * nData, 0.0,
                     2.0 * nData);

    RooRealVar n_bkg("n_bkg", "Fitted background yield", 0.5 * nData, 0.0,
                     2.0 * nData);

    RooAddPdf model("model", "Signal + Background",
                    RooArgList(pdf_sig, pdf_bkg), RooArgList(n_sig, n_bkg));

    // ---------------------------------------------------------
    // 7. Fit
    // ---------------------------------------------------------
    std::cout << "Starting binned extended template fit..." << std::endl;

    RooFitResult *fitResult =
        model.fitTo(data, Extended(true), Save(true), PrintLevel(-1));

    if (fitResult) {
        fitResult->Print("v");
    }

    // ---------------------------------------------------------
    // 8. Plot
    // ---------------------------------------------------------
    TCanvas *canvas = new TCanvas("canvas", "Jet Multiplicity Fit", 800, 600);

    RooPlot *frame = nJets.frame(Title("Template fit: N_valid_jets_tot"));

    data.plotOn(frame, Name("Data"), MarkerStyle(20));

    model.plotOn(frame, Components(pdf_bkg), Name("Background"),
                 LineColor(kRed), FillColor(kRed - 9), DrawOption("F"));

    model.plotOn(frame, Name("FullModel"), LineColor(kBlue), LineWidth(2));

    model.plotOn(frame, Components(pdf_sig), Name("Signal"), LineStyle(kDashed),
                 LineColor(kGreen + 2), LineWidth(2));

    data.plotOn(frame, Name("DataAgain"), MarkerStyle(20));

    frame->GetXaxis()->SetTitle("N_valid_jets_tot");
    frame->GetYaxis()->SetTitle("Events / bin");

    frame->Draw();

    canvas->SaveAs("fit_n_valid_jets_tot.pdf");
    canvas->SaveAs("fit_n_valid_jets_tot.png");

    // ---------------------------------------------------------
    // 9. Print result
    // ---------------------------------------------------------
    std::cout << "\n========================================================"
              << std::endl;
    std::cout << "                 TEMPLATE FIT RESULTS                   "
              << std::endl;
    std::cout << "========================================================"
              << std::endl;
    std::cout << "Data events          : " << nData << std::endl;
    std::cout << "--------------------------------------------------------"
              << std::endl;
    std::cout << "Fitted signal n_sig  : " << n_sig.getVal() << " +/- "
              << n_sig.getError() << std::endl;
    std::cout << "Fitted bkg n_bkg     : " << n_bkg.getVal() << " +/- "
              << n_bkg.getError() << std::endl;
    std::cout << "--------------------------------------------------------"
              << std::endl;

    double purity = n_sig.getVal() / (n_sig.getVal() + n_bkg.getVal());

    std::cout << "Fitted purity        : " << purity << std::endl;
    std::cout << "========================================================"
              << std::endl;

    file_data->Close();
}
