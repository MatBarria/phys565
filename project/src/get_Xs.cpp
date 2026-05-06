
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

void get_Xs() {
    // ---------------------------------------------------------
    // 1. Define observable
    // bins: 1, 2, 3, 4, >=5 jets
    // ---------------------------------------------------------
    RooRealVar nJets("nJets", "Jet multiplicity", 0.5, 5.5);
    nJets.setBins(5);

    RooRealVar nBTags("nBTags", "b-tag multiplicity", 0.5, 2.5);
    nBTags.setBins(2);

    // ---------------------------------------------------------
    // 2. Create histograms
    // ---------------------------------------------------------
    TH2D *h_mc_sig = new TH2D(
        "h_mc_sig", "Signal template; n_valid_jets_tot; n_valid_bjets_tot", 5,
        0.5, 5.5, 2, 0.5, 2.5);

    TH2D *h_mc_bkg = new TH2D(
        "h_mc_bkg", "Background template; n_valid_jets_tot; n_valid_bjets_tot",
        5, 0.5, 5.5, 2, 0.5, 2.5);
    TH2D *h_data =
        new TH2D("h_data", "Data; n_valid_jets_tot; n_valid_bjets_tot", 5, 0.5,
                 5.5, 2, 0.5, 2.5);

    h_mc_sig->Sumw2();
    h_mc_bkg->Sumw2();
    h_data->Sumw2();

    // ---------------------------------------------------------
    // 3. Fill MC templates using your fractions
    // ---------------------------------------------------------

    double sig_fractions[5][2] = {
        {0.060687332, 0},          // 1 jet
        {0.17761989, 0.076765424}, // 2 jets
        {0.19453496, 0.13924528},  // 3 jets
        {0.12662257, 0.10835013},  // 4 jets
        {0.056224331, 0.059950072} // >=5 jets
    };

    double bkg_fractions[5][2] = {
        {0.61658524, 0},             // 1 jet
        {0.25721943, 0.013623876},   // 2 jets
        {0.075947519, 0.009301414},  // 3 jets
        {0.018317789, 0.0021972499}, // 4 jets
        {0.0047723009, 0.0020351877} // >=5 jets
    };

    // double sig_fractions[5] = {0.06348451, 0.2794353, 0.32464346, 0.20753953,
    // 0.12489715};
    // double bkg_fractions[5] = {6.2237692e-01, 2.8062779e-01, 6.7370966e-02,
    // 2.7822454e-02, 0.00180189614};

    // Arbitrary normalization, because RooHistPdf will normalize the shapes.
    // Only the shape matters here.
    double template_norm = 100.0;

    for (int iJet = 1; iJet <= 5; ++iJet) {
        for (int iB = 1; iB <= 2; ++iB) {
            double sig_content =
                sig_fractions[iJet - 1][iB - 1] * template_norm;
            double bkg_content =
                bkg_fractions[iJet - 1][iB - 1] * template_norm;

            h_mc_sig->SetBinContent(iJet, iB, sig_content);
            h_mc_bkg->SetBinContent(iJet, iB, bkg_content);

            h_mc_sig->SetBinError(iJet, iB, std::sqrt(sig_content));
            h_mc_bkg->SetBinError(iJet, iB, std::sqrt(bkg_content));
        }
    }
    // ---------------------------------------------------------
    // 4. Read data from ROOT file
    // ---------------------------------------------------------
    TFile *file_data =
        TFile::Open("../data/proccess_tuples/no_cuts/data_tuples_JES1p00.root");

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

        int jet_bin = N_valid_jets_tot;
        int btag_bin = N_valid_b_jets;

        // 5 means 5 or more jets
        if (jet_bin >= 5) {
            jet_bin = 5;
        }

        // 2 means 2 or more b-tags
        if (btag_bin >= 2) {
            btag_bin = 2;
        }

        h_data->Fill(jet_bin, btag_bin);
    }

    std::cout << "Loaded data events: " << h_data->Integral() << std::endl;

    double nData = h_data->Integral();

    RooArgList observables(nJets, nBTags);

    RooDataHist rd_sig("rd_sig", "Signal template", observables,
                       Import(*h_mc_sig));

    RooDataHist rd_bkg("rd_bkg", "Background template", observables,
                       Import(*h_mc_bkg));

    RooDataHist data("data", "Observed data", observables, Import(*h_data));

    RooHistPdf pdf_sig("pdf_sig", "Signal PDF", RooArgSet(nJets, nBTags),
                       rd_sig);

    RooHistPdf pdf_bkg("pdf_bkg", "Background PDF", RooArgSet(nJets, nBTags),
                       rd_bkg);

    // ---------------------------------------------------------
    // 6. Build extended model
    // ---------------------------------------------------------
    RooRealVar n_sig("n_sig", "Fitted signal yield", 0.5 * nData, 0.0,
                     2.0 * nData);

    RooRealVar n_bkg("n_bkg", "Fitted background yield", 0.5 * nData, 0.0,
                     2.0 * nData);

    RooAddPdf model("model", "Signal + Background",
                    RooArgList(pdf_sig, pdf_bkg), RooArgList(n_sig, n_bkg));

    // ---------------------------------------------------------
    // 7. Fit
    // ---------------------------------------------------------
    RooFitResult *fitResult =
        model.fitTo(data, Extended(true), Save(true), PrintLevel(-1));

    if (fitResult) {
        fitResult->Print("v");
    }

    // ---------------------------------------------------------
    // 8. Print result
    // ---------------------------------------------------------
    std::cout << "\n========================================================"
              << std::endl;
    std::cout << "              2D TEMPLATE FIT RESULTS                  "
              << std::endl;
    std::cout << "========================================================"
              << std::endl;
    std::cout << "Data events          : " << nData << std::endl;
    std::cout << "Fitted signal n_sig  : " << n_sig.getVal() << " +/- "
              << n_sig.getError() << std::endl;
    std::cout << "Fitted bkg n_bkg     : " << n_bkg.getVal() << " +/- "
              << n_bkg.getError() << std::endl;

    double purity = n_sig.getVal() / (n_sig.getVal() + n_bkg.getVal());

    std::cout << "Fitted purity        : " << purity << std::endl;
    std::cout << "========================================================"
              << std::endl;

    double chi2 = 0.0;
    int nBinsUsed = 0;

    double fitted_sig = n_sig.getVal();
    double fitted_bkg = n_bkg.getVal();

    double sig_total = h_mc_sig->Integral();
    double bkg_total = h_mc_bkg->Integral();

    std::cout << "\n========================================================"
              << std::endl;
    std::cout << "              BIN-BY-BIN FIT CHECK                     "
              << std::endl;
    std::cout << "========================================================"
              << std::endl;

    std::cout << "jet_bin  btag_bin   data    expected    pull" << std::endl;

    for (int iJet = 1; iJet <= h_data->GetNbinsX(); ++iJet) {
        for (int iB = 1; iB <= h_data->GetNbinsY(); ++iB) {

            double data_bin = h_data->GetBinContent(iJet, iB);

            double f_sig = h_mc_sig->GetBinContent(iJet, iB) / sig_total;
            double f_bkg = h_mc_bkg->GetBinContent(iJet, iB) / bkg_total;

            double expected = fitted_sig * f_sig + fitted_bkg * f_bkg;

            if (expected <= 0)
                continue;

            double pull = (data_bin - expected) / std::sqrt(expected);

            chi2 += (data_bin - expected) * (data_bin - expected) / expected;
            nBinsUsed++;

            std::cout << iJet << "        " << iB << "          " << data_bin
                      << "    " << expected << "    " << pull << std::endl;
        }
    }

    int nFitPars = 2;
    int ndf = nBinsUsed - nFitPars;

    std::cout << "--------------------------------------------------------"
              << std::endl;
    std::cout << "chi2      = " << chi2 << std::endl;
    std::cout << "ndf       = " << ndf << std::endl;
    std::cout << "chi2/ndf  = " << chi2 / ndf << std::endl;
    std::cout << "========================================================"
              << std::endl;

    // ---------------------------------------------------------
    // 9. Save input histograms for checking
    // ---------------------------------------------------------
    TFile *out = new TFile("fit_2D_templates.root", "RECREATE");
    h_mc_sig->Write();
    h_mc_bkg->Write();
    h_data->Write();
    if (fitResult)
        fitResult->Write("fitResult");
    out->Close();

    file_data->Close();
}
