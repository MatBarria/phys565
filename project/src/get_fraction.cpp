#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

#include <iomanip>
#include <iostream>
#include <string>

void fillHistograms(const std::string &filename, TH1D *h_jets, TH1D *h_btags,
                    TH2D *h_2d, const std::string &treeName = "tree_output") {
    TFile *file = TFile::Open(filename.c_str());

    if (!file || file->IsZombie()) {
        std::cerr << "ERROR: could not open file " << filename << std::endl;
        return;
    }

    TTree *tree = nullptr;
    file->GetObject(treeName.c_str(), tree);

    if (!tree) {
        std::cerr << "ERROR: could not find tree " << treeName << " in file "
                  << filename << std::endl;
        file->Close();
        return;
    }

    UInt_t N_valid_jets_tot = 0;
    UInt_t N_valid_b_jets = 0;
    UInt_t NMuon_valid = 0;
    UChar_t triggerIsoMu24 = true;
    float weight = 0;

    if (!tree->GetBranch("N_valid_jets_tot")) {
        std::cerr << "ERROR: branch N_valid_jets_tot not found in " << filename
                  << std::endl;
        file->Close();
        return;
    }

    if (!tree->GetBranch("N_valid_b_jets")) {
        std::cerr << "ERROR: branch N_valid_b_jets not found in " << filename
                  << std::endl;
        file->Close();
        return;
    }

    tree->SetBranchAddress("N_valid_jets_tot", &N_valid_jets_tot);
    tree->SetBranchAddress("N_valid_b_jets", &N_valid_b_jets);
    tree->SetBranchAddress("NMuon_valid", &NMuon_valid);
    tree->SetBranchAddress("triggerIsoMu24", &triggerIsoMu24);
    tree->SetBranchAddress("weight", &weight);

    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (N_valid_jets_tot == 0)
            continue;
        if (N_valid_b_jets == 0)
            continue;
        if (NMuon_valid == 0)
            continue;
        if (!triggerIsoMu24) {
            // std::cout << filename << std::endl;
            continue;
        }

        int jet_bin = N_valid_jets_tot;
        int btag_bin = N_valid_b_jets;

        // 5 means 5 or more jets
        if (jet_bin >= 5) {
            jet_bin = 5;
        }

        // 2 means 2 or more b-tagged jets
        if (btag_bin >= 2) {
            btag_bin = 2;
        }

        h_jets->Fill(jet_bin, weight);
        h_btags->Fill(btag_bin, weight);
        h_2d->Fill(jet_bin, btag_bin, weight);
    }

    file->Close();
}

void print1DFractions(TH1D *h, const std::string &label) {
    double total = h->Integral();

    std::cout << "\n" << label << " 1D fractions " << total << ":" << std::endl;

    if (total <= 0) {
        std::cout << "Histogram is empty." << std::endl;
        return;
    }

    std::cout << "{";
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double frac = h->GetBinContent(i) / total;
        std::cout << std::setprecision(8) << frac;

        if (i != h->GetNbinsX()) {
            std::cout << ", ";
        }
    }
    std::cout << "}" << std::endl;
}

void print2DFractions(TH2D *h, const std::string &label) {
    double total = h->Integral();

    std::cout << "\n"
              << label << " 2D fractions [5][2] " << total << " :" << std::endl;

    if (total <= 0) {
        std::cout << "Histogram is empty." << std::endl;
        return;
    }

    std::cout << "{" << std::endl;

    for (int iJet = 1; iJet <= h->GetNbinsX(); ++iJet) {
        std::cout << "    {";

        for (int iB = 1; iB <= h->GetNbinsY(); ++iB) {
            double frac = h->GetBinContent(iJet, iB) / total;
            std::cout << std::setprecision(8) << frac;

            if (iB != h->GetNbinsY()) {
                std::cout << ", ";
            }
        }

        std::cout << "}";

        if (iJet != h->GetNbinsX()) {
            std::cout << ",";
        }

        if (iJet == 1)
            std::cout << "  // 1 jet";
        if (iJet == 2)
            std::cout << "  // 2 jets";
        if (iJet == 3)
            std::cout << "  // 3 jets";
        if (iJet == 4)
            std::cout << "  // 4 jets";
        if (iJet == 5)
            std::cout << "  // >=5 jets";

        std::cout << std::endl;
    }

    std::cout << "};" << std::endl;
}

void get_fraction() {
    // ---------------------------------------------------------
    // Histograms for signal
    // ---------------------------------------------------------
    TH1D *h_sig_jets = new TH1D(
        "h_sig_jets", "Signal jet multiplicity; N_valid_jets_tot; Events", 5,
        0.5, 5.5);

    TH1D *h_sig_btags = new TH1D(
        "h_sig_btags", "Signal b-tag multiplicity; N_valid_b_jets; Events", 2,
        0.5, 2.5);

    TH2D *h_sig_2d = new TH2D(
        "h_sig_2d", "Signal 2D template; N_valid_jets_tot; N_valid_b_jets", 5,
        0.5, 5.5, 2, 0.5, 2.5);

    // ---------------------------------------------------------
    // Histograms for background
    // ---------------------------------------------------------
    TH1D *h_bkg_jets = new TH1D(
        "h_bkg_jets", "Background jet multiplicity; N_valid_jets_tot; Events",
        5, 0.5, 5.5);

    TH1D *h_bkg_btags = new TH1D(
        "h_bkg_btags", "Background b-tag multiplicity; N_valid_b_jets; Events",
        2, 0.5, 2.5);

    TH2D *h_bkg_2d = new TH2D(
        "h_bkg_2d", "Background 2D template; N_valid_jets_tot; N_valid_b_jets",
        5, 0.5, 5.5, 2, 0.5, 2.5);

    // ---------------------------------------------------------
    // Fill from ROOT files
    // ---------------------------------------------------------
    fillHistograms("../data/proccess_tuples/no_cuts/signal.root", h_sig_jets,
                   h_sig_btags, h_sig_2d);
    fillHistograms("../data/proccess_tuples/no_cuts/background.root",
                   h_bkg_jets, h_bkg_btags, h_bkg_2d);

    // ---------------------------------------------------------
    // Print total events
    // ---------------------------------------------------------
    std::cout << "\nSignal entries after selection: " << h_sig_2d->Integral()
              << std::endl;

    std::cout << "Background entries after selection: " << h_bkg_2d->Integral()
              << std::endl;

    // ---------------------------------------------------------
    // Print 1D fractions
    // ---------------------------------------------------------
    print1DFractions(h_sig_jets, "Signal jets");
    print1DFractions(h_bkg_jets, "Background jets");

    print1DFractions(h_sig_btags, "Signal b-tags");
    print1DFractions(h_bkg_btags, "Background b-tags");

    // ---------------------------------------------------------
    // Print 2D fractions
    // ---------------------------------------------------------
    print2DFractions(h_sig_2d, "Signal");
    print2DFractions(h_bkg_2d, "Background");

    // ---------------------------------------------------------
    // Save histograms for checking
    // ---------------------------------------------------------
    TFile *out = new TFile("fractions_check.root", "RECREATE");

    h_sig_jets->Write();
    h_sig_btags->Write();
    h_sig_2d->Write();

    h_bkg_jets->Write();
    h_bkg_btags->Write();
    h_bkg_2d->Write();

    out->Close();

    std::cout << "\nSaved histograms to fractions_check.root" << std::endl;
}
