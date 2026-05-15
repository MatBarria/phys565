#include "../lib/CreateTuples.h"
#include "../lib/kinematics.h"
#include <iostream>
#include <ostream>
CreateTuple::CreateTuple(TString input, TString output, TString channel,
                         float JES_input) {

    tree_input = new TChain("events");
    tree_output = new TTree("tree_output", "tree_output");
    input_name = input;
    output_directory = output;
    // output_name = channel + "_tuples.root";
    JES = JES_input;
    TString JES_string = Form("%.2f", JES);
    JES_string.ReplaceAll(".", "p");

    output_name =
        Form("%s_tuples_JES%s.root", channel.Data(), JES_string.Data());

    fillChain();
}

CreateTuple::~CreateTuple() {
    delete tree_input->GetCurrentFile();
    delete tree_output;
    delete output_file;
}

void CreateTuple::fillChain() {
    tree_input->Add(input_name);
    if (!tree_input) {
        throw std::runtime_error("Tuple file not found");
        return;
    }
    std::cout << "Tuple loaded" << std::endl;
}

void CreateTuple::setBranchesAddressesOutput() {

    tree_output->Branch("NMuon", &NMuon, "NMuon/i");
    tree_output->Branch("NJet", &NJet, "NJet/i");
    tree_output->Branch("NPhoton", &NPhoton, "NPhoton/i");
    tree_output->Branch("NElectron", &NElectron, "NElectron/i");
    tree_output->Branch("MET_px", &MET_px, "MET_px/f");
    tree_output->Branch("MET_py", &MET_py, "MET_py/f");
    tree_output->Branch("MET_pt", &MET_pt, "MET_pt/f");
    tree_output->Branch("MChadronicBottom_px", &MChadronicBottom_px,
                        "MChadronicBottom_px/f");
    tree_output->Branch("MChadronicBottom_py", &MChadronicBottom_py,
                        "MChadronicBottom_py/f");
    tree_output->Branch("MChadronicBottom_pz", &MChadronicBottom_pz,
                        "MChadronicBottom_pz/f");
    tree_output->Branch("MCleptonicBottom_px", &MCleptonicBottom_px,
                        "MCleptonicBottom_px/f");
    tree_output->Branch("MCleptonicBottom_py", &MCleptonicBottom_py,
                        "MCleptonicBottom_py/f");
    tree_output->Branch("MCleptonicBottom_pz", &MCleptonicBottom_pz,
                        "MCleptonicBottom_pz/f");
    tree_output->Branch("MChadronicWDecayQuark_px", &MChadronicWDecayQuark_px,
                        "MChadronicWDecayQuark_px/f");
    tree_output->Branch("MChadronicWDecayQuark_py", &MChadronicWDecayQuark_py,
                        "MChadronicWDecayQuark_py/f");
    tree_output->Branch("MChadronicWDecayQuark_pz", &MChadronicWDecayQuark_pz,
                        "MChadronicWDecayQuark_pz/f");
    tree_output->Branch("MChadronicWDecayQuarkBar_px",
                        &MChadronicWDecayQuarkBar_px,
                        "MChadronicWDecayQuarkBar_px/f");
    tree_output->Branch("MChadronicWDecayQuarkBar_py",
                        &MChadronicWDecayQuarkBar_py,
                        "MChadronicWDecayQuarkBar_py/f");
    tree_output->Branch("MChadronicWDecayQuarkBar_pz",
                        &MChadronicWDecayQuarkBar_pz,
                        "MChadronicWDecayQuarkBar_pz/f");
    tree_output->Branch("MClepton_px", &MClepton_px, "MClepton_px/f");
    tree_output->Branch("MClepton_py", &MClepton_py, "MClepton_py/f");
    tree_output->Branch("MClepton_pz", &MClepton_pz, "MClepton_pz/f");
    tree_output->Branch("MCleptonPDGid", &MCleptonPDGid, "MCleptonPDGid/I");
    tree_output->Branch("MCneutrino_px", &MCneutrino_px, "MCneutrino_px/f");
    tree_output->Branch("MCneutrino_py", &MCneutrino_py, "MCneutrino_py/f");
    tree_output->Branch("MCneutrino_pz", &MCneutrino_pz, "MCneutrino_pz/f");
    tree_output->Branch("NPrimaryVertices", &NPrimaryVertices,
                        "NPrimaryVertices/i");
    tree_output->Branch("triggerIsoMu24", &triggerIsoMu24, "triggerIsoMu24/b");

    tree_output->Branch("MCtop_mass_leptonic", &MCtop_mass_leptonic,
                        "MCtop_mass_leptonic/f");
    tree_output->Branch("MCtop_mass_hadronic", &MCtop_mass_hadronic,
                        "MCtop_mass_hadronic/f");
    tree_output->Branch("N_valid_b_jets", &N_valid_b_jets, "N_valid_b_jets/i");
    tree_output->Branch("N_valid_jets", &N_valid_jets, "N_valid_jets/i");
    tree_output->Branch("N_valid_jets_tot", &N_valid_jets_tot,
                        "N_valid_jets_tot/i");
    tree_output->Branch("diMuon_mass", &diMuon_mass, "diMuon_mass/f");
    tree_output->Branch("weight", &weight, "weight/f");
    tree_output->Branch("permutation_weight_sum", &permutation_weight_sum,
                        "permutation_weight_sum/f");
    tree_output->Branch("total_weight", &total_weight, "total_weight/f");
    tree_output->Branch("total_weight_norm", &total_weight_norm,
                        "total_weight_norm/f");
    tree_output->Branch("chi2", &chi2, "chi2/f");
    tree_output->Branch("permutation_weight", &permutation_weight,
                        "permutation_weight/f");

    tree_output->Branch("W_leptonic_mass", &W_leptonic_mass,
                        "W_leptonic_mass/f");
    tree_output->Branch("W_hadronic_mass", &W_hadronic_mass,
                        "W_hadronic_mass/f");
    tree_output->Branch("W_leptonic_mass_reco", &W_leptonic_mass_reco,
                        "W_leptonic_mass_reco/f");
    tree_output->Branch("W_hadronic_mass_reco", &W_hadronic_mass_reco,
                        "W_hadronic_mass_reco/f");
    tree_output->Branch("top_hadronic_mass_1", &top_hadronic_mass_1,
                        "top_hadronic_mass_1/f");
    tree_output->Branch("top_hadronic_mass_2", &top_hadronic_mass_2,
                        "top_hadronic_mass_2/f");
    tree_output->Branch("top_leptoninc_mass_1", &top_leptoninc_mass_1,
                        "top_leptoninc_mass_1/f");
    tree_output->Branch("top_leptoninc_mass_2", &top_leptoninc_mass_2,
                        "top_leptoninc_mass_2/f");
    tree_output->Branch("top_hadronic_mass", &top_hadronic_mass,
                        "top_hadronic_mass/f");
    tree_output->Branch("top_hadronic_mass_reco", &top_hadronic_mass_reco,
                        "top_hadronic_mass_reco/f");
    tree_output->Branch("top_leptoninc_mass", &top_leptoninc_mass,
                        "top_leptoninc_mass/f");
    tree_output->Branch("mu1_Px", &mu1_Px, "mu1_Px/f");
    tree_output->Branch("mu1_Py", &mu1_Py, "mu1_Py/f");
    tree_output->Branch("mu1_Pz", &mu1_Pz, "mu1_Pz/f");
    tree_output->Branch("mu1_E", &mu1_E, "mu1_E/f");
    tree_output->Branch("mu1_Iso", &mu1_Iso, "mu1_Iso/f");
    tree_output->Branch("mu1_Pt", &mu1_Pt, "mu1_Pt/f");
    tree_output->Branch("mu2_Px", &mu2_Px, "mu2_Px/f");
    tree_output->Branch("mu2_Py", &mu2_Py, "mu2_Py/f");
    tree_output->Branch("mu2_Pz", &mu2_Pz, "mu2_Pz/f");
    tree_output->Branch("mu2_E", &mu2_E, "mu2_E/f");
    tree_output->Branch("mu2_Iso", &mu2_Iso, "mu2_Iso/f");
    tree_output->Branch("mu2_Pt", &mu2_Pt, "mu2_Pt/f");
    tree_output->Branch("mu3_Px", &mu3_Px, "mu3_Px/f");
    tree_output->Branch("mu3_Py", &mu3_Py, "mu3_Py/f");
    tree_output->Branch("mu3_Pz", &mu3_Pz, "mu3_Pz/f");
    tree_output->Branch("mu3_E", &mu3_E, "mu3_E/f");
    tree_output->Branch("mu3_Iso", &mu3_Iso, "mu3_Iso/f");
    tree_output->Branch("mu3_Pt", &mu3_Pt, "mu3_Pt/f");
    tree_output->Branch("e1_Px", &e1_Px, "e1_Px/f");
    tree_output->Branch("e1_Py", &e1_Py, "e1_Py/f");
    tree_output->Branch("e1_Pz", &e1_Pz, "e1_Pz/f");
    tree_output->Branch("e1_E", &e1_E, "e1_E/f");
    tree_output->Branch("e1_Iso", &e1_Iso, "e1_Iso/f");
    tree_output->Branch("e1_Pt", &e1_Pt, "e1_Pt/f");
    tree_output->Branch("e2_Px", &e2_Px, "e2_Px/f");
    tree_output->Branch("e2_Py", &e2_Py, "e2_Py/f");
    tree_output->Branch("e2_Pz", &e2_Pz, "e2_Pz/f");
    tree_output->Branch("e2_E", &e2_E, "e2_E/f");
    tree_output->Branch("NMuon_valid", &NMuon_valid, "NMuon_valid/i");
    tree_output->Branch("Ne_valid", &Ne_valid, "Ne_valid/i");
    tree_output->Branch("Nlep_valid", &Nlep_valid, "Nlep_valid/i");
    tree_output->Branch("NMuon_valid_mc", &NMuon_valid_mc, "NMuon_valid_mc/i");
    tree_output->Branch("Ne_valid_mc", &Ne_valid_mc, "Ne_valid_mc/i");
    tree_output->Branch("Nlep_valid_mc", &Nlep_valid_mc, "Nlep_valid_mc/i");

    tree_output->Branch("jet1_Px", &jet_Px[0], "jet1_Px/F");
    tree_output->Branch("jet1_Py", &jet_Py[0], "jet1_Py/F");
    tree_output->Branch("jet1_Pz", &jet_Pz[0], "jet1_Pz/F");
    tree_output->Branch("jet1_E", &jet_E[0], "jet1_E/F");
    tree_output->Branch("jet1_Pt", &jet_Pt[0], "jet1_Pt/F");

    tree_output->Branch("jet2_Px", &jet_Px[1], "jet2_Px/F");
    tree_output->Branch("jet2_Py", &jet_Py[1], "jet2_Py/F");
    tree_output->Branch("jet2_Pz", &jet_Pz[1], "jet2_Pz/F");
    tree_output->Branch("jet2_E", &jet_E[1], "jet2_E/F");
    tree_output->Branch("jet2_Pt", &jet_Pt[1], "jet2_Pt/F");

    tree_output->Branch("jet3_Px", &jet_Px[2], "jet3_Px/F");
    tree_output->Branch("jet3_Py", &jet_Py[2], "jet3_Py/F");
    tree_output->Branch("jet3_Pz", &jet_Pz[2], "jet3_Pz/F");
    tree_output->Branch("jet3_E", &jet_E[2], "jet3_E/F");
    tree_output->Branch("jet3_Pt", &jet_Pt[2], "jet3_Pt/F");

    tree_output->Branch("jet4_Px", &jet_Px[3], "jet4_Px/F");
    tree_output->Branch("jet4_Py", &jet_Py[3], "jet4_Py/F");
    tree_output->Branch("jet4_Pz", &jet_Pz[3], "jet4_Pz/F");
    tree_output->Branch("jet4_E", &jet_E[3], "jet4_E/F");
    tree_output->Branch("jet4_Pt", &jet_Pt[3], "jet4_Pt/F");

    tree_output->Branch("jet5_Px", &jet_Px[4], "jet5_Px/F");
    tree_output->Branch("jet5_Py", &jet_Py[4], "jet5_Py/F");
    tree_output->Branch("jet5_Pz", &jet_Pz[4], "jet5_Pz/F");
    tree_output->Branch("jet5_E", &jet_E[4], "jet5_E/F");
    tree_output->Branch("jet5_Pt", &jet_Pt[4], "jet5_Pt/F");
    tree_output->Branch("jet6_Px", &jet_Px[5], "jet6_Px/F");
    tree_output->Branch("jet6_Py", &jet_Py[5], "jet6_Py/F");
    tree_output->Branch("jet6_Pz", &jet_Pz[5], "jet6_Pz/F");
    tree_output->Branch("jet6_E", &jet_E[5], "jet6_E/F");
    tree_output->Branch("jet6_Pt", &jet_Pt[5], "jet6_Pt/F");

    tree_output->Branch("bjet1_Px", &bjet_Px[0], "bjet1_Px/F");
    tree_output->Branch("bjet1_Py", &bjet_Py[0], "bjet1_Py/F");
    tree_output->Branch("bjet1_Pz", &bjet_Pz[0], "bjet1_Pz/F");
    tree_output->Branch("bjet1_E", &bjet_E[0], "bjet1_E/F");
    tree_output->Branch("bjet1_Pt", &bjet_Pt[0], "bjet1_Pt/F");

    tree_output->Branch("bjet2_Px", &bjet_Px[1], "bjet2_Px/F");
    tree_output->Branch("bjet2_Py", &bjet_Py[1], "bjet2_Py/F");
    tree_output->Branch("bjet2_Pz", &bjet_Pz[1], "bjet2_Pz/F");
    tree_output->Branch("bjet2_E", &bjet_E[1], "bjet2_E/F");
    tree_output->Branch("bjet2_Pt", &bjet_Pt[1], "bjet2_Pt/F");

    tree_output->Branch("bjet3_Px", &bjet_Px[2], "bjet3_Px/F");
    tree_output->Branch("bjet3_Py", &bjet_Py[2], "bjet3_Py/F");
    tree_output->Branch("bjet3_Pz", &bjet_Pz[2], "bjet3_Pz/F");
    tree_output->Branch("bjet3_E", &bjet_E[2], "bjet3_E/F");
    tree_output->Branch("bjet3_Pt", &bjet_Pt[2], "bjet3_Pt/F");

    tree_output->Branch("bjet4_Px", &bjet_Px[3], "bjet4_Px/F");
    tree_output->Branch("bjet4_Py", &bjet_Py[3], "bjet4_Py/F");
    tree_output->Branch("bjet4_Pz", &bjet_Pz[3], "bjet4_Pz/F");
    tree_output->Branch("bjet4_E", &bjet_E[3], "bjet4_E/F");
    tree_output->Branch("bjet4_Pt", &bjet_Pt[3], "bjet4_Pt/F");

    tree_output->Branch("bjet1_Btag", &bjet_Btag[0], "bjet1_Btag/f");
    tree_output->Branch("bjet2_Btag", &bjet_Btag[1], "bjet2_Btag/f");
    tree_output->Branch("bjet3_Btag", &bjet_Btag[2], "bjet3_Btag/f");
    tree_output->Branch("bjet4_Btag", &bjet_Btag[3], "bjet4_Btag/f");
    std::cout << "Output Branches addressed setted" << std::endl;
}

void CreateTuple::setBranchesAddressesInput() {

    tree_input->SetBranchAddress("NJet", &NJet);
    tree_input->SetBranchAddress("Jet_Px", &Jet_Px);
    tree_input->SetBranchAddress("Jet_Py", &Jet_Py);
    tree_input->SetBranchAddress("Jet_Pz", &Jet_Pz);
    tree_input->SetBranchAddress("Jet_E", &Jet_E);
    tree_input->SetBranchAddress("Jet_btag", &Jet_btag);
    tree_input->SetBranchAddress("Jet_ID", &Jet_ID);
    tree_input->SetBranchAddress("NMuon", &NMuon);
    tree_input->SetBranchAddress("Muon_Px", &Muon_Px);
    tree_input->SetBranchAddress("Muon_Py", &Muon_Py);
    tree_input->SetBranchAddress("Muon_Pz", &Muon_Pz);
    tree_input->SetBranchAddress("Muon_E", &Muon_E);
    tree_input->SetBranchAddress("Muon_Charge", &Muon_Charge);
    tree_input->SetBranchAddress("Muon_Iso", &Muon_Iso);
    tree_input->SetBranchAddress("NElectron", &NElectron);
    tree_input->SetBranchAddress("Electron_Px", &Electron_Px);
    tree_input->SetBranchAddress("Electron_Py", &Electron_Py);
    tree_input->SetBranchAddress("Electron_Pz", &Electron_Pz);
    tree_input->SetBranchAddress("Electron_E", &Electron_E);
    tree_input->SetBranchAddress("Electron_Charge", &Electron_Charge);
    tree_input->SetBranchAddress("Electron_Iso", &Electron_Iso);
    tree_input->SetBranchAddress("NPhoton", &NPhoton);
    tree_input->SetBranchAddress("Photon_Px", &Photon_Px);
    tree_input->SetBranchAddress("Photon_Py", &Photon_Py);
    tree_input->SetBranchAddress("Photon_Pz", &Photon_Pz);
    tree_input->SetBranchAddress("Photon_E", &Photon_E);
    tree_input->SetBranchAddress("Photon_Iso", &Photon_Iso);
    tree_input->SetBranchAddress("MET_px", &MET_px);
    tree_input->SetBranchAddress("MET_py", &MET_py);
    tree_input->SetBranchAddress("MChadronicBottom_px", &MChadronicBottom_px);
    tree_input->SetBranchAddress("MChadronicBottom_py", &MChadronicBottom_py);
    tree_input->SetBranchAddress("MChadronicBottom_pz", &MChadronicBottom_pz);
    tree_input->SetBranchAddress("MCleptonicBottom_px", &MCleptonicBottom_px);
    tree_input->SetBranchAddress("MCleptonicBottom_py", &MCleptonicBottom_py);
    tree_input->SetBranchAddress("MCleptonicBottom_pz", &MCleptonicBottom_pz);
    tree_input->SetBranchAddress("MChadronicWDecayQuark_px",
                                 &MChadronicWDecayQuark_px);
    tree_input->SetBranchAddress("MChadronicWDecayQuark_py",
                                 &MChadronicWDecayQuark_py);
    tree_input->SetBranchAddress("MChadronicWDecayQuark_pz",
                                 &MChadronicWDecayQuark_pz);
    tree_input->SetBranchAddress("MChadronicWDecayQuarkBar_px",
                                 &MChadronicWDecayQuarkBar_px);
    tree_input->SetBranchAddress("MChadronicWDecayQuarkBar_py",
                                 &MChadronicWDecayQuarkBar_py);
    tree_input->SetBranchAddress("MChadronicWDecayQuarkBar_pz",
                                 &MChadronicWDecayQuarkBar_pz);
    tree_input->SetBranchAddress("MClepton_px", &MClepton_px);
    tree_input->SetBranchAddress("MClepton_py", &MClepton_py);
    tree_input->SetBranchAddress("MClepton_pz", &MClepton_pz);
    tree_input->SetBranchAddress("MCleptonPDGid", &MCleptonPDGid);
    tree_input->SetBranchAddress("MCneutrino_px", &MCneutrino_px);
    tree_input->SetBranchAddress("MCneutrino_py", &MCneutrino_py);
    tree_input->SetBranchAddress("MCneutrino_pz", &MCneutrino_pz);
    tree_input->SetBranchAddress("NPrimaryVertices", &NPrimaryVertices);
    tree_input->SetBranchAddress("triggerIsoMu24", &triggerIsoMu24);
    tree_input->SetBranchAddress("EventWeight", &EventWeight);

    std::cout << "Input Branches addressed setted" << std::endl;
}

void CreateTuple::fillOutputTree(TString channel) {

    std::cout << "Filling Output Tree" << std::endl;
    Long64_t total_entries = tree_input->GetEntries();
    std::cout << "Entries: " << total_entries << std::endl;

    float MC_lepton_mass;

    for (int event_index = 0; event_index < total_entries; event_index++) {
        tree_input->GetEntry(event_index);

        weight = EventWeight;

        // if (triggerIsoMu24 == 0) {
        // continue;
        //}

        MET_pt = std::hypot(MET_px, MET_py);

        diMuon_mass = -1.0;
        if (NMuon >= 2) {
            bool found_pair = false;
            for (int i = 0; i < NMuon && !found_pair; ++i) {
                for (int j = i + 1; j < NMuon && !found_pair; ++j) {
                    if (Muon_Charge[i] * Muon_Charge[j] == -1) {
                        ROOT::Math::PxPyPzEVector mu1(Muon_Px[i], Muon_Py[i],
                                                      Muon_Pz[i], Muon_E[i]);
                        ROOT::Math::PxPyPzEVector mu2(Muon_Px[j], Muon_Py[j],
                                                      Muon_Pz[j], Muon_E[j]);
                        diMuon_mass = (mu1 + mu2).M();
                        found_pair = true;
                    }
                }
            }
        }

        mu1_Px = mu1_Py = mu1_Pz = mu1_E = mu1_Iso = mu1_Pt = -1;
        mu2_Px = mu2_Py = mu2_Pz = mu2_E = mu2_Iso = mu2_Pt = -1;
        mu3_Px = mu3_Py = mu3_Pz = mu3_E = mu3_Iso = mu3_Pt = -1;
        mu1_Iso = -2;
        e1_Px = e1_Py = e1_Pz = e1_E = e1_Iso = e1_Pt = -1;
        e2_Px = e2_Py = e2_Pz = e2_E = e2_Iso = e2_Pt = -1;
        NMuon_valid = Ne_valid = Nlep_valid = 0;

        // auto fill_muon = [&](int i, float &px, float &py, float &pz, float
        // &e, float &iso, float &pt, int &n_valid) {
        // if (NMuon > i) {
        // if (std::hypot(Muon_Px[i], Muon_Py[i]) > 25 &&
        // Muon_Iso[i] < 2.625) {
        //// if (pt > 25) {
        // n_valid = n_valid + 1;
        // px = Muon_Px[i];
        // py = Muon_Py[i];
        // pz = Muon_Pz[i];
        // e = Muon_E[i];
        // iso = Muon_Iso[i];
        // pt = std::hypot(px, py);
        //}
        //}
        //};

        auto fill_electron = [&](int i, float &px, float &py, float &pz,
                                 float &e, float &iso, float &pt,
                                 int &n_valid) {
            if (NElectron > i) {
                px = Electron_Px[i];
                py = Electron_Py[i];
                pz = Electron_Pz[i];
                e = Electron_E[i];
                iso = Electron_Iso[i];
                pt = std::hypot(px, py);
                if (pt > 25 && iso < 2.6) {
                    // if (pt > 25) {
                    n_valid = n_valid + 1;
                }
            }
        };

        auto save_muon = [&](int out_index, int in_index) {
            if (out_index == 0) {
                mu1_Px = Muon_Px[in_index];
                mu1_Py = Muon_Py[in_index];
                mu1_Pz = Muon_Pz[in_index];
                mu1_E = Muon_E[in_index];
                mu1_Iso = Muon_Iso[in_index];
                mu1_Pt = std::hypot(mu1_Px, mu1_Py);
            } else if (out_index == 1) {
                mu2_Px = Muon_Px[in_index];
                mu2_Py = Muon_Py[in_index];
                mu2_Pz = Muon_Pz[in_index];
                mu2_E = Muon_E[in_index];
                mu2_Iso = Muon_Iso[in_index];
                mu2_Pt = std::hypot(mu2_Px, mu2_Py);
            } else if (out_index == 2) {
                mu3_Px = Muon_Px[in_index];
                mu3_Py = Muon_Py[in_index];
                mu3_Pz = Muon_Pz[in_index];
                mu3_E = Muon_E[in_index];
                mu3_Iso = Muon_Iso[in_index];
                mu3_Pt = std::hypot(mu3_Px, mu3_Py);
            }
        };

        for (int i = 0; i < NMuon; i++) {
            float pt = std::hypot(Muon_Px[i], Muon_Py[i]);
            float iso = Muon_Iso[i];

            if (pt > 25 && iso < 2.625) {
                if (NMuon_valid < 3) {
                    save_muon(NMuon_valid, i);
                }

                NMuon_valid++;
            }
        }

        // fill_muon(0, mu1_Px, mu1_Py, mu1_Pz, mu1_E, mu1_Iso, mu1_Pt,
        // NMuon_valid);
        // fill_muon(1, mu2_Px, mu2_Py, mu2_Pz, mu2_E, mu2_Iso, mu2_Pt,
        // NMuon_valid);
        // fill_muon(2, mu3_Px, mu3_Py, mu3_Pz, mu3_E, mu3_Iso, mu3_Pt,
        // NMuon_valid);
        fill_electron(0, e1_Px, e1_Py, e1_Pz, e1_E, e1_Iso, e1_Pt, Ne_valid);
        fill_electron(1, e2_Px, e2_Py, e2_Pz, e2_E, e2_Iso, e2_Pt, Ne_valid);

        bool ALL_EVENTS = false;
        MC_lepton_mass = -1;
        Nlep_valid = NMuon_valid + Ne_valid;
        if (NMuon_valid == 0 && !ALL_EVENTS)
            continue;
        MCtop_mass_hadronic = MCtop_mass_leptonic = -1;
        NMuon_valid_mc = Ne_valid_mc = Nlep_valid_mc = -1;

        if (channel == "ttbar") {
            NMuon_valid_mc = Ne_valid_mc = Nlep_valid_mc = 0;
            ROOT::Math::PxPyPzMVector MC_B_hadronic_vector(
                MChadronicBottom_px, MChadronicBottom_py, MChadronicBottom_pz,
                BOTTOM_MASS);
            ROOT::Math::PxPyPzMVector MC_B_leptonic_vector(
                MCleptonicBottom_px, MCleptonicBottom_py, MCleptonicBottom_pz,
                BOTTOM_MASS);

            if (std::abs(MCleptonPDGid) == 11) {
                MC_lepton_mass = E_MASS;
                Ne_valid_mc++;
                Nlep_valid_mc++;
            }
            if (std::abs(MCleptonPDGid) == 13) {
                MC_lepton_mass = MU_MASS;
                NMuon_valid_mc++;
                Nlep_valid_mc++;
            }
            if (std::abs(MCleptonPDGid) == 15) {
                MC_lepton_mass = TAU_MASS;
                Nlep_valid_mc++;
            }

            if (Nlep_valid_mc != 0) {

                ROOT::Math::PxPyPzMVector MC_lepton_vector(
                    MClepton_px, MClepton_py, MClepton_pz, MC_lepton_mass);
                ROOT::Math::PxPyPzMVector MC_nu_vector(
                    MCneutrino_px, MCneutrino_py, MCneutrino_pz, NU_MASS);

                MCtop_mass_leptonic =
                    (MC_B_leptonic_vector + MC_nu_vector + MC_lepton_vector)
                        .M();

                ROOT::Math::PxPyPzMVector MC_q_vector(
                    MChadronicWDecayQuark_px, MChadronicWDecayQuark_py,
                    MChadronicWDecayQuark_pz, Q_MASS);
                ROOT::Math::PxPyPzMVector MC_qbar_vector(
                    MChadronicWDecayQuarkBar_px, MChadronicWDecayQuarkBar_py,
                    MChadronicWDecayQuarkBar_pz, Q_MASS);

                MCtop_mass_hadronic =
                    (MC_qbar_vector + MC_q_vector + MC_B_hadronic_vector).M();
            }
        }
        N_valid_jets = N_valid_b_jets = 0;
        const float BTAG_CUT = 1.511;
        //const float BTAG_CUT = 1.7625899280575539;
        std::vector<int> valid_jets_idx;
        std::vector<int> valid_b_jets_idx;
        std::vector<int> valid_jets_tot_idx;

        // std::cout << "---------------------" << std::endl;
        for (int jet = 0; jet < NJet; ++jet) {
            if (!Jet_ID[jet])
                continue;

            valid_jets_tot_idx.push_back(jet);

            if (Jet_btag[jet] > BTAG_CUT) {
                valid_b_jets_idx.push_back(jet);
                ++N_valid_b_jets;
            } else {
                valid_jets_idx.push_back(jet);
                ++N_valid_jets;
            }
        }

        N_valid_jets_tot = N_valid_jets + N_valid_b_jets;
        for (int i = 0; i < 6; ++i) {
            jet_Px[i] = -1;
            jet_Py[i] = -1;
            jet_Pz[i] = -1;
            jet_E[i] = -1;
            jet_Pt[i] = -1;
        }

        for (int i = 0; i < 4; ++i) {
            bjet_Px[i] = -1;
            bjet_Py[i] = -1;
            bjet_Pz[i] = -1;
            bjet_E[i] = -1;
            bjet_Pt[i] = -1;
            bjet_Btag[i] = -1;
        }

        for (int i = 0; i < std::min((int)valid_jets_idx.size(), 6); ++i) {
            int j = valid_jets_idx[i];
            jet_Px[i] = Jet_Px[j];
            jet_Py[i] = Jet_Py[j];
            jet_Pz[i] = Jet_Pz[j];
            jet_E[i] = Jet_E[j];
            jet_Pt[i] = std::hypot(Jet_Px[j], Jet_Py[j]);
        }

        for (int i = 0; i < std::min((int)valid_b_jets_idx.size(), 4); ++i) {
            int j = valid_b_jets_idx[i];
            bjet_Px[i] = Jet_Px[j];
            bjet_Py[i] = Jet_Py[j];
            bjet_Pz[i] = Jet_Pz[j];
            bjet_E[i] = Jet_E[j];
            bjet_Pt[i] = std::hypot(Jet_Px[j], Jet_Py[j]);
            bjet_Btag[i] = Jet_btag[j];
        }

        W_leptonic_mass = W_hadronic_mass = top_hadronic_mass_1 =
            top_hadronic_mass_2 = top_leptoninc_mass_1 = top_leptoninc_mass_2 =
                -1;
        // if (N_valid_b_jets >= 2) {
        // N_valid_b_jets = 2;
        //}
        // if (N_valid_jets_tot >= 5) {
        // N_valid_jets_tot = 5;
        //}

        chi2 = 0;
        permutation_weight = 0;

        if (ALL_EVENTS) {
            tree_output->Fill();
            continue;
        }

        if (N_valid_jets_tot < 4) {
            // tree_output->Fill();
            continue;
        }

        if (N_valid_b_jets == 0) {
            // tree_output->Fill();
            continue;
        }
        if (MET_pt < 7) {
            continue;
        }

        if (channel == "data") {
            weight = weight * .9;
        }
        // std::cout <<"I guess we are here now " << N_valid_jets
        // <<std::endl;
        double dijet_mass_diff = 1e9;
        ROOT::Math::PxPyPzEVector W_hadronic_vector;
        if (N_valid_b_jets == 1) {
            valid_b_jets_idx.push_back(std::move(valid_jets_idx.front()));
            valid_jets_idx.erase(valid_jets_idx.begin());
        }
        if (N_valid_b_jets == 0) {
            valid_b_jets_idx.push_back(std::move(valid_jets_idx.front()));
            valid_jets_idx.erase(valid_jets_idx.begin());
            valid_b_jets_idx.push_back(std::move(valid_jets_idx.front()));
            valid_jets_idx.erase(valid_jets_idx.begin());
        }

        // ROOT::Math::PxPyPzEVector jet1();
        // ROOT::Math::PxPyPzEVector jet2();
        // int jet_1, jet_2;
        for (int idx1 = 0; idx1 < valid_jets_idx.size(); ++idx1) {
            for (int idx2 = idx1 + 1; idx2 < valid_jets_idx.size(); ++idx2) {
                int j1 = valid_jets_idx[idx1];
                int j2 = valid_jets_idx[idx2];

                ROOT::Math::PxPyPzEVector jet1_vector(Jet_Px[j1], Jet_Py[j1],
                                                      Jet_Pz[j1], Jet_E[j1]);

                ROOT::Math::PxPyPzEVector jet2_vector(Jet_Px[j2], Jet_Py[j2],
                                                      Jet_Pz[j2], Jet_E[j2]);

                ROOT::Math::PxPyPzEVector W_vector = jet1_vector + jet2_vector;
                float W_vector_mass = W_vector.M();
                float mass_diff = std::abs(W_vector_mass - W_MASS);
                if (mass_diff < dijet_mass_diff) {
                    dijet_mass_diff = mass_diff;
                    W_hadronic_mass = W_vector_mass;
                    // jet_1 = j1;
                    // jet_2 = j2;
                    W_hadronic_vector = W_vector;
                }
            }
        }
        // top_hadronic_mass
        dijet_mass_diff = 1e9;

        // if (valid_b_jets_idx.size() < 3) {
        int bj1 = valid_b_jets_idx[0];
        ROOT::Math::PxPyPzEVector b1_vector(Jet_Px[bj1], Jet_Py[bj1],
                                            Jet_Pz[bj1], Jet_E[bj1]);

        top_hadronic_mass_1 = (W_hadronic_vector + b1_vector).M();
        int bj2 = valid_b_jets_idx[1];
        ROOT::Math::PxPyPzEVector b2_vector(Jet_Px[bj2], Jet_Py[bj2],
                                            Jet_Pz[bj2], Jet_E[bj2]);
        //}

        top_hadronic_mass_1 = (W_hadronic_vector + b1_vector).M();
        top_hadronic_mass_2 = (W_hadronic_vector + b2_vector).M();

        ROOT::Math::PxPyPzEVector mu_vector(mu1_Px, mu1_Py, mu1_Pz, mu1_E);

        nu_pz = CalculateNeutrinoPz(mu1_Px, mu1_Py, mu1_Pz, MET_px, MET_py);

        ROOT::Math::PxPyPzMVector nu_vector(MET_px, MET_py, nu_pz, 0);
        W_leptonic_mass = (mu_vector + nu_vector).M();
        top_leptoninc_mass_1 = (mu_vector + nu_vector + b1_vector).M();
        top_leptoninc_mass_2 = (mu_vector + nu_vector + b2_vector).M();

        // if (std::abs(top_hadronic_mass_1 - top_leptoninc_mass_2) <
        // std::abs(top_hadronic_mass_2 - top_leptoninc_mass_2)) {
        // top_hadronic_mass = top_hadronic_mass_1;
        // top_leptoninc_mass = top_leptoninc_mass_2;
        //} else {
        // top_hadronic_mass = top_hadronic_mass_2;
        // top_leptoninc_mass = top_leptoninc_mass_1;
        //}
        // int N_permu = 0;
        permutation_weight_sum = 0;
        std::cout << "Entering the zoine" << std::endl;
        for (size_t a = 0; a < 4; ++a) {
            for (size_t b = a + 1; b < 4; ++b) {

                // These two jets are the hadronic W jets.
                // a < b avoids double-counting W jet order.

                if (valid_jets_tot_idx[a] == valid_b_jets_idx[0] ||
                    valid_jets_tot_idx[b] == valid_b_jets_idx[0]) {
                    continue;
                }
                if (N_valid_b_jets > 1) {
                    if (valid_jets_tot_idx[a] == valid_b_jets_idx[1] ||
                        valid_jets_tot_idx[b] == valid_b_jets_idx[1]) {
                        continue;
                    }
                }

                for (size_t c = 0; c < 4; ++c) {
                    if (c == a || c == b)
                        continue;

                    // This jet is b_had.

                    for (size_t d = 0; d < 4; ++d) {
                        if (d == a || d == b || d == c)
                            continue;

                        // N_permu++;
                        //  This jet is b_lep.
                        chi2 = 0;
                        permutation_weight = 0;
                        std::array<int, 4> jet_perm = {
                            valid_jets_tot_idx[a], // J1: W jet 1
                            valid_jets_tot_idx[b], // J2: W jet 2
                            valid_jets_tot_idx[c], // J3: b hadronic
                            valid_jets_tot_idx[d]  // J4: b leptonic
                        };

                        std::array<double, N_PAR> p{};
                        std::array<double, N_MEAS> x_meas{};
                        std::array<double, N_MEAS> sigma{};

                        // fillFitArrays(jet_perm, mu1_Px, mu1_Py, mu1_Pz,
                        // MET_px, MET_py, p, x_meas, sigma);
                        fillFitArrays(jet_perm, p, x_meas, sigma);

                        FitResult fit = minimizeEvent(p, x_meas, sigma);

                        if (!fit.valid)
                            continue;

                        // Use the minimized total function for ranking.
                        // This includes detector chi2 + constraint penalty.

                        // Rebuild fitted objects
                        auto mu_fit = muonFromPxPyPz(
                            fit.pfit[MU_PX], fit.pfit[MU_PY], fit.pfit[MU_PZ]);

                        auto j1_fit =
                            jetFromEetaPhi(fit.pfit[J1_E], fit.pfit[J1_ETA],
                                           fit.pfit[J1_PHI], 0.);

                        double jet2_e_fit = solveJ2E(
                            fit.pfit[J1_E], fit.pfit[J1_ETA], fit.pfit[J1_PHI],
                            fit.pfit[J2_ETA], fit.pfit[J2_PHI]);
                        if (jet2_e_fit < 0.0)
                            continue;
                        auto j2_fit = jetFromEetaPhi(
                            jet2_e_fit, fit.pfit[J2_ETA], fit.pfit[J2_PHI], 0.);

                        auto j3_fit =
                            jetFromEetaPhi(fit.pfit[J3_E], fit.pfit[J3_ETA],
                                           fit.pfit[J3_PHI], BOTTOM_MASS);

                        // auto j4_fit =
                        // jetFromEetaPhi(fit.pfit[J4_E], fit.pfit[J4_ETA],
                        // fit.pfit[J4_PHI], W_MASS);

                        double nu_pz_fit = CalculateNeutrinoPz(
                            fit.pfit[MU_PX], fit.pfit[MU_PY], fit.pfit[MU_PZ],
                            fit.pfit[MET_X], fit.pfit[MET_Y]);
                        auto nu_fit = neutrinoFromMet(
                            fit.pfit[MET_X], fit.pfit[MET_Y], nu_pz_fit);

                        auto W_lep_fit = mu_fit + nu_fit;
                        auto W_had_fit = j1_fit + j2_fit;
                        auto top_had_fit = W_had_fit + j3_fit;
                        double jet4_e_fit = solveJ4EFromTopEquality(
                            W_lep_fit, top_had_fit.M(), fit.pfit[J4_ETA],
                            fit.pfit[J4_PHI], x_meas[X_J4_E], BOTTOM_MASS);
                        if (jet4_e_fit < BOTTOM_MASS)
                            continue;
                        auto j4_fit =
                            jetFromEetaPhi(jet4_e_fit, fit.pfit[J4_ETA],
                                           fit.pfit[J4_PHI], BOTTOM_MASS);

                        auto top_lep_fit = W_lep_fit + j4_fit;
                        ROOT::Math::PxPyPzEVector j1_reco(
                            Jet_Px[jet_perm[0]], Jet_Py[jet_perm[0]],
                            Jet_Pz[jet_perm[0]], Jet_E[jet_perm[0]]);

                        ROOT::Math::PxPyPzEVector j2_reco(
                            Jet_Px[jet_perm[1]], Jet_Py[jet_perm[1]],
                            Jet_Pz[jet_perm[1]], Jet_E[jet_perm[1]]);

                        // auto mu_reco = muonFromPxPyPz(mu1_Px, mu1_Py,
                        // mu1_Pz); double nu_pz_reco = CalculateNeutrinoPz(
                        // mu1_Px, mu1_Py, mu1_Pz, MET_px, MET_py);
                        // auto nu_reco =
                        // neutrinoFromMet(MET_px, MET_py, nu_pz_reco);

                        chi2 = chi2Diagonal(fit.pfit.data(), x_meas, sigma);
                        if (std::abs(top_had_fit.M() - top_lep_fit.M()) > 1) {
                            continue;
                        }
                        if (chi2 > 15) {
                            // std::cout << "are we here?" << std::endl;
                            continue;
                        }
                        permutation_weight = std::exp(-0.5 * chi2);
                        if (permutation_weight < 0.01)
                            continue;

                        permutation_weight_sum =
                            permutation_weight_sum + permutation_weight;
                        // total_weight = weight * permutation_weight;
                        //  auto W_had_reco = W_had_fit;
                        // auto W_had_reco = j1_reco + j2_reco;
                        // auto W_lep_reco = mu_reco + nu_reco;

                        // ROOT::Math::PxPyPzEVector j3_reco(
                        // Jet_Px[jet_perm[2]], Jet_Py[jet_perm[2]],
                        // Jet_Pz[jet_perm[2]], Jet_E[jet_perm[2]]);
                        // auto top_reco = W_had_reco + j3_reco;

                        // W_hadronic_mass_reco = W_had_reco.M();
                        // W_leptonic_mass_reco = W_lep_reco.M();
                        // top_hadronic_mass_reco = top_reco.M();
                        //  fit.mW_had_fit = W_had_fit.M();
                        //  std::cout
                        //<< "target top mass = " << top_had_fit.M()
                        //<< " | solved j4E = " << jet4_e_fit
                        //<< " | top_lep mass = " << top_lep_fit.M()
                        //<< " | diff = " << top_had_fit.M() -
                        //  top_lep_fit.M()
                        //<< std::endl; // fit.mW_lep_fit = W_lep_fit.M();

                        // top_hadronic_mass = top_had_fit.M();
                        // top_leptoninc_mass = top_lep_fit.M();
                        //  fit.mt_fit = 0.5 * (fit.mt_had_fit +
                        //  fit.mt_lep_fit);

                        // tree_output->Fill();
                        // allFits.push_back(fit);
                    }
                }
            }
        } // std::cout << "xi = " << chi2 << std::endl;

        for (size_t a = 0; a < 4; ++a) {
            for (size_t b = a + 1; b < 4; ++b) {

                // These two jets are the hadronic W jets.
                // a < b avoids double-counting W jet order.
                if (valid_jets_tot_idx[a] == valid_b_jets_idx[0] ||
                    valid_jets_tot_idx[b] == valid_b_jets_idx[0]) {
                    continue;
                }
                if (N_valid_b_jets > 1) {
                    if (valid_jets_tot_idx[a] == valid_b_jets_idx[1] ||
                        valid_jets_tot_idx[b] == valid_b_jets_idx[1]) {
                        continue;
                    }
                }

                for (size_t c = 0; c < 4; ++c) {
                    if (c == a || c == b)
                        continue;

                    // This jet is b_had.

                    for (size_t d = 0; d < 4; ++d) {
                        if (d == a || d == b || d == c)
                            continue;

                        // N_permu++;
                        //  This jet is b_lep.
                        chi2 = 0;
                        permutation_weight = 0;
                        std::array<int, 4> jet_perm = {
                            valid_jets_tot_idx[a], // J1: W jet 1
                            valid_jets_tot_idx[b], // J2: W jet 2
                            valid_jets_tot_idx[c], // J3: b hadronic
                            valid_jets_tot_idx[d]  // J4: b leptonic
                        };

                        std::array<double, N_PAR> p{};
                        std::array<double, N_MEAS> x_meas{};
                        std::array<double, N_MEAS> sigma{};

                        // fillFitArrays(jet_perm, mu1_Px, mu1_Py, mu1_Pz,
                        // MET_px, MET_py, p, x_meas, sigma);
                        fillFitArrays(jet_perm, p, x_meas, sigma);

                        FitResult fit = minimizeEvent(p, x_meas, sigma);

                        if (!fit.valid)
                            continue;

                        // Use the minimized total function for ranking.
                        // This includes detector chi2 + constraint penalty.

                        // Rebuild fitted objects
                        auto mu_fit = muonFromPxPyPz(
                            fit.pfit[MU_PX], fit.pfit[MU_PY], fit.pfit[MU_PZ]);

                        auto j1_fit =
                            jetFromEetaPhi(fit.pfit[J1_E], fit.pfit[J1_ETA],
                                           fit.pfit[J1_PHI], 0.);

                        double jet2_e_fit = solveJ2E(
                            fit.pfit[J1_E], fit.pfit[J1_ETA], fit.pfit[J1_PHI],
                            fit.pfit[J2_ETA], fit.pfit[J2_PHI]);
                        if (jet2_e_fit < 0.0)
                            continue;
                        auto j2_fit = jetFromEetaPhi(
                            jet2_e_fit, fit.pfit[J2_ETA], fit.pfit[J2_PHI], 0.);

                        auto j3_fit =
                            jetFromEetaPhi(fit.pfit[J3_E], fit.pfit[J3_ETA],
                                           fit.pfit[J3_PHI], BOTTOM_MASS);

                        // auto j4_fit =
                        // jetFromEetaPhi(fit.pfit[J4_E], fit.pfit[J4_ETA],
                        // fit.pfit[J4_PHI], W_MASS);

                        double nu_pz_fit = CalculateNeutrinoPz(
                            fit.pfit[MU_PX], fit.pfit[MU_PY], fit.pfit[MU_PZ],
                            fit.pfit[MET_X], fit.pfit[MET_Y]);
                        auto nu_fit = neutrinoFromMet(
                            fit.pfit[MET_X], fit.pfit[MET_Y], nu_pz_fit);

                        auto W_lep_fit = mu_fit + nu_fit;
                        auto W_had_fit = j1_fit + j2_fit;
                        auto top_had_fit = W_had_fit + j3_fit;
                        double jet4_e_fit = solveJ4EFromTopEquality(
                            W_lep_fit, top_had_fit.M(), fit.pfit[J4_ETA],
                            fit.pfit[J4_PHI], x_meas[X_J4_E], BOTTOM_MASS);
                        if (jet4_e_fit < BOTTOM_MASS)
                            continue;
                        auto j4_fit =
                            jetFromEetaPhi(jet4_e_fit, fit.pfit[J4_ETA],
                                           fit.pfit[J4_PHI], BOTTOM_MASS);

                        auto top_lep_fit = W_lep_fit + j4_fit;
                        ROOT::Math::PxPyPzEVector j1_reco(
                            Jet_Px[jet_perm[0]], Jet_Py[jet_perm[0]],
                            Jet_Pz[jet_perm[0]], Jet_E[jet_perm[0]]);

                        ROOT::Math::PxPyPzEVector j2_reco(
                            Jet_Px[jet_perm[1]], Jet_Py[jet_perm[1]],
                            Jet_Pz[jet_perm[1]], Jet_E[jet_perm[1]]);

                        auto mu_reco = muonFromPxPyPz(mu1_Px, mu1_Py, mu1_Pz);
                        double nu_pz_reco = CalculateNeutrinoPz(
                            mu1_Px, mu1_Py, mu1_Pz, MET_px, MET_py);
                        auto nu_reco =
                            neutrinoFromMet(MET_px, MET_py, nu_pz_reco);

                        chi2 = chi2Diagonal(fit.pfit.data(), x_meas, sigma);
                        if (std::abs(top_had_fit.M() - top_lep_fit.M()) > 1) {
                            continue;
                        }
                        if (chi2 > 10) {
                            // std::cout << "are we here?" << std::endl;
                            continue;
                        }
                        permutation_weight = std::exp(-0.5 * chi2);
                        if (permutation_weight < 0.01)
                            continue;

                        total_weight = weight * permutation_weight;
                        total_weight_norm =
                            total_weight / permutation_weight_sum;
                        // auto W_had_reco = W_had_fit;
                        auto W_had_reco = j1_reco + j2_reco;
                        auto W_lep_reco = mu_reco + nu_reco;

                        ROOT::Math::PxPyPzEVector j3_reco(
                            Jet_Px[jet_perm[2]], Jet_Py[jet_perm[2]],
                            Jet_Pz[jet_perm[2]], Jet_E[jet_perm[2]]);

                        auto top_reco = W_had_reco + j3_reco;

                        W_hadronic_mass_reco = W_had_reco.M();
                        W_leptonic_mass_reco = W_lep_reco.M();
                        top_hadronic_mass_reco = top_reco.M();
                        // fit.mW_had_fit = W_had_fit.M();
                        // std::cout
                        //<< "target top mass = " << top_had_fit.M()
                        //<< " | solved j4E = " << jet4_e_fit
                        //<< " | top_lep mass = " << top_lep_fit.M()
                        //<< " | diff = " << top_had_fit.M() -
                        // top_lep_fit.M()
                        //<< std::endl; // fit.mW_lep_fit = W_lep_fit.M();

                        top_hadronic_mass = top_had_fit.M();
                        top_leptoninc_mass = top_lep_fit.M();
                        // fit.mt_fit = 0.5 * (fit.mt_had_fit +
                        // fit.mt_lep_fit);

                        tree_output->Fill();
                        //  allFits.push_back(fit);
                    }
                }
            }
        } // std::cout << "xi = " << chi2 << std::endl;
    }
}

void CreateTuple::saveTree() {
    std::cout << "Events in the final tuples: " << tree_output->GetEntries()
              << std::endl;
    output_file = new TFile(output_directory + output_name, "RECREATE");
    output_file->cd();

    tree_output->Write();

    output_file->Close();
    std::cout << "Tuple saved" << std::endl;
}

void CreateTuple::fillFitArrays(const std::array<int, 4> &jet_perm,

                                // double mu1_Px, double mu1_Py, double
                                // mu1_Pz, double MET_px, double MET_py,

                                std::array<double, N_PAR> &p,
                                std::array<double, N_MEAS> &x_meas,
                                std::array<double, N_MEAS> &sigma) {
    double p_mu =
        std::sqrt(mu1_Px * mu1_Px + mu1_Py * mu1_Py + mu1_Pz * mu1_Pz);

    double sigmaMu = muSigma(p_mu);

    // --------------------
    // Muon
    // --------------------
    p[MU_PX] = mu1_Px;
    p[MU_PY] = mu1_Py;
    p[MU_PZ] = mu1_Pz;

    x_meas[X_MU_PX] = mu1_Px;
    x_meas[X_MU_PY] = mu1_Py;
    x_meas[X_MU_PZ] = mu1_Pz;

    sigma[X_MU_PX] = sigmaMu;
    sigma[X_MU_PY] = sigmaMu;
    sigma[X_MU_PZ] = sigmaMu;

    // --------------------
    // Jets
    // jet_perm[0], jet_perm[1] -> W hadronic jets
    // jet_perm[2] -> b hadronic
    // jet_perm[3] -> b leptonic
    double E[4];
    double eta[4];
    double phi[4];
    double sigmaE[4]; // --------------------
    for (size_t i = 0; i < 4; ++i) {
        int idx = jet_perm[i];

        ROOT::Math::PxPyPzEVector jet_vector(
            Jet_Px[idx] * JES, Jet_Py[idx] * JES, Jet_Pz[idx] * JES,
            Jet_E[idx] * JES);

        E[i] = Jet_E[idx] * JES;
        eta[i] = jet_vector.Eta();
        phi[i] = jet_vector.Phi();

        // std::cout << "JES: " << JES << std::endl;
        sigmaE[i] = jetSigmaE(E[i], eta[i]);
    }

    double sigmaEta = jetSigmaEta();
    double sigmaPhi = jetSigmaPhi();

    // --------------------
    // Fit parameters
    // --------------------

    // J1: W light jet 1, energy is free
    p[J1_E] = E[0];
    p[J1_ETA] = eta[0];
    p[J1_PHI] = phi[0];

    // J2: W light jet 2, energy is NOT free
    // J2_E will be solved from M(j1+j2)=MW
    p[J2_ETA] = eta[1];
    p[J2_PHI] = phi[1];

    // J3: b hadronic
    p[J3_E] = E[2];
    p[J3_ETA] = eta[2];
    p[J3_PHI] = phi[2];

    // J4: b leptonic
    // p[J4_E] = E[3];
    p[J4_ETA] = eta[3];
    p[J4_PHI] = phi[3];

    // --------------------
    // Measured values
    // --------------------

    x_meas[X_J1_E] = E[0];
    x_meas[X_J1_ETA] = eta[0];
    x_meas[X_J1_PHI] = phi[0];

    x_meas[X_J2_E] = E[1];
    x_meas[X_J2_ETA] = eta[1];
    x_meas[X_J2_PHI] = phi[1];

    x_meas[X_J3_E] = E[2];
    x_meas[X_J3_ETA] = eta[2];
    x_meas[X_J3_PHI] = phi[2];

    x_meas[X_J4_E] = E[3];
    x_meas[X_J4_ETA] = eta[3];
    x_meas[X_J4_PHI] = phi[3];

    // --------------------
    // Jet uncertainties
    // --------------------

    sigma[X_J1_E] = sigmaE[0];
    sigma[X_J1_ETA] = sigmaEta;
    sigma[X_J1_PHI] = sigmaPhi;

    sigma[X_J2_E] = sigmaE[1];
    sigma[X_J2_ETA] = sigmaEta;
    sigma[X_J2_PHI] = sigmaPhi;

    sigma[X_J3_E] = sigmaE[2];
    sigma[X_J3_ETA] = sigmaEta;
    sigma[X_J3_PHI] = sigmaPhi;

    sigma[X_J4_E] = sigmaE[3];
    sigma[X_J4_ETA] = sigmaEta;
    sigma[X_J4_PHI] = sigmaPhi;
    // --------------------
    // MET
    // --------------------
    p[MET_X] = MET_px;
    p[MET_Y] = MET_py;

    x_meas[X_MET_X] = MET_px;
    x_meas[X_MET_Y] = MET_py;

    sigma[X_MET_X] = 12.0;
    sigma[X_MET_Y] = 12.0;

    // --------------------
    // Neutrino pz
    // --------------------
    // p[NU_PZ] = 0.0;
}
