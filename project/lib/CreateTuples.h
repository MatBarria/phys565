#ifndef LIB_CreateTuple_H
#define LIB_CreateTuple_H

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TTree.h>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <math.h>
#include <ostream>

class CreateTuple {
  public:
    CreateTuple(TString input, TString output, TString channel);

    virtual ~CreateTuple();

    void fillChain();
    double CalculateNeutrinoPz(double mu_E, double mu_px, double mu_py,
                               double mu_pz, double nu_px, double nu_py);

    void setBranchesAddressesInput();

    void setBranchesAddressesOutput();

    void fillOutputTree(TString channel);

    /**
     * @brief Saves the output tree to a ROOT file in the specified output
     * directory.
     */
    void saveTree();
    const float BOTTOM_MASS = 4.183;
    const float TOP_MASS = 173;
    const float W_MASS = 80.36;
    const float NU_MASS = 0.;
    const float E_MASS = 0.000511;
    const float MU_MASS = 0.10566;
    const float TAU_MASS = 1.77693;
    const float Q_MASS = 0;

  private:
    TChain *tree_input;       /**< Pointer to input TChain. */
    TTree *tree_output;       /**< Pointer to output TTree. */
    TString input_name;       /**< Path of the input file. */
    TString output_name;      /**< Name of the output file. */
    TString output_directory; /**< Directory for saving the output file. */
    TFile *output_file;       /**< Pointer to output TFile. */

    Int_t NJet;
    Float_t Jet_Px[20];   //[NJet]
    Float_t Jet_Py[20];   //[NJet]
    Float_t Jet_Pz[20];   //[NJet]
    Float_t Jet_E[20];    //[NJet]
    Float_t Jet_btag[20]; //[NJet]
    Bool_t Jet_ID[20];    //[NJet]
    Int_t NMuon;
    Float_t Muon_Px[3];   //[NMuon]
    Float_t Muon_Py[3];   //[NMuon]
    Float_t Muon_Pz[3];   //[NMuon]
    Float_t Muon_E[3];    //[NMuon]
    Int_t Muon_Charge[3]; //[NMuon]
    Float_t Muon_Iso[3];  //[NMuon]
    Int_t NElectron;
    Float_t Electron_Px[2];   //[NElectron]
    Float_t Electron_Py[2];   //[NElectron]
    Float_t Electron_Pz[2];   //[NElectron]
    Float_t Electron_E[2];    //[NElectron]
    Int_t Electron_Charge[2]; //[NElectron]
    Float_t Electron_Iso[2];  //[NElectron]
    Int_t NPhoton;
    Float_t Photon_Px[2];  //[NPhoton]
    Float_t Photon_Py[2];  //[NPhoton]
    Float_t Photon_Pz[2];  //[NPhoton]
    Float_t Photon_E[2];   //[NPhoton]
    Float_t Photon_Iso[2]; //[NPhoton]
    Float_t MET_px;
    Float_t MET_py;
    Float_t MChadronicBottom_px;
    Float_t MChadronicBottom_py;
    Float_t MChadronicBottom_pz;
    Float_t MCleptonicBottom_px;
    Float_t MCleptonicBottom_py;
    Float_t MCleptonicBottom_pz;
    Float_t MChadronicWDecayQuark_px;
    Float_t MChadronicWDecayQuark_py;
    Float_t MChadronicWDecayQuark_pz;
    Float_t MChadronicWDecayQuarkBar_px;
    Float_t MChadronicWDecayQuarkBar_py;
    Float_t MChadronicWDecayQuarkBar_pz;
    Float_t MClepton_px;
    Float_t MClepton_py;
    Float_t MClepton_pz;
    Int_t MCleptonPDGid;
    Float_t MCneutrino_px;
    Float_t MCneutrino_py;
    Float_t MCneutrino_pz;
    Int_t NPrimaryVertices;
    Bool_t triggerIsoMu24;
    Float_t EventWeight;

    /** New Event variables */
    float diMuon_mass;

    float mu1_Px, mu1_Py, mu1_Pz, mu1_E, mu1_Iso, mu1_Pt;
    float mu2_Px, mu2_Py, mu2_Pz, mu2_E, mu2_Iso, mu2_Pt;
    float mu3_Px, mu3_Py, mu3_Pz, mu3_E, mu3_Iso, mu3_Pt;
    float e1_Px, e1_Py, e1_Pz, e1_E, e1_Iso, e1_Pt;
    float e2_Px, e2_Py, e2_Pz, e2_E, e2_Iso, e2_Pt;
    int NMuon_valid, Ne_valid, Nlep_valid;
    int NMuon_valid_mc, Ne_valid_mc, Nlep_valid_mc;

    float jet_Px[6] = {-1}, jet_Py[6] = {-1}, jet_Pz[6] = {-1}, jet_E[6] = {0},
          jet_Pt[6] = {-1};
    float bjet_Px[4] = {-1}, bjet_Py[4] = {-1}, bjet_Pz[4] = {-1},
          bjet_E[4] = {-1}, bjet_Pt[4] = {-1}, bjet_Btag[4] = {-1};
    int N_valid_jets, N_valid_b_jets, N_valid_jets_tot;

    float MET_pt;

    float weight;
    float MCtop_mass_hadronic, MCtop_mass_leptonic;
    float W_leptonic_mass, W_hadronic_mass, top_hadronic_mass_1,
        top_hadronic_mass_2, top_leptoninc_mass_1, top_leptoninc_mass_2;
};

CreateTuple::CreateTuple(TString input, TString output, TString channel) {

    tree_input = new TChain("events");
    tree_output = new TTree("tree_output", "tree_output");
    input_name = input;
    output_directory = output;
    output_name = channel + "_tuples.root";

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

    tree_output->Branch("W_leptonic_mass", &W_leptonic_mass,
                        "W_leptonic_mass/f");
    tree_output->Branch("W_hadronic_mass", &W_hadronic_mass,
                        "W_hadronic_mass/f");
    tree_output->Branch("top_hadronic_mass_1", &top_hadronic_mass_1,
                        "top_hadronic_mass_1/f");
    tree_output->Branch("top_hadronic_mass_2", &top_hadronic_mass_2,
                        "top_hadronic_mass_2/f");
    tree_output->Branch("top_leptoninc_mass_1", &top_leptoninc_mass_1,
                        "top_leptoninc_mass_1/f");
    tree_output->Branch("top_leptoninc_mass_2", &top_leptoninc_mass_2,
                        "top_leptoninc_mass_2/f");
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

        //if (triggerIsoMu24 == 0) {
            //continue;
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
        e1_Px = e1_Py = e1_Pz = e1_E = e1_Iso = e1_Pt = -1;
        e2_Px = e2_Py = e2_Pz = e2_E = e2_Iso = e2_Pt = -1;
        NMuon_valid = Ne_valid = Nlep_valid = 0;

        auto fill_muon = [&](int i, float &px, float &py, float &pz, float &e,
                             float &iso, float &pt, int &n_valid) {
            if (NMuon > i) {
                px = Muon_Px[i];
                py = Muon_Py[i];
                pz = Muon_Pz[i];
                e = Muon_E[i];
                iso = Muon_Iso[i];
                pt = std::hypot(px, py);
                if (pt > 25 && iso < 0.1) {
                    // if (pt > 25) {
                    n_valid = n_valid + 1;
                }
            }
        };

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
                if (pt > 25 && iso < 0.1) {
                    // if (pt > 25) {
                    n_valid = n_valid + 1;
                }
            }
        };

        fill_muon(0, mu1_Px, mu1_Py, mu1_Pz, mu1_E, mu1_Iso, mu1_Pt,
                  NMuon_valid);
        fill_muon(1, mu2_Px, mu2_Py, mu2_Pz, mu2_E, mu2_Iso, mu2_Pt,
                  NMuon_valid);
        fill_muon(2, mu3_Px, mu3_Py, mu3_Pz, mu3_E, mu3_Iso, mu3_Pt,
                  NMuon_valid);
        fill_electron(0, e1_Px, e1_Py, e1_Pz, e1_E, e1_Iso, e1_Pt, Ne_valid);
        fill_electron(1, e2_Px, e2_Py, e2_Pz, e2_E, e2_Iso, e2_Pt, Ne_valid);

        MC_lepton_mass = -1;
        Nlep_valid = NMuon_valid + Ne_valid;
        if (NMuon_valid != 1)
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
        const float BTAG_CUT = 1.;
        std::vector<int> valid_jets_idx;
        std::vector<int> valid_b_jets_idx;

        // std::cout << "---------------------" << std::endl;
        for (int jet = 0; jet < NJet; ++jet) {
            if (!Jet_ID[jet])
                continue;

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
        // if (N_valid_b_jets < 2 || N_valid_jets < 2) {
        //// tree_output->Fill();
        // continue;
        //}
        if (N_valid_jets_tot < 4) {
            //tree_output->Fill();
            continue;
        }

        // if (N_valid_b_jets == 0) {
        //// tree_output->Fill();
        // continue;
        //}
        // std::cout <<"I guess we are here now " << N_valid_jets <<std::endl;
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
        int jet_1, jet_2;
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
                    jet_1 = j1;
                    jet_2 = j2;
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

        float nu_pz =
            CalculateNeutrinoPz(mu1_E, mu1_Px, mu1_Py, mu1_Pz, MET_px, MET_py);

        ROOT::Math::PxPyPzMVector nu_vector(MET_px, MET_py, nu_pz, 0);
        W_leptonic_mass = (mu_vector + nu_vector).M();
        top_leptoninc_mass_1 = (mu_vector + nu_vector + b1_vector).M();
        top_leptoninc_mass_2 = (mu_vector + nu_vector + b1_vector).M();
        // for (size_t idx1 = 0; idx1 < valid_b_jets_idx.size(); ++idx1) {
        // int bj1 = valid_b_jets_idx[idx1];
        // ROOT::Math::PxPyPzEVector b_vector(Jet_Px[bj1], Jet_Py[bj1],
        // Jet_Pz[bj1], Jet_E[bj1]);

        // float top_vector_mass = (W_hadronic_vector + b_vector).M();
        // float mass_diff = std::abs(top_vector_mass - TOP_MASS);
        // if (mass_diff < hadronic_top_mass_diff) {
        // hadronic_top_mass_diff = mass_diff;
        // top_hadronic_mass = top_vector_mass;
        //}
        //}

        tree_output->Fill();
    }

    std::cout << "Output Tuple filled" << std::endl;
}

double CreateTuple::CalculateNeutrinoPz(double mu_E, double mu_px, double mu_py,
                                        double mu_pz, double nu_px,
                                        double nu_py) {

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
    // Scenario B: Discriminant is Negative (Complex roots, mismeasurement)
    else {
        // Set D = 0 and just take the real part
        pz_nu = -B / (2.0 * A);
    }

    return pz_nu;
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
// for (size_t idx1 = 0; idx1 < 2; ++idx1) {
#endif // if LIB_CreateTuple_H
