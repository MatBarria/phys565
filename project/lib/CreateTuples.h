#ifndef LIB_CreateTuple_H
#define LIB_CreateTuple_H

#include "kinematics.h"
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TTree.h>
// #include <algorithm>
#include <cmath>
// #include <cstddef>
// #include <iostream>
// #include <math.h>
// #include <ostream>

class CreateTuple {
  public:
    CreateTuple(TString input, TString output, TString channel);

    virtual ~CreateTuple();

    void fillChain();

    void setBranchesAddressesInput();

    void setBranchesAddressesOutput();

    void fillOutputTree(TString channel);

    void fillFitArrays(const std::array<int, 4> &jet_perm,

                       // double mu1_Px, double mu1_Py, double mu1_Pz,
                       // double MET_px, double MET_py,

                       std::array<double, N_PAR> &p,
                       std::array<double, N_MEAS> &x_meas,
                       std::array<double, N_MEAS> &sigma);
    /**
     * @brief Saves the output tree to a ROOT file in the specified output
     * directory.
     */
    void saveTree();

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
    float nu_pz;
    int NMuon_valid, Ne_valid, Nlep_valid;
    int NMuon_valid_mc, Ne_valid_mc, Nlep_valid_mc;
    float jet_Px[6] = {-1}, jet_Py[6] = {-1}, jet_Pz[6] = {-1}, jet_E[6] = {0},
          jet_Pt[6] = {-1};
    float bjet_Px[4] = {-1}, bjet_Py[4] = {-1}, bjet_Pz[4] = {-1},
          bjet_E[4] = {-1}, bjet_Pt[4] = {-1}, bjet_Btag[4] = {-1};
    int N_valid_jets, N_valid_b_jets, N_valid_jets_tot;

    float MET_pt;

    float weight, chi2, permutation_weight;
    float MCtop_mass_hadronic, MCtop_mass_leptonic;
    float W_leptonic_mass, W_hadronic_mass, W_leptonic_mass_reco,
        W_hadronic_mass_reco, top_hadronic_mass_1, top_hadronic_mass_2,
        top_leptoninc_mass_1, top_leptoninc_mass_2, top_hadronic_mass,
        top_leptoninc_mass;
};

#endif // if LIB_CreateTuple_H
