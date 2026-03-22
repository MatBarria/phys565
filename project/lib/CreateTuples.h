#ifndef LIB_CreateTuple_H
#define LIB_CreateTuple_H

#include "Constants.h"
#include "Run3Constants.h"
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <TString.h>
#include <TTree.h>
#include <cmath>
#include <iostream>
#include <math.h>

/**
 * @class CreateTuple
 * @brief This Class read outputs from HmmAnalyzer and return a smaller
 * tree_input with the variables necessary for the categorization
 */
class CreateTuple {
public:
  /**
   * @brief Constructor that initializes the CreateTuple class with specified
   * input and output settings.
   *
   * @param input The full path of the input ROOT file containing the data.
   * @param output Directory path for the output ROOT file.
   */
  CreateTuple(TString input, TString output, TString era, TString channel,
              bool is_data, bool is_signal_);
  /**
   * @brief Destructor for CreateTuple, cleaning up dynamically allocated
   * objects.
   */
  virtual ~CreateTuple();

  /**
   * @brief Loads the input file and adds its data to the input TChain.
   *
   * @throws std::runtime_error if the input file is not found.
   */
  void fillChain();

  /**
   * @brief Sets up branch addresses for reading variables from the input
   * tree.
   */
  void setBranchesAddressesInput();

  /**
   * @brief Sets up branch addresses for writing variables to the output tree.
   */
  void setBranchesAddressesOutput();

  /**
   * @brief Processes events and fills the output tree with calculated
   * variables.
   */
  void fillOutputTree();

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

  /** Read event variables */
  float NPrimaryVertices, triggerIsoMu24, EventWeight;

  /** Read MET variables **/
  float MET_Px, MET_Py;

  /** Read Muon variables */
  int NMuon;
  std::vector<float> *Muon_Px, *Muon_Py, *Muon_Pz, *Muon_E, *Muon_Charge,
      *Muon_Iso;

  /** Read electron variables */
  int NElectron;
  std::vector<float> *Electron_Px, *Electron_Py, *Electron_Pz, *Electron_E,
      *Electron_Charge, *Electron_Iso;

  /** Read Jet variables */
  int NJet;
  std::vector<float> *NJet, *Jet_Px, *Jet_Py, *Jet_Pz, *Jet_E, *Jet_btag,
      *Jet_ID;

  /** Read Photon variables */
  int NPhoton;
  std::vector<float> *Photon_Px, *Photon_Py, *Photon_Pz, *Photon_E, *Photon_Iso;

  /** Read MC variables */
  float *MChadronicBottom_px, *MChadronicBottom_py, *MChadronicBottom_pz,
      *MCleptonicBottom_px, *MCleptonicBottom_py, *MCleptonicBottom_pz,
      *MChadronicWDecayQuark_px, *MChadronicWDecayQuark_py,
      *MChadronicWDecayQuark_pz, *MChadronicWDecayQuarkBar_px,
      *MChadronicWDecayQuarkBar_py, *MChadronicWDecayQuarkBar_pz, *MClepton_px,
      *MClepton_py, *MClepton_pz, *MCleptonPDGid, *MCneutrino_px,
      *MCneutrino_py, *MCneutrino_pz;

  /** New Event variables */
  float diMuon_mass;
};

CreateTuple::CreateTuple(TString input, TString output) {

  tree_input = new TChain("events");
  tree_output = new TTree("tree_output", "tree_output");
  input_name = input;
  output_directory = output;
  output_name = channel + "tuples.root";
  Muon_Px = nullptr;
  Muon_Py = nullptr;
  Muon_Pz = nullptr;
  Muon_E = nullptr;
  Muon_Charge = nullptr;
  Muon_Iso = nullptr;

  Electron_Px = nullptr;
  Electron_Py = nullptr;
  Electron_Pz = nullptr;
  Electron_E = nullptr;
  Electron_Iso = nullptr;

  Jet_Px = nullptr;
  Jet_Py = nullptr;
  Jet_Pz = nullptr;
  Jet_E = nullptr;
  Jet_btag = nullptr;
  Jet_ID = nullptr;

  Photon_Px = nullptr;
  Photon_Py = nullptr;
  Photon_Pz = nullptr;
  Photon_E = nullptr;
  Photon_Iso = nullptr;

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

  tree_output->Branch("gen_weight", &gen_weight, "gen_weight/f");
  tree_output->Branch("pileup_weight", &pileup_weight, "pileup_weight/f");
  tree_output->Branch("pileup_weight_up", &pileup_weight_up,
                      "pileup_weight_up/f");
  tree_output->Branch("pileup_weight_down", &pileup_weight_down,
                      "pileup_weight_down/f");
  tree_output->Branch("scale_factor", &scale_factor, "scale_factor/d");
  tree_output->Branch("weight", &weight, "weight/D");
  tree_output->Branch("weight_no_lumi", &weight_no_lumi, "weight_no_lumi/d");
  // tree_output->Branch("is_data", &is_data_int, "is_data/i");
  // tree_output->Branch("is_signal", &is_signal_int, "signal/i");
  tree_output->Branch("rho", &rho, "rho/f");
  tree_output->Branch("PV", &pv, "PV/i");
  tree_output->Branch("n_SoftJet_pt2", &n_SoftJet_pt2, "n_SoftJet_pt2/i");
  tree_output->Branch("n_SoftJet_pt5", &n_SoftJet_pt5, "n_SoftJet_pt5/i");
  tree_output->Branch("n_SoftJet_pt10", &n_SoftJet_pt10, "n_SoftJet_pt10/i");
  tree_output->Branch("HT", &HT, "HT/f");
  tree_output->Branch("HT_pt2", &HT_pt2, "HT_pt2/f");
  tree_output->Branch("HT_pt5", &HT_pt5, "HT_pt5/f");
  tree_output->Branch("HT_pt10", &HT_pt10, "HT_pt10/f");

  tree_output->Branch("is_ggH_category", &is_ggH_category, "is_ggH_category/i");
  tree_output->Branch("is_VBF_category", &is_VBF_category, "is_VBF_category/i");

  // DiMuon variables
  tree_output->Branch("diMuon_mass", &diMuon_mass, "diMuon_mass/f");
  tree_output->Branch("diMuon_bsConstrainedMass", &diMuon_bsConstrainedMass,
                      "diMuon_bsConstrainedMass/f");
  tree_output->Branch("diMuon_pt", &diMuon_pt, "diMuon_pt/f");
  tree_output->Branch("diMuon_bsConstrainedPt", &diMuon_bsConstrainedPt,
                      "diMuon_bsConstrainedPt/f");
  tree_output->Branch("diMuon_phi", &diMuon_phi, "diMuon_phi/f");
  tree_output->Branch("diMuon_eta", &diMuon_eta, "diMuon_eta/f");
  tree_output->Branch("diMuon_rapidity", &diMuon_rapidity, "diMuon_rapidity/f");

  // MET variables
  tree_output->Branch("MET_phi", &MET_phi, "MET_phi/f");
  tree_output->Branch("MET_pt", &MET_pt, "MET_pt/f");
  tree_output->Branch("MET_sumEt", &MET_sumEt, "MET_sumEt/f");

  tree_output->Branch("ChsMET_phi", &ChsMET_phi, "ChsMET_phi/f");
  tree_output->Branch("ChsMET_pt", &ChsMET_pt, "ChsMET_pt/f");
  tree_output->Branch("ChsMET_sumEt", &ChsMET_sumEt, "ChsMET_sumEt/f");

  tree_output->Branch("PuppiMET_phi", &PuppiMET_phi, "PuppiMET_phi/f");
  tree_output->Branch("PuppiMET_pt", &PuppiMET_pt, "PuppiMET_pt/f");
  tree_output->Branch("PuppiMET_sumEt", &PuppiMET_sumEt, "PuppiMET_sumEt/f");

  // Muon variables
  tree_output->Branch("mu1_pt_mass_ratio", &mu1_pt_mass_ratio,
                      "mu1_pt_mass_ratio/f");
  tree_output->Branch("mu2_pt_mass_ratio", &mu2_pt_mass_ratio,
                      "mu2_pt_mass_ratio/f");
  tree_output->Branch("mu1_bsConstrainedPt_mass_ratio",
                      &mu1_bsConstrainedPt_mass_ratio,
                      "mu1_bsConstrainedPt_mass_ratio/f");
  tree_output->Branch("mu2_bsConstrainedPt_mass_ratio",
                      &mu2_bsConstrainedPt_mass_ratio,
                      "mu2_bsConstrainedPt_mass_ratio/f");
  tree_output->Branch("mu1_eta", &mu1_eta, "mu1_eta/f");
  tree_output->Branch("mu2_eta", &mu2_eta, "mu2_eta/f");
  tree_output->Branch("mu1_pt", &mu1_pt, "mu1_pt/f");
  tree_output->Branch("mu2_pt", &mu2_pt, "mu2_pt/f");
  tree_output->Branch("mu1_bsConstrainedPt", &mu1_bsConstrainedPt,
                      "mu1_bsConstrainedPt/f");
  tree_output->Branch("mu2_bsConstrainedPt", &mu2_bsConstrainedPt,
                      "mu2_bsConstrainedPt/f");
  tree_output->Branch("mu1_ptErr", &mu1_ptErr, "mu1_ptErr/f");
  tree_output->Branch("mu2_ptErr", &mu2_ptErr, "mu2_ptErr/f");
  tree_output->Branch("mu1_bsConstrainedPtErr", &mu1_bsConstrainedPtErr,
                      "mu1_bsConstrainedPtErr/f");
  tree_output->Branch("mu2_bsConstrainedPtErr", &mu2_bsConstrainedPtErr,
                      "mu2_bsConstrainedPtErr/f");
  tree_output->Branch("phi_CS", &phi_CS, "phi_CS/f");
  tree_output->Branch("cos_theta_CS", &cos_theta_CS, "cos_theta_CS/f");

  // Error variables
  tree_output->Branch("relative_diMuon_mass_error", &relative_diMuon_mass_error,
                      "relative_diMuon_mass_error/f");
  tree_output->Branch("relative_diMuon_bsConstrainedMass_error",
                      &relative_diMuon_bsConstrainedMass_error,
                      "relative_diMuon_bsConstrainedMass_error/f");

  // Jet variables
  tree_output->Branch("n_jet", &n_jet, "n_jet/i");
  tree_output->Branch("n_bjet", &n_bjet, "n_bjet/i");
  tree_output->Branch("n_bjet_Loose", &n_bjet_Loose, "n_bjet_Loose/i");
  tree_output->Branch("leading_jet_pt", &leading_jet_pt, "leading_jet_pt/f");
  tree_output->Branch("leading_jet_eta", &leading_jet_eta, "leading_jet_eta/f");
  tree_output->Branch("subleading_jet_pt", &subleading_jet_pt,
                      "subleading_jet_pt/f");
  tree_output->Branch("subleading_jet_eta", &subleading_jet_eta,
                      "subleading_jet_eta/f");
  //
  // diJet variables
  tree_output->Branch("diJet_mass", &diJet_mass, "diJet_mass/f");
  tree_output->Branch("delta_eta_diJet", &delta_eta_diJet, "delta_eta_diJet/f");
  tree_output->Branch("delta_phi_diJet", &delta_phi_diJet, "delta_phi_diJet/f");
  tree_output->Branch("z_zeppenfeld", &z_zeppenfeld, "z_zeppenfeld/f");
  tree_output->Branch("pt_balance", &pt_balance, "pt_balance/f");
  tree_output->Branch("pt_centrality", &pt_centrality, "pt_centrality/f");
  tree_output->Branch("min_delta_eta_diMuon_jet", &min_delta_eta_diMuon_jet,
                      "min_delta_eta_diMuon_jet/f");
  tree_output->Branch("min_delta_phi_diMuon_jet", &min_delta_phi_diMuon_jet,
                      "min_delta_phi_diMuon_jet/f");
  // DiJet variables
  std::cout << "Output Branches addressed setted" << std::endl;
}

void CreateTuple::setBranchesAddressesInput() {

  tree_input->SetBranchStatus("*", 0);

  tree_input->SetBranchStatus("t_puWeight", 1);
  tree_input->SetBranchAddress("t_puWeight", &pileup_weight);
  tree_input->SetBranchStatus("t_puWeightUp", 1);
  tree_input->SetBranchAddress("t_puWeightUp", &pileup_weight_up);
  tree_input->SetBranchStatus("t_puWeightDown", 1);
  tree_input->SetBranchAddress("t_puWeightDown", &pileup_weight_down);
  tree_input->SetBranchStatus("t_genWeight", 1);
  tree_input->SetBranchAddress("t_genWeight", &gen_weight);
  tree_input->SetBranchStatus("t_Rho", 1);
  tree_input->SetBranchAddress("t_Rho", &rho);
  tree_input->SetBranchStatus("t_PV_npvsGood", 1);
  tree_input->SetBranchAddress("t_PV_npvsGood", &pv);
  tree_input->SetBranchStatus("t_SoftActivityJetNjets2", 1);
  tree_input->SetBranchAddress("t_SoftActivityJetNjets2", &n_SoftJet_pt2);
  tree_input->SetBranchStatus("t_SoftActivityJetNjets5", 1);
  tree_input->SetBranchAddress("t_SoftActivityJetNjets5", &n_SoftJet_pt5);
  tree_input->SetBranchStatus("t_SoftActivityJetNjets10", 1);
  tree_input->SetBranchAddress("t_SoftActivityJetNjets10", &n_SoftJet_pt10);

  tree_input->SetBranchStatus("t_SoftActivityJetHT", 1);
  tree_input->SetBranchAddress("t_SoftActivityJetHT", &HT);
  tree_input->SetBranchStatus("t_SoftActivityJetHT2", 1);
  tree_input->SetBranchAddress("t_SoftActivityJetHT2", &HT_pt2);
  tree_input->SetBranchStatus("t_SoftActivityJetHT5", 1);
  tree_input->SetBranchAddress("t_SoftActivityJetHT5", &HT_pt5);
  tree_input->SetBranchStatus("t_SoftActivityJetHT10", 1);
  tree_input->SetBranchAddress("t_SoftActivityJetHT10", &HT_pt10);

  tree_input->SetBranchStatus("t_genWeight", 1);
  tree_input->SetBranchAddress("t_genWeight", &gen_weight);
  // DiMuon variables
  tree_input->SetBranchStatus("t_diMuon_bsConstrainedMass", 1);
  tree_input->SetBranchAddress("t_diMuon_bsConstrainedMass",
                               &diMuon_bsConstrainedMass);
  tree_input->SetBranchStatus("t_diMuon_bsConstrainedPt", 1);
  tree_input->SetBranchAddress("t_diMuon_bsConstrainedPt",
                               &diMuon_bsConstrainedPt);
  tree_input->SetBranchStatus("t_diMuon_mass", 1);
  tree_input->SetBranchAddress("t_diMuon_mass", &diMuon_mass);
  tree_input->SetBranchStatus("t_diMuon_pt", 1);
  tree_input->SetBranchAddress("t_diMuon_pt", &diMuon_pt);
  tree_input->SetBranchStatus("t_diMuon_phi", 1);
  tree_input->SetBranchAddress("t_diMuon_phi", &diMuon_phi);
  tree_input->SetBranchStatus("t_diMuon_eta", 1);
  tree_input->SetBranchAddress("t_diMuon_eta", &diMuon_eta);

  // MET varaibles
  tree_input->SetBranchStatus("t_MET_phi", 1);
  tree_input->SetBranchAddress("t_MET_phi", &MET_phi);
  tree_input->SetBranchStatus("t_MET_pt", 1);
  tree_input->SetBranchAddress("t_MET_pt", &MET_pt);
  tree_input->SetBranchStatus("t_MET_sumEt", 1);
  tree_input->SetBranchAddress("t_MET_sumEt", &MET_sumEt);

  tree_input->SetBranchStatus("t_ChsMET_phi", 1);
  tree_input->SetBranchAddress("t_ChsMET_phi", &ChsMET_phi);
  tree_input->SetBranchStatus("t_ChsMET_pt", 1);
  tree_input->SetBranchAddress("t_ChsMET_pt", &ChsMET_pt);
  tree_input->SetBranchStatus("t_ChsMET_sumEt", 1);
  tree_input->SetBranchAddress("t_ChsMET_sumEt", &ChsMET_sumEt);

  tree_input->SetBranchStatus("t_PuppiMET_phi", 1);
  tree_input->SetBranchAddress("t_PuppiMET_phi", &PuppiMET_phi);
  tree_input->SetBranchStatus("t_PuppiMET_pt", 1);
  tree_input->SetBranchAddress("t_PuppiMET_pt", &PuppiMET_pt);
  tree_input->SetBranchStatus("t_PuppiMET_sumEt", 1);
  tree_input->SetBranchAddress("t_PuppiMET_sumEt", &PuppiMET_sumEt);

  // Muon variables
  tree_input->SetBranchStatus("t_mu1", 1);
  tree_input->SetBranchAddress("t_mu1", &mu1_index);
  tree_input->SetBranchStatus("t_mu2", 1);
  tree_input->SetBranchAddress("t_mu2", &mu2_index);
  tree_input->SetBranchStatus("t_Mu_charge", 1);
  tree_input->SetBranchAddress("t_Mu_charge", &mu_charge);
  tree_input->SetBranchStatus("t_Mu_pt", 1);
  tree_input->SetBranchAddress("t_Mu_pt", &mu_pt);
  tree_input->SetBranchStatus("t_Mu_bsConstrainedPt", 1);
  tree_input->SetBranchAddress("t_Mu_bsConstrainedPt", &mu_bsConstrainedPt);
  tree_input->SetBranchStatus("t_Mu_ptErr", 1);
  tree_input->SetBranchAddress("t_Mu_ptErr", &mu_ptErr);
  tree_input->SetBranchStatus("t_Mu_bsConstrainedPtErr", 1);
  tree_input->SetBranchAddress("t_Mu_bsConstrainedPtErr",
                               &mu_bsConstrainedPtErr);
  tree_input->SetBranchStatus("t_Mu_phi", 1);
  tree_input->SetBranchAddress("t_Mu_phi", &mu_phi);
  tree_input->SetBranchStatus("t_Mu_eta", 1);
  tree_input->SetBranchAddress("t_Mu_eta", &mu_eta);

  // Electron variables
  tree_input->SetBranchStatus("t_El_pt", 1);
  tree_input->SetBranchAddress("t_El_pt", &elec_pt);

  // Jet variables
  tree_input->SetBranchStatus("t_nbJet", 1);
  tree_input->SetBranchAddress("t_nbJet", &n_bjet);
  tree_input->SetBranchStatus("t_nbJet_Loose", 1);
  tree_input->SetBranchAddress("t_nbJet_Loose", &n_bjet_Loose);
  tree_input->SetBranchStatus("t_nJet", 1);
  tree_input->SetBranchAddress("t_nJet", &n_jet);
  tree_input->SetBranchStatus("t_Jet_mass", 1);
  tree_input->SetBranchAddress("t_Jet_mass", &jet_mass);
  tree_input->SetBranchStatus("t_Jet_pt", 1);
  tree_input->SetBranchAddress("t_Jet_pt", &jet_pt);
  tree_input->SetBranchStatus("t_Jet_phi", 1);
  tree_input->SetBranchAddress("t_Jet_phi", &jet_phi);
  tree_input->SetBranchStatus("t_Jet_eta", 1);
  tree_input->SetBranchAddress("t_Jet_eta", &jet_eta);

  // DiJet variables
  tree_input->SetBranchStatus("t_diJet_mass", 1);
  tree_input->SetBranchAddress("t_diJet_mass", &diJet_mass);
  tree_input->SetBranchStatus("t_diJet_mass_mo", 1);
  tree_input->SetBranchAddress("t_diJet_mass_mo", &diJet_mass_mo);
  tree_input->SetBranchStatus("t_diJet_pt", 1);
  tree_input->SetBranchAddress("t_diJet_pt", &diJet_pt);
  tree_input->SetBranchStatus("t_diJet_phi", 1);
  tree_input->SetBranchAddress("t_diJet_phi", &diJet_phi);
  tree_input->SetBranchStatus("t_diJet_eta", 1);
  tree_input->SetBranchAddress("t_diJet_eta", &diJet_eta);

  std::cout << "Input Branches addressed setted" << std::endl;
}

void CreateTuple::fillOutputTree() {

  std::cout << "Filling Output Tree" << std::endl;
  double total_entries = tree_input->GetEntries();
  std::cout << "Entries: " << total_entries << std::endl;
  TLorentzVector mu1_vector;
  TLorentzVector mu2_vector;
  TLorentzVector mu1BSC_vector;
  TLorentzVector mu2BSC_vector;
  std::pair<float, float> angles_CS;
  for (int event_index = 0; event_index < total_entries; event_index++) {
    tree_input->GetEntry(event_index);
    // if (diMuon_mass < 110 || diMuon_mass > 150)
    // if (diMuon_mass < 100 || diMuon_mass > 180)
    if (diMuon_mass < 70 || diMuon_mass > 180)
      continue;

    weight = GetEventWeight(gen_weight, pileup_weight, scale_factor);
    weight_no_lumi = weight / luminosity;

    // DiMuon variables
    diMuon_rapidity = (mu1_vector + mu2_vector).Rapidity();

    // Muon variables
    mu1_vector.SetPtEtaPhiM((*mu_pt)[mu1_index], (*mu_eta)[mu1_index],
                            (*mu_phi)[mu1_index], MUON_MASS);
    mu2_vector.SetPtEtaPhiM((*mu_pt)[mu2_index], (*mu_eta)[mu2_index],
                            (*mu_phi)[mu2_index], MUON_MASS);
    mu1BSC_vector.SetPtEtaPhiM((*mu_bsConstrainedPt)[mu1_index],
                               (*mu_eta)[mu1_index], (*mu_phi)[mu1_index],
                               MUON_MASS);
    mu2BSC_vector.SetPtEtaPhiM((*mu_bsConstrainedPt)[mu2_index],
                               (*mu_eta)[mu2_index], (*mu_phi)[mu2_index],
                               MUON_MASS);

    mu1_pt = (*mu_pt)[mu1_index];
    mu2_pt = (*mu_pt)[mu2_index];
    mu1_bsConstrainedPt = (*mu_bsConstrainedPt)[mu1_index];
    mu2_bsConstrainedPt = (*mu_bsConstrainedPt)[mu2_index];
    mu1_ptErr = (*mu_ptErr)[mu1_index];
    mu2_ptErr = (*mu_ptErr)[mu2_index];
    mu1_bsConstrainedPtErr = (*mu_bsConstrainedPtErr)[mu1_index];
    mu2_bsConstrainedPtErr = (*mu_bsConstrainedPtErr)[mu2_index];
    relative_diMuon_mass_error = std::sqrt(std::pow(mu1_ptErr / mu1_pt, 2) +
                                           std::pow(mu2_ptErr / mu2_pt, 2));
    relative_diMuon_bsConstrainedMass_error =
        std::sqrt(std::pow(mu1_bsConstrainedPtErr / mu1_bsConstrainedPt, 2) +
                  std::pow(mu2_bsConstrainedPtErr / mu2_bsConstrainedPt, 2));

    angles_CS = CSAngles(mu1_vector, mu2_vector, (*mu_charge)[mu1_index]);

    mu1_pt_mass_ratio = (*mu_pt)[mu1_index] / diMuon_mass;
    mu2_pt_mass_ratio = (*mu_pt)[mu2_index] / diMuon_mass;
    mu1_bsConstrainedPt_mass_ratio =
        (*mu_bsConstrainedPt)[mu1_index] / diMuon_bsConstrainedMass;
    mu2_bsConstrainedPt_mass_ratio =
        (*mu_bsConstrainedPt)[mu2_index] / diMuon_bsConstrainedMass;
    mu1_eta = (*mu_eta)[mu1_index];
    mu2_eta = (*mu_eta)[mu2_index];
    phi_CS = angles_CS.second;
    cos_theta_CS = angles_CS.first;

    if (n_jet == 0) {
      leading_jet_pt = 0;
      leading_jet_eta = 0;
      diJet_mass = 0;
      subleading_jet_pt = 0;
      subleading_jet_eta = 0;
      delta_eta_diJet = 0;
      delta_phi_diJet = -1;
      z_zeppenfeld = 0;
      pt_balance = -1;
      pt_centrality = -1;
      min_delta_eta_diMuon_jet = 0;
      min_delta_phi_diMuon_jet = 0;
    } else if (n_jet == 1) {
      leading_jet_pt = jet_pt->at(0);
      leading_jet_eta = jet_eta->at(0);
      diJet_mass = 0;
      subleading_jet_pt = 0;
      subleading_jet_eta = 0;
      delta_eta_diJet = 0;
      delta_phi_diJet = -1;
      z_zeppenfeld = 0;
      pt_balance = -1;
      pt_centrality = -1;
      min_delta_eta_diMuon_jet = DeltaEta(diMuon_eta, jet_eta->at(0));
      min_delta_phi_diMuon_jet = DeltaPhi(diMuon_phi, jet_phi->at(0));
    } else {
      leading_jet_pt = jet_pt->at(0);
      leading_jet_eta = jet_eta->at(0);
      subleading_jet_pt = jet_pt->at(1);
      subleading_jet_eta = jet_eta->at(1);
      delta_eta_diJet = DeltaEta(jet_eta->at(0), jet_eta->at(1));
      delta_phi_diJet = DeltaPhi(jet_phi->at(0), jet_phi->at(1));
      z_zeppenfeld = GetZZeppenfeldVariable(diMuon_rapidity, jet_pt, jet_phi,
                                            jet_eta, jet_mass);
      pt_balance = GetPtBalanceVariable(mu1_vector + mu2_vector, jet_pt,
                                        jet_phi, jet_eta, jet_mass);
      pt_centrality = GetPtCentralityVariable(diMuon_pt, jet_pt, jet_phi,
                                              jet_eta, jet_mass);
      min_delta_eta_diMuon_jet =
          TMath::Min(DeltaEta(diMuon_eta, jet_eta->at(0)),
                     DeltaEta(diMuon_eta, jet_eta->at(1)));
      min_delta_phi_diMuon_jet =
          TMath::Min(TMath::Abs(DeltaPhi(diMuon_phi, jet_phi->at(0))),
                     TMath::Abs(DeltaPhi(diMuon_phi, jet_phi->at(1))));
    }

    // Choose category

    tree_output->Fill();
  }

  std::cout << "Output Tuple filled" << std::endl;
}

void CreateTuple::saveTree() {

  output_file = new TFile(output_directory + output_name, "RECREATE");
  output_file->cd();

  tree_output->Write();

  output_file->Close();
  std::cout << "Tuple saved" << std::endl;
}

#endif // if LIB_CreateTuple_H
