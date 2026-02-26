#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include <iostream>

static TStyle *tdrStyle = 0;

void setTDRStyle() {
  // delete old style if present
  if (tdrStyle)
    delete tdrStyle;
  tdrStyle = new TStyle("tdrStyle", "CMS TDR style");

  // use plain white canvas
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetPadBorderMode(0);
  tdrStyle->SetPadColor(kWhite);

  // white frame
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameFillColor(kWhite);

  // no statistics box
  tdrStyle->SetOptStat(0);
  tdrStyle->SetOptTitle(0);

  // tick marks on left & right
  tdrStyle->SetPadTickX(1);
  tdrStyle->SetPadTickY(1);

  // axis labels & titles
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  tdrStyle->SetTitleOffset(1.2, "Y");

  // margins
  tdrStyle->SetPadTopMargin(0.07);
  tdrStyle->SetPadRightMargin(0.19);
  tdrStyle->SetPadBottomMargin(0.16);
  tdrStyle->SetPadLeftMargin(0.16);

  // lines, markers
  tdrStyle->SetLineWidth(2);
  tdrStyle->SetHistLineWidth(2);
  tdrStyle->SetMarkerSize(1.1);

  // apply
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
}

// ---- auto-run on load ----
__attribute__((constructor)) static void _auto_apply_style() { setTDRStyle(); }

void draw_dalitz() {

  const float PI_MASS = .13957;
  const float P_MASS = .93827; // proton mass
  const float N_MASS = .93956; // neutron mass
  const float E_BEAM = .650;   // neutron mass

  TLorentzVector target(0.0, 0.0, 0.0, P_MASS);
  TLorentzVector beam(0.0, 0.0,
                      TMath::Sqrt(E_BEAM * E_BEAM - PI_MASS * PI_MASS), E_BEAM);
  TLorentzVector W =
      beam +
      target; // Boost from CM(W-rest) -> LAB is +beta, from LAB -> CM is -beta

  TVector3 boost_vector = W.BoostVector(); // CM->lab

  //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[3] = {N_MASS, PI_MASS, PI_MASS};

  TLorentzVector beam_CM = beam;

  // Cheking if the boost make sense
  // std::cout << "P CM Before = ("
  //<< beam_CM.X() << ", "
  //<< beam_CM.Y() << ", "
  //<< beam_CM.Z() << ")  "
  //<< "E = " << beam_CM.T()
  //<< "  M = " << beam_CM.M()
  //<< std::endl;

  beam_CM.Boost(-boost_vector); // minus beacuse we want lab-> CM
  // TVector3 beam_CM_norm = beam_CM.Vect().Unit();
  // beam_CM_norm.Print();

  // Cheking if the boost make sense
  // std::cout << "P CM After = ("
  //<< beam_CM.X() << ", "
  //<< beam_CM.Y() << ", "
  //<< beam_CM.Z() << ")  "
  //<< "E = " << beam_CM.T()
  //<< "  M = " << beam_CM.M()
  //<< std::endl;

  TGenPhaseSpace event;
  event.SetDecay(W, 3, masses);

  TH1F *h_theta_lab = new TH1F("", "", 90, 0, 180);

  TH1F *h_theta_cm = new TH1F("", "", 90, 0, 180);
  TH2F *dalitz = new TH2F("", "", 80, -3.1, 3.8, 80, -1.1, -1.8);

  TCanvas *canvas_1 = new TCanvas("", "", 1600, 1400);
  TCanvas *canvas_2 = new TCanvas("", "", 1600, 1400);
  TCanvas *canvas_3 = new TCanvas("", "", 1600, 1400);

  const int TOTAL_EVENTS = 12000;
  for (int n = 0; n < TOTAL_EVENTS; n++) {

    Double_t weight = event.Generate();

    // Dalitz
    TLorentzVector *p_neutron = event.GetDecay(0);

    TLorentzVector *p_pi_plus = event.GetDecay(1);
    TLorentzVector *p_pi_minus = event.GetDecay(2);

    TLorentzVector p_n_pi_plus = *p_neutron + *p_pi_plus;
    TLorentzVector p_n_pi_minus = *p_neutron + *p_pi_minus;

    dalitz->Fill(p_n_pi_minus.M2(), p_n_pi_plus.M2(), weight);

    // CM angle
    float theta_lab =
        p_neutron->Vect().Angle(beam.Vect()) * 180.0 / TMath::Pi();
    h_theta_lab->Fill(theta_lab, weight);

    // LAB angle: boost neutron from lab -> cm
    TLorentzVector p_neutron_cm = *p_neutron;
    p_neutron_cm.Boost(-boost_vector);
    float theta_cm =
        p_neutron_cm.Vect().Angle(beam_CM.Vect()) * 180.0 / TMath::Pi();
    h_theta_cm->Fill(theta_cm, weight);
  }

  dalitz->GetZaxis()->SetTitle("Events/bin");
  dalitz->GetXaxis()->SetTitle("$m^2_{n\\pi^-}$");
  dalitz->GetYaxis()->SetTitle("$m^2_{n\\pi^+}$");

  h_theta_cm->GetYaxis()->SetTitle("Events/bin");
  h_theta_cm->GetXaxis()->SetTitle("$\\theta_{cm} $");

  h_theta_lab->GetYaxis()->SetTitle("Events/bin");
  h_theta_lab->GetXaxis()->SetTitle("$\\theta_{lab} $");

  canvas_1->cd();
  dalitz->Draw("COLZ");
  canvas_1->SaveAs("dalitz.pdf");
  canvas_1->SaveAs("dalitz.png");
  delete canvas_1;
  canvas_2->cd();
  h_theta_cm->Draw("");
  canvas_2->SaveAs("theta_cm.pdf");
  canvas_2->SaveAs("theta_cm.png");
  delete canvas_2;

  canvas_3->cd();
  h_theta_lab->Draw("");
  canvas_3->SaveAs("theta_lab.pdf");
  canvas_3->SaveAs("theta_lab.png");
  delete canvas_3;
}
