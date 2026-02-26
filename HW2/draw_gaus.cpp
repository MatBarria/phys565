#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
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
  // tdrStyle->SetFrameBorderMode(0);
  // tdrStyle->SetFrameFillColor(kWhite);

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
  tdrStyle->SetPadRightMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.12);
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

void draw_gaus() {

  TCanvas *c1 = new TCanvas("", "", 800, 600);

  TH1D *h = new TH1D("", "", 100, -5, 5);

  TRandom3 rng(0);

  double mean = 0.0;
  double sigma = 1.0;

  int N = 20000;
  for (int i = 0; i < N; i++) {
    double value = rng.Gaus(mean, sigma);
    h->Fill(value);
  }

  h->SetLineColor(kBlack);
  h->Draw("E");

  TF1 *fit = new TF1("fit", "gaus", -5, 5);
  fit->SetLineColor(kRed);
  h->Fit(fit, "R");

  // Save output
  c1->SaveAs("Gaussian_fit.pdf");
  c1->SaveAs("Gaussian_fit.png");
}
