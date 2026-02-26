#include "TColor.h"
#include "TROOT.h"
#include "TStyle.h"
#include <iostream>

static TStyle *tdrStyle = 0;

void setTDRStyle() {
  // delete old style if present
  if (tdrStyle)
    delete tdrStyle;
  tdrStyle = new TStyle("tdrStyle", "TDR style");

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

__attribute__((constructor)) static void _auto_apply_style() { setTDRStyle(); }
void draw_sin() {

  TCanvas *canvas_1 = new TCanvas("", "", 1175, 1000);
  TF2 *f2 = new TF2("f2", "sin(x)*sin(y)/(x*y)", 0, 5, 0, 5);
  f2->GetXaxis()->SetTitle("x");
  f2->GetYaxis()->SetTitle("y");

  f2->SetNpx(400);
  f2->SetNpy(400);

  f2->Draw("COLZ");
  canvas_1->SaveAs("Trig_func.pdf");
  canvas_1->SaveAs("Trig_func.png");
}
