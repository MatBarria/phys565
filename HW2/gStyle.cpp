// tdrstyle.C
// CMS TDR style for ROOT plots

#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"

static TStyle *tdrStyle = 0;

void setTDRStyle() {
  // delete old style if present
  if (tdrStyle) delete tdrStyle;
  tdrStyle = new TStyle("tdrStyle","CMS TDR style");

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
  tdrStyle->SetTitleFont(42,"XYZ");
  tdrStyle->SetLabelFont(42,"XYZ");
  tdrStyle->SetLabelSize(0.05,"XYZ");
  tdrStyle->SetTitleSize(0.06,"XYZ");
  tdrStyle->SetTitleOffset(1.2,"Y");

  // margins
  tdrStyle->SetPadTopMargin(0.07);
  tdrStyle->SetPadRightMargin(0.05);
  tdrStyle->SetPadBottomMargin(0.12);
  tdrStyle->SetPadLeftMargin(0.16);

  // lines, markers
  tdrStyle->SetLineWidth(2);
  tdrStyle->SetHistLineWidth(2);
  tdrStyle->SetMarkerSize(1.1);

  // colors
  tdrStyle->SetPalette(1);

  // apply
  gROOT->SetStyle("tdrStyle");
  gROOT->ForceStyle();
}
