from ROOT import TASImage, TGaxis, TLatex, TPad, TStyle, gPad, kBlack, kWhite

print("Updated style  4")


def tdrGrid(gridOn, tdrStyle=None):
    if tdrStyle is None:
        tdrStyle = setTDRStyle()
    tdrStyle.SetPadGridX(gridOn)
    tdrStyle.SetPadGridY(gridOn)


def fixOverlay():
    gPad.RedrawAxis()


def setTDRStyle():
    tdrStyle = TStyle("tdrStyle", "Style for P-TDR")
    tdrStyle.SetImageScaling(3.);

    tdrStyle.SetCanvasBorderMode(0)
    tdrStyle.SetCanvasColor(kWhite)
    tdrStyle.SetCanvasDefH(600)
    tdrStyle.SetCanvasDefW(600)
    tdrStyle.SetCanvasDefX(0)
    tdrStyle.SetCanvasDefY(0)

    tdrStyle.SetPadBorderMode(0)
    tdrStyle.SetPadColor(kWhite)
    tdrStyle.SetPadGridX(False)
    tdrStyle.SetPadGridY(False)
    tdrStyle.SetGridColor(0)
    tdrStyle.SetGridStyle(3)
    tdrStyle.SetGridWidth(1)

    tdrStyle.SetFrameBorderMode(0)
    tdrStyle.SetFrameBorderSize(1)
    tdrStyle.SetFrameFillColor(0)
    tdrStyle.SetFrameFillStyle(0)
    tdrStyle.SetFrameLineColor(1)
    tdrStyle.SetFrameLineStyle(1)
    tdrStyle.SetFrameLineWidth(1)

    tdrStyle.SetHistLineColor(2)
    tdrStyle.SetHistLineStyle(1)
    tdrStyle.SetHistLineWidth(1)

    tdrStyle.SetEndErrorSize(2)
    tdrStyle.SetErrorX(0)
    tdrStyle.SetMarkerStyle(20)

    tdrStyle.SetOptFit(1)
    tdrStyle.SetFitFormat("5.4g")
    tdrStyle.SetFuncColor(2)
    tdrStyle.SetFuncStyle(1)
    tdrStyle.SetFuncWidth(1)

    tdrStyle.SetOptDate(0)

    tdrStyle.SetOptFile(0)
    tdrStyle.SetOptStat(0)
    tdrStyle.SetStatColor(kWhite)
    tdrStyle.SetStatFont(42)
    tdrStyle.SetStatFontSize(0.025)
    tdrStyle.SetStatTextColor(1)
    tdrStyle.SetStatFormat("6.4g")
    tdrStyle.SetStatBorderSize(1)
    tdrStyle.SetStatH(0.1)
    tdrStyle.SetStatW(0.15)

    tdrStyle.SetPadTopMargin(0.05)
    tdrStyle.SetPadBottomMargin(0.13)
    tdrStyle.SetPadLeftMargin(0.16)
    tdrStyle.SetPadRightMargin(0.02)

    # tdrStyle.SetOptTitle(1)
    # tdrStyle.SetTitleFont(42)
    # tdrStyle.SetTitleColor(1)
    # tdrStyle.SetTitleTextColor(1)
    # tdrStyle.SetTitleFillColor(10)
    # tdrStyle.SetTitleFontSize(0.05)

    # tdrStyle.SetTitleColor(1, 'XYZ')
    # tdrStyle.SetTitleFont(42, 'XYZ')
    # tdrStyle.SetTitleSize(0.06, 'XYZ')
    # tdrStyle.SetTitleXOffset(0.9)
    # tdrStyle.SetTitleYOffset(1.25)
    tdrStyle.SetOptTitle(1)
    tdrStyle.SetTitleBorderSize(0)
    tdrStyle.SetTitleFillColor(0)
    tdrStyle.SetTitleStyle(0)
    tdrStyle.SetTitleAlign(23)  # centered
    tdrStyle.SetTitleX(0.5)  # center x
    tdrStyle.SetTitleY(0.995)  # near top
    tdrStyle.SetTitleW(0.0)  # let ROOT size it
    tdrStyle.SetTitleH(0.0)
    tdrStyle.SetLabelColor(1, "XYZ")
    tdrStyle.SetLabelFont(42, "XYZ")
    tdrStyle.SetLabelOffset(0.007, "XYZ")
    tdrStyle.SetLabelSize(0.05, "XYZ")

    tdrStyle.SetAxisColor(1, "XYZ")
    tdrStyle.SetStripDecimals(True)
    tdrStyle.SetTickLength(0.03, "XYZ")
    tdrStyle.SetNdivisions(510, "XYZ")
    tdrStyle.SetPadTickX(1)
    tdrStyle.SetPadTickY(1)

    tdrStyle.SetOptLogx(0)
    tdrStyle.SetOptLogy(0)
    tdrStyle.SetOptLogz(0)

    tdrStyle.SetPaperSize(20.0, 20.0)
    tdrStyle.SetHatchesLineWidth(5)
    tdrStyle.SetHatchesSpacing(0.05)

    tdrStyle.cd()
    return tdrStyle
