// =============================================================================
//  topmass_measurement.C
//
//  Option 2 top mass measurement using a parametric Novosibirsk PDF.
//
//  Strategy:
//    1. Fit MC fitted-mass distribution → extract peak_0, width, tail
//    2. Build parametric PDF: peak(Mt) = Mt + (peak_0 - Mt_gen)
//       The peak shifts 1:1 with Mt; width and tail are fixed from MC.
//    3. Fit data with Mt as the single free parameter → measurement
//    4. Calibration check using shifted pseudo-data → calibration line
//    5. Apply calibration to get final Mt_true
//
//  Usage (in ROOT):
//    .L topmass_measurement.C
//    runAll("mc.root", "data.root")   // two separate files
//    runAll("mc.root", "mc.root")     // same file: MC = pseudo-data
//
// =============================================================================

#include "RooDataSet.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooFormulaVar.h"
#include "RooNovosibirsk.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <vector>

// =============================================================================
// USER CONFIGURATION — change these to match your TTree
// =============================================================================
namespace cfg {
const char *TREE_NAME = "tree_output";         // TTree name
const char *BRANCH_MASS = "top_hadronic_mass"; // fitted top mass branch
const char *BRANCH_WEIGHT = "total_weight";    // event weight branch
//const char *BRANCH_WEIGHT_SUM = "permuation_weight_sum"; // event weight branch
const double MASS_MIN = 50.0;                            // fit range low  [GeV]
const double MASS_MAX = 400.0;                           // fit range high [GeV]
const double MT_GEN = 172.5;    // generated top mass [GeV]
const double MT_FIT_LO = 160.0; // Mt search range low
const double MT_FIT_HI = 185.0; // Mt search range high
} // namespace cfg

// =============================================================================
// Helpers
// =============================================================================

// Load a weighted RooDataSet from a TTree
// Events are weighted by cfg::BRANCH_WEIGHT — chi2 cut already applied upstream
RooDataSet *loadDataset(const char *filename, const char *name, RooRealVar &x) {
    TFile *f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return nullptr;
    }
    TTree *tree = (TTree *)f->Get(cfg::TREE_NAME);
    if (!tree) {
        std::cerr << "ERROR: TTree '" << cfg::TREE_NAME << "' not found in "
                  << filename << "\n";
        return nullptr;
    }

    float mass = 0., weight = 1., weight_sum = 1;
    tree->SetBranchAddress(cfg::BRANCH_MASS, &mass);
    tree->SetBranchAddress(cfg::BRANCH_WEIGHT, &weight);
    //tree->SetBranchAddress(cfg::BRANCH_WEIGHT_SUM, &weight_sum);

    // RooDataSet needs a weight variable declared in the arg set
    RooRealVar w("w", "event weight", 1., -1e15, 1e10);
    RooDataSet *ds =
        new RooDataSet(name, name, RooArgSet(x, w), RooFit::WeightVar(w));

    Long64_t nEntries = tree->GetEntries();
    Long64_t nPass = 0;
    double sumW = 0.;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        // if (mass < cfg::MASS_MIN || mass > cfg::MASS_MAX) continue;
        x.setVal(mass);
        w.setVal(weight);
        ds->add(RooArgSet(x, w),
                weight ); // second arg = weight value
        sumW += weight ;
        ++nPass;
    }

    std::cout << "[loadDataset] " << name << ": " << nPass
              << " events, sum(w) = " << sumW << "\n";

    return ds;
}

// =============================================================================
// STEP 1: Fit MC shape — returns {peak_0, width, tail}
// =============================================================================
struct ShapeParams {
    double peak0, width, tail;
};

ShapeParams fitMCShape(const char *mcFile, RooRealVar &x) {
    std::cout << "\n========== STEP 1: Fit MC shape ==========\n";

    RooDataSet *mcData = loadDataset(mcFile, "mcData", x);
    if (!mcData || mcData->numEntries() == 0) {
        std::cerr << "ERROR: empty MC dataset\n";
        return {cfg::MT_GEN, 15.0, 0.1};
    }

    // Novosibirsk parameters — start near physically motivated values
    RooRealVar peak("peak", "peak", cfg::MT_GEN, cfg::MT_GEN - 20.,
                    cfg::MT_GEN + 20.);
    RooRealVar width("width", "width", 15, 1., 60.);
    // RooRealVar tail("tail", "tail", 0.01, -0.25, 0.25);
    RooRealVar tail("tail", "tail", 0.05, -4.55, 0.35);

    RooNovosibirsk pdf("mcPdf", "MC Novosibirsk", x, peak, width, tail);

    RooFitResult *r =
        pdf.fitTo(*mcData, RooFit::Save(), RooFit::SumW2Error(true),
                  RooFit::PrintLevel(-1), RooFit::Warnings(false));

    ShapeParams sp{peak.getVal(), width.getVal(), tail.getVal()};

    std::cout << "  peak_0 = " << sp.peak0
              << "  (offset from Mt_gen: " << sp.peak0 - cfg::MT_GEN
              << " GeV)\n"
              << "  width  = " << sp.width << " GeV\n"
              << "  tail   = " << sp.tail << "\n";

    // Plot
    TCanvas *c = new TCanvas("cMC", "MC shape fit", 800, 600);
    RooPlot *frame = x.frame(RooFit::Title("MC fitted top mass"));
    mcData->plotOn(frame, RooFit::Name("mcpoints"));
    pdf.plotOn(frame, RooFit::LineColor(kBlue), RooFit::Name("mcfit"));
    // Draw chi2/ndf on plot
    frame->GetXaxis()->SetTitle("Fitted top mass [GeV]");
    frame->GetYaxis()->SetTitle("Events / bin");
    frame->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.DrawLatex(0.15, 0.85, Form("peak = %.2f GeV", sp.peak0));
    latex.DrawLatex(0.15, 0.80, Form("width = %.2f GeV", sp.width));
    latex.DrawLatex(0.15, 0.75, Form("tail = %.3f", sp.tail));
    c->SaveAs("step1_mc_shape_fit.pdf");
    std::cout << "  Saved: step1_mc_shape_fit.pdf\n";

    return sp;
}

// =============================================================================
// STEP 2: Fit data — Mt is the single free parameter
// =============================================================================
double fitData(const char *dataFile, RooRealVar &x, const ShapeParams &sp) {
    std::cout << "\n========== STEP 2: Fit data ==========\n";

    RooDataSet *obsData = loadDataset(dataFile, "obsData", x);
    if (!obsData || obsData->numEntries() == 0) {
        std::cerr << "ERROR: empty data dataset\n";
        return -1.;
    }

    // The only free parameter: the top mass we are measuring
    RooRealVar Mt("Mt", "top quark mass [GeV]", cfg::MT_GEN, cfg::MT_FIT_LO,
                  cfg::MT_FIT_HI);

    // Fixed shape parameters from MC
    RooRealVar width("width", "width", sp.width); // fixed
    RooRealVar tail("tail", "tail", sp.tail);     // fixed

    // peak(Mt) = Mt + (peak_0 - Mt_gen)
    // This encodes the 1:1 shift: if Mt changes by δ, peak changes by δ
    double offsetVal = sp.peak0 - cfg::MT_GEN;
    RooRealVar offset("offset", "peak offset [GeV]", offsetVal); // fixed
    RooFormulaVar peak("peak", "Mt + offset", RooArgList(Mt, offset));

    RooNovosibirsk pdf("dataPdf", "parametric Novosibirsk", x, peak, width,
                       tail);

    // Fit — use Minos for correct asymmetric uncertainties
    // RooFitResult *r = pdf.fitTo(
    //*obsData, RooFit::Save(), RooFit::SumW2Error(true), RooFit::Minos(true),
    // RooFit::PrintLevel(-1), RooFit::Warnings(false));
    RooFitResult *r =
        pdf.fitTo(*obsData, RooFit::Save(), RooFit::Minos(true),
                  RooFit::PrintLevel(-1), RooFit::Warnings(false));

    double Mt_fit = Mt.getVal();
    double err_hi = Mt.getAsymErrorHi();
    double err_lo = Mt.getAsymErrorLo(); // negative by convention

    std::cout << "\n  *** Pre-calibration result ***\n"
              << "  Mt = " << Mt_fit << " +" << err_hi << " " << err_lo
              << " GeV  (stat)\n";

    // Plot
    TCanvas *c = new TCanvas("cData", "Data mass fit", 800, 600);
    RooPlot *frame = x.frame(RooFit::Title("Data: fitted top mass"));
    obsData->plotOn(frame, RooFit::Name("datapoints"));
    pdf.plotOn(frame, RooFit::LineColor(kRed), RooFit::Name("datafit"));
    frame->GetXaxis()->SetTitle("Fitted top mass [GeV]");
    frame->GetYaxis()->SetTitle("Events / bin");
    frame->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.038);
    latex.DrawLatex(0.15, 0.88, "CMS top mass measurement (parametric fit)");
    latex.DrawLatex(0.15, 0.82,
                    Form("M_{t} = %.2f ^{+%.2f}_{%.2f} GeV (stat)", Mt_fit,
                         err_hi, err_lo));
    c->SaveAs("step2_data_fit.pdf");
    std::cout << "  Saved: step2_data_fit.pdf\n";

    return Mt_fit;
}

// =============================================================================
// STEP 3: Calibration check using shifted pseudo-data
//
//   Shifts the MC dataset by δ = {-6,-3,0,+3,+6} GeV,
//   fits each shifted sample, plots Mt_fit vs Mt_true,
//   fits a calibration line Mt_fit = a + b*Mt_true.
// =============================================================================
struct CalibLine {
    double a, b;
}; // Mt_fit = a + b*Mt_true

CalibLine calibrationCheck(const char *mcFile, RooRealVar &x,
                           const ShapeParams &sp) {
    std::cout << "\n========== STEP 3: Calibration check ==========\n";

    RooDataSet *mcData = loadDataset(mcFile, "mcDataCalib", x);
    if (!mcData || mcData->numEntries() == 0) {
        std::cerr << "ERROR: empty MC dataset for calibration\n";
        return {0., 1.};
    }

    std::vector<double> deltas = {-6., -3., 0., +3., +6.};
    std::vector<double> Mt_true_pts, Mt_fit_pts, Mt_fit_err;

    for (double delta : deltas) {
        double Mt_true = cfg::MT_GEN + delta;

        // Build shifted weighted dataset
        RooRealVar ww("w", "weight", 1., -1e6, 1e6);
        RooDataSet shifted("shifted", "shifted", RooArgSet(x, ww),
                           RooFit::WeightVar(ww));
        for (int i = 0; i < mcData->numEntries(); ++i) {
            const RooArgSet *row = mcData->get(i); // also loads weight
            double val = ((RooRealVar *)row->find("x"))->getVal() + delta;
            if (val < cfg::MASS_MIN || val > cfg::MASS_MAX)
                continue;
            double eventW = mcData->weight(); // weight of current entry
            x.setVal(val);
            ww.setVal(eventW);
            shifted.add(RooArgSet(x, ww), eventW);
        }

        // Fit the shifted sample
        RooRealVar Mt("Mt", "Mt", Mt_true, cfg::MT_FIT_LO, cfg::MT_FIT_HI);
        RooRealVar width("width", "width", sp.width);
        RooRealVar tail("tail", "tail", sp.tail);
        double offsetVal = sp.peak0 - cfg::MT_GEN;
        RooRealVar offset("offset", "offset", offsetVal);
        RooFormulaVar peak("peak", "Mt + offset", RooArgList(Mt, offset));
        RooNovosibirsk pdf("pdf", "pdf", x, peak, width, tail);
        pdf.fitTo(shifted, RooFit::SumW2Error(true), RooFit::PrintLevel(-1),
                  RooFit::Warnings(false));

        std::cout << "  delta = " << std::showpos << delta << std::noshowpos
                  << "  Mt_true = " << Mt_true << "  Mt_fit = " << Mt.getVal()
                  << " +/- " << Mt.getError() << "\n";

        Mt_true_pts.push_back(Mt_true);
        Mt_fit_pts.push_back(Mt.getVal());
        Mt_fit_err.push_back(Mt.getError());
    }

    // Fit calibration line: Mt_fit = a + b * Mt_true
    int N = deltas.size();
    TGraphErrors *cal = new TGraphErrors(N);
    for (int i = 0; i < N; ++i) {
        cal->SetPoint(i, Mt_true_pts[i], Mt_fit_pts[i]);
        cal->SetPointError(i, 0., Mt_fit_err[i]);
    }
    cal->SetMarkerStyle(20);
    cal->SetMarkerColor(kBlue);

    TF1 *line =
        new TF1("calLine", "[0] + [1]*x", cfg::MT_FIT_LO, cfg::MT_FIT_HI);
    line->SetParameters(0., 1.); // start near identity
    cal->Fit(line, "Q");         // Q = quiet

    CalibLine cl{line->GetParameter(0), line->GetParameter(1)};

    std::cout << "\n  Calibration line: Mt_fit = " << cl.a << " + " << cl.b
              << " * Mt_true\n"
              << "  Inversion:        Mt_true = (Mt_fit - " << cl.a << ") / "
              << cl.b << "\n";

    // Also draw the ideal diagonal for reference
    TF1 *diag = new TF1("diag", "x", cfg::MT_FIT_LO, cfg::MT_FIT_HI);
    diag->SetLineStyle(kDashed);
    diag->SetLineColor(kGray + 1);

    TCanvas *c = new TCanvas("cCalib", "Calibration", 800, 600);
    cal->SetTitle("Calibration: Mt_{fit} vs Mt_{true};"
                  "Mt_{true} [GeV];Mt_{fit} [GeV]");
    cal->Draw("AP");
    line->Draw("same");
    diag->Draw("same");

    TLegend *leg = new TLegend(0.15, 0.70, 0.55, 0.88);
    leg->AddEntry(cal, "pseudo-data points", "p");
    leg->AddEntry(line, Form("fit: %.3f + %.4f #cdot M_{t}^{true}", cl.a, cl.b),
                  "l");
    leg->AddEntry(diag, "ideal (slope=1, intercept=0)", "l");
    leg->Draw();
    c->SaveAs("step3_calibration.pdf");
    std::cout << "  Saved: step3_calibration.pdf\n";

    return cl;
}

// =============================================================================
// MAIN entry point
// =============================================================================
void topmass_measurement() {
    const char *mcFile = "../data/proccess_tuples/signal.root"; // your MC file
    const char *dataFile =
        "../data/proccess_tuples/data_tuples.root"; // your data file
    // const char *mcFile = "../data/proccess_tuples/signal_2.root"; // your MC
    // file const char *dataFile =
    //"../data/proccess_tuples/signal_1.root"; // your data file
    // Suppress RooFit banner
    RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

    // Observable: fitted top mass
    RooRealVar x("x", "fitted top mass [GeV]", cfg::MASS_MIN, cfg::MASS_MAX);

    // ── Step 1: extract shape from MC ────────────────────────────────────
    ShapeParams sp = fitMCShape(mcFile, x);

    // ── Step 2: fit data → pre-calibration Mt ────────────────────────────
    double Mt_precalib = fitData(dataFile, x, sp);
    if (Mt_precalib < 0.)
        return;

    // ── Step 3: calibration check ────────────────────────────────────────
    CalibLine cl = calibrationCheck(mcFile, x, sp);

    // ── Apply calibration ────────────────────────────────────────────────
    double Mt_calibrated = (Mt_precalib - cl.a) / cl.b;

    std::cout << "\n========================================\n"
              << "  FINAL RESULT\n"
              << "  Mt (pre-calib)  = " << Mt_precalib << " GeV\n"
              << "  Mt (calibrated) = " << Mt_calibrated << " GeV\n"
              << "  (stat uncertainty from Step 2 Minos errors,\n"
              << "   also divide by calib slope b = " << cl.b << ")\n"
              << "========================================\n";
}

// runAll("../data/proccess_tuples/signal.root","../data/proccess_tuples/signal.root")
