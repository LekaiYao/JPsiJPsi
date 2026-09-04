#include "TFile.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "TParameter.h"
#include "TSystem.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

namespace {

string validationEdgeLabel(double value) {
    ostringstream label;
    label << setprecision(12) << value;
    string result = label.str();
    for(char &character : result) {
        if(character == '-') character = 'm';
        else if(character == '.') character = 'p';
    }
    return result;
}

string differentialResultPath(
    const string &resultTag,
    const string &variable,
    double minimum,
    double maximum,
    bool isReference
) {
    return "fit_results/" + resultTag + "/" + variable + "_" +
        validationEdgeLabel(minimum) + "_" + validationEdgeLabel(maximum) +
        (isReference ? "_reference.root" : "_nominal.root");
}

bool acceptedFit(const RooFitResult *result) {
    return result && result->status() == 0 && result->covQual() == 3 &&
        std::isfinite(result->edm()) && result->edm() < 0.01;
}

void failValidation(const string &message) {
    cerr << message << endl;
    gSystem->Exit(1);
}

double approvedCorrectionRelative(
    const string &variable,
    double minimum,
    double maximum,
    bool isTotal
) {
    if(isTotal) return 0.08805293089840492;
    struct CorrectionBin {
        const char *variable;
        double minimum;
        double maximum;
        double relative;
    };
    static const CorrectionBin bins[] = {
        {"delta_y", 0.0, 0.5, 0.07580713063},
        {"delta_y", 0.5, 1.0, 0.1199520085},
        {"delta_y", 1.0, 1.5, 0.2144261991},
        {"delta_y", 1.5, 2.0, 0.2540822226},
        {"delta_y", 2.0, 2.5, 0.1081856788},
        {"delta_y", 2.5, 4.0, 0.1286451326},
        {"delta_phi", 0.0, 0.3927, 0.1091414027},
        {"delta_phi", 0.3927, 0.7854, 0.07074449717},
        {"delta_phi", 0.7854, 1.1781, 0.07730933945},
        {"delta_phi", 1.1781, 1.5708, 0.1180495287},
        {"delta_phi", 1.5708, 1.9635, 0.3488942865},
        {"delta_phi", 1.9635, 2.3562, 0.5390212072},
        {"delta_phi", 2.3562, 2.7489, 0.109167718},
        {"delta_phi", 2.7489, 3.1416, 0.1680099942},
        {"evt_mass", 7.5, 17.5, 0.05468848022},
        {"evt_mass", 17.5, 27.5, 0.07334206956},
        {"evt_mass", 27.5, 37.5, 0.2080534422},
        {"evt_mass", 37.5, 47.5, 0.2677139347},
        {"evt_mass", 47.5, 57.5, 0.3871274033},
        {"evt_mass", 57.5, 67.5, 0.1243808473},
        {"evt_mass", 67.5, 107.5, 0.06072998757},
        {"evt_y", 0.0, 0.4, 0.09960871584},
        {"evt_y", 0.4, 0.8, 0.07809099214},
        {"evt_y", 0.8, 1.2, 0.1116045116},
        {"evt_y", 1.2, 1.6, 0.09404106481},
        {"evt_y", 1.6, 2.0, 0.1042348012},
        {"evt_pt", 0.0, 5.0, 0.09408537837},
        {"evt_pt", 5.0, 10.0, 0.3691573596},
        {"evt_pt", 10.0, 15.0, 0.2704994845},
        {"evt_pt", 15.0, 20.0, 0.1183484058},
        {"evt_pt", 20.0, 25.0, 0.08509011811},
        {"evt_pt", 25.0, 30.0, 0.07788903208},
        {"evt_pt", 30.0, 35.0, 0.08022597139},
        {"evt_pt", 35.0, 40.0, 0.09077947752},
        {"evt_pt", 40.0, 80.0, 0.11180612}
    };
    for(const CorrectionBin &bin : bins) {
        if(variable == bin.variable &&
           std::fabs(minimum - bin.minimum) < 1.0e-9 &&
           std::fabs(maximum - bin.maximum) < 1.0e-9) {
            return bin.relative;
        }
    }
    return -1.0;
}

const RooRealVar *signalYield(const RooFitResult *result) {
    return result ? dynamic_cast<const RooRealVar *>(
        result->floatParsFinal().find("n_P_P")) : 0;
}

int boundaryFallbackFlag(TFile &file, bool isTotal, const char *name) {
    if(isTotal) return 0;
    TParameter<int> *flag = dynamic_cast<TParameter<int> *>(
        file.Get(name));
    return flag ? flag->GetVal() : -1;
}

} // namespace

void Fit_valid(
    string var="total",
    double vmin=0.0,
    double vmax=0.0,
    string resultTag="corrected_error_nominal_root640_unseeded",
    string totalResultDirectory=""
) {
    const bool isTotal = var == "total";
    const string totalPrefix = totalResultDirectory.empty() ?
        "" : totalResultDirectory + "/";
    const string nominalPath = isTotal ?
        totalPrefix + "Fit_4D_tot_native_asymptotic_unseeded.root" :
        differentialResultPath(resultTag, var, vmin, vmax, false);
    const string referencePath = isTotal ?
        totalPrefix + "Fit_4D_tot_ref_native_asymptotic_unseeded.root" :
        differentialResultPath(resultTag, var, vmin, vmax, true);

    TFile nominalFile(nominalPath.c_str(), "READ");
    TFile referenceFile(referencePath.c_str(), "READ");
    if(nominalFile.IsZombie() || referenceFile.IsZombie()) {
        failValidation("Cannot open nominal/reference fit result for " + var);
        return;
    }
    RooFitResult *nominal = dynamic_cast<RooFitResult *>(
        nominalFile.Get("fit_result_native_asymptotic"));
    RooFitResult *reference = dynamic_cast<RooFitResult *>(
        referenceFile.Get("fit_result_native_asymptotic"));
    if(!acceptedFit(nominal) || !acceptedFit(reference)) {
        failValidation("Nominal/reference fit result fails the quality gate for " + var);
        return;
    }
    const RooRealVar *nominalYield = signalYield(nominal);
    const RooRealVar *referenceYield = signalYield(reference);
    if(!nominalYield || !referenceYield || nominalYield->getVal() <= 0.0) {
        failValidation("Cannot read a positive n_P_P result for " + var);
        return;
    }
    const char *sigCombFallbackName =
        "boundary_fallback_n_Sig_Comb_and_n_Comb_Sig_shared_fixed_zero";
    const char *combCombFallbackName =
        "boundary_fallback_n_Comb_Comb_fixed_zero";
    const int nominalSigCombFallback = boundaryFallbackFlag(
        nominalFile, isTotal, sigCombFallbackName);
    const int nominalCombCombFallback = boundaryFallbackFlag(
        nominalFile, isTotal, combCombFallbackName);
    const int referenceSigCombFallback = boundaryFallbackFlag(
        referenceFile, isTotal, sigCombFallbackName);
    const int referenceCombCombFallback = boundaryFallbackFlag(
        referenceFile, isTotal, combCombFallbackName);
    if(nominalSigCombFallback < 0 || nominalCombCombFallback < 0 ||
       referenceSigCombFallback < 0 || referenceCombCombFallback < 0) {
        failValidation("Missing boundary-fallback metadata for " + var);
        return;
    }

    const double luminosityFb = 36.684;
    const double branchingFraction = 0.05961;
    const double branchingRelative = 0.011;
    const double luminosityRelative = 0.012;
    const double lifetimeRelative = 0.003;
    const double correctionRelative = approvedCorrectionRelative(
        var, vmin, vmax, isTotal);
    if(correctionRelative < 0.0) {
        failValidation("No approved correction systematic for " + var);
        return;
    }
    const double width = isTotal ? 1.0 : vmax - vmin;
    if(width <= 0.0) {
        failValidation("Invalid bin width for " + var);
        return;
    }
    const double conversion = 1.0e-3 /
        (luminosityFb * branchingFraction * branchingFraction * width);
    const double yield = nominalYield->getVal();
    const double yieldError = nominalYield->getError();
    const double referenceValue = referenceYield->getVal();
    const double fitterRelative = std::fabs(referenceValue - yield) / yield;
    const double totalSystematicRelative = std::sqrt(
        branchingRelative * branchingRelative +
        luminosityRelative * luminosityRelative +
        correctionRelative * correctionRelative +
        lifetimeRelative * lifetimeRelative +
        fitterRelative * fitterRelative);

    cout << setprecision(15)
         << (isTotal ? "total" : "differential") << ","
         << var << "," << vmin << "," << vmax << ","
         << yield << "," << yieldError << "," << referenceValue << ","
         << fitterRelative << ","
         << yield * conversion << "," << yieldError * conversion << ","
         << yield * conversion * fitterRelative << ","
         << branchingRelative << "," << luminosityRelative << ","
         << correctionRelative << "," << lifetimeRelative << ","
         << totalSystematicRelative << ","
         << yield * conversion * totalSystematicRelative << ","
         << nominalSigCombFallback << ","
         << nominalCombCombFallback << ","
         << referenceSigCombFallback << ","
         << referenceCombCombFallback << ","
         << nominal->status() << "," << nominal->covQual() << ","
         << nominal->edm() << ","
         << reference->status() << "," << reference->covQual() << ","
         << reference->edm() << endl;
}
