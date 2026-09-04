#include "Plot_4D.hpp"
#include "TObjString.h"
#include "TParameter.h"
#include "TSystem.h"
#include <cmath>
#include <exception>
#include <iomanip>
#include <sstream>

namespace {

std::string edgeLabel(double value) {
    std::ostringstream label;
    label << std::setprecision(12) << value;
    std::string result = label.str();
    for(char &character : result) {
        if(character == '-') character = 'm';
        else if(character == '.') character = 'p';
    }
    return result;
}

std::string differentialBinTag(
    const std::string &variable,
    double minimum,
    double maximum,
    bool isReference
) {
    return variable + "_" + edgeLabel(minimum) + "_" + edgeLabel(maximum) +
        (isReference ? "_reference" : "_nominal");
}

bool acceptedDifferentialFit(const RooFitResult *result) {
    return result && result->status() == 0 && result->covQual() == 3 &&
        std::isfinite(result->edm()) && result->edm() < 0.01;
}

void failDifferentialFit(const std::string &message, int exitCode=1) {
    std::cerr << message << std::endl;
    gSystem->Exit(exitCode);
}

} // namespace

void Fit_4D_diff(
    string var,
    double vmin,
    double vmax,
    bool isRef=false,
    bool useNativeAsymptotic=true,
    string resultTag="corrected_error_nominal_root640_unseeded",
    string weightDataFileName="WeightData.root",
    string totalModelDirectory="",
    bool writeGenericWorkspace=true,
    bool allowStrategy2=true
) {
    if(!useNativeAsymptotic) {
        failDifferentialFit(
            "The corrected-error nominal differential fit requires native AsymptoticError(true)",
            2);
        return;
    }
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 40, 0)
    failDifferentialFit(
        "The corrected-error nominal differential fit requires ROOT >= 6.40",
        2);
    return;
#endif
    if(vmax <= vmin || resultTag.empty()) {
        failDifferentialFit("Invalid differential-bin range or empty result tag", 2);
        return;
    }

    const string ref = isRef ? "_ref" : "";
    const string binTag = differentialBinTag(var, vmin, vmax, isRef);
    const string resultDirectory = "fit_results/" + resultTag;
    const string resultFileName = resultDirectory + "/" + binTag + ".root";

    // Define variables.
    RooRealVar Jpsi_mass1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25);
    RooRealVar Jpsi_mass2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25);
    RooRealVar Jpsi_ctau1("Jpsi_ctau1", "Jpsi_ctau1", -0.03, 0.16);
    RooRealVar Jpsi_ctau2("Jpsi_ctau2", "Jpsi_ctau2", -0.03, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 1000);
    RooRealVar kine(var.c_str(), var.c_str(), vmin, vmax);
    RooArgSet variables;
    variables.add(Jpsi_mass1);
    variables.add(Jpsi_mass2);
    variables.add(Jpsi_ctau1);
    variables.add(Jpsi_ctau2);
    variables.add(evt_weight);
    variables.add(kine);

    // Preserve the historical open-bin selection used by the AN chain.
    const string sel = var + " > " + to_string(vmin) +
        " && " + var + " < " + to_string(vmax);
    TFile dataFile(weightDataFileName.c_str(), "READ");
    if(dataFile.IsZombie()) {
        failDifferentialFit("Cannot open " + weightDataFileName);
        return;
    }
    TTree *dataTree = dynamic_cast<TTree *>(dataFile.Get("data"));
    if(!dataTree) {
        failDifferentialFit("Cannot find tree 'data' in WeightData.root");
        return;
    }
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 40, 0)
    RooDataSet data(
        "data", "data", variables, Import(*dataTree), Cut(sel.c_str()),
        WeightVar("evt_weight"));
#endif
    if(data.numEntries() == 0) {
        failDifferentialFit("The selected differential bin contains no events");
        return;
    }

    // Shape parameters are fixed to the matching no-seed total-fit workspace.
    const string totalModelBaseName = "Model_4D_tot" + ref + ".root";
    const string totalModelFileName = totalModelDirectory.empty() ?
        totalModelBaseName : totalModelDirectory + "/" + totalModelBaseName;
    TFile totalModelFile(totalModelFileName.c_str(), "READ");
    if(totalModelFile.IsZombie()) {
        failDifferentialFit("Cannot open " + totalModelFileName);
        return;
    }
    RooWorkspace *wsp = dynamic_cast<RooWorkspace *>(totalModelFile.Get("wsp"));
    RooAddPdf *pdfAll = wsp ? dynamic_cast<RooAddPdf *>(wsp->pdf("pdf_all")) : 0;
    if(!wsp || !pdfAll) {
        failDifferentialFit("Cannot load workspace/pdf_all from " + totalModelFileName);
        return;
    }
    RooAbsPdf *pdf_P_P = wsp->pdf("pdf_P_P");
    RooAbsPdf *pdf_P_NP = wsp->pdf("pdf_P_NP");
    RooAbsPdf *pdf_NP_P = wsp->pdf("pdf_NP_P");
    RooAbsPdf *pdf_NP_NP = wsp->pdf("pdf_NP_NP");
    RooAbsPdf *pdf_Sig_Comb = wsp->pdf("pdf_Sig_Comb");
    RooAbsPdf *pdf_Comb_Sig = wsp->pdf("pdf_Comb_Sig");
    RooAbsPdf *pdf_Comb_Comb = wsp->pdf("pdf_Comb_Comb");
    if(!pdf_P_P || !pdf_P_NP || !pdf_NP_P || !pdf_NP_NP ||
       !pdf_Sig_Comb || !pdf_Comb_Sig || !pdf_Comb_Comb) {
        failDifferentialFit("The total-fit workspace is missing component PDFs");
        return;
    }

    // Keep the original differential-fit yield initial values and ranges,
    // except for the user-approved physical combinatorial-yield boundaries at 0.
    // pdf_Sig_Comb and pdf_Comb_Sig intentionally share n_Sig_Comb in pdf_all.
    RooRealVar *n_P_P = wsp->var("n_P_P");
    RooRealVar *n_P_NP = wsp->var("n_P_NP");
    RooRealVar *n_NP_NP = wsp->var("n_NP_NP");
    RooRealVar *n_Sig_Comb = wsp->var("n_Sig_Comb");
    RooRealVar *n_Comb_Comb = wsp->var("n_Comb_Comb");
    if(!n_P_P || !n_P_NP || !n_NP_NP || !n_Sig_Comb || !n_Comb_Comb) {
        failDifferentialFit("The total-fit workspace is missing yield parameters");
        return;
    }
    n_P_P->setVal(1e3);
    n_P_NP->setVal(1e2);
    n_NP_NP->setVal(1e2);
    n_Sig_Comb->setVal(1e3);
    n_Comb_Comb->setVal(1e2);
    n_P_P->setMin(1);
    n_P_NP->setMin(0);
    n_NP_NP->setMin(0);
    n_Sig_Comb->setMin(0);
    n_Comb_Comb->setMin(0);

    RooAbsData *fitData = data.reduce(
        RooArgSet(Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2));
    RooFitResult *result = 0;
    bool sigCombBoundaryFallback = false;
    bool combCombBoundaryFallback = false;
    int acceptedStrategy = -1;
    bool acceptedOffset = false;
    const int maxFitAttempts = allowStrategy2 ? 4 : 1;
    for(int attempt = 0; attempt < maxFitAttempts; ++attempt) {
        delete result;
        result = 0;
        const int minimizerStrategy = attempt == 0 ? 1 : 2;
        const bool useOffset = attempt > 0;
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 40, 0)
        try {
            result = pdfAll->fitTo(
                *fitData, Save(), Extended(kTRUE), AsymptoticError(kTRUE),
                Strategy(minimizerStrategy), Offset(useOffset));
        } catch(const std::exception &error) {
            cerr << "Fit attempt " << attempt + 1
                 << " threw an exception: " << error.what() << endl;
            const std::string errorMessage = error.what();
            const bool outsideDefaultRange =
                errorMessage.find("outside the default range") != std::string::npos;
            const bool explicitSigCombBoundaryException = outsideDefaultRange &&
                errorMessage.find("n_Sig_Comb") != std::string::npos;
            const bool explicitCombCombBoundaryException = outsideDefaultRange &&
                errorMessage.find("n_Comb_Comb") != std::string::npos;
            const double sigCombBoundaryTolerance = 1.0e-4 * std::max(
                1.0, std::fabs(n_Sig_Comb->getError()));
            const double combCombBoundaryTolerance = 1.0e-4 * std::max(
                1.0, std::fabs(n_Comb_Comb->getError()));
            RooRealVar *fallbackYield = 0;
            bool *fallbackFlag = 0;
            if(!sigCombBoundaryFallback && explicitSigCombBoundaryException) {
                fallbackYield = n_Sig_Comb;
                fallbackFlag = &sigCombBoundaryFallback;
            } else if(!combCombBoundaryFallback && explicitCombCombBoundaryException) {
                fallbackYield = n_Comb_Comb;
                fallbackFlag = &combCombBoundaryFallback;
            } else if(!sigCombBoundaryFallback &&
                      n_Sig_Comb->getMin() <= sigCombBoundaryTolerance &&
                      n_Sig_Comb->getVal() <= sigCombBoundaryTolerance) {
                fallbackYield = n_Sig_Comb;
                fallbackFlag = &sigCombBoundaryFallback;
            } else if(!combCombBoundaryFallback &&
                      n_Comb_Comb->getMin() <= combCombBoundaryTolerance &&
                      n_Comb_Comb->getVal() <= combCombBoundaryTolerance) {
                fallbackYield = n_Comb_Comb;
                fallbackFlag = &combCombBoundaryFallback;
            }
            if(fallbackYield && fallbackFlag) {
                fallbackYield->setVal(0.0);
                fallbackYield->setConstant(kTRUE);
                *fallbackFlag = true;
                cerr << "Activating approved boundary fallback: fixing "
                     << fallbackYield->GetName()
                     << "=0 before the corrected-covariance refit" << endl;
            } else {
                // Preserve the bounded retry for non-boundary exceptions.
                RooRealVar *yields[] = {
                    n_P_P, n_P_NP, n_NP_NP, n_Sig_Comb, n_Comb_Comb
                };
                for(RooRealVar *yield : yields) {
                    if(yield->getVal() < yield->getMin()) {
                        yield->setVal(yield->getMin());
                    }
                    if(yield->getVal() > yield->getMax()) {
                        yield->setVal(yield->getMax());
                    }
                }
            }
        }
#endif
        if(result) {
            cout << "Fit attempt " << attempt + 1
                 << ": strategy=" << minimizerStrategy
                 << ", offset=" << (useOffset ? 1 : 0)
                 << ": status=" << result->status()
                 << ", covQual=" << result->covQual()
                 << ", EDM=" << result->edm() << endl;
        }
        if(acceptedDifferentialFit(result)) {
            acceptedStrategy = minimizerStrategy;
            acceptedOffset = useOffset;
            break;
        }
        const auto atLowerBoundary = [](const RooRealVar *parameter) {
            if(!parameter || !parameter->hasMin()) return false;
            const double tolerance = 1.0e-4 * std::max(
                1.0, std::fabs(parameter->getError()));
            return std::fabs(parameter->getMin()) <= tolerance &&
                parameter->getVal() <= parameter->getMin() + tolerance;
        };
        if(!sigCombBoundaryFallback && atLowerBoundary(n_Sig_Comb)) {
            n_Sig_Comb->setVal(0.0);
            n_Sig_Comb->setConstant(kTRUE);
            sigCombBoundaryFallback = true;
            cerr << "Activating approved boundary fallback after failed fit: "
                 << "fixing n_Sig_Comb=0 (shared with n_Comb_Sig)" << endl;
        }
        if(!combCombBoundaryFallback && atLowerBoundary(n_Comb_Comb)) {
            n_Comb_Comb->setVal(0.0);
            n_Comb_Comb->setConstant(kTRUE);
            combCombBoundaryFallback = true;
            cerr << "Activating approved boundary fallback after failed fit: "
                 << "fixing n_Comb_Comb=0" << endl;
        }
    }
    if(!acceptedDifferentialFit(result)) {
        std::ostringstream failure;
        failure << "Differential 4D fit failed the status=0, covQual=3, EDM<0.01 gate "
                << "after " << maxFitAttempts << " attempts for " << var
                << " in (" << vmin << ", " << vmax << ")"
                << (isRef ? " [reference]" : " [nominal]");
        delete result;
        delete fitData;
        failDifferentialFit(failure.str());
        return;
    }

    const string plotPrefix = "fig/native_asymptotic/diff/" + resultTag + "/" +
        (isRef ? "reference/" : "nominal/") + binTag + "/";
    gSystem->mkdir(plotPrefix.c_str(), kTRUE);
    Plot_4D(
        &data, *pdfAll, *pdf_P_P, *pdf_P_NP, *pdf_NP_P, *pdf_NP_NP,
        *pdf_Sig_Comb, *pdf_Comb_Sig, *pdf_Comb_Comb, plotPrefix,
        Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2, "", "pdf");

    // Retain generic workspaces for compatibility and save a unique result
    // artifact so later bins cannot overwrite the fit used for the table.
    if(writeGenericWorkspace) {
        wsp->writeToFile(("Model_4D_diff" + ref + ".root").c_str());
    }
    gSystem->mkdir(resultDirectory.c_str(), kTRUE);
    TFile resultFile(resultFileName.c_str(), "RECREATE");
    result->Write("fit_result_native_asymptotic");
    TObjString variable(var.c_str());
    TObjString selection(sel.c_str());
    TObjString errorConvention("ROOT native AsymptoticError(true), no external seed");
    TObjString totalModel(totalModelFileName.c_str());
    variable.Write("variable");
    selection.Write("selection");
    errorConvention.Write("error_convention");
    totalModel.Write("total_model_workspace");
    TParameter<double>("bin_min", vmin).Write();
    TParameter<double>("bin_max", vmax).Write();
    TParameter<int>("selected_entries", data.numEntries()).Write();
    TParameter<int>("boundary_fallback_n_Comb_Comb_fixed_zero",
        combCombBoundaryFallback ? 1 : 0).Write();
    TParameter<int>("boundary_fallback_n_Sig_Comb_and_n_Comb_Sig_shared_fixed_zero",
        sigCombBoundaryFallback ? 1 : 0).Write();
    TParameter<int>("minimizer_strategy", acceptedStrategy).Write();
    TParameter<int>("offset_enabled", acceptedOffset ? 1 : 0).Write();
    resultFile.Close();

    cout << std::setprecision(15)
         << "Accepted differential result: " << resultFileName
         << ", n_P_P=" << n_P_P->getVal()
         << " +/- " << n_P_P->getError()
         << ", n_Sig_Comb_n_Comb_Sig_shared_fallback="
         << (sigCombBoundaryFallback ? 1 : 0)
         << ", n_Comb_Comb_fallback="
         << (combCombBoundaryFallback ? 1 : 0)
         << ", strategy=" << acceptedStrategy
         << ", offset=" << (acceptedOffset ? 1 : 0) << endl;
    delete result;
    delete fitData;
}
