#include "Plot_4D.hpp"
#include "RooDerivative.h"
#include "RooFormulaVar.h"
#include "TDecompChol.h"
#include "TMatrixDSymEigen.h"
#include "TObjString.h"
#include "RVersion.h"
#include "TSystem.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <vector>

namespace {

struct CovarianceDiagnostics {
    bool finitePositiveDiagonal;
    bool correlationCholesky;
    double minCorrelationEigenvalue;
    double maxCorrelationEigenvalue;
};

class MutableRooFitResult : public RooFitResult {
public:
    explicit MutableRooFitResult(const RooFitResult &other) : RooFitResult(other) {}
    void replaceCovariance(TMatrixDSym &covariance, int covarianceQuality) {
        setCovarianceMatrix(covariance);
        setCovQual(covarianceQuality);
    }
};

bool loadSeedValues(RooAbsPdf &pdf, const RooAbsData &data, const char *seedFileName) {
    TFile seedFile(seedFileName, "READ");
    if(seedFile.IsZombie()) return false;
    RooWorkspace *seedWorkspace = dynamic_cast<RooWorkspace *>(seedFile.Get("wsp"));
    if(!seedWorkspace) return false;
    std::unique_ptr<RooArgSet> targetParameters(pdf.getParameters(data));
    if(!targetParameters.get()) return false;
    targetParameters->assignValueOnly(seedWorkspace->allVars());
    cout << "Initial values loaded from " << seedFileName << endl;
    return true;
}

bool calculateAsymptoticCovariance(
    RooAbsPdf &pdf,
    RooAbsData &data,
    const RooFitResult &fitResult,
    double relativeStep,
    TMatrixDSym &covariance
) {
    RooFormulaVar logpdf("logpdf", "log(pdf)", "log(@0)", RooArgList(pdf));
    std::unique_ptr<RooArgSet> observables(logpdf.getObservables(data));
    std::unique_ptr<RooArgSet> allParameters(logpdf.getParameters(data));
    if(!observables.get() || !allParameters.get()) return false;

    const RooArgList &floated = fitResult.floatParsFinal();
    const int nParameters = floated.getSize();
    if(nParameters <= 0 || fitResult.covarianceMatrix().GetNrows() != nParameters) return false;

    TMatrixDSym scoreOuterProduct(nParameters);
    scoreOuterProduct.Zero();
    std::vector<std::unique_ptr<RooDerivative> > derivatives;
    std::vector<RooRealVar *> internalParameters;
    derivatives.reserve(nParameters);
    internalParameters.reserve(nParameters);

    for(int k = 0; k < nParameters; ++k) {
        RooRealVar *fitParameter = dynamic_cast<RooRealVar *>(floated.at(k));
        RooRealVar *internalParameter = fitParameter ?
            dynamic_cast<RooRealVar *>(allParameters->find(*fitParameter)) : 0;
        if(!fitParameter || !internalParameter ||
           !std::isfinite(fitParameter->getError()) || fitParameter->getError() <= 0.0) {
            cerr << "Cannot construct asymptotic derivative for parameter index " << k << endl;
            return false;
        }
        internalParameter->setVal(fitParameter->getVal());
        const double step = relativeStep * fitParameter->getError();
        derivatives.emplace_back(logpdf.derivative(*internalParameter, *observables, 1, step));
        internalParameters.push_back(internalParameter);
    }

    // The extended likelihood score also contains d(log(N_expected))/d(theta).
    std::vector<double> expectedDerivatives(nParameters, 0.0);
    const double expectedEvents = pdf.expectedEvents(*observables);
    if(std::isfinite(expectedEvents) && expectedEvents > 0.0) {
        for(int k = 0; k < nParameters; ++k) {
            RooRealVar *fitParameter = dynamic_cast<RooRealVar *>(floated.at(k));
            RooRealVar *internalParameter = internalParameters[k];
            const double value = fitParameter->getVal();
            const double step = relativeStep * fitParameter->getError();
            internalParameter->setVal(value + step);
            const double expectedPlus = pdf.expectedEvents(*observables);
            internalParameter->setVal(value - step);
            const double expectedMinus = pdf.expectedEvents(*observables);
            internalParameter->setVal(value);
            if(!std::isfinite(expectedPlus) || !std::isfinite(expectedMinus) ||
               expectedPlus <= 0.0 || expectedMinus <= 0.0) return false;
            expectedDerivatives[k] =
                (std::log(expectedPlus) - std::log(expectedMinus)) / (2.0 * step);
        }
    }

    std::vector<double> scores(nParameters, 0.0);
    for(int event = 0; event < data.numEntries(); ++event) {
        const RooArgSet *row = data.get(event);
        if(!row) return false;
        observables->assignValueOnly(*row);
        for(int k = 0; k < nParameters; ++k) {
            const double derivative = derivatives[k]->getVal();
            internalParameters[k]->setVal(
                dynamic_cast<RooRealVar *>(floated.at(k))->getVal());
            if(!std::isfinite(derivative)) return false;
            scores[k] = derivative + expectedDerivatives[k];
        }
        const double weightSquared = data.weightSquared();
        if(!std::isfinite(weightSquared) || weightSquared < 0.0) return false;
        for(int k = 0; k < nParameters; ++k) {
            for(int l = 0; l <= k; ++l) {
                scoreOuterProduct(k, l) += weightSquared * scores[k] * scores[l];
                if(k != l) scoreOuterProduct(l, k) = scoreOuterProduct(k, l);
            }
        }
    }

    covariance = scoreOuterProduct;
    covariance.Similarity(fitResult.covarianceMatrix());
    return true;
}

CovarianceDiagnostics diagnoseCovariance(const TMatrixDSym &covariance) {
    CovarianceDiagnostics diagnostics = {true, false, 0.0, 0.0};
    const int nParameters = covariance.GetNrows();
    TMatrixDSym correlation(nParameters);
    for(int k = 0; k < nParameters; ++k) {
        const double diagonal = covariance(k, k);
        if(!std::isfinite(diagonal) || diagonal <= 0.0) diagnostics.finitePositiveDiagonal = false;
    }
    if(!diagnostics.finitePositiveDiagonal) return diagnostics;
    for(int k = 0; k < nParameters; ++k) {
        for(int l = 0; l < nParameters; ++l) {
            correlation(k, l) = covariance(k, l) /
                std::sqrt(covariance(k, k) * covariance(l, l));
        }
    }
    TDecompChol cholesky(correlation);
    diagnostics.correlationCholesky = cholesky.Decompose();
    TMatrixDSymEigen eigenSystem(correlation);
    const TVectorD &eigenvalues = eigenSystem.GetEigenValues();
    diagnostics.minCorrelationEigenvalue = eigenvalues[0];
    diagnostics.maxCorrelationEigenvalue = eigenvalues[0];
    for(int k = 1; k < eigenvalues.GetNrows(); ++k) {
        diagnostics.minCorrelationEigenvalue =
            std::min(diagnostics.minCorrelationEigenvalue, eigenvalues[k]);
        diagnostics.maxCorrelationEigenvalue =
            std::max(diagnostics.maxCorrelationEigenvalue, eigenvalues[k]);
    }
    return diagnostics;
}

double maximumRelativeErrorShift(const TMatrixDSym &reference, const TMatrixDSym &variation) {
    double maximumShift = 0.0;
    for(int k = 0; k < reference.GetNrows(); ++k) {
        const double referenceError = std::sqrt(reference(k, k));
        const double variationError = std::sqrt(variation(k, k));
        maximumShift = std::max(maximumShift,
            std::fabs(variationError / referenceError - 1.0));
    }
    return maximumShift;
}

int parameterIndex(const RooArgList &parameters, const char *name) {
    for(int k = 0; k < parameters.getSize(); ++k) {
        if(std::string(parameters.at(k)->GetName()) == name) return k;
    }
    return -1;
}

void applyCovarianceToResultAndModel(
    MutableRooFitResult &fitResult,
    RooAbsPdf &pdf,
    const RooAbsData &data,
    TMatrixDSym &covariance,
    int covarianceQuality
) {
    fitResult.replaceCovariance(covariance, covarianceQuality);
    std::unique_ptr<RooArgSet> modelParameters(pdf.getParameters(data));
    for(int k = 0; k < fitResult.floatParsFinal().getSize(); ++k) {
        RooRealVar *resultParameter =
            dynamic_cast<RooRealVar *>(fitResult.floatParsFinal().at(k));
        RooRealVar *modelParameter = resultParameter ?
            dynamic_cast<RooRealVar *>(modelParameters->find(*resultParameter)) : 0;
        const double correctedError = std::sqrt(covariance(k, k));
        if(resultParameter) resultParameter->setError(correctedError);
        if(modelParameter) modelParameter->setError(correctedError);
    }
}

std::string parameterOrder(const RooArgList &parameters) {
    std::ostringstream order;
    for(int k = 0; k < parameters.getSize(); ++k) {
        if(k) order << ",";
        order << parameters.at(k)->GetName();
    }
    return order.str();
}

} // namespace

void Fit_4D_tot(
    bool isRef=false,
    bool useSumW2=false,
    bool useAsymptotic=false,
    bool useNativeAsymptotic=true,
    bool useSeed=false,
    int minimizerStrategy=2
) {
    if(useNativeAsymptotic && (useSumW2 || useAsymptotic)) {
        cerr << "Native AsymptoticError cannot be combined with another error branch" << endl;
        return;
    }
    if(minimizerStrategy < 0 || minimizerStrategy > 2) {
        cerr << "Minimizer strategy must be 0, 1, or 2" << endl;
        return;
    }
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 40, 0)
    if(useNativeAsymptotic) {
        cerr << "Native AsymptoticError migration requires ROOT >= 6.40" << endl;
        return;
    }
#endif
    // Define variables and read input file
    RooRealVar Jpsi_mass1("Jpsi_mass1", "Jpsi_mass1", 2.95, 3.25);
    RooRealVar Jpsi_mass2("Jpsi_mass2", "Jpsi_mass2", 2.95, 3.25);
    RooRealVar Jpsi_ctau1("Jpsi_ctau1", "Jpsi_ctau1", -0.03, 0.16);
    RooRealVar Jpsi_ctau2("Jpsi_ctau2", "Jpsi_ctau2", -0.03, 0.16);
    RooRealVar evt_weight("evt_weight", "evt_weight", 0, 1000);
    RooArgSet variables;
    variables.add(Jpsi_mass1);
    variables.add(Jpsi_mass2);
    variables.add(Jpsi_ctau1);
    variables.add(Jpsi_ctau2);
    variables.add(evt_weight);
    TFile *dataFile = new TFile("WeightData.root", "READ");
    TTree *dataTree = (TTree*)dataFile->Get("data");
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 40, 0)
    RooDataSet *data = new RooDataSet(
        "data", "data", variables, Import(*dataTree), WeightVar("evt_weight"));
#else
    RooDataSet *data = new RooDataSet("data", "data", dataTree, variables, "", "evt_weight");
#endif
    
    // Define J/psi Mass p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mean("Jpsi_mean", "Jpsi_mean", 3.0969, 3.05, 3.15);
    RooRealVar Jpsi_devia1("Jpsi_devia1", "Jpsi_devia1", 0.01, 0, 0.03);
    RooRealVar Jpsi_alpha1("Jpsi_alpha1", "Jpsi_alpha1", 1, 0.1, 5);
    RooRealVar Jpsi_nx1("Jpsi_nx1", "Jpsi_nx1", 1, 0, 50);
    RooRealVar Jpsi_devia2("Jpsi_devia2", "Jpsi_devia2", 0.05, 0.02, 0.08);
    RooRealVar Jpsi_ratio("Jpsi_ratio", "Jpsi_ratio", 0.6, 0, 1);
    RooGaussian Jpsi_gaussian_1("Jpsi_gaussian_1", "Jpsi_gaussian_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia2);
    RooGaussian Jpsi_gaussian_2("Jpsi_gaussian_2", "Jpsi_gaussian_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia2);
    RooCBShape Jpsi_crysBall_1("Jpsi_crysBall_1", "Jpsi_crysBall_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    RooCBShape Jpsi_crysBall_2("Jpsi_crysBall_2", "Jpsi_crysBall_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia1, Jpsi_alpha1, Jpsi_nx1);
    RooGaussian Jpsi_gaussRef_1("Jpsi_gaussRef_1", "Jpsi_gaussRef_1", Jpsi_mass1, Jpsi_mean, Jpsi_devia1);
    RooGaussian Jpsi_gaussRef_2("Jpsi_gaussRef_2", "Jpsi_gaussRef_2", Jpsi_mass2, Jpsi_mean, Jpsi_devia1);
    RooAbsPdf &Jpsi_core_1 = isRef ? static_cast<RooAbsPdf &>(Jpsi_gaussRef_1) : static_cast<RooAbsPdf &>(Jpsi_crysBall_1);
    RooAbsPdf &Jpsi_core_2 = isRef ? static_cast<RooAbsPdf &>(Jpsi_gaussRef_2) : static_cast<RooAbsPdf &>(Jpsi_crysBall_2);
    RooAddPdf JpsiMassSig1("JpsiMassSig1", "JpsiMassSig1", RooArgList(Jpsi_core_1, Jpsi_gaussian_1), Jpsi_ratio);
    RooAddPdf JpsiMassSig2("JpsiMassSig2", "JpsiMassSig2", RooArgList(Jpsi_core_2, Jpsi_gaussian_2), Jpsi_ratio);
    // Background p.d.f.
    RooChebychev JpsiMassComb1("JpsiMassComb1", "JpsiMassComb1", Jpsi_mass1, RooArgList());
    RooChebychev JpsiMassComb2("JpsiMassComb2", "JpsiMassComb2", Jpsi_mass2, RooArgList());

    // Define J/psi Ctau p.d.f.
    // Signal p.d.f.
    RooRealVar Jpsi_mu1("Jpsi_mu1", "Jpsi_mu1", 0, -0.005, 0.005);
    RooRealVar Jpsi_sigma1("Jpsi_sigma1", "Jpsi_sigma1", 0.001, 0, 0.003);
    RooRealVar Jpsi_sigma2("Jpsi_sigma2", "Jpsi_sigma2", 0.004, 0.002, 0.008);
    RooRealVar Jpsi_prop1("Jpsi_prop1", "Jpsi_prop1", 0.5, 0, 1);
    RooGaussian Jpsi_gauss1_1("Jpsi_gauss1_1", "Jpsi_gauss1_1", Jpsi_ctau1, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2_1("Jpsi_gauss2_1", "Jpsi_gauss2_1", Jpsi_ctau1, Jpsi_mu1, Jpsi_sigma2);
    RooGaussian Jpsi_gauss1_2("Jpsi_gauss1_2", "Jpsi_gauss1_2", Jpsi_ctau2, Jpsi_mu1, Jpsi_sigma1);
    RooGaussian Jpsi_gauss2_2("Jpsi_gauss2_2", "Jpsi_gauss2_2", Jpsi_ctau2, Jpsi_mu1, Jpsi_sigma2);
    RooAddPdf JpsiCtauSig1("JpsiCtauSig1", "JpsiCtauSig1", RooArgList(Jpsi_gauss1_1, Jpsi_gauss2_1), Jpsi_prop1);
    RooAddPdf JpsiCtauSig2("JpsiCtauSig2", "JpsiCtauSig2", RooArgList(Jpsi_gauss1_2, Jpsi_gauss2_2), Jpsi_prop1);
    // Background p.d.f.
    RooRealVar Jpsi_sigma3("Jpsi_sigma3", "Jpsi_sigma3", 0.001, 0, 0.01);
    RooRealVar Jpsi_coef1("Jpsi_coef1", "Jpsi_coef1", 0.06, 0.01, 0.1);
    RooGExpModel JpsiCtauBkg1("JpsiCtauBkg1", "JpsiCtauBkg1", Jpsi_ctau1, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    RooGExpModel JpsiCtauBkg2("JpsiCtauBkg2", "JpsiCtauBkg2", Jpsi_ctau2, Jpsi_sigma3, Jpsi_coef1, false, RooGExpModel::Type::Flipped);
    // Combinatorial signal p.d.f.
    RooRealVar Jpsi_prop2("Jpsi_prop2", "Jpsi_prop2", 0.5, 0, 1);
    RooAddPdf JpsiCtauCombSig1("JpsiCtauCombSig1", "JpsiCtauCombSig1", RooArgList(JpsiCtauSig1, JpsiCtauBkg1), Jpsi_prop2);
    RooAddPdf JpsiCtauCombSig2("JpsiCtauCombSig2", "JpsiCtauCombSig2", RooArgList(JpsiCtauSig2, JpsiCtauBkg2), Jpsi_prop2);
    // Combinatorial background p.d.f.
    RooRealVar Jpsi_sigma7("Jpsi_sigma7", "Jpsi_sigma7", 0.003, 0.0008, 0.02);
    RooRealVar Jpsi_coef3("Jpsi_coef3", "Jpsi_coef3", 0.06, 0.01, 0.2);
    RooGExpModel JpsiCtauCombBkg1("JpsiCtauCombBkg1", "JJpsiCtauCombBkg1", Jpsi_ctau1, Jpsi_sigma7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    RooGExpModel JpsiCtauCombBkg2("JpsiCtauCombBkg2", "JJpsiCtauCombBkg2", Jpsi_ctau2, Jpsi_sigma7, Jpsi_coef3, false, RooGExpModel::Type::Flipped);
    
    // Form 4-dim p.d.f.
    Int_t N = 100000;
    // J/psi1 Mass + J/psi2 Mass dimension
    RooProdPdf pdf_mass_SigSig("pdf_mass_SigSig", "pdf_mass_SigSig", JpsiMassSig1, JpsiMassSig2);
    RooProdPdf pdf_mass_SigComb("pdf_mass_SigComb", "pdf_mass_SigComb", JpsiMassSig1, JpsiMassComb2);
    RooProdPdf pdf_mass_CombSig("pdf_mass_CombSig", "pdf_mass_CombSig", JpsiMassComb1, JpsiMassSig2);
    RooProdPdf pdf_mass_CombComb("pdf_mass_CombComb", "pdf_mass_CombComb", JpsiMassComb1, JpsiMassComb2);
    // J/psi1 Ctau + J/psi2 Ctau dimension
    RooProdPdf pdf_ctau_PP("pdf_ctau_PP", "pdf_ctau_PP", JpsiCtauSig1, JpsiCtauSig2);
    RooProdPdf pdf_ctau_PNP("pdf_ctau_PNP", "pdf_ctau_PNP", JpsiCtauSig1, JpsiCtauBkg2);
    RooProdPdf pdf_ctau_NPP("pdf_ctau_NPP", "pdf_ctau_NPP", JpsiCtauBkg1, JpsiCtauSig2);
    RooProdPdf pdf_ctau_NPNP("pdf_ctau_NPNP", "pdf_ctau_NPNP", JpsiCtauBkg1, JpsiCtauBkg2);
    RooProdPdf pdf_ctau_SigBkg("pdf_ctau_SigBkg", "pdf_ctau_SigBkg", JpsiCtauCombSig1, JpsiCtauCombBkg2);
    RooProdPdf pdf_ctau_BkgSig("pdf_ctau_BkgSig", "pdf_ctau_BkgSig", JpsiCtauCombBkg1, JpsiCtauCombSig2);
    RooProdPdf pdf_ctau_BkgBkg("pdf_ctau_BkgBkg", "pdf_ctau_BkgBkg", JpsiCtauCombBkg1, JpsiCtauCombBkg2);
    // Combine Mass and Ctau dimensions
    RooProdPdf pdf_P_P("pdf_P_P", "pdf_P_P", pdf_mass_SigSig, pdf_ctau_PP);
    RooProdPdf pdf_P_NP("pdf_P_NP", "pdf_P_NP", pdf_mass_SigSig, pdf_ctau_PNP);
    RooProdPdf pdf_NP_P("pdf_NP_P", "pdf_NP_P", pdf_mass_SigSig, pdf_ctau_NPP);
    RooProdPdf pdf_NP_NP("pdf_NP_NP", "pdf_NP_NP", pdf_mass_SigSig, pdf_ctau_NPNP);
    RooProdPdf pdf_Sig_Comb("pdf_Sig_Comb", "pdf_Sig_Comb", pdf_mass_SigComb, pdf_ctau_SigBkg);
    RooProdPdf pdf_Comb_Sig("pdf_Comb_Sig", "pdf_Comb_Sig", pdf_mass_CombSig, pdf_ctau_BkgSig);
    RooProdPdf pdf_Comb_Comb("pdf_Comb_Comb", "pdf_Comb_Comb", pdf_mass_CombComb, pdf_ctau_BkgBkg);
    RooRealVar n_P_P("n_P_P", "n_P_P", 1e3, 1e1, N);
    RooRealVar n_P_NP("n_P_NP", "n_P_NP", 5e2, 1, N);
    RooRealVar n_NP_NP("n_NP_NP", "n_NP_NP", 1e3, 1, N);
    RooRealVar n_Sig_Comb("n_Sig_Comb", "n_Sig_Comb", 1e4, 1, N);
    RooRealVar n_Comb_Comb("n_Comb_Comb", "n_Comb_Comb", 1e2, 1, N);
    RooAddPdf pdf_all("pdf_all", "pdf_all",
        RooArgList(pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb),
        RooArgList(n_P_P, n_P_NP, n_P_NP, n_NP_NP, n_Sig_Comb, n_Sig_Comb, n_Comb_Comb)
    );
    RooAbsData *fitData = data->reduce(
        RooArgSet(Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2));
    const char *seedFile = isRef ?
        "Model_4D_tot_ref_previous_method.root" : "Model_4D_tot_previous_method.root";
    if(useAsymptotic) loadSeedValues(pdf_all, *fitData, seedFile);
    if(useNativeAsymptotic && useSeed &&
       !loadSeedValues(pdf_all, *fitData, seedFile)) {
        cerr << "Cannot load native-AsymptoticError seed workspace " << seedFile << endl;
        delete fitData;
        return;
    }
    if(useNativeAsymptotic) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 40, 0)
        RooFitResult *prefitResult = 0;
        const int maximumPrefitAttempts = 6;
        for(int attempt = 0; attempt < maximumPrefitAttempts; ++attempt) {
            delete prefitResult;
            prefitResult = pdf_all.fitTo(
                *fitData, Save(), Extended(kTRUE), SumW2Error(kFALSE),
                Strategy(minimizerStrategy));
            if(prefitResult) {
                cout << "Central-value prefit attempt " << attempt + 1
                     << ": strategy=" << minimizerStrategy
                     << ", status=" << prefitResult->status()
                     << ", covQual=" << prefitResult->covQual()
                     << ", EDM=" << prefitResult->edm() << endl;
            }
            if(prefitResult && !prefitResult->status() &&
               prefitResult->edm() < 0.01) break;
        }
        if(!prefitResult || prefitResult->status() ||
           prefitResult->edm() >= 0.01) {
            cerr << "Total 4D central-value prefit failed" << endl;
            delete prefitResult;
            delete fitData;
            return;
        }
        delete prefitResult;
#endif
    }
    RooFitResult *res = 0;
    // evt_weight is an inverse acceptance/efficiency correction.  The
    // weighted likelihood determines the central values; SumW2Error(kTRUE)
    // applies RooFit's squared-weight Hessian correction to the covariance.
    // This covariance correction is deterministic at a fitted parameter
    // point, so do not repeat a failed Hessian on the same state.
    const int maxFitAttempts = (useAsymptotic || useNativeAsymptotic) ? 6 :
        (useSumW2 ? 1 : 50);
    for(int attempt = 0; attempt < maxFitAttempts; ++attempt) {
        delete res;
        if(useNativeAsymptotic) {
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 40, 0)
            res = pdf_all.fitTo(
                *fitData, Save(), Extended(kTRUE), AsymptoticError(kTRUE),
                Strategy(minimizerStrategy));
#endif
        } else if(useAsymptotic) {
            // Keep the original weighted extended likelihood and obtain its
            // inverse Hessian.  The sandwich correction is calculated below.
            res = pdf_all.fitTo(
                *fitData, Save(), Extended(kTRUE), SumW2Error(kFALSE),
                Strategy(minimizerStrategy));
        } else if(useSumW2) {
            res = pdf_all.fitTo(
                *fitData, Save(), Extended(kTRUE), SumW2Error(kTRUE),
                Strategy(minimizerStrategy));
        } else {
            // Reproduce the previous AN fit convention exactly: weighted
            // central likelihood with RooFit's default SumW2Error(false).
            res = pdf_all.fitTo(*fitData, Save(), Strategy(minimizerStrategy));
        }
        if(res) {
            cout << "Fit attempt " << attempt + 1
                 << ": strategy=" << minimizerStrategy
                 << ": status=" << res->status()
                 << ", covQual=" << res->covQual()
                 << ", EDM=" << res->edm() << endl;
        }
        if(res && !res->status() && res->edm()<0.01) break;
    }
    if(!res || res->status() || res->edm()>=0.01) {
        cerr << "Total 4D fit failed after " << maxFitAttempts << " attempts" << endl;
        delete res;
        delete fitData;
        return;
    }

    RooFitResult *correctedResult = 0;
    TMatrixDSym asymptoticCovariance(res->floatParsFinal().getSize());
    TMatrixDSym covarianceHalfStep(res->floatParsFinal().getSize());
    TMatrixDSym covarianceDoubleStep(res->floatParsFinal().getSize());
    CovarianceDiagnostics covarianceDiagnostics = {false, false, 0.0, 0.0};
    double halfStepMaximumShift = 0.0;
    double doubleStepMaximumShift = 0.0;
    double nPPErrorHalfStep = 0.0;
    double nPPErrorNominalStep = 0.0;
    double nPPErrorDoubleStep = 0.0;
    if(useAsymptotic) {
        const bool nominalOK = calculateAsymptoticCovariance(
            pdf_all, *fitData, *res, 1.0e-4, asymptoticCovariance);
        const bool halfStepOK = calculateAsymptoticCovariance(
            pdf_all, *fitData, *res, 0.5e-4, covarianceHalfStep);
        const bool doubleStepOK = calculateAsymptoticCovariance(
            pdf_all, *fitData, *res, 2.0e-4, covarianceDoubleStep);
        if(!nominalOK || !halfStepOK || !doubleStepOK) {
            cerr << "Asymptotic covariance calculation failed" << endl;
            delete res;
            delete fitData;
            return;
        }
        covarianceDiagnostics = diagnoseCovariance(asymptoticCovariance);
        const CovarianceDiagnostics halfDiagnostics = diagnoseCovariance(covarianceHalfStep);
        const CovarianceDiagnostics doubleDiagnostics = diagnoseCovariance(covarianceDoubleStep);
        if(!covarianceDiagnostics.finitePositiveDiagonal ||
           !covarianceDiagnostics.correlationCholesky ||
           !halfDiagnostics.finitePositiveDiagonal || !halfDiagnostics.correlationCholesky ||
           !doubleDiagnostics.finitePositiveDiagonal || !doubleDiagnostics.correlationCholesky) {
            cerr << "Asymptotic covariance failed positive-definiteness checks" << endl;
            delete res;
            delete fitData;
            return;
        }
        halfStepMaximumShift = maximumRelativeErrorShift(asymptoticCovariance, covarianceHalfStep);
        doubleStepMaximumShift = maximumRelativeErrorShift(asymptoticCovariance, covarianceDoubleStep);
        const int nPPIndex = parameterIndex(res->floatParsFinal(), "n_P_P");
        if(nPPIndex < 0) {
            cerr << "n_P_P is absent from the floating-parameter list" << endl;
            delete res;
            delete fitData;
            return;
        }
        nPPErrorHalfStep = std::sqrt(covarianceHalfStep(nPPIndex, nPPIndex));
        nPPErrorNominalStep = std::sqrt(asymptoticCovariance(nPPIndex, nPPIndex));
        nPPErrorDoubleStep = std::sqrt(covarianceDoubleStep(nPPIndex, nPPIndex));
        MutableRooFitResult mutableResult(*res);
        applyCovarianceToResultAndModel(
            mutableResult, pdf_all, *fitData, asymptoticCovariance, res->covQual());
        correctedResult = new RooFitResult(mutableResult);
        correctedResult->SetName("fit_result_asymptotic");
        correctedResult->SetTitle("Weighted 4D fit with asymptotic sandwich covariance");
    }

    // Draw data point and p.d.f. curve
    string prefix;
    if(useNativeAsymptotic) {
        const string seedLabel = useSeed ? "seeded/" : "unseeded/";
        prefix = isRef ? "fig/native_asymptotic/ref/" + seedLabel :
            "fig/native_asymptotic/" + seedLabel;
    } else if(useAsymptotic) prefix = isRef ? "fig/asymptotic/ref/" : "fig/asymptotic/";
    else if(!useSumW2) prefix = isRef ? "fig/previous_method/ref/" : "fig/previous_method/";
    else prefix = isRef ? "fig/ref/" : "fig/";
    gSystem->mkdir(prefix.c_str(), kTRUE);
    Plot_4D(data, pdf_all, pdf_P_P, pdf_P_NP, pdf_NP_P, pdf_NP_NP, pdf_Sig_Comb, pdf_Comb_Sig, pdf_Comb_Comb,
        prefix, Jpsi_mass1, Jpsi_mass2, Jpsi_ctau1, Jpsi_ctau2, "", "pdf");
    pdf_all.getVariables()->Print("v");

    // Set parameters to constant save model to file
    Jpsi_mean.setConstant(kTRUE);
    Jpsi_devia1.setConstant(kTRUE);
    Jpsi_devia2.setConstant(kTRUE);
    Jpsi_alpha1.setConstant(kTRUE);
    Jpsi_nx1.setConstant(kTRUE);
    Jpsi_ratio.setConstant(kTRUE);

    Jpsi_mu1.setConstant(kTRUE);
    Jpsi_sigma1.setConstant(kTRUE);
    Jpsi_sigma2.setConstant(kTRUE);
    Jpsi_prop1.setConstant(kTRUE);
    Jpsi_sigma3.setConstant(kTRUE);
    Jpsi_coef1.setConstant(kTRUE);
    Jpsi_prop2.setConstant(kTRUE);
    Jpsi_sigma7.setConstant(kTRUE);
    Jpsi_coef3.setConstant(kTRUE);

    RooWorkspace *wsp = new RooWorkspace("wsp", "wsp");
    wsp->import(pdf_all);
    const char *modelFile = useNativeAsymptotic ?
        (useSeed ?
        (isRef ? "Model_4D_tot_ref_native_asymptotic_seeded.root" :
                 "Model_4D_tot_native_asymptotic_seeded.root") :
        (isRef ? "Model_4D_tot_ref.root" : "Model_4D_tot.root")) :
        (useAsymptotic ?
        (isRef ? "Model_4D_tot_ref_asymptotic.root" : "Model_4D_tot_asymptotic.root") :
        (!useSumW2 ?
        (isRef ? "Model_4D_tot_ref_previous_method.root" : "Model_4D_tot_previous_method.root") :
        (isRef ? "Model_4D_tot_ref.root" : "Model_4D_tot.root")));
    wsp->writeToFile(modelFile);
    if(useNativeAsymptotic) {
        const char *fitResultFile = isRef ?
            (useSeed ? "Fit_4D_tot_ref_native_asymptotic_seeded.root" :
                       "Fit_4D_tot_ref_native_asymptotic_unseeded.root") :
            (useSeed ? "Fit_4D_tot_native_asymptotic_seeded.root" :
                       "Fit_4D_tot_native_asymptotic_unseeded.root");
        TFile outputFile(fitResultFile, "RECREATE");
        res->Write("fit_result_native_asymptotic");
        TObjString order(parameterOrder(res->floatParsFinal()).c_str());
        order.Write("parameter_order");
        outputFile.Close();
    } else if(useAsymptotic) {
        const char *fitResultFile = isRef ?
            "Fit_4D_tot_ref_asymptotic.root" : "Fit_4D_tot_asymptotic.root";
        TFile outputFile(fitResultFile, "RECREATE");
        correctedResult->Write("fit_result_asymptotic");
        res->Write("fit_result_weighted_hessian");
        asymptoticCovariance.Write("asymptotic_covariance");
        covarianceHalfStep.Write("asymptotic_covariance_half_step");
        covarianceDoubleStep.Write("asymptotic_covariance_double_step");
        TObjString order(parameterOrder(res->floatParsFinal()).c_str());
        order.Write("parameter_order");
        outputFile.Close();
    }
    // Log message to screen
    cout<<"Error convention: "<<(useNativeAsymptotic ? "native AsymptoticError(true)" :
        (useAsymptotic ? "asymptotic sandwich covariance" :
        (useSumW2 ? "SumW2Error(true)" : "previous SumW2Error(false)")))<<endl;
    cout<<"Minimizer strategy: "<<minimizerStrategy<<endl;
    cout<<"Status: "<<res->status()<<endl;
    cout<<"Event yield: "<<n_P_P.getVal()<<" +/- "<<n_P_P.getError()<<endl;
    if(useAsymptotic) {
        cout<<"n_P_P asymptotic error (half/nominal/double derivative step): "
            <<nPPErrorHalfStep<<" / "<<nPPErrorNominalStep<<" / "<<nPPErrorDoubleStep<<endl;
        cout<<"Maximum relative parameter-error shift (half/double step): "
            <<halfStepMaximumShift<<" / "<<doubleStepMaximumShift<<endl;
        cout<<"Correlation eigenvalue range: "
            <<covarianceDiagnostics.minCorrelationEigenvalue<<" to "
            <<covarianceDiagnostics.maxCorrelationEigenvalue<<endl;
    }
    delete correctedResult;
    delete res;
    delete fitData;
    return;
}
