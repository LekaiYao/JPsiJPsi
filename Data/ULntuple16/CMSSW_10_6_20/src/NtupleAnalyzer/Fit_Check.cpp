#include "Plot_4D.hpp"
#include "RVersion.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TPad.h"
#include "TROOT.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TTree.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooBinning.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RooFit;
using namespace std;

namespace {

const int kMaximumFitAttempts = 3;
const double kBoundarySigmaTolerance = 0.1;

struct ObservedEvent {
    double values[4];
    double weight;
};

struct ProjectionDefinition {
    string name;
    RooRealVar *variable;
    vector<double> edges;
};

struct Chi2Result {
    double chi2;
    double minimumEffectiveEntries;
    vector<double> sumWeights;
    vector<double> sumWeightsSquared;
    vector<double> expected;
};

struct NumericalCdf {
    unique_ptr<RooAbsPdf> marginal;
    RooRealVar *variable;
    vector<double> grid;
    vector<double> cumulative;

    double value(double coordinate) const {
        if(coordinate <= grid.front()) return 0.0;
        if(coordinate >= grid.back()) return 1.0;
        auto upper = std::upper_bound(grid.begin(), grid.end(), coordinate);
        const int highIndex = static_cast<int>(upper - grid.begin());
        const int lowIndex = highIndex - 1;
        const double fraction = (coordinate - grid[lowIndex]) /
            (grid[highIndex] - grid[lowIndex]);
        return cumulative[lowIndex] + fraction *
            (cumulative[highIndex] - cumulative[lowIndex]);
    }
};

bool validResultTag(const string &tag) {
    if(tag.empty()) return false;
    for(char character : tag) {
        if(!(std::isalnum(static_cast<unsigned char>(character)) ||
             character == '_' || character == '-')) return false;
    }
    return true;
}

const RooRealVar *fittedParameter(
    const RooFitResult *result,
    const string &name
) {
    return result ? dynamic_cast<const RooRealVar *>(
        result->floatParsFinal().find(name.c_str())) : nullptr;
}

bool hasFinitePositiveError(
    const RooFitResult *result,
    const string &name
) {
    const RooRealVar *parameter = fittedParameter(result, name);
    return parameter && std::isfinite(parameter->getError()) &&
        parameter->getError() > 0.0;
}

bool acceptedToyFit(const RooFitResult *result) {
    return result && result->status() == 0 && result->covQual() == 3 &&
        std::isfinite(result->edm()) && result->edm() < 0.01 &&
        hasFinitePositiveError(result, "n_P_P");
}

string csvField(string value) {
    string escaped;
    escaped.reserve(value.size() + 2);
    escaped.push_back('"');
    for(char character : value) {
        if(character == '"') escaped.push_back('"');
        if(character == '\n' || character == '\r') escaped.push_back(' ');
        else escaped.push_back(character);
    }
    escaped.push_back('"');
    return escaped;
}

string nearBoundaryParameters(const RooFitResult *result, int &count) {
    count = 0;
    if(!result) return "";
    ostringstream names;
    const RooArgList &parameters = result->floatParsFinal();
    for(int index = 0; index < parameters.getSize(); ++index) {
        const RooRealVar *parameter = dynamic_cast<const RooRealVar *>(
            parameters.at(index));
        if(!parameter || !parameter->hasMin() || !parameter->hasMax()) continue;
        const double range = parameter->getMax() - parameter->getMin();
        if(!std::isfinite(range) || range <= 0.0) continue;
        const double tolerance = std::max(
            1e-9, kBoundarySigmaTolerance * std::fabs(parameter->getError()));
        const bool nearLower =
            std::fabs(parameter->getVal() - parameter->getMin()) <= tolerance;
        const bool nearUpper =
            std::fabs(parameter->getMax() - parameter->getVal()) <= tolerance;
        if(!nearLower && !nearUpper) continue;
        if(count++) names << ";";
        names << parameter->GetName() << (nearLower ? ":lower" : ":upper");
    }
    return names.str();
}

uint32_t toySeed(uint32_t baseSeed, int toyId) {
    uint64_t mixed = static_cast<uint64_t>(baseSeed) +
        104729ULL * static_cast<uint64_t>(toyId + 1);
    mixed ^= mixed >> 16;
    mixed *= 0x7feb352dULL;
    mixed ^= mixed >> 15;
    uint32_t seed = static_cast<uint32_t>(mixed % 900000000ULL);
    return seed ? seed : 1U;
}

vector<ObservedEvent> readObservedEvents(TTree &tree) {
    double mass1 = 0.0;
    double mass2 = 0.0;
    double ctau1 = 0.0;
    double ctau2 = 0.0;
    double weight = 0.0;
    tree.SetBranchAddress("Jpsi_mass1", &mass1);
    tree.SetBranchAddress("Jpsi_mass2", &mass2);
    tree.SetBranchAddress("Jpsi_ctau1", &ctau1);
    tree.SetBranchAddress("Jpsi_ctau2", &ctau2);
    tree.SetBranchAddress("evt_weight", &weight);

    vector<ObservedEvent> events;
    for(Long64_t entry = 0; entry < tree.GetEntries(); ++entry) {
        tree.GetEntry(entry);
        if(mass1 < 2.95 || mass1 > 3.25 ||
           mass2 < 2.95 || mass2 > 3.25 ||
           ctau1 < -0.03 || ctau1 > 0.16 ||
           ctau2 < -0.03 || ctau2 > 0.16) continue;
        if(!std::isfinite(weight) || weight <= 0.0) {
            throw runtime_error("Observed fit-range data contain an invalid weight");
        }
        ObservedEvent event = {{mass1, mass2, ctau1, ctau2}, weight};
        events.push_back(event);
    }
    if(events.empty()) throw runtime_error("No observed events are inside the 4D fit range");
    return events;
}

void resetFloatingParameters(
    RooWorkspace &workspace,
    const RooFitResult &truthResult,
    bool makeConstant
) {
    const RooArgList &truthParameters = truthResult.floatParsFinal();
    for(int index = 0; index < truthParameters.getSize(); ++index) {
        const RooRealVar *truth = dynamic_cast<const RooRealVar *>(
            truthParameters.at(index));
        RooRealVar *parameter = truth ? workspace.var(truth->GetName()) : nullptr;
        if(!truth || !parameter) {
            throw runtime_error("Cannot map a nominal floating parameter into the workspace");
        }
        parameter->setVal(truth->getVal());
        parameter->setError(truth->getError());
        parameter->setConstant(makeConstant);
    }
}

void resetYieldOnlyToyFit(
    RooWorkspace &workspace,
    const RooFitResult &truthResult
) {
    resetFloatingParameters(workspace, truthResult, true);
    const vector<string> yieldNames = {
        "n_P_P", "n_P_NP", "n_NP_NP", "n_Sig_Comb", "n_Comb_Comb"
    };
    for(const string &name : yieldNames) {
        RooRealVar *yield = workspace.var(name.c_str());
        if(!yield) {
            throw runtime_error("Cannot find a required yield in the nominal workspace");
        }
        if(name == "n_Comb_Comb") yield->setMin(0.0);
        yield->setConstant(false);
    }
}

NumericalCdf numericalCdf(
    RooAbsPdf &pdf,
    RooRealVar &variable,
    const vector<RooRealVar *> &observables,
    int gridIntervals=4096
) {
    RooArgSet integratedVariables;
    for(RooRealVar *observable : observables) {
        if(observable != &variable) integratedVariables.add(*observable);
    }
    unique_ptr<RooAbsPdf> marginal(pdf.createProjection(integratedVariables));
    if(!marginal) throw runtime_error("Cannot construct a one-dimensional PDF projection");
    const double previous = variable.getVal();
    RooArgSet normalization(variable);
    vector<double> grid(gridIntervals + 1, 0.0);
    vector<double> density(gridIntervals + 1, 0.0);
    vector<double> cumulative(gridIntervals + 1, 0.0);
    const double step = (variable.getMax() - variable.getMin()) / gridIntervals;
    for(int index = 0; index <= gridIntervals; ++index) {
        grid[index] = variable.getMin() + index * step;
        variable.setVal(grid[index]);
        density[index] = marginal->getVal(&normalization);
        if(!std::isfinite(density[index]) || density[index] < 0.0) {
            throw runtime_error("A numerical marginal-PDF integration point is invalid");
        }
        if(index > 0) {
            cumulative[index] = cumulative[index - 1] +
                0.5 * step * (density[index - 1] + density[index]);
        }
    }
    variable.setVal(previous);
    const double integral = cumulative.back();
    if(!std::isfinite(integral) || integral <= 0.0) {
        throw runtime_error("Numerical marginal-PDF normalization failed");
    }
    for(double &value : cumulative) value /= integral;
    return {std::move(marginal), &variable, std::move(grid), std::move(cumulative)};
}

vector<double> quantileEdges(
    RooAbsPdf &pdf,
    RooRealVar &variable,
    const vector<RooRealVar *> &observables,
    int numberOfBins
) {
    NumericalCdf cdf = numericalCdf(pdf, variable, observables);
    vector<double> edges(numberOfBins + 1, 0.0);
    edges.front() = variable.getMin();
    edges.back() = variable.getMax();
    for(int edgeIndex = 1; edgeIndex < numberOfBins; ++edgeIndex) {
        const double target = static_cast<double>(edgeIndex) / numberOfBins;
        double low = variable.getMin();
        double high = variable.getMax();
        for(int iteration = 0; iteration < 80; ++iteration) {
            const double middle = 0.5 * (low + high);
            if(cdf.value(middle) < target) low = middle;
            else high = middle;
        }
        edges[edgeIndex] = 0.5 * (low + high);
        if(!(edges[edgeIndex] > edges[edgeIndex - 1])) {
            throw runtime_error("Nominal PDF quantile edges are not strictly increasing");
        }
    }
    return edges;
}

int binIndex(double value, const vector<double> &edges) {
    if(value < edges.front() || value > edges.back()) return -1;
    if(value == edges.back()) return static_cast<int>(edges.size()) - 2;
    auto upper = std::upper_bound(edges.begin(), edges.end(), value);
    const int index = static_cast<int>(upper - edges.begin()) - 1;
    return index >= 0 && index + 1 < static_cast<int>(edges.size()) ? index : -1;
}

Chi2Result projectionChi2(
    RooAbsPdf &pdf,
    RooRealVar &variable,
    const vector<RooRealVar *> &observables,
    const RooArgSet &normalizationSet,
    const vector<double> &edges,
    const vector<ObservedEvent> *observedEvents,
    const RooAbsData *toyData,
    int valueIndex
) {
    const int numberOfBins = static_cast<int>(edges.size()) - 1;
    Chi2Result result;
    result.chi2 = 0.0;
    result.minimumEffectiveEntries = std::numeric_limits<double>::infinity();
    result.sumWeights.assign(numberOfBins, 0.0);
    result.sumWeightsSquared.assign(numberOfBins, 0.0);
    result.expected.assign(numberOfBins, 0.0);

    if(observedEvents) {
        for(const ObservedEvent &event : *observedEvents) {
            const int index = binIndex(event.values[valueIndex], edges);
            if(index < 0) continue;
            result.sumWeights[index] += event.weight;
            result.sumWeightsSquared[index] += event.weight * event.weight;
        }
    } else if(toyData) {
        for(int entry = 0; entry < toyData->numEntries(); ++entry) {
            const RooArgSet *row = toyData->get(entry);
            const double value = row->getRealValue(variable.GetName());
            const double weight = toyData->weight();
            const int index = binIndex(value, edges);
            if(index < 0) continue;
            result.sumWeights[index] += weight;
            result.sumWeightsSquared[index] += weight * weight;
        }
    } else {
        throw runtime_error("Projection chi2 received neither observed nor toy data");
    }

    NumericalCdf cdf = numericalCdf(pdf, variable, observables);
    const double expectedTotal = pdf.expectedEvents(&normalizationSet);
    if(!std::isfinite(expectedTotal) || expectedTotal <= 0.0) {
        throw runtime_error("The fitted model has an invalid expected event yield");
    }
    for(int index = 0; index < numberOfBins; ++index) {
        const double probability = cdf.value(edges[index + 1]) -
            cdf.value(edges[index]);
        result.expected[index] = expectedTotal * probability;
        if(result.sumWeightsSquared[index] <= 0.0 || probability <= 0.0) {
            throw runtime_error("A nominal-quantile chi2 bin has zero variance or expectation");
        }
        const double difference = result.sumWeights[index] - result.expected[index];
        result.chi2 += difference * difference / result.sumWeightsSquared[index];
        const double effectiveEntries = result.sumWeights[index] * result.sumWeights[index] /
            result.sumWeightsSquared[index];
        result.minimumEffectiveEntries = std::min(
            result.minimumEffectiveEntries, effectiveEntries);
    }
    return result;
}

unique_ptr<RooDataSet> makeObservedPlotData(
    const vector<ObservedEvent> &events,
    const vector<RooRealVar *> &observables
) {
    RooRealVar eventWeight("evt_weight", "evt_weight", 0.0, 1000.0);
    RooArgSet plotVariables;
    for(RooRealVar *observable : observables) plotVariables.add(*observable);
    plotVariables.add(eventWeight);
    unique_ptr<RooDataSet> data(new RooDataSet(
        "observed_plot_data", "observed_plot_data", plotVariables,
        WeightVar(eventWeight)));
    for(const ObservedEvent &event : events) {
        for(int index = 0; index < 4; ++index) {
            observables[index]->setVal(event.values[index]);
        }
        eventWeight.setVal(event.weight);
        data->add(plotVariables, event.weight);
    }
    if(data->numEntries() != static_cast<int>(events.size())) {
        throw runtime_error("Observed plotting dataset has an unexpected size");
    }
    return data;
}

void drawProjectionPlot(
    const string &plotDirectory,
    const RooAbsData &data,
    const RooAddPdf &pdf,
    const RooAbsPdf &pdfPP,
    const RooAbsPdf &pdfPNP,
    const RooAbsPdf &pdfNPP,
    const RooAbsPdf &pdfNPNP,
    const RooAbsPdf &pdfSigComb,
    const RooAbsPdf &pdfCombSig,
    const RooAbsPdf &pdfCombComb,
    const ProjectionDefinition &projection,
    const Chi2Result &chi2,
    const string &axisTitle,
    bool logarithmicY
) {
    const int numberOfBins = static_cast<int>(projection.edges.size()) - 1;
    const string objectSuffix = projection.name + "_" + std::to_string(numberOfBins);
    const string binningName = "chi2_binning_" + objectSuffix;
    RooBinning binning(numberOfBins, projection.edges.data(), binningName.c_str());
    unique_ptr<RooPlot> frame(projection.variable->frame(Bins(numberOfBins)));
    data.plotOn(frame.get(), Binning(binning), DataError(RooAbsData::SumW2),
        MarkerStyle(20), MarkerSize(1.0), Name("Data"));
    plot_on(frame.get(), pdf, pdfPP, pdfPNP, pdfNPP, pdfNPNP,
        pdfSigComb, pdfCombSig, pdfCombComb, "");
    data.plotOn(frame.get(), Binning(binning), DataError(RooAbsData::SumW2),
        MarkerStyle(20), MarkerSize(1.0), Name("DataTop"));

    TCanvas canvas(("canvas_" + objectSuffix).c_str(),
        ("canvas_" + objectSuffix).c_str(), 1500, 1500);
    canvas.Divide(1, 2);
    TPad *upperPad = dynamic_cast<TPad *>(canvas.cd(1));
    upperPad->SetPad(0.01, 0.20, 0.99, 0.99);
    upperPad->SetLeftMargin(0.15);
    upperPad->SetTicks(1, 1);
    frame->SetTitle("");
    if(logarithmicY) {
        frame->SetMinimum(0.5);
        frame->SetMaximum(8.0 * frame->GetMaximum());
        upperPad->SetLogy();
    } else {
        frame->SetMaximum(1.25 * frame->GetMaximum());
    }
    frame->Draw();
    TLegend legend(.65, .50, .90, .90);
    add_entry(frame.get(), &legend);
    legend.Draw();
    TLatex experimentLabel;
    experimentLabel.SetNDC();
    add_latex(experimentLabel);
    TLatex chi2Label;
    chi2Label.SetNDC();
    chi2Label.SetTextFont(42);
    chi2Label.SetTextSize(0.032);
    ostringstream chi2Text;
    chi2Text << fixed << setprecision(2) << "#chi^{2}/ndf = ";
    if(logarithmicY) {
        chi2Text << chi2.chi2 / (numberOfBins - 1);
    } else {
        chi2Text << chi2.chi2 << "/" << numberOfBins - 1 << " = "
            << chi2.chi2 / (numberOfBins - 1);
    }
    const double annotationX = logarithmicY ? 0.42 : 0.18;
    const double annotationY = 0.78;
    chi2Label.DrawLatex(annotationX, annotationY, chi2Text.str().c_str());
    ostringstream binText;
    binText << numberOfBins << (logarithmicY ? " quantile bins" :
        " nominal-PDF quantile bins");
    chi2Label.DrawLatex(annotationX, annotationY - 0.045, binText.str().c_str());

    TPad *pullPad = dynamic_cast<TPad *>(canvas.cd(2));
    pullPad->SetPad(0.01, 0.03, 0.99, 0.25);
    pullPad->SetLeftMargin(0.15);
    pullPad->SetBottomMargin(0.25);
    pullPad->SetTicks(1, 1);
    TH1D pullAxis(("pull_axis_" + objectSuffix).c_str(), "",
        numberOfBins, projection.edges.data());
    pullAxis.SetDirectory(nullptr);
    pullAxis.SetStats(false);
    pullAxis.GetYaxis()->SetTitle("Pull");
    pullAxis.GetYaxis()->SetNdivisions(505);
    pullAxis.GetYaxis()->SetTitleSize(0.10);
    pullAxis.GetYaxis()->SetTitleOffset(0.55);
    pullAxis.GetYaxis()->SetLabelSize(0.10);
    pullAxis.GetXaxis()->SetTitle(axisTitle.c_str());
    pullAxis.GetXaxis()->SetTitleSize(0.10);
    pullAxis.GetXaxis()->SetLabelSize(0.10);
    double maximumAbsolutePull = 0.0;
    TGraphErrors pullGraph(numberOfBins);
    pullGraph.SetName(("Pull_" + objectSuffix).c_str());
    pullGraph.SetMarkerStyle(20);
    pullGraph.SetMarkerSize(1.0);
    pullGraph.SetLineWidth(1);
    for(int index = 0; index < numberOfBins; ++index) {
        const double pull = (chi2.sumWeights[index] - chi2.expected[index]) /
            std::sqrt(chi2.sumWeightsSquared[index]);
        maximumAbsolutePull = std::max(maximumAbsolutePull, std::fabs(pull));
        pullGraph.SetPoint(index,
            0.5 * (projection.edges[index] + projection.edges[index + 1]), pull);
        pullGraph.SetPointError(index, 0.0, 0.0);
    }
    const double pullRange = std::max(5.0, std::ceil(maximumAbsolutePull + 0.5));
    pullAxis.SetMinimum(-pullRange);
    pullAxis.SetMaximum(pullRange);
    pullAxis.Draw("AXIS");
    TLine zeroLine(projection.edges.front(), 0.0, projection.edges.back(), 0.0);
    set_line_style(&zeroLine);
    zeroLine.Draw("same");
    pullGraph.Draw("P same");

    const string outputBase = plotDirectory + "/projection_" + projection.name;
    canvas.SaveAs((outputBase + ".pdf").c_str());
    canvas.SaveAs((outputBase + ".png").c_str());
}

void drawToyFitProjectionPlot(
    const string &plotDirectory,
    const RooAbsData &data,
    const RooAddPdf &pdf,
    const RooAbsPdf &pdfPP,
    const RooAbsPdf &pdfPNP,
    const RooAbsPdf &pdfNPP,
    const RooAbsPdf &pdfNPNP,
    const RooAbsPdf &pdfSigComb,
    const RooAbsPdf &pdfCombSig,
    const RooAbsPdf &pdfCombComb,
    RooRealVar &variable,
    const string &outputName,
    const string &axisTitle,
    bool logarithmicY,
    int toyId,
    uint32_t seed,
    const RooFitResult &fitResult
) {
    const int numberOfBins = 100;
    const string objectSuffix = outputName + "_" + std::to_string(toyId);
    unique_ptr<RooPlot> frame(variable.frame(Bins(numberOfBins)));
    unique_ptr<RooPlot> pullFrame(variable.frame(Bins(numberOfBins)));
    data.plotOn(frame.get(), DataError(RooAbsData::SumW2),
        MarkerStyle(20), MarkerSize(1.0), Name("Data"));
    plot_on(frame.get(), pdf, pdfPP, pdfPNP, pdfNPP, pdfNPNP,
        pdfSigComb, pdfCombSig, pdfCombComb, "");
    data.plotOn(frame.get(), DataError(RooAbsData::SumW2),
        MarkerStyle(20), MarkerSize(1.0), Name("DataTop"));
    RooHist *pull = frame->pullHist("Data", "All");
    if(!pull) throw runtime_error("Cannot construct a toy-fit projection pull");
    pullFrame->addPlotable(pull, "P");

    TCanvas canvas(("toy_canvas_" + objectSuffix).c_str(),
        ("toy_canvas_" + objectSuffix).c_str(), 1500, 1500);
    canvas.Divide(1, 2);
    TPad *upperPad = dynamic_cast<TPad *>(canvas.cd(1));
    upperPad->SetPad(0.01, 0.20, 0.99, 0.99);
    upperPad->SetLeftMargin(0.15);
    upperPad->SetTicks(1, 1);
    frame->SetTitle("");
    if(logarithmicY) {
        frame->SetMinimum(1.0);
        frame->SetMaximum(5.0 * frame->GetMaximum());
        upperPad->SetLogy();
    } else {
        frame->SetMaximum(1.25 * frame->GetMaximum());
    }
    frame->Draw();

    TLegend legend(.65, .50, .90, .90);
    const string dataLabel = "Weighted toy " + std::to_string(toyId);
    legend.AddEntry(frame->findObject("Data"), dataLabel.c_str(), "LEP");
    legend.AddEntry(frame->findObject("All"), "Total p.d.f.", "L");
    legend.AddEntry(frame->findObject("P_P"), "prompt, prompt", "L");
    legend.AddEntry(frame->findObject("P_NP"), "prompt, non-prompt", "L");
    legend.AddEntry(frame->findObject("NP_P"), "non-prompt, prompt", "L");
    legend.AddEntry(frame->findObject("NP_NP"), "non-prompt, non-prompt", "L");
    legend.AddEntry(frame->findObject("Sig_Comb"), "J/#psi, #mu^{+}#mu^{-}", "L");
    legend.AddEntry(frame->findObject("Comb_Sig"), "#mu^{+}#mu^{-}, J/#psi", "L");
    legend.AddEntry(frame->findObject("Comb_Comb"),
        "#mu^{+}#mu^{-}, #mu^{+}#mu^{-}", "L");
    legend.Draw();

    TLatex label;
    label.SetNDC();
    label.SetTextFont(62);
    label.SetTextSize(0.05);
    label.DrawLatex(0.18, 0.85, "CMS");
    label.SetTextFont(52);
    label.SetTextSize(0.04);
    const string toyLabel = "Toy " + std::to_string(toyId);
    label.DrawLatex(0.27, 0.85, toyLabel.c_str());
    label.SetTextFont(42);
    label.SetTextSize(0.024);
    ostringstream seedText;
    seedText << "seed " << seed;
    label.DrawLatex(0.43, 0.83, seedText.str().c_str());
    ostringstream statusText;
    statusText << "st=" << fitResult.status() << ", EDM="
        << scientific << setprecision(1) << fitResult.edm();
    label.DrawLatex(0.43, 0.79, statusText.str().c_str());

    TPad *pullPad = dynamic_cast<TPad *>(canvas.cd(2));
    pullPad->SetPad(0.01, 0.03, 0.99, 0.25);
    pullPad->SetLeftMargin(0.15);
    pullPad->SetBottomMargin(0.25);
    pullPad->SetTicks(1, 1);
    set_pull_style(pullFrame.get());
    pullFrame->GetXaxis()->SetTitle(axisTitle.c_str());
    double maximumAbsolutePull = 0.0;
    for(int index = 0; index < pull->GetN(); ++index) {
        maximumAbsolutePull = std::max(
            maximumAbsolutePull, std::fabs(pull->GetY()[index]));
    }
    const double pullRange = std::max(5.0, std::ceil(maximumAbsolutePull + 0.5));
    pullFrame->SetMinimum(-pullRange);
    pullFrame->SetMaximum(pullRange);
    pullFrame->Draw();
    TLine zeroLine(variable.getMin(), 0.0, variable.getMax(), 0.0);
    set_line_style(&zeroLine);
    zeroLine.Draw("same");

    const string outputBase = plotDirectory + "/" + outputName;
    canvas.SaveAs((outputBase + ".pdf").c_str());
    canvas.SaveAs((outputBase + ".png").c_str());
}

int existingToyRows(const string &path) {
    ifstream input(path);
    if(!input) return -1;
    string line;
    int lines = -1;
    while(std::getline(input, line)) {
        if(!line.empty()) ++lines;
    }
    return std::max(lines, 0);
}

} // namespace

void Fit_Check(
    int startToy=0,
    int numberOfToys=1,
    string resultTag="gof_root640_v1",
    unsigned int baseSeed=20260901,
    int projectionBins=20,
    bool makeToy0Diagnostics=true
) {
#if ROOT_VERSION_CODE < ROOT_VERSION(6, 40, 0)
    cerr << "The corrected-error toy validation requires ROOT >= 6.40" << endl;
    gSystem->Exit(2);
    return;
#endif
    if(startToy < 0 || numberOfToys < 0 ||
       (numberOfToys == 0 && startToy != 0) || projectionBins < 2 ||
       !validResultTag(resultTag)) {
        cerr << "Invalid toy range, result tag, or projection bin count" << endl;
        gSystem->Exit(2);
        return;
    }

    const string resultDirectory = "fit_results/" + resultTag;
    const string resultsPath = resultDirectory + "/toy_results.csv";
    const string fitAttemptsPath = resultDirectory + "/fit_attempts.csv";
    gSystem->mkdir(resultDirectory.c_str(), kTRUE);
    try {
        TFile dataFile("WeightData.root", "READ");
        TTree *dataTree = dynamic_cast<TTree *>(dataFile.Get("data"));
        if(dataFile.IsZombie() || !dataTree) {
            throw runtime_error("Cannot read WeightData.root:data");
        }
        const vector<ObservedEvent> observedEvents = readObservedEvents(*dataTree);
        double observedSumWeights = 0.0;
        double observedSumWeightsSquared = 0.0;
        for(const ObservedEvent &event : observedEvents) {
            observedSumWeights += event.weight;
            observedSumWeightsSquared += event.weight * event.weight;
        }

        TFile modelFile("Model_4D_tot.root", "READ");
        RooWorkspace *workspace = dynamic_cast<RooWorkspace *>(modelFile.Get("wsp"));
        RooAddPdf *pdf = workspace ? dynamic_cast<RooAddPdf *>(workspace->pdf("pdf_all")) : nullptr;
        RooAbsPdf *pdfPP = workspace ? workspace->pdf("pdf_P_P") : nullptr;
        RooAbsPdf *pdfPNP = workspace ? workspace->pdf("pdf_P_NP") : nullptr;
        RooAbsPdf *pdfNPP = workspace ? workspace->pdf("pdf_NP_P") : nullptr;
        RooAbsPdf *pdfNPNP = workspace ? workspace->pdf("pdf_NP_NP") : nullptr;
        RooAbsPdf *pdfSigComb = workspace ? workspace->pdf("pdf_Sig_Comb") : nullptr;
        RooAbsPdf *pdfCombSig = workspace ? workspace->pdf("pdf_Comb_Sig") : nullptr;
        RooAbsPdf *pdfCombComb = workspace ? workspace->pdf("pdf_Comb_Comb") : nullptr;
        if(modelFile.IsZombie() || !workspace || !pdf || !pdfPP || !pdfPNP ||
           !pdfNPP || !pdfNPNP || !pdfSigComb || !pdfCombSig || !pdfCombComb) {
            throw runtime_error("Cannot read the nominal 4D model and all components");
        }
        TFile fitResultFile("Fit_4D_tot_native_asymptotic_unseeded.root", "READ");
        RooFitResult *truthResult = dynamic_cast<RooFitResult *>(
            fitResultFile.Get("fit_result_native_asymptotic"));
        if(fitResultFile.IsZombie() || !acceptedToyFit(truthResult)) {
            throw runtime_error("The nominal ROOT 6.40 fit result is missing or invalid");
        }

        RooRealVar *mass1 = workspace->var("Jpsi_mass1");
        RooRealVar *mass2 = workspace->var("Jpsi_mass2");
        RooRealVar *ctau1 = workspace->var("Jpsi_ctau1");
        RooRealVar *ctau2 = workspace->var("Jpsi_ctau2");
        RooRealVar *nPP = workspace->var("n_P_P");
        if(!mass1 || !mass2 || !ctau1 || !ctau2 || !nPP) {
            throw runtime_error("The nominal workspace is missing required variables");
        }
        vector<RooRealVar *> observables = {mass1, mass2, ctau1, ctau2};
        RooArgSet observableSet;
        for(RooRealVar *observable : observables) observableSet.add(*observable);
        resetFloatingParameters(*workspace, *truthResult, true);
        const RooRealVar *truthNPP = dynamic_cast<const RooRealVar *>(
            truthResult->floatParsFinal().find("n_P_P"));
        if(!truthNPP) throw runtime_error("n_P_P is absent from the nominal fit result");
        const double nPPTruth = truthNPP->getVal();
        const double expectedTruth = pdf->expectedEvents(&observableSet);
        const double poissonRate = expectedTruth / observedSumWeights;
        if(!std::isfinite(poissonRate) || poissonRate <= 0.0 ||
           std::fabs(poissonRate - 1.0) > 0.01) {
            throw runtime_error("Nominal expectedEvents and observed sumw differ by more than 1%");
        }

        vector<ProjectionDefinition> projections = {
            {"Jpsi_mass1", mass1, {}},
            {"Jpsi_mass2", mass2, {}},
            {"Jpsi_ctau1", ctau1, {}},
            {"Jpsi_ctau2", ctau2, {}}
        };
        for(ProjectionDefinition &projection : projections) {
            projection.edges = quantileEdges(
                *pdf, *projection.variable, observables, projectionBins);
        }

        vector<Chi2Result> observedChi2;
        for(int index = 0; index < 4; ++index) {
            observedChi2.push_back(projectionChi2(
                *pdf, *projections[index].variable, observables, observableSet,
                projections[index].edges, &observedEvents, nullptr, index));
        }

        if(startToy == 0) {
            if(!gSystem->AccessPathName(resultsPath.c_str())) {
                throw runtime_error("Refusing to overwrite an existing toy_results.csv");
            }
            ofstream configuration(resultDirectory + "/configuration.txt");
            configuration << setprecision(17)
                << "root_version=" << gROOT->GetVersion() << "\n"
                << "error_convention=ROOT native AsymptoticError(true)\n"
                << "external_seed=false; toy fits start from nominal truth as in the historical toy method\n"
                << "fit_parameters=n_P_P,n_P_NP,n_NP_NP,n_Sig_Comb,n_Comb_Comb floating; all nominal shape parameters fixed\n"
                << "n_Comb_Comb_range=physical lower boundary 0; upper boundary inherited from the nominal workspace\n"
                << "fit_gate=status=0,covQual=3,EDM<0.01\n"
                << "fit_error_gate=n_P_P error finite and positive\n"
                << "fit_attempt_1=baseline yield-only configuration\n"
                << "fit_attempt_2=continue from attempt 1 endpoint with Strategy(2) and Offset(true); use n_Comb_Comb fixed at 0 if attempt 1 failed at its lower boundary\n"
                << "fit_attempt_3=only if attempt 2 first fails at the n_Comb_Comb lower boundary; refit with n_Comb_Comb fixed at 0\n"
                << "n_Comb_Comb_boundary_fallback=after a failed fit only, fix n_Comb_Comb=0 when the result is at the lower boundary or ROOT reports its out-of-range covariance evaluation\n"
                << "near_boundary_definition=distance <= max(1e-9,0.1*fitted_error)\n"
                << "observed_fit_entries=" << observedEvents.size() << "\n"
                << "observed_sumw=" << observedSumWeights << "\n"
                << "observed_sumw2=" << observedSumWeightsSquared << "\n"
                << "nominal_expected_events=" << expectedTruth << "\n"
                << "poisson_weight_rate=" << poissonRate << "\n"
                << "n_P_P_truth=" << nPPTruth << "\n"
                << "projection_bins=" << projectionBins << " nominal-PDF quantiles\n"
                << "projection_integration=numerical trapezoid of the nominal marginal PDF on 4096 intervals; ROOT analytical ctau CDF is not used\n"
                << "chi2_definition=sum((sumw-expected)^2/sumw2); no binned refit\n"
                << "projection_plot=weighted data and all nominal 4D PDF components in the exact chi2 quantile bins\n"
                << "projection_pull=(sumw-expected)/sqrt(sumw2), identical to the chi2 bin contribution\n"
                << "toy0_fit_diagnostics=" << (makeToy0Diagnostics ?
                    "enabled: 100 uniform-bin nominal-style projections plus ROOT artifact" :
                    "disabled: no toy-0 projections or ROOT artifact") << "\n"
                << "base_seed=" << baseSeed << "\n";
            configuration.close();

            ofstream binsFile(resultDirectory + "/projection_bins.csv");
            binsFile << "projection,bin,low,high,data_sumw,data_sumw2,effective_entries,expected\n";
            binsFile << setprecision(17);
            for(int projectionIndex = 0; projectionIndex < 4; ++projectionIndex) {
                for(int bin = 0; bin < projectionBins; ++bin) {
                    const Chi2Result &chi2 = observedChi2[projectionIndex];
                    const double effectiveEntries = chi2.sumWeights[bin] * chi2.sumWeights[bin] /
                        chi2.sumWeightsSquared[bin];
                    binsFile << projections[projectionIndex].name << "," << bin << ","
                        << projections[projectionIndex].edges[bin] << ","
                        << projections[projectionIndex].edges[bin + 1] << ","
                        << chi2.sumWeights[bin] << ","
                        << chi2.sumWeightsSquared[bin] << ","
                        << effectiveEntries << ","
                        << chi2.expected[bin] << "\n";
                }
            }
            binsFile.close();

            ofstream observedFile(resultDirectory + "/observed_chi2.csv");
            observedFile << "projection,chi2,bins,minimum_effective_entries\n";
            observedFile << setprecision(17);
            for(int index = 0; index < 4; ++index) {
                observedFile << projections[index].name << ","
                    << observedChi2[index].chi2 << "," << projectionBins << ","
                    << observedChi2[index].minimumEffectiveEntries << "\n";
            }
            observedFile.close();

            const string plotDirectory = resultDirectory + "/plots";
            gSystem->mkdir(plotDirectory.c_str(), kTRUE);
            unique_ptr<RooDataSet> observedPlotData = makeObservedPlotData(
                observedEvents, observables);
            const vector<string> axisTitles = {
                "M(J/#psi_{1}) [GeV]", "M(J/#psi_{2}) [GeV]",
                "c#tau(J/#psi_{1}) [cm]", "c#tau(J/#psi_{2}) [cm]"
            };
            for(int index = 0; index < 4; ++index) {
                drawProjectionPlot(
                    plotDirectory, *observedPlotData, *pdf, *pdfPP, *pdfPNP,
                    *pdfNPP, *pdfNPNP, *pdfSigComb, *pdfCombSig, *pdfCombComb,
                    projections[index], observedChi2[index], axisTitles[index],
                    index >= 2);
            }

            ofstream results(resultsPath);
            results << "toy_id,seed,n_raw,sumw,sumw2,attempts,accepted,status,covQual,edm,n_P_P_true,n_P_P_fit,n_P_P_error,pull,minNll,chi2_mass1,chi2_mass2,chi2_ctau1,chi2_ctau2\n";
            results.close();

            ofstream fitAttempts(fitAttemptsPath);
            fitAttempts << "toy_id,seed,attempt,mode,status,covQual,edm,minNll,n_P_P_fit,n_P_P_error,near_boundary_count,near_boundary_parameters,exception\n";
            fitAttempts.close();
        } else {
            const int rows = existingToyRows(resultsPath);
            if(rows != startToy) {
                ostringstream message;
                message << "toy_results.csv has " << rows
                    << " rows but the requested startToy is " << startToy;
                throw runtime_error(message.str());
            }
        }

        ofstream results(resultsPath, ios::app);
        if(!results) throw runtime_error("Cannot append toy_results.csv");
        results << setprecision(17);
        ofstream fitAttempts(fitAttemptsPath, ios::app);
        if(!fitAttempts) throw runtime_error("Cannot append fit_attempts.csv");
        fitAttempts << setprecision(17);

        for(int toyId = startToy; toyId < startToy + numberOfToys; ++toyId) {
            const uint32_t seed = toySeed(baseSeed, toyId);
            TRandom3 poissonRandom(seed ^ 0x5bd1e995U);
            vector<double> toyWeights;
            toyWeights.reserve(observedEvents.size() + 256);
            double toySumWeights = 0.0;
            double toySumWeightsSquared = 0.0;
            for(const ObservedEvent &event : observedEvents) {
                const int multiplicity = poissonRandom.Poisson(poissonRate);
                for(int copy = 0; copy < multiplicity; ++copy) {
                    toyWeights.push_back(event.weight);
                    toySumWeights += event.weight;
                    toySumWeightsSquared += event.weight * event.weight;
                }
            }
            if(toyWeights.empty()) throw runtime_error("A toy has zero Poisson-resampled events");
            std::mt19937 shuffleRandom(seed ^ 0x9e3779b9U);
            std::shuffle(toyWeights.begin(), toyWeights.end(), shuffleRandom);

            resetFloatingParameters(*workspace, *truthResult, true);
            RooRandom::randomGenerator()->SetSeed(seed);
            unique_ptr<RooDataSet> generated(pdf->generate(
                observableSet, NumEvents(static_cast<int>(toyWeights.size()))));
            if(!generated || generated->numEntries() != static_cast<int>(toyWeights.size())) {
                throw runtime_error("Nominal PDF generation returned an unexpected event count");
            }

            TTree toyTree("toy_tree", "toy_tree");
            double toyMass1 = 0.0;
            double toyMass2 = 0.0;
            double toyCtau1 = 0.0;
            double toyCtau2 = 0.0;
            double toyWeight = 0.0;
            toyTree.Branch("Jpsi_mass1", &toyMass1);
            toyTree.Branch("Jpsi_mass2", &toyMass2);
            toyTree.Branch("Jpsi_ctau1", &toyCtau1);
            toyTree.Branch("Jpsi_ctau2", &toyCtau2);
            toyTree.Branch("evt_weight", &toyWeight);
            for(int entry = 0; entry < generated->numEntries(); ++entry) {
                const RooArgSet *row = generated->get(entry);
                toyMass1 = row->getRealValue("Jpsi_mass1");
                toyMass2 = row->getRealValue("Jpsi_mass2");
                toyCtau1 = row->getRealValue("Jpsi_ctau1");
                toyCtau2 = row->getRealValue("Jpsi_ctau2");
                toyWeight = toyWeights[entry];
                toyTree.Fill();
            }

            RooRealVar fitWeight("evt_weight", "evt_weight", 0.0, 1000.0);
            RooArgSet toyVariables;
            toyVariables.add(observableSet);
            toyVariables.add(fitWeight);
            RooDataSet toyData(
                "toy_data", "toy_data", toyVariables, Import(toyTree),
                WeightVar("evt_weight"));

            resetYieldOnlyToyFit(*workspace, *truthResult);
            RooRealVar *toyNCombComb = workspace->var("n_Comb_Comb");
            if(!toyNCombComb) {
                throw runtime_error("Cannot find n_Comb_Comb for the toy boundary fallback");
            }
            const auto atCombCombLowerBoundary = [](const RooRealVar *parameter) {
                if(!parameter || !parameter->hasMin()) return false;
                const double tolerance = 1e-4 * std::max(
                    1.0, std::fabs(parameter->getError()));
                return std::fabs(parameter->getMin()) <= tolerance &&
                    parameter->getVal() <= parameter->getMin() + tolerance;
            };
            bool combCombBoundaryFallbackActive = false;
            bool combCombBoundaryFallbackUsed = false;
            const auto activateCombCombBoundaryFallback = [&]() {
                if(combCombBoundaryFallbackActive) return;
                toyNCombComb->setVal(0.0);
                toyNCombComb->setConstant(true);
                combCombBoundaryFallbackActive = true;
                cerr << "Toy " << toyId
                    << ": fixing n_Comb_Comb=0 after a failed boundary fit"
                    << endl;
            };
            unique_ptr<RooFitResult> fitResult;
            int attempts = 0;
            string lastException;
            for(int attempt = 0; attempt < kMaximumFitAttempts; ++attempt) {
                if(attempt >= 2 &&
                   (!combCombBoundaryFallbackActive ||
                    combCombBoundaryFallbackUsed)) break;
                attempts = attempt + 1;
                const bool usingCombCombBoundaryFallback =
                    combCombBoundaryFallbackActive;
                if(usingCombCombBoundaryFallback) {
                    combCombBoundaryFallbackUsed = true;
                }
                const string fitMode = attempt == 0 ?
                    "baseline_yield_only" :
                    (usingCombCombBoundaryFallback ?
                        "strategy2_offset_n_Comb_Comb_fixed_zero" :
                        "strategy2_offset_yield_only");
                string attemptException;
                unique_ptr<RooFitResult> attemptResult;
                try {
                    if(attempt == 0) {
                        attemptResult.reset(pdf->fitTo(
                            toyData, Save(), Extended(kTRUE), AsymptoticError(kTRUE),
                            PrintLevel(-1), Warnings(kFALSE), Verbose(kFALSE)));
                    } else {
                        attemptResult.reset(pdf->fitTo(
                            toyData, Save(), Extended(kTRUE), AsymptoticError(kTRUE),
                            Strategy(2), Offset(true), PrintLevel(-1),
                            Warnings(kFALSE), Verbose(kFALSE)));
                    }
                } catch(const std::exception &error) {
                    attemptException = error.what();
                    lastException = attemptException;
                    cerr << "Toy " << toyId << " attempt " << attempts
                        << " threw: " << attemptException << endl;
                    const bool explicitCombCombRangeException =
                        attemptException.find("outside the default range") !=
                            string::npos &&
                        attemptException.find("n_Comb_Comb") != string::npos;
                    if(!combCombBoundaryFallbackActive &&
                       (explicitCombCombRangeException ||
                        atCombCombLowerBoundary(toyNCombComb))) {
                        activateCombCombBoundaryFallback();
                    }
                    const RooArgList &truthParameters = truthResult->floatParsFinal();
                    for(int parameterIndex = 0;
                        parameterIndex < truthParameters.getSize(); ++parameterIndex) {
                        const RooRealVar *truth = dynamic_cast<const RooRealVar *>(
                            truthParameters.at(parameterIndex));
                        RooRealVar *parameter = truth ? workspace->var(truth->GetName()) : nullptr;
                        if(!parameter) continue;
                        if(parameter->getVal() < parameter->getMin()) {
                            parameter->setVal(parameter->getMin());
                        }
                        if(parameter->getVal() > parameter->getMax()) {
                            parameter->setVal(parameter->getMax());
                        }
                    }
                }
                const bool acceptedAttempt =
                    acceptedToyFit(attemptResult.get()) &&
                    attemptException.empty();
                if(!acceptedAttempt && !combCombBoundaryFallbackActive) {
                    const RooRealVar *attemptNCombComb = fittedParameter(
                        attemptResult.get(), "n_Comb_Comb");
                    if(atCombCombLowerBoundary(attemptNCombComb)) {
                        activateCombCombBoundaryFallback();
                    }
                }
                int nearBoundaryCount = -1;
                const string boundaryParameters = nearBoundaryParameters(
                    attemptResult.get(), nearBoundaryCount);
                const RooRealVar *attemptNPP = attemptResult ?
                    dynamic_cast<const RooRealVar *>(
                        attemptResult->floatParsFinal().find("n_P_P")) : nullptr;
                fitAttempts << toyId << "," << seed << "," << attempts << ","
                    << fitMode << ","
                    << (attemptResult ? attemptResult->status() : -999) << ","
                    << (attemptResult ? attemptResult->covQual() : -999) << ","
                    << (attemptResult ? attemptResult->edm() :
                        std::numeric_limits<double>::quiet_NaN()) << ","
                    << (attemptResult ? attemptResult->minNll() :
                        std::numeric_limits<double>::quiet_NaN()) << ","
                    << (attemptNPP ? attemptNPP->getVal() :
                        std::numeric_limits<double>::quiet_NaN()) << ","
                    << (attemptNPP ? attemptNPP->getError() :
                        std::numeric_limits<double>::quiet_NaN()) << ","
                    << nearBoundaryCount << "," << csvField(boundaryParameters)
                    << "," << csvField(attemptException) << "\n";
                fitAttempts.flush();
                if(attemptResult) fitResult = std::move(attemptResult);
                if(acceptedAttempt) break;
            }

            if(makeToy0Diagnostics && toyId == 0 && fitResult) {
                const bool accepted = acceptedToyFit(fitResult.get());
                const string fitPlotDirectory = resultDirectory +
                    (accepted ? "/toy_0_fit_plots" : "/failed_fit_plots");
                gSystem->mkdir(fitPlotDirectory.c_str(), kTRUE);
                const vector<string> outputNames = {
                    "4D_JpsiMass1", "4D_JpsiMass2",
                    "4D_JpsiCtau1", "4D_JpsiCtau2"
                };
                const vector<string> axisTitles = {
                    "M(J/#psi_{1}) [GeV]", "M(J/#psi_{2}) [GeV]",
                    "c#tau(J/#psi_{1}) [cm]", "c#tau(J/#psi_{2}) [cm]"
                };
                for(int index = 0; index < 4; ++index) {
                    drawToyFitProjectionPlot(
                        fitPlotDirectory, toyData, *pdf, *pdfPP, *pdfPNP,
                        *pdfNPP, *pdfNPNP, *pdfSigComb, *pdfCombSig, *pdfCombComb,
                        *observables[index], outputNames[index], axisTitles[index],
                        index >= 2, toyId, seed, *fitResult);
                }
                if(accepted) {
                    TFile fitArtifact(
                        (resultDirectory + "/toy_0_fit.root").c_str(), "RECREATE");
                    if(fitArtifact.IsZombie()) {
                        throw runtime_error("Cannot create the accepted toy-0 ROOT artifact");
                    }
                    toyData.Write("toy_data");
                    fitResult->Write("fit_result");
                    fitArtifact.Close();
                }
            }

            if(!acceptedToyFit(fitResult.get())) {
                const int status = fitResult ? fitResult->status() : -999;
                const int covarianceQuality = fitResult ? fitResult->covQual() : -999;
                const double edm = fitResult ? fitResult->edm() :
                    std::numeric_limits<double>::quiet_NaN();
                results << toyId << "," << seed << "," << toyWeights.size() << ","
                    << toySumWeights << "," << toySumWeightsSquared << ","
                    << attempts << ",0," << status << "," << covarianceQuality << ","
                    << edm << "," << nPPTruth
                    << ",nan,nan,nan,nan,nan,nan,nan,nan\n";
                results.flush();
                ofstream failure(resultDirectory + "/failure.txt");
                failure << "toy_id=" << toyId << "\nseed=" << seed
                    << "\nattempts=" << attempts << "\nstatus=" << status
                    << "\ncovQual=" << covarianceQuality << "\nedm=" << edm
                    << "\nexception=" << lastException
                    << "\nfailure_artifact=failure_toy_" << toyId << ".root\n";
                failure.close();
                TFile failureArtifact(
                    (resultDirectory + "/failure_toy_" +
                        std::to_string(toyId) + ".root").c_str(), "RECREATE");
                if(failureArtifact.IsZombie()) {
                    throw runtime_error("Cannot create the failed-toy ROOT artifact");
                }
                toyData.Write("toy_data");
                if(fitResult) fitResult->Write("fit_result_last");
                failureArtifact.Close();
                cerr << "Toy " << toyId << " failed after " << attempts
                    << " attempts; stopping without replacement" << endl;
                gSystem->Exit(1);
                return;
            }

            vector<double> toyChi2;
            for(int index = 0; index < 4; ++index) {
                toyChi2.push_back(projectionChi2(
                    *pdf, *projections[index].variable, observables, observableSet,
                    projections[index].edges, nullptr, &toyData, index).chi2);
            }
            const double fittedNPP = nPP->getVal();
            const double fittedNPPError = nPP->getError();
            if(!std::isfinite(fittedNPP) || !std::isfinite(fittedNPPError) ||
               fittedNPPError <= 0.0) {
                throw runtime_error("An accepted toy has an invalid n_P_P result");
            }
            const double pull = (fittedNPP - nPPTruth) / fittedNPPError;
            results << toyId << "," << seed << "," << toyWeights.size() << ","
                << toySumWeights << "," << toySumWeightsSquared << ","
                << attempts << ",1," << fitResult->status() << ","
                << fitResult->covQual() << "," << fitResult->edm() << ","
                << nPPTruth << "," << fittedNPP << "," << fittedNPPError << ","
                << pull << "," << fitResult->minNll();
            for(double chi2 : toyChi2) results << "," << chi2;
            results << "\n";
            results.flush();
            cout << "Accepted toy " << toyId << ": n_P_P=" << fittedNPP
                << " +/- " << fittedNPPError << ", pull=" << pull << endl;
        }
    } catch(const std::exception &error) {
        cerr << "Fit_Check failed: " << error.what() << endl;
        gSystem->Exit(1);
        return;
    }
}
