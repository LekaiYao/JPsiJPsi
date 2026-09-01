#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
#include "TLegend.h"
#include "TObject.h"
#include "TSystem.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace RooFit;
using namespace std;

namespace {

vector<string> splitCsv(const string &line) {
    vector<string> fields;
    string field;
    stringstream stream(line);
    while(std::getline(stream, field, ',')) fields.push_back(field);
    return fields;
}

double sampleMean(const vector<double> &values) {
    double sum = 0.0;
    for(double value : values) sum += value;
    return sum / values.size();
}

double sampleStandardDeviation(const vector<double> &values, double mean) {
    double sum = 0.0;
    for(double value : values) {
        const double difference = value - mean;
        sum += difference * difference;
    }
    return std::sqrt(sum / (values.size() - 1));
}

struct ObservedChi2 {
    double value;
    int bins;
    double minimumEffectiveEntries;
};

} // namespace

void Fit_pull(string resultTag="gof_root640_v1", int expectedToys=5000) {
    const string resultDirectory = "fit_results/" + resultTag;
    try {
        if(expectedToys < 2) throw runtime_error("At least two toys are required");
        ifstream input(resultDirectory + "/toy_results.csv");
        if(!input) throw runtime_error("Cannot open toy_results.csv");
        string line;
        if(!std::getline(input, line)) throw runtime_error("toy_results.csv is empty");
        vector<double> pulls;
        vector<double> fittedYields;
        vector<double> fittedErrors;
        vector<vector<double>> toyChi2(4);
        set<int> toyIds;
        double truth = std::numeric_limits<double>::quiet_NaN();
        int failedRows = 0;
        while(std::getline(input, line)) {
            if(line.empty()) continue;
            const vector<string> fields = splitCsv(line);
            if(fields.size() != 19) throw runtime_error("Malformed toy_results.csv row");
            const int toyId = std::stoi(fields[0]);
            if(!toyIds.insert(toyId).second) throw runtime_error("Duplicate toy ID");
            const int accepted = std::stoi(fields[6]);
            if(!accepted) {
                ++failedRows;
                continue;
            }
            if(std::stoi(fields[7]) != 0 || std::stoi(fields[8]) != 3 ||
               std::stod(fields[9]) >= 0.01) {
                throw runtime_error("An accepted toy violates the fit-quality gate");
            }
            truth = std::stod(fields[10]);
            fittedYields.push_back(std::stod(fields[11]));
            fittedErrors.push_back(std::stod(fields[12]));
            pulls.push_back(std::stod(fields[13]));
            for(int index = 0; index < 4; ++index) {
                toyChi2[index].push_back(std::stod(fields[15 + index]));
            }
        }
        input.close();
        if(failedRows != 0) throw runtime_error("toy_results.csv contains failed toys");
        if(static_cast<int>(pulls.size()) != expectedToys ||
           static_cast<int>(toyIds.size()) != expectedToys) {
            throw runtime_error("The number of successful toys does not match the request");
        }
        for(int toyId = 0; toyId < expectedToys; ++toyId) {
            if(!toyIds.count(toyId)) throw runtime_error("Toy IDs are not contiguous");
        }

        map<string, ObservedChi2> observed;
        ifstream observedInput(resultDirectory + "/observed_chi2.csv");
        if(!observedInput || !std::getline(observedInput, line)) {
            throw runtime_error("Cannot read observed_chi2.csv");
        }
        while(std::getline(observedInput, line)) {
            if(line.empty()) continue;
            const vector<string> fields = splitCsv(line);
            if(fields.size() != 4) throw runtime_error("Malformed observed_chi2.csv row");
            observed[fields[0]] = {
                std::stod(fields[1]), std::stoi(fields[2]), std::stod(fields[3])
            };
        }
        const vector<string> projectionNames = {
            "Jpsi_mass1", "Jpsi_mass2", "Jpsi_ctau1", "Jpsi_ctau2"
        };
        if(observed.size() != projectionNames.size()) {
            throw runtime_error("observed_chi2.csv does not contain four projections");
        }

        const double pullMean = sampleMean(pulls);
        const double pullRms = sampleStandardDeviation(pulls, pullMean);
        const double yieldMean = sampleMean(fittedYields);
        const double yieldRms = sampleStandardDeviation(fittedYields, yieldMean);
        const double errorMean = sampleMean(fittedErrors);
        const double minimumPull = *std::min_element(pulls.begin(), pulls.end());
        const double maximumPull = *std::max_element(pulls.begin(), pulls.end());
        const double plotMinimum = -5.0;
        const double plotMaximum = 5.0;
        const int plotBins = 50;
        if(minimumPull <= plotMinimum || maximumPull >= plotMaximum) {
            throw runtime_error("A pull is outside the fixed [-5,5] formal plot range");
        }

        TTree pullTree("pull_tree", "pull_tree");
        double pullValue = 0.0;
        pullTree.Branch("pull", &pullValue);
        for(double value : pulls) {
            pullValue = value;
            pullTree.Fill();
        }
        RooRealVar pullVariable("pull", "pull", plotMinimum, plotMaximum);
        RooArgSet pullVariables(pullVariable);
        RooDataSet pullData(
            "pull_data", "pull_data", pullVariables, Import(pullTree));
        RooRealVar gaussianMean(
            "gaussian_mean", "gaussian_mean", pullMean, -1.0, 1.0);
        RooRealVar gaussianSigma(
            "gaussian_sigma", "gaussian_sigma", pullRms, 0.2, 5.0);
        RooGaussian gaussian("gaussian", "gaussian", pullVariable, gaussianMean, gaussianSigma);
        unique_ptr<RooFitResult> gaussianResult(gaussian.fitTo(
            pullData, Save(), PrintLevel(-1), Warnings(kFALSE), Verbose(kFALSE)));
        if(!gaussianResult || gaussianResult->status() != 0 ||
           gaussianResult->covQual() != 3) {
            throw runtime_error("The full-range Gaussian pull fit failed");
        }

        TCanvas canvas("pull_canvas", "pull_canvas", 1200, 1000);
        canvas.SetLeftMargin(0.12);
        canvas.SetBottomMargin(0.14);
        unique_ptr<RooPlot> frame(pullVariable.frame(
            Title("n_{P,P} pull distribution"), Bins(plotBins)));
        pullData.plotOn(
            frame.get(), DataError(RooAbsData::SumW2), Name("Data"));
        gaussian.plotOn(frame.get(), LineColor(kBlue), LineWidth(2), Name("Gaussian"));
        frame->GetXaxis()->SetTitle(
            "Pull = (N_{P,P}^{fit}-N_{P,P}^{true})/#sigma_{P,P}^{native}");
        frame->GetXaxis()->CenterTitle(kTRUE);
        frame->GetXaxis()->SetTitleOffset(1.25);
        frame->GetYaxis()->SetTitle("Toys / 0.2");
        frame->Draw();
        ostringstream meanLabel;
        meanLabel << fixed << setprecision(3) << "#mu = "
            << gaussianMean.getVal() << " #pm " << gaussianMean.getError();
        ostringstream sigmaLabel;
        sigmaLabel << fixed << setprecision(3) << "#sigma = "
            << gaussianSigma.getVal() << " #pm " << gaussianSigma.getError();
        TLegend legend(0.60, 0.66, 0.88, 0.88);
        legend.AddEntry(frame->findObject("Data"), "Toy fits", "lep");
        legend.AddEntry(frame->findObject("Gaussian"), "Gaussian fit", "l");
        legend.AddEntry((TObject *)nullptr, meanLabel.str().c_str(), "");
        legend.AddEntry((TObject *)nullptr, sigmaLabel.str().c_str(), "");
        legend.Draw();
        canvas.SaveAs((resultDirectory + "/pull_distribution.pdf").c_str());
        canvas.SaveAs((resultDirectory + "/pull_distribution.png").c_str());

        ofstream summary(resultDirectory + "/gof_summary.txt");
        summary << setprecision(17)
            << "status=complete\n"
            << "toys_requested=" << expectedToys << "\n"
            << "toys_accepted=" << pulls.size() << "\n"
            << "toys_failed=" << failedRows << "\n"
            << "n_P_P_truth=" << truth << "\n"
            << "n_P_P_fit_mean=" << yieldMean << "\n"
            << "n_P_P_fit_sample_sd=" << yieldRms << "\n"
            << "n_P_P_bias=" << yieldMean - truth << "\n"
            << "n_P_P_error_mean=" << errorMean << "\n"
            << "pull_sample_mean=" << pullMean << "\n"
            << "pull_sample_sd=" << pullRms << "\n"
            << "pull_min=" << minimumPull << "\n"
            << "pull_max=" << maximumPull << "\n"
            << "pull_plot_min=" << plotMinimum << "\n"
            << "pull_plot_max=" << plotMaximum << "\n"
            << "pull_plot_bins=" << plotBins << "\n"
            << "pull_fit_method=unbinned Gaussian fit over the fixed full range [-5,5]\n"
            << "pull_gaussian_mean=" << gaussianMean.getVal() << "\n"
            << "pull_gaussian_mean_error=" << gaussianMean.getError() << "\n"
            << "pull_gaussian_sigma=" << gaussianSigma.getVal() << "\n"
            << "pull_gaussian_sigma_error=" << gaussianSigma.getError() << "\n"
            << "pull_gaussian_fit_status=" << gaussianResult->status() << "\n"
            << "pull_gaussian_fit_covQual=" << gaussianResult->covQual() << "\n"
            << "pull_gaussian_fit_edm=" << gaussianResult->edm() << "\n";
        for(int index = 0; index < 4; ++index) {
            const ObservedChi2 &dataChi2 = observed.at(projectionNames[index]);
            int tailCount = 0;
            for(double toyValue : toyChi2[index]) {
                if(toyValue >= dataChi2.value) ++tailCount;
            }
            const double pValue = static_cast<double>(tailCount + 1) /
                static_cast<double>(expectedToys + 1);
            const double pError = std::sqrt(
                pValue * (1.0 - pValue) / static_cast<double>(expectedToys + 1));
            summary << projectionNames[index] << "_chi2=" << dataChi2.value << "\n"
                << projectionNames[index] << "_bins=" << dataChi2.bins << "\n"
                << projectionNames[index] << "_minimum_effective_entries="
                << dataChi2.minimumEffectiveEntries << "\n"
                << projectionNames[index] << "_toy_tail_count=" << tailCount << "\n"
                << projectionNames[index] << "_empirical_p=" << pValue << "\n"
                << projectionNames[index] << "_empirical_p_mc_error=" << pError << "\n";
        }
        summary.close();

        TFile output((resultDirectory + "/pull_fit.root").c_str(), "RECREATE");
        gaussianResult->Write("pull_gaussian_fit_result");
        pullTree.Write();
        canvas.Write("pull_canvas");
        output.Close();
        cout << setprecision(10)
            << "Completed " << expectedToys << " toys: pull mean=" << pullMean
            << ", pull SD=" << pullRms << ", Gaussian mean="
            << gaussianMean.getVal() << ", Gaussian sigma="
            << gaussianSigma.getVal() << endl;
    } catch(const std::exception &error) {
        cerr << "Fit_pull failed: " << error.what() << endl;
        gSystem->Exit(1);
        return;
    }
}
