#include <string>

using namespace std;

#define REPHRASE_IO_CONFIGURED
string inputManifest = "";
int expectedInputFiles = 65;
string acceptanceFileName = "";
string efficiencyFileName = "";
string outFile =
    "/eos/home-l/leyao/26JJ/JPsiJPsi/Data_driven/results/"
    "route9p1_nominal_dps_mc/WeightDPS_current.root";

#include "../Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/rephrase.cpp"

void build_dps_reference(
    const char *output =
        "Data_driven/results/dps_reference/WeightDPS_current.root",
    const char *manifest = "",
    const char *acceptance = "",
    const char *efficiency = "",
    bool symmetrizeJpsiLabels = false,
    unsigned int symmetrizationSeed = 50,
    int expectedFiles = 60) {
    outFile = output;
    inputManifest = manifest;
    acceptanceFileName = acceptance;
    efficiencyFileName = efficiency;
    expectedInputFiles = expectedFiles;
    rephrase(symmetrizeJpsiLabels, symmetrizationSeed);
}
