#include <string>

using namespace std;

#define REPHRASE_IO_CONFIGURED
#define N_DIR 1
string prefix[N_DIR] = {""};
string infix =
    "/eos/home-l/leyao/26JJ/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/"
    "NtupleAnalyzer/DPS_ntuple/Ntuple_2016_DPS";
int suffix[N_DIR] = {65};
string outFile =
    "/eos/home-l/leyao/26JJ/JPsiJPsi/Data_driven/results/"
    "route9p1_nominal_dps_mc/WeightDPS_current.root";

#include "../Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/rephrase.cpp"

void build_dps_reference(
    const char *output =
        "Data_driven/results/dps_reference/WeightDPS_current.root") {
    outFile = output;
    rephrase();
}
