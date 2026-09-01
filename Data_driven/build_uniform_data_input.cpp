#include <string>

using namespace std;

#define REPHRASE_IO_CONFIGURED
#define N_DIR 7

string prefix[N_DIR] = {
    "B/Ntuple_2016_B", "C/Ntuple_2016_C", "D/Ntuple_2016_D",
    "E/Ntuple_2016_E", "F/Ntuple_2016_F", "G/Ntuple_2016_G",
    "H/Ntuple_2016_H"};
string infix =
    "/eos/user/c/chensh/JPsiJPsi/Data/ULntuple16/CMSSW_10_6_20/src/"
    "NtupleAnalyzer/";
int suffix[N_DIR] = {20, 9, 14, 3, 8, 29, 36};
string outFile =
    "Data_driven/results/unified_input/"
    "WeightData_mJJ7p5.root";

#include "../Data/ULntuple16/CMSSW_10_6_20/src/NtupleAnalyzer/rephrase.cpp"

void build_uniform_data_input(
    const char *output =
        "Data_driven/results/unified_input/WeightData_mJJ7p5.root") {
    outFile = output;
    rephrase();
}
