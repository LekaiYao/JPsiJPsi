#include "TFile.h"
#include "RooWorkspace.h"

using namespace std;

void Fit_valid() {
    // Read fit model from file
    TFile *fData = new TFile("Model_4D_diff.root");
    RooWorkspace *wk = (RooWorkspace *)fData->Get("wsp");
    double sig = wk->var("n_P_P")->getVal(), sigErr = wk->var("n_P_P")->getError();
    TFile *fRef = new TFile("Model_4D_diff_ref.root");
    RooWorkspace *wr = (RooWorkspace *)fRef->Get("wsp");
    double ref = wr->var("n_P_P")->getVal();
    // Log out put to screen
    cout<<"Event yield: "<<sig<<" +/- "<<sigErr<<"(Sta.) +/- "<<fabs(ref-sig)/sig*100<<"\%(Sys.4)"<<endl;
    return;
}