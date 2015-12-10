#ifndef FakeFunctions_h
#define FakeFunctions_h

#include <iostream>
using namespace std;
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TLorentzVector.h"


double FailFail(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate=0, Double_t eventWeight=0.);

double FailFailToSingle(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate=0, Double_t eventWeight=0.);

double FailPass(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate=0, Double_t eventWeight=0.);

double PassFail(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate=0, Double_t eventWeight=0.);

double PassPass(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate=0, Double_t eventWeight=0.);
 
#endif
