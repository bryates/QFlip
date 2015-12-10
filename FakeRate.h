#ifndef FakeFunctions_h
#define FakeFunctions_h

#include <iostream>
using namespace std;
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"


double DoublebackGround(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate=0);

double SinglebackGround(TH2F* fakerate, std::vector<Lepton>& leptonColl, UInt_t &ilep, Double_t ***fakeN, UInt_t &type, Double_t weight);

double DoubleTOSinglebkg(TH2F* fakerate, std::vector<Lepton>& leptonColl, UInt_t &ilep, UInt_t &jlep);

void DoubleANDSinglebkg(std::vector<Lepton>& leptonColli, UInt_t &ilep, std::vector<Lepton>& leptonCollj, UInt_t &jlep, Double_t *****fakeN, UInt_t &type);

static const double arrayeta[] = {0.0,0.8,1.479,2.0,2.5};
static const double arraypT [] = {10.,15.,20.,25.,30.,35.,45.,60.,80.,100.};
static const Int_t nintpT=7;
static const Int_t ninteta=4;
 
#endif
