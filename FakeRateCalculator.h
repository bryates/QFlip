#ifndef FakeRateCalculator_h
#define FakeRateCalculator_h

#include "Analyzer.h"

class FRCalculator : public Analyzer {

  static const Bool_t debug = false;
  static const Int_t nintpT=9;
  static const Double_t minbin = 35.0;
  static const Double_t binwidh = 5.0;
  Double_t *arraypT;
  static const Int_t ninteta=4;
  //Double_t *arrayeta;
  
//  static const Double_t integratedlumi = 1.927196301; HLT_Mu5,8
  static const Double_t effectivelumiMu8Mu12 = 4.391111; //HLT_Mu8,12
  static const Double_t effectivelumiMu17 = 25.014000; //HLT_Mu17
  static const Double_t effectivelumiMu24 = 93.985000; //HLT_Mu24
  static const Double_t effectivelumiMu40 = 18796.0; //HLT_Mu40
//  static const Double_t integratedlumi = 123.9391;
//  static const Double_t integratedlumi = 3.711478; //Period A

  TH1F *h_MT, *h_HT, *h_METPhi, *h_dr, *h_dPhi, *h_ptMuOverptJet;
  TH2F *h_nEvents, *h_nEventsFO, *h_FOrate;
  UInt_t index;
  Double_t HT;
  MuonPlots *h_TLnum, *h_TLden;
  Double_t dPhi, ptMuOverptJet;

 public:

  FRCalculator();
  ~FRCalculator();

  void LoopFR();
};

#endif
