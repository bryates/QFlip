#ifndef SignalPlots_h
#define SignalPlots_h

#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "Lepton.h"
#include "Jet.h"
#include "StdPlots.h"
#include <iostream>
class SignalPlots : public StdPlots {
 
  Double_t dijetmass_tmp, dijetmass;
  static const Double_t Mass_W = 80.398;
  static const Double_t pi = 3.1415926535;
 
 public:
  TH1F *h_jjmass, *h_llmass, *h_l1jjmass, *h_l2jjmass, *h_lljjmass, *h_MET, *h_bTag;
  TH2F *h_WandNmass;
  TH3F *h_3Dparm, *h_3DparmJ;
  TH1F *h_leadingMuonPt, *h_secondMuonPt, *h_leadingJetPt, *h_secondJetPt;
  TH1F *h_leadingMuonIso, *h_secondMuonIso, *h_paircharge, *h_muonseta, *h_jetseta;
  TH1F *h_cosTheta1, *h_cosTheta2;
  TH1F *h_DeltaPhi;
  TH1F *h_dMETphilead, *h_cosdMETphilead;
  TH1F *h_dMETphisecond, *h_cosdMETphisecond;
  TH1F *h_Njets;
  TH1F *h_leadingMuonEta, *h_secondMuonEta;
  TH1F *h_leadingMuonPhi, *h_secondMuonPhi;
  TH1F *h_inv_leadingMuonPt, *h_inv_secondMuonPt;
  TH1F *h_invpt;
 
  SignalPlots(TString name);
  ~SignalPlots();
  void Fill(Double_t MET, std::vector<Lepton>& muons, std::vector<Jet>& jets, Double_t weight, Bool_t ptok, Bool_t ssok);
  void Fill(Double_t MET, std::vector<Lepton>& muons, std::vector<Lepton>& muonsloose, std::vector<Jet>& jets, Double_t weight, Bool_t ptok, Bool_t ssok);
  void Fill(Double_t MET, Double_t METphi, std::vector<Lepton>& muons, Double_t weight, Bool_t ptok, Bool_t ssok);
  void Fill(Double_t MET, Double_t METphi, std::vector<Lepton>& muons, std::vector<Lepton>& muonsloose, Double_t weight, Bool_t ptok, Bool_t ssok);
  void Write();

};


#endif
