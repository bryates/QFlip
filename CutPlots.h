#ifndef CUT_PLOTS_H
#define CUT_PLOTS_H

#include "StdPlots.h"
#include "Data.h"
#include "Lepton.h"
#include "SelectionFunctions.h"

class CutPlots : public StdPlots {

 public:
  TH1F *h_muonCharge, *h_mass, *h_MET, *h_PFSumET, *h_PFSumETMinusMu, *h_HT;
  TH1F *h_muonCharge_nw, *h_mass_nw, *h_MET_nw, *h_PFSumET_nw, *h_PFSumETMinusMu_nw, *h_HT_nw;
  TH1F *h_muonCharge_barrel, *h_muonCharge_disk, *h_muonCharge_OS, *h_muonCharge_SS;
  vector<TH1F*> h_charge; //[all mass,>50,Z+/-10]
  TH1F *h_nvtx;
  TH1F *h_nvtx_nw;
  TH1F *h_mass_50, *h_mass_Z, *h_mass_weight, *h_mass_weight_Z;
  TH2F *h_invpt_flip;
  TH1F *h_probe_pt, *h_probe_eta, *h_probe_phi;
  TH1F *h_invpt;
  Double_t noweight;
  Bool_t debug;
 
  CutPlots(TString name);
  ~CutPlots();
  void Fill(Double_t weight, Int_t muonCharge, Double_t MET, Double_t PFSumET, Double_t PFSumETMinusMu, Double_t HT, Double_t eta);
  void Fill(Double_t weight, Double_t mass, Int_t charge);
  void Fill(Double_t weight, std::vector<Lepton> &Coll, Int_t tag = -1, Int_t pdgid = 0);
  void SetVertex(Double_t weight, std::vector<double> &VertexNDF, std::vector<bool> &VertexIsFake, std::vector<double> &VertexX, std::vector<double> &VertexY, std::vector<double> &VertexZ);
  //void SetVertex(Double_t weight, Int_t goodSize);
  void NoWeight(Double_t weight);
  void Write();

};

#endif
