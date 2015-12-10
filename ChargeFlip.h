#ifndef CharegeFlip_h
#define ChargeFLip_h

#include "StdPlots.h"
#include "Data.h"
#include "Lepton.h"
#include "SelectionFunctions.h"
#include "TH2.h"

class ChargeFlip : public StdPlots {

 public:
  TH1F *h_chargeFakes, *h_charge_T, *h_charge_SS;
  TH1F *h_charge_pt, *h_charge_pt_SS, *h_charge_pt_T;
  TH1F *h_charge_pt_m50, *h_charge_pt_T_m50, *h_charge_pt_SS_m50;
  TH1F *h_charge_pt_mZ, *h_charge_pt_T_mZ, *h_charge_pt_SS_mZ;
  TH1F *h_deltaphi, *h_mass_weight, *h_mass_SS, *h_mass_T;
  TH2F *h_weight, *h_pt_eta, *h_pt_eta_T, *h_pt_eta_SS;
  TH1F *h_cosdphi;
  TH1F *h_probe_pt, *h_probe_eta, *h_probe_phi;
  TH1F *h_invpt;
  TVector3 MET, v;
  //TH1F *h_METPhi, *h_muPhi;

  ChargeFlip(TString name);
  ~ChargeFlip();
  void Fill(Double_t weight, Double_t mass, Int_t charge, std::vector<Lepton> &Coll, Int_t tag = -1); // kinematics
  void Fill(Double_t weight, Double_t MET, Double_t METphi, Double_t vec, Double_t phi); // deltaphi
  void Write();
};

#endif
