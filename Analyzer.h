#ifndef Analyzer_h
#define Analyzer_h

#include <set>
#include "Data.h"
#include "ElectronSelection.h"
#include "MuonSelection.h"
#include "JetSelection.h"
#include "GenSelection.h"
#include "SelectionFunctions.h"
#include "OtherFunctions.h"
#include "ElectronPlots.h"
#include "MuonPlots.h"
#include "JetPlots.h"
#include "SignalPlots.h"
#include "Reweight.cc"
#include "BTagSFUtil.h"

#include "CutPlots.h"

class Analyzer : public Data {

  TH2F *h_nEvents, *h_nEventsFO, *h_FOrate;

  static const Bool_t debug = false; 
  static const Double_t integratedlumi = 19711.2;//19762.501;
//  static const Double_t integratedlumi = 1.927196301; HLT_Mu5,8
//  static const Double_t integratedlumi = 4.391111; //HLT_Mu8,12
//  static const Double_t integratedlumi = 22.945019; HLT_Mu17
//  static const Double_t integratedlumi = 83.483; HLT_Mu24
//  static const Double_t integratedlumi = 123.9391;
//  static const Double_t integratedlumi = 1.0; //for Fakes
  static const Double_t Mass_Z = 91.1876;
  static const Double_t Mass_W = 80.398;
  static const Double_t trigeff = 0.94;
  static const Double_t mu1scale = 0.927;
  static const Double_t mu2scale = 0.992;

  Double_t *****doubleFake; Double_t ***singleFake; Double_t *****doubleANDsingleFake;
  Double_t *finalbkg1, *finalbkgerror1, *finalbkg2, *finalbkgerror2, *realsingle, *realsingleerror, *realdouble, *realtotal, *doubletosingle, *errdoubletosingle;
  Double_t jets2mass, triggerweight;
  Double_t ptMu0, ptMu1, ID_weight_0, ID_weight_1;
  Int_t tempCharge, index, triggerMu1, triggerMu2;
  UInt_t dataType;
  Double_t Wcand_tmp, Wcand, METcut, METcontrol;

  Bool_t VETO, SINGLEFAKE, DOUBLEFAKE, b_found, muonbad;
  Bool_t triggerMatched[2];
    
  BTagSFUtil *fBTagSF;
  Int_t jetFlavour;
  Double_t METEnDown, METEnUp, METResDown, METResUp, METEUnDown, METEUnUp;
  Bool_t b_foundEffDown, b_foundEffUp, b_foundMissDown, b_foundMissUp;

 public:
  static const Bool_t MC_pu = true;

  ReweightPU *reweightPU;
  TH1F *h_nvtx_norw, *h_nvtx_rw;
  UInt_t numberVertices;
  TString completename;

  Bool_t *goodVerticies;
  Bool_t noTriggerMatched;
  TDirectory *Dir;
  TH1F *h_prova, *h_RelIsoFR;
  TH1F *h_nVertex, *h_nVertex0, *h_nVertex1, *h_nVertex2;
  TH1F *h_nsignal, *h_cutflow;
  TH2F *h_singlefake, *h_doublefake;
  TH1F *h_MET, *h_METsign, *h_MuonMissCharge, *h_EventFakeType;
  TH1F *h_muonCharge, *h_PFSumET, *h_HT;
  TH2F *FRhisto, *h_dRvsbTag;
  TH2F *Mu10_STAT, *Mu10_SYS, *Mu20_STAT, *Mu20_SYS;
  TH2F *ID_Iso;
  TH2I *h_LeptvsVert;

  TFile *outfile;

  Long64_t entrieslimit;
  Double_t METx, METy, MET, dr, MCweight, weight;
  UInt_t VertexN;
  Int_t prescaler, HT;

  MuonSel MuonPOG;//MuonTight, MuonLooseButNOTight, MuonLoose, MuonVeto;
  //GenSel GenTight;
  ElectronSel Electron;//Tight, ElectronLoose;
  JJ  Jets; 
  //  std::vector<Lepton> lepton;
  //ElectronPlots *h_electrons, *h_electronsLoose;
  MuonPlots *h_muonsPOG, *h_muonsPOG2;//*h_muons, *h_muonsLoose, *h_LnotT;// *h_muonCharge;
  JetPlots *h_jets, *h_jets_2muons;
  //SignalPlots *h_signal3;
  /*SignalPlots *h_signal, *h_signalMET50, *h_signalbTag, *h_signalTOT, *h_WZcontrol;
  SignalPlots *h_singlefakes, *h_doublefakes, *h_totalfakes;
  SignalPlots *h_singlefakesMET50, *h_doublefakesMET50, *h_totalfakesMET50;
  SignalPlots *h_singlefakesbTag, *h_doublefakesbTag, *h_totalfakesbTag;
  SignalPlots *h_singlefakesTOT, *h_doublefakesTOT, *h_totalfakesTOT;*/

  CutPlots *h_muons, *h_twoMu, *h_singleIso, *h_pt, *h_PFRange;

  Analyzer();
  ~Analyzer();
  void Loop();
  void SetWeight(Double_t CrossSection, Double_t nevents);
  void SetName(TString name, Int_t version);
  void SetEvtN(Long64_t events);

};
#endif
