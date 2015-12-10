#include "AnalyzerFR.h"

AnalyzerFR::AnalyzerFR() {

  h_nEvents = new TH2F ("h_nEvents", "Number of Events",ninteta,arrayeta,nintpT,arraypT);
  h_muonPos = new TH1F ("h_muonPos", "Muon Charge;#eta;",50,-2.5,2.5);
  h_muonNeg = new TH1F ("h_muonNeg", "Muon Charge;#eta;",50,-2.5,2.5);
  h_muonMiss = new TH1F ("h_muonMiss", "Muon Charge;#eta;",50,-2.5,2.5);
  h_muonPos->Sumw2();
  h_muonNeg->Sumw2();
  h_muonMiss->Sumw2();
  //h_muonCharge->Sumw2();
  h_MET = new TH1F ("h_MET", "Missing E_{T};E_{T} (GeV);Events",300,0.0,150.0);//300
  h_METsign = new TH1F ("h_METsign", "Missing E_{T} for two muons;E_{T} (GeV);Events",300,0.0,150.0);//300
  h_PFSumET = new TH1F ("h_PFSumET", "Sum E_{T};E_{T} (GeV};Events",500,0.0,3000.0);//5000
  h_HT = new TH1F ("h_HT", "Sum P_{T} of jets;P_{T} (GeV/c);Events", 100,0.0, 1000.0);//3000

  h_nvtx_norw = new TH1F("h_nvtx_norw","Nvtx per bunch crossing at BX = 0 noreweight",60,0.0,60.0);
  h_nvtx_rw = new TH1F("h_nvtx_rw","Nvtx per bunch crossing at BX = 0 reweight",60,0.0,60.0);

  h_DY_pt_SS = new TH1F("h_DY_pt_SS",";P_{T} (GeV/C);Events", 100, 0.0, 500.);
  h_DY_pt_T = new TH1F("h_DY_pt_T",";P_{T} (GeV/C);Events", 100, 0.0, 500.);
  h_DY_pt_Flip = new TH1F("h_DY_pt_Flip",";P_{T} (GeV/C);Q Flip Rate", 100, 0.0, 500.);

  h_DY_pt_SS->Sumw2();
  h_DY_pt_T->Sumw2();
  h_DY_pt_Flip->Sumw2();

  if (debug) cout<<"inizio"<<endl;

  h_muonsPOG = new MuonPlots("POG_muons");
  h_muonsPOG2 = new MuonPlots("POG_two_muons");
  h_jets = new JetPlots("jets");
  h_jets_2muons = new JetPlots("jets_two_muons");

 
  h_noCuts = new StdPlots("noCuts"); 
  h_muons = new CutPlots("muons");
  h_twoMu = new CutPlots("TwoMuons");
  h_singleIso = new CutPlots("singleIso");
  h_pt = new CutPlots("pt");
  h_PFRange = new CutPlots("PFRange");
  h_METRange = new CutPlots("METRange");
  //h_HTRange = new CutPlots("HTRange");
  //h_Njets = new CutPlots("NJetsRange");
  h_NoJets = new CutPlots("NoJets");
  h_NoJets_SS = new CutPlots("NoJets_SS");
  h_NoJets_OS = new CutPlots("NoJets_OS");
  h_NoJets_Flip_SS = new Fakes("NoJets_Flip_SS");
  h_NoJets_Flip_OS = new Fakes("NoJets_Flip_OS");
  h_ChargeFlip = new Fakes("ChargeFlip");
  //h_ChargeFlipCuts = new Fakes("ChargeFlipCuts");
  h_TagFlip = new Fakes("TagFlip");
  h_DYFlip = new Fakes("DYFlip");
  //h_DYFlipCut = new Fakes("DYFlipCut");
  h_DYFlipTag = new Fakes("DYFlipTag");


  // Fake Rates
  h_nVertex= new TH1F("h_nVertex","number of verteces",60,0,60);
  h_nVertex0= new TH1F("h_nVertex0","number of verteces t0",60,0,60);
  h_nVertex1= new TH1F("h_nVertex1","number of verteces t1",60,0,60);
  h_nVertex2= new TH1F("h_nVertex2","number of verteces t2",60,0,60);

  h_electrons = new ElectronPlots("electrons");
  h_electronsLoose = new ElectronPlots("loose electrons");
  h_muonsFR = new MuonPlots("tight_muons");
  h_muonsLoose = new MuonPlots("loose_muons");
  h_LnotT = new MuonPlots("loose_not_tight");
  //h_muonCharge = new MuonPlots("misscharge_muon");
  //h_jetsFR = new JetPlots("jets_FR");
  h_jets_veto = new JetPlots("jets_w_veto");
  h_signal = new SignalPlots("signal");
  h_WZcontrol = new SignalPlots("WZcontrol");
  h_signalMET50 = new SignalPlots("signal_MET50");
  h_signalbTag = new SignalPlots("signal_bTag");
  h_signalTOT = new SignalPlots("signal_TOT");
  h_singlefakes = new SignalPlots("sf");
  h_doublefakes = new SignalPlots("df");
  h_totalfakes = new SignalPlots("tf");
  h_singlefakesMET50 = new SignalPlots("sf_MET50");
  h_doublefakesMET50 = new SignalPlots("df_MET50");
  h_totalfakesMET50 = new SignalPlots("tf_MET50");
  h_singlefakesbTag = new SignalPlots("sf_bTag");
  h_doublefakesbTag = new SignalPlots("df_bTag");
  h_totalfakesbTag = new SignalPlots("tf_bTag");
  h_singlefakesTOT = new SignalPlots("sf_TOT");
  h_doublefakesTOT = new SignalPlots("df_TOT");
  h_totalfakesTOT = new SignalPlots("tf_TOT");
  h_nsignal = new TH1F("h_signal","number of signal events ",20,-1,19);
  h_cutflow = new TH1F("h_cutflow","number of signal events in cut flow",40,0,40);
  h_singlefake = new TH2F("h_singlefake","number of single fakes ",4,0,4,4,-1,3);
  h_doublefake = new TH2F("h_doublefake","number of double fakes ",4,0,4,4,-1,3);

  //TFile *infile = new TFile("/uscms_data/d2/fgior8/LQntuple_18b/CMSSW_5_3_14_patch2_LQ/src/code/Total_FRcorr40_130.root"); //Total_FRcorr40_PoG.root
  TFile *infile = new TFile("/uscms_data/d3/byates/CMSSW_5_3_14_patch2/src/muon/Total_FRcorr40_PoG.root");

  infile->cd();
  TDirectory *dir=gDirectory;
  dir->GetObject("h_FOrate3",FRhisto);
//  infile->Close();

  TFile *HLT_SF_file = new TFile("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/HLT_DoubleMu_Efficiencies_Run_2012ABCD_53X_L3Match0p1.root");

  HLT_SF_file->cd();
  TDirectory *dirHLT=gDirectory;
  dirHLT->GetObject("DATA_over_MC_Mu17TkMu8_Tight_Mu1_10To20_and_Mu2_20ToInfty_with_STAT_uncrt",Mu10_STAT);
  dirHLT->GetObject("DATA_over_MC_Mu17TkMu8_Tight_Mu1_10To20_and_Mu2_20ToInfty_with_SYST_uncrt",Mu10_SYS);
  dirHLT->GetObject("DATA_over_MC_Mu17TkMu8_Tight_Mu1_20ToInfty_and_Mu2_20ToInfty_with_STAT_uncrt",Mu20_STAT);
  dirHLT->GetObject("DATA_over_MC_Mu17TkMu8_Tight_Mu1_20ToInfty_and_Mu2_20ToInfty_with_SYST_uncrt",Mu20_SYS);
//  HLT_SF_file->Close();

  TFile *ID_Iso_file = new TFile("/uscms_data/d2/fgior8/LQntuple_18b/CMSSW_5_3_14_patch2_LQ/src/code/FinalSF_NoJetRequirement.root");
  TDirectory *dirIDIso=gDirectory;
  dirIDIso->GetObject("etavspt",ID_Iso);


  if (debug) cout<<"fine"<<endl;
}

AnalyzerFR::~AnalyzerFR() { }

void AnalyzerFR::SetName(TString name, Int_t version) {
  completename = name + "_";
  completename += version;
  completename += ".root";
  outfile = new TFile(completename,"RECREATE");
}

void AnalyzerFR::SetWeight(Double_t CrossSection, Double_t nevents) {

  MCweight = integratedlumi * CrossSection / nevents;
  cout<<"MCweight = "<<MCweight<<endl;
 
}

void AnalyzerFR::SetEvtN(Long64_t events) {
  events ? entrieslimit=events :  entrieslimit=-1;
  cout<<"events "<<events<<endl<<"entrieslimit "<<entrieslimit<<endl;
}

void AnalyzerFR::Loop() {

  cout << "total number of entries " <<nentries<<endl;

  if (debug) cout<< "loop begins" <<endl;

  if (debug) cout<< "PU histos loaded" <<endl;

  if(!MCweight) MCweight=1; 

  weight=MCweight;

  if (fChain == 0) 
    cout << "Ciao!" << endl;

  h_muons->NoWeight(weight); //FIXME wet to 1?
  h_twoMu->NoWeight(weight); 
  h_singleIso->NoWeight(weight);
  h_pt->NoWeight(weight);
  h_PFRange->NoWeight(weight);
  reweightPU = new ReweightPU("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/MyDataPileupHistogram_69400.root");

  // Fake Rates
  fBTagSF = new BTagSFUtil("CSVM");

  Double_t SingleFake=0; Double_t DoubleFake=0; Double_t Single_Double=0;
  Int_t nSingleFake=0; Int_t nDoubleFake=0;

  UInt_t nbinX=FRhisto->GetNbinsX(); UInt_t nbinY=FRhisto->GetNbinsY(); UInt_t nSplit=4;

  if (debug) cout<< "Something wrong reading the FR histo" <<endl;

  singleFake=new Double_t**[nSplit];
  doubleFake=new Double_t****[nSplit];
  doubleANDsingleFake=new Double_t ****[nSplit];
  finalbkg1=new Double_t[nSplit];
  finalbkgerror1=new Double_t[nSplit];
  finalbkg2=new Double_t[nSplit];
  finalbkgerror2=new Double_t[nSplit];
  realsingle=new Double_t[nSplit];
  realsingleerror=new Double_t[nSplit];
  realdouble=new Double_t[nSplit];
  realtotal=new Double_t[nSplit];
  doubletosingle=new Double_t[nSplit];
  errdoubletosingle=new Double_t[nSplit];
  for (UInt_t z=0; z<nSplit; z++) {
    singleFake[z]=new Double_t*[nbinX];
    doubleFake[z]=new Double_t***[nbinX];
    doubleANDsingleFake[z]=new Double_t***[nbinX];
    finalbkg1[z]=0;
    finalbkgerror1[z]=0;
    finalbkg2[z]=0;
    finalbkgerror2[z]=0;
    realsingle[z]=0;
    realsingleerror[z]=0;
    realdouble[z]=0;
    realtotal[z]=0;
    doubletosingle[z]=0;
    errdoubletosingle[z]=0;
  }
  for (UInt_t z=0; z<nSplit; z++)
    for (UInt_t i=0; i<nbinX; i++) {
      singleFake[z][i]=new Double_t[nbinY];
      doubleFake[z][i]=new Double_t**[nbinY];
      doubleANDsingleFake[z][i]=new Double_t**[nbinY];
    }
  for (UInt_t z=0; z<nSplit; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++) {
        singleFake[z][i][j]=0;
        doubleFake[z][i][j]=new Double_t*[nbinX];
        doubleANDsingleFake[z][i][j]=new Double_t*[nbinX];
      }
  for (UInt_t z=0; z<nSplit; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++)
        for (UInt_t m=0; m<nbinX; m++) {
          doubleFake[z][i][j][m]=new Double_t[nbinY];
          doubleANDsingleFake[z][i][j][m]=new Double_t[nbinY];
        }
  for (UInt_t z=0; z<nSplit; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++)
        for (UInt_t m=0; m<nbinX; m++)
          for (UInt_t n=0; n<nbinY; n++) {
            doubleFake[z][i][j][m][n]=0;
            doubleANDsingleFake[z][i][j][m][n]=0;
          }



  if (debug) cout<< "at the loop" <<endl;
  std::set<int> runs;
  for (Long64_t jentry = 0; jentry < nentries; jentry++ ) {

    //    watch_getentry.Start(false);
    if (debug) cout<< "Event number " <<jentry<<endl;
    if (debug) cout<<"begin loop"<<endl;
    if (!(jentry % 10000))
      cout << jentry << endl;

    if (!fChain) cout<<"porcaccia"<<endl;
    fChain->GetEntry(jentry);

    if (!MC_pu) {
      if (isTrackingFailure || passTrackingFailureFilter) continue;
      if (!passBeamHaloFilterLoose) continue;
      if (passBadEESupercrystalFilter || passEcalDeadCellBoundaryEnergyFilter || passEcalDeadCellTriggerPrimitiveFilter || passEcalLaserCorrFilter) continue;
      if (!passHBHENoiseFilter) continue; // || passHcalLaserEventFilter) continue;
    }

    std::vector<TString> triggerslist;
    triggerslist.push_back("HLT_IsoMu24_eta2p1_v");
    // Fakes Rates
    triggerslist.push_back("HLT_Mu17_TkMu8_v");

    if ( !TriggerSelector(triggerslist, *HLTInsideDatasetTriggerNames, *HLTInsideDatasetTriggerDecisions, *HLTInsideDatasetTriggerPrescales, prescaler) ) continue;
  
    if (debug) cout<<"trigger passed"<<endl;

    if (MC_pu) {
      /// ***PU reweghting*** ///
      h_nvtx_norw->Fill(PileUpInteractionsTrue->at(0), MCweight);
      weight = reweightPU->GetWeight(PileUpInteractionsTrue->at(0))*MCweight;
      h_nvtx_rw->Fill(PileUpInteractionsTrue->at(0), weight);
    }

    if (debug) cout<<"pileup reweghting applied"<<endl<<"weight: "<<weight<<endl;

    // Vertex Select
    numberVertices = VertexNDF->size();
    goodVerticies = new Bool_t [numberVertices];
    if ( !isGoodEvent(numberVertices, *VertexIsFake, *VertexNDF, *VertexX, *VertexY, *VertexZ, goodVerticies) ) continue;

    for(UInt_t vv=0; vv<VertexNDF->size(); vv++) {
      if(goodVerticies[vv]) {
        VertexN=vv;
        break;
      }
    }
 
    ///// STARTING WITH PHYSICS OBJECTS COLLECTIONS /////

    if (debug) cout<< "Starting physics object collections " <<jentry<<endl;

    // Fake Rates
    std::vector<GenParticle> genTightColl;
    if (MC_pu) {
      GenTight.SetPt(10);
      GenTight.SetEta(3.0);
      GenTight.SetBSdxy(0.20);
      GenTight.GenSelection(*GenParticleEta, *GenParticlePt, *GenParticlePx, *GenParticlePy, *GenParticlePz, *GenParticleEnergy, *GenParticleVX, *GenParticleVY, *GenParticleVZ, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *GenParticlePdgId, *GenParticleStatus, *GenParticleNumDaught, *GenParticleMotherIndex, genTightColl);
    }

    
    std::vector<Lepton> muonPOGColl;
    MuonPOG.SetPt(15); 
    MuonPOG.SetEta(2.4);
    MuonPOG.SetRelIso(0.12);
    MuonPOG.SetChiNdof(10); 
    MuonPOG.SetBSdxy(0.20);
    MuonPOG.SetBSdz(0.50);
    MuonPOG.SetDeposits(400.0,600.0);
    MuonPOG.MuonSelection(*MuonIsPF, *MuonIsGlobal, *MuonEta, *MuonPhi, *MuonPt, *MuonPtError, *MuonEnergy, *MuonPFIsoR04ChargedHadron, *MuonPFIsoR04NeutralHadron, *MuonPFIsoR04Photon, *MuonEcalVetoIso, *MuonHcalVetoIso, *MuonCharge, *MuonGlobalTrkValidHits, *MuonTrkPixelHits, *MuonStationMatches, *MuonTrackLayersWithMeasurement, *MuonGlobalChi2, *MuonTrkVx, *MuonTrkVy, *MuonTrkVz, *MuonTrkD0, *MuonTrkD0Error, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *MuonPFIsoR04PU, muonPOGColl);

    // Fake Rates
    std::vector<Lepton> muonLooseButNOTightColl;
    MuonLooseButNOTight.SetPt(15);
    MuonLooseButNOTight.SetEta(2.4);
    MuonLooseButNOTight.SetRelIso(0.12,0.40);
    MuonLooseButNOTight.SetChiNdof(10,50);
    MuonLooseButNOTight.SetBSdxy(0.20);
    MuonLooseButNOTight.SetBSdz(0.50);
    MuonLooseButNOTight.SetDeposits(400.0,600.0);
    MuonLooseButNOTight.MuonSelection(*MuonIsPF, *MuonIsGlobal, *MuonEta, *MuonPhi, *MuonPt, *MuonPtError, *MuonEnergy, *MuonPFIsoR03ChargedHadron, *MuonPFIsoR03NeutralHadron, *MuonPFIsoR03Photon, *MuonEcalVetoIso, *MuonHcalVetoIso, *MuonCharge, *MuonGlobalTrkValidHits, *MuonTrkPixelHits, *MuonStationMatches, *MuonTrackLayersWithMeasurement, *MuonGlobalChi2, *MuonTrkVx, *MuonTrkVy, *MuonTrkVz, *MuonTrkD0, *MuonTrkD0Error, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *MuonPFIsoR03PU, muonLooseButNOTightColl);

    std::vector<Lepton> muonLooseColl;
    MuonLoose.SetPt(15);
    MuonLoose.SetEta(2.4);
    MuonLoose.SetRelIso(0.40);
    MuonLoose.SetChiNdof(50);
    MuonLoose.SetBSdxy(0.20);
    MuonLoose.SetBSdz(0.50);
    MuonLoose.SetDeposits(400.0,600.0);
    MuonLoose.MuonSelection(*MuonIsPF, *MuonIsGlobal, *MuonEta, *MuonPhi, *MuonPt, *MuonPtError, *MuonEnergy, *MuonPFIsoR03ChargedHadron, *MuonPFIsoR03NeutralHadron, *MuonPFIsoR03Photon, *MuonEcalVetoIso, *MuonHcalVetoIso, *MuonCharge, *MuonGlobalTrkValidHits, *MuonTrkPixelHits, *MuonStationMatches, *MuonTrackLayersWithMeasurement, *MuonGlobalChi2, *MuonTrkVx, *MuonTrkVy, *MuonTrkVz, *MuonTrkD0, *MuonTrkD0Error, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *MuonPFIsoR03PU, muonLooseColl);

    if (muonLooseColl.size() < 2) continue;

    if (debug) cout<<"matching trigger"<<endl;

    //muon Trigger marching and SF
    triggerMatched[0]=false;
    triggerMatched[1]=false;
    triggerweight=0.;
    triggerMu1=triggerMu2=-1.;
    if (muonLooseColl.size() > 0) {
      for (UInt_t i=0; i<muonLooseColl.size(); i++) {
        index=muonLooseColl[i].ilepton();
        if(triggerMu1==-1 && triggerMu2==-1 && MuonHLTDoubleMuonMatched->at(index)) {
          triggerMu1=i;
          triggerMatched[0]=true;
        }
        if(triggerMu1!=-1 && triggerMu2==-1 && MuonHLTDoubleMuonMatched->at(index)) {
          triggerMu2=i;
          triggerMatched[1]=true;
        }
        if(triggerMu1!=-1 && triggerMu2!=-1)
          break;
      }
    }
    if (debug) cout<<"first matching loop"<<endl;
    if(triggerMatched[0] && triggerMatched[1])
      if(muonLooseColl[1].lorentzVec().Pt()<20) {
        triggerweight=Mu10_STAT->GetBinContent(Mu10_STAT->GetXaxis()->FindBin(fabs(muonLooseColl[triggerMu2].eta())),Mu10_STAT->GetYaxis()->FindBin(fabs(muonLooseColl[triggerMu1].eta())));
      }
      else {
        if (fabs(muonLooseColl[triggerMu2].eta())>fabs(muonLooseColl[triggerMu1].eta())) {
          triggerweight=Mu20_STAT->GetBinContent(Mu20_STAT->GetXaxis()->FindBin(fabs(muonLooseColl[triggerMu2].eta())),Mu20_STAT->GetYaxis()->FindBin(fabs(muonLooseColl[triggerMu1].eta())));
        }
        else {
          triggerweight=Mu20_STAT->GetBinContent(Mu20_STAT->GetXaxis()->FindBin(fabs(muonLooseColl[triggerMu1].eta())),Mu20_STAT->GetYaxis()->FindBin(fabs(muonLooseColl[triggerMu2].eta())));
        }
      }
    else
      continue;
    if (MC_pu)
      weight*=triggerweight;

    if (debug) cout<<"trigger matched"<<endl;

    muonLooseColl[0].lorentzVec().Pt()<300. ? ptMu0=muonLooseColl[0].lorentzVec().Pt() : ptMu0=299.;
    muonLooseColl[1].lorentzVec().Pt()<300. ? ptMu1=muonLooseColl[1].lorentzVec().Pt() : ptMu1=299.;
    ID_weight_0 = ID_Iso->GetBinContent(ID_Iso->GetXaxis()->FindBin(fabs(muonLooseColl[0].eta())),ID_Iso->GetYaxis()->FindBin(ptMu0));
    ID_weight_1 = ID_Iso->GetBinContent(ID_Iso->GetXaxis()->FindBin(fabs(muonLooseColl[1].eta())),ID_Iso->GetYaxis()->FindBin(ptMu1));

    if (MC_pu)
      weight*=ID_weight_0*ID_weight_1;

    if (debug) cout<<"Iso and ID weights applied"<<endl;

    std::vector<Lepton> muonVetoColl;
    MuonVeto.SetPt(10);
    MuonVeto.SetEta(2.4);
    MuonVeto.SetRelIso(0.60);
    MuonVeto.SetChiNdof(500);
    MuonVeto.SetBSdxy(20.00);
    MuonVeto.SetBSdz(100.00);
    MuonVeto.SetDeposits(400.0,600.0);
    MuonVeto.LooseMuonSelection(*MuonIsPF, *MuonIsTracker, *MuonIsGlobal, *MuonEta, *MuonPhi, *MuonPt, *MuonPtError, *MuonEnergy, *MuonPFIsoR03ChargedHadron, *MuonPFIsoR03NeutralHadron, *MuonPFIsoR03Photon, *MuonEcalVetoIso, *MuonHcalVetoIso, *MuonCharge, *MuonGlobalTrkValidHits, *MuonTrkPixelHits, *MuonStationMatches, *MuonTrackLayersWithMeasurement, *MuonGlobalChi2, *MuonTrkVx, *MuonTrkVy, *MuonTrkVz, *MuonTrkD0, *MuonTrkD0Error, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *MuonPFIsoR03PU, muonVetoColl);

    std::vector<Lepton> electronTightColl;
    ElectronTight.SetPt(15);
    ElectronTight.SetEta(2.5);
    ElectronTight.SetRelIso(0.15);
    ElectronTight.SetBSdxy(0.01);
    ElectronTight.SetBSdz(0.10);
    ElectronTight.ElectronSelection(*ElectronIsEB, *ElectronIsEE, *ElectronHasTrackerDrivenSeed, *ElectronHasEcalDrivenSeed, *ElectronEta, *ElectronPhi, *ElectronPt, *ElectronEnergy, *ElectronPFPhotonIso03, *ElectronPFNeutralHadronIso03, *ElectronPFChargedHadronIso03, *ElectronCharge, *ElectronGsfCtfScPixCharge, *ElectronMissingHitsEG, *ElectronHasMatchedConvPhot, *ElectronDeltaEtaTrkSC, *ElectronDeltaPhiTrkSC, *ElectronSigmaIEtaIEta, *ElectronHoE, *ElectronCaloEnergy, *ElectronESuperClusterOverP, *ElectronTrackVx, *ElectronTrackVy, *ElectronTrackVz, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), rhoJets, electronTightColl);

    std::vector<Lepton> electronLooseColl;
    ElectronTight.SetPt(10);
    ElectronTight.SetEta(2.5);
    ElectronTight.SetRelIso(0.60);
    ElectronTight.SetBSdxy(20.0);
    ElectronTight.SetBSdz(100.0);
    ElectronTight.ElectronSelectionLoose(*ElectronPassEGammaIDVeto, *ElectronIsEB, *ElectronIsEE, *ElectronHasTrackerDrivenSeed, *ElectronHasEcalDrivenSeed, *ElectronEta, *ElectronPhi, *ElectronPt, *ElectronEnergy, *ElectronPFPhotonIso03, *ElectronPFNeutralHadronIso03, *ElectronPFChargedHadronIso03, *ElectronCharge, *ElectronGsfCtfScPixCharge, *ElectronMissingHitsEG, *ElectronHasMatchedConvPhot, *ElectronDeltaEtaTrkSC, *ElectronDeltaPhiTrkSC, *ElectronSigmaIEtaIEta, *ElectronHoE, *ElectronCaloEnergy, *ElectronESuperClusterOverP, *ElectronTrackVx, *ElectronTrackVy, *ElectronTrackVz, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), rhoJets, electronLooseColl);

    std::vector<Jet> jetVetoColl;
    JetsVeto.SetPt(20);
    JetsVeto.SetEta(2.5);
    JetsVeto.JetSelectionLeptonVeto(*PFJetPileupjetIDpassLooseWP, *PFJetEta, *PFJetPhi, *PFJetPt, *PFJetEnergy, *PFJetNeutralEmEnergyFraction, *PFJetNeutralHadronEnergyFraction, *PFJetChargedEmEnergyFraction, *PFJetChargedHadronEnergyFraction, *PFJetChargedMultiplicity, *PFJetNConstituents, *PFJetCombinedSecondaryVertexBTag, *PFJetClosestVertexWeighted3DSeparation, electronTightColl, muonLooseColl, jetVetoColl);




    //checking if we got a good muon and some trigger
    if(muonPOGColl.size()<1) continue;
    for(UInt_t i=0; i<muonPOGColl.size(); i++)
      h_noCuts->Fill(weight, (Int_t) muonPOGColl.size(), muonPOGColl[i].lorentzVec().Pt(), muonPOGColl[i].eta(), muonPOGColl[i].lorentzVec().Phi());
    //Matching for singleIso trigger
    Bool_t singleIso = false;
    for (UInt_t i=0; i<1; i++) {
      if (MC_pu) {
        if (MuonHLTSingleIsoMuonMatched->at(index) && muonPOGColl[i].lorentzVec().Pt()>30.) {
  	  singleIso = true;
          break;
	}
      }
      else {
        index=muonPOGColl[i].ilepton();
        if (MuonHLTSingleIsoMuonMatched->at(index) && muonPOGColl[i].lorentzVec().Pt()>30.) {
          singleIso = true;
          break;
        }
      }
    }
    if (!singleIso) continue;

    std::vector<Lepton> electronColl;
    Electron.SetPt(10);
    Electron.SetEta(2.5);
    Electron.SetRelIso(0.15);
    //Electron.SetChiNdof(10);
    Electron.SetBSdxy(0.02);
    Electron.SetBSdz(0.10);
    //Electron.SetDeposits(4.0,6.0);
    Electron.ElectronSelection(*ElectronIsEB, *ElectronIsEE, *ElectronHasTrackerDrivenSeed, *ElectronHasEcalDrivenSeed, *ElectronEta, *ElectronPhi, *ElectronPt, *ElectronEnergy, *ElectronPFPhotonIso03, *ElectronPFNeutralHadronIso03, *ElectronPFChargedHadronIso03, *ElectronCharge, *ElectronGsfCtfScPixCharge, *ElectronMissingHitsEG, *ElectronHasMatchedConvPhot, *ElectronDeltaEtaTrkSC, *ElectronDeltaPhiTrkSC, *ElectronSigmaIEtaIEta, *ElectronHoE, *ElectronCaloEnergy, *ElectronESuperClusterOverP, *ElectronTrackVx, *ElectronTrackVy, *ElectronTrackVz, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), rhoJets, electronColl);

    std::vector<Jet> jetColl;
    Jets.SetPt(30);
    Jets.SetEta(2.4);
    Jets.JetSelectionLeptonVeto(*PFJetPileupjetIDpassLooseWP, *PFJetEta, *PFJetPhi, *PFJetPt, *PFJetEnergy, *PFJetNeutralEmEnergyFraction, *PFJetNeutralHadronEnergyFraction, *PFJetChargedEmEnergyFraction, *PFJetChargedHadronEnergyFraction, *PFJetChargedMultiplicity, *PFJetNConstituents, *PFJetCombinedSecondaryVertexBTag, *PFJetClosestVertexWeighted3DSeparation, electronColl, muonPOGColl, jetColl);

    // filling standard plots for muons, electrons and jets
    if (muonPOGColl.size() > 0) {
      for (UInt_t i=0; i<muonPOGColl.size(); i++) {
        index=muonPOGColl[i].ilepton();
        h_muonsPOG->Fill(weight, (Int_t) muonPOGColl.size(), muonPOGColl[i].lorentzVec().Pt(), muonPOGColl[i].eta(), muonPOGColl[i].lorentzVec().Phi(), muonPOGColl[i].charge(), MuonTrkIso->at(index), MuonEcalIso->at(index), MuonHcalIso->at(index), MuonEcalVetoIso->at(index), MuonHcalVetoIso->at(index), MuonPFIsoR03Photon->at(index), MuonPFIsoR03ChargedHadron->at(index), MuonPFIsoR03NeutralHadron->at(index), muonPOGColl[i].chiNdof(), muonPOGColl[i].dxy_BS(), muonPOGColl[i].dz_BS(), MuonPFIsoR03PU->at(index), rhoJets);
      }
    }

    ///// SOME STANDARD PLOTS /////

    std::vector<GenParticle> genBColl;
    if (MC_pu) {
      GenTight.SetPt(10);
      GenTight.SetEta(3.0);
      GenTight.SetBSdxy(0.20);
      GenTight.GenSelectionB(*GenParticleEta, *GenParticlePt, *GenParticlePx, *GenParticlePy, *GenParticlePz, *GenParticleEnergy, *GenParticleVX, *GenParticleVY, *GenParticleVZ, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *GenParticlePdgId, *GenParticleStatus, *GenParticleNumDaught, *GenParticleMotherIndex, genBColl);
    }

    if (muonLooseButNOTightColl.size() > 0) {
      for (UInt_t i=0; i<muonLooseButNOTightColl.size(); i++) {
        index=muonLooseButNOTightColl[i].ilepton();
        h_LnotT->Fill(weight, (Int_t) muonLooseButNOTightColl.size(), muonLooseButNOTightColl[i].lorentzVec().Pt(), muonLooseButNOTightColl[i].eta(), muonLooseButNOTightColl[i].lorentzVec().Phi(), muonLooseButNOTightColl[i].charge(), MuonTrkIso->at(index), MuonEcalIso->at(index), MuonHcalIso->at(index), MuonEcalVetoIso->at(index), MuonHcalVetoIso->at(index), MuonPFIsoR03Photon->at(index), MuonPFIsoR03ChargedHadron->at(index), MuonPFIsoR03NeutralHadron->at(index), muonLooseButNOTightColl[i].chiNdof(), muonLooseButNOTightColl[i].dxy_BS(), muonLooseButNOTightColl[i].dz_BS(), MuonPFIsoR03PU->at(index), rhoJets);
        //      cout << " GEN size " << genBColl.size() << endl;
        for (UInt_t g=0; g<genBColl.size(); g++) {
          if ( genBColl[g].lorentzVec().DeltaR(muonLooseButNOTightColl[i].lorentzVec()) < 0.3 ) {
            h_cutflow->Fill(fabs(genBColl[g].pdgId()));
          }
        }
      }
    }

    if (muonLooseColl.size() > 0) {
      for (UInt_t i=0; i<muonLooseColl.size(); i++) {
        index=muonLooseColl[i].ilepton();
        h_muonsLoose->Fill(weight, (Int_t) muonLooseColl.size(), muonLooseColl[i].lorentzVec().Pt(), muonLooseColl[i].eta(), muonLooseColl[i].lorentzVec().Phi(), muonLooseColl[i].charge(), MuonTrkIso->at(index), MuonEcalIso->at(index), MuonHcalIso->at(index), MuonEcalVetoIso->at(index), MuonHcalVetoIso->at(index), MuonPFIsoR03Photon->at(index), MuonPFIsoR03ChargedHadron->at(index), MuonPFIsoR03NeutralHadron->at(index), muonLooseColl[i].chiNdof(), muonLooseColl[i].dxy_BS(), muonLooseColl[i].dz_BS(), MuonPFIsoR03PU->at(index), rhoJets);
      }
    }

    if (electronTightColl.size() > 0) {
      for (UInt_t i=0; i<electronTightColl.size(); i++) {
        index=electronTightColl[i].ilepton();
        h_electrons->Fill(weight, (Int_t) electronTightColl.size(), electronTightColl[i].lorentzVec().Pt(), electronTightColl[i].eta(), electronTightColl[i].lorentzVec().Phi(), electronTightColl[i].charge(), ElectronTrkIsoDR03->at(index), ElectronEcalIsoDR03->at(index), ElectronHcalIsoDR03->at(index), electronTightColl[i].dxy_BS(), electronTightColl[i].dz_BS(), rhoJets);
      }
    }

    b_found = false;
    if (jetVetoColl.size() > 0) {
      for (UInt_t i=0; i<jetVetoColl.size(); i++) {
        index=jetVetoColl[i].ijet();
        h_jets_veto->Fill( weight, (Int_t) jetVetoColl.size(), jetVetoColl[i].lorentzVec().Pt(), jetVetoColl[i].eta(), jetVetoColl[i].lorentzVec().Phi(), PFJetTrackCountingHighPurBTag->at(index), PFJetJetProbabilityBTag->at(index), jetVetoColl[i].btag_disc(), PFJetClosestVertexWeightedXYSeparation->at(index), PFJetClosestVertexWeightedZSeparation->at(index), PFJetClosestVertexWeighted3DSeparation->at(index) );
        if(jetVetoColl[i].btag_disc()>0.679)
          b_found = true;
      }
    }


    HT = 0.;
    if (jetColl.size() > 0) {
      for (UInt_t i=0; i<jetColl.size(); i++) {
        HT += jetColl[i].lorentzVec().Pt();
        index=jetColl[i].ijet();
        h_jets->Fill( weight, (Int_t) jetColl.size(), jetColl[i].lorentzVec().Pt(), jetColl[i].eta(), jetColl[i].lorentzVec().Phi(), PFJetTrackCountingHighPurBTag->at(index), PFJetJetProbabilityBTag->at(index), jetColl[i].btag_disc(), PFJetClosestVertexWeightedXYSeparation->at(index), PFJetClosestVertexWeightedZSeparation->at(index), PFJetClosestVertexWeighted3DSeparation->at(index) );
      }
    }
    else
      h_jets->Fill( weight, 0, 0, -99.0, -99.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

    h_muons->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    //h_muons->StdPlots::Fill(weight, muonPOGcoll.size(), muonPOGColl[i].lorentzVec().Pt(), muonPOGColl[i].eta(), muonPOGColl[i].lorentzVec().Phi());
    h_muons->Fill(weight, muonPOGColl);
    h_muons->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    //h_muons->Fill(weight, goodVerticies.size());

    //if (HT > 0)
      //h_HT->Fill(HT, weight);
    h_MET->Fill(PFMETType01XYCor->at(0), weight);
    h_PFSumET->Fill(PFSumETType01XYCor->at(0), weight);

    ///// WZ control region /////
    Double_t tmpZcand=0;
    Double_t Zcand=0;
    Bool_t oneOK = false;
    Bool_t twoOK = false;
    Bool_t threeOK = false;
    if (!b_found && MET>40 && muonVetoColl.size()==3 && jetVetoColl.size() >= 2) {
      if (muonPOGColl.size() == 3 && muonPOGColl[0].lorentzVec().Pt()>20.) {
        for (UInt_t l1=0; l1<muonPOGColl.size()-1; l1++)
          for (UInt_t l2=1; l2<muonPOGColl.size(); l2++) {
            if (muonPOGColl[l1].charge()!=muonPOGColl[l2].charge())
              tmpZcand = (muonPOGColl[l1].lorentzVec()+muonPOGColl[l2].lorentzVec()).M();
            if (fabs(tmpZcand-Mass_Z) < fabs(Zcand-Mass_Z))
              Zcand=tmpZcand;
          }
        if (fabs(Zcand-Mass_Z) < 20) {
          if (MC_pu) {
            for(UInt_t g=0; g<genTightColl.size(); g++) {
              if ( fabs( genTightColl[g].lorentzVec().Pt()-MuonMatchedGenParticlePt->at( muonPOGColl[0].ilepton() ) )<1.0 )
                oneOK = true;
              if ( fabs( genTightColl[g].lorentzVec().Pt()-MuonMatchedGenParticlePt->at( muonPOGColl[1].ilepton() ) )<1.0 )
                twoOK = true;
              if ( fabs( genTightColl[g].lorentzVec().Pt()-MuonMatchedGenParticlePt->at( muonPOGColl[2].ilepton() ) )<1.0 )
                threeOK = true;
            }
            if (oneOK && twoOK && threeOK) {
                h_WZcontrol->Fill(MET, muonPOGColl, jetVetoColl, weight, false, false);
                h_nsignal->Fill(-1.,weight);
              }
          }
          else {
            h_WZcontrol->Fill(MET, muonPOGColl, jetVetoColl, weight, false, false);
            h_nsignal->Fill(-1.,weight);
          }
        }
      }
    }
    if (debug) cout << "WZ control done" << endl;
    Double_t Wpair=999.9;
    Double_t temp_Wpair=999.9;
    Double_t lljj=999.9;
    Double_t temp_lljj=999.9;
    if (jetVetoColl.size() >= 2 && muonLooseColl.size()==2) {
      for (UInt_t i=0; i<jetVetoColl.size()-1; i++) {
        for (UInt_t j=1; j<jetVetoColl.size(); j++) {
          temp_Wpair = (jetVetoColl[i].lorentzVec() + jetVetoColl[j].lorentzVec()).M();
          temp_lljj = (muonLooseColl[0].lorentzVec() + muonLooseColl[1].lorentzVec() + jetVetoColl[i].lorentzVec() + jetVetoColl[j].lorentzVec()).M();
          //if ( fabs(temp_lljj-Mass_W) < fabs(lljj-Mass_W) ) {
          if ( fabs(temp_Wpair-Mass_W) < fabs(Wpair-Mass_W) ) {
            lljj = temp_lljj;
            Wpair = temp_Wpair;
          }
        }
      }
    }
    /// ***simple check for double muon invariant mass and 3rd lepton Z veto*** ///
    Double_t masslow=999.9;
    if (muonLooseColl.size() >= 2) {
      masslow = (muonLooseColl[0].lorentzVec() + muonLooseColl[1].lorentzVec()).M();
      //h_prova->Fill(masslow, weight);
    }
    if (masslow < 10.0) continue;
    Double_t mass3rd=999.9;
    Double_t temp_mass3rd=999.9;
    if (muonPOGColl.size() > 2) {
      for(UInt_t i=0; i<muonPOGColl.size()-1; i++)
        for(UInt_t j=i+1; j<muonPOGColl.size(); j++) {
          if ( muonPOGColl[i].charge() != muonPOGColl[j].charge() ) {
            temp_mass3rd = (muonPOGColl[i].lorentzVec() + muonPOGColl[j].lorentzVec()).M();
            if ( fabs(temp_mass3rd-Mass_Z) < fabs(mass3rd-Mass_Z) )
              mass3rd=temp_mass3rd;
          }
        }
    }
    if (mass3rd > (Mass_Z-15) && mass3rd < (Mass_Z+15) ) continue;
    if ( muonVetoColl.size()>2 ) continue;
    if ( electronLooseColl.size()>0 ) continue;
    //if ( jetVetoColl.size()>4 ) continue;
    METcut = 50.;
    METcontrol = 50.;


    ///// BACKGROUND /////
    if (debug) cout<<"Background selection"<<endl;

    DOUBLEFAKE=false;
    Wcand_tmp=Wcand=0;
    if (muonLooseButNOTightColl.size() == 2 && muonPOGColl.size() == 0 && jetVetoColl.size() >= 2) {
      for(UInt_t i=0; i<muonLooseButNOTightColl.size()-1; i++)
        for(UInt_t j=i+1; j<muonLooseButNOTightColl.size(); j++) {
          if (muonLooseButNOTightColl[i].charge()*muonLooseButNOTightColl[j].charge()>0)
            if (muonLooseButNOTightColl[i].lorentzVec().Pt() >=20) {
              if(b_found)
                dataType=2;
              else if (MET>METcontrol)
                dataType=1;
              else if (MET<=METcut)
                dataType=3;
              DOUBLEFAKE=true;
              DoubleFake=DoublebackGround(FRhisto, muonLooseButNOTightColl, i, j, doubleFake, dataType, 1);
              Single_Double=DoubleTOSinglebkg(FRhisto, muonLooseButNOTightColl, i, j);
            }
        }
    }
    if (debug) cout<<"Double done"<<endl;

    if (DOUBLEFAKE) {
      if (debug) cout<<"        Double found"<<endl;
      nDoubleFake++;
      h_nVertex0->Fill(numberVertices, weight*DoubleFake);
      h_doublefakesTOT->Fill(MET, muonLooseButNOTightColl, jetVetoColl, DoubleFake*weight, false, false);
      h_totalfakesTOT->Fill(MET, muonLooseButNOTightColl, jetVetoColl, (DoubleFake+Single_Double)*weight, false, false);
      if (b_found) {
        h_doublefakesbTag->Fill(MET, muonLooseButNOTightColl, jetVetoColl, DoubleFake*weight, false, false);
        h_totalfakesbTag->Fill(MET, muonLooseButNOTightColl, jetVetoColl, (DoubleFake+Single_Double)*weight, false, false);
      }
      else {
        if (MET>METcontrol) {
          h_doublefakesMET50->Fill(MET, muonLooseButNOTightColl, jetVetoColl, DoubleFake*weight, false, false);
          h_totalfakesMET50->Fill(MET, muonLooseButNOTightColl, jetVetoColl, (DoubleFake+Single_Double)*weight, false, false);
        }
        else {
          if (MET<=METcut) {
//          h_nVertex0->Fill(numberVertices, weight);
            h_doublefakes->Fill(MET, muonLooseButNOTightColl, jetVetoColl, DoubleFake*weight, false, false);
            h_totalfakes->Fill(MET, muonLooseButNOTightColl, jetVetoColl, (DoubleFake+Single_Double)*weight, false, false);
          }
        }
      }
    }
    SINGLEFAKE=false;
    Wcand_tmp=Wcand=0;
    if (muonLooseButNOTightColl.size() == 1 && muonPOGColl.size() == 1 && jetVetoColl.size() >= 2) {
      for(UInt_t i=0; i<muonPOGColl.size(); i++)
        for(UInt_t j=0; j<muonLooseButNOTightColl.size(); j++) {
          if (muonLooseButNOTightColl[j].charge()*muonPOGColl[i].charge()>0)
            if (muonLooseButNOTightColl[j].lorentzVec().Pt() >=20 || muonPOGColl[i].lorentzVec().Pt() >=20) {
              if (debug) cout<<"             Single found"<<endl;
              if(b_found)
                dataType=2;
              else if (MET>METcontrol)
                dataType=1;
              else if (MET<=METcut)
                dataType=3;
              SINGLEFAKE=true;
              SingleFake=SinglebackGround(FRhisto, muonLooseButNOTightColl, j, singleFake, dataType, 1);
              DoubleANDSinglebkg(muonPOGColl, i, muonLooseButNOTightColl, j, doubleANDsingleFake, dataType);
              goto endSingle;
            }
        }
    }
  endSingle:

    if (debug) cout<<"Single done, single = "<<SingleFake<<endl;

    if (SINGLEFAKE) {
      nSingleFake++;
      h_nVertex1->Fill(numberVertices, weight*SingleFake);
      h_singlefakesTOT->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
      h_totalfakesTOT->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
      if (b_found) {
        h_singlefakesbTag->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
        h_totalfakesbTag->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
      }
      else {
        if (MET>METcontrol) {
          h_singlefakesMET50->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
          h_totalfakesMET50->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
        }
        else {
          if (MET<=METcut) {
            h_singlefakes->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
            h_totalfakes->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
            //h_nVertex1->Fill(numberVertices, weight);
          }
        }
      }
    }
    /// BACKGROUND END ///


    //now we require events to have two muons

    if (electronColl.size()>0 || muonPOGColl.size()!=2) continue;
    if ( (muonPOGColl[0].lorentzVec()+muonPOGColl[1].lorentzVec()).M() < 20) continue;

    for (UInt_t i=0; i<muonPOGColl.size(); i++) {
      index=muonPOGColl[i].ilepton();
      h_muonsPOG2->Fill(weight, (Int_t) muonPOGColl.size(), muonPOGColl[i].lorentzVec().Pt(), muonPOGColl[0].eta(), muonPOGColl[i].lorentzVec().Phi(), muonPOGColl[i].charge(), MuonTrkIso->at(index), MuonEcalIso->at(index), MuonHcalIso->at(index), MuonEcalVetoIso->at(index), MuonHcalVetoIso->at(index), MuonPFIsoR03Photon->at(index), MuonPFIsoR03ChargedHadron->at(index), MuonPFIsoR03NeutralHadron->at(index), muonPOGColl[i].chiNdof(), muonPOGColl[i].dxy_BS(), muonPOGColl[i].dz_BS(), MuonPFIsoR03PU->at(index), rhoJets);
    }
    
    HT = 0;
    if (jetColl.size() > 0) {
      for (UInt_t i=0; i<jetColl.size(); i++) {
        index=jetColl[i].ijet();
        HT += jetColl[i].lorentzVec().Pt();
        h_jets_2muons->Fill( weight, (Int_t) jetColl.size(), jetColl[i].lorentzVec().Pt(), jetColl[0].eta(), jetColl[i].lorentzVec().Phi(), PFJetTrackCountingHighPurBTag->at(index), PFJetJetProbabilityBTag->at(index), jetColl[i].btag_disc(), PFJetClosestVertexWeightedXYSeparation->at(index), PFJetClosestVertexWeightedZSeparation->at(index), PFJetClosestVertexWeighted3DSeparation->at(index) );
      }
    }
    else
      h_jets_2muons->Fill( weight, 0, 0, -99.0, -99.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 );

    h_twoMu->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    TLorentzVector s = muonPOGColl[0].lorentzVec() + muonPOGColl[1].lorentzVec();
    h_twoMu->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge());
    h_twoMu->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_twoMu->Fill(weight, muonPOGColl);
   
    h_ChargeFlip->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge(), muonPOGColl);
    //h_ChargeFlip->Fill(weight, PFMETPhiType01XYCor->at(0),  muonPOGColl[0].lorentzVec().Phi());
    //h_ChargeFlip->Fill(weight, PFMETPhiType01XYCor->at(0),  muonPOGColl[1].lorentzVec().Phi());
    h_ChargeFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[0].lorentzVec().Pt(), muonPOGColl[0].lorentzVec().Phi());
    h_ChargeFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[1].lorentzVec().Pt(), muonPOGColl[1].lorentzVec().Phi());

    tag = -1;
    if ( fabs(muonPOGColl[0].lorentzVec().Pt()-48.) <= 10. )
      tag = 0;
    else if ( fabs(muonPOGColl[1].lorentzVec().Pt()-48.) <= 10. )
      tag = 1;

/*
    if (((muonPOGColl[0].lorentzVec().Pt() > 20 && tag == 1) || (muonPOGColl[1].lorentzVec().Pt() > 20 && tag == 0)) && PFMETType01XYCor->at(0) < 30) { //muonPOGColl[0].lorentzVec().GetFirstMother == 23 || muonPOGColl[0].lorentzVec().GetFirstMother == 22
      h_ChargeFlipCuts->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge(), muonPOGColl);
      h_ChargeFlipCuts->Fill(weight, PFMETPhiType01XYCor->at(0),  muonPOGColl[0].lorentzVec().Phi());
      h_ChargeFlipCuts->Fill(weight, PFMETPhiType01XYCor->at(0),  muonPOGColl[1].lorentzVec().Phi());
    }
*/

    /*
    if(muonPOGColl[0].charge() > 0)
      h_muonPos->Fill(muonPOGColl[0].eta(),muonPOGColl[0].lorentzVec().Pt(),1);
    else
      h_muonNeg->Fill(muonPOGColl[0].eta(),muonPOGColl[0].lorentzVec().Pt(),1);
    if(muonPOGColl[1].charge() > 0)
      h_muonPos->Fill(muonPOGColl[1].eta(),muonPOGColl[1].lorentzVec().Pt(),1);
    else
      h_muonNeg->Fill(muonPOGColl[1].eta(),muonPOGColl[1].lorentzVec().Pt(),1);
    */
    if(muonPOGColl[0].charge()*muonPOGColl[1].charge() > 0) {
      h_muonPos->Fill(muonPOGColl[0].eta(),1);
      h_muonPos->Fill(muonPOGColl[1].eta(),1);
    }
    else {
      h_muonNeg->Fill(muonPOGColl[0].eta(),1);
      h_muonNeg->Fill(muonPOGColl[1].eta(),1);
    }

    h_METsign->Fill(PFMETType01XYCor->at(0), weight);
    //h_PFSumET_two->Fill(PFSumETType01XYCor->at(0), weight);

    //Check pT = 48+/-10
    bool ptRange = false;
    //if( fabs(muonPOGColl[0].lorentzVec().Pt()-48.) <= 10. ) ptRange = true;
    //else if( fabs(muonPOGColl[1].lorentzVec().Pt()-48.) <= 10. ) ptRange = true;
    if(tag != -1) ptRange = true;
    
    if (!ptRange) continue;

    HT = 0;
    if (jetColl.size() > 0) {
      for (UInt_t i=0; i<jetColl.size(); i++) 
        HT += jetColl[i].lorentzVec().Pt();
    }

    h_pt->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    h_pt->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_pt->Fill(weight, muonPOGColl);
    h_pt->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge());
    /*
    Int_t tag = -1;
    if ( fabs(muonPOGColl[0].lorentzVec().Pt()-48.) <= 10. )
      tag = 0;
    else if ( fabs(muonPOGColl[1].lorentzVec().Pt()-48.) <= 10. )
      tag = 1;
    */
    h_TagFlip->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge(), muonPOGColl);
    if(tag == 1)
    //h_TagFlip->Fill(weight, PFMETPhiType01XYCor->at(0),  muonPOGColl[0].lorentzVec().Phi());
    h_TagFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[0].lorentzVec().Pt(), muonPOGColl[0].lorentzVec().Phi());
    if(tag == 0)
    //h_TagFlip->Fill(weight, PFMETPhiType01XYCor->at(0),  muonPOGColl[1].lorentzVec().Phi());
    h_TagFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[1].lorentzVec().Pt(), muonPOGColl[1].lorentzVec().Phi());

    //Check PFSumET
    bool PFSumETRange = false;
    if ( (PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et())) < 400) PFSumETRange = true;

    //if(!PFSumETRange) continue;

    HT = 0;
    if (jetColl.size() > 0) {
      for (UInt_t i=0; i<jetColl.size(); i++) 
        HT += jetColl[i].lorentzVec().Pt();
    }

    //h_muonCharge->Fill(muonPOGColl[0].charge()*muonPOGColl[1].charge(), weight);
    h_nEvents->Fill(fabs(muonPOGColl[1].eta()),muonPOGColl[1].lorentzVec().Pt(), weight);
    h_PFRange->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    h_PFRange->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_PFRange->Fill(weight, muonPOGColl);

    //MET cuts
    bool METRange = false;
    if( PFMETType01XYCor->at(0) < 54 ) METRange = true; //FoM for pure DY
    //FIXME 99

    if(!METRange) continue;

    h_METRange->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    h_METRange->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_METRange->Fill(weight, muonPOGColl);

    //HT cuts
    /*bool HTRange = false;
    //if( HT < 30) HTRange = true;
    if( (HT < 26) && (jetColl.size() == 0)) HTRange = true;

    if(!HTRange) continue;

    h_HTRange->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    s = muonPOGColl[0].lorentzVec() + muonPOGColl[1].lorentzVec();
    h_HTRange->Fill(weight, s.M());
    h_HTRange->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);

    //N jets cuts
    bool NJetsRange = false;
    if(jetColl.size() <=2) NJetsRange = true;

    if(!NJetsRange) continue;

    h_Njets->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    s = muonPOGColl[0].lorentzVec() + muonPOGColl[1].lorentzVec();
    h_Njets->Fill(weight, s.M());
    h_Njets->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    */

    //No jets
    bool NoJets = false;
    if( jetColl.size() == 0 ) NoJets = true;

    if(!NoJets) continue;

    h_NoJets->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    s = muonPOGColl[0].lorentzVec() + muonPOGColl[1].lorentzVec();
    h_NoJets->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge());
    h_NoJets->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_NoJets->Fill(weight, muonPOGColl);
    
    // Plot same sign
    if(muonPOGColl[0].charge()*muonPOGColl[1].charge() > 0) {
      h_NoJets_SS->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
      h_NoJets_SS->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge());
      h_NoJets_SS->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
      h_NoJets_SS->Fill(weight, muonPOGColl);
      if(tag == 1)
        h_NoJets_Flip_SS->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[0].lorentzVec().Pt(), muonPOGColl[0].lorentzVec().Phi());
      if(tag == 0)
        h_NoJets_Flip_SS->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[1].lorentzVec().Pt(), muonPOGColl[1].lorentzVec().Phi());
    }
  
    s = muonPOGColl[0].lorentzVec() + muonPOGColl[1].lorentzVec();
    // Plot opposite sign
    if(muonPOGColl[0].charge()*muonPOGColl[1].charge() < 0) {
      h_NoJets_OS->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
      h_NoJets_OS->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge());
      h_NoJets_OS->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
      h_NoJets_OS->Fill(weight, muonPOGColl);
      if(tag == 1)
        h_NoJets_Flip_OS->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[0].lorentzVec().Pt(), muonPOGColl[0].lorentzVec().Phi());
      if(tag == 0)
        h_NoJets_Flip_OS->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[1].lorentzVec().Pt(), muonPOGColl[1].lorentzVec().Phi());
    }

    s = muonPOGColl[0].lorentzVec() + muonPOGColl[1].lorentzVec();
    h_DYFlip->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge(), muonPOGColl);
    h_DYFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[0].lorentzVec().Pt(), muonPOGColl[0].lorentzVec().Phi());
    h_DYFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[1].lorentzVec().Pt(), muonPOGColl[1].lorentzVec().Phi());
      if(tag == 1) {
        h_DYFlipTag->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge(), muonPOGColl, tag);
        h_DYFlipTag->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[0].lorentzVec().Pt(), muonPOGColl[0].lorentzVec().Phi());
      }
      if(tag == 0) {
        h_DYFlipTag->Fill(weight, s.M(), muonPOGColl[0].charge()*muonPOGColl[1].charge(), muonPOGColl, tag);
        h_DYFlipTag->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonPOGColl[1].lorentzVec().Pt(), muonPOGColl[1].lorentzVec().Phi());
      }

    for(int num = 0; num <= 1; num++) {
      if(tag == num) continue;
      h_DY_pt_T->Fill(muonPOGColl[num].lorentzVec().Pt());
      if(muonPOGColl[0].charge()*muonPOGColl[1].charge() > 0)
        h_DY_pt_SS->Fill(muonPOGColl[num].lorentzVec().Pt());
    }

    //}

    //h_ChargeFlip->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), muonPOGColl);

    ///// SIGNAL and CONTROL region/////
    if (debug) cout<<"Signal selection"<<endl;

    VETO=false;
    Wcand_tmp=Wcand=0;
    if ( muonPOGColl.size() == 2 && jetVetoColl.size() >= 2 ) {
      for(UInt_t i=0; i<muonPOGColl.size()-1; i++)
        for(UInt_t j=i+1; j<muonPOGColl.size(); j++) {
          if ( muonPOGColl[i].charge()*muonPOGColl[j].charge() > 0 )
            if ( muonPOGColl[i].lorentzVec().Pt() >= 20) {
              if (MC_pu) {
                for(UInt_t g=0; g<genTightColl.size(); g++)
                  if ( fabs( genTightColl[g].lorentzVec().Pt()-MuonMatchedGenParticlePt->at( muonPOGColl[i].ilepton() ) )<1.0 ) {
                    VETO=true;
                    //cout << "UNO" <<endl;
                  }
              }
              else {
                VETO=true;
              //  cout << "UNO" <<endl;
              }
            }
        }
    }

    if (debug) cout<<"Filling signal histos"<<endl;
    if(VETO) {
      h_nVertex2->Fill(numberVertices, weight);
      h_nsignal->Fill(0.,weight);
      h_signalTOT->Fill(MET, muonPOGColl, jetVetoColl, weight, false, false);
      if (b_found) {
        h_nsignal->Fill(2., weight);
        h_signalbTag->Fill(MET, muonPOGColl, jetVetoColl, weight, false, false);
      }
      else
        if (MET>METcontrol) {
          h_nsignal->Fill(1., weight);
          h_signalMET50->Fill(MET, muonPOGColl, jetVetoColl, weight, false, false);
        }
        else {
          if (MET<=METcut) {
            h_nsignal->Fill(3., weight);
            h_signal->Fill(MET, muonPOGColl, jetVetoColl, weight, false, false);
//          h_nVertex2->Fill(numberVertices, weight);
          }
        }
    }

    if (debug) cout<<"cleaning"<<endl;
    muonPOGColl.clear();  muonLooseButNOTightColl.clear(); muonLooseColl.clear(); muonVetoColl.clear();
    electronTightColl.clear(); electronLooseColl.clear(); jetVetoColl.clear();
    genTightColl.clear(); //genBColl.clear();
    if (debug) cout<<"exiting loop"<<endl;
  }
  if (debug) cout<< "out of the loop" <<endl;


  BackGroundEstimate(FRhisto, singleFake, doubleANDsingleFake, doubleFake, finalbkg1, finalbkgerror1, finalbkg2, finalbkgerror2, realsingle, realsingleerror, realdouble, realtotal, doubletosingle, errdoubletosingle);

  cout<<"Single Fake n "<<nSingleFake<<" value "<<SingleFake<<endl;
  for (UInt_t z=0; z<nSplit; z++) {
    h_singlefake->SetBinContent(z+1,3,finalbkg1[z]);
    h_singlefake->SetBinError(z+1,3,finalbkgerror1[z]);
    h_singlefake->SetBinContent(z+1,1,realsingle[z]);
    h_singlefake->SetBinError(z+1,1,realsingleerror[z]);
    h_doublefake->SetBinContent(z+1,3,finalbkg2[z]);
    h_doublefake->SetBinError(z+1,3,finalbkgerror2[z]);
    h_doublefake->SetBinContent(z+1,1,realdouble[z]);
    h_doublefake->SetBinContent(z+1,2,doubletosingle[z]);
    h_doublefake->SetBinError(z+1,2,errdoubletosingle[z]);
  }
  cout<<"Double Fake n "<<nDoubleFake<<" value "<<DoubleFake<<endl;
  cout<<"Single_Double "<<Single_Double<<endl;
  cout<<"totale "<<realtotal[0]<<", of which : "<<realsingle[0]<<" single and "<<realdouble[0]<<" double"<<endl;


  /*}
  if (debug) cout<< "out of the loop" <<endl;*/
  
  outfile->cd();
  //h_muonMiss = (TH2F*)h_muonPos->Clone();
  //h_muonMiss->SetName("h_muonMiss");
  //h_muonMiss->Divide(h_muonNeg);
  h_muonMiss->Divide(h_muonPos,h_muonNeg);
  h_muonMiss->Write();
  h_muonPos->Write();
  h_muonNeg->Write();
  //h_muonCharge->Write();

  h_nvtx_norw->Write();
  h_nvtx_rw->Write();
  //h_pt->Write();
  //h_PFRange->Write();

  h_DY_pt_T->Write();
  h_DY_pt_SS->Write();
  h_DY_pt_Flip->Write();

  Dir = outfile->mkdir("noCuts");
  outfile->cd( Dir->GetName() );
  h_noCuts->Write();
  outfile->cd();

  Dir = outfile->mkdir("Jets");
  outfile->cd( Dir->GetName() );
  h_jets->Write();
  h_jets_2muons->Write();
  outfile->cd();

  Dir = outfile->mkdir("Muons");
  outfile->cd( Dir->GetName() );
  h_muonsPOG->Write();
  //h_MET->Write();
  //h_PFSumET->Write();
  //h_HT->Write();
  h_muons->Write();
  outfile->cd();

  Dir = outfile->mkdir("TwoMuons");
  outfile->cd( Dir->GetName() );
  h_muonsPOG2->Write();
  h_twoMu->Write();
  h_METsign->Write();
  outfile->cd();

  /*Dir = outfile->mkdir("singleIso");
  outfile->cd( Dir->getName() );
  h_singleIso->Write()
  outfile->cd();*/

  Dir = outfile->mkdir("pt");
  outfile->cd( Dir->GetName() );
  h_pt->Write();
  outfile->cd();

/*
  Dir = outfile->mkdir("PFRange");
  outfile->cd( Dir->GetName() );
  h_PFRange->Write();
  outfile->cd();
*/

  Dir = outfile->mkdir("METRange");
  outfile->cd( Dir->GetName() );
  h_METRange->Write();
  outfile->cd();

/*
  Dir = outfile->mkdir("HTRange");
  outfile->cd( Dir->GetName() );
  h_HTRange->Write();
  outfile->cd();
*/
  
  Dir = outfile->mkdir("NoJets");
  outfile->cd( Dir->GetName() );
  h_NoJets->Write();
  outfile->cd();

  Dir = outfile->mkdir("NoJets_SS");
  outfile->cd( Dir->GetName() );
  h_NoJets_SS->Write();
  outfile->cd();
  
  Dir = outfile->mkdir("NoJets_OS");
  outfile->cd( Dir->GetName() );
  h_NoJets_OS->Write();
  outfile->cd();
  
  Dir = outfile->mkdir("NoJets_Flip_SS");
  outfile->cd( Dir->GetName() );
  h_NoJets_Flip_SS->Write();
  outfile->cd();
  
  Dir = outfile->mkdir("NoJets_Flip_OS");
  outfile->cd( Dir->GetName() );
  h_NoJets_Flip_OS->Write();
  outfile->cd();

/*
  Dir = outfile->mkdir("NJetsRange");
  outfile->cd( Dir->GetName() );
  h_Njets->Write();
  outfile->cd();
*/

  Dir = outfile->mkdir("ChargeFlip");
  outfile->cd( Dir->GetName() );
  h_ChargeFlip->Write();
  outfile->cd();

/*
  Dir = outfile->mkdir("ChargeFlipCuts");
  outfile->cd( Dir->GetName() );
  h_ChargeFlipCuts->Write();
  outfile->cd();
*/

  Dir = outfile->mkdir("TagFlip");
  outfile->cd( Dir->GetName() );
  h_TagFlip->Write();
  outfile->cd();

  Dir = outfile->mkdir("DYFlip");
  outfile->cd( Dir->GetName() );
  h_DYFlip->Write();
  outfile->cd();

/*
  Dir = outfile->mkdir("DYFlipCut");
  outfile->cd( Dir->GetName() );
  h_DYFlipCut->Write();
  outfile->cd();
*/

  Dir = outfile->mkdir("DYFlipTag");
  outfile->cd( Dir->GetName() );
  h_DYFlipTag->Write();
  outfile->cd();

  Dir = outfile->mkdir("SingleFakes");
  outfile->cd( Dir->GetName() );
  h_singlefakes->Write();
  h_singlefakesMET50->Write();
  h_singlefakesbTag->Write();
  h_singlefakesTOT->Write();
  outfile->cd();
  Dir = outfile->mkdir("DoubleFakes");
  outfile->cd( Dir->GetName() );
  h_doublefakes->Write();
  h_doublefakesMET50->Write();
  h_doublefakesbTag->Write();
  h_doublefakesTOT->Write();
  outfile->cd();
  Dir = outfile->mkdir("TotalFakes");
  outfile->cd( Dir->GetName() );
  h_totalfakes->Write();
  h_totalfakesMET50->Write();
  h_totalfakesbTag->Write();
  h_totalfakesTOT->Write();
  outfile->cd();

  Dir = outfile->mkdir("debug");
  outfile->cd( Dir->GetName() );
  FRhisto->Write();
  h_muonsLoose->Write();
  h_LnotT->Write();
  h_jets_veto->Write();


  h_nEvents->Write();

  outfile->Close();

}
