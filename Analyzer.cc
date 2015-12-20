#include "Analyzer.h"

Analyzer::Analyzer() {

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

  h_DY_pt_SS = new TH1F("h_DY_pt_SS",";P_{T} (GeV/C);Events", 60, 0.0, 300.);
  h_DY_pt_T = new TH1F("h_DY_pt_T",";P_{T} (GeV/C);Events", 60, 0.0, 300.);
  h_DY_pt_Flip = new TH1F("h_DY_pt_Flip",";P_{T} (GeV/C);Q Flip Rate", 60, 0.0, 300.);
  h_genFlip = new TH2F("h_genFlip",";P_{T} (GeV/C);Gen vs Reco flip", 60, 0.0, 300.,5,-2,2);

  h_DY_pt_SS->Sumw2();
  h_DY_pt_T->Sumw2();
  h_DY_pt_Flip->Sumw2();

  h_TV = new TH1F("h_TV",";P_{T} (GeV/C);Events", 50, 0.0, 50.);
  h_TV2 = new TH1F("h_TV2",";P_{T} (GeV/C);Events", 50, 0.0, 50.);

  if (debug) cout<<"inizio"<<endl;

  h_muonsPOG = new MuonPlots("POG_muons");
  h_muonsPOG2 = new MuonPlots("POG_two_muons");
  h_jets = new JetPlots("jets");
  h_jets_2muons = new JetPlots("jets_two_muons");

  h_muonsRej = new StdPlots("Rej_muons");
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
  h_NoJets_Flip_SS = new ChargeFlip("NoJets_Flip_SS");
  h_NoJets_Flip_OS = new ChargeFlip("NoJets_Flip_OS");
  h_ChargeFlip = new ChargeFlip("ChargeFlip");
  //h_ChargeFlipCuts = new ChargeFlip("ChargeFlipCuts");
  h_TagFlip = new ChargeFlip("TagFlip");
  h_DYFlip = new ChargeFlip("DYFlip");
  //h_DYFlipCut = new ChargeFlip("DYFlipCut");
  h_DYFlipTag = new ChargeFlip("DYFlipTag");


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
  h_signalTOTDY = new SignalPlots("signal_TOTDY");
  h_singlefakes = new SignalPlots("sf");
  h_doublefakes = new SignalPlots("df");
  h_singlefakesDY = new SignalPlots("sfDY");
  h_doublefakesDY = new SignalPlots("dfDY");
  h_totalfakes = new SignalPlots("tf");
  h_totalfakesDY = new SignalPlots("tfDY");
  h_singlefakesMET50 = new SignalPlots("sf_MET50");
  h_doublefakesMET50 = new SignalPlots("df_MET50");
  h_totalfakesMET50 = new SignalPlots("tf_MET50");
  h_singlefakesbTag = new SignalPlots("sf_bTag");
  h_doublefakesbTag = new SignalPlots("df_bTag");
  h_totalfakesbTag = new SignalPlots("tf_bTag");
  h_singlefakesTOT = new SignalPlots("sf_TOT");
  h_doublefakesTOT = new SignalPlots("df_TOT");
  h_totalfakesTOT = new SignalPlots("tf_TOT");
  h_singlefakesTOTDY = new SignalPlots("sf_TOTDY");
  h_doublefakesTOTDY = new SignalPlots("df_TOTDY");
  h_totalfakesTOTDY = new SignalPlots("tf_TOTDY");
  h_nsignal = new TH1F("h_signal","number of signal events ",20,-1,19);
  h_cutflow = new TH1F("h_cutflow","number of signal events in cut flow",40,0,40);
  h_singlefake = new TH2F("h_singlefake","number of single fakes ",4,0,4,4,-1,3);
  h_doublefake = new TH2F("h_doublefake","number of double fakes ",4,0,4,4,-1,3);
  h_singlefakeDY = new TH2F("h_singlefakeDY","number of single fakes DY ",4,0,4,4,-1,3);
  h_doublefakeDY = new TH2F("h_doublefakeDY","number of double fakes DY ",4,0,4,4,-1,3);

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

Analyzer::~Analyzer() { }

void Analyzer::SetName(TString name, Int_t version) {
  completename = name + "_";
  completename += version;
  completename += ".root";
  outfile = new TFile(completename,"RECREATE");
}

void Analyzer::SetWeight(Double_t CrossSection, Double_t nevents) {

  MCweight = integratedlumi * CrossSection / nevents;
  cout<<"MCweight = "<<MCweight<<endl;
 
}

void Analyzer::SetEvtN(Long64_t events) {
  events ? entrieslimit=events :  entrieslimit=-1;
  cout<<"events "<<events<<endl<<"entrieslimit "<<entrieslimit<<endl;
}

void Analyzer::Loop() {

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
  //fBTagSF = new BTagSFUtil("CSVM");

  Double_t SingleFake=0; Double_t DoubleFake=0; Double_t Single_Double=0;
  Double_t SingleFakeDY=0; Double_t DoubleFakeDY=0; Double_t Single_DoubleDY=0;
  Int_t nSingleFake=0; Int_t nDoubleFake=0;
  Int_t nSingleFakeDY=0; Int_t nDoubleFakeDY=0;

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

  UInt_t nSplitDY=2;

  singleFakeDY=new Double_t**[nSplitDY];
  doubleFakeDY=new Double_t****[nSplitDY];
  doubleANDsingleFakeDY=new Double_t ****[nSplitDY];
  finalbkg1DY=new Double_t[nSplitDY];
  finalbkgerror1DY=new Double_t[nSplitDY];
  finalbkg2DY=new Double_t[nSplitDY];
  finalbkgerror2DY=new Double_t[nSplitDY];
  realsingleDY=new Double_t[nSplitDY];
  realsingleerrorDY=new Double_t[nSplitDY];
  realdoubleDY=new Double_t[nSplitDY];
  realtotalDY=new Double_t[nSplitDY];
  doubletosingleDY=new Double_t[nSplitDY];
  errdoubletosingleDY=new Double_t[nSplitDY];
  for (UInt_t z=0; z<nSplitDY; z++) {
    singleFakeDY[z]=new Double_t*[nbinX];
    doubleFakeDY[z]=new Double_t***[nbinX];
    doubleANDsingleFakeDY[z]=new Double_t***[nbinX];
    finalbkg1DY[z]=0;
    finalbkgerror1DY[z]=0;
    finalbkg2DY[z]=0;
    finalbkgerror2DY[z]=0;
    realsingleDY[z]=0;
    realsingleerrorDY[z]=0;
    realdoubleDY[z]=0;
    realtotalDY[z]=0;
    doubletosingleDY[z]=0;
    errdoubletosingleDY[z]=0;
  }
  for (UInt_t z=0; z<nSplitDY; z++)
    for (UInt_t i=0; i<nbinX; i++) {
      singleFakeDY[z][i]=new Double_t[nbinY];
      doubleFakeDY[z][i]=new Double_t**[nbinY];
      doubleANDsingleFakeDY[z][i]=new Double_t**[nbinY];
    }
  for (UInt_t z=0; z<nSplitDY; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++) {
        singleFakeDY[z][i][j]=0;
        doubleFakeDY[z][i][j]=new Double_t*[nbinX];
        doubleANDsingleFakeDY[z][i][j]=new Double_t*[nbinX];
      }
  for (UInt_t z=0; z<nSplitDY; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++)
        for (UInt_t m=0; m<nbinX; m++) {
          doubleFakeDY[z][i][j][m]=new Double_t[nbinY];
          doubleANDsingleFakeDY[z][i][j][m]=new Double_t[nbinY];
        }
  for (UInt_t z=0; z<nSplitDY; z++)
    for (UInt_t i=0; i<nbinX; i++)
      for (UInt_t j=0; j<nbinY; j++)
        for (UInt_t m=0; m<nbinX; m++)
          for (UInt_t n=0; n<nbinY; n++) {
            doubleFakeDY[z][i][j][m][n]=0;
            doubleANDsingleFakeDY[z][i][j][m][n]=0;
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
    //triggerslist.push_back("HLT_Mu17_TkMu8_v");

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

    // Gen Matching
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
    //MuonPOG.SetRelIso(0.05);
    MuonPOG.SetChiNdof(10); 
    MuonPOG.SetBSdxy(0.20);
    //MuonPOG.SetBSdxy(0.005);
    MuonPOG.SetBSdz(0.50);
    //MuonPOG.SetBSdz(0.10);
    MuonPOG.SetDeposits(400.0,600.0);
    MuonPOG.MuonSelection(*MuonIsPF, *MuonIsGlobal, *MuonEta, *MuonPhi, *MuonPt, *MuonPtError, *MuonEnergy, *MuonPFIsoR04ChargedHadron, *MuonPFIsoR04NeutralHadron, *MuonPFIsoR04Photon, *MuonEcalVetoIso, *MuonHcalVetoIso, *MuonCharge, *MuonGlobalTrkValidHits, *MuonTrkPixelHits, *MuonStationMatches, *MuonTrackLayersWithMeasurement, *MuonGlobalChi2, *MuonTrkVx, *MuonTrkVy, *MuonTrkVz, *MuonTrkD0, *MuonTrkD0Error, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *MuonPFIsoR04PU, muonPOGColl);

  
    for(UInt_t i = 0; i < muonPOGColl.size(); i++) {
      h_TV->Fill((Int_t) muonPOGColl.size(),weight);
    }
    if(muonPOGColl.size()==2 && muonPOGColl[0].charge()*muonPOGColl[1].charge() < 0)
      for(UInt_t i = 0; i < muonPOGColl.size(); i++) {
        h_TV2->Fill((Int_t) muonPOGColl.size(),weight);
      }

    std::vector<Lepton> muonGenColl;
    std::vector<Lepton> muonRejColl;
    if(MC_pu && Gen) {
      bool match = false;
      double dRtmp = 1;
      double deltaR = 1;
      std::vector<UInt_t> genMatch;
      UInt_t gen = -1;
      // Loop over POG collection
      for(UInt_t i = 0; i < muonPOGColl.size(); i++) {
        // Find potential matches with all gen particles
        for(UInt_t j = 0; j < genTightColl.size(); j++) {
        //for(UInt_t j = 0; j < GenParticlePt->size(); j++) {
          // Skip already matched gen particles
          if(std::find(genMatch.begin(), genMatch.end(), j) != genMatch.end()) continue;
          //if(! (fabs(GenParticlePdgId->at(j))==13 || fabs(GenParticlePdgId->at(j))==15)) continue;
          //if(fabs(GenParticlePdgId->at(j)) == 2212) continue;
          //if(! (fabs(GenParticlePdgId->at(j))==13 && fabs(GenParticlePdgId->at(GenParticleMotherIndex->at(j)))==23 && fabs(GenParticlePdgId->at(GenParticleMotherIndex->at(j)))==24)) continue;
          //int ilep = muonPOGColl[i].ilepton();
          //TLorentzVector Gen;
          //Gen.SetPtEtaPhiE(GenParticlePt->at(j), GenParticleEta->at(j), GenParticlePhi->at(j), GenParticleEnergy->at(j));
          //dRtmp = muonPOGColl[i].lorentzVec().DeltaR(Gen);
          dRtmp = muonPOGColl[i].lorentzVec().DeltaR(genTightColl[j].lorentzVec());
          // Find match with smallest DeltaR (always < 0.3)
          if(dRtmp < 0.3 && dRtmp < deltaR) {
            match = true;
            deltaR = dRtmp;
            //deltaR = muonPOGColl[i].lorentzVec().DeltaR(genTightColl[j].lorentzVec());
            gen = j; // Index of best match
          }
        }
        // Build vector of matches and rejects
        if(match) {
          //muonPOGColl[i].SetMatchedId(GenParticlePdgId->at(gen));
          muonPOGColl[i].SetMatchedId(GenParticlePdgId->at(genTightColl[gen].ilepton()));
          muonGenColl.push_back(muonPOGColl[i]);
          genMatch.push_back(gen);
          //if(muonPOGColl[i].charge()*genTightColl[gen].pdgId() > 0)
            h_genFlip->Fill(muonPOGColl[i].lorentzVec().Pt(), muonPOGColl[i].charge()*muonPOGColl[i].MatchedCharge());

          //Debugging for charge flip
          if(muonPOGColl.size() != 2) continue;
          //if(muonPOGColl[0].charge()*muonPOGColl[1].charge() < 0) continue;
          if(fabs((muonPOGColl[0].lorentzVec()+muonPOGColl[1].lorentzVec()).M() - 91) > 10) continue;
          cout << "PDGID, charge, pT, phi, eta, status, mother" << endl;
          cout << "Muon: " << endl;
          cout << muonPOGColl[i].MatchedId() << " " << muonPOGColl[i].charge() << " " << muonPOGColl[i].lorentzVec().Pt() << " " << muonPOGColl[i].lorentzVec().Phi() << " " << muonPOGColl[i].eta() << endl;
          
          cout << "Gen:" << endl;
/*
          for(int igen = 0; igen < GenParticlePt->size(); igen++) {
            //if(! (fabs(GenParticlePdgId->at(igen))==13 || fabs(GenParticlePdgId->at(igen))==15)) continue;
            if(fabs(GenParticlePdgId->at(j)) == 2212) continue;
            if(! (fabs(GenParticlePdgId->at(igen))==13 || fabs(GenParticlePdgId->at(igen))==15)) cout << "**** PARTICLE ****" << endl;
            cout << GenParticlePdgId->at(igen) << " G " << GenParticlePt->at(igen) << " " << GenParticlePhi->at(igen) << " " << GenParticleEta->at(igen) << " " << GenParticleStatus->at(igen) << " " <<  GenParticlePdgId->at(GenParticleMotherIndex->at(igen)) << endl;
          }
*/
          for(UInt_t it = 0; it < genTightColl.size(); it++) {
            int igen = genTightColl[it].ilepton();
            cout << GenParticlePdgId->at(igen) << " G " << GenParticlePt->at(igen) << " " << GenParticlePhi->at(igen) << " " << GenParticleEta->at(igen) << " " << GenParticleStatus->at(igen) << " " <<  GenParticlePdgId->at(GenParticleMotherIndex->at(igen)) << endl;
          }
          //End debugging

        }
        else muonRejColl.push_back(muonPOGColl[i]);
        match = false;
        dRtmp = 1;
        deltaR = 1;
        gen = -1;
      }

      if((muonGenColl.size() + muonRejColl.size()) != muonPOGColl.size()) {
        cout << "Not all particles accounted for!" << endl;
        return;
      }
      for(UInt_t i = 0; i < muonRejColl.size(); i++) {
        h_muonsRej->Fill(weight, (Int_t) muonRejColl.size(), muonRejColl[i].lorentzVec().Pt(), muonRejColl[i].eta(), muonRejColl[i].lorentzVec().Phi());
      }
    }
    else
      muonGenColl = muonPOGColl;

    // Fake Rates
    std::vector<Lepton> muonLooseButNOTightColl;
    MuonLooseButNOTight.SetPt(15);
    MuonLooseButNOTight.SetEta(2.4);
    MuonLooseButNOTight.SetRelIso(0.12,0.80);
    //MuonLooseButNOTight.SetRelIso(0.05,0.40);
    MuonLooseButNOTight.SetChiNdof(10,50);
    MuonLooseButNOTight.SetBSdxy(0.20, 0.20);
    //MuonLooseButNOTight.SetBSdxy(0.005,0.20);
    MuonLooseButNOTight.SetBSdz(0.50);
    //MuonLooseButNOTight.SetBSdz(0.10);
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
    Jets.JetSelectionLeptonVeto(*PFJetPileupjetIDpassLooseWP, *PFJetEta, *PFJetPhi, *PFJetPt, *PFJetEnergy, *PFJetNeutralEmEnergyFraction, *PFJetNeutralHadronEnergyFraction, *PFJetChargedEmEnergyFraction, *PFJetChargedHadronEnergyFraction, *PFJetChargedMultiplicity, *PFJetNConstituents, *PFJetCombinedSecondaryVertexBTag, *PFJetClosestVertexWeighted3DSeparation, electronColl, muonGenColl, jetColl);

    if (debug) cout<<"matching trigger"<<endl;

/*
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
      if(muonLooseColl[triggerMu1].lorentzVec().Pt()<20) {
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
    //else
      //continue;
    if (MC_pu)
      weight*=triggerweight;
*/

    if (debug) cout<<"trigger matched"<<endl;

    muonLooseColl[0].lorentzVec().Pt()<300. ? ptMu0=muonLooseColl[0].lorentzVec().Pt() : ptMu0=299.;
    muonLooseColl[1].lorentzVec().Pt()<300. ? ptMu1=muonLooseColl[1].lorentzVec().Pt() : ptMu1=299.;
    ID_weight_0 = ID_Iso->GetBinContent(ID_Iso->GetXaxis()->FindBin(fabs(muonLooseColl[0].eta())),ID_Iso->GetYaxis()->FindBin(ptMu0));
    ID_weight_1 = ID_Iso->GetBinContent(ID_Iso->GetXaxis()->FindBin(fabs(muonLooseColl[1].eta())),ID_Iso->GetYaxis()->FindBin(ptMu1));

    if (MC_pu && 0) //FIXME
      weight*=ID_weight_0*ID_weight_1;

    if (debug) cout<<"Iso and ID weights applied"<<endl;

    std::vector<Lepton> electronVetoColl;
    ElectronVeto.SetPt(15);
    ElectronVeto.SetEta(2.5);
    ElectronVeto.SetRelIso(0.15);
    ElectronVeto.SetBSdxy(0.01);
    ElectronVeto.SetBSdz(0.10);
    ElectronVeto.ElectronSelection(*ElectronIsEB, *ElectronIsEE, *ElectronHasTrackerDrivenSeed, *ElectronHasEcalDrivenSeed, *ElectronEta, *ElectronPhi, *ElectronPt, *ElectronEnergy, *ElectronPFPhotonIso03, *ElectronPFNeutralHadronIso03, *ElectronPFChargedHadronIso03, *ElectronCharge, *ElectronGsfCtfScPixCharge, *ElectronMissingHitsEG, *ElectronHasMatchedConvPhot, *ElectronDeltaEtaTrkSC, *ElectronDeltaPhiTrkSC, *ElectronSigmaIEtaIEta, *ElectronHoE, *ElectronCaloEnergy, *ElectronESuperClusterOverP, *ElectronTrackVx, *ElectronTrackVy, *ElectronTrackVz, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), rhoJets, electronVetoColl);


    //checking if we got a good muon and some trigger
    
    //if(muonLooseColl.size()<1) continue;
    for(UInt_t i=0; i<muonGenColl.size(); i++)
      h_noCuts->Fill(weight, (Int_t) muonGenColl.size(), muonGenColl[i].lorentzVec().Pt(), muonGenColl[i].eta(), muonGenColl[i].lorentzVec().Phi());
    //Matching for singleIso trigger
    bool singleIso = false;
    if (debug) cout << "starting single Iso trigger matching" << endl;
    for (UInt_t i=0; i<1; i++) {
      index=muonLooseColl[i].ilepton();
      if (MuonHLTSingleIsoMuonMatched->at(index) && muonLooseColl[i].lorentzVec().Pt()>30.) {
        singleIso = true;
        break;
      }
    }

    if (debug) cout << "single Iso trigger done" << endl;
    if (!singleIso) continue;
    //bool singleIso = (triggerweight > 0);

    // filling standard plots for muons, electrons and jets
    if (muonGenColl.size() > 0) {
      for (UInt_t i=0; i<muonGenColl.size(); i++) {
        index=muonGenColl[i].ilepton();
        h_muonsPOG->Fill(weight, (Int_t) muonGenColl.size(), muonGenColl[i].lorentzVec().Pt(), muonGenColl[i].eta(), muonGenColl[i].lorentzVec().Phi(), muonGenColl[i].charge(), MuonTrkIso->at(index), MuonEcalIso->at(index), MuonHcalIso->at(index), MuonEcalVetoIso->at(index), MuonHcalVetoIso->at(index), MuonPFIsoR03Photon->at(index), MuonPFIsoR03ChargedHadron->at(index), MuonPFIsoR03NeutralHadron->at(index), muonGenColl[i].chiNdof(), muonGenColl[i].dxy_BS(), muonGenColl[i].dz_BS(), MuonPFIsoR03PU->at(index), rhoJets);
      }
    }

    std::vector<Jet> jetVetoColl;
    JetsVeto.SetPt(20);
    JetsVeto.SetEta(2.5);
    JetsVeto.JetSelectionLeptonVeto(*PFJetPileupjetIDpassLooseWP, *PFJetEta, *PFJetPhi, *PFJetPt, *PFJetEnergy, *PFJetNeutralEmEnergyFraction, *PFJetNeutralHadronEnergyFraction, *PFJetChargedEmEnergyFraction, *PFJetChargedHadronEnergyFraction, *PFJetChargedMultiplicity, *PFJetNConstituents, *PFJetCombinedSecondaryVertexBTag, *PFJetClosestVertexWeighted3DSeparation, electronColl, muonLooseColl, jetVetoColl);

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

    if (electronVetoColl.size() > 0) {
      for (UInt_t i=0; i<electronVetoColl.size(); i++) {
        index=electronVetoColl[i].ilepton();
        h_electrons->Fill(weight, (Int_t) electronVetoColl.size(), electronVetoColl[i].lorentzVec().Pt(), electronVetoColl[i].eta(), electronVetoColl[i].lorentzVec().Phi(), electronVetoColl[i].charge(), ElectronTrkIsoDR03->at(index), ElectronEcalIsoDR03->at(index), ElectronHcalIsoDR03->at(index), electronVetoColl[i].dxy_BS(), electronVetoColl[i].dz_BS(), rhoJets);
      }
    }

    b_found = false;

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

    bool POG = ((muonGenColl.size() > 0) && singleIso);
    if(POG) {
    h_muons->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
    //h_muons->StdPlots::Fill(weight, muonPOGcoll.size(), muonPOGColl[i].lorentzVec().Pt(), muonPOGColl[i].eta(), muonPOGColl[i].lorentzVec().Phi());
    h_muons->Fill(weight, muonGenColl);
    h_muons->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    //h_muons->Fill(weight, goodVerticies.size());
    }

    //if (HT > 0)
      //h_HT->Fill(HT, weight);
    h_MET->Fill(PFMETType01XYCor->at(0), weight);
    h_PFSumET->Fill(PFSumETType01XYCor->at(0), weight);

    MET = PFMETType01XYCor->at(0);

    // Background WAS here

    //now we require events to have two muons

    if (electronVetoColl.size()>0 || muonLooseColl.size()!=2) continue;
    if ( (muonLooseColl[0].lorentzVec()+muonLooseColl[1].lorentzVec()).M() < 20) continue;
    METcut = 50.;
    METcontrol = 50.;

    for (UInt_t i=0; i<muonGenColl.size(); i++) {
      index=muonGenColl[i].ilepton();
      h_muonsPOG2->Fill(weight, (Int_t) muonGenColl.size(), muonGenColl[i].lorentzVec().Pt(), muonGenColl[0].eta(), muonGenColl[i].lorentzVec().Phi(), muonGenColl[i].charge(), MuonTrkIso->at(index), MuonEcalIso->at(index), MuonHcalIso->at(index), MuonEcalVetoIso->at(index), MuonHcalVetoIso->at(index), MuonPFIsoR03Photon->at(index), MuonPFIsoR03ChargedHadron->at(index), MuonPFIsoR03NeutralHadron->at(index), muonGenColl[i].chiNdof(), muonGenColl[i].dxy_BS(), muonGenColl[i].dz_BS(), MuonPFIsoR03PU->at(index), rhoJets);
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
  
    bool twoPOG = ((muonGenColl.size() == 2) && POG);
    if(twoPOG) {
    h_twoMu->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
    TLorentzVector s = muonGenColl[0].lorentzVec() + muonGenColl[1].lorentzVec();
    h_twoMu->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge());
    h_twoMu->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_twoMu->Fill(weight, muonGenColl);
   
    h_ChargeFlip->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge(), muonGenColl);
    //h_ChargeFlip->Fill(weight, PFMETPhiType01XYCor->at(0),  muonPOGColl[0].lorentzVec().Phi());
    //h_ChargeFlip->Fill(weight, PFMETPhiType01XYCor->at(0),  muonPOGColl[1].lorentzVec().Phi());
    h_ChargeFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[0].lorentzVec().Pt(), muonGenColl[0].lorentzVec().Phi());
    h_ChargeFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[1].lorentzVec().Pt(), muonGenColl[1].lorentzVec().Phi());
    }

    tag = -1;
    double pt0 = fabs(muonLooseColl[0].lorentzVec().Pt()-48.);
    double pt1 = fabs(muonLooseColl[1].lorentzVec().Pt()-48.);
    if ( pt0 <= 10. && pt1 <= 10.) {
      if(pt0 > pt1) tag = 1;
      if(pt0 < pt1) tag = 0;
    }
    else if ( fabs(muonLooseColl[0].lorentzVec().Pt()-48.) <= 10. )
      tag = 0;
    else if ( fabs(muonLooseColl[1].lorentzVec().Pt()-48.) <= 10. )
      tag = 1;

    POGtag = -1;
    if(twoPOG) {
      pt0 = fabs(muonGenColl[0].lorentzVec().Pt()-48.);
      pt1 = fabs(muonGenColl[1].lorentzVec().Pt()-48.);
      if ( pt0 <= 10. && pt1 <= 10.) {
        if(pt0 > pt1) POGtag = 1;
        if(pt0 < pt1) POGtag = 0;
      }
      else if ( fabs(muonGenColl[0].lorentzVec().Pt()-48.) <= 10. )
        POGtag = 0;
      else if ( fabs(muonGenColl[1].lorentzVec().Pt()-48.) <= 10. )
        POGtag = 1;
    }
/*
    if (((muonGenColl[0].lorentzVec().Pt() > 20 && tag == 1) || (muonGenColl[1].lorentzVec().Pt() > 20 && tag == 0)) && PFMETType01XYCor->at(0) < 30) { //muonGenColl[0].lorentzVec().GetFirstMother == 23 || muonGenColl[0].lorentzVec().GetFirstMother == 22
      h_ChargeFlipCuts->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge(), muonGenColl);
      h_ChargeFlipCuts->Fill(weight, PFMETPhiType01XYCor->at(0),  muonGenColl[0].lorentzVec().Phi());
      h_ChargeFlipCuts->Fill(weight, PFMETPhiType01XYCor->at(0),  muonGenColl[1].lorentzVec().Phi());
    }
*/

    /*
    if(muonGenColl[0].charge() > 0)
      h_muonPos->Fill(muonGenColl[0].eta(),muonGenColl[0].lorentzVec().Pt(),1);
    else
      h_muonNeg->Fill(muonGenColl[0].eta(),muonGenColl[0].lorentzVec().Pt(),1);
    if(muonGenColl[1].charge() > 0)
      h_muonPos->Fill(muonGenColl[1].eta(),muonGenColl[1].lorentzVec().Pt(),1);
    else
      h_muonNeg->Fill(muonGenColl[1].eta(),muonGenColl[1].lorentzVec().Pt(),1);
    */
    if(twoPOG) {
    if(muonGenColl[0].charge()*muonGenColl[1].charge() > 0) {
      h_muonPos->Fill(muonGenColl[0].eta(),1);
      h_muonPos->Fill(muonGenColl[1].eta(),1);
    }
    else {
      h_muonNeg->Fill(muonGenColl[0].eta(),1);
      h_muonNeg->Fill(muonGenColl[1].eta(),1);
    }
    }

    h_METsign->Fill(PFMETType01XYCor->at(0), weight);
    //h_PFSumET_two->Fill(PFSumETType01XYCor->at(0), weight);

    //Check pT = 48+/-10
    //bool ptRange = false;
    //if( fabs(muonGenColl[0].lorentzVec().Pt()-48.) <= 10. ) ptRange = true;
    //else if( fabs(muonGenColl[1].lorentzVec().Pt()-48.) <= 10. ) ptRange = true;
    bool ptRange = ((POGtag != -1) && twoPOG);
    
    //if (!ptRange) continue;

    HT = 0;
    if (jetColl.size() > 0 && ptRange) {
      for (UInt_t i=0; i<jetColl.size(); i++) 
        HT += jetColl[i].lorentzVec().Pt();
    }

    if(ptRange) {
    TLorentzVector s = muonGenColl[0].lorentzVec() + muonGenColl[1].lorentzVec();
    h_pt->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
    h_pt->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_pt->Fill(weight, muonGenColl);
    h_pt->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge());
    /*
    Int_t tag = -1;
    if ( fabs(muonGenColl[0].lorentzVec().Pt()-48.) <= 10. )
      tag = 0;
    else if ( fabs(muonGenColl[1].lorentzVec().Pt()-48.) <= 10. )
      tag = 1;
    */
    h_TagFlip->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge(), muonGenColl);
    if(POGtag == 1)
    //h_TagFlip->Fill(weight, PFMETPhiType01XYCor->at(0),  muonGenColl[0].lorentzVec().Phi());
    h_TagFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[0].lorentzVec().Pt(), muonGenColl[0].lorentzVec().Phi());
    if(POGtag == 0)
    //h_TagFlip->Fill(weight, PFMETPhiType01XYCor->at(0),  muonGenColl[1].lorentzVec().Phi());
    h_TagFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[1].lorentzVec().Pt(), muonGenColl[1].lorentzVec().Phi());
    }

    //Check PFSumET
    //bool PFSumETRange = false;
    //if ( (PFSumETType01XYCor->at(0)-(muonLooseColl[0].lorentzVec().Et()+muonLooseColl[1].lorentzVec().Et())) < 400) PFSumETRange = true;

    //if(!PFSumETRange) continue;

    HT = 0;
    if (jetColl.size() > 0 && ptRange) {
      for (UInt_t i=0; i<jetColl.size(); i++) 
        HT += jetColl[i].lorentzVec().Pt();
    }


    if(ptRange) {
    //h_muonCharge->Fill(muonGenColl[0].charge()*muonGenColl[1].charge(), weight);
    h_nEvents->Fill(fabs(muonGenColl[1].eta()),muonGenColl[1].lorentzVec().Pt(), weight);
    h_PFRange->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
    h_PFRange->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_PFRange->Fill(weight, muonGenColl);
    }

    //MET cuts
    //bool METRange = false;
    //if( PFMETType01XYCor->at(0) < 54 ) METRange = true; //FoM for pure DY
    bool METRange = ((PFMETType01XYCor->at(0) < 54) && ptRange);
    //FIXME 99

    //if(!METRange) continue;

    if(METRange) {
    h_METRange->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
    h_METRange->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_METRange->Fill(weight, muonGenColl, POGtag);
    }

    //HT cuts
    /*bool HTRange = false;
    //if( HT < 30) HTRange = true;
    if( (HT < 26) && (jetColl.size() == 0)) HTRange = true;

    if(!HTRange) continue;

    h_HTRange->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
    s = muonGenColl[0].lorentzVec() + muonGenColl[1].lorentzVec();
    h_HTRange->Fill(weight, s.M());
    h_HTRange->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);

    //N jets cuts
    bool NJetsRange = false;
    if(jetColl.size() <=2) NJetsRange = true;

    if(!NJetsRange) continue;

    h_Njets->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
    s = muonGenColl[0].lorentzVec() + muonGenColl[1].lorentzVec();
    h_Njets->Fill(weight, s.M());
    h_Njets->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    */

    //No jets
    //bool NoJets = false;
    //if( jetColl.size() == 0 ) NoJets = true;
    bool NoJets = ((jetColl.size() == 0) && METRange);
    //bool NoJets = (METRange);

    //if(!NoJets) continue;

    if(NoJets) {
    TLorentzVector s = muonGenColl[0].lorentzVec() + muonGenColl[1].lorentzVec();
    h_NoJets->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
    s = muonGenColl[0].lorentzVec() + muonGenColl[1].lorentzVec();
    h_NoJets->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge());
    h_NoJets->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
    h_NoJets->Fill(weight, muonGenColl);
    
    // Plot same sign
    if(muonGenColl[0].charge()*muonGenColl[1].charge() > 0) {
      h_NoJets_SS->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
      h_NoJets_SS->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge());
      h_NoJets_SS->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
      h_NoJets_SS->Fill(weight, muonGenColl, POGtag);
      if(POGtag == 1)
        h_NoJets_Flip_SS->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[0].lorentzVec().Pt(), muonGenColl[0].lorentzVec().Phi());
      if(POGtag == 0)
        h_NoJets_Flip_SS->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[1].lorentzVec().Pt(), muonGenColl[1].lorentzVec().Phi());
    }
  
    s = muonGenColl[0].lorentzVec() + muonGenColl[1].lorentzVec();
    // Plot opposite sign
    if(muonGenColl[0].charge()*muonGenColl[1].charge() < 0) {
      h_NoJets_OS->Fill(weight, muonGenColl[0].charge()*muonGenColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonGenColl[0].lorentzVec().Et()+muonGenColl[1].lorentzVec().Et()), HT, muonGenColl[0].eta());
      h_NoJets_OS->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge());
      h_NoJets_OS->SetVertex(weight, *VertexNDF, *VertexIsFake, *VertexX, *VertexY, *VertexZ);
      h_NoJets_OS->Fill(weight, muonGenColl, POGtag);
      if(POGtag == 1)
        h_NoJets_Flip_OS->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[0].lorentzVec().Pt(), muonGenColl[0].lorentzVec().Phi());
      if(POGtag == 0)
        h_NoJets_Flip_OS->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[1].lorentzVec().Pt(), muonGenColl[1].lorentzVec().Phi());
    }

    s = muonGenColl[0].lorentzVec() + muonGenColl[1].lorentzVec();
    h_DYFlip->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge(), muonGenColl, POGtag);
    h_DYFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[0].lorentzVec().Pt(), muonGenColl[0].lorentzVec().Phi());
    h_DYFlip->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[1].lorentzVec().Pt(), muonGenColl[1].lorentzVec().Phi());
      if(POGtag == 1) {
        h_DYFlipTag->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge(), muonGenColl, POGtag);
        h_DYFlipTag->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[0].lorentzVec().Pt(), muonGenColl[0].lorentzVec().Phi());
      }
      if(POGtag == 0) {
        h_DYFlipTag->Fill(weight, s.M(), muonGenColl[0].charge()*muonGenColl[1].charge(), muonGenColl, POGtag);
        h_DYFlipTag->Fill(weight, PFMETType01XYCor->at(0), PFMETPhiType01XYCor->at(0), muonGenColl[1].lorentzVec().Pt(), muonGenColl[1].lorentzVec().Phi());
      }

/*
      if(POGtag == 1) {
        h_DY_pt_T->Fill(muonGenColl[0].lorentzVec().Pt());
        if(muonGenColl[0].charge()*muonGenColl[1].charge() > 0)
          h_DY_pt_SS->Fill(muonGenColl[0].lorentzVec().Pt());
      }
      if(POGtag == 0) {
        h_DY_pt_T->Fill(muonGenColl[1].lorentzVec().Pt());
        if(muonGenColl[0].charge()*muonGenColl[1].charge() > 0)
          h_DY_pt_SS->Fill(muonGenColl[1].lorentzVec().Pt());
      }
*/
    for(int num = 0; num <= 1; num++) {
      if(POGtag == num) continue;
      h_DY_pt_T->Fill(muonGenColl[num].lorentzVec().Pt(), weight);
      if(muonGenColl[0].charge()*muonGenColl[1].charge() > 0)
        h_DY_pt_SS->Fill(muonGenColl[num].lorentzVec().Pt(), weight);
    }
    }

    ///// BACKGROUND /////
    if (debug) cout<<"Background selection"<<endl;

    DOUBLEFAKE=false;
    Wcand_tmp=Wcand=0;
    if (muonLooseButNOTightColl.size() == 2 && muonPOGColl.size() == 0 && jetVetoColl.size() >= 2) {
      for(UInt_t i=0; i<muonLooseButNOTightColl.size()-1; i++)
        for(UInt_t j=i+1; j<muonLooseButNOTightColl.size(); j++) {
          if (muonLooseButNOTightColl[i].charge()*muonLooseButNOTightColl[j].charge()>0)
            if (muonLooseButNOTightColl[i].lorentzVec().Pt() >=20) {
              if (debug) cout << "DOUBLE" << endl;
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

    bool DYcuts = singleIso && (tag != -1) && (PFMETType01XYCor->at(0) < 54) && (jetColl.size() == 0);
    if (muonLooseButNOTightColl.size() == 2 && muonPOGColl.size() == 0 && 0) {
      if(singleIso) cout << "Single Iso Double passed" << endl;
      if(tag !=-1) cout << "tag Double passed" << endl;
      if(PFMETType01XYCor->at(0) < 54) cout << "MET Double passed" << endl;
      if(jetColl.size() == 0) cout << "NoJets Double passed" << endl;
    }
    DOUBLEFAKEDY=false;
    if (muonLooseButNOTightColl.size() == 2 && muonPOGColl.size() == 0 && DYcuts) {
      //cout << "Enter Double DY" << endl;
      for(UInt_t i=0; i<muonLooseButNOTightColl.size()-1; i++)
        for(UInt_t j=i+1; j<muonLooseButNOTightColl.size(); j++) {
          if (muonLooseButNOTightColl[i].charge()*muonLooseButNOTightColl[j].charge()>0) {
            //cout << "Same sign Double DY" << endl;
            if (muonLooseButNOTightColl[i].lorentzVec().Pt() >=20) {
              //cout << "DOUBLE DY" << endl;
              dataType=1;
              DOUBLEFAKEDY=true;
              DoubleFakeDY=DoublebackGround(FRhisto, muonLooseButNOTightColl, i, j, doubleFakeDY, dataType, 1);
              Single_DoubleDY=DoubleTOSinglebkg(FRhisto, muonLooseButNOTightColl, i, j);
            }
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
    if (DOUBLEFAKEDY) {
      if (debug) cout<<"        Double found"<<endl;
      nDoubleFakeDY++;
      //h_nVertex0->Fill(numberVertices, weight*DoubleFake);
      h_doublefakesTOTDY->Fill(MET, PFMETPhiType01XYCor->at(0), muonLooseButNOTightColl, DoubleFakeDY*weight, false, false);
      h_totalfakesTOTDY->Fill(MET, PFMETPhiType01XYCor->at(0), muonLooseButNOTightColl, (DoubleFakeDY+Single_DoubleDY)*weight, false, false);
      /*
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
          */
//          h_nVertex0->Fill(numberVertices, weight);
            h_doublefakesDY->Fill(MET, PFMETPhiType01XYCor->at(0), muonLooseButNOTightColl, DoubleFakeDY*weight, false, false);
            h_totalfakesDY->Fill(MET, PFMETPhiType01XYCor->at(0), muonLooseButNOTightColl, (DoubleFakeDY+Single_DoubleDY)*weight, false, false);
          /*
          }
        }
      }
      */
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

    if (muonLooseButNOTightColl.size() == 1 && muonPOGColl.size() == 1 && 0) {
      if(singleIso) cout << "Single Iso Single passed" << endl;
      if(tag !=-1) cout << "tag Single passed" << endl;
      if(PFMETType01XYCor->at(0) < 54) cout << "MET Single passed" << endl;
      if(jetColl.size() == 0) cout << "NoJets Single passed" << endl;
  }
    SINGLEFAKEDY=false;
    Wcand_tmp=Wcand=0;
    if (muonLooseButNOTightColl.size() == 1 && muonPOGColl.size() == 1 && DYcuts) {
      //cout << "Enter Single DY" << endl;
      for(UInt_t i=0; i<muonPOGColl.size(); i++)
        for(UInt_t j=0; j<muonLooseButNOTightColl.size(); j++) {
          if (muonLooseButNOTightColl[j].charge()*muonPOGColl[i].charge()>0) {
            //cout << "Same sign Single DY" << endl;
            if (muonLooseButNOTightColl[j].lorentzVec().Pt() >=20 || muonPOGColl[i].lorentzVec().Pt() >=20) {
              if(muonLooseButNOTightColl[j].ilepton() == muonPOGColl[i].ilepton()) continue;
              //cout << "SINGLE DY" << endl;
              if (debug) cout<<"             Single found"<<endl;
              dataType=1;
              SINGLEFAKEDY=true;
              SingleFakeDY=SinglebackGround(FRhisto, muonLooseButNOTightColl, j, singleFakeDY, dataType, 1);
              DoubleANDSinglebkg(muonPOGColl, i, muonLooseButNOTightColl, j, doubleANDsingleFakeDY, dataType);
              goto endSingleDY;
            }
          }
        }
    }
  endSingleDY:


    if (debug) cout<<"Single done, single = "<<SingleFake<<endl;

    if (SINGLEFAKE) {
      if (debug) cout << "Singlefake" << endl;
      nSingleFake++;
      h_nVertex1->Fill(numberVertices, weight*SingleFake);
      h_singlefakesTOT->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
      h_totalfakesTOT->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
      if (debug) cout << "Plots filled" << endl;
      if (b_found) {
        if (debug) cout << "b Found" << endl;
        h_singlefakesbTag->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
        h_totalfakesbTag->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
      }
      else {
        if (debug) cout << "Check MET" << endl;
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
    if (debug) cout<<"Single done, single = "<<SingleFakeDY<<endl;

    if (SINGLEFAKEDY) {
      if (debug) cout << "Singlefake" << endl;
      nSingleFakeDY++;
      //h_nVertex1->Fill(numberVertices, weight*SingleFake);
      h_singlefakesTOTDY->Fill(MET, PFMETPhiType01XYCor->at(0), muonPOGColl, muonLooseButNOTightColl, SingleFakeDY*weight, false, false);
      h_totalfakesTOTDY->Fill(MET, PFMETPhiType01XYCor->at(0), muonPOGColl, muonLooseButNOTightColl, SingleFakeDY*weight, false, false);
      if (debug) cout << "Plots filled" << endl;
      /*
      if (b_found) {
        if (debug) cout << "b Found" << endl;
        h_singlefakesbTag->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
        h_totalfakesbTag->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
      }
      else {
        if (debug) cout << "Check MET" << endl;
        if (MET>METcontrol) {
          h_singlefakesMET50->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
          h_totalfakesMET50->Fill(MET, muonPOGColl, muonLooseButNOTightColl, jetVetoColl, SingleFake*weight, false, false);
        }
        else {
          if (MET<=METcut) {
          */
            h_singlefakesDY->Fill(MET, PFMETPhiType01XYCor->at(0), muonPOGColl, muonLooseButNOTightColl, SingleFakeDY*weight, false, false);
            h_totalfakes->Fill(MET, PFMETPhiType01XYCor->at(0), muonPOGColl, muonLooseButNOTightColl, SingleFakeDY*weight, false, false);
            //h_nVertex1->Fill(numberVertices, weight);
          /*
          }
        }
      }
      */
    }

    /// BACKGROUND END ///

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
    muonPOGColl.clear();  muonLooseButNOTightColl.clear(); muonLooseColl.clear();
    electronColl.clear(); electronVetoColl.clear(); jetColl.clear();
    genTightColl.clear(); //genBColl.clear();
    if (debug) cout<<"exiting loop"<<endl;
  }
  if (debug) cout<< "out of the loop" <<endl;


  BackGroundEstimate(FRhisto, singleFake, doubleANDsingleFake, doubleFake, finalbkg1, finalbkgerror1, finalbkg2, finalbkgerror2, realsingle, realsingleerror, realdouble, realtotal, doubletosingle, errdoubletosingle, nSplit);

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
  cout<<"Double Fake "<<nDoubleFake<<" value "<<DoubleFake<<endl;
  cout<<"Single_Double "<<Single_Double<<endl;
  cout<<"totale "<<realtotal[0]<<", of which : "<<realsingle[0]<<" single and "<<realdouble[0]<<" double"<<endl;

  BackGroundEstimate(FRhisto, singleFakeDY, doubleANDsingleFakeDY, doubleFakeDY, finalbkg1DY, finalbkgerror1DY, finalbkg2DY, finalbkgerror2DY, realsingleDY, realsingleerrorDY, realdoubleDY, realtotalDY, doubletosingleDY, errdoubletosingleDY, nSplitDY);

  cout<<"Single Fake DY n "<<nSingleFakeDY<<" value "<<SingleFakeDY<<endl;
  for (UInt_t z=0; z<nSplitDY; z++) {
    h_singlefakeDY->SetBinContent(z+1,3,finalbkg1DY[z]);
    h_singlefakeDY->SetBinError(z+1,3,finalbkgerror1DY[z]);
    h_singlefakeDY->SetBinContent(z+1,1,realsingleDY[z]);
    h_singlefakeDY->SetBinError(z+1,1,realsingleerrorDY[z]);
    h_doublefakeDY->SetBinContent(z+1,3,finalbkg2DY[z]);
    h_doublefakeDY->SetBinError(z+1,3,finalbkgerror2DY[z]);
    h_doublefakeDY->SetBinContent(z+1,1,realdoubleDY[z]);
    h_doublefakeDY->SetBinContent(z+1,2,doubletosingleDY[z]);
    h_doublefakeDY->SetBinError(z+1,2,errdoubletosingleDY[z]);
  }
  cout<<"Double Fake DY "<<nDoubleFakeDY<<" value "<<DoubleFakeDY<<endl;
  cout<<"Single_Double DY "<<Single_DoubleDY<<endl;
  cout<<"totale DY "<<realtotalDY[0]<<", of which : "<<realsingleDY[0]<<" single and "<<realdoubleDY[0]<<" double"<<endl;

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
  h_genFlip->Write();
  h_TV->Write();
  h_TV2->Write();

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
  h_singlefakesDY->Write();
  h_singlefakesMET50->Write();
  h_singlefakesbTag->Write();
  h_singlefakesTOT->Write();
  h_singlefakesTOTDY->Write();
  outfile->cd();
  Dir = outfile->mkdir("DoubleFakes");
  outfile->cd( Dir->GetName() );
  h_doublefakes->Write();
  h_doublefakesDY->Write();
  h_doublefakesMET50->Write();
  h_doublefakesbTag->Write();
  h_doublefakesTOT->Write();
  h_doublefakesTOTDY->Write();
  outfile->cd();
  Dir = outfile->mkdir("TotalFakes");
  outfile->cd( Dir->GetName() );
  h_totalfakes->Write();
  h_totalfakesDY->Write();
  h_totalfakesMET50->Write();
  h_totalfakesbTag->Write();
  h_totalfakesTOT->Write();
  h_totalfakesTOTDY->Write();
  outfile->cd();

  Dir = outfile->mkdir("debug");
  outfile->cd( Dir->GetName() );
  FRhisto->Write();
  h_muonsLoose->Write();
  h_LnotT->Write();
  h_jets_veto->Write();
  h_signal->Write();
  h_muonsRej->Write();


  h_nEvents->Write();

  outfile->Close();

}
