#include "AnalyzerCosmic.h"

AnalyzerCosmic::AnalyzerCosmic() {

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

  h_pt_SS = new TH1F("h_pt_SS",";P_{T} (GeV/C);Events", 100, 0.0, 500.);
  h_pt_T = new TH1F("h_pt_T",";P_{T} (GeV/C);Events", 100, 0.0, 500.);
  h_pt_Flip = new TH1F("h_pt_Flip",";P_{T} (GeV/C);Q Flip Rate", 100, 0.0, 500.);
  h_DY_pt_SS = new TH1F("h_DY_pt_SS",";P_{T} (GeV/C);Events", 100, 0.0, 500.);
  h_DY_pt_T = new TH1F("h_DY_pt_T",";P_{T} (GeV/C);Events", 100, 0.0, 500.);
  h_DY_pt_Flip = new TH1F("h_DY_pt_Flip",";P_{T} (GeV/C);Q Flip Rate", 100, 0.0, 500.);

  h_pt_Flip->Sumw2();
  h_DY_pt_SS->Sumw2();
  h_DY_pt_T->Sumw2();
  h_DY_pt_Flip->Sumw2();

  if (debug) cout<<"inizio"<<endl;

  h_muonsPOG = new MuonPlots("POG_muons");
  h_muonsPOG2 = new MuonPlots("POG_two_muons");
  h_jets = new JetPlots("jets");
  h_jets_2muons = new JetPlots("jets_two_muons");

 
  h_noCuts = new StdPlots("noCuts"); 
  h_muons = new CutPlots("Muons");
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

  if (debug) cout<<"fine"<<endl;
}

AnalyzerCosmic::~AnalyzerCosmic() { }

void AnalyzerCosmic::SetName(TString name, Int_t version) {
  completename = name + "_";
  completename += version;
  completename += ".root";
  outfile = new TFile(completename,"RECREATE");
}

void AnalyzerCosmic::SetWeight(Double_t CrossSection, Double_t nevents) {

  MCweight = integratedlumi * CrossSection / nevents;
  cout<<"MCweight = "<<MCweight<<endl;
 
}

void AnalyzerCosmic::SetEvtN(Long64_t events) {
  events ? entrieslimit=events :  entrieslimit=-1;
  cout<<"events "<<events<<endl<<"entrieslimit "<<entrieslimit<<endl;
}

void AnalyzerCosmic::Loop() {

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

    if ( !TriggerSelector(triggerslist, *HLTInsideDatasetTriggerNames, *HLTInsideDatasetTriggerDecisions, *HLTInsideDatasetTriggerPrescales, prescaler) ) continue;
  
    if (debug) cout<<"trigger passed"<<endl;

    if (MC_pu) {
      /// ***PU reweghting*** ///
      h_nvtx_norw->Fill(PileUpInteractionsTrue->at(0), MCweight);
      weight = reweightPU->GetWeight(PileUpInteractionsTrue->at(0))*MCweight;
      h_nvtx_rw->Fill(PileUpInteractionsTrue->at(0), weight);
    }

    if (debug) cout<<"pileup reweghting applied"<<endl;

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
    
    std::vector<Lepton> muonPOGColl;
    MuonPOG.SetPt(10); 
    //MuonPOG.SetEta(2.4);
    MuonPOG.SetEta(0.8);
    MuonPOG.SetRelIso(0.12);
    MuonPOG.SetChiNdof(10); 
    MuonPOG.SetBSdxy(0.20);
    MuonPOG.SetBSdz(0.50);
    MuonPOG.SetDeposits(400.0,600.0);
    MuonPOG.MuonSelection(*MuonIsPF, *MuonIsGlobal, *MuonEta, *MuonPhi, *MuonPt, *MuonPtError, *MuonEnergy, *MuonPFIsoR04ChargedHadron, *MuonPFIsoR04NeutralHadron, *MuonPFIsoR04Photon, *MuonEcalVetoIso, *MuonHcalVetoIso, *MuonCharge, *MuonGlobalTrkValidHits, *MuonTrkPixelHits, *MuonStationMatches, *MuonTrackLayersWithMeasurement, *MuonGlobalChi2, *MuonTrkVx, *MuonTrkVy, *MuonTrkVz, *MuonTrkD0, *MuonTrkD0Error, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *MuonPFIsoR04PU, muonPOGColl);

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

    for(int num = 0; num <= 1; num++) {
      if(tag == num) continue;
      h_pt_T->Fill(muonPOGColl[num].lorentzVec().Pt());
      if(muonPOGColl[0].charge()*muonPOGColl[1].charge() > 0)
        h_pt_SS->Fill(muonPOGColl[num].lorentzVec().Pt());
    }

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
    //}

    for(int num = 0; num <= 1; num++) {
      if(tag == num) continue;
      h_DY_pt_T->Fill(muonPOGColl[num].lorentzVec().Pt());
      if(muonPOGColl[0].charge()*muonPOGColl[1].charge() > 0)
        h_DY_pt_SS->Fill(muonPOGColl[num].lorentzVec().Pt());
    }

    //h_ChargeFlip->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), muonPOGColl);


  }
  if (debug) cout<< "out of the loop" <<endl;
  
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

   h_pt_Flip->Divide(h_pt_SS,h_pt_T);
   h_pt_Flip->Write();
   h_DY_pt_Flip->Divide(h_DY_pt_SS,h_DY_pt_T);
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

  h_nEvents->Write();

  outfile->Close();

}
