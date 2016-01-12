#include "Analyzer.h"

Analyzer::Analyzer() {

  h_nEvents = new TH2F ("h_nEvents", "Number of Events",ninteta,arrayeta,nintpT,arraypT);
  h_muonCharge = new TH1F ("h_muonCharge", "Muon Charge;Charge;Events",50,-3.0,3.0);
  //h_muonCharge->Sumw2();
  h_MET = new TH1F ("h_MET", "Missing E_{T};E_{T} (GeV);Events",300,0.0,150.0);//300
  h_METsign = new TH1F ("h_METsign", "Missing E_{T} for two muons;E_{T} (GeV);Events",300,0.0,150.0);//300
  h_PFSumET = new TH1F ("h_PFSumET", "Sum E_{T};E_{T} (GeV};Events",500,0.0,3000.0);//5000
  h_HT = new TH1F ("h_HT", "Sum P_{T} of jets;P_{T} (GeV/c);Events", 100,0.0, 500.0);//3000

  h_nvtx_norw = new TH1F("h_nvtx_norw","Nvtx per bunch crossing at BX = 0 noreweight",60,0.0,60.0);
  h_nvtx_rw = new TH1F("h_nvtx_rw","Nvtx per bunch crossing at BX = 0 reweight",60,0.0,60.0);

  if (debug) cout<<"inizio"<<endl;

  h_muonsPOG = new MuonPlots("POG_muons");
  h_muonsPOG2 = new MuonPlots("POG_two_muons");
  h_jets = new JetPlots("jets");
  h_jets_2muons = new JetPlots("jets_two_muons");

  
  h_muons = new CutPlots("muons");
  h_twoMu = new CutPlots("TwoMuons");
  h_singleIso = new CutPlots("singleIso");
  h_pt = new CutPlots("pt");
  h_PFRange = new CutPlots("PFRange");
  h_METRange = new CutPlots("METRange");
  h_NoJets = new CutPlots("NoJets");
  h_GoldenZ = new CutPlots("GoldenZ");

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

    h_nvtx_norw->Fill(VertexNDF->size(), MCweight);

    if (MC_pu) {
      /// ***PU reweghting*** ///
//      h_nvtx_norw->Fill(PileUpInteractionsTrue->at(0), MCweight);
      weight = reweightPU->GetWeight(PileUpInteractionsTrue->at(0))*MCweight;
//      h_nvtx_rw->Fill(PileUpInteractionsTrue->at(0), weight);
    }
    h_nvtx_rw->Fill(VertexNDF->size(), weight);
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
    //MuonPOG.SetRelIso(0.05);
    MuonPOG.SetChiNdof(10); 
    MuonPOG.SetBSdxy(0.20);
    //MuonPOG.SetBSdxy(0.005);
    MuonPOG.SetBSdz(0.50);
    //MuonPOG.SetBSdz(0.10);
    MuonPOG.SetDeposits(400.0,600.0);
    //MuonPOG.SetDeposits(4.0,6.0);
    MuonPOG.MuonSelection(*MuonIsPF, *MuonIsGlobal, *MuonEta, *MuonPhi, *MuonPt, *MuonPtError, *MuonEnergy, *MuonPFIsoR04ChargedHadron, *MuonPFIsoR04NeutralHadron, *MuonPFIsoR04Photon, *MuonEcalVetoIso, *MuonHcalVetoIso, *MuonCharge, *MuonGlobalTrkValidHits, *MuonTrkPixelHits, *MuonStationMatches, *MuonTrackLayersWithMeasurement, *MuonGlobalChi2, *MuonTrkVx, *MuonTrkVy, *MuonTrkVz, *MuonTrkD0, *MuonTrkD0Error, VertexX->at(VertexN), VertexY->at(VertexN), VertexZ->at(VertexN), *MuonPFIsoR04PU, muonPOGColl);

    //checking if we got a good muon and some trigger
    if(muonPOGColl.size()<1) continue;
    //Matching for singleIso trigger
    Bool_t singleIso = false;
    for (UInt_t i=0; i<1; i++) {
      if (MC_pu) {
        index=muonPOGColl[i].ilepton();
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

    if (HT > 0)
      h_HT->Fill(HT, weight);
    h_MET->Fill(PFMETType01XYCor->at(0), weight);
    h_PFSumET->Fill(PFSumETType01XYCor->at(0), weight);

    //now we require events to have two muons

    if (electronColl.size()>0 || muonPOGColl.size()!=2) continue;
    if ( (muonPOGColl[0].lorentzVec()+muonPOGColl[1].lorentzVec()).M() < 20) continue;

    for (UInt_t i=0; i<muonPOGColl.size(); i++) {
      if (jetColl.size() > 0) break;
      index=muonPOGColl[i].ilepton();
      h_muonsPOG2->Fill(weight, (Int_t) muonPOGColl.size(), muonPOGColl[i].lorentzVec().Pt(), muonPOGColl[i].eta(), muonPOGColl[i].lorentzVec().Phi(), muonPOGColl[i].charge(), MuonTrkIso->at(index), MuonEcalIso->at(index), MuonHcalIso->at(index), MuonEcalVetoIso->at(index), MuonHcalVetoIso->at(index), MuonPFIsoR03Photon->at(index), MuonPFIsoR03ChargedHadron->at(index), MuonPFIsoR03NeutralHadron->at(index), muonPOGColl[i].chiNdof(), muonPOGColl[i].dxy_BS(), muonPOGColl[i].dz_BS(), MuonPFIsoR03PU->at(index), rhoJets);
    }
    
    if (jetColl.size() > 0) {
      for (UInt_t i=0; i<jetColl.size(); i++) {
        index=jetColl[i].ijet();
        h_jets_2muons->Fill( weight, (Int_t) jetColl.size(), jetColl[i].lorentzVec().Pt(), jetColl[i].eta(), jetColl[i].lorentzVec().Phi(), PFJetTrackCountingHighPurBTag->at(index), PFJetJetProbabilityBTag->at(index), jetColl[i].btag_disc(), PFJetClosestVertexWeightedXYSeparation->at(index), PFJetClosestVertexWeightedZSeparation->at(index), PFJetClosestVertexWeighted3DSeparation->at(index) );
      }
    }

    h_twoMu->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    //TLorentzVector s = muonPOGColl[0].lorentzVec() + muonPOGColl[1].lorentzVec();
    //h_twoMu->Fill(weight, s.M());
    h_twoMu->Fill(weight, muonPOGColl);

    h_METsign->Fill(PFMETType01XYCor->at(0), weight);
    //h_PFSumET_two->Fill(PFSumETType01XYCor->at(0), weight);

    //Check pT = 44+/-10
    bool ptRange = false;
    //if( fabs(muonPOGColl[0].lorentzVec().Pt()-44.) <= 10. ) ptRange = true;
    //else if ( fabs(muonPOGColl[1].lorentzVec().Pt()-44.) <= 10. ) ptRange = true;
    int tag = -1;
    double pt0 = fabs(muonPOGColl[0].lorentzVec().Pt()-44.);
    double pt1 = fabs(muonPOGColl[1].lorentzVec().Pt()-44.);
    if ( pt0 <= 10. && pt1 <= 10.) {
      if(pt0 > pt1) tag = 1;
      if(pt0 < pt1) tag = 0;
    }
    else if ( pt0 <= 10. )
      tag = 0;
    else if ( pt1 <= 10. )
      tag = 1;

    if (tag != -1) ptRange = true;
    
    if (!ptRange) continue;
    h_pt->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    h_pt->Fill(weight, muonPOGColl, tag);

    //Check PFSumET
/*
    bool PFSumETRange = false;
    if ( (PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et())) < 400) PFSumETRange = true;

    if(!PFSumETRange) continue;
*/

    h_muonCharge->Fill(muonPOGColl[0].charge()+muonPOGColl[1].charge(), weight);
    h_nEvents->Fill(fabs(muonPOGColl[1].eta()),muonPOGColl[1].lorentzVec().Pt(), weight);
    //h_PFRange->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());

    bool METRange = false;
    if (PFMETType01XYCor->at(0) < 26) METRange = true;

    if (!METRange) continue;
    h_METRange->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    h_METRange->Fill(weight, muonPOGColl);

    bool NoJets = false;
    if (jetColl.size() == 0) NoJets = true;

    if (!NoJets) continue;
    h_NoJets->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    h_NoJets->Fill(weight, muonPOGColl, tag);

    if ( fabs((muonPOGColl[0].lorentzVec()+muonPOGColl[1].lorentzVec()).M() - 90) > 10) continue;
    h_GoldenZ->Fill(weight, muonPOGColl[0].charge()*muonPOGColl[1].charge(), PFMETType01XYCor->at(0), PFSumETType01XYCor->at(0), PFSumETType01XYCor->at(0)-(muonPOGColl[0].lorentzVec().Et()+muonPOGColl[1].lorentzVec().Et()), HT, muonPOGColl[0].eta());
    h_GoldenZ->Fill(weight, muonPOGColl, tag);

  }
  if (debug) cout<< "out of the loop" <<endl;
  
  outfile->cd();
  //h_muonCharge->Write();

  h_nvtx_norw->Write();
  h_nvtx_rw->Write();

  Dir = outfile->mkdir("Jets");
  outfile->cd( Dir->GetName() );
  h_jets->Write();
  h_jets_2muons->Write();
  outfile->cd();

  Dir = outfile->mkdir("Muons");
  outfile->cd( Dir->GetName() );
  h_muonsPOG->Write();
  h_MET->Write();
  h_PFSumET->Write();
  h_HT->Write();
  outfile->cd();

  Dir = outfile->mkdir("TwoMuons");
  outfile->cd( Dir->GetName() );
  h_muonsPOG2->Write();
  h_twoMu->Write();
  h_METsign->Write();
  outfile->cd();

  Dir = outfile->mkdir("pt");
  outfile->cd( Dir->GetName() );
  h_pt->Write();
  outfile->cd();

  Dir = outfile->mkdir("METRange");
  outfile->cd( Dir->GetName() );
  h_METRange->Write();
  outfile->cd();

  Dir = outfile->mkdir("NoJets");
  outfile->cd( Dir->GetName() );
  h_NoJets->Write();
  outfile->cd();

  Dir = outfile->mkdir("GoldenZ");
  outfile->cd( Dir->GetName() );
  h_GoldenZ->Write();
  outfile->cd();

  h_nEvents->Write();

  Dir = outfile->mkdir("debug");
  outfile->cd( Dir->GetName() );
  h_PFRange->Write();
  outfile->cd();

  outfile->Close();

}
