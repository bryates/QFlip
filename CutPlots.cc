#include "CutPlots.h"
#include "TMath.h"
#include "WeightChargeFlip.C"

CutPlots::CutPlots(TString name) : StdPlots(name) {
  debug = true;
  noweight = 1;

  h_muonCharge     = new TH1F ("h_"+name+"_muonCharge", "Muon Charge "+name+";Charge;Events",5,-3.0,3.0);
  h_charge.push_back( new TH1F ("h_"+name+"_muonCharge_m0", "Muon Charge "+name+";Charge;Events",5,-3.0,3.0) );
  h_charge.push_back( new TH1F ("h_"+name+"_muonCharge_m50", "Muon Charge "+name+" (M > 50);Charge;Events",5,-3.0,3.0) );
  h_charge.push_back( new TH1F ("h_"+name+"_muonCharge_mZ", "Muon Charge "+name+" (M = M_{Z}#pm10);Charge;Events",5,-3.0,3.0) );
  h_muonCharge_OS  = new TH1F ("h_"+name+"_muonCharge_OS", "Muon Charge "+name+"(Opposite Sign);Charge;Events",5,-3.0,3.0);
  h_muonCharge_SS  = new TH1F ("h_"+name+"_muonCharge_SS", "Muon Charge "+name+"(Same Sign);Charge;Events",5,-3.0,3.0);
  h_muonCharge_barrel     = new TH1F ("h_"+name+"_muonCharge_barrel", "Muon Charge "+name+" |#eta| < 1.5;Charge;Events",50,-3.0,3.0);
  h_muonCharge_disk     = new TH1F ("h_"+name+"_muonCharge_disk", "Muon Charge "+name+" 1.5 < |#eta| < 2.5;Charge;Events",50,-3.0,3.0);
  h_mass           = new TH1F ("h_"+name+"_mass", "Invariant mass "+name+";Mass (GeV/c^{2});",100,0,500);
  h_MET            = new TH1F ("h_"+name+"_MET", "Missing E_{T} "+name+";E_{T} (GeV);Events",300,0.0,150.0);//300
  h_PFSumET        = new TH1F ("h_"+name+"_PFSumET", "Sum E_{T} "+name+";E_{T} (GeV};Events",500,0.0,3000.0);//5000
  h_PFSumETMinusMu = new TH1F ("h_"+name+"_PFSumETMinusMu", "Sum E_{T} - Muon E_{T} "+name+";E_{T} (GeV};Events",500,0.0,3000.0);//5000
  h_HT = new TH1F ("h_"+name+"_HT", "Sum P_{T} of jets;P_{T} (GeV/c);Events", 100,0.0, 1000.0);//5000
  h_nvtx           = new TH1F ("h_"+name+"_nvtx", "Nvtx "+name+";Number of vertices;",60,0.0,60.0);
  h_muonCharge_nw     = new TH1F ("h_"+name+"_muonCharge_nw", "Muon Charge "+name+" (no weight);Charge;Events",50,-3.0,3.0);
  h_mass_nw           = new TH1F ("h_"+name+"_mass_nw", "Invariant mass "+name+" (no weight);Mass (GeV/c^{2});",100,0,500);
  h_MET_nw            = new TH1F ("h_"+name+"_MET_nw", "Missing E_{T} "+name+" (no weight);E_{T} (GeV);Events",300,0.0,150.0);//300
  h_PFSumET_nw        = new TH1F ("h_"+name+"_PFSumET_nw", "Sum E_{T} "+name+" (no weight);E_{T} (GeV};Events",500,0.0,3000.0);//5000
  h_PFSumETMinusMu_nw = new TH1F ("h_"+name+"_PFSumETMinusMu_nw", "Sum E_{T} - Muon E_{T} "+name+" (no weight);E_{T} (GeV};Events",500,0.0,3000.0);//5000
  h_HT_nw             = new TH1F ("h_"+name+"_HT_nw", "Sum P_{T} of jets;P_{T} (GeV/c) (no weight);Events", 100,0.0, 1000.0);//5000
  h_nvtx_nw           = new TH1F ("h_"+name+"_nvtx_nw", "Nvtx "+name+" (no weight);Number of vertices;",60,0.0,60.0);
  h_mass_50           = new TH1F("h_"+name+"_mass_50", "Mass_{50} "+name+";Mass (GeV/c^{2});",100,0,500);
  h_mass_Z           = new TH1F("h_"+name+"_mass_Z", "Mass_{Z} "+name+";Mass (GeV/c^{2});",100,0,500);
  h_mass_weight           = new TH1F("h_"+name+"_mass_weight", "Mass_{prediction} "+name+";Mass (GeV/c^{2});",100,0,500);
  h_mass_weight_Z           = new TH1F("h_"+name+"_mass_weight_Z", "Mass_{Z prediction} "+name+";Mass (GeV/c^{2});",100,0,500);
  h_invpt_flip        = new TH2F("h_"+name+"_invpt_flip",  name+" 1/p_{T} (GeV)", 15,0,0.05,5,-3.0,3.0);

if(0) {

  h_pt           = new TH1F ("h_"+name+"_pt", name+" p_{T}"";Events;",300,0.0,500.0);
  //h_invpt           = new TH1F ("h_"+name+"_invpt", name+" 1/p_{T}"";Events;",15,0,0.05);
  h_eta           = new TH1F ("h_"+name+"_eta", name+" #eta"";Events;",5,0.0,2.4);
  h_phi           = new TH1F ("h_"+name+"_phi", name+" #phi"";Events;",6,0.0,6.28);
  }
}

CutPlots::~CutPlots() {
  delete h_mass;
  delete h_muonCharge;
  delete h_muonCharge_OS;
  delete h_muonCharge_SS;
  delete h_muonCharge_barrel;
  delete h_muonCharge_disk;
  delete h_MET;
  delete h_PFSumET;
  delete h_PFSumETMinusMu;
  delete h_HT;
  delete h_nvtx;
  delete h_mass_nw;
  delete h_muonCharge_nw;
  delete h_MET_nw;
  delete h_PFSumET_nw;
  delete h_PFSumETMinusMu_nw;
  delete h_HT_nw;
  delete h_mass_50;
  delete h_mass_Z;
  delete h_mass_weight;
  delete h_mass_weight_Z;
  //delete h_invpt;
  delete h_invpt_flip;
  for(unsigned int i = 0; i < h_charge.size(); i++)
    delete h_charge[i];
  if(debug && 0) {
    delete h_pt;
    delete h_eta;
    delete h_phi;
    delete h_probe_pt;
    delete h_probe_eta;
    delete h_probe_phi;
  }
}

void CutPlots::Fill(Double_t weight, Int_t muonCharge, Double_t MET, Double_t PFSumET, Double_t PFSumETMinusMu, Double_t HT, Double_t eta) { // Main cuts
  h_muonCharge->Fill(muonCharge, weight);
  if(fabs(eta) < 1.5)
    h_muonCharge_barrel->Fill(muonCharge, weight);
  else if(fabs(eta) < 2.5)
    h_muonCharge_disk->Fill(muonCharge, weight);
  if(muonCharge > 0)
    h_muonCharge_SS->Fill(muonCharge, weight);
  if(muonCharge < 0)
    h_muonCharge_OS->Fill(muonCharge, weight);
  h_MET->Fill(MET, weight);
  h_PFSumET->Fill(PFSumET, weight);
  h_PFSumETMinusMu->Fill(PFSumETMinusMu, weight);
  h_HT->Fill(HT, weight);
  h_muonCharge_nw->Fill(muonCharge, noweight);
  h_MET_nw->Fill(MET, noweight);
  h_PFSumET_nw->Fill(PFSumET, noweight);
  h_PFSumETMinusMu_nw->Fill(PFSumETMinusMu, noweight);
  h_HT_nw->Fill(HT, noweight);
}

void CutPlots::NoWeight(Double_t weight) {
  noweight = weight;
}

void CutPlots::Fill(Double_t weight, Double_t mass, Int_t charge) { // Mass
  h_mass->Fill(mass, weight);
  h_mass_nw->Fill(mass, noweight);
  if (mass >= 50) {
    h_mass_50->Fill(mass, weight);
    h_charge[1]->Fill(charge, weight);
  }
  if (fabs(mass - 91) < 10) {
    h_mass_Z->Fill(mass, weight);
    h_charge[2]->Fill(charge, weight);
  }
    h_charge[0]->Fill(charge, weight);
}

void CutPlots::Fill(Double_t weight, std::vector<Lepton> &Coll, Int_t tag, Int_t pdgid) { // kinematics 
  h_particles->Fill((Int_t) Coll.size(), weight);
  Double_t W = WeightChargeFlip(Coll[0].lorentzVec().Pt(), Coll[1].lorentzVec().Pt());
  if(isinf(W)) return;
  Double_t M = (Coll[0].lorentzVec()+Coll[1].lorentzVec()).M();
  if(fabs(Coll[0].eta()) < 0.8 && fabs(Coll[1].eta()) < 0.8) {
    //h_mass_Z->Fill(M, weight);
    h_mass_weight->Fill(M, W);
    if (fabs(M - 91) < 10)
      h_mass_weight->Fill(M, W);
  }
  if(debug) {
    for (UInt_t i=0; i<Coll.size(); i++) {
      //if(fabs(M - 91) < 10 && fabs(Coll[i].eta()) < 0.8)
        h_invpt_flip->Fill(1/(Coll[i].lorentzVec().Pt()),Coll[i].charge()*Coll[i].MatchedCharge(),weight);
      if(tag == i) {
        StdPlots::h_tag_pt->Fill(Coll[i].lorentzVec().Pt(),weight);
        StdPlots::h_tag_eta->Fill(Coll[i].eta(),weight);
        StdPlots::h_tag_phi->Fill(Coll[i].lorentzVec().Phi(),weight);  
        if(fabs(Coll[i].eta()) < 0.8)
          StdPlots::h_CRAFT_tag->Fill(Coll[i].lorentzVec().Pt(),weight);
      }
      else {
        StdPlots::h_pt->Fill(Coll[i].lorentzVec().Pt(),weight);
        StdPlots::h_invpt->Fill(1/(Coll[i].lorentzVec().Pt()),weight);
        StdPlots::h_invpT->Fill(1/(Coll[i].lorentzVec().Pt()),Coll[i].eta(),weight);
        StdPlots::h_eta->Fill(Coll[i].eta(),weight);
        StdPlots::h_phi->Fill(Coll[i].lorentzVec().Phi(),weight);
        if(fabs(Coll[i].eta()) < 0.8)
          StdPlots::h_CRAFT->Fill(Coll[i].lorentzVec().Pt(),weight);
      }
    }
  }
}

void CutPlots::SetVertex(Double_t weight, std::vector<double> &VertexNDF, std::vector<bool> &VertexIsFake, std::vector<double> &VertexX, std::vector<double> &VertexY, std::vector<double> &VertexZ) {
  Bool_t *goodVerticies = new Bool_t[VertexNDF.size()];
  unsigned int nvtx = 0;
  while(nvtx<VertexNDF.size()) {
    if(!isGoodEvent((Int_t)VertexNDF.size(), VertexIsFake, VertexNDF, VertexX, VertexY, VertexZ, goodVerticies)) continue;
    nvtx++;
  }
  h_nvtx->Fill(nvtx, weight);
  h_nvtx_nw->Fill(nvtx, noweight);
}

void CutPlots:: Write() {
  h_mass->Write();
  h_muonCharge->Write();
  h_muonCharge_disk->Write();
  h_muonCharge_barrel->Write();
  h_MET->Write();
  h_PFSumET->Write();
  h_PFSumETMinusMu->Write();
  h_HT->Write();
  h_mass_nw->Write();
  h_muonCharge_nw->Write();
  h_MET_nw->Write();
  h_PFSumET_nw->Write();
  h_PFSumETMinusMu_nw->Write();
  h_nvtx->Write();
  h_nvtx_nw->Write();
  h_particles->Write();
  h_mass_50->Write();
  h_mass_Z->Write();
  h_mass_weight->Write();
  for(unsigned int i = 0; i < h_charge.size(); i++)
    h_charge[i]->Write();
  StdPlots::Write();
  h_invpt_flip->Write();
/*
  if(debug) {
    h_pt->Write();
    h_eta->Write();
    h_phi->Write();
  }
*/
}
