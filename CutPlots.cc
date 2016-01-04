#include "CutPlots.h"
#include "TMath.h"

CutPlots::CutPlots(TString name) : StdPlots(name) {
  debug = true;
  noweight = 1;

  h_muonCharge     = new TH1F ("h_"+name+"_muonCharge", "Muon Charge "+name+";Charge;Events",50,-3.0,3.0);
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

if(0) {

  h_pt           = new TH1F ("h_"+name+"_pt", name+" p_{T}"";Events;",300,0.0,500.0);
  h_eta           = new TH1F ("h_"+name+"_eta", name+" #eta"";Events;",5,0.0,2.4);
  h_phi           = new TH1F ("h_"+name+"_phi", name+" #phi"";Events;",6,0.0,6.28);
  }
}

CutPlots::~CutPlots() {
  delete h_mass;
  delete h_muonCharge;
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
  if(debug) {
    delete h_pt;
    delete h_eta;
    delete h_phi;
  }
}

void CutPlots::Fill(Double_t weight, Int_t muonCharge, Double_t MET, Double_t PFSumET, Double_t PFSumETMinusMu, Double_t HT, Double_t eta) { // Main cuts
  h_muonCharge->Fill(muonCharge, weight);
  if(fabs(eta) < 1.5)
    h_muonCharge_barrel->Fill(muonCharge, weight);
  else if(fabs(eta) < 2.5)
    h_muonCharge_disk->Fill(muonCharge, weight);
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

void CutPlots::Fill(Double_t weight, Double_t mass) { // Mass
  h_mass->Fill(mass, weight);
  h_mass_nw->Fill(mass, noweight);
}

void CutPlots::Fill(Double_t weight, std::vector<Lepton> &Coll, Int_t tag) { // kinematics 
  h_particles->Fill((Int_t) Coll.size(), weight);
  if (Coll.size() == 2)
    h_mass->Fill((Coll[0].lorentzVec()+Coll[1].lorentzVec()).M(), weight);
  for (UInt_t i=0; i<Coll.size(); i++) {
    if (i != tag) {
      h_pt->Fill(Coll[i].lorentzVec().Pt(),weight);
      h_eta->Fill(Coll[i].eta(),weight);
      h_phi->Fill(Coll[i].lorentzVec().Phi(),weight);
    }
    else {
      h_tag_pt->Fill(Coll[i].lorentzVec().Pt(),weight);
      h_tag_eta->Fill(Coll[i].eta(),weight);
      h_tag_phi->Fill(Coll[i].lorentzVec().Phi(),weight);
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
  StdPlots::Write();
  h_mass->Write();
  h_muonCharge->Write();
  h_muonCharge_disk->Write();
  h_muonCharge_barrel->Write();
  h_MET->Write();
  h_PFSumET->Write();
  h_PFSumETMinusMu->Write();
  h_HT->Write();
  //h_mass_nw->Write();
  h_muonCharge_nw->Write();
  //h_MET_nw->Write();
  //h_PFSumET_nw->Write();
  //h_PFSumETMinusMu_nw->Write();
  h_nvtx->Write();
  h_nvtx_nw->Write();
}
