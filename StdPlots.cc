#include "StdPlots.h"

StdPlots::StdPlots(TString name) {  

  h_particles = new TH1F("h_N_"+name,      "Number of "+name,   50,0,50);
  h_pt        = new TH1F("h_"+name+"_pt",  name+" p_{T} (GeV)", 60,0,300);
  h_invpt        = new TH1F("h_"+name+"_invpt",  name+" 1/p_{T} (GeV)", 15,0,0.05);
  h_eta       = new TH1F("h_"+name+"_eta", name+" #eta",        100,-5,5);
  h_phi       = new TH1F("h_"+name+"_phi", name+" #phi",        100,-3.1415926535,3.1415926535);
  h_tag_pt        = new TH1F("h_"+name+"_tag_pt",  name+" 1/p_{T} (GeV)", 60,0,300);
  h_tag_eta       = new TH1F("h_"+name+"_tag_eta", name+" #eta",        100,-5,5);
  h_tag_phi       = new TH1F("h_"+name+"_tag_phi", name+" #phi",        100,-3.1415926535,3.1415926535);
  Double_t edges[7] = {10,20,30,50,100,120,1000};
  h_CRAFT        = new TH1F("h_"+name+"_CRAFT",  name+" p_{T} (GeV)", 6, edges);
  h_CRAFT_tag        = new TH1F("h_"+name+"_CRAFT_tag",  name+" p_{T} (GeV)", 6, edges);
  h_invpT        = new TH2F("h_"+name+"_invpT",  name+" 1/p_{T} (GeV)", 25,0,0.05,100,-5,5);

}

StdPlots::~StdPlots() {
  delete h_particles;
  delete h_pt;
  delete h_invpt;
  delete h_invpT;
  delete h_eta;
  delete h_phi;
  delete h_tag_pt;
  delete h_tag_eta;
  delete h_tag_phi;
  delete h_CRAFT;
  delete h_CRAFT_tag;
}

void StdPlots::Fill(Double_t weight, Int_t N, Double_t pt, Double_t eta, Double_t phi) {
  if ( N > 0 )
    h_particles->Fill(N,weight/N);
  else
    h_particles->Fill(N,weight);
  h_pt->Fill(pt,weight);
  h_invpt->Fill(1/pt,weight);
  h_eta->Fill(eta,weight);
  h_phi->Fill(phi,weight);
  h_invpT->Fill(1/pt,eta,weight);
}

void StdPlots::Write() {
 h_particles->Write();
 h_pt->Write();
 h_invpt->Write();
 h_eta->Write();
 h_phi->Write();
 h_tag_pt->Write();
 h_tag_eta->Write();
 h_tag_phi->Write();
 h_CRAFT->Write();
 h_CRAFT_tag->Write();
 h_invpT->Write();
}

