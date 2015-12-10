#include "ChargeFlip.h"
#include "WeightChargeFlip.C"

static const Double_t pi = 3.1415926535;

ChargeFlip::ChargeFlip(TString name) : StdPlots(name) {

  h_charge_SS = new TH1F("h_"+name+"_charge_SS","Charge (Same Sign);Events",5,-3.0,3.0);
  h_charge_T = new TH1F("h_"+name+"_charge_T","Charge (Opposite Sign);Events",5,-3.0,3.0);
  h_chargeFakes = new TH1F("h_"+name+"_chargeFlip","CF;Charge;Q Flip Fraction",5,-3.0,3.0);
  h_charge_pt_SS = new TH1F("h_"+name+"_charge_pt_SS","CF_{SS} 1/P_{T};P_{T} (1/GeV);#mu^{#pm} #mu^{#pm}",15,0,0.05);
  h_charge_pt_T = new TH1F("h_"+name+"_charge_pt_T","CF_{T} 1/P_{T};P_{T} (1/GeV);Total",15,0,0.05);
  h_charge_pt = new TH1F("h_"+name+"_charge_pt","CF P_{T};1/P_{T} (1/GeV);Q Flip Fraction",15,0,0.05);
  h_charge_pt_SS_m50 = new TH1F("h_"+name+"_charge_pt_SS_m50","CF_{SS} 1/P_{T} (M > 50);1/P_{T} (1/GeV);#mu^{#pm} #mu^{#pm}",15,0,0.05);
  h_charge_pt_T_m50 = new TH1F("h_"+name+"_charge_pt_T_m50","CF_{T} P_{T} (M > 50;1/P_{T} (1/GeV);Total",15,0,0.05);
  h_charge_pt_m50 = new TH1F("h_"+name+"_charge_pt_m50","CF P_{T} (M > 50);1/P_{T} (1/GeV);Q Flip Fraction",15,0,0.05);
  h_charge_pt_SS_mZ = new TH1F("h_"+name+"_charge_pt_SS_mZ","CF_{SS} P_{T} (M_{Z} #pm 10);1/P_{T} (1/GeV);#mu^{#pm} #mu^{#pm}",15,0,0.05);
  h_charge_pt_T_mZ = new TH1F("h_"+name+"_charge_pt_T_mZ","CF_{T} P_{T} (M_{Z} #pm 10);1/P_{T} (1/GeV);Total",15,0,0.05);
  h_charge_pt_mZ = new TH1F("h_"+name+"_charge_pt_mZ","CF P_{T} (M_{Z} #pm 10);1/P_{T} (1/GeV);Q Flip Fraction",15,0,0.05);
  h_deltaphi = new TH1F("h_"+name+"_deltaphi","#phi_{MET} - #phi_{#mu};#phi (rad);Events",100,0,pi);//-3.1415926535,3.1415926535);
  h_mass_weight = new TH1F("h_"+name+"_mass_weight",";M#mu#mu (GeV);Events",100,0,200);
  h_mass_SS = new TH1F("h_"+name+"_mass_SS",";M#mu#mu (GeV);Events",100,0,200);
  h_mass_T = new TH1F("h_"+name+"_mass_T",";M#mu#mu (GeV);Events",100,0,200);
  h_weight = new TH2F("h_"+name+"_weight","CF P_{T};1/P_{T} (1/GeV);Q Flip Fraction",15,0,0.05,200,0,1);
  h_pt_eta = new TH2F("h_"+name+"_pt_eta",";1/P_{T} (1/GeV);#eta;Q Flip Fraction",15,0,0.05,50,-2.5,2.5);
  h_pt_eta_T = new TH2F("h_"+name+"_pt_eta_T",";1/P_{T} (1/GeV);#eta;Total",15,0,0.05,50,-2.5,2.5);
  h_pt_eta_SS = new TH2F("h_"+name+"_pt_eta_SS",";1/P_{T} (1/GeV);#eta;#mu^{#pm} #mu^{#pm}",15,0,0.05,50,-2.5,2.5);
  //h_METPhi = new TH1F("h_"+name+"_MET_phi",";#phi (rad)",50,-6.283185307,6.283185307);
  //h_muPhi = new TH1F("h_"+name+"_mu_phi",";#phi (rad)",50,-6.283185307,6.283185307);
  h_cosdphi = new TH1F("h_"+name+"_cosdphi","cos(#phi_{MET} - #phi_{#mu});cos(#phi);Events",100,-1,1);//-3.1415926535,3.1415926535);
  //h_invpt           = new TH1F ("h_"+name+"_invpt", name+" 1/p_{T}"";Events;",15,0,0.05);


  h_charge_SS->Sumw2();
  h_charge_T->Sumw2();
  h_charge_pt_SS->Sumw2();
  h_charge_pt_T->Sumw2();
  h_charge_pt_SS_m50->Sumw2();
  h_charge_pt_T_m50->Sumw2();
  h_charge_pt_SS_mZ->Sumw2();
  h_charge_pt_T_mZ->Sumw2();
  h_pt_eta_T->Sumw2();
  h_pt_eta_SS->Sumw2();

}

ChargeFlip::~ChargeFlip() {
  delete h_charge_SS;
  delete h_charge_T;
  delete h_chargeFakes;
  delete h_charge_pt;
  delete h_charge_pt_SS_m50;
  delete h_charge_pt_T_m50;
  delete h_charge_pt_m50;
  delete h_charge_pt_SS_mZ;
  delete h_charge_pt_T_mZ;
  delete h_charge_pt_mZ;
  delete h_mass_weight;
  delete h_mass_SS;
  delete h_mass_T;
  delete h_weight;
  delete h_pt_eta;
  //delete h_METPhi;
  //delete h_muPhi;
  delete h_cosdphi;
  //delete h_invpt;
}

void ChargeFlip::Fill(Double_t weight, Double_t mass, Int_t charge, std::vector<Lepton> &Coll, Int_t tag) { // kinematics
  h_particles->Fill((Int_t) Coll.size(), weight);
  Double_t W = WeightChargeFlip(Coll[0].lorentzVec().Pt(), Coll[1].lorentzVec().Pt());
  if(isinf(W)) return;
  Double_t M = (Coll[0].lorentzVec()+Coll[1].lorentzVec()).M();
  if(fabs(M - 91) <= 20) {
    h_mass_T->Fill(M, weight);
    if(charge < 0) {
      h_mass_weight->Fill(M, W);
    }
    else if(charge > 0) {
      h_mass_SS->Fill(M, weight);
      //h_mass_weight->Fill(M, weight);
    }
  }
  for(Int_t i = 0; i < (Int_t) Coll.size(); i++) {
    if(tag == i) continue;
    h_pt_eta_T->Fill(1/Coll[i].lorentzVec().Pt(), Coll[i].eta(), weight);
    h_weight->Fill(1/Coll[i].lorentzVec().Pt(), weight);
    if(charge > 0) {
      h_charge_pt_SS->Fill(1/Coll[i].lorentzVec().Pt(), weight);
      h_pt_eta_SS->Fill(1/Coll[i].lorentzVec().Pt(), Coll[i].eta(), weight);
    }
    h_charge_pt_T->Fill(1/Coll[i].lorentzVec().Pt(), weight);
    if(mass < 50) continue;
    if(charge > 0)
      h_charge_pt_SS_m50->Fill(1/Coll[i].lorentzVec().Pt(), weight);
    h_charge_pt_T_m50->Fill(1/Coll[i].lorentzVec().Pt(), weight);
    if(fabs(mass - 91) > 10) continue;
    if(charge > 0)
      h_charge_pt_SS_mZ->Fill(1/Coll[i].lorentzVec().Pt(), weight);
    h_charge_pt_T_mZ->Fill(1/Coll[i].lorentzVec().Pt(), weight);
  }

  for (Int_t i=0; i < (Int_t) Coll.size(); i++) {
    if(tag == i) {
      StdPlots::h_tag_pt->Fill(Coll[i].lorentzVec().Pt(), weight);
      StdPlots::h_tag_eta->Fill(Coll[i].eta(), weight);
      StdPlots::h_tag_pt->Fill(Coll[i].lorentzVec().Phi(), weight);
    }
    else {
      StdPlots::h_pt->Fill(Coll[i].lorentzVec().Pt(), weight);
      StdPlots::h_invpt->Fill(1/(Coll[i].lorentzVec().Pt()), weight);
      StdPlots::h_eta->Fill(Coll[i].eta(),weight);
      StdPlots::h_phi->Fill(Coll[i].lorentzVec().Phi(), weight);
    }
  }
}

void ChargeFlip::Fill(Double_t weight, Double_t m, Double_t METphi, Double_t vec, Double_t phi) {
  MET = TVector3(m * cos(METphi),m * sin(METphi), 0);
  v = TVector3(vec * cos(phi), vec * sin(phi), 0);
  h_deltaphi->Fill(MET.Angle(v), weight);
  h_cosdphi->Fill(cos(MET.Angle(v)), weight);
  //h_phi->Fill(phi, weight);
  //h_deltaphi->Fill((METphi - phi), weight);
  //h_cosdphi->Fill(cos(METphi - phi), weight);
  //h_METPhi->Fill(METphi, weight);
  //h_muPhi->Fill(phi, weight);
}

void ChargeFlip::Write() {
  h_particles->Write();
  h_charge_SS->Write();
  h_charge_T->Write();
  h_chargeFakes->Divide(h_charge_SS, h_charge_T);
  h_chargeFakes->Write();
  h_charge_pt_SS->Write();
  h_charge_pt_T->Write();
  h_charge_pt->Divide(h_charge_pt_SS, h_charge_pt_T);
  h_charge_pt->Write();
  h_charge_pt_SS_m50->Write();
  h_charge_pt_T_m50->Write();
  h_charge_pt_SS_mZ->Write();
  h_charge_pt_T_mZ->Write();
  h_charge_pt_m50->Divide(h_charge_pt_SS_m50, h_charge_pt_T_m50);
  h_charge_pt_mZ->Divide(h_charge_pt_SS_mZ, h_charge_pt_T_mZ);
  h_charge_pt_m50->Write();
  h_charge_pt_mZ->Write();
  h_pt_eta->Divide(h_pt_eta_SS, h_pt_eta_T);
  h_pt->Write();
  //h_invpt->Write();
  h_phi->Write();
  h_deltaphi->Write();
  h_mass_weight->Write();
  h_mass_SS->Write();
  h_mass_T->Write();
  h_weight->Write();
  h_pt_eta_T->Write();
  h_pt_eta_SS->Write();
  h_pt_eta->Write();
  //h_METPhi->Write();
  //h_muPhi->Write();
  h_cosdphi->Write();
  StdPlots::Write();
/*
  h_tag_pt->Write();
  h_tag_eta->Write();
  h_tag_phi->Write();
*/
}
