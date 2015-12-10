#include "TFile.h"
#include "TString.h"

void fakeClone(TString name="FR") {

TFile *mu = new TFile("SingleMu_"+name+".root", "READ");
TFile *fake = new TFile("Fakes_"+name+".root", "RECREATE");
TH1F *h_NoJets_SS_pt = (TH1F*)mu->Get("TotalFakes/h_secondMuonPt_tf_TOTDY");
TH1F *h_NoJets_SS_tag_pt = (TH1F*)mu->Get("TotalFakes/h_leadingMuonPt_tf_TOTDY");
TH1F *h_NoJets_SS_phi = (TH1F*)mu->Get("TotalFakes/h_secondMuonPhi_tf_TOTDY");
TH1F *h_NoJets_SS_tag_phi = (TH1F*)mu->Get("TotalFakes/h_leadingMuonPhi_tf_TOTDY");
TH1F *h_NoJets_SS_eta = (TH1F*)mu->Get("TotalFakes/h_secondMuonEta_tf_TOTDY");
TH1F *h_NoJets_SS_tag_eta = (TH1F*)mu->Get("TotalFakes/h_leadingMuonEta_tf_TOTDY");
TH1F *h_NoJets_SS_mass = (TH1F*)mu->Get("TotalFakes/h_llmass_tf_TOTDY");
TH1F *h_NoJets_SS_MET = (TH1F*)mu->Get("TotalFakes/h_MET_tf_TOTDY");
TH1F *h_NoJets_SS_invpt = (TH1F*)mu->Get("TotalFakes/h_inv_leadingMuonPt_tf_TOTDY");
TH1F *h_NoJets_SS_invpT = (TH1F*)mu->Get("TotalFakes/h_tf_TOTDY_invpT");
TH1F *h_NoJets_SS_CRAFT = (TH1F*)mu->Get("TotalFakes/h_tf_TOTDY_CRAFT");
TH1F *h_NoJets_SS_CRAFT_tag = (TH1F*)mu->Get("TotalFakes/h_tf_TOTDY_CRAFT_tag");
TH1F *h_NoJets_OS_pt = (TH1F*)mu->Get("TotalFakes/h_secondMuonPt_tf_TOTDY");
TH1F *h_NoJets_OS_tag_pt = (TH1F*)mu->Get("TotalFakes/h_leadingMuonPt_tf_TOTDY");
TH1F *h_NoJets_OS_phi = (TH1F*)mu->Get("TotalFakes/h_secondMuonPhi_tf_TOTDY");
TH1F *h_NoJets_OS_tag_phi = (TH1F*)mu->Get("TotalFakes/h_leadingMuonPhi_tf_TOTDY");
TH1F *h_NoJets_OS_eta = (TH1F*)mu->Get("TotalFakes/h_secondMuonEta_tf_TOTDY");
TH1F *h_NoJets_OS_tag_eta = (TH1F*)mu->Get("TotalFakes/h_leadingMuonEta_tf_TOTDY");
TH1F *h_NoJets_OS_mass = (TH1F*)mu->Get("TotalFakes/h_llmass_tf_TOTDY");
TH1F *h_NoJets_OS_MET = (TH1F*)mu->Get("TotalFakes/h_MET_tf_TOTDY");
TH1F *h_NoJets_SS_invpt = (TH1F*)mu->Get("TotalFakes/h_inv_secondMuonPt_tf_TOTDY");
TH1F *h_NoJets_OS_CRAFT = (TH1F*)mu->Get("TotalFakes/h_tf_TOTDY_CRAFT");
TH1F *h_NoJets_OS_CRAFT_tag = (TH1F*)mu->Get("TotalFakes/h_tf_TOTDY_CRAFT_tag");

TH1F *h_NoJets_Flip_SS_deltaphi = (TH1F*)mu->Get("TotalFakes/h_dMETphisecond_tf_TOTDY");
TH1F *h_NoJets_Flip_SS_cosdphi = (TH1F*)mu->Get("TotalFakes/h_cosdMETphisecond_tf_TOTDY");
TH1F *h_NoJets_Flip_OS_deltaphi = (TH1F*)mu->Get("TotalFakes/h_dMETphilead_tf_TOTDY");
TH1F *h_NoJets_Flip_OS_cosdphi = (TH1F*)mu->Get("TotalFakes/h_cosdMETphilead_tf_TOTDY");

fake->mkdir("NoJets_SS");
fake->cd("NoJets_SS");
h_NoJets_SS_pt->SetName("h_NoJets_SS_pt");
h_NoJets_SS_pt->Write();
h_NoJets_SS_CRAFT->SetName("h_NoJets_SS_CRAFT");
h_NoJets_SS_CRAFT->Write();
h_NoJets_SS_CRAFT_tag->SetName("h_NoJets_SS_CRAFT_tag");
h_NoJets_SS_CRAFT_tag->Write();
h_NoJets_SS_tag_pt->SetName("h_NoJets_SS_tag_pt");
h_NoJets_SS_tag_pt->Write();
h_NoJets_SS_mass->SetName("h_NoJets_SS_mass");
h_NoJets_SS_mass->Write();
h_NoJets_SS_MET->SetName("h_NoJets_SS_MET");
h_NoJets_SS_MET->Write();
h_NoJets_SS_eta->SetName("h_NoJets_SS_eta");
h_NoJets_SS_eta->Write();
h_NoJets_SS_phi->SetName("h_NoJets_SS_phi");
h_NoJets_SS_phi->Write();
h_NoJets_SS_tag_eta->SetName("h_NoJets_SS_tag_eta");
h_NoJets_SS_tag_eta->Write();
h_NoJets_SS_tag_phi->SetName("h_NoJets_SS_tag_phi");
h_NoJets_SS_tag_phi->Write();
h_NoJets_SS_invpT->SetName("h_NoJets_SS_invpT");
h_NoJets_SS_invpT->Write();

fake->mkdir("NoJets_OS");
fake->cd("NoJets_OS");
h_NoJets_OS_pt->SetName("h_NoJets_OS_pt");
h_NoJets_OS_pt->Write();
h_NoJets_OS_CRAFT->SetName("h_NoJets_OS_CRAFT");
h_NoJets_OS_CRAFT->Write();
h_NoJets_OS_CRAFT_tag->SetName("h_NoJets_OS_CRAFT_tag");
h_NoJets_OS_CRAFT_tag->Write();
h_NoJets_OS_tag_pt->SetName("h_NoJets_OS_tag_pt");
h_NoJets_OS_tag_pt->Write();
h_NoJets_OS_mass->SetName("h_NoJets_OS_mass");
h_NoJets_OS_mass->Write();
h_NoJets_OS_MET->SetName("h_NoJets_OS_MET");
h_NoJets_OS_MET->Write();
h_NoJets_OS_eta->SetName("h_NoJets_OS_eta");
h_NoJets_OS_eta->Write();
h_NoJets_OS_phi->SetName("h_NoJets_OS_phi");
h_NoJets_OS_phi->Write();
h_NoJets_OS_tag_eta->SetName("h_NoJets_OS_tag_eta");
h_NoJets_OS_tag_eta->Write();
h_NoJets_OS_tag_phi->SetName("h_NoJets_OS_tag_phi");
h_NoJets_OS_tag_phi->Write();

fake->mkdir("NoJets");
fake->cd("NoJets");
TH1F *h_NoJets_pt = (TH1F*)h_NoJets_OS_pt->Clone("h_NoJets_pt");
h_NoJets_pt->Add(h_NoJets_SS_pt);
h_NoJets_pt->Write();
TH1F *h_NoJets_CRAFT = (TH1F*)h_NoJets_OS_CRAFT->Clone("h_NoJets_CRAFT");
h_NoJets_CRAFT->Add(h_NoJets_SS_CRAFT);
h_NoJets_CRAFT->Write();
TH1F *h_NoJets_CRAFT_tag = (TH1F*)h_NoJets_OS_CRAFT_tag->Clone("h_NoJets_CRAFT_tag");
h_NoJets_CRAFT_tag->Add(h_NoJets_SS_CRAFT_tag);
h_NoJets_CRAFT_tag->Write();
TH1F *h_NoJets_tag_pt = (TH1F*)h_NoJets_OS_tag_pt->Clone("h_NoJets_tag_pt");
h_NoJets_tag_pt->Add(h_NoJets_SS_tag_pt);
h_NoJets_tag_pt->Write();
TH1F *h_NoJets_mass = (TH1F*)h_NoJets_OS_mass->Clone("h_NoJets_mass");
h_NoJets_mass->Add(h_NoJets_SS_mass);
h_NoJets_mass->Write();
TH1F *h_NoJets_MET = (TH1F*)h_NoJets_OS_MET->Clone("h_NoJets_MET");
h_NoJets_MET->Add(h_NoJets_SS_MET);
h_NoJets_MET->Write();
TH1F *h_NoJets_eta = (TH1F*)h_NoJets_OS_eta->Clone("h_NoJets_eta");
h_NoJets_eta->Add(h_NoJets_SS_eta);
h_NoJets_eta->Write();
TH1F *h_NoJets_phi = (TH1F*)h_NoJets_OS_phi->Clone("h_NoJets_phi");
h_NoJets_phi->Add(h_NoJets_SS_phi);
h_NoJets_phi->Write();
TH1F *h_NoJets_tag_eta = (TH1F*)h_NoJets_OS_tag_eta->Clone("h_NoJets_tag_eta");
h_NoJets_tag_eta->Add(h_NoJets_SS_tag_eta);
h_NoJets_tag_eta->Write();
TH1F *h_NoJets_tag_phi = (TH1F*)h_NoJets_OS_tag_phi->Clone("h_NoJets_tag_phi");
h_NoJets_tag_phi->Add(h_NoJets_SS_tag_phi);
h_NoJets_tag_phi->Write();

fake->mkdir("NoJets_Flip_SS");
fake->cd("NoJets_Flip_SS");
h_NoJets_Flip_SS_deltaphi->SetName("h_NoJets_Flip_SS_deltaphi");
h_NoJets_Flip_SS_deltaphi->Write();
h_NoJets_Flip_SS_cosdphi->SetName("h_NoJets_Flip_SS_cosdphi");
h_NoJets_Flip_SS_cosdphi->Write();

fake->mkdir("NoJets_Flip_OS");
fake->cd("NoJets_Flip_OS");
h_NoJets_Flip_OS_deltaphi->SetName("h_NoJets_Flip_OS_deltaphi");
h_NoJets_Flip_OS_deltaphi->Write();
h_NoJets_Flip_OS_cosdphi->SetName("h_NoJets_Flip_OS_cosdphi");
h_NoJets_Flip_OS_cosdphi->Write();

}
