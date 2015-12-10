#include <vector>
#include "TString.h"
#include "TFile.h"
#include "THStack.h"

void fakeCRAFT(TString version="FR"){

vector<TString> sigSS_str;
vector<TFile*> sigSS;
vector<TString> MC_str;
vector<TFile*> MC;
TH1F *h_sigSS;
TH1F *h_sigOS;
TH1F *h_sigSSBG;
TH1F *h_sigSSFR;
TH1F *h_sigSSTOT;
TH1F *h_bgSS;
TH1F *h_bgOS;
TH1F *h_fake;
TH1F *h_flip;
THStack *hstack;


sigSS_str.push_back("SingleMu_"+version+".root");
//MC_str.push_back("DY_10-50_"+version+".root");
//MC_str.push_back("DY_50_"+version+".root");
//MC_str.push_back("ttbar_"+version+".root");
MC_str.push_back("ZZ_inclusive_"+version+".root");
MC_str.push_back("WZ_inclusive_"+version+".root");
//MC_str.push_back("WW_inclusive_"+version+".root");
MC_str.push_back("WpWp_"+version+".root");
MC_str.push_back("WmWm_"+version+".root");
/*MC_str.push_back("ttZ_"+version+".root");
MC_str.push_back("ttW_"+version+".root");
MC_str.push_back("ttWW_"+version+".root");*/
MC_str.push_back("ttV_"+version+".root");
/*MC_str.push_back("WWW_"+version+".root");
MC_str.push_back("ZZZ_"+version+".root");
MC_str.push_back("WWZ_"+version+".root");
MC_str.push_back("WZZ_"+version+".root");*/
MC_str.push_back("VVV_"+version+".root");
/*#MC_str.push_back("W_jets_"+version+".root");
#MC_str.push_back("QCD_mu15_"+version+".root");*/

for(int i = 0; i < sigSS_str.size(); i++) {
  TFile *f_tmp = new TFile(sigSS_str.at(i), "READ");
  sigSS.push_back(f_tmp);
}

for(int i = 0; i < MC_str.size(); i++) {
  TFile *f_tmp = new TFile(MC_str.at(i), "READ");
  MC.push_back(f_tmp);
}

h_sigSS = (TH1F*)sigSS.at(0)->Get("NoJets_SS/h_NoJets_SS_CRAFT");
h_sigOS = (TH1F*)sigSS.at(0)->Get("NoJets_OS/h_NoJets_OS_CRAFT");
h_fake = (TH1F*)sigSS.at(0)->Get("TotalFakes/h_tf_TOTDY_CRAFT");
h_bgSS = (TH1F*)MC.at(0)->Get("NoJets_SS/h_NoJets_SS_CRAFT");
h_bgOS = (TH1F*)MC.at(0)->Get("NoJets_OS/h_NoJets_OS_CRAFT");

hstack = new THStack("hstack", "Signal and Background plots");

for(int i = 1; i < MC.size(); i++) {
  TH1F *h_tmp = (TH1F*)MC.at(i)->Get("NoJets_SS/h_NoJets_SS_CRAFT");
  h_bgSS->Add(h_tmp,1);
  h_tmp = (TH1F*)MC.at(i)->Get("NoJets_OS/h_NoJets_OS_CRAFT");
  h_bgOS->Add(h_tmp,1);
}

TCanvas c;

//h_sigSS->SetFillColor(kBlue);
h_sigSS->SetLineColor(kBlack);
h_sigSS->SetMarkerColor(kBlack);
h_sigSS->SetMarkerColor(kBlack);

h_sigSSBG = (TH1F*)h_sigSS->Clone();
h_sigSSBG->Add(h_bgSS,-1);
h_sigOS->Add(h_bgOS,-1);
h_sigSS->Sumw2();
h_sigOS->Sumw2();
h_sigSSBG->Sumw2();

h_sigSSBG->Draw();
h_bgSS->SetFillColor(kRed);
h_bgSS->SetLineColor(kBlack);
gStyle->SetOptStat(0);
c.SaveAs("Sig_minus_BG.pdf");
c.SaveAs("Sig_minus_BG.eps");

h_sigSSFR = (TH1F*)h_sigSS->Clone();
h_sigSSFR->Add(h_fake,-1);
h_fake->SetFillColor(kGreen+2);
h_fake->SetLineColor(kBlack);

h_sigSSFR->Draw();
gStyle->SetOptStat(0);
c.SaveAs("Sig_minus_Fake.pdf");
c.SaveAs("Sig_minus_Fake.eps");

h_sigSSTOT = (TH1F*)h_sigSSBG->Clone();
h_sigSSTOT->Sumw2();
h_sigSSTOT->Add(h_fake,-1);
h_sigSSTOT->Draw();
h_sigSSTOT->SetFillColor(kBlue);

//hstack->Add(h_sigSS);
hstack->Add(h_fake);
hstack->Add(h_bgSS);
//hstack->Add(h_sigSSTOT);

gStyle->SetOptStat(0);
c.SaveAs("Sig_minus_All.pdf");
c.SaveAs("Sig_minus_All.eps");
h_sigSS->Draw();

TLegend *legend = new TLegend(.65, 0.945 - 4*0.065, .93, 0.91);
legend->SetFillColor(0);
legend->AddEntry(h_sigSS, "Signal", "ple");
legend->AddEntry(h_bgSS, "MC", "fl");
legend->AddEntry(h_fake, "Fakes", "fl");
//legend->AddEntry(h_sigSSTOT, "Signal - MC - fakes", "fl");

gPad->SetLogy();
h_sigSS->Draw();
hstack->Draw("9histsame");
h_sigSS->Draw("same");
legend->Draw("same");
c.SaveAs("hstack.pdf");
c.SaveAs("hstack.eps");

h_flip = (TH1F*)h_sigSSTOT->Clone();
h_flip->Sumw2();
h_flip->Divide(h_sigOS);
h_flip->Draw();
gStyle->SetOptStat(0);
c.SaveAs("Charge_Filp_Rate_CRAFT.pdf");
h_flip->SetAxisRange(10,10e3, "X");
gPad->SetLogx();
c.SaveAs("Charge_Filp_Rate_log_CRAFT.pdf");
c.SaveAs("Charge_Filp_Rate_log_CRAFT.eps");

/*
TH1F *h_flip_ptinv = new TH1F("h_flip_ptinv", "Q Flip vs 1/P_{T}", h_flip->GetNbinsX(), 0, 0.05);
for(int i = 1; i < h_flip->GetNbinsX(); i++) {
  double num = h_flip->GetBinContent(i);
  h_flip_ptinv->SetBinContent(1./i, num);
}

gPad->SetLogx(0);
h_flip_ptinv->Draw();
gPad->SetLogy();
TGaxis::SetMaxDigits(2);
c.SaveAs("Charge_Filp_Rate_ptinv_log_"+version+".pdf");
c.SaveAs("Charge_Filp_Rate_ptinv_log_"+version+".eps");*/


}
