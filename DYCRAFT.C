#include <vector>
#include "TString.h"
#include "TFile.h"
#include "THStack.h"

void DYCRAFT(TString version="FR"){

vector<TString> name_str;
vector<TFile*> Gen;
vector<TFile*> FR;
TH1F *h_GenSS;
TH1F *h_GenOS;
TH1F *h_GenSSBG;
TH1F *h_GenSSFR;
TH1F *h_GenSSTOT;
TH1F *h_FRSS;
TH1F *h_FROS;
TH1F *h_fake;
TH1F *h_flip;
THStack *hstack;


name_str.push_back("SingleMu_");
name_str.push_back("DY_10-50_");
name_str.push_back("DY_50_");

for(int i = 0; i < name_str.size(); i++) {
  TFile *f_tmp = new TFile(name_str.at(i)+"Gen.root", "READ");
  Gen.push_back(f_tmp);
}

for(int i = 0; i < name_str.size(); i++) {
  TFile *f_tmp = new TFile(name_str.at(i)+"FR.root", "READ");
  FR.push_back(f_tmp);
}

h_GenSS = (TH1F*)Gen.at(0)->Get("NoJets_SS/h_NoJets_SS_CRAFT");
h_GenOS = (TH1F*)Gen.at(0)->Get("NoJets_OS/h_NoJets_OS_CRAFT");
h_FRSS = (TH1F*)FR.at(0)->Get("NoJets_SS/h_NoJets_SS_CRAFT");
h_FROS = (TH1F*)FR.at(0)->Get("NoJets_OS/h_NoJets_OS_CRAFT");

hstack = new THStack("hstack", "Signal and Background plots");

for(int i = 1; i < Gen.size(); i++) {
  TH1F *h_tmp = (TH1F*)Gen.at(i)->Get("NoJets_SS/h_NoJets_SS_CRAFT");
  h_GenSS->Add(h_tmp,1);
  h_tmp = (TH1F*)FR.at(i)->Get("NoJets_OS/h_NoJets_OS_CRAFT");
  h_FROS->Add(h_tmp,1);
}

TCanvas c;

//h_GenSS->SetFillColor(kBlue);
h_GenSS->SetLineColor(kBlack);
h_GenSS->SetMarkerColor(kBlack);
h_GenSS->SetMarkerColor(kBlack);

h_flip = (TH1F*)h_GenSS->Clone();
h_flip->Sumw2();
h_flip->Divide(h_FROS);
h_flip->Draw();
gStyle->SetOptStat(0);
gPad->SetLogy();
c.SaveAs("Charge_Filp_Rate_DY_CRAFT.pdf");
h_flip->SetAxisRange(10,10e3, "X");
gPad->SetLogx();
c.SaveAs("Charge_Filp_Rate_DY_log_CRAFT.pdf");
c.SaveAs("Charge_Filp_Rate_DY_log_CRAFT.eps");

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
