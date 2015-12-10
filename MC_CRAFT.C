#include <vector>
#include "TString.h"
#include "TFile.h"
#include "THStack.h"

void MC_CRAFT(TString version="FR"){

vector<TString> MC_str;
vector<TFile*> MC;
TH1F *h_bgSS;
TH1F *h_bgOS;
TH1F *h_fake;
TH1F *h_flip;
THStack *hstack;


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

for(int i = 0; i < MC_str.size(); i++) {
  TFile *f_tmp = new TFile(MC_str.at(i), "READ");
  MC.push_back(f_tmp);
}

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

//h_bgSS->SetFillColor(kBlue);
h_bgSS->SetLineColor(kBlack);
h_bgSS->SetMarkerColor(kBlack);
h_bgSS->SetMarkerColor(kBlack);
h_bgSS->Draw();
c.SaveAs("BG_MC.pdf");
c.SaveAs("BG_MC.eps");

//hstack->Add(h_bgSS);
hstack->Add(h_bgSS);
//hstack->Add(h_bgSSTOT);

TLegend *legend = new TLegend(.65, 0.945 - 4*0.065, .93, 0.91);
legend->SetFillColor(0);
legend->AddEntry(h_bgSS, "MC", "fl");
//legend->AddEntry(h_bgSSTOT, "Signal - MC - fakes", "fl");

gPad->SetLogy();
h_bgSS->Draw();
hstack->Draw("9histsame");
h_bgSS->Draw("same");
legend->Draw("same");
c.SaveAs("hstack_MC.pdf");
c.SaveAs("hstack_MC.eps");

h_flip = (TH1F*)h_bgSS->Clone();
h_flip->Sumw2();
h_flip->Divide(h_bgOS);
h_flip->Draw();
gStyle->SetOptStat(0);
c.SaveAs("Charge_Filp_Rate_MC_CRAFT.pdf");
h_flip->SetAxisRange(10,10e3, "X");
gPad->SetLogx();
c.SaveAs("Charge_Filp_Rate_log_MC_CRAFT.pdf");
c.SaveAs("Charge_Filp_Rate_log_MC_CRAFT.eps");

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
