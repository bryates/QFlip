#include <vector>
#include "TString.h"
#include "TFile.h"
#include "THStack.h"

void MCinv(TString version="FR"){

vector<TFile*> bgSS;
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

TH2F *tmp_SS = (TH2F*)MC.at(0)->Get("NoJets_SS/h_NoJets_SS_invpT");
h_bgSS = (TH1F*)tmp_SS->ProjectionX("h_NoJets_SS_invpt",43,59);
TH2F *tmp_OS = (TH2F*)MC.at(0)->Get("NoJets_OS/h_NoJets_OS_invpT");
h_bgOS = (TH1F*)tmp_OS->ProjectionX("h_NoJets_OS_invpt",43,59);

//hstack = new THStack("hstack", "Signal and Background plots");

for(int i = 1; i < MC.size(); i++) {
  TH2F *h_tmp = (TH2F*)MC.at(i)->Get("NoJets_SS/h_NoJets_SS_invpT");
  h_bgSS->Add(h_tmp->ProjectionX("h_NoJets_SS_invpT",43,59),1);
  h_tmp = (TH2F*)MC.at(i)->Get("NoJets_OS/h_NoJets_OS_invpT");
  h_bgOS->Add(h_tmp->ProjectionX("h_NoJets_OS_invpT",43,59),1);
}

TCanvas c;

//h_bgSS->SetFillColor(kBlue);
h_bgSS->SetLineColor(kBlack);
h_bgSS->SetMarkerColor(kBlack);
h_bgSS->SetMarkerColor(kBlack);

h_bgSS->Sumw2();
h_bgOS->Sumw2();

h_flip = (TH1F*)h_bgSS->Clone();
h_flip->Divide(h_bgSS,h_bgOS,1,1,"B");
h_flip->Draw();
gStyle->SetOptStat(0);
c.SaveAs("Charge_MC_Filp_Rate_invpt_"+version+".pdf");
h_flip->SetAxisRange(10,10e3, "X");
h_flip->SetMinimum(0);
//gPad->SetLogx();
gPad->SetLogy(0);
TGaxis::SetMaxDigits(2);
c.SaveAs("Charge_MC_Filp_Rate_invpt_log_"+version+".pdf");
c.SaveAs("Charge_MC_Filp_Rate_invpt_log_"+version+".eps");

TFile *out = new TFile("SingleMu_"+version+".root","UPDATE");
h_flip->Write("h_MC_flip_rage_invpt",TObject::kOverwrite);
out->Close();


}
