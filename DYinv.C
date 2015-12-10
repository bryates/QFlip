#include <vector>
#include "TString.h"
#include "TFile.h"
#include "THStack.h"

void DYinv(){

vector<TString> name_str;
vector<TFile*> FR;
vector<TFile*> Gen;
TH1F *h_FR;
TH1F *h_Gen;
TH1F *h_flip;


name_str.push_back("DY_10-50_");
//name_str.push_back("DY_50_");
//MC_str.push_back("ttbar_"+version+".root");
/*
MC_str.push_back("ZZ_inclusive_"+version+".root");
MC_str.push_back("WZ_inclusive_"+version+".root");
//MC_str.push_back("WW_inclusive_"+version+".root");
MC_str.push_back("WpWp_"+version+".root");
MC_str.push_back("WmWm_"+version+".root");
MC_str.push_back("ttZ_"+version+".root");
MC_str.push_back("ttW_"+version+".root");
MC_str.push_back("ttWW_"+version+".root");
MC_str.push_back("ttV_"+version+".root");
MC_str.push_back("WWW_"+version+".root");
MC_str.push_back("ZZZ_"+version+".root");
MC_str.push_back("WWZ_"+version+".root");
MC_str.push_back("WZZ_"+version+".root");
MC_str.push_back("VVV_"+version+".root");
#MC_str.push_back("W_jets_"+version+".root");
#MC_str.push_back("QCD_mu15_"+version+".root");
*/

for(int i = 0; i < name_str.size(); i++) {
  TFile *f_tmp = new TFile(name_str.at(i)+"FR.root", "READ");
  FR.push_back(f_tmp);
}

for(int i = 0; i < name_str.size(); i++) {
  TFile *f_tmp = new TFile(name_str.at(i)+"Gen.root", "READ");
  Gen.push_back(f_tmp);
}

h_FR = (TH1F*)FR.at(0)->Get("NoJets_OS/h_NoJets_OS_invpt_flip");
h_Gen = (TH1F*)Gen.at(0)->Get("NoJets_OS/h_NoJets_OS_invpt_flip");

for(int i = 1; i < Gen.size(); i++) {
  TH1F *h_tmp = (TH1F*)Gen.at(i)->Get("NoJets_OS/h_NoJets_OS_invpt_flip");
  h_Gen->Add(h_tmp,1);
  
  h_FR->Add((TH1F*)FR.at(i)->Get("NoJets_OS/h_NoJets_OS_invpt_flip"),1);
}

TCanvas c;

h_flip = (TH1F*)h_Gen->Clone();
h_flip->Divide(h_FR);
h_flip->SetMinimum(0);
h_flip->Draw();
gStyle->SetOptStat(0);
c.SaveAs("Charge_Filp_Rate_DY_invpt.pdf");
h_flip->SetAxisRange(10,10e3, "X");
//gPad->SetLogx();
gPad->SetLogy(0);
TGaxis::SetMaxDigits(2);
c.SaveAs("Charge_Filp_Rate_DY_invpt_log.pdf");
c.SaveAs("Charge_Filp_Rate_DY_invpt_log.eps");

TFile *out = new TFile("SingleMu_FR.root","UPDATE");
h_flip->Write("h_Flip_Rate_DY_invpt",TObject::kOverwrite);
out->Close();


}
