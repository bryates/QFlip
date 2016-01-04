#include <vector>
std::vector<TString> data;
std::vector<TString> bg;
TFile *files,*fileb;
TH1F *h_tmp, *h_s, *h_b;
TString cut = "TwoMuons";
//TString cut = "PFRange";
TString type = "MET";
Float_t xmax = 0.0;
Int_t xbins = 0;
Double_t s = 0.0;
Double_t b = 0.0;

double FoMmax[4]; // (x,s,b,FoM)

void FoM() {
data.push_back("DY_10-50_FR.root");
data.push_back("DY_50_FR.root");
bg.push_back("ttbar_FR.root");
bg.push_back("ZZ_inclusive_FR.root");
bg.push_back("WZ_inclusive_FR.root");
bg.push_back("WW_inclusive_FR.root");
bg.push_back("WpWp_FR.root");
bg.push_back("WmWm_FR.root");
bg.push_back("ttZ_FR.root");
bg.push_back("ttW_FR.root");
bg.push_back("ttWW_FR.root");
bg.push_back("WWW_FR.root");
bg.push_back("ZZZ_FR.root");
bg.push_back("WWZ_FR.root");
bg.push_back("WZZ_FR.root");
//bg.push_back("Wjets_FR.root");
//bg.push_back("QCD_mu15_v3_FR.root");

for(int i=0; i < data.size(); i++) {
  files = new TFile(data[i],"READ");
  if (!files->IsOpen())
    continue;
  files->cd();
  h_tmp = (TH1F*)files->Get( cut + "/h_" + cut + "_" + type );
  //h_tmp = (TH1F*)files->Get( cut + "/h_" + type );
  if (h_tmp->GetNbinsX() > xbins) xbins = h_tmp->GetNbinsX();
  if (h_tmp->GetBinWidth(1)*xbins > xmax) xmax = h_tmp->GetBinWidth(1)*xbins;
  files->Close();
}

//TProfile *h_FoM = new TProfile("h_FoM", cut+"/h_FoM_"+type,xbins,0,xmax);
TH1F *h_FoM = new TH1F("h_FoM", cut+"/h_FoM_"+type,xbins,0,xmax);

h_s = new TH1F("h_s","h_s",xbins,0,xmax);
h_b = new TH1F("h_b","h_b",xbins,0,xmax);
for(int idata=0; idata < data.size(); idata++) {
  files = new TFile(data[idata],"READ");
  if (!files->IsOpen())
    continue;
  files->cd();
  h_tmp = (TH1F*)files->Get( cut + "/h_" + cut + "_" + type );
  //h_tmp = (TH1F*)files->Get( cut + "/h_" + type );
  h_s->Add(h_tmp,1);
  files->Close();
}
for(int ibg=0; ibg < bg.size(); ibg++) {
  fileb = new TFile(bg[ibg],"READ");
  if (!fileb->IsOpen())
    continue;
  fileb->cd();
  h_tmp = (TH1F*)fileb->Get( cut + "/h_" + cut + "_" + type );
  //h_tmp = (TH1F*)fileb->Get( cut + "/h_" + type );
  h_b->Add(h_tmp,1);
  fileb->Close();
}
for(Int_t x=1; x <= xbins; x++) {
  if(1) {
  s = h_s->Integral(1,x);
  b = h_b->Integral(1,x);
  }
  else {
  //hdevMaj->SetBinContent(i,h_th1f[ifile][iplot]->Integral(i,hdevMaj->GetNbinsX())/sqrt(h_th1f[ifile][iplot]->Integral(i,hdevMaj->GetNbinsX())+htotal->Integral(i,hdevMaj->GetNbinsX()) ) );

  s = h_s->Integral(x,xbins);
  b = h_b->Integral(x,xbins);
  //h_FoM->SetBinContent(x,h_s->Integral(x,h_FoM->GetNbinsX())/sqrt(h_b->Integral(x,h_FoM->GetNbinsX())));
  //h_FoM->SetBinContent(x,h_s->Integral(x,h_FoM->GetNbinsX())/sqrt(h_b->Integral(x,h_FoM->GetNbinsX())+h_s->Integral(x,h_FoM->GetNbinsX())));
  }
  if(b == 0) // Avoid zeroes/infinities
    continue;
  double FoM = s / sqrt(b);
  //cout << s << ", " << b << ", " << FoM << endl;
  //h_FoM->Fill(FoM,1);
  h_FoM->SetBinContent(x,FoM);
  //h_FoM->SetBinError(x,0);
  if(FoM > FoMmax[3]) {
    FoMmax[0] = x;
    FoMmax[1] = s;
    FoMmax[2] = b;
    FoMmax[3] = FoM;
  }
}
  h_b->Draw();
  h_s->Draw("same");

cout << "x=" << FoMmax[0] << " ";
cout << "s=" << FoMmax[1] << " ";
cout << "b=" << FoMmax[2] << " ";
cout << "FoM=" <<FoMmax[3] << " " << std::endl;
//h_FoM->Rebin(2);
h_FoM->Draw("L");
}
