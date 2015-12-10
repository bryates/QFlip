#ifndef StdPlots_h
#define StdPlots_h

#include <string>
#include "TH1F.h"
#include "TH2F.h"

class StdPlots {

 public:
  TH1F *h_particles, *h_pt, *h_invpt, *h_eta, *h_phi;
  TH1F *h_tag_pt, *h_tag_eta, *h_tag_phi;
  TH1F *h_CRAFT, *h_CRAFT_tag;
  TH2F *h_invpT;
 
  StdPlots(TString name);
  ~StdPlots();
  void Fill(Double_t weight, Int_t N, Double_t pt, Double_t eta, Double_t phi);
  void Write();

};


#endif
