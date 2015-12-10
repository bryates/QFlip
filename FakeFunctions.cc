#include "FakeFunctions.h"

double FailFailToSingle(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate, Double_t eventWeight) {
    
    UInt_t FRbinLep1, FRbinLep2;
    
    if (0==eventWeight) eventWeight=1.0;
    
    FRbinLep1=fakerate->FindBin(Lep1.Eta(),Lep1.Pt());
    FRbinLep2=fakerate->FindBin(Lep2.Eta(),Lep2.Pt());
    
    Double_t Lep1FR, Lep2FR;
    Double_t Lep1FRerr, Lep2FRerr;
    
    Lep1FR=fakerate->GetBinContent(FRbinLep1);
    Lep1FRerr=fakerate->GetBinError(FRbinLep1);
    Lep2FR=fakerate->GetBinContent(FRbinLep2);
    Lep2FRerr=fakerate->GetBinError(FRbinLep2);
    
    weight = (Lep1FR*(1-Lep2FR)+Lep2FR*(1-Lep1FR)) / ((1-Lep1FR)*(1-Lep2FR));
    return weight;
    
    /*
     Double_t epsilon1, epsilon2;
     Double_t epsilon1err, epsilon2err;
     
     epsilon1 = Lep1FR / (1-Lep1FR);
     epsilon2 = Lep2FR / (1-Lep2FR);
     */
}


double FailFail(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate, Double_t eventWeight) {
    
    UInt_t FRbinLep1, FRbinLep2;
    
    if (0==eventWeight) eventWeight=1.0;
  
    FRbinLep1=fakerate->FindBin(fabs(Lep1.Eta()),Lep1.Pt());
cout<< "Eta " << fabs(Lep1.Eta()) << " and pt " << Lep1.Pt() << endl;
    FRbinLep2=fakerate->FindBin(fabs(Lep2.Eta()),Lep2.Pt());
    
    Double_t Lep1FR, Lep2FR;
    Double_t Lep1FRerr, Lep2FRerr;

    Lep1FR=fakerate->GetBinContent(FRbinLep1);
    Lep1FRerr=fakerate->GetBinError(FRbinLep1);
    Lep2FR=fakerate->GetBinContent(FRbinLep2);
    Lep2FRerr=fakerate->GetBinError(FRbinLep2);
    
    weight = - Lep1FR*Lep2FR / ((1-Lep1FR)*(1-Lep2FR));
    return weight;
    
    /*
    Double_t epsilon1, epsilon2;
    Double_t epsilon1err, epsilon2err;
    
    epsilon1 = Lep1FR / (1-Lep1FR);
    epsilon2 = Lep2FR / (1-Lep2FR);
    */
}

double FailPass(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate, Double_t eventWeight) {
    
    UInt_t FRbinLep1, FRbinLep2;
    
    if (0==eventWeight) eventWeight=1.0;
    
    FRbinLep1=fakerate->FindBin(Lep1.Eta(),Lep1.Pt());
    FRbinLep2=fakerate->FindBin(Lep2.Eta(),Lep2.Pt());
    
    Double_t Lep1FR, Lep2FR;
    Double_t Lep1FRerr, Lep2FRerr;
    
    Lep1FR=fakerate->GetBinContent(FRbinLep1);
    Lep1FRerr=fakerate->GetBinError(FRbinLep1);
    Lep2FR=fakerate->GetBinContent(FRbinLep2);
    Lep2FRerr=fakerate->GetBinError(FRbinLep2);
    
    weight = Lep1FR / (1-Lep1FR);
    return weight;
    
    /*
     Double_t epsilon1, epsilon2;
     Double_t epsilon1err, epsilon2err;
     
     epsilon1 = Lep1FR / (1-Lep1FR);
     epsilon2 = Lep2FR / (1-Lep2FR);
     */
}

double PassFail(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate, Double_t eventWeight) {
    
    UInt_t FRbinLep1, FRbinLep2;
    
    if (0==eventWeight) eventWeight=1.0;
    
    FRbinLep1=fakerate->FindBin(Lep1.Eta(),Lep1.Pt());
    FRbinLep2=fakerate->FindBin(Lep2.Eta(),Lep2.Pt());
    
    Double_t Lep1FR, Lep2FR;
    Double_t Lep1FRerr, Lep2FRerr;
    
    Lep1FR=fakerate->GetBinContent(FRbinLep1);
    Lep1FRerr=fakerate->GetBinError(FRbinLep1);
    Lep2FR=fakerate->GetBinContent(FRbinLep2);
    Lep2FRerr=fakerate->GetBinError(FRbinLep2);
    
    weight = Lep2FR / (1-Lep2FR);
    return weight;
    
    /*
     Double_t epsilon1, epsilon2;
     Double_t epsilon1err, epsilon2err;
     
     epsilon1 = Lep1FR / (1-Lep1FR);
     epsilon2 = Lep2FR / (1-Lep2FR);
     */
}

double PassPass(Double_t &weight, Double_t &error, TLorentzVector Lep1, TLorentzVector Lep2, TH2F* fakerate, TH2F* promptrate, Double_t eventWeight) {
    
    UInt_t FRbinLep1, FRbinLep2;
    
    if (0==eventWeight) eventWeight=1.0;
    
    FRbinLep1=fakerate->FindBin(Lep1.Eta(),Lep1.Pt());
    FRbinLep2=fakerate->FindBin(Lep2.Eta(),Lep2.Pt());
    
    Double_t Lep1FR, Lep2FR;
    Double_t Lep1FRerr, Lep2FRerr;
    
    Lep1FR=fakerate->GetBinContent(FRbinLep1);
    Lep1FRerr=fakerate->GetBinError(FRbinLep1);
    Lep2FR=fakerate->GetBinContent(FRbinLep2);
    Lep2FRerr=fakerate->GetBinError(FRbinLep2);
    
    weight = Lep1FR*Lep2FR / ((1-Lep1FR)*(1-Lep2FR));
    return weight;
    
    /*
     Double_t epsilon1, epsilon2;
     Double_t epsilon1err, epsilon2err;
     
     epsilon1 = Lep1FR / (1-Lep1FR);
     epsilon2 = Lep2FR / (1-Lep2FR);
     */
}


