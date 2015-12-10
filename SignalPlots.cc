#include "SignalPlots.h"

SignalPlots::SignalPlots(TString name) : StdPlots(name) {
  h_jjmass =         new TH1F("h_dijetsmass_"    + name,"Invariant mass of the two leading jets",100,0,1000);
  h_llmass =         new TH1F("h_llmass_"        + name,"Invariant mass of the two leading muons",200,0,1000);
  h_l1jjmass =       new TH1F("h_l1jjmass_"      + name,"Invariant mass of the two leading jets and leading muon",100,0,1000);
  h_l2jjmass =       new TH1F("h_l2jjmass_"      + name,"Invariant mass of the two leading jets and second muon",100,0,1000);
  h_lljjmass =       new TH1F("h_lljjmass_"      + name,"Invariant mass of the four particles",200,0,2000);
  h_WandNmass =      new TH2F("h_WandNmass_"     + name,"Invariant mass of the W and the N",200,0,2000,200,0,2000);
  h_3Dparm =         new TH3F("h_3Dpar_"         + name,"m(lljj) and muon p_{T}_{1} and muon p_{T}_{2}",200,0,2000,60,0,300,40,0,200);
  h_3DparmJ =        new TH3F("h_3DparJ_"        + name,"leading jet p_{T} and muon p_{T}_{1} and muon p_{T}_{2}",60,0,300,60,0,300,40,0,200);
  h_leadingMuonPt =  new TH1F("h_leadingMuonPt_" + name,"leading muon pt",60,0,300);
  h_secondMuonPt =   new TH1F("h_secondMuonPt_"  + name,"secondary muon pt",60,0,300);
  h_leadingJetPt =   new TH1F("h_leadingJetPt_"  + name,"leading jet pt",60,0,300);
  h_secondJetPt =    new TH1F("h_secondJetPt_"   + name,"secondary jet pt",60,0,300);
  h_leadingMuonIso = new TH1F("h_leadingMuonIso_"+ name,"leading muon relIso",40,0,0.4);
  h_secondMuonIso =  new TH1F("h_secondMuonIso_" + name,"secondary muon relIso",40,0,0.4);
  h_MET =            new TH1F("h_MET_"           + name,"Missing Et",300,0.0,150.0);
  h_paircharge =     new TH1F("h_paircharge_"    + name,"Charge of the muon pair",5,-2,3);
  h_muonseta =       new TH1F("h_muonseta_"      + name,"#eta distribution of the two muons",50,-3,3);
  h_jetseta =        new TH1F("h_jetseta_"       + name,"#eta distribution of the two jets",50,-3,3);
  h_bTag =           new TH1F("h_bTag_"          + name,"bTag discrimant",100,-1,3);
  h_cosTheta1 =      new TH1F("h_cosTheta1_"     + name,"cos#theta first muon",100,-1,1);
  h_cosTheta2 =      new TH1F("h_cosTheta2_"     + name,"cos#theta second muon",100,-1,1);
  h_DeltaPhi =       new TH1F("h_DeltaPhi_"      + name,"#Delta#phi between two muons",100,0.0,3.15);
  h_Njets =          new TH1F("h_Njets_"         + name,"number of jets",10,0,10);
  h_dMETphilead = new TH1F("h_dMETphilead_"+name,"#phi_{MET} - #phi_{leading #mu};#phi (rad);Events",100,0,pi);
  h_dMETphisecond = new TH1F("h_dMETphisecond_"+name,"#phi_{MET} - #phi_{second #mu};#phi (rad);Events",100,0,pi);
  h_cosdMETphilead = new TH1F("h_cosdMETphilead_"+name,"cos(#phi_{MET} - #phi_{leading #mu});cos(#phi);Events",100,-1,1);
  h_cosdMETphisecond = new TH1F("h_cosdMETphisecond_"+name,"cos(#phi_{MET} - #phi_{second #mu});cos(#phi);Events",100,-1,1);
  //h_invPt =   new TH1F("h_invPt_"  + name,"secondary muon 1/pt",60,0,300);
  h_leadingMuonEta = new TH1F("h_leadingMuonEta_"+name, " #eta", 100,-5,5);
  h_leadingMuonPhi = new TH1F("h_leadingMuonPhi_"+name, " #phi", 100,-3.1415926535,3.1415926535);
  h_secondMuonEta = new TH1F("h_secondMuonEta_"+name, " #eta", 100,-5,5);
  h_secondMuonPhi = new TH1F("h_secondMuonPhi_"+name, " #phi", 100,-3.1415926535,3.1415926535);
  h_inv_leadingMuonPt = new TH1F("h_inv_leadingMuonPt_"+name,"1/P_{T}",15,0,0.05);
  h_inv_secondMuonPt = new TH1F("h_inv_secondMuonPt_"+name,"1/P_{T}",15,0,0.05);
}

SignalPlots::~SignalPlots() {
  delete h_jjmass;
  delete h_llmass;
  delete h_l1jjmass;
  delete h_l2jjmass;
  delete h_lljjmass;
  delete h_WandNmass;
  delete h_3Dparm;
  delete h_3DparmJ;
  delete h_leadingMuonPt;
  delete h_secondMuonPt;
  delete h_leadingJetPt;
  delete h_secondJetPt;
  delete h_leadingMuonIso;
  delete h_secondMuonIso;
  delete h_MET;
  delete h_paircharge;
  delete h_muonseta;
  delete h_jetseta;
  delete h_bTag;
  delete h_cosTheta1;
  delete h_cosTheta2;
  delete h_DeltaPhi;
  delete h_Njets;
  delete h_dMETphilead;
  delete h_dMETphisecond;
  delete h_cosdMETphilead;
  delete h_cosdMETphisecond;
  //h_invPt->Write();
  //delete h_invPt;
  delete h_leadingMuonPhi;
  delete h_leadingMuonEta;
  delete h_secondMuonPhi;
  delete h_secondMuonEta;
  delete h_inv_leadingMuonPt;
  delete h_inv_secondMuonPt;
}

void SignalPlots::Fill(Double_t MET, std::vector<Lepton>& muons, std::vector<Jet>& jets, Double_t weight, Bool_t ptok, Bool_t ssok) {
  dijetmass_tmp=dijetmass=9999.9;
  UInt_t m,n;
  for(UInt_t i=0; i<muons.size()-1; i++)
    for(UInt_t j=i+1; j<muons.size(); j++) {
      if (muons[i].charge()*muons[j].charge()>0 || ssok)
	if (muons[i].lorentzVec().Pt()>=20 || ptok) {
          for(UInt_t emme=0; emme<jets.size()-1; emme++)
	    for(UInt_t enne=1; enne<jets.size(); enne++) {
              //dijetmass_tmp = (jets[emme].lorentzVec()+jets[enne].lorentzVec()).M();
              dijetmass_tmp = (muons[i].lorentzVec()+muons[j].lorentzVec()+jets[emme].lorentzVec()+jets[enne].lorentzVec()).M();
              if ( fabs(dijetmass_tmp-Mass_W) < fabs(dijetmass-Mass_W) ) {
                dijetmass = dijetmass_tmp;
                m = emme;
                n = enne;
	      }
	    }
          //if (dijetmass > 200 || (jets[m].lorentzVec()+jets[n].lorentzVec()).M()>120) goto nogoodW;
          //if (dijetmass<50 || dijetmass>110) goto nogoodW;
          h_MET->Fill(MET, weight);
	  //m=0; n=1;
	  h_jjmass->Fill( (jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
          h_llmass->Fill( (muons[i].lorentzVec()+muons[j].lorentzVec()).M(),weight);
	  h_l1jjmass->Fill( (muons[i].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	  h_l2jjmass->Fill( (muons[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	  h_lljjmass->Fill( (muons[i].lorentzVec()+muons[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	  h_WandNmass->Fill( (muons[i].lorentzVec()+muons[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M() , (muons[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
          h_3Dparm->Fill( (muons[i].lorentzVec()+muons[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(), muons[i].lorentzVec().Pt(), muons[j].lorentzVec().Pt(), weight); 
          h_3DparmJ->Fill( jets[m].lorentzVec().Pt(), muons[i].lorentzVec().Pt(), muons[j].lorentzVec().Pt(), weight);
	  h_leadingMuonPt->Fill( muons[i].lorentzVec().Pt(),weight);
	  h_secondMuonPt->Fill( muons[j].lorentzVec().Pt(),weight);
	  h_leadingJetPt->Fill( jets[m].lorentzVec().Pt(),weight);
	  h_secondJetPt->Fill( jets[n].lorentzVec().Pt(),weight);
	  h_leadingMuonIso->Fill( muons[i].relIso(),weight);
	  h_secondMuonIso->Fill( muons[j].relIso(),weight);
          h_paircharge->Fill(muons[i].charge(),weight);
          h_muonseta->Fill(muons[i].eta(),weight);
          h_muonseta->Fill(muons[j].eta(),weight);
          h_jetseta->Fill(jets[m].eta(),weight);
          h_jetseta->Fill(jets[n].eta(),weight);
          h_bTag->Fill(jets[m].btag_disc(),weight);
          h_bTag->Fill(jets[n].btag_disc(),weight);
          h_cosTheta1->Fill(cos(muons[i].lorentzVec().Theta()),weight);
          h_cosTheta2->Fill(cos(muons[j].lorentzVec().Theta()),weight);
          h_DeltaPhi->Fill(fabs(muons[i].lorentzVec().DeltaPhi(muons[j].lorentzVec())),weight);
          h_Njets->Fill(jets.size(), weight);
	}
   }
   nogoodW:
   ;
}

void SignalPlots::Fill(Double_t MET, Double_t METphi, std::vector<Lepton>& muons, Double_t weight, Bool_t ptok, Bool_t ssok) {
  UInt_t m,n;
  for(UInt_t i=0; i<muons.size()-1; i++)
    for(UInt_t j=i+1; j<muons.size(); j++) {
      if (muons[i].charge()*muons[j].charge()>0 || ssok)
        if (muons[i].lorentzVec().Pt()>=20 || ptok) {
          //if (dijetmass > 200 || (jets[m].lorentzVec()+jets[n].lorentzVec()).M()>120) goto nogoodW;
          //if (dijetmass<50 || dijetmass>110) goto nogoodW;
          TVector3 vMET = TVector3(MET * cos(METphi), MET * sin(METphi), 0);
          TVector3 vlead = TVector3(muons[i].lorentzVec().Pt() * cos(muons[i].lorentzVec().Phi()), muons[i].lorentzVec().Pt() * sin(muons[i].lorentzVec().Phi()), 0);
          TVector3 vsec = TVector3(muons[j].lorentzVec().Pt() * cos(muons[j].lorentzVec().Phi()), muons[j].lorentzVec().Pt() * sin(muons[j].lorentzVec().Phi()), 0);
          h_MET->Fill(MET, weight);
          h_dMETphilead->Fill(vMET.Angle(vlead), weight);
          h_cosdMETphilead->Fill(cos(vMET.Angle(vlead)), weight);
          h_dMETphisecond->Fill(vMET.Angle(vsec), weight);
          h_cosdMETphisecond->Fill(cos(vMET.Angle(vsec)), weight);
          //m=0; n=1;
          h_llmass->Fill( (muons[i].lorentzVec()+muons[j].lorentzVec()).M(),weight);
          h_leadingMuonPt->Fill( muons[i].lorentzVec().Pt(),weight);
          h_secondMuonPt->Fill( muons[j].lorentzVec().Pt(),weight);
          StdPlots::h_invpT->Fill( muons[j].lorentzVec().Pt(),muons[j].eta(),weight);
          if(fabs(muons[j].eta()) < 0.8) {
            StdPlots::h_CRAFT->Fill( muons[j].lorentzVec().Pt(),weight);
            StdPlots::h_CRAFT_tag->Fill( muons[i].lorentzVec().Pt(),weight);
          }
          h_inv_leadingMuonPt->Fill( 1./muons[i].lorentzVec().Pt(),weight);
          h_inv_secondMuonPt->Fill( 1./muons[j].lorentzVec().Pt(),weight);
          h_leadingMuonPhi->Fill( muons[i].lorentzVec().Phi(),weight);
          h_secondMuonPhi->Fill( muons[j].lorentzVec().Phi(),weight);
          h_leadingMuonEta->Fill( muons[i].eta(),weight);
          h_secondMuonEta->Fill( muons[j].eta(),weight);
          //h_invPt->Fill( 1/(muons[j].lorentzVec().Pt()),weight);
          h_leadingMuonIso->Fill( muons[i].relIso(),weight);
          h_secondMuonIso->Fill( muons[j].relIso(),weight);
          h_paircharge->Fill(muons[i].charge(),weight);
          h_muonseta->Fill(muons[i].eta(),weight);
          h_muonseta->Fill(muons[j].eta(),weight);
          h_cosTheta1->Fill(cos(muons[i].lorentzVec().Theta()),weight);
          h_cosTheta2->Fill(cos(muons[j].lorentzVec().Theta()),weight);
          h_DeltaPhi->Fill(fabs(muons[i].lorentzVec().DeltaPhi(muons[j].lorentzVec())),weight);
        }
   }
   nogoodW:
   ;
}

void SignalPlots::Fill(Double_t MET, std::vector<Lepton>& muons, std::vector<Lepton>& muonsloose, std::vector<Jet>& jets, Double_t weight, Bool_t ptok, Bool_t ssok) {
  dijetmass_tmp=dijetmass=9999.9;
  UInt_t m,n;
  for(UInt_t i=0; i<muons.size(); i++)
    for(UInt_t j=0; j<muonsloose.size(); j++) {
      if (muons[i].charge()*muonsloose[j].charge()>0 || ssok)
	if (muons[i].lorentzVec().Pt()>=20 || muonsloose[j].lorentzVec().Pt()>=20 || ptok) {
	  if (muons[i].lorentzVec().Pt()>=muonsloose[j].lorentzVec().Pt()) {
	    for(UInt_t emme=0; emme<jets.size()-1; emme++)
              for(UInt_t enne=1; enne<jets.size(); enne++) {
                //dijetmass_tmp = (jets[emme].lorentzVec()+jets[enne].lorentzVec()).M();
                dijetmass_tmp = (muons[i].lorentzVec()+muonsloose[j].lorentzVec()+jets[emme].lorentzVec()+jets[enne].lorentzVec()).M();
                if ( fabs(dijetmass_tmp-Mass_W) < fabs(dijetmass-Mass_W) ) {
                  dijetmass = dijetmass_tmp;
                  m = emme;
                  n = enne;
                }
              }
            //if (dijetmass > 200 || (jets[m].lorentzVec()+jets[n].lorentzVec()).M()>120) goto nogoodW;
            //if (dijetmass<50 || dijetmass>110) goto nogoodW; 
            //m=0; n=1;
	    h_MET->Fill(MET, weight);
	    h_jjmass->Fill( (jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	    h_llmass->Fill( (muons[i].lorentzVec()+muonsloose[j].lorentzVec()).M(),weight);
	    h_l1jjmass->Fill( (muons[i].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	    h_l2jjmass->Fill( (muonsloose[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	    h_lljjmass->Fill( (muons[i].lorentzVec()+muonsloose[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	    h_WandNmass->Fill( (muons[i].lorentzVec()+muonsloose[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M() , (muons[i].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
            h_3Dparm->Fill( (muons[i].lorentzVec()+muons[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(), muons[i].lorentzVec().Pt(), muonsloose[j].lorentzVec().Pt(), weight);
	    h_3DparmJ->Fill( jets[m].lorentzVec().Pt(), muons[i].lorentzVec().Pt(), muons[j].lorentzVec().Pt(), weight);
            h_leadingMuonPt->Fill( muons[i].lorentzVec().Pt(),weight);
	    h_secondMuonPt->Fill( muonsloose[j].lorentzVec().Pt(),weight);
	    h_leadingJetPt->Fill( jets[m].lorentzVec().Pt(),weight);
	    h_secondJetPt->Fill( jets[n].lorentzVec().Pt(),weight);
	    h_leadingMuonIso->Fill( muons[i].relIso(),weight);
	    h_secondMuonIso->Fill( muonsloose[j].relIso(),weight);
	    h_paircharge->Fill(muons[i].charge(),weight);
	    h_muonseta->Fill(muons[i].eta(),weight);
	    h_muonseta->Fill(muons[j].eta(),weight);
	    h_jetseta->Fill(jets[m].eta(),weight);
	    h_jetseta->Fill(jets[n].eta(),weight);
	    h_bTag->Fill(jets[m].btag_disc(),weight);
	    h_bTag->Fill(jets[n].btag_disc(),weight);
	    h_cosTheta1->Fill(cos(muons[i].lorentzVec().Theta()),weight);
	    h_cosTheta2->Fill(cos(muonsloose[j].lorentzVec().Theta()),weight);
            h_DeltaPhi->Fill(fabs(muons[i].lorentzVec().DeltaPhi(muonsloose[j].lorentzVec())),weight); 
            h_Njets->Fill(jets.size(), weight);
	  }
	  else {
            for(UInt_t emme=0; emme<jets.size()-1; emme++)
              for(UInt_t enne=1; enne<jets.size(); enne++) {
                //dijetmass_tmp = (jets[emme].lorentzVec()+jets[enne].lorentzVec()).M();
                dijetmass_tmp = (muons[i].lorentzVec()+muonsloose[j].lorentzVec()+jets[emme].lorentzVec()+jets[enne].lorentzVec()).M();
                if ( fabs(dijetmass_tmp-Mass_W) < fabs(dijetmass-Mass_W) ) {
                  dijetmass = dijetmass_tmp;
                  m = emme;
                  n = enne;
                }
              }
            //if (dijetmass > 200 || (jets[m].lorentzVec()+jets[n].lorentzVec()).M()>120) goto nogoodW;
            //if (dijetmass<50 || dijetmass>110) goto nogoodW;
            //m=0; n=1;	      
	    h_MET->Fill(MET, weight);
	    h_jjmass->Fill( (jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	    h_llmass->Fill( (muonsloose[j].lorentzVec()+muons[i].lorentzVec()).M(),weight);
	    h_l1jjmass->Fill( (muonsloose[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	    h_l2jjmass->Fill( (muons[i].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	    h_lljjmass->Fill( (muonsloose[j].lorentzVec()+muons[i].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
	    h_WandNmass->Fill( (muonsloose[j].lorentzVec()+muons[i].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M() , (muonsloose[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(),weight);
            h_3Dparm->Fill( (muons[i].lorentzVec()+muons[j].lorentzVec()+jets[m].lorentzVec()+jets[n].lorentzVec()).M(), muonsloose[j].lorentzVec().Pt(), muons[i].lorentzVec().Pt(), weight);
            h_3DparmJ->Fill( jets[m].lorentzVec().Pt(), muons[i].lorentzVec().Pt(), muons[j].lorentzVec().Pt(), weight);
	    h_leadingMuonPt->Fill( muonsloose[j].lorentzVec().Pt(),weight);
	    h_secondMuonPt->Fill( muons[i].lorentzVec().Pt(),weight);
	    h_leadingJetPt->Fill( jets[m].lorentzVec().Pt(),weight);
	    h_secondJetPt->Fill( jets[n].lorentzVec().Pt(),weight);
	    h_leadingMuonIso->Fill( muonsloose[j].relIso(),weight);
	    h_secondMuonIso->Fill( muons[i].relIso(),weight);
	    h_paircharge->Fill(muons[i].charge(),weight);
	    h_muonseta->Fill(muons[i].eta(),weight);
	    h_muonseta->Fill(muons[j].eta(),weight);
	    h_jetseta->Fill(jets[m].eta(),weight);
	    h_jetseta->Fill(jets[n].eta(),weight);
	    h_bTag->Fill(jets[m].btag_disc(),weight);
	    h_bTag->Fill(jets[n].btag_disc(),weight);
	    h_cosTheta1->Fill(cos(muons[i].lorentzVec().Theta()),weight);
	    h_cosTheta2->Fill(cos(muonsloose[j].lorentzVec().Theta()),weight); 
            h_DeltaPhi->Fill(fabs(muons[i].lorentzVec().DeltaPhi(muonsloose[j].lorentzVec())),weight);
            h_Njets->Fill(jets.size(), weight);
	  }
	}
    }
  nogoodW:
  ;
}

void SignalPlots::Fill(Double_t MET, Double_t METphi, std::vector<Lepton>& muons, std::vector<Lepton>& muonsloose, Double_t weight, Bool_t ptok, Bool_t ssok) {
  dijetmass_tmp=dijetmass=9999.9;
  for(UInt_t i=0; i<muons.size(); i++)
    for(UInt_t j=0; j<muonsloose.size(); j++) {
      if (muons[i].charge()*muonsloose[j].charge()>0 || ssok)
        if (muons[i].lorentzVec().Pt()>=20 || muonsloose[j].lorentzVec().Pt()>=20 || ptok) {
          if (muons[i].lorentzVec().Pt()>=muonsloose[j].lorentzVec().Pt()) {
            //if (dijetmass > 200 || (jets[m].lorentzVec()+jets[n].lorentzVec()).M()>120) goto nogoodW;
            //if (dijetmass<50 || dijetmass>110) goto nogoodW;
            //m=0; n=1;
            TVector3 vMET = TVector3(MET * cos(METphi), MET * sin(METphi), 0);
            TVector3 vlead = TVector3(muons[i].lorentzVec().Pt() * cos(muons[i].lorentzVec().Phi()), muons[i].lorentzVec().Pt() * sin(muons[i].lorentzVec().Phi()), 0);
            TVector3 vsec = TVector3(muonsloose[i].lorentzVec().Pt() * cos(muonsloose[i].lorentzVec().Phi()), muonsloose[i].lorentzVec().Pt() * sin(muonsloose[i].lorentzVec().Phi()), 0);
            h_MET->Fill(MET, weight);
            h_dMETphilead->Fill(vMET.Angle(vlead), weight);
            h_cosdMETphilead->Fill(cos(vMET.Angle(vlead)), weight);
            h_dMETphisecond->Fill(vMET.Angle(vsec), weight);
            h_cosdMETphisecond->Fill(cos(vMET.Angle(vsec)), weight);
            h_llmass->Fill( (muons[i].lorentzVec()+muonsloose[j].lorentzVec()).M(),weight);
            if((muons[i].lorentzVec()+muonsloose[j].lorentzVec()).M() < 20)
              cout << "POG=" << muons[i].ilepton() << " Loose=" << muonsloose[j].ilepton() << endl;
            h_leadingMuonPt->Fill( muons[i].lorentzVec().Pt(),weight);
            h_secondMuonPt->Fill( muonsloose[j].lorentzVec().Pt(),weight);
            StdPlots::h_invpT->Fill( muons[j].lorentzVec().Pt(),muons[j].eta(),weight);
            if(fabs(muons[j].eta()) < 0.8) {
              StdPlots::h_CRAFT->Fill( muons[j].lorentzVec().Pt(),weight);
              StdPlots::h_CRAFT_tag->Fill( muons[i].lorentzVec().Pt(),weight);
            }
            h_inv_leadingMuonPt->Fill( 1./muons[i].lorentzVec().Pt(),weight);
            h_inv_secondMuonPt->Fill( 1./muonsloose[j].lorentzVec().Pt(),weight);
            h_leadingMuonPhi->Fill( muons[i].lorentzVec().Phi(),weight);
            h_secondMuonPhi->Fill( muonsloose[j].lorentzVec().Phi(),weight);
            h_leadingMuonEta->Fill( muons[i].eta(),weight);
            h_secondMuonEta->Fill( muonsloose[j].eta(),weight);
            //h_invPt->Fill( 1/(muonsloose[j].lorentzVec().Pt()),weight);
            h_leadingMuonIso->Fill( muons[i].relIso(),weight);
            h_secondMuonIso->Fill( muonsloose[j].relIso(),weight);
            h_paircharge->Fill(muons[i].charge(),weight);
            h_muonseta->Fill(muons[i].eta(),weight);
            h_muonseta->Fill(muons[j].eta(),weight);
            h_cosTheta1->Fill(cos(muons[i].lorentzVec().Theta()),weight);
            h_cosTheta2->Fill(cos(muonsloose[j].lorentzVec().Theta()),weight);
            h_DeltaPhi->Fill(fabs(muons[i].lorentzVec().DeltaPhi(muonsloose[j].lorentzVec())),weight);
          }
          else {
            //if (dijetmass > 200 || (jets[m].lorentzVec()+jets[n].lorentzVec()).M()>120) goto nogoodW;
            //if (dijetmass<50 || dijetmass>110) goto nogoodW;
            //m=0; n=1;
            TVector3 vMET = TVector3(MET * cos(METphi), MET * sin(METphi), 0);
            TVector3 vsec = TVector3(muons[i].lorentzVec().Pt() * cos(muons[i].lorentzVec().Phi()), muons[i].lorentzVec().Pt() * sin(muons[i].lorentzVec().Phi()), 0);
            TVector3 vlead = TVector3(muonsloose[i].lorentzVec().Pt() * cos(muonsloose[i].lorentzVec().Phi()), muonsloose[i].lorentzVec().Pt() * sin(muonsloose[i].lorentzVec().Phi()), 0);
            h_MET->Fill(MET, weight);
            h_dMETphisecond->Fill(vMET.Angle(vsec), weight);
            h_cosdMETphisecond->Fill(cos(vMET.Angle(vsec)), weight);
            h_dMETphilead->Fill(vMET.Angle(vlead), weight);
            h_cosdMETphilead->Fill(cos(vMET.Angle(vlead)), weight);
            h_llmass->Fill( (muonsloose[j].lorentzVec()+muons[i].lorentzVec()).M(),weight);
            h_leadingMuonPt->Fill( muonsloose[j].lorentzVec().Pt(),weight);
            h_secondMuonPt->Fill( muons[i].lorentzVec().Pt(),weight);
            StdPlots::h_invpT->Fill( muons[i].lorentzVec().Pt(),muons[i].eta(),weight);
            if(fabs(muons[i].eta()) < 0.8) {
              StdPlots::h_CRAFT->Fill( muons[i].lorentzVec().Pt(),weight);
              StdPlots::h_CRAFT_tag->Fill( muons[j].lorentzVec().Pt(),weight);
            }
            h_inv_leadingMuonPt->Fill( 1./muonsloose[j].lorentzVec().Pt(),weight);
            h_inv_secondMuonPt->Fill( 1./muons[i].lorentzVec().Pt(),weight);
            h_leadingMuonPhi->Fill( muonsloose[j].lorentzVec().Phi(),weight);
            h_secondMuonPhi->Fill( muons[i].lorentzVec().Phi(),weight);
            h_leadingMuonEta->Fill( muonsloose[j].eta(),weight);
            h_secondMuonEta->Fill( muons[i].eta(),weight);
            //h_invPt->Fill( 1/(muons[i].lorentzVec().Pt()),weight);
            h_leadingMuonIso->Fill( muonsloose[j].relIso(),weight);
            h_secondMuonIso->Fill( muons[i].relIso(),weight);
            h_paircharge->Fill(muons[i].charge(),weight);
            h_muonseta->Fill(muons[i].eta(),weight);
            h_muonseta->Fill(muons[j].eta(),weight);
            h_cosTheta1->Fill(cos(muons[i].lorentzVec().Theta()),weight);
            h_cosTheta2->Fill(cos(muonsloose[j].lorentzVec().Theta()),weight);
            h_DeltaPhi->Fill(fabs(muons[i].lorentzVec().DeltaPhi(muonsloose[j].lorentzVec())),weight);
          }
        }
    }
  nogoodW:
  ;
}

void SignalPlots::Write() {
  h_jjmass->Write();
  h_llmass->Write();
  h_l1jjmass->Write();
  h_l2jjmass->Write();
  h_lljjmass->Write();
  h_WandNmass->Write();
  h_3Dparm->Write();
  h_3DparmJ->Write();
  h_leadingMuonPt->Write();
  h_inv_leadingMuonPt->Write();
  h_secondMuonPhi->Write();
  h_leadingMuonPhi->Write();
  h_secondMuonEta->Write();
  h_leadingMuonEta->Write();
  h_secondMuonPt->Write();
  h_inv_secondMuonPt->Write();
  h_leadingJetPt->Write();
  h_secondJetPt->Write();
  h_leadingMuonIso->Write();
  h_secondMuonIso->Write();
  h_MET->Write();
  h_paircharge->Write();
  h_muonseta->Write();
  h_jetseta->Write();
  h_bTag->Write();
  h_cosTheta1->Write();
  h_cosTheta2->Write();
  h_DeltaPhi->Write();
  h_Njets->Write();
  h_dMETphilead->Write();
  h_dMETphisecond->Write();
  h_cosdMETphilead->Write();
  h_cosdMETphisecond->Write();
  StdPlots::Write();
  //h_invPt->Write();
}

