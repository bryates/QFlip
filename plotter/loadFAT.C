#include <vector>

#include "SinglePlot.h"
#include "Single2dPlot.h"

using namespace std;

void loadCFO(std::vector<TString>& filename, std::vector<TString>& legendname, std::vector<TString>& plotlabel, std::vector<int>& color, std::vector<int>& linecol, std::vector<std::string>& type, std::vector<SinglePlot>& hist1d, std::vector<Single2dPlot>& hist2d, std::vector<double>& weight, std::vector<bool>& legend) {

  const double luminosity = 1.;
  const double EWK_scaling = 1.0;
  const double Z_scaling = 2.00463;
  const double W_scaling = 1.;

  //const TString directory = "/Volumes/Documents/Analyses/Majorana_Neutrino/PlotMakerMCOnly/files300/";
  const TString directory = "/uscms_data/d3/byates/CMSSW_5_3_14_patch2/src/muon/";
  std::vector<TString> classe;
  std::vector<TString> cuts;
  enum logbool {      nolog,         log };
  enum normbool {     nonorm,        norm };
  enum normfirstbool {nonormToFirst, normToFirst };
  enum stackbool {    nostack,       stack };
  enum overflowbool { nooverflow,    overflow };

  Bool_t ttbar=true; Bool_t DY_10=true; Bool_t DY_50=true;  Bool_t W_jets=true; Bool_t QCD=false;
  Bool_t WW = true; Bool_t WZ=true; Bool_t ZZ=true; Bool_t WpWp = true; Bool_t WmWm = true;
  //Bool_t WWW = true; Bool_t WWZ = true; Bool_t WZZ = true;
  Bool_t VVV = true;
  //Bool_t ttZ = true; Bool_t ttW = true; Bool_t ttWW = true;
  Bool_t ttV = true;
  Bool_t data=true;
  Bool_t plot1=true;
  Bool_t jet30=true;
  Bool_t fakes=false;
  
  if(plot1) {
    //classe.push_back("TL_denominator");
    classe.push_back("POG_muons");
    classe.push_back("POG_two_muons");
    //classe.push_back("TL_numerator");

    //cuts.push_back("Muons");
    cuts.push_back("TwoMuons");
    //cuts.push_back("MuJets");
    cuts.push_back("pt");
    //cuts.push_back("TagFlip");
    //cuts.push_back("PFRange");
    cuts.push_back("METRange");
    //cuts.push_back("HTRange");
    cuts.push_back("NoJets");
    cuts.push_back("NoJets_SS");
    cuts.push_back("NoJets_OS");
    //cuts.push_back("NJetsRange");

    for(UInt_t iii=0; iii<1; iii++) {
      cout<<classe[iii]<<endl;
      if (true) {
	hist1d.push_back( SinglePlot("Muons/h_N_"+classe[iii], 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_Detector_RelIso_rho", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_pt", 5,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "","muon p_{T} (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_eta", 2,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "","muon #eta","Events") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_phi", 2,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_charge", 1,  log, nonorm, nonormToFirst, 2.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_HCalIso", 1,  nolog, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_ECalIso", 1,  nolog, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_TrkIso", 1,  nolog, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_Detector_RelIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_HCalIsoDeposit", 1,  nolog, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_ECalIsoDeposit", 1,  nolog, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_GlbChi2", 1,  nolog, nonorm, nonormToFirst, 50.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_PF_RelIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_PF_RelIso_beta", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_photonIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_chargedHadronIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_neutralHadronIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_PUpt", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_dxy", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Muons/h_"+classe[iii]+"_dz", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("Muons/h_HT", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("Muons/h_MET", 2,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","#slash{E}_{T} (GeV)","Events/2 GeV") );
        //hist1d.push_back( SinglePlot("Muons/h_PFSumET", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("h_MET", 2,  nolog, nonorm, nonormToFirst, 200.0, nooverflow, stack, "","#slash{E}_{T} (GeV)","Events/2 GeV") );
    
	hist1d.push_back( SinglePlot("TwoMuons/h_N_"+classe[iii+1], 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_Detector_RelIso_rho", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_pt", 5,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "","muon p_{T} (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_eta", 2,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "","muon #eta","Events") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_phi", 2,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_charge", 1,  log, nonorm, nonormToFirst, 2.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_HCalIsoDeposit", 1,  nolog, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_ECalIsoDeposit", 1,  nolog, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_GlbChi2", 1,  nolog, nonorm, nonormToFirst, 50.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_PF_RelIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_PF_RelIso_beta", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_photonIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_chargedHadronIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_neutralHadronIso", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_PUpt", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_dxy", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_"+classe[iii+1]+"_dz", 1,  nolog, nonorm, nonormToFirst, .5, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("TwoMuons/h_twoMu_mass", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("TwoMuons/h_METsign", 2,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","#slash{E}_{T} (GeV)","Events/2 GeV") );
        //hist1d.push_back( SinglePlot("TwoMuons/h_twoMu_muonCharge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	//hist1d.push_back( SinglePlot("TwoMuons/h_twoMu_PFSumET", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","PF Sum {E}_{T} (GeV)","Events GeV") );

/*
	hist1d.push_back( SinglePlot("Jets/h_N_jets", 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Jets/h_jets_pt", 5,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "","muon p_{T} (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Jets/h_jets_eta", 2,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "","muon #eta","Events") );
	hist1d.push_back( SinglePlot("Jets/h_jets_phi", 2,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","") );
	
        hist1d.push_back( SinglePlot("Jets/h_N_jets_two_muons", 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
        hist1d.push_back( SinglePlot("Jets/h_jets_two_muons_pt", 5,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "","muon p_{T} (GeV)","Events/10 GeV") );
        hist1d.push_back( SinglePlot("Jets/h_jets_two_muons_eta", 2,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "","muon #eta","Events") );
        hist1d.push_back( SinglePlot("Jets/h_jets_two_muons_phi", 2,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","") );

        hist1d.push_back( SinglePlot("h_pt_mass", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
        hist1d.push_back( SinglePlot("h_pt_muonCharge", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","muon p_{T} (GeV)","Events/10 GeV") );
        hist1d.push_back( SinglePlot("h_pt_MET", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","muon #eta","Events") );
        hist1d.push_back( SinglePlot("h_pt_PFSumET", 1,  nolog, nonorm, nonormToFirst, 1.0 , nooverflow, stack, "","","") );
	
*/
	
      }
    }


    for(UInt_t i=0; i < cuts.size(); i++) {
        if (true) break;
/*
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_muonCharge", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Chage","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_muonCharge_barrel", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Chage","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_muonCharge_disk", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Chage","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_muonCharge_nw", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Chage (no weight)","Events") );
*/
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_mass", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Mass (GeV/c^{2})","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_mass", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Mass (GeV/c^{2})","Events") );
//        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_mass_nw", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Mass (GeV/c^{2}) (no weight)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_MET", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Missing E_{T} (GeV)","Events") );
/*
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_MET_nw", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Missing E_{T} (GeV) (no weight)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_PFSumET", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Sum E_{T} (GeV)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_PFSumET_nw", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Sum E_{T} (GeV) (no weight)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_PFSumETMinusMu", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Sum E_{T} - Muon E_{T} (GeV)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_PFSumETMinusMu_nw", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Sum E_{T} - Muon E_{T} (GeV) (no weight)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_HT", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Sum P_{T} of jets (GeV/c)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_HT", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Sum P_{T} of jets (GeV/c)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_nvtx", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Nvtx of "+cuts[i]+"","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_nvtx_nw", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Nvtx of "+cuts[i]+" (no weight)","Events") );
*/
        hist1d.push_back( SinglePlot(cuts[i]+"/h_N_"+cuts[i], 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Number of "+cuts[i]+"","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_pt", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" p_{T} (GeV/C)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_tag_pt", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" p_{T} (GeV/C)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_eta", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" #eta","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_tag_eta", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" #eta","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_phi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" #phi (phi)","Events") );
        hist1d.push_back( SinglePlot(cuts[i]+"/h_"+cuts[i]+"_tag_phi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" #phi (rad)","Events") );
/*
        hist1d.push_back( SinglePlot("DYFlipTag/h_DYFlipTag_deltaphi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","#phi_{MET} - #phi_{#mu} of DYFlipTag","Events") );
        hist1d.push_back( SinglePlot("DYFlipTag/h_DYFlipTag_cosdphi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","cos(#phi_{MET} - #phi_{#mu}) of DYFlipTag","Events") );
*/
        hist1d.push_back( SinglePlot("NoJets_Flip_SS/h_NoJets_Flip_SS_cosdphi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","cos(#phi_{MET} - #phi_{#mu}) of "+cuts[i]+"","Events") );
        hist1d.push_back( SinglePlot("NoJets_Flip_SS/h_NoJets_Flip_SS_deltaphi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","#phi_{MET} - #phi_{#mu} of "+cuts[i]+"","Events") );
        hist1d.push_back( SinglePlot("NoJets_Flip_OS/h_NoJets_Flip_OS_cosdphi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","cos(#phi_{MET} - #phi_{#mu}) of "+cuts[i]+"","Events") );
        hist1d.push_back( SinglePlot("NoJets_Flip_OS/h_NoJets_Flip_OS_deltaphi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","#phi_{MET} - #phi_{#mu} of "+cuts[i]+"","Events") );
/*
        hist1d.push_back( SinglePlot("debug/h_N_Rej_muons", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Number of muonRej","Events") );
        hist1d.push_back( SinglePlot("debug/h_Rej_muons_pt", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" p_{T} (GeV/C)","Events") );
        hist1d.push_back( SinglePlot("debug/h_Rej_muons_tag_pt", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" p_{T} (GeV/C)","Events") );
        hist1d.push_back( SinglePlot("debug/h_Rej_muons_eta", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" #eta","Events") );
        hist1d.push_back( SinglePlot("debug/h_Rej_muons_tag_eta", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" #eta","Events") );
        hist1d.push_back( SinglePlot("debug/h_Rej_muons_phi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" #phi (phi)","Events") );
        hist1d.push_back( SinglePlot("debug/h_Rej_muons_tag_phi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "",cuts[i]+" #phi (phi)","Events") );
*/
        //hist1d.push_back( SinglePlot("h_genFlip", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","Gen vs. Reco Q Flip of "+cuts[i]+"","Events") );
/*
	hist1d.push_back( SinglePlot("Jets/h_N_jets", 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Jets/h_jets_pt", 5,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "","muon p_{T} (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Jets/h_jets_eta", 2,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "","muon #eta","Events") );
	hist1d.push_back( SinglePlot("Jets/h_jets_phi", 2,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Jets/h_N_jets_two_muons", 1,  log, nonorm, nonormToFirst, 10.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("Jets/h_jets_two_muons_pt", 5,  log, nonorm, nonormToFirst, 300.0, nooverflow, stack, "","muon p_{T} (GeV)","Events/10 GeV") );
	hist1d.push_back( SinglePlot("Jets/h_jets_two_muons_eta", 2,  nolog, nonorm, nonormToFirst, 5.0, nooverflow, stack, "","muon #eta","Events") );
	hist1d.push_back( SinglePlot("Jets/h_jets_two_muons_phi", 2,  nolog, nonorm, nonormToFirst, 3.15, nooverflow, stack, "","","") );
*/
/*
	hist1d.push_back( SinglePlot("TotalFakes/h_lljjmass_tf", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("debug/h_lljjmass_signal", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
	hist1d.push_back( SinglePlot("debug/h_lljjmass_signal", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","") );
*/
/*
	hist1d.push_back( SinglePlot("h_TV", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","N before gen matching") );
	hist1d.push_back( SinglePlot("h_TVgen", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "","","N after gen matching") );
*/
    }
  }

 if (fakes) {
   filename.push_back(directory+"Fakes_FR.root");
   legendname.push_back("Fakes");
   plotlabel.push_back("Fakes");
   color.push_back(kCyan); linecol.push_back(kBlack);
   legend.push_back(true);
   type.push_back("mc");
   weight.push_back(1);
 }

 if (ZZ) {
    //if ( jet30 )
      filename.push_back(directory+"ZZ_inclusive_FR.root");
    //else
      //filename.push_back(directory+"ZZ_inclusive_21_5-22-15.root");
    legendname.push_back("ZZ");
    plotlabel.push_back("ZZ");
    color.push_back(kYellow+2); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double ZZ_scale_factor = 17.654/979998;
    //weight.push_back(ZZ_scale_factor * luminosity );
    weight.push_back(1);
  }
  
  if (WpWp) {
    filename.push_back(directory+"WpWp_FR.root");
    legendname.push_back("WpWp");
    plotlabel.push_back("WpWp");
    color.push_back(kGreen-1); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    weight.push_back(1);
  }

  if (WmWm) {
    filename.push_back(directory+"WmWm_FR.root");
    legendname.push_back("WmWm");
    plotlabel.push_back("WmWm");
    color.push_back(kGreen-2); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    weight.push_back(1);
  }

  if (VVV) {
    filename.push_back(directory+"VVV_FR.root");
    legendname.push_back("VVV");
    plotlabel.push_back("VVV");
    color.push_back(kRed); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    weight.push_back(1);
  }

  if (WW) {
    //if ( jet30 )
      filename.push_back(directory+"WW_inclusive_FR.root");
    //else
      //filename.push_back(directory+"WW_inclusive_21_5-22-15.root");
    legendname.push_back("WW");
    plotlabel.push_back("WW");
    color.push_back(kGreen+2); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double WW_scale_factor = 54.838/10000431;
    //weight.push_back(WW_scale_factor * luminosity );
    weight.push_back(1);
  }

 if (WZ) {
    //if ( jet30 )
      filename.push_back(directory+"WZ_inclusive_FR.root");
    //else
      //filename.push_back(directory+"WZ_inclusive_21_5-22-15.root");
    legendname.push_back("WZ");
    plotlabel.push_back("WZ");
    color.push_back(kGreen+3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double WZ_scale_factor = 33.21/10000293;
    //weight.push_back(WZ_scale_factor * luminosity );
    weight.push_back(1);
  }

  if (ttbar) {
    //if ( jet30 )
      filename.push_back(directory+"ttbar_FR.root");
    //else
      //filename.push_back(directory+"ttbar_21_5-22-15.root");
    legendname.push_back("t#bar{t}");
    plotlabel.push_back("ttbar");
    color.push_back(kBlue); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double ttbar_scale_factor = 248.9 / 234.;
    weight.push_back(1);
  }

  if (ttV) {
    filename.push_back(directory+"ttV_FR.root");
    legendname.push_back("t#bar{t}V");
    plotlabel.push_back("ttV");
    color.push_back(kMagenta); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double ttbar_scale_factor = 248.9 / 234.;
    weight.push_back(1);
  }

  if (DY_50) { // Signal only
    //if ( jet30 )
      filename.push_back(directory+"DY_50_FR.root");
    //else
      //filename.push_back(directory+"DY_50_21_5-22-15.root");
    legendname.push_back("DY m(ll)>50");
    plotlabel.push_back("DY");
    color.push_back(kYellow-3); linecol.push_back(kBlack);
    legend.push_back(false);
    type.push_back("singal_mc");
    const double DY50_scale_factor = 1177.3*3 / 3503.71 * Z_scaling;
    //weight.push_back(DY50_scale_factor);
    weight.push_back(1);
    //weight.push_back(Z_scaling);
  }

  if (DY_10) { // Signal only
    //if ( jet30 )
      filename.push_back(directory+"DY_10-50_FR.root");
    //else
      //filename.push_back(directory+"DY_10-50_21_5-22-15.root");
    legendname.push_back("DY 10<m(ll)<50");
    plotlabel.push_back("DY");
    color.push_back(kYellow); linecol.push_back(kBlack);
    legend.push_back(false);
    type.push_back("signal_mc");
    const double DY10_scale_factor = EWK_scaling*11050 / 37835275. * 1529/1413;
    //weight.push_back(DY10_scale_factor * luminosity );
    //weight.push_back(Z_scaling);
    weight.push_back(1);
  }

  if (DY_10) {
    //if ( jet30 )
      filename.push_back(directory+"DY_10-50_FR.root");
    //else
      //filename.push_back(directory+"DY_10-50_21_5-22-15.root");
    legendname.push_back("DY 10<m(ll)<50");
    plotlabel.push_back("DY");
    color.push_back(kYellow); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double DY10_scale_factor = EWK_scaling*11050 / 37835275. * 1529/1413;
    //weight.push_back(DY10_scale_factor * luminosity );
    //weight.push_back(Z_scaling);
    weight.push_back(1);
  }

  if (DY_50) {
    //if ( jet30 )
      filename.push_back(directory+"DY_50_FR.root");
    //else
      //filename.push_back(directory+"DY_50_21_5-22-15.root");
    legendname.push_back("DY m(ll)>50");
    plotlabel.push_back("DY");
    color.push_back(kYellow-3); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double DY50_scale_factor = 1177.3*3 / 3503.71 * Z_scaling;
    //weight.push_back(DY50_scale_factor);
    weight.push_back(1);
    //weight.push_back(Z_scaling);
  }

 if (W_jets) {
    //if ( jet30 )
      filename.push_back(directory+"Wjets_FR.root");
    //else
      //filename.push_back(directory+"Wjets_21_5-22-15.root");
    legendname.push_back("W_jets");
    plotlabel.push_back("W_jets");
    color.push_back(kGreen); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double W_jets_scale_factor = 12234.4*3 / 37509.0 * W_scaling;
    //weight.push_back(W_jets_scale_factor);
    weight.push_back(1);
    //weight.push_back(W_scaling);
  }

  if (QCD) {
    //if ( jet30 )
      filename.push_back(directory+"QCD_mu15_FR.root");
    //else
      //filename.push_back(directory+"QCD_mu15_v3_21_5-22-15.root");
    legendname.push_back("QCD");
    plotlabel.push_back("QCD");
    color.push_back(kOrange); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("mc");
    const double W_jets_scale_factor = 12234.4*3 / 37509.0;
    weight.push_back(1.);
    //weight.push_back(W_scaling);
  }

  if (data) {
    //if ( jet30 )
      //filename.push_back(directory+"SingleMu_10-6-15.root");
      filename.push_back(directory+"SingleMu_FR.root");
    //else
      //filename.push_back(directory+"SingleMu_21_5-22-15.root");
    legendname.push_back("Data");
    plotlabel.push_back("Data");
    color.push_back(kBlack); linecol.push_back(kBlack);
    legend.push_back(true);
    type.push_back("data");
    weight.push_back(1);
  }  
  

  /////////////////////
  //   Formatting    //
  /////////////////////

  // SinglePlot(std::string name, unsigned int rebin, bool log, bool normalize, bool normToFirst, double scaleXmax,
  //            bool overflowbin, bool stacked, TString title)

  // Single2dPlot(std::string name, std::string title, std::string drawOption, unsigned int rebinX, unsigned int rebinY)

} 

