#include <vector>

#include "SinglePlot.h"
#include "Single2dPlot.h"

void loadCFO(std::vector<TString>& filename, std::vector<TString>& legendname, std::vector<TString>& plotlabel, std::vector<int>& color, std::vector<int>& linecol, std::vector<std::string>& type, std::vector<SinglePlot>& hist1d, std::vector<Single2dPlot>& hist2d, std::vector<double>& weight, std::vector<bool>& legend) {

//  const double luminosity = 36.09;
  const double luminosity = 80.582679;

  //const TString directory = "/Volumes/Documents/Analyses/SUSY_FATjets/files/1/";
  const TString directory = "/uscms_data/d3/byates/CMSSW_5_3_14_patch2/src/muon/";
  std::vector<TString> classe;
  enum logbool {      nolog,         log };
  enum normbool {     nonorm,        norm };
  enum normfirstbool {nonormToFirst, normToFirst };
  enum stackbool {    nostack,       stack };
  enum overflowbool { nooverflow,    overflow };

  Bool_t ttbar=true; Bool_t Z_jets=false; Bool_t W_jets=true; Bool_t QCD=false;
  Bool_t T1tttt=false;
  Bool_t data=true;
  Bool_t signal=true;
  Bool_t DY=true;

  
  if(signal) {
    //control
    /*hist1d.push_back( SinglePlot("h_leadingJetPt_control", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetPtLog_control", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetMass_control", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetMassLog_control", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    //smoothing
    hist1d.push_back( SinglePlot("h_leadingJetPt_smoothing", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetPtLog_smoothing", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetMass_smoothing", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetMassLog_smoothing", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    //dressed
    hist1d.push_back( SinglePlot("h_leadingJetPt_dressed", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetPtLog_dressed", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetMass_dressed", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetMassLog_dressed", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    //signal
    hist1d.push_back( SinglePlot("h_leadingJetPt_signal", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetPtLog_signal", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetMass_signal", 1,  log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );
    hist1d.push_back( SinglePlot("h_leadingJetMassLog_signal", 1,  nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "", "leading jet p_{T} (GeV)","Events/10 GeV") );

    hist2d.push_back( Single2dPlot("h_leadingJetPtvsMassLog_control", "h_leadingJetPtvsMassLog_control", "COLZ", 1, 1) );
    hist2d.push_back( Single2dPlot("h_leadingJetPtvsMassLog_smoothing", "h_leadingJetPtvsMassLog_smoothing", "COLZ", 1, 1) );*/
    /*hist1d.push_back( SinglePlot("h_muonCharge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Muon Charge", "", "") );
    hist1d.push_back( SinglePlot("h_MET", 1, nolog, nonorm,normToFirst, 1.0, nooverflow,stack, "Missing E_{T}", "E_{T} (GeV)", "Events") );
    hist1d.push_back( SinglePlot("h_PFSumET", 1, nolog, nonorm,nonormToFirst, 1.0, nooverflow,stack, "Sum E_{T}", "E_{T} (GeV)", "Events") );
    hist1d.push_back( SinglePlot("h_PFSumETMinusMu", 1, nolog, nonorm,nonormToFirst, 1.0, nooverflow,stack, "Sum E_{T} - Muon E_{T}", "E_{T} (GeV)", "Events") );
    hist1d.push_back( SinglePlot("h_HT", 1, nolog, nonorm,nonormToFirst, 1.0, nooverflow,stack, "Sum P_{T} of jets", "P_{T} (GeV/c)", "Events") );*/
    hist2d.push_back( Single2dPlot("h_nEvents", "", "COLZ", 1, 1) );

//Charge
    hist1d.push_back( SinglePlot("h_muons_muonCharge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Muon Charge", "", "") );
    hist1d.push_back( SinglePlot("h_twoMu_muonCharge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Muon Charge", "", "") );
    hist1d.push_back( SinglePlot("h_singleIso_muonCharge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Muon Charge", "", "") );
    hist1d.push_back( SinglePlot("h_pt_muonCharge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Muon Charge", "", "") );
    hist1d.push_back( SinglePlot("h_PFRange_muonCharge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Muon Charge", "", "") );

//MET    
    hist1d.push_back( SinglePlot("h_muons_MET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "MET", "", "") );
    hist1d.push_back( SinglePlot("h_twoMu_MET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "MET", "", "") );
    hist1d.push_back( SinglePlot("h_singleIso_MET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "MET", "", "") );
    hist1d.push_back( SinglePlot("h_pt_MET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "MET", "", "") );
    hist1d.push_back( SinglePlot("h_PFRange_MET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "MET", "", "") );

//PFSumET
    hist1d.push_back( SinglePlot("h_muons_PFSumET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET", "", "") );
    hist1d.push_back( SinglePlot("h_twoMu_PFSumET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET", "", "") );
    hist1d.push_back( SinglePlot("h_singleIso_PFSumET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET", "", "") );
    hist1d.push_back( SinglePlot("h_pt_PFSumET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET", "", "") );
    hist1d.push_back( SinglePlot("h_PFRange_PFSumET", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET", "", "") );

//PFSumETMinusMu
    hist1d.push_back( SinglePlot("h_muons_PFSumETMinusMu", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET - Muon E_{T}", "", "") );
    hist1d.push_back( SinglePlot("h_twoMu_PFSumETMinusMu", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET - Muon E_{T}", "", "") );
    hist1d.push_back( SinglePlot("h_singleIso_PFSumETMinusMu", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET - Muon E_{T}", "", "") );
    hist1d.push_back( SinglePlot("h_pt_PFSumETMinusMu", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "PFSumET - Muon E_{T}", "", "") );
    hist1d.push_back( SinglePlot("h_PFRange_PFSumETMinusMu", 1, nolog, nonorm, nonormToFirst, 500.0, nooverflow, stack, "PFSumET - Muon E_{T}", "", "") );

//Muons
    hist1d.push_back( SinglePlot("Muons/h_N_POG_muons", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Number of POG_muons", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_pt", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons p_{t} (GeV)", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_eta", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons #eta", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_phi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons #phi", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_charge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Charge of POG_muons", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_HCalIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons HCal Iso", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_ECalIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons ECal Iso", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_TrkIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Tracker Iso", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_Detector_RelIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Detector_RelIso", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_Detector_RelIso_rho", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Detector_RelIso #rho corrected", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_HCalIsoDeposit", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons HCal Iso Deposit", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_ECalIsoDeposit", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons ECal Iso Deposit", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_photonIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons photon Iso", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_chargedHadronIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons charged Hadron Iso", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_neutralHadronIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons neutral Hadron Iso", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_PF_RelIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons neutral Hadron Iso", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_PF_RelIso_beta", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Particle Flow RelIso #beta corrected", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_PUpt", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Particle Flow pileup correction", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_GlbChi2", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Global #chi^{2} per #DoF", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_dxy", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons transverse IP", "", "") );
    hist1d.push_back( SinglePlot("Muons/h_POG_muons_dz", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons longitudinal IP", "", "") );

//TwoMuons
    hist1d.push_back( SinglePlot("h_twoMu_mass", 1, nolog, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Invariant xass of two_muons", "Mass (GeV/c^{2})", "") );
    hist1d.push_back( SinglePlot("h_twoMu_mass", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Invariant xass of two_muons", "Mass (GeV/c^{2})", "") );
    hist1d.push_back( SinglePlot("h_twoMu_MET", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "MET of two_muons", "MET (GeV)", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_N_POG_two_muons", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Number of POG_two_muons", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_pt", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons p_{t} (GeV)", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_eta", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons #eta", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_phi", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons #phi", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_charge", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "Charge of POG_muons", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_HCalIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons HCal Iso", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_ECalIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons ECal Iso", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_TrkIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Tracker Iso", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_Detector_RelIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Detector_RelIso", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_Detector_RelIso_rho", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Detector_RelIso #rho corrected", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_HCalIsoDeposit", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons HCal Iso Deposit", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_ECalIsoDeposit", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons ECal Iso Deposit", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_photonIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons photon Iso", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_chargedHadronIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons charged Hadron Iso", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_neutralHadronIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons neutral Hadron Iso", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_PF_RelIso", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons neutral Hadron Iso", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_PF_RelIso_beta", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Particle Flow RelIso #beta corrected", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_PUpt", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Particle Flow pileup correction", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_GlbChi2", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons Global #chi^{2} per #DoF", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_dxy", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons transverse IP", "", "") );
    hist1d.push_back( SinglePlot("TwoMuons/h_POG_two_muons_dz", 1, log, nonorm, nonormToFirst, 1.0, nooverflow, stack, "POG_muons longitudinal IP", "", "") );

  }

  if (signal && ttbar) {
    filename.push_back(directory+"ttbar_1-21-15.root");
    legendname.push_back("t#bar{t}");
    plotlabel.push_back("t#bar{t}");
    color.push_back(kRed); linecol.push_back(kRed);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.0816);
    weight.push_back(1);
  }

  if (signal && W_jets) {
    filename.push_back(directory+"Wjets_1-21-15.root");
    legendname.push_back("W");
    plotlabel.push_back("W");
    color.push_back(kGreen); linecol.push_back(kGreen);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.1043);
    weight.push_back(1);
  }
  
  if (signal && Z_jets) {
    filename.push_back(directory+"ZJets_1.root");
    legendname.push_back("Z");
    plotlabel.push_back("Z");
    color.push_back(kYellow); linecol.push_back(kYellow);
    legend.push_back(true);
    type.push_back("mc");
    weight.push_back(0.1237);
  }
  
  if (signal && QCD) {
    filename.push_back(directory+"QCD_1000_1.root");
    legendname.push_back("QCD 1000");
    plotlabel.push_back("QCD 1000");
    color.push_back(kOrange); linecol.push_back(kOrange);
    legend.push_back(true);
    type.push_back("mc");
    weight.push_back(0.3965);
  }
  
  if (signal && QCD) {
    filename.push_back(directory+"QCD_500HT1000_1.root");
    legendname.push_back("QCD 500-1000");
    plotlabel.push_back("QCD 500-1000");
    color.push_back(kOrange+2); linecol.push_back(kOrange+2);
    legend.push_back(true);
    type.push_back("mc");
    weight.push_back(7.1869);
  }

  if (signal && T1tttt) {
    filename.push_back(directory+"SMS_T1tttt_mGo1400_mLSP25_1.root");
    legendname.push_back("T1tttt");
    plotlabel.push_back("T1tttt");
    color.push_back(kMagenta); linecol.push_back(kMagenta);
    legend.push_back(true);
    type.push_back("signal");
    weight.push_back(0.00013);
  }

  if (signal && DY) {
    filename.push_back(directory+"DY_50_1-21-15.root");
    legendname.push_back("DY 50");
    plotlabel.push_back("DY 50");
    color.push_back(kMagenta); linecol.push_back(kMagenta);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.4411);
    weight.push_back(1);
  }
  
  if (signal && DY) {
    filename.push_back(directory+"DY_10-50_1-21-15.root");
    legendname.push_back("DY 10-50");
    plotlabel.push_back("DY 10-50");
    color.push_back(kMagenta+2); linecol.push_back(kMagenta+2);
    legend.push_back(true);
    type.push_back("mc");
    //weight.push_back(0.4411);
    weight.push_back(1);
  }

  if (signal && data) {
    filename.push_back(directory+"SingleMu_1-21-15.root");
    legendname.push_back("Data");
    plotlabel.push_back("Data");
    color.push_back(kBlue); linecol.push_back(kBlue);
    legend.push_back(true);
    type.push_back("data");
    weight.push_back(1);
  }

  /////////////////////
  //   Formatting    //
  /////////////////////

  // SinglePlot(std::string name, unsigned int rebin, bool nolog, bool normalize, bool normToFirst, double scaleXmax,
  //            bool overflowbin, bool stacked, TString title)

  // Single2dPlot(std::string name, std::string title, std::string drawOption, unsigned int rebinX, unsigned int rebinY)

} 

