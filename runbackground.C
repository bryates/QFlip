{
  gROOT->ProcessLine(".L Reweight.cc+g");
  gROOT->ProcessLine(".L Data.cc+g");
  gROOT->ProcessLine(".L OtherFunctions.cc+g");
  gROOT->ProcessLine(".L SelectionFunctions.cc+g");
  gROOT->ProcessLine(".L LeptonSelection.cc+g");
  gROOT->ProcessLine(".L ElectronSelection.cc+g");
  gROOT->ProcessLine(".L MuonSelection.cc+g");
  //gROOT->ProcessLine(".L MuonSelectionProbe.cc+g");
  gROOT->ProcessLine(".L JetSelection.cc+g");
  gROOT->ProcessLine(".L BTagSFUtil.C+g");
  gROOT->ProcessLine(".L GenSelection.cc+g");
  gROOT->ProcessLine(".L StdPlots.cc+g");
  gROOT->ProcessLine(".L ElectronPlots.cc+g");
  gROOT->ProcessLine(".L MuonPlots.cc+g");
  gROOT->ProcessLine(".L JetPlots.cc+g");
  gROOT->ProcessLine(".L SignalPlots.cc+g");
  gROOT->ProcessLine(".L CutPlots.cc+g");
  gROOT->ProcessLine(".L Analyzer.cc+g");

  gROOT->ProcessLine(".L ChainMaker.C+g");

  //////////////////////////////////////////////////////////
  //                          MC                          //
  //////////////////////////////////////////////////////////

  if (0) {
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Wgamma_leptons.txt");
    FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Wgamma_FR40",31);
    Pippo.SetWeight(461.6, 4802358);
    std::cout << "Wgamma\n";  Pippo.LoopFR();
  }

  if (0) {
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/QCD_mu15.txt");
    //TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/QCD_mu15.txt");
    FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("QCD_mu15_FR",8); //Pippo.SetWeight(10.04, 2194752);
    std::cout << "QCD_mu\n";  Pippo.LoopFR();
  }

  if (0) {
    //TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/QCD_600-800_mu.txt");
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ggHZZ.txt");
    //FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("QCD_300-470_mu_FR_delete",8); //Pippo.SetWeight(10.04, 2194752);
    FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("ggHZZ",8); //Pippo.SetWeight(10.04, 2194752);
    std::cout << "QCD_mu\n";  Pippo.LoopFR();
  }

  if (0) {
    TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_0-15.txt");
    EffCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_0-15_Eff",9); //Pippo.SetWeight(10.04, 2194752);
    std::cout << "QCD_mu\n";  Pippo.LoopEff();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_11/CMSSW_5_3_8_LQ/src/code/DataSetList/ttbar.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("ttbar_FR",50); 
     Pippo.SetWeight(234., 6736135);
     std::cout << "ttbar\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Wjets.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Wjets_PLOT",8); Pippo.SetWeight(36257.2, 57709905*2316/1818);
     std::cout << "Wjets\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_0-15.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_0-15_FR40_30",24); 
     Pippo.SetWeight(4935.578, 200448);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_15-20.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_15-20_FR40_30",24);
     Pippo.SetWeight(172.6416, 200236);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_20-30.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_20-30_FR40_30",24); 
     Pippo.SetWeight(156.8129, 150120);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_30-50.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_30-50_FR40_30",24); 
     Pippo.SetWeight(102.5517, 150165);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_50-80.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_50-80_FR40_30",24); 
     Pippo.SetWeight(40.15132, 100200);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_80-120.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_80-120_FR40_30",24);
     Pippo.SetWeight(12.69172, 100172);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_120-170.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_120-170_FR40_30",24); 
     Pippo.SetWeight(3.56208, 100160);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_170-230.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_170-230_FR40_30",24);
     Pippo.SetWeight(0.9677884, 100172);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_230-300.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_230-300_FR40_30",24);
     Pippo.SetWeight(0.2684074, 70500);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/Zjet_300.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("Zjet_300_FR40_30",24); 
     Pippo.SetWeight(0.9677884, 100172);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/WW_inclusive.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("WW_inclusive_FR",9); //Pippo.SetWeight(0.9677884, 100172);
     std::cout << "Zjet\n";  Pippo.LoopFR();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/QCD_mumu.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_mumu_40",130); 
     //Pippo.SetWeight(154.0, 3701947);
     std::cout << "QCD_mumu\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ttbar.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ttbar",300); 
     Pippo.SetWeight(234.0, 6923652);
     std::cout << "ttbar\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_09/CMSSW_5_3_4_LQ/src/code/DataSetList/QCD_1000_mu.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_1000_mu_del",8); //Pippo.SetWeight(154.0, 3701947);
     std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_40.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana40_syn",114);
     //Pippo.SetWeight(1516, 99390);
     //Pippo.SetWeight(1, 99390);
     std::cout << "Majorana40\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_50.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana50",114); 
     Pippo.SetWeight(1071.1, 49996);
     std::cout << "Majorana50\n";  Pippo.Loop();
  }
 
  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_60.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana60",114);
     Pippo.SetWeight(607.7, 99494);
     //Pippo.SetWeight(1, 99494);
     std::cout << "Majorana60\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_70.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana70",114); 
     Pippo.SetWeight(211.96, 49994);
     std::cout << "Majorana70\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_80.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana80_sync",115);
     //Pippo.SetWeight(19.07, 98492);
     //Pippo.SetWeight(1, 98492);
     std::cout << "Majorana80\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_90.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana90",115); 
     Pippo.SetWeight(7.1047, 49996);
     std::cout << "Majorana90\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_100.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana100",115); 
     Pippo.SetWeight(3.5618, 49994);
     //Pippo.SetWeight(1, 49994);
     std::cout << "Majorana100\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_125.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana125",115);
     Pippo.SetWeight(1.0767, 49995);
     std::cout << "Majorana125\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_150.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana150",115);
     Pippo.SetWeight(0.4594, 49995);
     //Pippo.SetWeight(1, 49995);
     std::cout << "Majorana150\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_175.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana175",115);
     Pippo.SetWeight(0.23266, 49995);
     std::cout << "Majorana175\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_200.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana200_sync",117);
     //Pippo.SetWeight(0.13127, 49992);
     std::cout << "Majorana200\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_250.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana250",117);
     Pippo.SetWeight(0.050928, 49997);
     std::cout << "Majorana250\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_300.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana300_syn",117);
     //Pippo.SetWeight(0.023214, 49996);
     std::cout << "Majorana300\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_350.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana350",117);
     Pippo.SetWeight(0.011705, 49995);
     std::cout << "Majorana350\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_400.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana400",117);
     Pippo.SetWeight(0.0063324, 49996);
     std::cout << "Majorana400\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_500.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana500",117); 
     Pippo.SetWeight(0.0021542, 49995);
     std::cout << "Majorana500\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_600.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana600_syn",117);
     //Pippo.SetWeight(8.545E-04, 49994);
     std::cout << "Majorana600\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_700.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana700_sync",117);
     //Pippo.SetWeight(3.83E-04, 49998);
     std::cout << "Majorana700\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ZZ_inclusive.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ZZ_inclusive_WZctrl",130); 
     Pippo.SetWeight(17.654, 9799908);
     std::cout << "ZZ\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WZ_inclusive.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WZ_inclusive_WZctrl",130); 
     Pippo.SetWeight(33.21, 10000283);
     std::cout << "WZ\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WpWp.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WpWp",117); 
     Pippo.SetWeight(0.2482, 99985);
     std::cout << "WpWp\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WmWm.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WmWm",117); 
     Pippo.SetWeight(0.08888, 96392);
     std::cout << "WmWm\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WW_dp.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WW_dp",114); 
     Pippo.SetWeight(0.5879, 834040);
     std::cout << "WW_dp\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WWW.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WWW",117);
     Pippo.SetWeight(0.08217, 220549);
     std::cout << "WWW\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WWZ.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WWZ",117);
     Pippo.SetWeight(0.05795, 222234);
     std::cout << "WWZ\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WZZ.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WZZ",117);
     Pippo.SetWeight(0.01968, 219835);
     std::cout << "WZZ\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ZZZ.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ZZZ",117);
     Pippo.SetWeight(0.005527, 224904);
     std::cout << "ZZZ\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WWG.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WWG",101);
     Pippo.SetWeight(0.528, 304285);
     std::cout << "WWG\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ttW.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ttW",117); 
     Pippo.SetWeight(0.232, 196046);
     std::cout << "ttW\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ttWW.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ttWW",117);
     Pippo.SetWeight(0.002037, 217820);
     std::cout << "ttWW\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ttZ.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ttZ",117); 
     Pippo.SetWeight(0.1743, 210160);
     std::cout << "ttZ\n";  Pippo.Loop();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/DY_10-50.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("WR700_MNu350",46); 
     Pippo.SetWeight(8.598, 500195);
     std::cout << "WR\n";  Pippo.LoopFR();
  }

  //////////////////////////////////////////////////////////
  //                        Data                          //
  //////////////////////////////////////////////////////////


  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/DoubleMu.txt");
     FRCalculator Pippo; Pippo.Init(chain); Pippo.SetName("SingleMuA_del",1); //Pippo.SetWeight(.02);
     std::cout << "mu\n";  Pippo.LoopFR();
     //Pippo.~Analyzer();
  }

  if (0) {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src//code/DataSetList/SingleMuA.txt");
     // 1076311113.674 (mb)^-1 
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("DoubleMu_prova",1); //Pippo.SetWeight(.02);
     std::cout << "mu\n";  Pippo.Loop();
  }

}

