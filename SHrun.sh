#!/bin/sh


dataset=$1
version=$2

cat > rootiamo.C<<EOF

{
  gROOT->ProcessLine(".L Reweight_cc.so");
  gROOT->ProcessLine(".L Data_cc.so");
  gROOT->ProcessLine(".L OtherFunctions_cc.so");
  gROOT->ProcessLine(".L SelectionFunctions_cc.so");
  gROOT->ProcessLine(".L LeptonSelection_cc.so");
  gROOT->ProcessLine(".L ElectronSelection_cc.so");
  gROOT->ProcessLine(".L MuonSelection_cc.so");
  gROOT->ProcessLine(".L MuonSelectionProbe_cc.so");
  gROOT->ProcessLine(".L JetSelection_cc.so");
  gROOT->ProcessLine(".L BTagSFUtil_C.so");
  gROOT->ProcessLine(".L GenSelection_cc.so");
  gROOT->ProcessLine(".L StdPlots_cc.so");
  gROOT->ProcessLine(".L ElectronPlots_cc.so");
  gROOT->ProcessLine(".L MuonPlots_cc.so");
  gROOT->ProcessLine(".L JetPlots_cc.so");
  gROOT->ProcessLine(".L SignalPlots_cc.so");
  gROOT->ProcessLine(".L Analyzer_cc.so");
  gROOT->ProcessLine(".L Analyzer_Ele_cc.so");
  gROOT->ProcessLine(".L FakeRateCalculator_cc.so");
  gROOT->ProcessLine(".L FakeRateCalculator_Ele_cc.so");
  gROOT->ProcessLine(".L EfficiencyCalculator_cc.so");
  gROOT->ProcessLine(".L ChainMaker_C.so");

  TString extralabel = "";

  //////////////////////////////////////////////////////////
  //                          MC                          //
  //////////////////////////////////////////////////////////
  if ("$dataset" == "N_40") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_40.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana40",$version);
     Pippo.SetWeight(1516, 99390);
     //Pippo.SetWeight(1, 99390);
     std::cout << "Majorana40\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_50") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_50.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana50",$version);
     Pippo.SetWeight(1071.1, 49996);
     std::cout << "Majorana50\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_60") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_60.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana60",$version);
     Pippo.SetWeight(607.7, 99494);
     //Pippo.SetWeight(1, 99494);
     std::cout << "Majorana60\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_70") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_70.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana70",$version);
     Pippo.SetWeight(211.96, 49994);
     std::cout << "Majorana70\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_80") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_80.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana80",$version);
     Pippo.SetWeight(19.07, 98492);
     //Pippo.SetWeight(1, 98492);
     std::cout << "Majorana80\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_90") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_90.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana90",$version);
     Pippo.SetWeight(7.1047, 49996);
     std::cout << "Majorana90\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_100") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_100.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana100",$version);
     Pippo.SetWeight(3.5618, 49994);
     //Pippo.SetWeight(1, 49994);
     std::cout << "Majorana100\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_125") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_125.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana125",$version);
     Pippo.SetWeight(1.0767, 49995);
     std::cout << "Majorana125\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_150") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_150.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana150",$version);
     Pippo.SetWeight(0.4594, 49995);
     //Pippo.SetWeight(1, 49995);
     std::cout << "Majorana150\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_175") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_175.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana175",$version);
     Pippo.SetWeight(0.23266, 49995);
     std::cout << "Majorana175\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_200") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_200.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana200",$version);
     Pippo.SetWeight(0.13127, 49992);
     std::cout << "Majorana200\n";  Pippo.Loop();
  }
    if ("$dataset" == "N_250") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_250.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana250",$version);
     Pippo.SetWeight(0.050928, 49997);
     std::cout << "Majorana250\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_300") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_300.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana300",$version);
     Pippo.SetWeight(0.023214, 49996);
     std::cout << "Majorana300\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_350") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_350.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana350",$version);
     Pippo.SetWeight(0.011705, 49995);
     std::cout << "Majorana350\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_400") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_400.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana400",$version);
     Pippo.SetWeight(0.0063324, 49996);
     std::cout << "Majorana400\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_500") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_500.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana500",$version);
     Pippo.SetWeight(0.0021542, 49995);
     std::cout << "Majorana500\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_600") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_600.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana600",$version);
     Pippo.SetWeight(8.545E-04, 49994);
     std::cout << "Majorana600\n";  Pippo.Loop();
  }
  if ("$dataset" == "N_700") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/N_700.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Majorana700",$version);
     Pippo.SetWeight(3.83E-04, 49998);
     std::cout << "Majorana700\n";  Pippo.Loop();
  }


  if ("$dataset" == "ZZ") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ZZ_inclusive.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ZZ_inclusive",$version);
     Pippo.SetWeight(17.654, 9799908);
     std::cout << "ZZ\n";  Pippo.Loop();
  }
  if ("$dataset" == "WZ") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WZ_inclusive.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WZ_inclusive",$version);
     Pippo.SetWeight(33.21, 10000283);
     std::cout << "WZ\n";  Pippo.Loop();
  }
  if ("$dataset" == "WpWp") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WpWp.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WpWp",$version);
     Pippo.SetWeight(0.2482, 99985);
     std::cout << "WpWp\n";  Pippo.Loop();
  }
  if ("$dataset" == "WmWm") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WmWm.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WmWm",$version);
     Pippo.SetWeight(0.08888, 96392);
     std::cout << "WmWm\n";  Pippo.Loop();
  }


  if ("$dataset" == "WWW") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WWW.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WWW",$version);
     Pippo.SetWeight(0.08217, 220549);
     std::cout << "WWW\n";  Pippo.Loop();
  }
  if ("$dataset" == "WWZ") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WWZ.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WWZ",$version);
     Pippo.SetWeight(0.05795, 222234);
     std::cout << "WWZ\n";  Pippo.Loop();
  }
  if ("$dataset" == "WZZ") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/WZZ.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WZZ",$version);
     Pippo.SetWeight(0.01968, 219835);
     std::cout << "WZZ\n";  Pippo.Loop();
  }
  if ("$dataset" == "ZZZ") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ZZZ.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ZZZ",$version);
     Pippo.SetWeight(0.005527, 224904);
     std::cout << "ZZZ\n";  Pippo.Loop();
  }

  if ("$dataset" == "ttW") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ttW.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ttW",$version);
     Pippo.SetWeight(0.232, 196046);
     std::cout << "ttW\n";  Pippo.Loop();
  }
  if ("$dataset" == "ttWW") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ttWW.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ttWW",$version);
     Pippo.SetWeight(0.002037, 217820);
     std::cout << "ttWW\n";  Pippo.Loop();
  }
  if ("$dataset" == "ttZ") {
     TChain* chain = ChainMaker("/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/ttZ.txt");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ttZ",$version);
     Pippo.SetWeight(0.1743, 210160);
     std::cout << "ttZ\n";  Pippo.Loop();
  }

}


EOF

root -q -b rootiamo.C

