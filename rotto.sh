#!/bin/sh

dir=/uscms_data/d3/byates/CMSSW_5_3_14_patch2/src/muon
#dir=/uscms_data/d2/fgior8/LumiCal/CMSSW_5_3_9/src/code
echo $dir
cd $dir

#source /uscmst1/prod/sw/cms/bashrc prod
echo `pwd`
eval `scramv1 runtime -sh`

clusterNumber=$1
list=$2
num=$3
dataset=$4
echo $dataset
cd ${_CONDOR_SCRATCH_DIR}
echo  `pwd`

#cp /uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/*so .
ls $dir/*.so
cp $dir/*.so .

echo "Seeds:"
echo $var0 $var1

export CMSSW_SEARCH_PATH=${CMSSW_SEARCH_PATH}:`pwd`


cat > rootiamo_job${clusterNumber}.C<<EOF

{
  gROOT->ProcessLine(".L Reweight_cc.so");
  gROOT->ProcessLine(".L Data_cc.so");
  gROOT->ProcessLine(".L OtherFunctions_cc.so");
  gROOT->ProcessLine(".L SelectionFunctions_cc.so");
  gROOT->ProcessLine(".L LeptonSelection_cc.so");
  gROOT->ProcessLine(".L ElectronSelection_cc.so");
  gROOT->ProcessLine(".L MuonSelection_cc.so");
  //gROOT->ProcessLine(".L MuonSelectionProbe.cc+g");
  gROOT->ProcessLine(".L JetSelection_cc.so");
  gROOT->ProcessLine(".L BTagSFUtil_C.so");
  gROOT->ProcessLine(".L GenSelection_cc.so");
  gROOT->ProcessLine(".L StdPlots_cc.so");
  gROOT->ProcessLine(".L ElectronPlots_cc.so");
  gROOT->ProcessLine(".L MuonPlots_cc.so");
  gROOT->ProcessLine(".L JetPlots_cc.so");
  gROOT->ProcessLine(".L SignalPlots_cc.so");
  gROOT->ProcessLine(".L CutPlots_cc.so");
  gROOT->ProcessLine(".L Analyzer_cc.so");
  gROOT->ProcessLine(".L ChainMaker_C.so");

  TString extralabel = "";

  //////////////////////////////////////////////////////////
  //                          MC                          //
  //////////////////////////////////////////////////////////

  if ("$dataset" == "QCD_15-20_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_15-20_mu_job$num",300); 
    Pippo.SetWeight(702200000*0.0039, 1722681);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_20-30_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_20-30_mu_job$num",300); 
    Pippo.SetWeight(287000000*0.0065, 8486904);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_30-50_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_30-50_mu_job$num",300); 
    Pippo.SetWeight(66090000*0.0122, 9560265);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_50-80_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_50-80_mu_job$num",300); 
    Pippo.SetWeight(8082000*0.0218, 10365230);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_80-120_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_80-120_mu_job$num",300); 
    Pippo.SetWeight(1024000*0.0395, 9238642);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_120-170_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_120-170_mu_job$num",300); 
    Pippo.SetWeight(157800*0.0473, 8501935);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_170-300_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_170-300_mu_job$num",300);
    Pippo.SetWeight(34020*0.0676, 7669947);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_300-470_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_300-470_mu_job$num",300);
    Pippo.SetWeight(1757*0.0864, 7832261);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_470-600_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_470-600_mu_job$num",300);
    Pippo.SetWeight(115.2*0.1024, 3783069);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_600-800_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_600-800_mu_job$num",300);
    Pippo.SetWeight(27.01*0.0996, 4119000);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_800-1000_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_800-1000_mu_job$num",300);
    Pippo.SetWeight(3.57*0.1033, 4107853);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_1000_mu") {
    TChain* chain = ChainMaker("$dir/$list");
    Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_1000_mu_job$num",300);
    Pippo.SetWeight(0.774*0.1097, 3873970);
    std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "ttbar") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ttbar_job$num",300); 
     Pippo.SetWeight(234., 6923652);
     std::cout << "ttbar\n";  Pippo.Loop();
  }

  if ("$dataset" == "W_jets") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("Wjets_job$num",300);
     Pippo.SetWeight(37509.0, 57709905);
     std::cout << "Wjets\n";  Pippo.Loop();
  }

  if ("$dataset" == "DY_10-50") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("DY_10-50_job$num",300); 
     Pippo.SetWeight(11050., 37835275);
     std::cout << "DY\n";  Pippo.Loop();
  }

  if ("$dataset" == "DY_50") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("DY_50_job$num",300); 
     Pippo.SetWeight(3503.71, 30459503);
     std::cout << "DY\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_mu15") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_mu15_v3_job$num",300); 
     Pippo.SetWeight(364000000*0.00037, 21484602);
     std::cout << "QCD_mu\n";  Pippo.Loop();
  }

  if ("$dataset" == "QCD_mumu") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("QCD_mumu_job$num",31); 
     Pippo.SetWeight(4.9590002E10*4.0E-7, 1074400);
     std::cout << "QCD_mumu\n";  Pippo.Loop();
  }

  if ("$dataset" == "ZZ_inclusive") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("ZZ_inclusive_job$num",300); 
     Pippo.SetWeight(17.654, 9799908);
     std::cout << "ZZ\n";  Pippo.Loop();
  }

  if ("$dataset" == "WZ_inclusive") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WZ_inclusive_job$num",300); 
     Pippo.SetWeight(33.21, 10000283);
     std::cout << "WZ\n";  Pippo.Loop();
  }

  if ("$dataset" == "WW_inclusive") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("WW_inclusive_job$num",300);
     Pippo.SetWeight(54.838, 10000431);
     std::cout << "WW\n";  Pippo.Loop();
  }

  //////////////////////////////////////////////////////////
  //                        Data                          //
  //////////////////////////////////////////////////////////

  if ("$dataset" == "SingleMu") {
     TChain* chain = ChainMaker("$dir/$list");
     Analyzer Pippo; Pippo.Init(chain); Pippo.SetName("SingleMu_job$num",300);
     std::cout << "mu\n";  Pippo.Loop();
  }
}


EOF

root -q -b rootiamo_job${clusterNumber}.C
