#!/bin/bash
if [ ! -z $1 ]
then
  FR="${1}_"
else
  FR=""
fi
month=`date +"%m-" | sed -e "s/^0*//g"`
day=`date +"%d-%y" | sed -e "s/^0*//g"`
now=$month$day
echo $now
#hadd -f SingleMu_$FR$now.root SingleMu_FR40_*
hadd -f Wjets_$FR$now.root Wjets_FR40_*
hadd -f DY_10-50_$FR$now.root DY_10-50_FR40_*
hadd -f DY_50_$FR$now.root DY_50_FR40_*
hadd -f ttbar_$FR$now.root ttbar_FR40_*
hadd -f ZZ_inclusive_$FR$now.root ZZ_inclusive_FR40_*
hadd -f WZ_inclusive_$FR$now.root WZ_inclusive_FR40_*
hadd -f WW_inclusive_$FR$now.root WW_inclusive_FR40_*
hadd -f QCD_mu15_v3_$FR$now.root QCD_mu15_v3_FR40_*
