#!/bin/sh
#hadd -f SingleMu
if [ ! -z $1 ]
then
  FR="${1}"
else
  FR="FR"
fi
hadd -f DY_10-50_$FR.root DY_10-50_FR40_*
hadd -f DY_50_$FR.root DY_50_FR40_*
hadd -f ttbar_$FR.root ttbar_FR40_*
hadd -f ZZ_inclusive_$FR.root ZZ_inclusive_FR40_*
hadd -f WZ_inclusive_$FR.root WZ_inclusive_FR40_*
hadd -f WW_inclusive_$FR.root WW_inclusive_FR40_*
hadd -f WpWp_$FR.root WpWp_FR40_*
hadd -f WmWm_$FR.root WmWm_FR40_*
hadd -f ttZ_$FR.root ttZ_FR40_*
hadd -f ttW_$FR.root ttW_FR40_*
hadd -f ttWW_$FR.root ttWW_FR40_*
hadd -f WWW_$FR.root WWW_FR40_*
hadd -f ZZZ_$FR.root ZZZ_FR40_*
hadd -f WWZ_$FR.root WWZ_FR40_*
hadd -f WZZ_$FR.root WZZ_FR40_*
#hadd -f W_jets_$FR.root W_jets_
#hadd -f QCD_mu15_$FR.root QCD_mu15_
./haddFR.sh $FR
