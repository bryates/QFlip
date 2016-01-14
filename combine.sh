#!/bin/sh
#hadd -f SingleMu
if [ ! -z $1 ]
then
  version="_${1}"
else
  version="_FR"
fi
hadd -f DY_10-50$version.root DY_10-50_job*
hadd -f DY_50$version.root DY_50_job*
hadd -f ttbar$version.root ttbar_job*
hadd -f ZZ_inclusive$version.root ZZ_inclusive_job*
hadd -f WZ_inclusive$version.root WZ_inclusive_job*
hadd -f WW_inclusive$version.root WW_inclusive_job*
hadd -f WpWp$version.root WpWp_job*
hadd -f WmWm$version.root WmWm_job*
hadd -f ttZ$version.root ttZ_job*
hadd -f ttW$version.root ttW_job*
hadd -f ttWW$version.root ttWW_job*
hadd -f WWW$version.root WWW_job*
hadd -f ZZZ$version.root ZZZ_job*
hadd -f WWZ$version.root WWZ_job*
hadd -f WZZ$version.root WZZ_job*
#hadd -f W_jets$version.root W_jets_
#hadd -f QCD_mu15$version.root QCD_mu15_
./haddFR.sh $version
