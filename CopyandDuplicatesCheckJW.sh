#!/bin/bash

directory=/pnfs/cms/WAX/11/store/user/jwwetzel/MAJ/QCD/$1
destination=/uscms_data/d2/fgior8/LQntuple_11/CMSSW_5_3_8_LQ/src/code/DataSetList
TEXT=QCD_X_mu.txt

for file in `ls -S $directory`
do
  a=`echo $file | awk -F_ '{print $5}'`
  copiedfiles[$a]=0
done  

for file in `ls -S $directory`
do
  b=`echo $file | awk -F_ '{print $5}'`
  echo $b
  if [ ${copiedfiles[$b]} -eq 1 ]
  then
    echo "Duplicate"
    continue
  fi
  echo "ok"
#  echo "$destination/$file" >> $destination/$TEXT
  echo "dcache:$directory/$file" >> $destination/$TEXT
#  dccp  $directory/$file $destination/$file
  copiedfiles[$b]=1
done                                  

exit 0
