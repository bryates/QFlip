#!/bin/bash

directory=/eos/uscms/store/user/etiras/MAJORANA/MAJ/QCD/ee_m$1/
destination=/uscms_data/d2/fgior8/LQntuple_11/CMSSW_5_3_8_LQ/src/code/DataSetList_moved
TEXT=Majorana_ee_$1.txt

for file in `ls -S $directory`
do
  a=`echo $file | awk -F_ '{print $3}'`
  copiedfiles[$a]=0
done  

for file in `ls -S $directory`
do
  b=`echo $file | awk -F_ '{print $3}'`
  echo $b
  if [ ${copiedfiles[$b]} -eq 1 ]
  then
    echo "Duplicate"
    continue
  fi
  echo "ok"
#  echo "$destination/$file" >> $destination/$TEXT
  echo "$directory/$file" >> $destination/$TEXT
#  dccp  $directory/$file $destination/$file
  copiedfiles[$b]=1
done                                  

exit 0
