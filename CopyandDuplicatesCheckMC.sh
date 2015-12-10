#!/bin/bash

directory=/eos/uscms/store/user/fgior8/Summer12_tag18/$1
destination=/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList
TEXT=$1.txt

for file in `ls -S $directory`
do
  a=`echo $file | awk -F_ '{print $2}'`
  copiedfiles[$a]=0
done  

for file in `ls -S $directory`
do
  b=`echo $file | awk -F_ '{print $2}'`
  echo $b
  if [ ${copiedfiles[$b]} -eq 1 ]
  then
    echo "Duplicate"
    rm $directory/$file
    continue
  fi
  echo "ok"
  echo "$directory/$file" >> $destination/$TEXT
#  echo "dcache:$directory/$file" >> $destination/$TEXT
#  dccp  $directory/$file $destination/$file
  copiedfiles[$b]=1
done                                  

for c in  {1..2}
do
  if [ -z "${copiedfiles[$c]}" ]
  then
    echo "Not copied"
    echo $c
  fi
done

exit 0
