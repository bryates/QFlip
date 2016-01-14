#!/bin/bash
dataset=$1

if [ $dataset == "QCD_mu15" ]
then
  files=`ls | grep QCD_mu15_v3_job`
  numfiles=`ls | grep QCD_mu15_v3_job | wc -l`
else
  files=`ls | grep ${dataset}_job`
  numfiles=`ls | grep ${dataset}_job | wc -l`
fi
#echo $files
num=`ls | grep $dataset.txt | wc -l`

if [ $num == 0 ]
then
  echo "${dataset} has not been run yet"
  exit 0
fi

find ${dataset}_job* ! -size +1k

if [ $num == $numfiles ]
then
  echo "No ${dataset} files missing"
  exit 0
fi
echo $num
dif=$((num-$numfiles))
echo $dif
if [ $num == 0 ]
then
  exit 0
fi
#for i in {0..$num}
for (( i=0; i<$num; i++ ))
do
  if [[ $files != *"job${i}_"* ]]
  then
    echo "${dataset}_job${i}"
  fi
done
