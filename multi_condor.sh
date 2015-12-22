#!/bin/bash
dataset=$1

#if [[ $dataset == *"DY_50"* ]]
#then
#  list=/uscms_data/d3/byates/CMSSW_5_3_14_patch2/src/muon/DY_50.txt
#else
  list=/uscms_data/d2/fgior8/LQntuple_18/CMSSW_5_3_14_patch2_LQ/src/code/DataSetList/$dataset.txt
#fi
a=0
little=little
txt=$dataset.txt
b=0
for line in `cat $list`
do
#  echo "ciao $line"
  echo "$line" >> $little$b$txt
  a=$(($a+1))
#  echo $a
    #request_memory = 200KB
    #request_disk = 1MB
    #Should_Transfer_Files = YES
    #WhenToTransferOutput = ON_EXIT
  if [ $a -eq 50 ]
  then
    echo "condor_submit condor $little$b$txt $1_FR_$b"

    cat > simulator<<EOF

    universe = vanilla
    Executable = rotto.sh
    Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
    request_memory = 400MB
    request_disk = 10MB
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Output = log/sleep_FR40_job\$(Cluster).stdout
    Error = log/sleep_FR40_job\$(Cluster).stderr
    Log = log/sleep_FR40_job\$(Cluster).log
    notify_user = brent.yates@email.ucr.edu
    notification= Never
    arguments = \$(Cluster) $little$b$txt $b $dataset
    Queue 1

EOF

    condor_submit simulator

    a=0
    #rm "$little$b$txt"
    chmod 664 $little$b$txt
    b=$(($b+1))
  fi
done
echo "condor_submit condor $little$b$txt $1_FR_$b"

    cat > simulator<<EOF

    universe = vanilla
    Executable = rotto.sh
    Requirements = OpSys == "LINUX"&& (Arch != "DUMMY" )
    request_memory = 400MB
    request_disk = 10MB
    Should_Transfer_Files = YES
    WhenToTransferOutput = ON_EXIT
    Output = log/sleep_FR40_job\$(Cluster).stdout
    Error = log/sleep_FR40_job\$(Cluster).stderr
    Log = log/sleep_FR40_job\$(Cluster).log
    notify_user = brent.yates@email.ucr.edu
    notification= Never
    arguments = \$(Cluster) $little$b$txt $b $dataset
    Queue 1

EOF

    condor_submit simulator

#rm "$little$b$txt"
#chmod 664 $little$b$txt
exit 0
