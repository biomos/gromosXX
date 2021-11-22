#!/usr/bin/bash 

#Load Gromos Settings
. ../../../gromos_settings.sh

cores_per_rep=(1); #2 
replicas=2

jobs_prefix="2DTHRepEx_${replicas}"
gromosBIN=${1}

for i in "${cores_per_rep[@]}"; 
do 
    echo "$i using nodes: $((i*replicas))";
    tmp_prefix="${jobs_prefix}_N${i}" 
    rm -rf ${tmp_prefix}
    cp -r template ${tmp_prefix}

    cd ${tmp_prefix}
    ./job.sh ${tmp_prefix} $((i*2)) ${gromosBIN}  || echo "Failed: ${x} " >> ../Failed.log
    cd ..

done


