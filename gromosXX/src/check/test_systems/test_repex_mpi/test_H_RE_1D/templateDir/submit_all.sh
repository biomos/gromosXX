#!/usr/bin/bash 

#Load Gromos Settings
. ../../../gromos_settings.sh

replicas=3
max_step=$((${maxCores}/${replicas}))
cores_per_rep=($(seq 1 ${stepCores} ${max_step}));

jobs_prefix="1DHRepEx_${replicas}"
gromosBIN=${1}

for i in "${cores_per_rep[@]}"; 
do 
    echo "$i using nodes: $((i*replicas))";
    tmp_prefix="${jobs_prefix}_N${i}" 
    rm -rf ${tmp_prefix}
    cp -r template ${tmp_prefix}

    cd ${tmp_prefix}
    ./job.sh ${tmp_prefix} $((${i}*${replicas})) ${gromosBIN}  || echo "Failed: ${x} " >> ../Failed.log
    cd ..

done


