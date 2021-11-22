#!/usr/bin/bash 

# load settings
. ../../gromos_settings.sh

#how many cares to be tested on 
cores_per_rep=(1 2 3 4 );

jobs_prefix="md_mpi"
gromosBIN=${1}

for i in "${cores_per_rep[@]}"; 
do 
    tmp_prefix="${jobs_prefix}_N${i}" 
    echo -e "\t Testing: $tmp_prefix using nodes: ${i}    ${gromosBIN}";


    rm -rf ${tmp_prefix}
    cp -r template ${tmp_prefix}

    cd ${tmp_prefix}
    ./job.sh ${tmp_prefix} ${i} ${gromosBIN} || echo "Failed: ${x} " >> ../Failed.log
    cd ..

done
