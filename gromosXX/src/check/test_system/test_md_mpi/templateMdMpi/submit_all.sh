#!/usr/bin/bash 
cores_per_rep=(1 2 3 4 );
#cores_per_rep=(8 10);
jobs_prefix="md_mpi"
gromosBIN=${1}

for i in "${cores_per_rep[@]}"; 
do 
echo "$i using nodes: ${i}";
tmp_prefix="${jobs_prefix}_N${i}" 
rm -r ${tmp_prefix}
cp -r template ${tmp_prefix}
cd ${tmp_prefix}
./job.sh ${tmp_prefix} ${i} ${gromosBIN}
cd ..

done
