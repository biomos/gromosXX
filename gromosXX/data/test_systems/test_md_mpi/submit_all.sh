#!/usr/bin/bash 
cores_per_rep=(1 2 3 4 6 8 10);
#cores_per_rep=(8 10);
jobs_prefix="md_mpi"

for i in "${cores_per_rep[@]}"; 
do 
echo "$i using nodes: ${i}";
tmp_prefix="${jobs_prefix}_N${i}" 
rm -r ${tmp_prefix}
cp -r template ${tmp_prefix}
cd ${tmp_prefix}
bsub -J${tmp_prefix} -n$((i*replicas)) -W 04:00 -R "select[model==XeonGold_5118]"  "./job.sh ${tmp_prefix}"
cd ..

done
