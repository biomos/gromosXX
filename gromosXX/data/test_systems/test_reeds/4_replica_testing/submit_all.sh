#!/usr/bin/bash 

# User specific aliases and functions
. ~/new_env.sh


cores_per_rep=(1 2 3 4 6 8 10);
cores_per_rep=(1 2 3 4 );
replicas=4
jobs_prefix="1DrepEx_${replicas}_newEnv_mpi2"

for i in "${cores_per_rep[@]}"; 
do 
echo "$i using nodes: $((i*replicas))";
tmp_prefix="${jobs_prefix}_N${i}" 
rm -r ${tmp_prefix}
cp -r template ${tmp_prefix}
cd ${tmp_prefix}
bsub -J${tmp_prefix} -n$((i*replicas)) -W 04:00   "./job.sh ${tmp_prefix}"
cd ..

done

#-R "select[model==XeonGold_5118]"
