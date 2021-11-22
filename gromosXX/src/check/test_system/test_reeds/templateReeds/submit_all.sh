#!/usr/bin/bash 

cores_per_rep=(2); #1
replicas=3
jobs_prefix="1DrepEx_${replicas}_newEnv_mpi2"
gromosBIN=${1}

for i in "${cores_per_rep[@]}"; 
do 
echo "$i using nodes: $((i*replicas))";
tmp_prefix="${jobs_prefix}_N${i}" 
rm -r ${tmp_prefix}
cp -r template ${tmp_prefix}
cd ${tmp_prefix}
./job.sh ${tmp_prefix} $((i*2)) ${gromosBIN}
cd ..

done

#-R "select[model==XeonGold_5118]"
