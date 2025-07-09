#!/bin/bash

SCRIPT=$PWD/gromos_dockerfile_generator.py

gitlab_project_tag=gromosx.boku.ac.at:5005/gromos/gromosxx/

#args[${#args[@]}]="--ubuntu 24.04"
#args[${#args[@]}]="--ubuntu 22.04"
#args[${#args[@]}]="--ubuntu 20.04"
#args[${#args[@]}]="--ubuntu 24.04 --llvm"
#args[${#args[@]}]="--ubuntu 22.04 --llvm"
#args[${#args[@]}]="--ubuntu 20.04 --llvm"
#args[${#args[@]}]="--ubuntu 24.04 --mpi"
#args[${#args[@]}]="--ubuntu 22.04 --mpi"
#args[${#args[@]}]="--ubuntu 20.04 --mpi"
#args[${#args[@]}]="--ubuntu 24.04 --cuda"
#args[${#args[@]}]="--ubuntu 22.04 --cuda"
#args[${#args[@]}]="--ubuntu 20.04 --cuda"
#args[${#args[@]}]="--debian 12"
#args[${#args[@]}]="--debian 11"
#args[${#args[@]}]="--debian 10"
#args[${#args[@]}]="--debian 12 --mpi"  #this one for MPI
#args[${#args[@]}]="--debian 11 --mpi"
#args[${#args[@]}]="--debian 10 --mpi"
args[${#args[@]}]="--debian 12 --cuda"  #this one for CUDA
#args[${#args[@]}]="--debian 11 --cuda"
#args[${#args[@]}]="--debian 10 --cuda"


for arg_string in "${args[@]}"; do
    tag=`echo $gitlab_project_tag"gromos-"$arg_string | sed 's/--//'g | sed 's/ /-/'g`
    dockerfile=`echo "Dockerfile_"$arg_string".txt" | sed 's/--//'g | sed 's/ /-/'g`
    echo -e "\033[32;7m$tag\033[0m"
    echo -e "\033[34;7mCreating Dockerfile backup file\033[0m"
    $(which python3) $SCRIPT $arg_string > $dockerfile
    echo -e "\033[34;7mGetting cached information\033[0m"
    $(which python3) $SCRIPT $arg_string | docker build -t $tag --cache-from $tag -
    echo -e "\033[34;7mBuilding docker image\033[0m"
    $(which python3) $SCRIPT $arg_string | docker build -t $tag -
    echo -e "\033[34;7mPushing docker imager\033[0m"
    docker push $tag
done

