#!/bin/bash

SCRIPT=$PWD/gromos_dockerfile_generator.py

gitlab_project_tag=gromosx.boku.ac.at:5005/gromos/gromosxx/

#args[${#args[@]}]="--ubuntu 22.04"
#args[${#args[@]}]="--ubuntu 20.04"
#args[${#args[@]}]="--ubuntu 18.04"
#args[${#args[@]}]="--ubuntu 16.04"
#args[${#args[@]}]="--ubuntu 22.04 --llvm"
#args[${#args[@]}]="--ubuntu 20.04 --llvm"
#args[${#args[@]}]="--ubuntu 18.04 --llvm"
#args[${#args[@]}]="--ubuntu 16.04 --llvm"
#args[${#args[@]}]="--ubuntu 22.04 --mpi"
args[${#args[@]}]="--ubuntu 20.04 --mpi"
#args[${#args[@]}]="--ubuntu 18.04 --mpi"
#args[${#args[@]}]="--ubuntu 16.04 --mpi"
args[${#args[@]}]="--ubuntu 20.04 --cuda 11.6.1"
#args[${#args[@]}]="--ubuntu 20.04 --cuda 11.4.1"
#args[${#args[@]}]="--ubuntu 16.04 --cuda 10.2"
#args[${#args[@]}]="--debian 11"
#args[${#args[@]}]="--debian 10"
#args[${#args[@]}]="--debian 9"
#args[${#args[@]}]="--debian 11 --mpi"
#args[${#args[@]}]="--debian 10 --mpi"
#args[${#args[@]}]="--debian 9 --mpi"


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

