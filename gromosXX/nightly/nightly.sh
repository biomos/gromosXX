#!/usr/bin/bash

export NIGHTHOME=~/projects/gromosXX/nightly

cd ${NIGHTHOME}
source ${NIGHTHOME}/options.all

usage(){
	    printf "\n"
	    printf "nightly\n"
	    printf "\tprepare\n"
	    printf "\tdistcheck\n"
	    printf "\tsolaris\n"
	    printf "\tdebian\n"
	    printf "\trpm\n"
	    printf "\n"
}
prepare(){
    echo "creating a tar.gz archive"
    source ${NIGHTHOME}/prepare.sh
}

distcheck(){
    echo "running distcheck"
    source ${NIGHTHOME}/distcheck.sh
}

pkg_solaris(){
    echo "creating a solaris package"
    source ${NIGHTHOME}/pkg_solaris.sh
}

pkg_debian(){
    echo "creating a debian package"
    source ${NIGHTHOME}/pkg_debian.sh
}

if [ -z "$1" ] ; then
    usage
fi

NIGHTLOG=${NIGHT}/night_${BUILD}.log

OS=`uname -s`"_"`uname -m`
source ${NIGHTHOME}/options.${OS}

OPTIONNAME=$1
export OPTIONNAME
source ${NIGHTHOME}/options.${OPTIONNAME}
shift

VERSION=`awk -F= 'BEGIN{x=0} {if (x==1) printf "."; \
			x=1; printf $2} END{printf "\n"}' \
	 ${VERSIONFILE}`

echo "${NAME}"

if [ -z "$1" ] ; then
    usage
fi

while [ ! -z "$1" ] ; do

    case $1 in
	"prepare") prepare
	;;
	"distcheck") distcheck
	;;
	"solaris") pkg_solaris
	;;
	"debian") pkg_debian
	;;
	*) usage
	;;
    esac

    shift

done
