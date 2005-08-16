#!/usr/bin/bash

BUILD=`date "+%d%m%y"`
VERSION=`awk -F= 'BEGIN{x=0} {if (x==1) printf "."; x=1; printf $2} END{printf "\n"}' ../VERSION`
NIGHT=~/RELEASE/nightly
NIGHTLOG=${NIGHT}/night_${BUILD}.log

OS=`uname -s`"_"`uname -m`
source options.${OS}

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
    source prepare.sh
}

distcheck(){
    echo "running distcheck"
    source distcheck.sh
}

pkg_solaris(){
    echo "creating a solaris package"
    source pkg_solaris.sh
}

pkg_debian(){
    echo "creating a debian package"
    source pkg_debian.sh
}

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
