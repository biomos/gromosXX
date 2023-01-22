#!/usr/bin/env python3

"""Dockerfile generator for CI testing images.

Authors:
    * Jakob Liu <koki@cloningcompany.at>
    * Yerko Escalona <yerko.escalona@boku.ac.at>

Usage:
    $ python3 gromos_dockerfile_generator.py --help
    $ python3 gromos_dockerfile_generator.py --ubuntu 16.04 --cuda 10.2
    $ python3 gromos_dockerfile_generator.py

"""

import argparse
import collections
import hpccm
from hpccm.primitives import baseimage

# list of packages to install
_common_packages = [
    'apt-utils',
    'automake',
    'autotools-dev',
    'build-essential',
    'bzip2',
    'ca-certificates',
    'curl',
    'git',
    'libfftw3-dev',
    'libgsl-dev',
    'libtool',
    'openssh-client',
    'software-properties-common',
    'wget',
    'zlib1g-dev',
    'ccache',
    'python3-pytest']

_pip = ['numpy', 'PyYAML']

_clang = ['clang',
          'libomp5',
          'libomp-dev']

_mpi = [
    'libopenmpi-dev',
    'libfftw3-3',
    'libfftw3-mpi3',
    'libfftw3-mpi-dev',
    'openmpi-bin',]

_cuda = [
    'nvidia-cuda-toolkit',
    'nvidia-driver-460',]

# parser

parser = argparse.ArgumentParser(description='GROMOS CI image options')

compiler_group = parser.add_mutually_exclusive_group()
compiler_group.add_argument('--gcc', type=int, nargs='?', const=7, default=None)
compiler_group.add_argument('--llvm', type=int, nargs='?', const=7, default=None)

distro_group = parser.add_mutually_exclusive_group()
distro_group.add_argument('--ubuntu', type=str, nargs='?', const='20.04', default=None,
                          choices=['16.04', '18.04', '20.04', '22.04'])
distro_group.add_argument('--debian', type=str, nargs='?', const='9', default=None,
                          choices=['9', '10', '11'])
distro_group.add_argument('--centos', type=str, nargs='?', const='9', default=None)

parser.add_argument('--mpi', type=str, nargs='?', const="openmpi", default=None)
parser.add_argument('--cuda', type=str, nargs='?', const="10.2", default=None)


parser.add_argument('--tag', type=str, nargs='?', default=None)

# functions

def set_linux_distro(args):
    """ Set the linux distribution according to hpccm. """
    if args.ubuntu is not None:
        ubuntu_distros = {'16.04': 'ubuntu16',
                          '18.04': 'ubuntu18',
                          '20.04': 'ubuntu20',
                          '22.04': 'ubuntu22'}
        linux_distro = ubuntu_distros[args.ubuntu]

    # ubuntu uses the debian's filenames 
    elif args.debian is not None:
        debian_distros = {'9': 'ubuntu16',
                          '10': 'ubuntu20',
                          '11': 'ubuntu22'}
        linux_distro = debian_distros[args.debian]

    return linux_distro

def get_image_tag(args):
    """ Get the correct image tag for looking at DockerHub. """

    if args.ubuntu is not None:
        if args.cuda is not None:
            image_tag = 'nvidia/cuda:' + args.cuda + '-devel' + '-ubuntu' + args.ubuntu
        else:
            image_tag = "ubuntu:" + args.ubuntu
    elif args.debian is not None:
        if args.cuda is not None:
            raise RuntimeError('debian + cuda is not implemented')
        else:
            image_tag = "debian:" + args.debian

    return image_tag

def get_compiler(args):
    """ Set parameters for gcc or clang. """

    if args.gcc is not None:
        compiler = hpccm.building_blocks.gnu(extra_repository=True,
                                             version=args.gcc,
                                             fortran=False)
    elif args.llvm is not None:
        compiler = hpccm.building_blocks.packages(ospackages=_clang)
        #compiler = hpccm.building_blocks.llvm(version=args.llvm, upstream=True)

    else:
        compiler = "# using default compiler" #print(compiler)
    return compiler


def get_mpi(args):
    """ Install MPI packages. """

    if args.mpi is not None:
        return hpccm.building_blocks.packages(ospackages=_mpi)
    else:
        return "# default MPI configuration"



def add_base_stage(name, input_args):
    """ Add base stage """

    output_stages = {}
    output_stages[name] = hpccm.Stage()
    output_stages[name] += hpccm.primitives.baseimage(image=get_image_tag(input_args), _distro=set_linux_distro(input_args))
    output_stages[name] += hpccm.building_blocks.packages(ospackages=_common_packages)
    output_stages[name] += hpccm.building_blocks.pip(packages=_pip, pip='pip3')
    output_stages[name] += get_compiler(input_args)
    output_stages[name] += get_mpi(input_args)
    return output_stages


def build_stages(args):
    """ Build stages """

    stages = {}

    stages.update(add_base_stage(name="base", input_args=args))

    #  space for new stages

    for build_stage in stages.values():
        if build_stage is not None:
            yield build_stage


if __name__ == "__main__":
    args = parser.parse_args()

    hpccm.config.set_linux_distro(set_linux_distro(args))
    recipe = build_stages(args)

    for stage in recipe:
        print(stage)

