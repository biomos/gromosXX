"""
This file is part of GROMOS.

Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
See <https://www.gromos.net> for details.

GROMOS is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

#!/usr/bin/env python3

"""Dockerfile generator for CI testing images.

Authors:
    * Jakob Liu <koki@cloningcompany.at>
    * Yerko Escalona <yerko.escalona@boku.ac.at>
    * Drazen Petrov <drazen.petrov@boku.ac.at>

Usage:
    $ python3 gromos_dockerfile_generator.py --help
    $ python3 gromos_dockerfile_generator.py --ubuntu 22.04 --cuda
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
    'cmake',
    'bzip2',
    'ca-certificates',
    'curl',
    'git',
    'gcc-11',
    'g++-11',
    'gnupg2',
    'libfftw3-dev',
    'libgsl-dev',
    'libtool',
    'openssh-client',
    'software-properties-common',
    'wget',
    'zlib1g-dev',
    'ccache',
    'python3-pytest',
    'python3-numpy',
    'python3-yaml']

_pip = ['numpy', 'PyYAML']

_clang = ['clang',
          'libomp5',
          'libomp-dev']

_mpi = [
    'libopenmpi-dev',
    'libfftw3-mpi3',
    'libfftw3-mpi-dev',
    'openmpi-bin',]
    #'libfftw3-3', #removed this package from the list because non-existent in debian12 anymore

_cuda = [
    'cuda-toolkit-11-7',
    'cuda-libraries-11-7',
    'cuda-libraries-dev-11-7',
    'cuda-cudart-dev-11-7',]

# parser

parser = argparse.ArgumentParser(description='GROMOS CI image options')

compiler_group = parser.add_mutually_exclusive_group()
compiler_group.add_argument('--gcc', type=int, nargs='?', const=7, default=None)
compiler_group.add_argument('--llvm', type=int, nargs='?', const=7, default=None)

distro_group = parser.add_mutually_exclusive_group()
distro_group.add_argument('--ubuntu', type=str, nargs='?', const='22.04', default=None,
                          choices=['20.04', '22.04', '24.04'])
distro_group.add_argument('--debian', type=str, nargs='?', const='12', default=None,
                          choices=['10', '11', '12'])
distro_group.add_argument('--centos', type=str, nargs='?', const='9', default=None)

parser.add_argument('--mpi', type=str, nargs='?', const="openmpi", default=None)
parser.add_argument('--cuda', type=str, nargs='?', const="11.7", default=None)


parser.add_argument('--tag', type=str, nargs='?', default=None)

# functions

def set_linux_distro(args):
    """ Set the linux distribution according to hpccm. """
    if args.ubuntu is not None:
        ubuntu_distros = {'20.04': 'ubuntu20',
                          '22.04': 'ubuntu22',
                          '24.04': 'ubuntu24'}
        linux_distro = ubuntu_distros[args.ubuntu]

    # ubuntu uses the debian's filenames 
    elif args.debian is not None:
        debian_distros = {'10': 'debian10',
                          '11': 'debian11',
                          '12': 'debian12'}
        linux_distro = debian_distros[args.debian]

    return linux_distro

def get_cuda_repo_distro(args):
    if args.ubuntu:
        return {
            '20.04': 'ubuntu2004',
            '22.04': 'ubuntu2204',
            '24.04': 'ubuntu2404'
        }[args.ubuntu]
    elif args.debian:
        return {
            '10': 'debian10',
            '11': 'debian11',
            '12': 'debian11' # NVIDIA does not provide CUDA 11.7 for deb12; deb11 works
        }[args.debian]
    else:
        raise RuntimeError("No CUDA repo distro for given args")

def get_image_tag(args):
    """ Get the correct image tag for looking at DockerHub. """

    if args.ubuntu is not None:
        image_tag = "ubuntu:" + args.ubuntu
    elif args.debian is not None:
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

    # make gcc-11 / g++-11 the default compiler
    output_stages[name] += hpccm.primitives.shell(
      commands=[
        'update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-11 100 '
          '--slave /usr/bin/g++ g++ /usr/bin/g++-11',
        'update-alternatives --set gcc /usr/bin/gcc-11'
      ]
    )

# pip installation doesn't work, complains about no venv available, numpy and PyYAML now installed through apt
#    output_stages[name] += hpccm.building_blocks.pip(packages=_pip, pip='pip3')
    output_stages[name] += get_compiler(input_args)
    output_stages[name] += get_mpi(input_args)

    if input_args.cuda:
        # 1) register NVIDIA CUDA apt repo
        cuda_repo_distro = get_cuda_repo_distro(input_args)
        output_stages[name] += hpccm.primitives.shell(commands=[
            f'curl -fsSL https://developer.download.nvidia.com/compute/cuda/repos/'
                f'{cuda_repo_distro}/x86_64/3bf863cc.pub | apt-key add -',
            f'echo "deb https://developer.download.nvidia.com/compute/cuda/repos/'
                f'{cuda_repo_distro}/x86_64/ /" > /etc/apt/sources.list.d/cuda.list',
            'apt-get update'
        ])
        # 2) install the CUDA toolkit & dev-libs
        output_stages[name] += hpccm.building_blocks.packages(
            ospackages=_cuda
        )
        # 3) set CUDA env vars so that nvcc, libcudart, etc. just work
        output_stages[name] += hpccm.primitives.environment(
            variables={
                'CUDA_HOME': f'/usr/local/cuda-{input_args.cuda}',
                'PATH':      f'/usr/local/cuda-{input_args.cuda}/bin:$PATH',
                'LD_LIBRARY_PATH': 
                     f'/usr/local/cuda-{input_args.cuda}/lib64:$LD_LIBRARY_PATH'
              }
        )

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

    # tell HPCCM which distro to target
    hpccm.config.set_linux_distro(set_linux_distro(args))

    # generate all stages
    recipe = list(build_stages(args))

    # render them all into one big string
    dockerfile_contents = "\n".join(str(stage) for stage in recipe)

    # print to stdout (so you still see it on the console)
    print(dockerfile_contents)
