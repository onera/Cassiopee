#!/bin/bash

export http_proxy=http://proxy.onera:80
export https_proxy=$http_proxy

version=2.4
folder=/scratch/jchiche/Anaconda/conda-bld/src_cache

mkdir -p $folder

address=http://elsa.onera.fr/Cassiopee/Download/KCore-${version}F.tar.gz

wget "http://elsa.onera.fr/Cassiopee/Download/KCore-2.4F.tar.gz"
