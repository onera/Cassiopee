#!/bin/bash
# Build Cassiopee wheels inside a manylinux2014 container.
# Usage (from the publish-pypi.yml workflow):
#   docker run --rm -v $GITHUB_WORKSPACE/..:/io \
#       quay.io/pypa/manylinux2014_x86_64:latest \
#       /bin/bash /io/Cassiopee/tools/build_wheels.sh
set -e

export PATH=/opt/python/cp311-cp311/bin:$PATH

yum update -y && yum install -y \
    gcc \
    gcc-c++ \
    gcc-gfortran \
    openmpi \
    openmpi-devel \
    hdf5 \
    hdf5-devel \
    tk \
    mesa-libGL \
    mesa-libGL-devel \
    mesa-libGLU \
    mesa-libGLU-devel \
    mesa-libOSMesa \
    mesa-libOSMesa-devel \
    libX11-devel \
    libXext-devel \
    libXmu-devel \
    libXi-devel \
    libXrender-devel \
    libpng-devel \
    zlib-devel \
    xorg-x11-server-devel \
    OCE-devel \
    OCE-draw \
    OCE-foundation \
    OCE-modeling \
    OCE-ocaf \
    OCE-visualization

pip3 install \
    wheel \
    auditwheel \
    setuptools \
    scons
pip3 install --only-binary=numpy "numpy>=2.0,<2.3"

export PATH=/usr/lib64/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
export CFLAGS="-std=c99"
pip3 install mpi4py
export CFLAGS="-std=c11"
python --version && python -m pip list

export CASSIOPEE=/io/Cassiopee
cd "$CASSIOPEE/Cassiopee"
sed -i 's/Modeler//g' MODULES
sed -i 's/OCC//g' MODULES

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Collect the built wheels from /tmp, repair them, and put the result
# in $CASSIOPEE/wheels/<int32|int64>/wheelhouse/
repairWheels() {
    local intSize="$1"
    rm -rf "$CASSIOPEE/tmp"
    mkdir -p "$CASSIOPEE/tmp/"
    find /tmp/ -type f -name 'cassiopee_*.whl' -exec mv {} "$CASSIOPEE/tmp/" \;
    mkdir -p "$CASSIOPEE/wheels/$intSize/wheelhouse/"
    cd "$CASSIOPEE/wheels/$intSize"
    for kwheel in $(find "$CASSIOPEE/tmp/" -name "*.whl"); do
        if [[ "$kwheel" == *-none-any.whl ]]; then
            cp "$kwheel" wheelhouse/
        else
            auditwheel repair "$kwheel" || exit 1
        fi
    done
    ls wheelhouse
}

# Build the meta-package
buildMetaWheel() {
    local intSize="$1"
    mkdir -p "$CASSIOPEE/wheels/$intSize/meta"
    python -m pip wheel \
        --no-deps \
        --no-build-isolation \
        --wheel-dir "$CASSIOPEE/wheels/$intSize/meta" \
        "$CASSIOPEE"
}

# Append -i8 to the project name and to every cassiopee* dependency in a
# pyproject.toml.
renamePackagesInt64() {
    local file="$1"
    sed -i -E \
        -e 's/^(name = "cassiopee[A-Za-z0-9_.-]*)"/\1-i8"/' \
        -e 's/^([[:space:]]*"cassiopee[A-Za-z0-9_.-]*)([^A-Za-z0-9_.-])/\1-i8\2/' \
        "$file"
    #grep -nE '^name|"cassiopee' "$file"
}

# ---------------------------------------------------------------------------
# Installation int32
# ---------------------------------------------------------------------------
export MACHINE=azure
. "$CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8"
cd "$CASSIOPEE"
./install || { echo "Cassiopee int32 install failed"; exit 1; }
buildMetaWheel int32
repairWheels int32

# ---------------------------------------------------------------------------
# Installation int64 (-i8)
# ---------------------------------------------------------------------------
cd "$CASSIOPEE/Cassiopee"
FREE_MODULES=$(head -1 MODULES | cut -d '=' -f 2 | tr -d "'")

renamePackagesInt64 "$CASSIOPEE/pyproject.toml"
for mod in $FREE_MODULES; do
    renamePackagesInt64 "$CASSIOPEE/Cassiopee/$mod/pyproject.toml"
done

export MACHINE=azure_i8
. "$CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8"
cd "$CASSIOPEE"
./install || { echo "Cassiopee int64 install failed"; exit 1; }
buildMetaWheel int64
repairWheels int64
