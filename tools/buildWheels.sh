#!/bin/bash
# Build Cassiopee wheels inside a manylinux2014 container.
# Usage (from the publish-pypi.yml workflow):
#   docker run --rm -v $GITHUB_WORKSPACE/..:/io \
#       quay.io/pypa/manylinux2014_x86_64:latest \
#       /bin/bash /io/Cassiopee/tools/build_wheels.sh
set -e

PYTHON_ABIS=(
    cp39-cp39
    cp310-cp310
    cp311-cp311
    cp312-cp312
    cp313-cp313
)

manylinux-interpreters ensure "${PYTHON_ABIS[@]}"

for py in "${PYTHON_ABIS[@]}"; do
    test -x "/opt/python/$py/bin/python" || {
        echo "ERROR: Missing Python interpreter: $py"
        exit 1
    }
done

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
    local pythonABI="$2"
    rm -rf "$CASSIOPEE/tmp"
    mkdir -p "$CASSIOPEE/tmp"

    find /tmp/ -type f \( \
        -name "cassiopee_*-${pythonABI}-*.whl" -o \
        -name "cassiopee_*-none-any.whl" \
    \) -exec mv {} "$CASSIOPEE/tmp/" \;

    if ! compgen -G "$CASSIOPEE/tmp/*-${pythonABI}-*.whl" >/dev/null; then
        echo "ERROR: No Cassiopee wheels found for $pythonABI"
        exit 1
    fi

    mkdir -p "$CASSIOPEE/wheels/$intSize/wheelhouse"
    cd "$CASSIOPEE/wheels/$intSize"

    for kwheel in "$CASSIOPEE"/tmp/*.whl; do
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

buildPythonWheels() {
    local intSize="$1"
    local py="$2"
    local python="/opt/python/$py/bin/python"

    echo "=========================================="
    echo "Building $intSize with $py"
    echo "=========================================="

    if [[ ! -x "$python" ]]; then
        echo "ERROR: Python interpreter not found: $python"
        exit 1
    fi

    "$python" --version
    "$python" -m pip --version

    "$python" -m pip install --upgrade pip
    "$python" -m pip install wheel auditwheel setuptools scons

    if [[ "$py" == "cp38-cp38" ]]; then
        "$python" -m pip install --only-binary=numpy "numpy>=1.23.3,<2"
    else
        "$python" -m pip install --only-binary=numpy "numpy>=2.0,<2.3"
    fi

    MPICC=mpicc CFLAGS="-std=c99" \
        "$python" -m pip install \
        --no-cache-dir --no-binary=mpi4py mpi4py

    export CFLAGS="-std=c11"

    "$python" -c \
        'import sys, numpy; print(sys.executable); print(numpy.__version__)'

    . "$CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8"

    # The environment script may modify PATH or PYTHONEXE.
    # Restore the selected interpreter for the Cassiopee installation.
    export PATH="/opt/python/$py/bin:$BASE_PATH"
    export PYTHONEXE="$python"
    hash -r

    cd "$CASSIOPEE"
    ./clean -f
    ./install || {
        echo "Cassiopee $intSize $py installation failed"
        exit 1
    }

    repairWheels "$intSize" "$py"
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

export PATH=/usr/lib64/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
BASE_PATH="$PATH"

for py in "${PYTHON_ABIS[@]}"; do
    buildPythonWheels int32 "$py"
done

buildMetaWheel int32

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
cd "$CASSIOPEE"

export PATH=/usr/lib64/openmpi/bin:$PATH
export LD_LIBRARY_PATH=/usr/lib64/openmpi/lib:$LD_LIBRARY_PATH
BASE_PATH="$PATH"

for py in "${PYTHON_ABIS[@]}"; do
    buildPythonWheels int64 "$py"
done

buildMetaWheel int64
