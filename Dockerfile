# Use the official Ubuntu base image
FROM ubuntu:latest

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Set the default Python version to 3.12.4, can be passed as an argument
ARG PYTHON_VERSION=3.12.4

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    scons \
    gcc \
    g++ \
    gfortran \
    libopenmpi-dev \
    libhdf5-openmpi-dev \
    mesa-common-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libosmesa6-dev \
    xorg-dev \
    time
    
# Install system dependencies required by pyenv
RUN apt-get update && apt-get install -y \
    git \
    libssl-dev \
    zlib1g-dev \
    libbz2-dev \
    libreadline-dev \
    libsqlite3-dev \
    llvm \
    libncurses5-dev \
    libncursesw5-dev \
    xz-utils \
    tk-dev \
    libffi-dev \
    liblzma-dev \
    python3-openssl \
    libyaml-dev \
    libgdbm-dev \
    libdb-dev \
    libpcap-dev \
    libmysqlclient-dev

# Set environment variables for compatibility with amd64 and arm64 architectures
ARG TARGETARCH
ENV MPI_HOME=/usr/lib/${TARGETARCH}-linux-gnu/openmpi
ENV CPATH=${MPI_HOME}/include:${CPATH:-}
ENV LIBRARY_PATH=${MPI_HOME}/lib:${LIBRARY_PATH:-}

# Set environment variables for PyEnv
ENV PYENV_ROOT=/root/.pyenv
ENV PATH=$PYENV_ROOT/shims:$PYENV_ROOT/bin:$PATH

# Install pyenv and the requested Python version
RUN curl https://pyenv.run | bash \
    && pyenv update \
    && pyenv install ${PYTHON_VERSION} \
    && pyenv global ${PYTHON_VERSION} \
    && pyenv rehash

# Shadow python3 with pyenv
RUN pyenv global ${PYTHON_VERSION} \
    && ln -sf $(pyenv which python) /usr/local/bin/python \
    && ln -sf $(pyenv which python) /usr/local/bin/python3 \
    && ln -sf $(pyenv which pip) /usr/local/bin/pip \
    && ln -sf $(pyenv which pip) /usr/local/bin/pip3

ENV PIP_ROOT_USER_ACTION=ignore

# Install Python dependencies using pip
RUN python -m pip install --upgrade pip setuptools wheel build
RUN python -m pip install \
        numpy>=1.23.3 \
        mpi4py>=3.1.3 \
        scons>=4.4.0 \
        matplotlib>=3.8.0

# Verify the version of Python and numpy
RUN python3 --version && pip3 --version
RUN python3 -m pip list \
    && python3 -c "import numpy; print(numpy.__version__)"

# Install opencascade and remove cache files after installing packages
RUN apt-get update && apt-get install -y \
    gedit \
    meld \
    cmake \
    git \
    libocct-foundation-dev \
    libocct-modeling-algorithms-dev \
    libocct-data-exchange-dev \
    libocct-modeling-data-dev \
    libocct-draw-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install latest versions of HDF5 and CGNS tools
# Adapted from James Clark's Dockerfile
# https://github.com/jameshclrk/cgns-hdf5-parmetis-docker/blob/master/Dockerfile
ENV HDF5_VERSION=1.10.5
ENV CGNS_VERSION=4.4.0
ENV HDF5_INSTALL_DIR="/opt/hdf5"
ENV CGNS_INSTALL_DIR="/opt/cgns"

RUN curl -fSL "https://support.hdfgroup.org/ftp/HDF5/current/src/hdf5-$HDF5_VERSION.tar.gz" -o hdf5.tar.gz \
	&& mkdir -p /usr/src/hdf5 \
	&& tar -xf hdf5.tar.gz -C /usr/src/hdf5 --strip-components=1 \
	&& rm hdf5.tar.gz* \
	&& cd /usr/src/hdf5 \
	&& CC=mpicc ./configure --prefix=${HDF5_INSTALL_DIR} --enable-parallel \
	&& make -j"$(nproc)" \
	&& make install \
	&& rm -rf /usr/src/hdf5
	
ENV HDF5_DIR="${HDF5_INSTALL_DIR}"

RUN curl -fSL "https://github.com/CGNS/CGNS/archive/refs/tags/v$CGNS_VERSION.tar.gz" -o cgns.tar.gz \
	&& mkdir -p /usr/src/cgns \
	&& tar -xf cgns.tar.gz -C /usr/src/cgns --strip-components=1 \
	&& rm cgns.tar.gz* \
	&& cd /usr/src/cgns \
	&& mkdir build && cd build \
	&& PATH="${HDF5_INSTALL_DIR}/bin:$PATH" cmake -DCGNS_ENABLE_HDF5=ON -DCMAKE_INSTALL_PREFIX:PATH=${CGNS_INSTALL_DIR} -DHDF5_NEED_MPI=ON .. \
	&& make -j"$(nproc)" \
	&& make install \
	&& rm -rf /usr/src/cgns

ENV LD_LIBRARY_PATH="${CGNS_INSTALL_DIR}/lib:${HDF5_INSTALL_DIR}/lib:${LD_LIBRARY_PATH}" PATH="${CGNS_INSTALL_DIR}/bin:${HDF5_INSTALL_DIR}/bin:$PATH"
ENV CGNS_DIR="${CGNS_INSTALL_DIR}"
    
# Set environment variables for Cassiopee
ENV CASSIOPEE=/Cassiopee
#ENV MACHINE=ubuntu
RUN case "${TARGETARCH}" in \
        amd64) echo "MACHINE=ubuntu" ;; \
        arm64) echo "MACHINE=ubuntu_arm64" ;; \
        *) echo "Unsupported architecture: ${TARGETARCH}" && exit 1 ;; \
    esac > /etc/cassiopee-machine.env

# Do not prevent mpirun to run as root
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Set the working directory in the container
WORKDIR $CASSIOPEE

# Copy the current directory contents into the container
COPY . $CASSIOPEE

# Change the default shell to be the bash shell
SHELL ["/bin/bash", "-c"]

# Source environment and run install script
RUN . /etc/cassiopee-machine.env && echo "Using MACHINE=${MACHINE}" \
    && . $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8 \
    && export PRODMODE=3 \
    && export OMP_NUM_THREADS=2 \
    && cd $CASSIOPEE/Cassiopee \
    && ./install

# Default command to run the application: start an interactive shell session
ENTRYPOINT ["/bin/bash", "-i", "-c", "source /etc/cassiopee-machine.env && source $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8 && exec /bin/bash -i"]
