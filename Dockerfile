# Use the official Ubuntu base image
FROM ubuntu:latest

# Set non-interactive mode for apt-get
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3-dev \
    python3-numpy \
    python3-distutils-extra \
    python3-pip \
    scons \
    gcc \
    g++ \
    gfortran \
    libopenmpi-dev \
    python3-mpi4py \
    libhdf5-openmpi-dev \
    python3-tk \
    mesa-common-dev \
    libgl1-mesa-dev \
    libglu1-mesa-dev \
    libosmesa6-dev \
    xorg-dev \
    time
    
ENV PIP_ROOT_USER_ACTION=ignore

# Install opencascade and remove cache files after installing packages
RUN apt-get update && apt-get install -y \
    build-essential \
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
ARG HDF5_VERSION=1.10.5
ARG CGNS_VERSION=4.4.0
ARG HDF5_INSTALL_DIR="/opt/hdf5"
ARG CGNS_INSTALL_DIR="/opt/cgns"

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
    
# Set environment variables
ENV CASSIOPEE=/Cassiopee
ENV MACHINE=ubuntu

# Do not prevent mpirun to run as root
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

# Set the working directory in the container
WORKDIR $CASSIOPEE

# Copy the current directory contents into the container
COPY . $CASSIOPEE

# Source environment and run install script
RUN . $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8 \
    && cd $CASSIOPEE/Cassiopee \
    && ./install

# Change the default shell to be the bash shell
SHELL ["/bin/bash", "-c"] 

# Define the default command to run the application: start an interactive shell
# session
ENTRYPOINT . $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8 && /bin/bash
