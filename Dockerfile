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
    cmake \
    git \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
#    tcl \
#    tk \
#    doxygen \
#    libfreetype6-dev \
#    libxmu-dev \
#    libxi-dev \
#    libxext-dev \
#    libocct-ocaf-7.6t64 \
    
# Set environment variables
#ENV PATH=/OCCT/bin:$PATH
#ENV LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/
ENV CASSIOPEE=/Cassiopee
ENV MACHINE=ubuntu

# Set the working directory in the container
WORKDIR $CASSIOPEE

# Copy the current directory contents into the container
COPY . $CASSIOPEE

# Exclude the OCC module temporarily
RUN echo -e "FREEMODULES='KCore XCore Converter Geom Transform Generator Post Initiator Connector Distributor2 Dist2Walls RigidMotion Compressor Modeler Intersector Apps CPlot'\nexport FREEMODULES\nFULLMODULES='KCore XCore Converter Geom Transform Generator Post Initiator Connector Distributor2 Dist2Walls RigidMotion Compressor Modeler Intersector Apps CPlot'\nexport FULLMODULES\nOTHERS=''" > $CASSIOPEE/Cassiopee/MODULES

# Source environment and run install script
RUN . $CASSIOPEE/Cassiopee/Envs/sh_Cassiopee_r8 \
    && cd $CASSIOPEE/Cassiopee \
    && ./install

# Define the default command to run the application: start an interactive shell
# session
CMD ["bash"]
