#!/usr/bin/env bash

# This file contains common settings and functions used by conda build.sh scripts.
# Another file contains the settings that may be changed by the user that builds
# the recipes: build-settings.sh
# This file is automatically sourced at the end of the current one.

###############################################################################
# Define some sane bash settings
# Exit on error. Append "|| true" if you expect an error.
set -o errexit
# Exit on error inside any functions or subshells.
set -o errtrace
# Do not allow use of undefined vars. Use ${VAR:-} to use an undefined VAR
set -o nounset
# Catch the error in case mysqldump fails (but gzip succeeds) in `mysqldump |gzip`
set -o pipefail
# Turn on traces, useful while debugging but commented out by default
# set -o xtrace

###############################################################################
# Define some installation paths
export PYTHON_INSTALL_PATH=${SP_DIR}
export LIB_INSTALL_PATH=${PREFIX}/lib
export INCLUDE_INSTALL_PATH=${PREFIX}/include
export RUNTIME_ENV_INSTALL_PATH=${PREFIX}/etc/conda
export DATA_INSTALL_PATH=${PREFIX}/share

###############################################################################
# Remove .svn directories
function remove_dot_svn_directories () {
  find . -type d -name .svn -exec rm -rf {} +
}

###############################################################################
# Copy files from a source to a target
# $1: path to installation directory
# $2, $3, ...: files to install
function copy_to () {
  install_with "cp -r" "${@}"
  return $?
}

###############################################################################
# Move files from a source to a target
# $1: path to installation directory
# $2, $3, ...: files to install
function move_to () {
  install_with "mv" "${@}"
  return $?
}

###############################################################################
# Install files from a source to a target with a command
# $1: command used puting the files (mv, cp -r, ...)
# $2: path to installation directory
# $3, $4, ...: files to install
function install_with () {
  local -r command="${1}"
  local -r target="${2}"
  local -r sources="${@:3}"

  mkdir -p ${target}
  ${command} ${sources} ${target}
}

###############################################################################
# Return whether a value is in a list of values
# $1: value to be checked
# $2, $3, ...: values to be checked against
function is_in () {
  local -r value="${1}"
  local -r items="${@:2}"

  for item in ${items}; do
    if [[ ${value} == ${item} ]]; then
      return 0
    fi
  done

  return 1
}

###############################################################################
# Replace placeholder in runtime environment scripts
# $1: placeholder (without @ padding) to be replaced
# $2: value to replace with
function replace_placeholder_in_runtime_env_scripts () {
  local -r placeholder="${1}"
  local -r value="${2}"

  sed -i \
    "s|@${placeholder}@|${value}|g" \
    $RUNTIME_ENV_INSTALL_PATH/activate.d/*
}

###############################################################################
# Set MPI backend compiler environment variable
# $1: kind of backend, one of: CC, CXX, F77, F90, FC
# $2: backend compiler executable
function set_mpi_back_end_compiler () {
  local -r kind=${1}
  local -r compiler=${2}

  # check kind
  local -r kinds="CC CXX F77 F90 FC"
  if ! is_in ${kind} ${kinds}; then
    echo 'Bad MPI compiler backend kind: '${kind}' is not one of '${kinds}
    exit 1
  fi

  # platform mpi
  export MPI_${kind}=${compiler}
  # open mpi
  export OMPI_${kind}=${compiler}
  # intel mpi
  export I_MPI_${kind}=${compiler}
  # mpich
  export MPICH_${kind}=${compiler}
}

###############################################################################
# Export a variable with compiler flags from predefined flags.
#
# The predefined flags are set in environment variables for different
# compilers and build types like:
#  * FLAGS_GNU_DBG, FLAGS_GNU_OPT
#  * FLAGS_INTEL_DBG, FLAGS_INTEL_OPT
#
# $1: kind of compiler, one of: GNU, INTEL
# $2: variable to be exported
function export_compiler_flags_to () {
  local -r variable=$1
  local -r kind=$2

  # check BUILD_TYPE
  local -r build_type=${BUILD_TYPE:-'none'} # deal with unset or empty variable
  local -r build_types='OPT DBG'
  if ! is_in ${build_type} ${build_types}; then
    echo 'Bad BUILD_TYPE: '${build_type}' is not one of '${build_types}
    exit 1
  fi

  local -r source_var_name=FLAGS_${kind}_${BUILD_TYPE}

  # check the predefined flags exist
  if [[ -z ${!source_var_name:-} ]]; then
    echo "${source_var_name} variable is not defined"
    exit 1
  fi

  export ${variable}="${!source_var_name}"
}

###############################################################################
function install_flowsim_plugin () {

  remove_dot_svn_directories
 
  export CC=mpicc
  export CXX=mpicc
  
  # use our build config file
  rm -rf config  
  cp -r $RECIPE_DIR/config .

  export FLOWSIM_PREFIX=$PREFIX
 
  rm -rf build
  python setup.py distclean
  python setup.py build
  python setup.py install

  # install python packages
  #move_to \
  #  $PYTHON_INSTALL_PATH  \
  #  $PREFIX/py/FS*
  move_to \
    $PYTHON_INSTALL_PATH  \
    $PREFIX/FS*

  # install headers
  echo "==================================================="
  module_rel_paths=$(find . -type d -path "./FS*" -prune)
  for module_rel_path in $module_rel_paths; do
    copy_to \
      $INCLUDE_INSTALL_PATH/$module_rel_path \
      $module_rel_path/include/*
  done
}
