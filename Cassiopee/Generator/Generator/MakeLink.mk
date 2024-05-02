#!/bin/sh
# ========================================================================
# Creation de liens dans le répertoire obj des fichiers f90 à compiler.
# ========================================================================
#
RSOURCE=$1
SOURCE=$2
DTARGET=$3

export SOURCE
export DTARGET

DSOURCE=`dirname ${SOURCE}`
FSOURCE=`basename ${SOURCE}`
FTARGET=${DTARGET}/${FSOURCE}

if [ -f "${DTARGET}/${FSOURCE}" ]
then
 :
else
  ln -sf "${RSOURCE}/${SOURCE}" "${DTARGET}/${FSOURCE}"
fi
