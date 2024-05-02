#!/bin/sh
# ========================================================================
# Enleve les répertoires pour ne garder que les noms de fichier
# ========================================================================
#
OBJST=$1

OBJS=`echo ${OBJST} | sed -e "s;[A-Za-z1-90_]*/;;g"`
echo ${OBJS}
