#  ==========================================================================
#  Project: elsA - Copyright (c) 1998-2005 by ONERA
#  Type   : Makefile     
#  File   : Make_obj.mk
#  Chrono : $Date: 2003/02/04 15:17:13 $
#           1.0 09/06/1998 elsA    0.1 
#  ==========================================================================
#
E_PYINCLUDE1 = -I$(E_PPREFIX)/include/python$(E_PYVERSION)
E_PYINCLUDES = $(E_PYINCLUDE1) -I$(E_PPREFIX1)/include/python$(E_PYVERSION)
E_PYLIBS     = $(E_PPREFIX)/lib/python$(E_PYVERSION)/config
E_PYINC      = -I$(E_PYLIBS) $(E_PYINCLUDES)
E_PYLIB      = -L$(E_PYLIBS) -lpython$(E_PYVERSION)

# Object files to put into library
# 
E_LIBOBJLIST=\
CDiff.o

E_OTHERLIBS=

OTHEROBJS=

E_EXECS=\
CDiff.x
#
# --------------------------------------------------------LAST LINE------
