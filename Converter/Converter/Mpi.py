# Interface pour MPI

import os
if 'MPIRUN' in os.environ: # si MPIRUN=0, force sequentiel
    if int(os.environ['MPIRUN'])>0:
        try:
            from .Mpi4py import *
        except: raise ImportError("Converter:Mpi: requires mpi4py module.")
    else:
       rank = 0; size = 1; KCOMM = None
else: # try import (may fail - core or hang)
    try:
        from .Mpi4py import *
    except:
        rank = 0; size = 1; KCOMM = None
        from .Distributed import setProc, _setProc, getProc, getProcDict, convertFile2SkeletonTree, computeGraph
        def barrier(): return
        print("Warning: Converter:Mpi: mpi4py is not available. Sequential behaviour.")

from .Distributed import _readZones, _convert2PartialTree, _convert2SkeletonTree, _readPyTreeFromPaths
from . import PyTree

#==============================================================================
# IN: t: full/loaded skel/partial
#==============================================================================
def center2Node(t, var=None, cellNType=0, graph=None):
    """Convert a zone or a field from centers to node."""
    if graph is None: graph = computeGraph(t, type='match')
    tl = addXZones(t, graph)
    tl = convert2PartialTree(tl)
    # print info
    #zones = Internal.getZones(tl)
    #print 'Rank %d has %d zones.'%(rank, len(zones))
    tl = PyTree.center2Node(tl, var, cellNType)
    tl = rmXZones(tl)
    return tl
