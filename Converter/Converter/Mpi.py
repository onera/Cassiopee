# Interface pour MPI

import os
if os.environ.has_key('MPIRUN'): # securise import mpi4py.MPI
    if int(os.environ['MPIRUN'])>0:
        try:
            from Mpi4py import *
        except: raise ImportError("Converter:Mpi: requires mpi4py module.")
    else:
       rank = 0; size = 1; KCOMM = None
else: # try import (may fail - core or hang)
    try:
        from Mpi4py import *
    except: raise ImportError("Converter:Mpi: requires mpi4py module.")

from Distributed import _readZones, _convert2PartialTree, _convert2SkeletonTree, _readPyTreeFromPaths
import PyTree

#==============================================================================
# IN: t: full/loaded skel/partial
#==============================================================================
def center2Node(t, var=None, cellNType=0, graph=None):
    """Convert a zone or a field from centers to node."""
    if graph is None: graph = computeGraph(t, type='match')
    tl = addXZones(t, graph)
    tl = convert2PartialTree(tl)
    # print info
    zones = Internal.getZones(tl)
    print 'Rank %d has %d zones.'%(rank, len(zones))
    tl = PyTree.center2Node(tl, var, cellNType)
    tl = rmXZones(tl)
    return tl
