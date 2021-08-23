# Interface pour MPI

import os
if 'MPIRUN' in os.environ: # si MPIRUN=0, force sequentiel
    if int(os.environ['MPIRUN'])>0:
        try: from .Mpi4py import *
        except: raise ImportError("Converter:Mpi: requires mpi4py module.")
    else:
       rank = 0; size = 1; KCOMM = None; COMM_WORLD = None
       from .Distributed import setProc, _setProc, getProc, getProcDict, convertFile2SkeletonTree, computeGraph, readZones, convert2PartialTree, convert2SkeletonTree, readPyTreeFromPaths
       def barrier(): return
       def bcast(a, root=0): return a
       def send(a, dest=0, tag=0): return None
       def recv(source=0, tag=0): return None # pb here
       def sendRecv(a, source=0, dest=0): return []
       def sendRecvC(a, source=0, dest=0): return []
       def seq(F, *args): F(*args)
       print("Warning: Converter:Mpi: Sequential behaviour is forced by MPIRUN=0.")
 
else: # try import (may fail - core or hang)
    try: from .Mpi4py import *
    except:
        rank = 0; size = 1; KCOMM = None; COMM_WORLD = None
        from .Distributed import setProc, _setProc, getProc, getProcDict, convertFile2SkeletonTree, computeGraph, readZones, convert2PartialTree, convert2SkeletonTree, readPyTreeFromPaths
        def barrier(): return
        def bcast(a, root=0): return a
        def send(a, dest=0, tag=0): return None
        def recv(source=0, tag=0): return None # pb here
        def sendRecv(a, source=0, dest=0): return []
        def sendRecvC(a, source=0, dest=0): return []
        def seq(F, *args): F(*args)
        print("Warning: Converter:Mpi: mpi4py is not available. Sequential behaviour.")

from .Distributed import _readZones, _convert2PartialTree, _convert2SkeletonTree, _readPyTreeFromPaths
from . import PyTree
from . import Internal

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

# addGhostCells parallel
# si modified=[], transferts tous les champs
# si modified=None, transfert aucun champ
def _addGhostCells(t, b, d, adaptBCs=0, modified=[], fillCorner=1):
    
    if modified == []: # all
        variables = PyTree.getVarNames(t, excludeXYZ=True, loc='nodes')[0]
        variables += PyTree.getVarNames(t, excludeXYZ=True, loc='centers')[0]
    elif modified is None: variables = []; modified = []
    else: variables = modified

    _addMXZones(t, depth=2, variables=variables, noCoordinates=False, 
                keepOldNodes=False)
    print("%d: addGC(max): Nblocs=%d, NPts(M)=%g"%(rank,len(Internal.getZones(t)), PyTree.getNPts(t)*1./1.e6), flush=True)
    Internal._addGhostCells(t, t, d, adaptBCs, modified, fillCorner)
    _rmMXZones(t)

    # ancienne version utilisant addXZones
    #graph = computeGraph(t, type='match', reduction=True)
    #_addXZones(t, graph, variables=variable, noCoordinates=False, 
    #           zoneGC=False, keepOldNodes=False)
    #print("%d: addGC(max): Nblocs=%d, NPts(M)=%g"%(rank,len(Internal.getZones(t)), PyTree.getNPts(t)*1./1.e6), flush=True)
    #Internal._addGhostCells(t, t, d, adaptBCs, modified, fillCorner)
    #_rmXZones(t)
    return None