# Interface pour MPI

import os, time, timeit, sys

if 'MPIRUN' in os.environ: # si MPIRUN=0, force sequentiel
    if int(os.environ['MPIRUN'])>0:
        try: from .Mpi4py import *
        except: raise ImportError("Converter:Mpi: requires mpi4py module.")
    else:
        rank = 0; size = 1; KCOMM = None; COMM_WORLD = None
        from .Distributed import setProc, _setProc, getProc, getProcDict, getProperty, getPropertyDict, convertFile2SkeletonTree, computeGraph, splitGraph, mergeGraph, readZones, convert2PartialTree, convert2SkeletonTree, readPyTreeFromPaths
        def barrier(): return
        def bcast(a, root=0): return a
        def Bcast(a, root=0): return a
        def gather(a, root=0): return a
        def Gather(a, root=0): return a
        def gatherZones(a, root=0): return a
        def allgatherZones(a, root=0): return a
        def send(a, dest=0, tag=0): return None
        def recv(source=0, tag=0): return None # pb here
        def sendRecv(a, source=0, dest=0): return []
        def sendRecvC(a, source=0, dest=0): return []
        def reduce(a, op=None, root=0): return a
        def Reduce(a, b, op=None, root=0): return a
        def allreduce(a, op=None): return a
        def Allreduce(a, b, op=None): return a
        def seq(F, *args): F(*args)
        print("Warning: Converter:Mpi: Sequential behaviour is forced by MPIRUN=0.")
 
else: # try import (may fail - core or hang)
    try: from .Mpi4py import *
    except:
        rank = 0; size = 1; KCOMM = None; COMM_WORLD = None
        from .Distributed import setProc, _setProc, getProc, getProcDict, getProperty, getPropertyDict, convertFile2SkeletonTree, computeGraph, splitGraph, mergeGraph, readZones, convert2PartialTree, convert2SkeletonTree, readPyTreeFromPaths
        def barrier(): return
        def bcast(a, root=0): return a
        def Bcast(a, root=0): return a
        def gather(a, root=0): return a
        def Gather(a, root=0): return a
        def gatherZones(a, root=0): return a
        def allgatherZones(a, root=0): return a
        def send(a, dest=0, tag=0): return None
        def recv(source=0, tag=0): return None # pb here
        def sendRecv(a, source=0, dest=0): return []
        def sendRecvC(a, source=0, dest=0): return []
        def reduce(a, op=None, root=0): return a
        def Reduce(a, b, op=None, root=0): return a
        def allreduce(a, op=None): return a
        def Allreduce(a, b, op=None): return a
        def seq(F, *args): F(*args)
        print("Warning: Converter:Mpi: mpi4py is not available. Sequential behaviour.")

from .Distributed import _readZones, _convert2PartialTree, _convert2SkeletonTree, _readPyTreeFromPaths, mergeGraph, splitGraph, isZoneSkeleton__
from . import PyTree as C
from . import Internal

# Previous times for CPU time measures
PREVFULLTIME = None # full

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
    #print('Rank %d has %d zones.'%(rank, len(zones)))
    tl = C.center2Node(tl, var, cellNType)
    tl = rmXZones(tl)
    return tl

# addGhostCells parallel
# si modified=[], transferts tous les champs
# si modified=None, transfert aucun champ
def _addGhostCells(t, b, d, adaptBCs=1, modified=[], fillCorner=1):
    """Add ghost cells to pyTree."""
    if modified == []: # all
        variables = C.getVarNames(t, excludeXYZ=True, loc='nodes')[0]
        variables += C.getVarNames(t, excludeXYZ=True, loc='centers')[0]
    elif modified is None: variables = []; modified = []
    else: variables = modified

    _addMXZones(t, depth=2, variables=variables, noCoordinates=False, 
                keepOldNodes=False)
    #print("%d: addGC(max): Nblocs=%d, NPts(M)=%g"%(rank,len(Internal.getZones(t)), C.getNPts(t)*1./1.e6), flush=True)
    Internal._addGhostCells(t, t, d, adaptBCs, modified, fillCorner)
    _rmMXZones(t)

    # ancienne version utilisant addXZones
    #graph = computeGraph(t, type='match', reduction=True)
    #_addXZones(t, graph, variables=variable, noCoordinates=False, 
    #           zoneGC=False, keepOldNodes=False)
    #print("%d: addGC(max): Nblocs=%d, NPts(M)=%g"%(rank,len(Internal.getZones(t)), C.getNPts(t)*1./1.e6), flush=True)
    #Internal._addGhostCells(t, t, d, adaptBCs, modified, fillCorner)
    #_rmXZones(t)
    return None

def getNPts(a):
    """Return the number of points in a."""
    npts = C.getNPts(a)
    npts = allreduce(npts, op=SUM)
    return npts

def getNCells(a):
    """Return the number of cells in t."""
    ncells = C.getNPts(a)
    ncells = allreduce(ncells, op=SUM)
    return ncells

def getMinValue(t, varName):
    """Get the minimum value of variable defined by var."""
    val = C.getMinValue(t, varName)
    val = allreduce(val, op=MIN)
    return val

def getMaxValue(t, varName):
    """Get the maximum value of variable defined by var."""
    val = C.getMaxValue(t, varName)
    val = allreduce(val, op=MAX)
    return val

def getMeanValue(t, varName):
    """Get the mean value of variable defined by var."""
    val = C.getMeanValue(t, varName)
    npts = C.getNPts(t)
    val = val*npts
    val = allreduce(val, op=SUM)
    npts = allreduce(npts, op=SUM)
    return val/npts

# Ecrit une trace dans un fichier proc0.out
# si cpu=True, ecrit le temps depuis le trace precedent
# si mem=True, ecrit l'etat de la mem du noeud
# si stdout=True, ecrit a l'ecran
def trace(text=">>> IN XXX: ", cpu=True, mem=True, stdout=False):
    """Write a trace of cpu and memory in a file or to stdout for current node."""
    global PREVFULLTIME
    msg = text
    if cpu:
        if PREVFULLTIME is None:
            dt = 0.
            PREVFULLTIME = timeit.default_timer()
        else:
            t = timeit.default_timer()
            dt = t - PREVFULLTIME
            PREVFULLTIME = t
        msg += ' [%g secs]'%dt
    if mem:
        pid = os.getpid()
        try: 
            f = open("/proc/{}/smaps".format(pid))
            s = f.readlines()
            f.close()
        except: s = []
        
        tot = 0.
        found = False
        for ts in s:
            if found:
                tot += int(ts[5:-3])
                found = False
            if ts.find("heap") >= 0: found = True
        if tot > 1.e6: msg += '[{} GB]'.format(tot/1.e6)
        elif tot > 1000.: msg += '[{} MB]'.format(tot/1000.)
        else: msg += '[{} kB]'.format(tot)
    msg += '\n'

    if stdout: # ecriture a l'ecran
        print('%d: %s'%(rank, msg)) 
        sys.stdout.flush()
    else: # dans des fichiers
        f = open('proc%03d.out'%rank, "a")
        f.write(msg)
        f.flush()
        f.close()
    return None
