# Interface pour MPI

import os, time, timeit, sys
from .Distributed import _readZones, _convert2PartialTree, _convert2SkeletonTree, _readPyTreeFromPaths, mergeGraph, splitGraph, isZoneSkeleton__
from . import PyTree as C
from . import Internal
from . import Distributed

if 'MPIRUN' in os.environ: # si MPIRUN=0, force sequentiel
    if int(os.environ['MPIRUN'])>0:
        try: from .Mpi4py import *
        except: raise ImportError("Converter:Mpi: requires mpi4py module.")
    else:
        rank = 0; size = 1; KCOMM = None; COMM_WORLD = None
        SUM = 0; MAX = 0; MIN = 0; LAND = 0
        from .Distributed import setProc, _setProc, getProc, getProcDict, getProperty, getPropertyDict, convertFile2SkeletonTree, computeGraph, splitGraph, mergeGraph, readZones, writeZones, convert2PartialTree, convert2SkeletonTree, readPyTreeFromPaths
        def barrier(): return
        def bcast(a, root=0): return a
        def Bcast(a, root=0): return a
        def gather(a, root=0): return a
        def Gather(a, root=0): return a
        def gatherZones(a, root=0): return a
        def allgather(a): return [a]
        def allgatherZones(a, root=0): return a
        def allgatherTree(a): return a
        def allgatherDict(a): return a
        def send(a, dest=0, tag=0): return None
        def recv(source=0, tag=0): return None # pb here
        def sendRecv(a, source=0, dest=0): return []
        def sendRecvC(a, source=0, dest=0): return []
        def reduce(a, op=None, root=0): return a
        def Reduce(a, b, op=None, root=0): return a
        def allreduce(a, op=None): return a
        def Allreduce(a, b, op=None): b[:] = a[:]; return None
        def seq(F, *args): F(*args)
        def convertFile2PyTree(fileName, format=None, proc=None): return C.convertFile2PyTree(fileName, format)
        def convertPyTree2File(t, fileName, format=None, links=[], ignoreProcNodes=False, merge=True): return C.convertPyTree2File(t, fileName, format, links)
        def addXZones(t, graph, variables=None, noCoordinates=False, cartesian=False, subr=True, keepOldNodes=True, zoneGC=True): return Internal.copyRef(t)
        def _addXZones(t, graph, variables=None, noCoordinates=False, cartesian=False, subr=True, keepOldNodes=True, zoneGC=True): return None
        def _addLXZones(t, graph, variables=None, cartesian=False, interDict=[], bboxDict={}, layers=2, subr=True): return None
        def _addBXZones(a, depth=2, allB=False): return None
        def rmXZones(t): return Internal.copyRef(t)
        def _rmXZones(t): return None
        def _rmBXZones(t): return None
        def createBboxDict(t):
            import Generator.PyTree as G
            bboxDict = {}
            for z in Internal.getZones(t): bboxDict[z[0]] = G.bbox(z)
            return bboxDict
        #print("Warning: Converter:Mpi: Sequential behaviour is forced by MPIRUN=0.")

else: # try import (may fail - core or hang)
    try: from .Mpi4py import *
    except:
        rank = 0; size = 1; KCOMM = None; COMM_WORLD = None
        SUM = 0; MAX = 0; MIN = 0; LAND = 0
        from .Distributed import setProc, _setProc, getProc, getProcDict, getProperty, getPropertyDict, convertFile2SkeletonTree, computeGraph, splitGraph, mergeGraph, readZones, convert2PartialTree, convert2SkeletonTree, readPyTreeFromPaths
        def barrier(): return
        def bcast(a, root=0): return a
        def Bcast(a, root=0): return a
        def gather(a, root=0): return a
        def Gather(a, root=0): return a
        def gatherZones(a, root=0): return a
        def allgather(a): return [a]
        def allgatherZones(a, root=0): return a
        def allgatherTree(a): return a
        def allgatherDict(a): return a
        def send(a, dest=0, tag=0): return None
        def recv(source=0, tag=0): return None # pb here
        def sendRecv(a, source=0, dest=0): return []
        def sendRecvC(a, source=0, dest=0): return []
        def reduce(a, op=None, root=0): return a
        def Reduce(a, b, op=None, root=0): return a
        def allreduce(a, op=None): return a
        def Allreduce(a, b, op=None): b[:] = a[:]; return None
        def seq(F, *args): F(*args)
        def convertFile2PyTree(fileName, format=None, proc=None): return C.convertFile2PyTree(fileName, format)
        def convertPyTree2File(t, fileName, format=None, links=[], ignoreProcNodes=False, merge=True): return C.convertPyTree2File(t, fileName, format, links)
        def addXZones(t, graph, variables=None, noCoordinates=False, cartesian=False, subr=True, keepOldNodes=True, zoneGC=True): return Internal.copyRef(t)
        def _addXZones(t, graph, variables=None, noCoordinates=False, cartesian=False, subr=True, keepOldNodes=True, zoneGC=True): return None
        def _addLXZones(t, graph, variables=None, cartesian=False, interDict=[], bboxDict={}, layers=2, subr=True): return None
        def _addBXZones(a, depth=2, allB=False): return None
        def rmXZones(t): return Internal.copyRef(t)
        def _rmXZones(t): return None
        def _rmBXZones(t): return None
        def createBboxDict(t):
            import Generator.PyTree as G
            bboxDict = {}
            for z in Internal.getZones(t): bboxDict[z[0]] = G.bbox(z)
            return bboxDict
        print("Warning: Converter:Mpi: mpi4py is not available. Sequential behaviour.")

# Previous times for CPU time measures (trace)
PREVFULLTIME = None # full

#==============================================================================
# IN: t: full/loaded skel/partial
#==============================================================================
def center2Node(t, var=None, cellNType=0, graph=None):
    allstructured = 1
    for z in Internal.getZones(t):
        type = Internal.getZoneType(z)
        if type != 1: allstructured = 0; break
    allstructured = allreduce(allstructured)
    if allstructured == size: # all zones are structured
        return center2Node1__(t, var, cellNType) # to be changed to 2
    else: # mixed or unstructured
        return center2Node1__(t, var, cellNType, graph)

# generic, use addXZones
def center2Node1__(t, var=None, cellNType=0, graph=None):
    """Convert a zone or a field from centers to node."""
    if graph is None: graph = computeGraph(t, type='match')
    tl = addXZones(t, graph)
    _convert2PartialTree(tl)
    #zones = Internal.getZones(tl)
    #print('Rank %d has %d zones.'%(rank, len(zones)))
    tl = C.center2Node(tl, var, cellNType)
    _rmXZones(tl)
    return tl

# only for all structured grids (use addMXZones)
def center2Node2__(t, var=None, cellNType=0):
    """Convert a zone or a field from centers to node (structured grids only)."""
    tl = Internal.copyRef(t)
    if var is None: noCoordinates = False
    else: noCoordinates = True
    _addMXZones(tl, variables=var, noCoordinates=noCoordinates, keepOldNodes=False)
    _convert2PartialTree(tl)
    #zones = Internal.getZones(tl)
    #print('Rank %d has %d zones.'%(rank, len(zones)))
    t2 = C.center2Node(tl, var, cellNType)
    _rmMXZones(tl)
    _rmMXZones(t2)
    return t2

def addGhostCells(t, b, d, adaptBCs=1, modified=[], fillCorner=1):
    """Add ghost cells to pyTree."""
    tp = Internal.copyRef(t)
    _addGhostCells(tp, b, d, adaptBCs=adaptBCs, modified=modified, fillCorner=fillCorner)
    return tp

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

def isFinite(t, var=None):
    """Return true if all fields in a have no NAN or INF values."""
    val = C.isFinite(t, var)
    val = allreduce(val, op=SUM)
    if val == size: return True
    else: return False

# Ecrit une trace dans un fichier proc0.out
# si cpu=True, ecrit le temps depuis le trace precedent
# si mem=True, ecrit l'etat de la memoire du noeud
# si stdout=True, ecrit a l'ecran
# si reset, vide le fichier log
# si filename="stdout", ecrit a l'ecran, sinon ecrit dans le fichier filename
def trace(text=">>> IN XXX: ", cpu=True, mem=True, stdout=False, reset=False, fileName="stdout"):
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
            f = open("/proc/%s/smaps"%(pid))
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
        if tot > 1.e6: msg += '[%f GB]'%(tot/1.e6)
        elif tot > 1000.: msg += '[%f MB]'%(tot/1000.)
        else: msg += '[%f kB]'%(tot)
    msg += '\n'

    if stdout:
        #print("Warning: trace: arg. stdout is obsolete. Please use fileName='stdout' instead.")
        fileName = "stdout"

    if fileName == "stdout":
        print('%d: %s'%(rank, msg))
        sys.stdout.flush()
    else: # dans des fichiers par processes
        fileName = fileName.split('.')
        if '%' in fileName[0]: fileName[0] = fileName[0]%rank
        else: fileName[0] += '%03d'%rank
        if len(fileName) == 1: fileName[0] += '.out'
        fileName = '.'.join(s for s in fileName)
        if reset: f = open(fileName, "w")
        else: f = open(fileName, "a")
        f.write(msg)
        f.flush()
        f.close()
    return None

#==============================================================================
# Construit un arbre de BBox a partir d'un arbre squelette charge
# ou d'un arbre partiel
# L'arbre des BBox final est identique sur tous les processeurs
# Nota Bene: Si method='OBB' et weighting=1, les zones de t doivent etre
# formees par des triangles. Pour cela, appliquer au prealable (par exemple):
# createBBoxTree(C.convertArray2Tetra(P.exteriorFaces(t)),isOBB=1,weighting=1)
#==============================================================================
def createBBoxTree(t, method='AABB', weighting=0, tol=0.):
    """Return a bbox tree of t."""
    try: import Generator.PyTree as G
    except: raise ImportError("createBBoxTree requires Generator module.")
    # bounding box des zones locales
    bases = Internal.getBases(t)
    zb = []
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            if not Distributed.isZoneSkeleton__(z):
                zbb = G.BB(z, method, weighting, tol=tol)
                # ajoute baseName/zoneName
                zbb[0] = b[0]+'/'+zbb[0]
                # Clean up (zoneSubRegion)
                Internal._rmNodesFromType(zbb, 'ZoneSubRegion_t')
                zb.append(zbb)

    # Echanges des zones locales de bounding box
    # (allgather serialise automatiquement les donnees)
    zones = allgather(zb)

    # On les remplace dans l'arbre (par noms)
    tb = Distributed.setZonesInTree(t, zones)
    return tb
