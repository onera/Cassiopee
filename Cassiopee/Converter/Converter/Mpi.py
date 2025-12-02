# Interface pour MPI

import os, timeit, sys
from .Distributed import (
    _readZones, _convert2PartialTree, _convert2SkeletonTree,
    _readPyTreeFromPaths, mergeGraph, splitGraph, isZoneSkeleton__
)
from . import PyTree as C
from . import Internal
from . import Distributed

if 'MPIRUN' in os.environ: # si MPIRUN=0, force sequentiel
    if int(os.environ['MPIRUN']) > 0:
        try: from .Mpi4py import *
        except: raise ImportError("Converter:Mpi: requires mpi4py module.")
    else:
        rank = 0; size = 1; KCOMM = None; COMM_WORLD = None
        master = True
        SUM = 0; MAX = 0; MIN = 0; LAND = 0
        from .Distributed import (
            setProc, _setProc, getProc, getProcDict,
            getProperty, getPropertyDict, convertFile2SkeletonTree,
            computeGraph, splitGraph, mergeGraph, readZones, writeZones,
            convert2PartialTree, convert2SkeletonTree, readPyTreeFromPaths
        )
        def abort(errorcode=0): os._exit(errorcode)
        def barrier(): return
        def bcast(a, root=0): return a
        def Bcast(a, root=0): return a
        def bcastZone(a, root=0, coord=True, variables=[]): return a
        def gather(a, root=0): return a
        def Gather(a, root=0): return a
        def gatherZones(a, root=0): return a
        def allgather(a): return [a]
        def allgatherZones(a, coord=True, variables=[]): return a
        def allgatherTree(a): return a
        def allgatherDict(a): return a
        def allgatherNext(a): return [a]
        def send(a, dest=0, tag=0): return None
        def isend(a, dest=0, tag=0): return None
        def recv(source=0, tag=0): return None # pb here
        def requestWaitall(reqs): return None
        def sendRecv(a, source=0, dest=0): return []
        def sendRecvC(a, source=0, dest=0): return []
        def reduce(a, op=None, root=0): return a
        def Reduce(a, b, op=None, root=0): return a
        def allreduce(a, op=None): return a
        def Allreduce(a, b, op=None): b[:] = a[:]; return None
        def getSizeOf(a): return Internal.getSizeOf(a)
        def passNext(a): return a
        def seq(F, *args): F(*args)
        def convertFile2PyTree(fileName, format=None, proc=None):
            return C.convertFile2PyTree(fileName, format)
        def convertPyTree2File(t, fileName, format=None, links=[], isize=8, rsize=8, ignoreProcNodes=False, merge=True): return C.convertPyTree2File(t, fileName, format=format, links=links, isize=isize, rsize=rsize)
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
        master = True
        SUM = 0; MAX = 0; MIN = 0; LAND = 0
        from .Distributed import (
            setProc, _setProc, getProc, getProcDict,
            getProperty, getPropertyDict, convertFile2SkeletonTree,
            computeGraph, splitGraph, mergeGraph, readZones, writeZones,
            convert2PartialTree, convert2SkeletonTree, readPyTreeFromPaths
        )
        def abort(errorcode=0): os._exit(errorcode)
        def barrier(): return
        def bcast(a, root=0): return a
        def Bcast(a, root=0): return a
        def bcastZone(a, root=0, coord=True, variables=[]): return a
        def gather(a, root=0): return a
        def Gather(a, root=0): return a
        def gatherZones(a, root=0): return a
        def allgather(a): return [a]
        def allgatherZones(a, coord=True, variables=[]): return a
        def allgatherTree(a): return a
        def allgatherDict(a): return a
        def allgatherNext(a): return [a]
        def send(a, dest=0, tag=0): return None
        def isend(a, dest=0, tag=0): return None
        def recv(source=0, tag=0): return None # pb here
        def requestWaitall(reqs): return None
        def sendRecv(a, source=0, dest=0): return []
        def sendRecvC(a, source=0, dest=0): return []
        def reduce(a, op=None, root=0): return a
        def Reduce(a, b, op=None, root=0): return a
        def allreduce(a, op=None): return a
        def Allreduce(a, b, op=None): b[:] = a[:]; return None
        def getSizeOf(a): return Internal.getSizeOf(a)
        def passNext(a): return a
        def seq(F, *args): F(*args)
        def convertFile2PyTree(fileName, format=None, proc=None):
            return C.convertFile2PyTree(fileName, format)
        def convertPyTree2File(t, fileName, format=None, links=[], isize=8, rsize=8, ignoreProcNodes=False, merge=True): return C.convertPyTree2File(t, fileName, format=format, links=links, isize=isize, rsize=rsize)
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

# Store state for trace
TRACESTATE = { 'prevFullTime': None, 'method': 0, 'fileName': 'stdout', 'mem': True, 'cpu': True, 'master': True }

#==============================================================================
# IN: t: full/loaded skel/partial
#==============================================================================
def center2Node(t, var=None, cellNType=0, graph=None):
    """Convert zone/fields defined at centers to nodes."""
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
    return None

def getNPts(a):
    """Return the number of points in a."""
    npts = C.getNPts(a)
    npts = allreduce(npts, op=SUM)
    return npts

def getNCells(a):
    """Return the number of cells in t."""
    ncells = C.getNCells(a)
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

# Trace memory and cpu time.
# if cpu=True, write cpu time ellapsed from previous trace call.
# if mem=True, write node used memory
# if reset is true, empty log
# if filename="stdout", write to screen, else write in file filename
# if master=True, only master write output and max is performed on mem and time.
# if method is present: Rss(0), tracemalloc(1), heap(2), VmRss(3), Resident(4)
# note: 1 counts only memory allocated by python and not the one allocated in pure C.
# note: trace default arguments are defined in TRACESTATE (line 96)
def trace(text=">>> IN XXX: ", cpu=None, mem=None, reset=False, fileName=None, method=None, master=None):
    """Write a trace of cpu and memory in a file or to stdout for current node."""
    global TRACESTATE
    msg = text
    if method is not None: TRACESTATE['method'] = method
    if fileName is not None:
        TRACESTATE['fileName'] = fileName
    if mem is not None: TRACESTATE['mem'] = mem
    if cpu is not None: TRACESTATE['cpu'] = cpu
    if master is not None: TRACESTATE['master'] = master

    dt = 0.; tot = 0.; peak = -1

    if TRACESTATE['mem'] and TRACESTATE['method'] == 1:
        import tracemalloc
        if TRACESTATE['prevFullTime'] is None: tracemalloc.start()
    if TRACESTATE['cpu']:
        if TRACESTATE['prevFullTime'] is None:
            dt = 0.
            TRACESTATE['prevFullTime'] = timeit.default_timer()
        else:
            t = timeit.default_timer()
            dt = t - TRACESTATE['prevFullTime']
            TRACESTATE['prevFullTime'] = t
        hours = int(dt // 3600)
        minutes = int((dt % 3600) // 60)
        seconds = (dt % 60)
        msg += ' [%g hr %g min %g sec]'%(hours, minutes, seconds)
    if TRACESTATE['mem']:
        tot = 0.; peak = -1
        if TRACESTATE['method'] == 0: # Rss in smaps
            pid = os.getpid()
            try:
                f = open("/proc/%s/smaps"%pid)
                s = f.readlines()
                f.close()
            except: s = []
            for ts in s:
                if ts.split()[0] == "Rss:":
                    tot += int(ts[4:-3])
        elif TRACESTATE['method'] == 1: # tracemalloc
            tot, peak = tracemalloc.get_traced_memory()
            tot = tot//1000 # in kb
            peak = peak//1000 # in kb
        elif TRACESTATE['method'] == 2: # heap in smaps
            pid = os.getpid()
            try:
                f = open("/proc/%s/smaps"%pid)
                s = f.readlines()
                f.close()
            except: s = []
            found = False
            for ts in s:
                if found:
                    tot += int(ts[5:-3])
                    found = False
                if ts.find("heap") >= 0: found = True
        elif TRACESTATE['method'] == 3: # VmRSS in status
            pid = os.getpid()
            try:
                f = open("/proc/%s/status"%pid)
                s = f.readlines()
                f.close()
            except: s = []
            for ts in s:
                if ts.split()[0] == "VmRSS:":
                    tot += int(ts[6:-3])
        elif TRACESTATE['method'] == 4: # Resident in statm
            pid = os.getpid()
            try:
                f = open("/proc/%s/statm"%pid)
                s = f.readlines()
                f.close()
            except: s = []
            for ts in s:
                tot += int(ts.split()[1])

    if TRACESTATE['master']: # reduction
        a = allgather((dt, tot))
        dt = 0.; tot = 0.
        for i in a:
            dt = max(dt, i[0]); tot = max(tot, i[1])

    # build msg
    if TRACESTATE['cpu']: msg += ' [%g secs]'%dt
    if TRACESTATE['mem']:
        if peak == -1: # peak is not known
            if tot > 1.e6: msg += '[%f GB]'%(tot/1.e6)
            elif tot > 1000.: msg += '[%f MB]'%(tot/1000.)
            else: msg += '[%f kB]'%(tot)
        else:
            if tot > 1.e6: msg += '[current: %f GB | peak: %f GB]'%(tot/1.e6, peak/1.e6)
            elif tot > 1000.: msg += '[current: %f MB | peak: %f MB]'%(tot/1000., peak/1000.)
            else: msg += '[current: %f kB | peak: %f kB]'%(tot, peak)
    msg += '\n'

    if TRACESTATE['fileName'] == "stdout":
        if TRACESTATE['master']:
            if rank == 0:
                print('%d: %s'%(rank, msg))
                sys.stdout.flush()
        else:
            print('%d: %s'%(rank, msg))
            sys.stdout.flush()

    else: # dans des fichiers par processes si pas master ou dans un seul si master
        if not TRACESTATE['master']:
            fileName = TRACESTATE['fileName'].split('.')
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
def createBBoxTree(t, method='AABB', weighting=0, tol=0., keepOldNodes=True):
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
                if keepOldNodes:
                    # Clean up (zoneSubRegion)
                    Internal._rmNodesFromType(zbb, 'ZoneSubRegion_t')
                else:
                    C._extractVars(zbb, None, keepOldNodes=False)
                    _setProc(zbb, getProc(z)) # keep proc
                zb.append(zbb)

    # Echanges des zones locales de bounding box
    # (allgather serialise automatiquement les donnees)
    zones = allgather(zb)

    # On les remplace dans l'arbre (par noms)
    tb = Distributed.setZonesInTree(t, zones)
    return tb
