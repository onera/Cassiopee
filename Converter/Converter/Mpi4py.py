# Interface de Cassiopee pour mpi
from . import PyTree as C
from . import Internal
from . import Distributed

# Acces a Distributed
from .Distributed import readZones, _readZones, convert2PartialTree, _convert2PartialTree, convert2SkeletonTree, readNodesFromPaths, readPyTreeFromPaths, writeNodesFromPaths

__all__ = ['rank', 'size', 'KCOMM', 'COMM_WORLD', 'setCommunicator', 'barrier', 'send', 'recv', 'sendRecv',
    'bcast', 'Bcast', 'bcastZone', 'allgatherZones',
    'allgather', 'readZones', 'writeZones', 'convert2PartialTree', 'convert2SkeletonTree', 'convertFile2DistributedPyTree', 
    'readNodesFromPaths', 'readPyTreeFromPaths', 'writeNodesFromPaths',
    'allgatherTree', 'convertFile2SkeletonTree', 'convertFile2PyTree', 'convertPyTree2File', 'seq', 'print0', 'printA',
    'createBBoxTree', 'createBboxDict', 'computeGraph', 'addXZones', '_addXZones', '_addMXZones', '_addBXZones', '_addLXZones',
    'rmXZones', '_rmXZones', '_rmMXZones', '_rmBXZones', 'getProcDict', 'getProc', 'setProc', '_setProc', 'COMM_WORLD']

from mpi4py import MPI
import numpy

try: range = xrange
except: pass

COMM_WORLD = MPI.COMM_WORLD
KCOMM = COMM_WORLD

rank = KCOMM.rank
size = KCOMM.size

# version collective
def writeZones(t, fileName, format=None, proc=None, zoneNames=None, links=None):
    """Write zones in parallel."""
    seq(Distributed.writeZones, t, fileName, format, proc, zoneNames, links)
    return None
        
#==============================================================================
# Change de communicateur
#==============================================================================
def setCommunicator(com):
    """Set MPI communicator to com."""
    global KCOMM, rank, size
    KCOMM = com
    rank = KCOMM.rank
    size = KCOMM.size

#==============================================================================
# Barrier
#==============================================================================
def barrier():
    KCOMM.barrier()

#==============================================================================
# Send
#==============================================================================
def send(obj, dest=None):
    KCOMM.send(obj, dest)

#==============================================================================
# Receive
#==============================================================================
def recv(source=None):
    return KCOMM.recv(source)

#==============================================================================
# Send and receive with a graph
# IN: datas: un dictionnaire des donnees a envoyer par proc de destination
# OUT: un dictionnaire des donnees recues par proc d'origine
#==============================================================================
def sendRecv(datas, graph):
    if graph == {}: return {}
    reqs = []
    
    if rank in graph:
        g = graph[rank] # graph du proc courant
        for oppNode in g:
            # Envoie les datas necessaires au noeud oppose
            #print('%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode]))
            s = KCOMM.isend(datas[oppNode], dest=oppNode)
            reqs.append(s)
    rcvDatas={}
    for node in graph:
        #print(rank, graph[node].keys())
        if rank in graph[node]:
            #print('%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank]))
            rcvDatas[node] = KCOMM.recv(source=node)
    MPI.Request.Waitall(reqs)
    return rcvDatas

#==============================================================================
# allgather
#==============================================================================
def allgather(data):
    ret = KCOMM.allgather(data)
    # Si dictionnaire de listes, on fusionne les listes
    # Si dictionnaire d'autre chose, on append dans des listes
    if isinstance(data, dict):
        out = {}
        for r in ret:
            for k in r:
                if k not in out: out[k] = r[k]
                else:
                    try: out[k] += r[k]
                    except: out[k].append(r[k])
        return out
    else: return ret

def allgatherTree(t):
    """Gather a distributed tree on all processors."""
    d = KCOMM.allgather(t)
    return Internal.merge(d)
    
#==============================================================================
# bcast from root
#==============================================================================
# data=All with pickle
def bcast(data, root=0):
    return KCOMM.bcast(data, root)

# data=numpy
def Bcast(data, root=0):
    return KCOMM.Bcast(data, root)

# data=tree
def bcastTree(t, root=0):
    if t is not None:
        zones = Internal.getZones(t)
        for z in zones: z = bcastZone(z)
    return t

# data=zone
# Envoie uniquement les coords de z
def bcastZone(z, root=0):
    # zp = squelette envoye en pickle
    if rank == root:
        zp = Internal.copyRef(z)
        Internal._rmNodesFromType(zp, 'GridCoordinates_t')
        Internal._rmNodesFromType(zp, 'FlowSolution_t')
        Internal._rmNodesFromType(zp, 'ZoneSubRegion_t')
    else: zp = None

    if rank == root:
        px = Internal.getNodeFromName2(z, 'CoordinateX')[1]
        py = Internal.getNodeFromName2(z, 'CoordinateY')[1]
        pz = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
        sx = px.shape; sy = py.shape; sz = pz.shape
    else:
        sx = None; sy = None; sz = None

    zp, sx, sy, sz = KCOMM.bcast((zp, sx, sy, sz), root)
        
    if rank != root:
        px = numpy.empty(sx, dtype=numpy.float64, order='F')
        py = numpy.empty(sy, dtype=numpy.float64, order='F')
        pz = numpy.empty(sz, dtype=numpy.float64, order='F')
    
    KCOMM.Bcast([px,MPI.DOUBLE], root)
    KCOMM.Bcast([py,MPI.DOUBLE], root)
    KCOMM.Bcast([pz,MPI.DOUBLE], root)
    
    if rank != root:
        # Reconstruction de la zone
        Internal._createUniqueChild(zp, Internal.__GridCoordinates__, 'GridCoordinates_t')
        n = Internal.getNodeFromName1(zp, Internal.__GridCoordinates__)
        Internal._createUniqueChild(n, 'CoordinateX', 'DataArray_t', value=px)
        Internal._createUniqueChild(n, 'CoordinateY', 'DataArray_t', value=py)
        Internal._createUniqueChild(n, 'CoordinateZ', 'DataArray_t', value=pz)
    else: zp = z
    return zp

# All gather une liste de zones, recuperation identique sur tous les procs
# dans une liste a plat
# Uniquement les coordonnees sont envoyees
def allgatherZones(zones):
    # Chaque processeur bcast ses zones vers les autres ranks
    zones = Internal.getZones(zones)
    lenZones = KCOMM.allgather(len(zones))
    allZones = []
    for i in range(size):
        for cz in range(lenZones[i]):
            if rank == i: zp = bcastZone(zones[cz], root=i)
            else: zp = bcastZone(None, root=i)
            allZones.append(zp)
    return allZones

#==============================================================================
# Lecture du squelette d'un arbre dans un fichier
# Lecture proc 0 + bcast
#==============================================================================
def convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5,
                             maxDepth=-1, links=None):
    """Read a file and return a skeleton tree."""
    if rank == 0: t = Distributed.convertFile2SkeletonTree(fileName, format, maxFloatSize, maxDepth, None, links)
    else: t = None
    t = KCOMM.bcast(t)
    if links is not None:
        lk = KCOMM.bcast(links)
        if rank > 0: links += lk
    return t
    
#==============================================================================
# Only for hdf
# Split a file on all processors
#==============================================================================
def convertFile2DistributedPyTree(fileName):
  import etc.toolbox.internal as tgi
  tmp = tgi.convertFile2DistributedPyTree(fileName)
  return tmp.getLocalTree()

#==============================================================================
# Lecture complete d'un arbre dans un fichier
# si proc=None, lecture proc 0 + bcast
# sinon lecture des zones correspondant a proc
#==============================================================================
def convertFile2PyTree(fileName, format=None, proc=None):
    """Read a file and return a full tree or partial tree."""
    if proc is None: # load full tree on all procs
        if rank == 0: t = C.convertFile2PyTree(fileName, format)
        else: t = None
        t = KCOMM.bcast(t)
    else:
        t = convertFile2SkeletonTree(fileName, format)
        _readZones(t, fileName, rank=proc)
        _convert2PartialTree(t, rank=proc)
    return t

# parallel merge on proc0 without zones
def _merge__(t):
    tp = Internal.copyRef(t)
    if not Internal.isTopTree(t): return tp
    if rank > 0:
        out = []
        for i in t[2]:
            if i[3] != 'Zone_t': out.append(i)
        tp[2] = out
        KCOMM.send(tp, dest=0)
    if rank == 0:
        for i in range(1, size):
            ret = KCOMM.recv(source=i)
            if ret is not None:
                tp = Internal.merge([tp, ret])
    return tp
    
#==============================================================================
# Ecriture sequentielle
# Avec recuperation de toutes les zones
# si ignoreProcNodes=True, toutes les zones non squelette sont ecrites,
# sinon ecrit seulement les zones procNode correspondant a rank
# si merge=True, effectue un merge prealable des arbres de tous les procs 
# sans zones (pour recuperer toutes les bases et les data des bases)
#==============================================================================
def convertPyTree2File(t, fileName, format=None, links=[], 
    ignoreProcNodes=False, merge=True):
    """Write a skeleton or partial tree."""
    tp = convert2PartialTree(t)
    tp = C.deleteEmptyZones(tp)
    Internal._adaptZoneNamesForSlash(tp)
    if merge: tp = _merge__(tp)
        
    nzones = len(Internal.getZones(tp))
    if rank == 0:
        if nzones > 0:
            C.convertPyTree2File(tp, fileName, format=format, links=links); go = 1
        else: go = 0
        if size > 1: KCOMM.send(go, dest=1)
    else:
        go = KCOMM.recv(source=rank-1)
        if go == 1:
            if ignoreProcNodes: Distributed.writeZones(tp, fileName, format=format, proc=-1, links=links)
            else: Distributed.writeZones(tp, fileName, format=format, proc=rank, links=links)
        else:
            if nzones > 0:
                C.convertPyTree2File(tp, fileName, format=format, links=links); go = 1
        if rank < size-1: KCOMM.send(go, dest=rank+1)
    barrier()
    
#==============================================================================
# Execute sequentiellement F sur tous les procs
#==============================================================================
def seq(F, *args):
    if rank == 0:
        F(*args)
        if size > 1: KCOMM.send(1, dest=1)
    else:
        go = KCOMM.recv(source=rank-1)
        F(*args)
        if rank < size-1: KCOMM.send(rank+1, dest=rank+1)
    barrier()
        
#==============================================================================
# Print uniquement du proc 0
#==============================================================================
def print0(a):
    if rank == 0: print(a)

#==============================================================================
# Print sur tous les procs sequentiellement
#==============================================================================
def printA(A):
    def fprint(A): print(A)
    seq(fprint, A)

#==============================================================================
# Construit un arbre de BBox a partir d'un arbre squelette charge
# ou d'un arbre partiel
# L'arbre des BBox final est identique sur tous les processeurs
# Nota Bene: Si method='OBB' et weighting=1, les zones de t doivent etre 
# formees par des triangles. Pour cela, appliquer au prealable (par exemple):
# createBBoxTree(C.convertArray2Tetra(P.exteriorFaces(t)),isOBB=1,weighting=1)
#==============================================================================
def createBBoxTree(t, method='AABB', weighting=0):
    """Return a bbox tree of t."""
    try: import Generator.PyTree as G
    except: raise ImportError("createBBoxTree requires Generator module.")
    # bounding box des zones locales
    tp = Internal.node2PyTree(t)
    bases = Internal.getBases(t)
    zb = []
    for b in bases:
        zones = Internal.getNodesFromType1(b, 'Zone_t')
        for z in zones:
            if not Distributed.isZoneSkeleton__(z):
                zbb = G.BB(z, method, weighting)
                # ajoute baseName/zoneName
                zbb[0] = b[0]+'/'+zbb[0]
                # Clean up (zoneSubRegion)
                Internal._rmNodesFromType(zbb, 'ZoneSubRegion_t')
                zb.append(zbb)
        
    # Echanges des zones locales de bounding box
    # (allgather serialise automatiquement les donnees)
    zones = KCOMM.allgather(zb)

    # On les remplace dans l'arbre (par noms)
    tb = Distributed.setZonesInTree(t, zones)
    return tb

#==============================================================================
# Calcule le dictionnaire des bbox de l'arbre complet
# Utile pour addXZones optimises
#==============================================================================
def createBboxDict(t):
    try: import Generator.PyTree as G
    except: raise ImportError("createBBoxTree requires Generator module.")
    bboxDict = {}
    zones = Internal.getZones(t)
    for z in zones:
        bboxDict[z[0]] = G.bbox(z)

    b = KCOMM.allgather(bboxDict)
    bboxDict = {}
    for i in range(len(b)):
        for j in b[i]:
            if not j in bboxDict: bboxDict[j] = b[i][j]
    return bboxDict

#==============================================================================
# Calcule le graph
# graph[proc1][proc2] est la liste des zones de proc1 intersectant au moins 
# une zone de proc2
# IN: type: type de graph
# IN: reduction: si True, on assure que le graph est le meme sur tous
# les procs. A utiliser si le graph est construit a partir d'arbres partiels
# ou pour type='proc'
# IN: procDict: le procDict de tous l'arbre, si None, on le calcul.
# IN: intersectionsDict: dictionnaire d'intersections. Consulter la doc de
# computeGraph dans Distributed.py pour plus de details.
#==============================================================================
def computeGraph(t, type='bbox', t2=None, procDict=None, reduction=True, 
                 intersectionsDict=None, exploc=False, procDict2=None, it=0):
    """Return the communication graph for different block relation types."""
    if not procDict: procDict = getProcDict(t)
    graph = Distributed.computeGraph(t, type, t2, procDict, rank, 
                                     intersectionsDict, exploc, procDict2, it)

    if reduction:
        # Assure que le graph est le meme pour tous les processeurs
        g = KCOMM.allgather(graph)
        graph = {}
        for i in g:
            for k in i:
                if not k in graph: graph[k] = {}
                for j in i[k]:
                    if not j in graph[k]: graph[k][j] = []
                    graph[k][j] += i[k][j]
                    graph[k][j] = list(set(graph[k][j]))

    return graph

#=============================================================================
# Calcule l'intersection de deux bbox
#=============================================================================
def GetIntersectionBbox(bbox1, bbox2):
   Ibbox = numpy.zeros(6)
   Ibbox[0:3] = [max(bbox1[i],bbox2[i]) for i in range(3)]
   Ibbox[3:6] = [min(bbox1[i],bbox2[i]) for i in range(3,6)]
   return Ibbox

#==============================================================================
# Recupere les zones specifiees dans le graph, les ajoute a l'arbre local t
# if subr=True, the ZoneSubRegions are kept during the exchange 
#==============================================================================
def addXZones(t, graph, variables=None, noCoordinates=False, cartesian=False, subr=True):
    """Add zones specified in graph on current proc."""
    tp = Internal.copyRef(t)
    _addXZones(tp, graph, variables, noCoordinates, cartesian, subr)
    return tp

def _addXZones(t, graph, variables=None, noCoordinates=False, cartesian=False, subr=True):
    if not graph: return t
    reqs = []
    if cartesian: import Compressor.PyTree as Compressor 
    if rank in graph:
        g = graph[rank] # graph du proc courant
        for oppNode in g:
            # Envoie les zones necessaires au noeud oppose
            #print('%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode]))
            names = g[oppNode]
            data = [] # data est une liste de zones
            for n in names:
                zone = Internal.getNodeFromName2(t, n)
                if variables is not None:
                    v = C.getVarNames(zone, excludeXYZ=True)[0]
                    for i in variables: v.remove(i)
                    if noCoordinates: v += ['CoordinateX', 'CoordinateY', 'CoordinateZ']
                    zonep = C.rmVars(zone, v)
                    if cartesian:
                        zonepc = Compressor.compressCartesian(zonep, subr=subr)
                        data.append(zonepc)
                    else: data.append(zonep)
                else:
                    if cartesian:
                        zonep = Compressor.compressCartesian(zone)
                        data.append(zonep)
                    else: data.append(zone)
            s = KCOMM.isend(data, dest=oppNode)
            reqs.append(s)

    # Reception
    for node in graph:
        #print(rank, graph[node].keys())
        if rank in graph[node]:
            #print('%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank]))
            data = KCOMM.recv(source=node)
            for z in data: # data est une liste de zones
                if cartesian:
                    import Compressor.PyTree as Compressor
                    Compressor._uncompressCartesian(z)
                #print('%d: recoit la zone %s.'%(rank,z[0]))
                Internal.createChild(z, 'XZone', 'UserDefinedData_t') 
                # Existe deja? 
                zone = Internal.getNodeFromName2(t, z[0])
                if zone is not None: # replace
                    bases = Internal.getBases(t)
                    for b in bases:
                        c = Internal.getNodePosition(zone, b)
                        if c != -1: b[2][c] = z
                else: # append to first base
                    bases = Internal.getBases(t)
                    bases[0][2].append(z)
    MPI.Request.Waitall(reqs)
    return t

#==============================================================================
# Recupere les zones specifiees dans le graph, les ajoute a l'arbre local t
# if layers not None, only communicate the desired number of layers
# bboxDict is dict with the zones of t as keys and their specific bboxes as key values, used when layers not None
# if subr=True, the ZoneSubRegions are kept during the exchange 
#==============================================================================
def addLXZones(t, graph, variables=None, cartesian=False, interDict=[], bboxDict={}, layers=2, subr=True):
    """Add zones specified in graph on current proc."""
    tp = Internal.copyRef(t)
    _addLXZones(tp, graph, variables, cartesian, interDict, bboxDict, layers, subr)
    return tp

def _addLXZones(t, graph, variables=None, cartesian=False, interDict=[], bboxDict={}, layers=2, subr=True):
    if not graph: return t
    # create bboXDict if necessary
    if not bboxDict: bboxDict = createBboxDict(t)
    reqs = []
    if cartesian: import Compressor.PyTree as Compressor 
    if rank in graph:
        g = graph[rank] # graph du proc courant
        for oppNode in g:
            # Envoie les zones necessaires au noeud oppose
            #print('%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode]))
            names = g[oppNode]
            data = [] # data est une liste de zones
            for n in names:
                zone = Internal.getNodeFromName2(t, n)
                if variables is not None:
                    v = C.getVarNames(zone, excludeXYZ=True)[0]
                    for i in variables: v.remove(i)
                    zonep = C.rmVars(zone, v)
                    if cartesian:
                        # Each block will be cut and compressed before being sent
                        # The information provided by bboxDict are crucial
                        for ni in interDict[n]:
                            if int(ni[-1]) == oppNode:
                                # bbox is the bbox of the overlapping between the two blocks
                                bbox = GetIntersectionBbox(bboxDict[n],bboxDict[ni])
                                # The copy is VERY important
                                zonepc = Internal.copyTree(zonep)
                                # call compressor with the right arguments
                                Compressor._compressCartesian(zonepc,bbox=bbox,layers=layers,subr=subr)
                                data.append(zonepc)
                    else: data.append(zonep)
                else:
                    if cartesian:
                        # Each block will be cut and compressed before being sent
                        # The information provided by bboxDict are crucial
                        for ni in interDict[n]:
                            if int(ni[-1]) == oppNode:
                                # bbox is the bbox of the overlapping between the two blocks
                                bbox = GetIntersectionBbox(bboxDict[n], bboxDict[ni])
                                # The copy is VERY important
                                zonep = Internal.copyTree(zone)
                                # call compressor with the right arguments
                                Compressor._compressCartesian(zonep, bbox=bbox, layers=layers)
                                data.append(zonep)
                    else: data.append(zone)
            s = KCOMM.isend(data, dest=oppNode)
            reqs.append(s)

    # Reception
    for node in graph:
        #print(rank, graph[node].keys())
        if rank in graph[node]:
            #print('%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank]))
            data = KCOMM.recv(source=node)
            for z in data: # data est une liste de zones
                if cartesian:
                    import Compressor.PyTree as Compressor
                    Compressor._uncompressCartesian(z)
                #print('%d: recoit la zone %s.'%(rank,z[0]))
                Internal.createChild(z, 'XZone', 'UserDefinedData_t') 
                # Existe deja? 
                zone = Internal.getNodeFromName2(t, z[0])
                if zone is not None: # replace
                    if cartesian: # a cartesian block might be splitted into several parts...
                    # Hence the counter, to avoid two identical zone names in the tree
                        cpt = 1
                        z[0] = z[0][:-2]+'#'+str(cpt)+z[0][-2:]
                        while Internal.getNodeFromName2(t, z[0]) is not None:
                            cpt += 1
                            z[0] = z[0][:-4]+'#'+str(cpt)+z[0][-2:]
                        bases = Internal.getBases(t)
                        bases[0][2].append(z)
                    else:
                        bases = Internal.getBases(t)
                        for b in bases:
                            c = Internal.getNodePosition(zone, b)
                            if c != -1: b[2][c] = z
                else: # append to first base
                    bases = Internal.getBases(t)
                    bases[0][2].append(z)
    MPI.Request.Waitall(reqs)
    return None

# Recupere les sous-zones de match de z correspondant a oppNode
def getMatchSubZones__(z, procDict, oppNode, depth):
    import Transform.PyTree as T
    dim = Internal.getZoneDim(z)
    out = []
    gcs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
    for g in gcs:
        nodes = Internal.getNodesFromType1(g, 'GridConnectivity1to1_t')
        for n in nodes:
            oppZoneName = Internal.getValue(n)
            oppNodeHere = procDict[oppZoneName]
            if oppNodeHere == oppNode:
                prange = Internal.getNodeFromName1(n, 'PointRange')
                prange = Internal.getValue(prange)
                wrange = Internal.range2Window(prange)
                imin = wrange[0]
                imax = wrange[1]
                jmin = wrange[2]
                jmax = wrange[3]
                kmin = wrange[4]
                kmax = wrange[5]
                if imin == imax and imin == 1: imax = 1+depth ; suffix = 'imin'
                elif imin == imax: imin = imax-depth ; suffix = 'imax'
                elif jmin == jmax and jmin == 1: jmax = 1+depth ; suffix = 'jmin'
                elif jmin == jmax: jmin = jmax-depth ; suffix = 'jmax'
                elif kmin == kmax and kmin == 1: kmax = 1+depth ; suffix = 'kmin'
                elif kmin == kmax: kmin = kmax-depth ; suffix = 'kmax'              
                oppZone = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax))
                oppZone[0] = z[0]+'_MX_'+oppZoneName+'-'+suffix 
                     
                z = Internal.rmNodesByName(z, 'ZoneGridConnectivity')
                Internal.createChild(oppZone, 'ZoneGridConnectivity_t', 'ZoneGridConnectivity_t')
                gcXZone = Internal.createNode('ZoneGridConnectivity_t', 'ZoneGridConnectivity_t')
                Internal._addChild(gcXZone, n)
                Internal._addChild(oppZone,gcXZone)
                
                Internal.createChild(oppZone, 'XZone', 'UserDefinedData_t')
                Internal._setLoc2Glob(oppZone, z[0], win=[imin,imax,jmin,jmax,kmin,kmax], sourceDim=[dim[1],dim[2],dim[3]])
                
                out.append(oppZone)
    return out

def updateGridConnectivity(a):
    # Update grid connectivities to be consistent with XZone (just after using addMXZones)
    zones     = Internal.getZones(a)
    zonesReal = []
    for z in zones:
        xz = Internal.getNodeFromName1(z, 'XZone')
        if xz is None: zonesReal.append(z)

    for z in zonesReal:
        gcs   = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
        for g in gcs:
            nodes = Internal.getNodesFromType1(g, 'GridConnectivity1to1_t')
            for n in nodes:
                # Recherche le nom de la bandelette en raccord 
                oppName = Internal.getValue(n)
                
                # suffix
                prange = Internal.getNodeFromName1(n, 'PointRangeDonor')
                prange = Internal.getValue(prange)
                wrange = Internal.range2Window(prange)
                imin = wrange[0] ; imax = wrange[1]
                jmin = wrange[2] ; jmax = wrange[3]
                kmin = wrange[4] ; kmax = wrange[5]
                if imin == imax and imin == 1:   suffix = 'imin'
                elif imin == imax: suffix = 'imax'
                elif jmin == jmax and jmin == 1: suffix = 'jmin'
                elif jmin == jmax: suffix = 'jmax'
                elif kmin == kmax and kmin == 1: suffix = 'kmin'
                elif kmin == kmax: suffix = 'kmax' 
                
                zopp = Internal.getNodeFromName(a, oppName+'_MX_'+z[0]+'-'+suffix)
                
                if zopp is not None:
                                           
                    Internal.setValue(n, zopp[0]) # renommage
            
                    src, loc2glob = Internal.getLoc2Glob(zopp)

                    # Update current zone
                    prd    = Internal.getNodeFromName1(n, 'PointRangeDonor')
                    p      = Internal.range2Window(prd[1])
                    p      = [p[0]-loc2glob[0]+1,p[1]-loc2glob[0]+1,p[2]-loc2glob[2]+1,p[3]-loc2glob[2]+1,p[4]-loc2glob[4]+1,p[5]-loc2glob[4]+1]
                    p      = Internal.window2Range(p)
                    Internal.setValue(prd, p)

                    # Update XZone 
                    gcopp  = Internal.getNodesFromType1(zopp, 'ZoneGridConnectivity_t')
                    match  = Internal.getNodesFromType1(gcopp, 'GridConnectivity1to1_t')[0] # 1 seul match dans les XZone
                
                    pr     = Internal.getNodeFromName1(match, 'PointRange')
                    Internal.setValue(pr, p) 
                                  
    return a 

def _revertMXGridConnectivity(a):
    # Restore grid connectivities with respect to real zone (after using addMXZones)
    zones     = Internal.getZones(a)
    zonesReal = []
    for z in zones:
        xz = Internal.getNodeFromName1(z, 'XZone')
        if xz is None: zonesReal.append(z)

    for z in zonesReal:
        gcs   = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
        for g in gcs:  
            nodes = Internal.getNodesFromType1(g, 'GridConnectivity1to1_t')
            for n in nodes:
                # Recherche le nom de la bandelette en raccord 
                oppName = Internal.getValue(n)
                zopp    = Internal.getNodeFromName(a, oppName)
                xzopp   = Internal.getNodeFromName1(zopp, 'XZone')

                if xzopp is not None:
                    newName = oppName.split('_MX_')[0]
                    Internal.setValue(n, newName)
      
                    src, loc2glob = Internal.getLoc2Glob(zopp)

                    # Update current zone
                    prd    = Internal.getNodeFromName1(n, 'PointRangeDonor')
                    p      = Internal.range2Window(prd[1])
                    p      = [p[0]+loc2glob[0]-1,p[1]+loc2glob[0]-1,p[2]+loc2glob[2]-1,p[3]+loc2glob[2]-1,p[4]+loc2glob[4]-1,p[5]+loc2glob[4]-1]
                    p      = Internal.window2Range(p)
                    Internal.setValue(prd, p)
                                  
    return None

def _revertBXGridConnectivity(a):
    # Restore grid connectivities with respect to real zone (after using addBXZones)
    zones     = Internal.getZones(a)
    zonesReal = []
    for z in zones:
        xz = Internal.getNodeFromName1(z, 'XZone')
        if xz is None: zonesReal.append(z)

    for z in zonesReal:
        gcs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
        for g in gcs:  
            nodes = Internal.getNodesFromType1(g, 'GridConnectivity1to1_t')
            for n in nodes:
                # Recherche le nom de la bandelette en raccord 
                oppName = Internal.getValue(n)
                zopp    = Internal.getNodeFromName(a, oppName)
                if zopp is not None:
                    xzopp = Internal.getNodeFromName1(zopp, 'XZone')
                    if xzopp is not None:
                        newName = oppName[:len(oppName)-2]
                        Internal.setValue(n, newName)
                        src, loc2glob = Internal.getLoc2Glob(zopp)

                        # Update current zone
                        prd    = Internal.getNodeFromName1(n, 'PointRangeDonor')
                        p      = Internal.range2Window(prd[1])
                        p      = [p[0]+loc2glob[0]-1,p[1]+loc2glob[0]-1,p[2]+loc2glob[2]-1,p[3]+loc2glob[2]-1,p[4]+loc2glob[4]-1,p[5]+loc2glob[4]-1]
                        #check if not same window - can happen if near match exists
                        remove=False
                        if z[0] == newName:
                            remove=True
                            prr = Internal.getNodeFromName1(n,'PointRange')
                            p2 = Internal.range2Window(prr[1])
                            for i in range(6):
                                if p2[i] != p[i]: 
                                    remove=False
                                    break
                        if not remove:
                            p = Internal.window2Range(p)
                            Internal.setValue(prd, p)                   
                        else:
                            Internal._rmNodesByName(z,n[0])

            nodes = Internal.getNodesFromType1(g, 'GridConnectivity_t')
            for n in nodes:
                gctype = Internal.getNodeFromType(n,'GridConnectivityType_t')
                gctype = Internal.getValue(gctype)
                if gctype=='Abutting':
                    # Recherche le nom de la bandelette en raccord 
                    oppName = Internal.getValue(n)
                    zopp    = Internal.getNodeFromName(a, oppName)
                    if zopp is not None:
                        xzopp = Internal.getNodeFromName1(zopp, 'XZone')
                        if xzopp is not None:
                            newName = oppName[:len(oppName)-2]
                            Internal.setValue(n, newName)

    return None
    
# Ajoute des sous-zones correspondant aux raccords sur un arbre distribue
def _addMXZones(a, depth=2):

    graph = computeGraph(a, type='match')
    procDict = getProcDict(a)
    reqs = []
    bases = Internal.getBases(a)
    for b in bases:
        zones = Internal.getZones(b)
        if rank in graph:
            g = graph[rank] # graph du proc courant
            for oppNode in g:
                data = []
                for z in zones:
                    zs = getMatchSubZones__(z, procDict, oppNode, depth)
                    for zl in zs: zl[0] = b[0]+'/'+zl[0]
                    data += zs
                s = KCOMM.isend(data, dest=oppNode)
                reqs.append(s)
    for node in graph:
        if rank in graph[node]:
            data = KCOMM.recv(source=node)
            for d in data:
                (baseName, zoneName) = d[0].split('/',1)
                d[0] = zoneName
                b = Internal.getNodeFromName1(a, baseName)
                if b is None: b = Internal.newCGNSBase(baseName, parent=a)
                b[2].append(d)
    MPI.Request.Waitall(reqs)

    a = updateGridConnectivity(a)
    return None
    
# IN: bb0 et bb1: [xmin,ymin,zmin,xmax,ymax,zmax]
# Retourne true si les bbox s'intersectent
def inters(bb0, bb1, tol=0.):
    if bb0[3] < bb1[0]-tol or bb0[0] > bb1[3]+tol: return False
    if bb0[4] < bb1[1]-tol or bb0[1] > bb1[4]+tol: return False
    if bb0[5] < bb1[2]-tol or bb0[2] > bb1[5]+tol: return False
    return True

# creation de la bandelette
def subzone(a, indMin, indMax, supp):
    import Transform.PyTree as T
    dim = Internal.getZoneDim(a)
    (imin,jmin,kmin) = (indMin[0],indMin[1],indMin[2])
    (imax,jmax,kmax) = (indMax[0],indMax[1],indMax[2])
    ap = T.subzone(a, indMin, indMax)
    ap[0] = a[0] + supp
    Internal.createChild(ap, 'XZone', 'UserDefinedData_t')
    Internal._setLoc2Glob(ap, a[0], win=[imin,imax,jmin,jmax,kmin,kmax], sourceDim=[dim[1],dim[2],dim[3]])
    return ap
    
# Ajoute les bandelettes des autres procs sur le procs locaux
# si allB=True, les 6 bandelettes de chaque zone voisine sont ramenees (necessaire pour connectMatchPeriodic).
def _addBXZones(a, depth=2, allB=False):
    import Generator.PyTree as G
    # Calcul des bbox des zones locales
    zones = Internal.getZones(a)
    bbz = {}
    for z in zones:
        bb = G.bbox(z)
        bbz[z[0]] = bb

    # Calcul des bandelettes et de leur bbox
    bbsz = {}; sz = {}

    bases = Internal.getBases(a)
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            dim = Internal.getZoneDim(z)
            ni = dim[1]; nj = dim[2]; nk = dim[3]
            # b1 = subzone(z, (1,1,1), (min(depth,ni),nj,nk), 'S1')
            # b2 = subzone(z, (max(ni-depth+1,1),1,1), (ni,nj,nk), 'S2')
            # b3 = subzone(z, (1,1,1), (ni,min(depth,nj),nk), 'S3')
            # b4 = subzone(z, (1,max(nj-depth+1,1),1), (ni,nj,nk), 'S4')
            # b5 = subzone(z, (1,1,1), (ni,nj,min(depth,nk)), 'S5')
            # b6 = subzone(z, (1,1,max(nk-depth+1,1)), (ni,nj,nk), 'S6')
 
            # Bandelettes non recouvrantes
            b1 = subzone(z, (1,1,1), (min(depth,ni),nj,nk), 'S1')
            b2 = subzone(z, (max(ni-depth+1,1),1,1), (ni,nj,nk), 'S2')
            b3 = subzone(z, (min(depth,ni),1,1), (max(ni-depth+1,1),min(depth,nj),nk), 'S3') # CW
            b4 = subzone(z, (min(depth,ni),max(nj-depth+1,1),1), (max(ni-depth+1,1),nj,nk), 'S4') # CW 
            b5 = subzone(z, (min(depth,ni),min(depth,nj),1), (max(ni-depth+1,1),max(nj-depth+1,1),min(depth,nk)), 'S5') # CW
            b6 = subzone(z, (min(depth,ni),min(depth,nj),max(nk-depth+1,1)), (max(ni-depth+1,1),max(nj-depth+1,1),nk), 'S6') # CW
            
            sz[b1[0]] = b1
            sz[b2[0]] = b2
            sz[b3[0]] = b3
            sz[b4[0]] = b4
            sz[b5[0]] = b5
            sz[b6[0]] = b6
            bbsz[b1[0]] = G.bbox(b1)
            bbsz[b2[0]] = G.bbox(b2)
            bbsz[b3[0]] = G.bbox(b3)
            bbsz[b4[0]] = G.bbox(b4)
            bbsz[b5[0]] = G.bbox(b5)
            bbsz[b6[0]] = G.bbox(b6)

            b1[0] = b[0]+'/'+b1[0]
            b2[0] = b[0]+'/'+b2[0]
            b3[0] = b[0]+'/'+b3[0]
            b4[0] = b[0]+'/'+b4[0]
            b5[0] = b[0]+'/'+b5[0]
            b6[0] = b[0]+'/'+b6[0]
            
    # allgather des bbox des bandelettes
    bboxes = KCOMM.allgather(bbz)

    # Echange par alltoall
    data   = []
    zone_i = []
    for i in range(size):
        data_i = {}
        if i != rank:
            for bb in bboxes[i]:
                for k in bbsz:
                    #SP : marge a 1e-6 lorsque les zones sont de meme z par exemple
                    if inters(bbsz[k], bboxes[i][bb], tol=1.e-6):
                        if not allB:
                            data_i[k] = sz[k]
                        else:
                            zname = k[:len(k)-2] 
                            if zname not in zone_i:
                                zone_i.append(zname) 
                                
            if allB:                    
                for z in zone_i:
                    data_i[z+'S1'] = sz[z+'S1']
                    data_i[z+'S2'] = sz[z+'S2']
                    data_i[z+'S3'] = sz[z+'S3']
                    data_i[z+'S4'] = sz[z+'S4']
                    data_i[z+'S5'] = sz[z+'S5']
                    data_i[z+'S6'] = sz[z+'S6']
           
        data.append(data_i)

    datar = COMM_WORLD.alltoall(data)
    for p in datar:
        for q in p:
            z = p[q]
            zoneName = z[0]
            (baseName, zoneName) = zoneName.split('/',1)
            b = Internal.getNodeFromName1(a, baseName)
            if b is None: b = Internal.newCGNSBase(baseName, parent=a)
            z[0] = zoneName
            b[2].append(z)
    
    return None

#==============================================================================
# Supprime les zones ajoutees par addXZones
#==============================================================================
def rmXZones(t):
    """Remove zones added by addXZones."""
    tp = Internal.copyRef(t)
    _rmXZones(tp)
    return tp

def _rmXZones(t):
    bases = Internal.getBases(t)
    if bases == []:
        zones = Internal.getZones(t)
        for z in zones:
            tag = Internal.getNodeFromName1(z, 'XZone')
            if tag is not None:
                (p, c) = Internal.getParentOfNode(t, z)
                del p[2][c]
    else:
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                tag = Internal.getNodeFromName1(z, 'XZone')
                if tag is not None:
                    c = Internal.getNodePosition(z, b)
                    del b[2][c]
    return None

def _rmMXZones(t):
    _revertMXGridConnectivity(t)
    _rmXZones(t)
    return None

def _rmBXZones(t):
    _revertBXGridConnectivity(t)
    _rmXZones(t)
    return None

#==============================================================================
# Retourne le dictionnaire proc pour chaque zoneName
# IN: t: full/S/LS/Partial
#==============================================================================
def getProcDict(t, reduction=True):
    """Return the dictionary proc['zoneName']."""
    d1 = Distributed.getProcDict(t)
    if reduction:
        d = KCOMM.allgather(d1)
        items = []
        for i in d: items += i.items()
        procDict = dict(items)
    else: procDict = d1
    return procDict

#==============================================================================
# IN: z: zone node
#==============================================================================
def getProc(z):
    """Return the proc where zone is affected to."""
    return Distributed.getProc(z)

#==============================================================================
# IN: t: zone node, tree, ...
# IN: rank: number of rank to set
#==============================================================================
def setProc(t, rank):
    """Set the proc number to a zone or a set of zones."""
    tp = Internal.copyRef(t)
    Distributed._setProc(tp, rank)
    return tp

def _setProc(t, rank):
    Distributed._setProc(t, rank)
    return None
