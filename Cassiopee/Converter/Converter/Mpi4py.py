# Interface de Cassiopee pour mpi
from . import PyTree as C
from . import Internal
from . import Distributed
from . import converter

# Acces a Distributed
from .Distributed import readZones, _readZones, convert2PartialTree, _convert2PartialTree, convert2SkeletonTree, readNodesFromPaths, readPyTreeFromPaths, writeNodesFromPaths, mergeGraph, splitGraph

__all__ = ['rank', 'size', 'KCOMM', 'COMM_WORLD', 'SUM', 'MIN', 'MAX', 'LAND',
           'setCommunicator', 'barrier', 'send', 'recv', 'sendRecv', 'sendRecvC',
           'bcast', 'Bcast', 'gather', 'Gather',
           'reduce', 'Reduce', 'allreduce', 'Allreduce',
           'bcastZone', 'gatherZones', 'allgatherZones',
           'createBBTree', 'intersect', 'intersect2', 'allgatherDict',
           'allgather', 'readZones', 'writeZones', 'convert2PartialTree',
           'convert2SkeletonTree',
           'readNodesFromPaths', 'readPyTreeFromPaths', 'writeNodesFromPaths',
           'allgatherTree', 'convertFile2SkeletonTree', 'convertFile2PyTree',
           'convertPyTree2File', 'seq', 'print0', 'printA',
           'createBboxDict', 'computeGraph', 'addXZones',
           '_addXZones', '_addMXZones', '_addBXZones', '_addLXZones',
           'rmXZones', '_rmXZones', '_rmMXZones', '_rmBXZones', 'getProcDict',
           'getProc', 'setProc', '_setProc', 'getPropertyDict', 'getProperty', 'COMM_WORLD']

from mpi4py import MPI
import numpy
import os

COMM_WORLD = MPI.COMM_WORLD
KCOMM = COMM_WORLD

rank = KCOMM.rank
size = KCOMM.size

SUM = MPI.SUM
MAX = MPI.MAX
MIN = MPI.MIN
LAND = MPI.LAND

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
# barrier
#==============================================================================
def barrier():
    KCOMM.barrier()

#=================================================================================
# send - blocking send data from a proc to another proc (using pickle, small data)
#=================================================================================
def send(obj, dest=None, tag=0):
    KCOMM.send(obj, dest, tag=tag)

#=================================================================================
# Send - blocking send data from a proc to another proc (for numpys)
#=================================================================================
def Send(obj, dest=None, tag=0):
    KCOMM.Send(obj, dest, tag=tag)

#=====================================================================================
# isend - non blocking send data from a proc to another proc (using pickle, small data)
#=====================================================================================
def isend(obj, dest=None, tag=0):
    KCOMM.isend(obj, dest, tag=tag)

#=====================================================================================
# iSend - non blocking send data from a proc to another proc (for numpys)
#=====================================================================================
def Isend(obj, dest=None, tag=0):
    KCOMM.Isend(obj, dest, tag=tag)

#==============================================================================
# receive - receive data from a proc
#==============================================================================
def recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG):
    return KCOMM.recv(source=source, tag=tag)

#==============================================================================
# Reduce to root (using pickle, small data)
#==============================================================================
def reduce(data, op=MPI.SUM, root=0):
    return KCOMM.reduce(data, op=op, root=root)

#==============================================================================
# Reduce to root (for numpy data)
# dataIn and dataOut must be numpys of same size and type
# or buffer specified [data,count,MPI_DOUBLE]
#==============================================================================
def Reduce(dataIn, dataOut, op=MPI.SUM):
    return KCOMM.Reduce(dataIn, dataOut, op=op)

#==============================================================================
# Reduce to all (using pickle, small data)
#==============================================================================
def allreduce(data, op=MPI.SUM):
    return KCOMM.allreduce(data, op=op)

#==============================================================================
# Reduce to all (for numpy data)
# dataIn and dataOut must be numpys of same size and type
# or buffer specified [data,count,MPI_DOUBLE]
#==============================================================================
def Allreduce(dataIn, dataOut, op=MPI.SUM):
    return KCOMM.Allreduce(dataIn, dataOut, op=op)

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
            if oppNode in datas: s = KCOMM.isend(datas[oppNode], dest=oppNode)
            else: s = KCOMM.isend(None, dest=oppNode)
            reqs.append(s)
    barrier()
    rcvDatas={}
    for node in graph:
        #print(rank, graph[node],graph[node].keys(),flush=True)
        if rank in graph[node]:
            #print('%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank]),flush=True)
            rec = KCOMM.recv(source=node)
            if rec is not None: rcvDatas[node] = rec
    MPI.Request.waitall(reqs)
    return rcvDatas

#==============================================================================
# Send and receive with a graph - C version - no pickle
# IN: datas: un dictionnaire des donnees a envoyer par proc de destination
# OUT: un dictionnaire des donnees recues par proc d'origine
# Attention: ne fonctionne que pour certaines datas (issues de transfer)
#==============================================================================
def sendRecvC(datas, graph):
    if graph == {}: return {}
    reqs = []

    if rank in graph:
        g = graph[rank] # graph du proc courant
        for oppNode in g:
            # Envoie les datas necessaires au noeud oppose
            #print('%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode]))
            if oppNode in datas:
                s = converter.iSend(datas[oppNode], oppNode, rank, KCOMM)
            else:
                s = converter.iSend(None, oppNode, rank, KCOMM)
            reqs.append(s)
    barrier()
    rcvDatas={}
    for node in graph:
        #print(rank, graph[node].keys())
        if rank in graph[node]:
            #print('%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank]))
            rec = converter.recv(node, rank, KCOMM)
            if rec is not None: rcvDatas[node] = rec

    a = converter.waitAll(reqs)
    return rcvDatas

#==============================================================================
# Construction d'un arbre de recherche pour des BBox
# IN : tBB arbre cgns de bbox
# OUT : objet C BBtree (hook)
#==============================================================================
def createBBTree(t):
    zones = Internal.getZones(t)
    minBBoxes = [] ; maxBBoxes = []
    for z in zones:
        # BBox de la zone
        gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        xCoords = Internal.getNodeFromName1(gc, 'CoordinateX')[1]
        yCoords = Internal.getNodeFromName1(gc, 'CoordinateY')[1]
        zCoords = Internal.getNodeFromName1(gc, 'CoordinateZ')[1]
        minBBoxes.append([numpy.min(xCoords), numpy.min(yCoords), numpy.min(zCoords)])
        maxBBoxes.append([numpy.max(xCoords), numpy.max(yCoords), numpy.max(zCoords)])

    return converter.createBBTree(minBBoxes, maxBBoxes)

#==============================================================================
# Recherche des intersections entre une bbox de zone et un arbre de BBox
# IN : Bbox d'une zone + arbre de recherche de BBox
# OUT : tableau de taille de nombre des BBox avec True or False
#==============================================================================
def intersect(zone, BBTree):
    gc = Internal.getNodeFromName1(zone, Internal.__GridCoordinates__)
    xCoords = Internal.getNodeFromName1(gc, 'CoordinateX')[1]
    yCoords = Internal.getNodeFromName1(gc, 'CoordinateY')[1]
    zCoords = Internal.getNodeFromName1(gc, 'CoordinateZ')[1]
    minBBox = [numpy.min(xCoords), numpy.min(yCoords), numpy.min(zCoords)]
    maxBBox = [numpy.max(xCoords), numpy.max(yCoords), numpy.max(zCoords)]
    return converter.intersect(minBBox, maxBBox, BBTree)

def intersect2(t, BBTree):
    zones = Internal.getZones(t)
    inBB = numpy.empty((6*len(zones)), dtype=numpy.float64)
    for c, z in enumerate(zones):
        gc = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
        xCoords = Internal.getNodeFromName1(gc, 'CoordinateX')[1]
        yCoords = Internal.getNodeFromName1(gc, 'CoordinateY')[1]
        zCoords = Internal.getNodeFromName1(gc, 'CoordinateZ')[1]
        inBB[6*c  ] = xCoords[0,0,0]
        inBB[6*c+1] = yCoords[0,0,0]
        inBB[6*c+2] = zCoords[0,0,0]
        inBB[6*c+3] = xCoords[1,0,0]
        inBB[6*c+4] = yCoords[0,1,0]
        inBB[6*c+5] = zCoords[0,0,1]

    return converter.intersect2(inBB, BBTree)

#==============================================================================
# Recherche des zones fixes non intersectees et ajout dans le dict
# IN : liste des zones fixes, dict des intersects
# OUT : dict a jour (reference ou pas?)
#==============================================================================
#def fillDict(zones, dictIntersect):
#    for zone in zones:
#        if zone[0] not in dictIntersect: dictIntersect[zone[0]] = []
#    return dictIntersect

#==============================================================================
# allGather dictionnaire (rejete les doublons de keys et values)
#==============================================================================
def allgatherDict(data):
    ret = KCOMM.allgather(data)
    if isinstance(data, dict):
        out = {}
        for r in ret:
            for k in r:
                # if rank==0: print('PROC%d : k=%s'%(rank, k), flush=True)
                if k not in out:
                    out[k] = []
                    for data in r[k]:
                        if data not in out[k]:
                            # if rank==0: print('PROC%d : data=%s ; out[k] = '%(rank, k), out[k], flush=True)
                            out[k].append(data)
                            # if rank==0: print('PROC%d : out[k] = '%(rank), out[k], flush=True)

                else:
                    for data in r[k]:
                        # if rank==0: print('PROC%d : data=%s ; out[k] = '%(rank, k), out[k], flush=True)
                        if data not in out[k]: out[k].append(data)
                        # if rank==0: print('PROC%d : out[k] = '%(rank), out[k], flush=True)
        return out
    else: return ret

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
# gather to root from all procs
#==============================================================================
# data=All with pickle
def gather(data, root=0):
    return KCOMM.gather(data, root)

# data=numpy
def Gather(data, root=0):
    return KCOMM.Gather(data, root)

#==============================================================================
# bcast from root to all procs
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
# if variables == [] and coord == True, envoie uniquement les coords de z
def bcastZone(z, root=0, coord=True, variables=[]):
    # zp = squelette envoye en pickle
    if rank == root:
        zp = Internal.copyRef(z)
        Internal._rmNodesFromType(zp, 'GridCoordinates_t')
        Internal._rmNodesFromType(zp, 'FlowSolution_t')
        Internal._rmNodesFromType(zp, 'ZoneSubRegion_t')
    else: zp = None

    zp = KCOMM.bcast((zp), root)

    if coord:
        if rank == root:
            px = Internal.getNodeFromName2(z, 'CoordinateX')[1]
            py = Internal.getNodeFromName2(z, 'CoordinateY')[1]
            pz = Internal.getNodeFromName2(z, 'CoordinateZ')[1]
            sx = px.shape; sy = py.shape; sz = pz.shape
        else:
            sx = None; sy = None; sz = None

        sx, sy, sz = KCOMM.bcast((sx, sy, sz), root)

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

    for cpt, var in enumerate(variables):
        loc = 'centers' if 'centers' in var else 'nodes'
        var = var.replace('centers:','')

        if rank == root:
            pv = Internal.getNodeFromName2(z, var)[1]
            sv = pv.shape
        else:
            sv = None

        sv = KCOMM.bcast((sv), root)

        if rank != root:
            pv = numpy.empty(sv, dtype=numpy.float64, order='F')

        KCOMM.Bcast([pv,MPI.DOUBLE], root)

        if rank != root:
            # Reconstruction de la zone
            flowSol = Internal.__FlowSolutionNodes__ if loc == 'nodes' else Internal.__FlowSolutionCenters__
            if cpt == 0:
                Internal._createUniqueChild(zp, flowSol, 'FlowSolution_t')
                n = Internal.getNodeFromName1(zp, flowSol)
                if loc == 'centers': Internal._createUniqueChild(n, 'GridLocation', 'GridLocation_t', 'CellCenter')
            else:
                n = Internal.getNodeFromName1(zp, flowSol)
            Internal._createUniqueChild(n, var, 'DataArray_t', value=pv)

    if rank == root:
        zp = Internal.copyRef(z)
        # suppression coord si coord=False et des autres champs
        if not coord: Internal._rmNodesFromType(zp, 'GridCoordinates_t')
        if variables == []: Internal._rmNodesFromType(zp, 'FlowSolution_t')
        else:
            varszp = C.getVarNames(zp, excludeXYZ=True, loc='both')[0]
            for var in varszp:
                if var not in variables:
                    Internal._rmNodesFromName(zp, var.replace('centers:',''))
    return zp

# All gather une liste de zones, recuperation identique sur tous les procs
# dans une liste a plat
# Partage les coordonnees si coord=True
# Partage les variables renseignees et supprime les autres
def allgatherZones(zones, coord=True, variables=[]):
    # Chaque processeur bcast ses zones vers les autres ranks
    zones = Internal.getZones(zones)
    lenZones = KCOMM.allgather(len(zones))
    allZones = []
    for i in range(size):
        for cz in range(lenZones[i]):
            if rank == i:
                zp = bcastZone(zones[cz], root=i, coord=coord, variables=variables)
            #    if variables == []:
            #        Internal._rmNodesFromType(zp, 'FlowSolution_t')
            #    else:
            #        varszp = C.getVarNames(zp, excludeXYZ=True, loc='both')[0]
            #        for var in varszp:
            #            if var not in variables: Internal._rmNodesFromName(zp, var.replace('centers:',''))
            else:
                zp = bcastZone(None, root=i, coord=coord, variables=variables)

            #if rank == i: zp = zones[cz]
            #else: zp = None
            #zp = bcastZone(zp, root=i, coord=coord, variables=variables)

            allZones.append(zp)
    return allZones

# Envoie les zones de chaque proc au proc root dans une liste a plat des zones
def gatherZones(zones, root=0):
    zones = Internal.getZones(zones)
    ret = KCOMM.gather(zones, root)
    out = []
    if ret is not None:
        for i in ret: out += i
    return out

#==============================================================================
# Lecture du squelette d'un arbre dans un fichier
# Lecture proc 0 + bcast
#==============================================================================
def convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5,
                             maxDepth=-1, links=None):
    """Read a file and return a skeleton tree."""
    if rank == 0: exists = int(os.path.exists(fileName))
    else: exists = 1
    exists = bcast(exists, root=0)
    if not exists: raise IOError("convertFile2SkeletonTree: file %s not found."%fileName)

    if rank == 0: t = Distributed.convertFile2SkeletonTree(fileName, format, maxFloatSize, maxDepth, None, links)
    else: t = None
    t = KCOMM.bcast(t)
    if links is not None:
        lk = KCOMM.bcast(links)
        if rank > 0: links += lk
    return t

#==============================================================================
# Lecture complete d'un arbre dans un fichier
# si proc=None, lecture proc 0 + bcast
# sinon lecture des zones correspondant a proc
#==============================================================================
def convertFile2PyTree(fileName, format=None, proc=None):
    """Read a file and return a full tree or partial tree."""
    if rank == 0: exists = int(os.path.exists(fileName))
    else: exists = 1
    exists = bcast(exists, root=0)
    if not exists: raise IOError("convertFile2PyTree: file %s not found."%fileName)

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
        for b in Internal.getBases(tp):
            out = []
            for i in b[2]:
                if i[3] != 'Zone_t': out.append(i)
            b[2] = out
        KCOMM.send(tp, dest=0)
        return t
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
# Calcule le dictionnaire des bbox de l'arbre complet
# Utile pour addXZones optimises
#==============================================================================
def createBboxDict(t):
    try: import Generator.PyTree as G
    except: raise ImportError("createBboxDict requires Generator module.")
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
                    graph[k][j] = sorted(list(set(graph[k][j])))

    return graph

#=============================================================================
# Calcule l'intersection de deux bbox
#=============================================================================
def GetIntersectionBbox(bbox1, bbox2):
    Ibbox = numpy.zeros(6, dtype=numpy.float64)
    Ibbox[0:3] = [max(bbox1[i],bbox2[i]) for i in range(3)]
    Ibbox[3:6] = [min(bbox1[i],bbox2[i]) for i in range(3,6)]
    return Ibbox

#==============================================================================
# Recupere les zones specifiees dans le graph, les ajoute a l'arbre local t
# IN: variables: None (send all vars) ou ['Density', ...] ou []
# IN: noCoordinates: False (Coordinates send), True (Coordinates not send)
# IN: keepOldNodes: send all other nodes than Coordinates and variables
# IN: cartesian: send compress coordinates for cartesian grids
# if subr=True, the ZoneSubRegions are sent
# if zoneGC=True, the ZoneGridConnectivity are sent
#==============================================================================
def addXZones(t, graph, variables=None, noCoordinates=False,
              cartesian=False, subr=True,
              keepOldNodes=True, zoneGC=True):
    """Add zones specified in graph on current proc."""
    tp = Internal.copyRef(t)
    _addXZones(tp, graph, variables, noCoordinates, cartesian, subr,
               keepOldNodes, zoneGC)
    return tp

def _addXZones(t, graph, variables=None, noCoordinates=False,
               cartesian=False, subr=True,
               keepOldNodes=True, zoneGC=True):
    """Add zones specified in graph on current proc."""
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
                (base,c) = Internal.getParentOfNode2(t, zone)

                if not keepOldNodes:
                    if variables is None: vars = C.getVarNames(zone, excludeXYZ=True)[0]
                    elif variables == []: vars = None
                    else: vars = variables
                    zonep = C.extractVars(zone, vars=variables, keepOldNodes=False)
                    if noCoordinates: C._rmVars(zonep, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
                    if cartesian: Compressor._compressCartesian(zonep)
                    if zoneGC:
                        zGC = Internal.getNodeFromType1(zone, 'ZoneGridConnectivity_t')
                        if zGC is not None: zonep[2].append(zGC)
                elif variables is not None: # all nodes but vars and coordinates
                    v = C.getVarNames(zone, excludeXYZ=True)[0]
                    for i in variables:
                        if i in v: v.remove(i)
                    if noCoordinates: v += ['CoordinateX', 'CoordinateY', 'CoordinateZ']
                    zonep = C.rmVars(zone, v)
                    if cartesian: Compressor._compressCartesian(zonep, subr=subr)
                else: # full zone
                    zonep = Internal.copyRef(zone)
                    if cartesian: Compressor._compressCartesian(zonep)
                if base is not None: zonep[0] = base[0]+'/'+zone[0]
                data.append(zonep)
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
                Internal.createChild(z, 'XZone', 'UserDefinedData_t')

                ret = z[0].split('/',1)
                if len(ret) == 2:
                    baseName = ret[0]; zoneName = ret[1]
                    z[0] = zoneName
                    base = Internal.getNodeFromName1(t, baseName)
                    if base is None:
                        if base is None: base = Internal.newCGNSBase(baseName, parent=t)
                    base[2].append(z)
                else:
                    #print('%d: recoit la zone %s.'%(rank,z[0]))
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
                    for i in variables:
                        if i in v: v.remove(i)
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
                imin = wrange[0] ; imax = wrange[1]
                jmin = wrange[2] ; jmax = wrange[3]
                kmin = wrange[4] ; kmax = wrange[5]
                if imin == imax and imin == 1: imax = 1+depth ; suffix = 'imin'+str(jmin)+str(kmin)
                elif imin == imax: imin = imax-depth ; suffix = 'imax'+str(jmin)+str(kmin)
                elif jmin == jmax and jmin == 1: jmax = 1+depth ; suffix = 'jmin'+str(imin)+str(kmin)
                elif jmin == jmax: jmin = jmax-depth ; suffix = 'jmax'+str(imin)+str(kmin)
                elif kmin == kmax and kmin == 1: kmax = 1+depth ; suffix = 'kmin'+str(imin)+str(jmin)
                elif kmin == kmax: kmin = kmax-depth ; suffix = 'kmax'+str(imin)+str(jmin)
                oppZone = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax))
                oppZone[0] = z[0]+'_MX_'+oppZoneName+'-'+suffix
                z = Internal.rmNodesByName(z, 'ZoneGridConnectivity')
                Internal.createChild(oppZone, 'ZoneGridConnectivity_t', 'ZoneGridConnectivity_t')
                gcXZone = Internal.createNode('ZoneGridConnectivity_t', 'ZoneGridConnectivity_t')
                Internal._addChild(gcXZone, n)
                Internal._addChild(oppZone, gcXZone)

                Internal.createChild(oppZone, 'XZone', 'UserDefinedData_t')
                Internal._setLoc2Glob(oppZone, z[0], win=[imin,imax,jmin,jmax,kmin,kmax], sourceDim=[dim[1],dim[2],dim[3]])
                out.append(oppZone)
    return out

def _updateGridConnectivity(a):
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
                if imin == imax and imin == 1: suffix = 'imin'+str(jmin)+str(kmin)
                elif imin == imax: suffix = 'imax'+str(jmin)+str(kmin)
                elif jmin == jmax and jmin == 1: suffix = 'jmin'+str(imin)+str(kmin)
                elif jmin == jmax: suffix = 'jmax'+str(imin)+str(kmin)
                elif kmin == kmax and kmin == 1: suffix = 'kmin'+str(imin)+str(jmin)
                elif kmin == kmax: suffix = 'kmax'+str(imin)+str(jmin)

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

    return None

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
                if gctype == 'Abutting':
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
def _addMXZones(a, depth=2, variables=None, noCoordinates=False, keepOldNodes=True):

    graph = computeGraph(a, type='match')
    bases = Internal.getBases(a)
    procDict = getProcDict(a)
    reqs = []
    if rank in graph:
        g = graph[rank] # graph du proc courant
        for oppNode in g:
            data = []
            for b in bases:
                zones = Internal.getZones(b)
                for z in zones:
                    zs = getMatchSubZones__(z, procDict, oppNode, depth)
                    #for zp in zs: print(z[0], 'found',zp[0])
                    if not keepOldNodes:
                        if variables is None: vars = C.getVarNames(zs, excludeXYZ=True)[0]
                        elif variables == []: vars = None
                        else: vars = variables
                        zsp = C.extractVars(zs, vars=variables, keepOldNodes=False)
                        if noCoordinates: C._rmVars(zsp, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
                        for c, z in enumerate(zsp):
                            ns = Internal.getNodeFromName1(zs[c], 'XZone')
                            if ns is not None: z[2].append(ns)
                            ns = Internal.getNodeFromName1(zs[c], '.Solver#ownData')
                            if ns is not None:
                                if 'Parameter_int' in ns[2]: ns[2].remove('Parameter_int')
                                if 'Parameter_real' in ns[2]: ns[2].remove('Parameter_real')
                                z[2].append(ns)
                            ns = Internal.getNodesFromType1(zs[c], 'ZoneGridConnectivity_t')
                            z[2] += ns
                        zs = zsp
                    elif variables is not None:
                        v = C.getVarNames(zs, excludeXYZ=True)[0]
                        for i in variables:
                            if i in v: v.remove(i)
                        if noCoordinates: v += ['CoordinateX', 'CoordinateY', 'CoordinateZ']
                        C._rmVars(zs, v)
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

    _updateGridConnectivity(a)

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
# IN: variables: None (all vars), ['Density'], []
def _addBXZones(a, depth=2, allB=False, variables=None):
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

            ri1 = min(depth,ni)
            ri2 = max(ni-depth+1,1)
            rj1 = min(depth,nj)
            rj2 = max(nj-depth+1,1)
            rk1 = min(depth,nk)
            rk2 = max(nk-depth+1,1)

            if ri2 - ri1 < 2: rip1 = int(ni/2)-1; rip2 = int(ni/2)+1
            else: rip1 = ri1; rip2 = ri2

            if rj2 - rj1 < 2: rjp1 = int(nj/2)-1; rjp2 = int(nj/2)+1
            else: rjp1 = rj1; rjp2 = rj2

            # Bandelettes non recouvrantes
            b1 = subzone(z, (1,1,1), (ri1,nj,nk), 'S1')
            b2 = subzone(z, (ri2,1,1), (ni,nj,nk), 'S2')
            b3 = subzone(z, (rip1,1,1), (rip2,rj1,nk), 'S3')
            b4 = subzone(z, (rip1,rj2,1), (rip2,nj,nk), 'S4')
            b5 = subzone(z, (rip1,rjp1,1), (rip2,rjp2,rk1), 'S5')
            b6 = subzone(z, (rip1,rjp1,rk2), (rip2,rjp2,nk), 'S6')

            if variables is not None: # no var
                b1 = C.extractVars(b1, vars=variables, keepOldNodes=False)
                b2 = C.extractVars(b2, vars=variables, keepOldNodes=False)
                b3 = C.extractVars(b3, vars=variables, keepOldNodes=False)
                b4 = C.extractVars(b4, vars=variables, keepOldNodes=False)
                b5 = C.extractVars(b5, vars=variables, keepOldNodes=False)
                b6 = C.extractVars(b6, vars=variables, keepOldNodes=False)

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


#==============================================================================
# Retourne le dictionnaire proc pour chaque zoneName
# IN: t: full/S/LS/Partial
#==============================================================================
def getPropertyDict(t, propName, reduction=True):
    """Return the dictionary proc['zoneName']."""
    d1 = Distributed.getPropertyDict(t, propName)
    #d1 = allgatherDict(d1)
    if reduction:
        d = KCOMM.allgather(d1)
        items = []
        for i in d: items += i.items()
        propDict = dict(items)
    else: propDict = d1
    return propDict

#==============================================================================
# IN: z: zone node
#==============================================================================
def getProperty(z, propName):
    """Return the property value for zone z."""
    return Distributed.getProperty(z, propName)
