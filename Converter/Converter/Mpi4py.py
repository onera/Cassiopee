# Interface de Cassiopee pour mpi
from . import PyTree as C
from . import Internal
from . import Distributed
import Compressor.PyTree as Compressor

# Acces a Distributed
from .Distributed import readZones, writeZones, convert2PartialTree, convert2SkeletonTree, readNodesFromPaths, readPyTreeFromPaths, writeNodesFromPaths

__all__ = ['rank', 'size', 'KCOMM', 'setCommunicator', 'barrier', 'send', 'recv', 'sendRecv', 'sendRecv2', 
    'allgather', 'readZones', 'writeZones', 'convert2PartialTree', 'convert2SkeletonTree', 'convertFile2DistributedPyTree', 
    'readNodesFromPaths', 'readPyTreeFromPaths', 'writeNodesFromPaths',
    'allgatherTree', 'convertFile2SkeletonTree', 'convertFile2PyTree', 'convertPyTree2File', 'seq', 'print0', 'printA',
    'createBBoxTree', 'computeGraph', 'addXZones', '_addXZones', 'rmXZones', '_rmXZones', 'getProcDict', 
    'getProc', 'setProc', '_setProc']

from mpi4py import MPI

COMM_WORLD = MPI.COMM_WORLD
KCOMM = COMM_WORLD

rank = KCOMM.rank
size = KCOMM.size

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
            #print '%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode])
            s = KCOMM.isend(datas[oppNode], dest=oppNode)
            reqs.append(s)
    rcvDatas={}
    for node in graph:
        #print rank, graph[node].keys()
        if rank in graph[node]:
            #print '%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank])
            rcvDatas[node] = KCOMM.recv(source=node)
    MPI.Request.Waitall(reqs)
    return rcvDatas

def sendRecv2(datas, graph):
    if graph == {}: return {}
    reqs = []
    if rank in graph:
        g = graph[rank] # graph du proc courant
        for oppNode in g:
            # Envoie les datas necessaires au noeud oppose
            print('%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode]))
            s = KCOMM.isend(datas[oppNode], dest=oppNode)
            reqs.append(s)
    rcvDatas={}
    for node in graph:
        #print rank, graph[node].keys()
        if rank in graph[node]:
            print('%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank]))
            rcvDatas[node] = KCOMM.recv(source=node)
    MPI.Request.Waitall(reqs)
    return rcvDatas

#==============================================================================
# allgather
#==============================================================================
def allgather(data):
    return KCOMM.allgather(data)

def allgatherTree(t):
    """Gather a distributed tree on all processors."""
    d = KCOMM.allgather(t)
    return Internal.merge(d)
    
#==============================================================================
# Lecture du squelette d'un arbre dans un fichier
# Lecture proc 0 + bcast
#==============================================================================
def convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5,
                             maxDepth=-1):
    """Read a file and return a skeleton tree."""
    if rank == 0: t = Distributed.convertFile2SkeletonTree(fileName, format, maxFloatSize, maxDepth)
    else: t = None
    t = KCOMM.bcast(t)
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
# Lecture proc 0 + bcast
#==============================================================================
def convertFile2PyTree(fileName, format=None):
    """Read a file and return a full tree."""
    if rank == 0: t = C.convertFile2PyTree(fileName, format)
    else: t = None
    t = KCOMM.bcast(t)
    return t

#==============================================================================
# Ecriture sequentielle
# Avec recuperation de toutes les zones
# Attention: les zones doivent avoir un procNode.
#==============================================================================
def convertPyTree2File(t, fileName, format=None, links=[]):
    """Write a skeleton or partial tree."""
    tp = convert2PartialTree(t)
    tp = C.deleteEmptyZones(tp)
    nzones = len(Internal.getZones(tp))
    if rank == 0:
        if nzones > 0:
            C.convertPyTree2File(tp, fileName, format=format, links=links); go = 1
        else: go = 0
        if size > 1: KCOMM.send(go, dest=1)
    else:
        go = KCOMM.recv(source=rank-1)
        if go == 1:
            Distributed.writeZones(tp, fileName, format=format, proc=rank, links=links)
        else:
            if nzones > 0:
                C.convertPyTree2File(tp, fileName, format=format, links=links); go = 1
        if rank < size-1: KCOMM.send(go, dest=rank+1)

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
                # Clean up (zoneSubRegion)
                for c, l in enumerate(zbb[2]):
                    if l[3] == 'ZoneSubRegion_t': del zbb[2][c]
                zb.append(zbb)
        
    # Echanges des zones locales de bounding box
    # (allgather serialise automatiquement les donnees)
    zones = KCOMM.allgather(zb)

    # On les remplace dans l'arbre (par noms)
    tb = Distributed.setZonesInTree(t, zones)
    return tb

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
                 intersectionsDict=None, exploc=False):
    """Return the communication graph for different block relation types."""
    if not procDict: procDict = getProcDict(t)

    graph = Distributed.computeGraph(t, type, t2, procDict, rank, 
                                     intersectionsDict, exploc)

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

#==============================================================================
# Recupere les zones specifiees dans le graph, les ajoute a l'arbre local t
#==============================================================================
def addXZones(t, graph, variables=None, cartesian=False):
    """Add zones specified in graph on current proc."""
    tp = Internal.copyRef(t)
    _addXZones(tp, graph, variables, cartesian)
    return tp

def _addXZones(t, graph, variables=None, cartesian=False):
    if graph == {}: return t
    reqs = []
    if rank in graph:
        g = graph[rank] # graph du proc courant
        for oppNode in g:
            # Envoie les zones necessaires au noeud oppose
            #print '%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode])
            names = g[oppNode]
            data = [] # data est une liste de zones
            for n in names:
                zone = Internal.getNodeFromName2(t, n)
                if variables is not None:
                    v = C.getVarNames(zone, excludeXYZ=True)[0]
                    for i in variables: v.remove(i)
                    zonep = C.rmVars(zone, v)
                    if cartesian: Compressor._compressCartesian(zonep)
                    data.append(zonep)
                else:
                    if cartesian: 
                        zonep = Compressor._compressCartesian(zonep)
                        data.append(zonep)
                    else: data.append(zone)
            s = KCOMM.isend(data, dest=oppNode)
            reqs.append(s)

    # Reception
    for node in graph:
        #print rank, graph[node].keys()
        if rank in graph[node]:
            #print '%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank])
            data = KCOMM.recv(source=node)
            
            for z in data: # data est une liste de zones
                if cartesian: Compressor._uncompressCartesian(z)
                #print '%d: recoit la zone %s.'%(rank,z[0])
                # tag z
                Internal.createChild(z, 'XZone', 'UserDefinedData_t') 
                # Existe deja? 
                zone = Internal.getNodeFromName2(t, z[0])
                if zone is not None: # replace
                    (p, c) = Internal.getParentOfNode(t, zone)
                    p[2][c] = z
                else: # append to first base
                    bases = Internal.getBases(t)
                    bases[0][2].append(z)
    MPI.Request.Waitall(reqs)
    return t

#==============================================================================
# Supprime les zones ajoutees par addXZones
#==============================================================================
def rmXZones(t):
    """Remove zones added by addXZones."""
    tp = Internal.copyRef(t)
    _rmXZones(tp)
    return tp

def _rmXZones(t):
    zones = Internal.getZones(t)
    for z in zones:
        tag = Internal.getNodeFromName1(z, 'XZone')
        if tag is not None:
            (p, c) = Internal.getParentOfNode(t, z)
            del p[2][c]
    return None

#==============================================================================
# Retourne le dictionnaire proc pour chaque zoneName
# IN: t: full/S/LS/Partial
#==============================================================================
def getProcDict(t):
    """Return the dictionary proc['zoneName']."""
    d1 = Distributed.getProcDict(t)
    d = KCOMM.allgather(d1)
    items = []
    for i in d: items += i.items()
    procDict = dict(items)
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
