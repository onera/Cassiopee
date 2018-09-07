# Interface de Cassiopee pour mpi
import PyTree as C
import Internal
import Distributed

# Acces a Distributed
from Distributed import readZones, writeZones, convert2PartialTree, convert2SkeletonTree, readNodesFromPaths, readPyTreeFromPaths, writeNodesFromPaths

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
    if graph.has_key(rank):
        g = graph[rank] # graph du proc courant
        for oppNode in g.keys():
            # Envoie les datas necessaires au noeud oppose
            #print '%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode])
            s = KCOMM.isend(datas[oppNode], dest=oppNode)
            reqs.append(s)
    rcvDatas={}
    for node in graph.keys():
        #print rank, graph[node].keys()
        if rank in graph[node].keys():
            #print '%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank])
            rcvDatas[node] = KCOMM.recv(source=node)
    MPI.Request.Waitall(reqs)
    return rcvDatas

def sendRecv2(datas, graph):
    if graph == {}: return {}
    reqs = []
    if graph.has_key(rank):
        g = graph[rank] # graph du proc courant
        for oppNode in g.keys():
            # Envoie les datas necessaires au noeud oppose
            print '%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode])
            s = KCOMM.isend(datas[oppNode], dest=oppNode)
            reqs.append(s)
    rcvDatas={}
    for node in graph.keys():
        #print rank, graph[node].keys()
        if rank in graph[node].keys():
            print '%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank])
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
def convertPyTree2File(t, fileName, format=None):
    """Write a skeleton or partial tree."""
    tp = convert2PartialTree(t)
    tp = C.deleteEmptyZones(tp)
    nzones = len(Internal.getZones(tp))
    if rank == 0:
        if nzones > 0:
            C.convertPyTree2File(tp, fileName, format=format); go = 1
        else: go = 0
        if size > 1: KCOMM.send(go, dest=1)
    else:
        go = KCOMM.recv(source=rank-1)
        if go == 1:
            Distributed.writeZones(tp, fileName, format=format, proc=rank)
        else:
            if nzones > 0:
                C.convertPyTree2File(tp, fileName, format=format); go = 1
        if rank < size-1: KCOMM.send(go, dest=rank+1)

#==============================================================================
# Execute sequentiellement F sur tous les procs
#==============================================================================
def seq(F):
    if rank == 0:
        F()
        KCOMM.send(1, dest=1)
    else:
        go = KCOMM.recv(source=rank-1)
        F()
        if rank < size-1: KCOMM.send(rank+1, dest=rank+1)
        
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
                zbb = G.BB(z, method, weighting); zb.append(zbb)
        
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
# IN: procDict: le procDict de tous l'arbre, si None, on le calcul.
# les procs. A utiliser si le graph est construit a partir d'arbres partiels
# ou pour type='proc'
# IN: intersectionsDict: dictionnaire d'intersections. Consulter la doc de
# computeGraph dans Distributed.py pour plus de details.
#==============================================================================
def computeGraph(t, type='bbox', t2=None, procDict=None, reduction=True, 
                 intersectionsDict=None):
    """Return the communication graph for different block relation types."""
    if not procDict: procDict = getProcDict(t)
    graph = Distributed.computeGraph(t, type, t2, procDict, rank, 
                                     intersectionsDict)

    if reduction:
        # Assure que le graph est le meme pour tous les processeurs
        g = KCOMM.allgather(graph)
        items = []
        for i in g: items += i.items()
        graph = dict(items)
        # a tester
        #d = {}
        #for i in g: d.update(i)

    return graph

#==============================================================================
# Recupere les zones specifiees dans le graph, les ajoute a l'arbre local t
#==============================================================================
def addXZones(t, graph):
    """Add zones specified in graph on current proc."""
    tp = Internal.copyRef(t)
    _addXZones(tp, graph)
    return tp

def _addXZones(t, graph):
    if graph == {}: return t
    reqs = []
    if graph.has_key(rank):
        g = graph[rank] # graph du proc courant
        for oppNode in g.keys():
            # Envoie les zones necessaires au noeud oppose
            #print '%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode])
            names = g[oppNode]
            data = [] # data est une liste de zones
            for n in names:
                zone = Internal.getNodesFromName2(t, n)[0]
                data.append(zone)
            s = KCOMM.isend(data, dest=oppNode)
            reqs.append(s)

    # Reception
    for node in graph.keys():
        #print rank, graph[node].keys()
        if rank in graph[node].keys():
            #print '%d: On doit recevoir de %d: %s'%(rank,node,graph[node][rank])
            data = KCOMM.recv(source=node)
            
            for z in data: # data est une liste de zones
                #print '%d: recoit la zone %s.'%(rank,z[0])
                # tag z
                Internal.createChild(z, 'XZone', 'UserDefinedData_t') 
                # Existe deja? 
                zone = Internal.getNodesFromName2(t, z[0])
                if len(zone) > 0: # replace
                    (p, c) = Internal.getParentOfNode(t, zone[0])
                    p[2][c] = z
                else: # append to first base
                    bases = Internal.getBases(t)
                    bases[0][2] += [z]
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
        tag = Internal.getNodesFromName1(z, 'XZone')
        if len(tag) > 0:
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
    Distributed._setProc(z, rank)
    return None

