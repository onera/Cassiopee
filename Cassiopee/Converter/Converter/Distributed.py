# Ce module permet de manipuler les SkeletonTree, les PartialTree et les
# BBoxTree
# - Skeleton tree: arbre complet avec les numpy float de taille > 100 non
# charge et valant None
# - Partial tree: arbre partiel ne contenant que les zones chargees
# - BBox tree: arbre identique a un pyTree mais ou les zones sont les
# BBox des zones du pyTree.

from . import Converter
from . import Internal
from . import PyTree
import numpy

#==============================================================================
# Lit un arbre squelette
# Warning: pour l'instant limite a hdf et adf
#==============================================================================
def convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5,
                             maxDepth=-1, dataShape=None, links=None):
    """Read a file and return a skeleton tree."""
    return PyTree.convertFile2PyTree(
        fileName, format, skeletonData=[maxFloatSize,maxDepth],
        dataShape=dataShape, links=links)

#==============================================================================
# Lit seulement un noeud de l'arbre ou ses enfants (suivant maxDepth)
#==============================================================================
def readNodesFromPaths(fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1,
                       dataShape=None, skipTypes=None, com=None, readIntMode=0):
    """Read nodes from file given their paths."""
    if format is None: format = Converter.convertExt2Format__(fileName)
    if not isinstance(paths, list): p = [paths]
    else: p = paths
    p = fixPaths__(p)
    if skipTypes is not None and isinstance(skipTypes, str): skipTypes = [skipTypes]
    if skipTypes is not None and isinstance(skipTypes, (str, tuple)): skipTypes = [skipTypes]

    ret = Converter.converter.readPyTreeFromPaths(fileName, p, format, maxFloatSize, maxDepth, readIntMode, dataShape, skipTypes, com)
    if not isinstance(paths, list): return ret[0]
    else: return ret

#==============================================================================
# Lit un noeud de l'arbre ou ses enfants (suivant maxDepth)
# et complete t
#==============================================================================
def readPyTreeFromPaths(t, fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1, setOnlyValue=True,
                        dataShape=None, skipTypes=None, com=None, readIntMode=0):
    """Read nodes from file given their path and complete t."""
    tp = Internal.copyRef(t)
    _readPyTreeFromPaths(tp, fileName, paths, format, maxFloatSize, maxDepth, setOnlyValue, dataShape, skipTypes, com, readIntMode)
    return tp

def _readPyTreeFromPaths(t, fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1, setOnlyValue=True,
                         dataShape=None, skipTypes=None, com=None, readIntMode=0):
    nodes = readNodesFromPaths(fileName, paths, format, maxFloatSize, maxDepth, dataShape, skipTypes, com, readIntMode)
    if not isinstance(paths, list): nodes = [nodes]; paths = [paths]
    if len(paths) != len(nodes):
        print("Warning: readPyTreeFromPaths: some paths can not be loaded. Nothing added to pyTree.")
        return None
    c = 0
    for p in paths:
        n = nodes[c]
        place = Internal.getNodeFromPath(t, p)
        if place is not None:
            place[0] = n[0]
            place[1] = n[1]
            if maxDepth == -1: place[2] = n[2]
            elif not setOnlyValue: place[2] = n[2]
            place[3] = n[3]
        else:
            parent = Internal.getPathAncestor(p)
            place = Internal.getNodeFromPath(t, parent)
            if place is not None:
                place[2].append(n)
            else:
                print("Warning: readPyTreeFromPaths: can not add node %s to t."%n[0])
        c += 1
    return None

#==============================================================================
# Ecrit seulement un noeud de l'arbre ou ses enfants (suivant maxDepth)
#==============================================================================
def writeNodesFromPaths(fileName, paths, nodes, format=None, maxDepth=-1, mode=0, isize=8, rsize=8):
    """Write nodes to file given their paths."""
    if format is None: format = Converter.convertExt2Format__(fileName)
    if not isinstance(paths, list): p = [paths]; n = [nodes]
    else: p = paths; n = nodes
    p = fixPaths__(p)
    Converter.converter.writePyTreePaths(fileName, n, p, format, maxDepth, mode, None, isize, rsize)
    return None

def writePyTreeFromPaths(t, fileName, paths, format=None, maxDepth=-1):
    """Write some nodes of the pyTree given their paths."""
    nodes = []
    if not isinstance(paths, list): p = [paths]
    else: p = paths
    paths = fixPaths__(p)
    opaths = []
    for p in paths:
        n = Internal.getNodeFromPath(t, p)
        if n is not None:
            nodes.append(n)
            opaths.append(p.rsplit('/',1)[0])
        else: print('Warning: write: path %s not found.'%p)
    writeNodesFromPaths(fileName, opaths, nodes, format, maxDepth, mode=0)
    return None

#========================================================================
# delete nodes in files (and all subsequent nodes) from paths
#========================================================================
def deletePaths(fileName, paths, format=None):
    """Delete nodes in file given their paths."""
    if format is None: format = Converter.convertExt2Format__(fileName)
    if format == 'bin_cgns' or format == 'unknown': format = Converter.checkFileType(fileName)
    if not isinstance(paths, list): p = [paths]
    else: p = paths
    p = fixPaths__(p)
    Converter.converter.deletePyTreePaths(fileName, p, format)
    return None

#=========================================================================
# Fix paths pour writePaths
# 1. Ne doit pas contenir CGNSTree en premier
# 2. Ne doit pas contenir nodeName en dernier
#=========================================================================
def fixPath__(path, nodeName=None):
    p = path; r = p
    if p[0:10] == '/CGNSTree/': r = p[9:]
    elif p[0:9] == 'CGNSTree/': r = p[8:]
    if nodeName is not None and nodeName == Internal.getPathLeaf(p): r = Internal.getPathAncestor(r)
    return r

def fixPaths__(paths, nodes=None):
    l = len(paths)
    out = []
    if nodes is not None:
        for i in range(l):
            r = fixPath__(paths[i], nodes[i][0])
            out.append(r)
    else:
        for i in range(l):
            r = fixPath__(paths[i])
            out.append(r)
    return out

#==============================================================================
# Determine si une zone est une zone squelette
# if ntype=0, se base sur les coord + FlowSolution
# if ntype=1, se base sur les coord uniquement
#==============================================================================
def isZoneSkeleton__(z, ntype=0):
    if len(z[2]) == 0: return True
    cx = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
    if cx is not None:
        t1 = Internal.getNodesFromType1(cx, 'DataArray_t')
        for d in t1:
            if d[1] is None: return True
        if len(cx[2]) == 0: return True
        return False
    if ntype == 0:
        cx = Internal.getNodesFromType1(z, 'FlowSolution_t')
        for x in cx:
            t1 = Internal.getNodesFromType1(x, 'DataArray_t')
            for d in t1:
                if d[1] is None: return True
    return False

#==============================================================================
# Converti un arbre en arbre squelette
# i.e supprime les noeuds DataArray_t et les remplace par None
#==============================================================================
def convert2SkeletonTree(t, maxSize=6):
    """Convert a tree to a skeleton tree."""
    tp = Internal.copyRef(t)
    _convert2SkeletonTree(tp, maxSize=maxSize)
    return tp

def _convert2SkeletonTree(t, maxSize=6):
    """Convert a tree to a skeleton tree."""
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromType(z, 'DataArray_t')
        for n in nodes:
            pt = n[1]
            if pt is not None and pt.size > maxSize: n[1] = None
    return None

#==============================================================================
# Converti un arbre squelette charge en arbre partiel
# rank=-1: enleve les zones squelettes (coord=None)
# rank>=0: enleve les zones de proc != rank
# Supprime aussi les champs squelettes
#==============================================================================
def convert2PartialTree(t, rank=-1):
    """Convert a tree to a partial tree."""
    tp = Internal.copyRef(t)
    _convert2PartialTree(tp, rank)
    return tp

def _convert2PartialTree(t, rank=-1):
    """Convert a tree to a partial tree."""
    zones = Internal.getZones(t)
    for z in zones:
        crit = False
        if rank == -1: crit = isZoneSkeleton__(z, ntype=1)
        else:
            proc = Internal.getNodeFromName2(z, 'proc')
            if proc is not None:
                proc = Internal.getValue(proc)
                crit = (rank != proc)
        if crit:
            # Supprime entierement la zone
            (p, c) = Internal.getParentOfNode(t, z)
            if Internal.isStdNode(t) == 0 and id(p) == id(t): del p[c]
            else: del p[2][c]
        else: # enleve les champs squelettes
            cx = Internal.getNodesFromType1(z, 'FlowSolution_t')
            for x in cx:
                t1 = Internal.getNodesFromType1(x, 'DataArray_t')
                for d in t1:
                    if d[1] is None:
                        p = Internal.getPath(x, d)
                        Internal._rmNodeByPath(z, p)
    return None

#==============================================================================
# Lit des zones et les remplace dans l'arbre
# IN: tree: arbre
# IN: rank: zones a charger
# IN: ou zoneNames: une liste des noms des zones a charger (['Base/Zone'])
# Warning: limite a adf et hdf
# Note: si un noeud Proc existe dans la zone remplacee, il est recopie
# dans la nouvelle zone
#==============================================================================
def readZones(t, fileName, format=None, rank=None, zoneNames=None):
    """Read some zones in a skeleton tree (by rank or name)."""
    tp = Internal.copyRef(t)
    _readZones(tp, fileName, format, rank, zoneNames)
    return tp

def _readZones(t, fileName, format=None, rank=None, zoneNames=None):
    """Read some zones in a skeleton tree (by rank or name)."""
    if zoneNames is None and rank is None: return None
    bases = Internal.getBases(t)
    if rank is not None: # load zones by rank
        # Chemins des zones a remplacer
        paths = []
        for b in bases:
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                nproc = Internal.getNodeFromName2(z, 'proc')
                if nproc is not None:
                    nproc = Internal.getValue(nproc)
                    if nproc == rank: paths.append('/'+b[0]+'/'+z[0])

    else: # by zone names
        paths = zoneNames[:]
        if isinstance(paths, str): paths = [paths]
        for c in range(len(paths)):
            paths[c] = paths[c].replace('CGNSTree', '')
            if paths[c][0] != '/': paths[c] = '/'+paths[c]

    #print('Reading '+fileName+' '+str(paths)+'...'),
    print('Reading %s [%d zones]...'%(fileName,len(paths))),
    if format is None: format = Converter.convertExt2Format__(fileName)
    if format == 'bin_cgns' or format == 'unknown': format = Converter.checkFileType(fileName)

    loadedZones = Converter.converter.readPyTreeFromPaths(fileName, paths, format, -1, -1, 0, None, None, None)

    import Compressor.PyTree as Compressor
    for z in loadedZones:
        Compressor._uncompressCartesian(z)
        Compressor._uncompressAll(z)

    # Replace/add now loaded zones
    m = 0
    for p in paths:
        z = Internal.getNodeFromPath(t, p)
        if z is not None:
            if rank is not None: # recopie les solveurs data
                nproc = Internal.getNodeFromName2(z, 'proc')
                if nproc is not None:
                    nproc = Internal.getValue(nproc)
                else: nproc = 0
            (p, c) = Internal.getParentOfNode(t, z)
            if Internal.isStdNode(t) == 0 and id(t) == id(p):
                p[c] = loadedZones[m]; zr = p[c]
            else: p[2][c] = loadedZones[m]; zr = p[2][c]
            if rank is not None:
                node = Internal.getNodeFromName1(zr, '.Solver#Param')
                if node is not None: param = node
                else:
                    param = ['.Solver#Param', None, [], 'UserDefinedData_t']
                    zr[2].append(param)
                v = numpy.zeros((1,1), Internal.E_NpyInt); v[0,0] = nproc
                node = Internal.getNodeFromName1(param, 'proc')
                if node is not None:
                    node[1] = v
                else:
                    a = ['proc', v, [], 'DataArray_t']
                    param[2].append(a)
            m += 1

    print('done.')
    return None

#==============================================================================
# Ecrit des zones dans un fichier deja cree
# si zoneNames != None: ecrit les zones specifiees
# si proc == rank: ecrit les zones dont le noeud proc correspond
# si proc == -1: ecrit tout les zones de t
# Warning: limite a adf et hdf
#==============================================================================
def writeZones(t, fileName, format=None, proc=None, zoneNames=None, links=None, isize=8, rsize=8):
    """Write some zones in an existing file (adf or hdf)."""
    if zoneNames is None and proc is None: return None
    tp, ntype = Internal.node2PyTree(t)
    bases = Internal.getBases(tp)
    if proc is not None: # write zones by proc
        # Chemins des noeuds a remplacer (Zone et IntegralData)
        paths = []; nodes = []
        for b in bases:
            zones = Internal.getNodesFromType1(b, 'Zone_t') + Internal.getNodesFromType1(b, 'IntegralData_t')
            for z in zones:
                if proc == -1:
                    paths.append('/%s'%b[0]); nodes.append(z)
                else:
                    nproc = Internal.getNodeFromName2(z, 'proc')
                    if nproc is not None:
                        nproc = Internal.getValue(nproc)
                        if nproc == proc:
                            paths.append('/%s'%b[0]); nodes.append(z)
                    else: # write nevertheless
                        paths.append('/%s'%b[0]); nodes.append(z)
    else: # by zone names
        paths = zoneNames[:]
        if isinstance(paths, str): paths = [paths]
        nodes = []
        for p in paths:
            n = Internal.getNodeFromPath(tp, p)
            nodes.append(n)
        for c in range(len(paths)):
            paths[c] = paths[c].replace('CGNSTree','')
            if paths[c][0] != '/': paths[c] = '/'+paths[c]
            paths[c] = Internal.getPathAncestor(paths[c])

    print('Writing %s [%d zones]...'%(fileName,len(paths))),
    if format is None: format = Converter.convertExt2Format__(fileName)
    Converter.converter.writePyTreePaths(fileName, nodes, paths, format, -1, 0, links, isize, rsize)
    print('done.')
    return None

#==============================================================================
# zones est une liste de zones, mais triees par procs
# Ex: [ [zonesDeProc0], [zonesDeProc1], ...]
# IN: liste des zones a remplacer dans t
# IN: size: nbre de processeurs
# Met les zones dans l'arbre par identification des noms des zones
# Si le nom des zones est base/zone, essaie de retrouver la base
#==============================================================================
def setZonesInTree(t, zones):
    tp = Internal.copyRef(t)
    _setZonesInTree(tp, zones)
    return tp

def _setZonesInTree(t, zones):
    size = len(zones)
    for i in range(size):
        for j in range(len(zones[i])):
            zone = zones[i][j]
            zoneName = zone[0]
            spl = zoneName.split('/', 1)
            if len(spl) >= 2: # baseName/zoneName
                baseName = spl[0]; zoneName = spl[1]
                zone[0] = zoneName
                base = Internal.getNodeFromName1(t, baseName)
                if base is not None:
                    z = Internal.getNodeFromName1(base, zoneName)
                    if z is not None:
                        (p, nb) = Internal.getParentOfNode(base, z)
                        p[2][nb] = zone
                    else:
                        base[2].append(zone)
                else:
                    base = Internal.newCGNSBase(baseName, parent=t)
                    base[2].append(zone)
            else:
                z = Internal.getNodeFromName2(t, zoneName)
                if z is not None:
                    (p, nb) = Internal.getParentOfNode(t, z)
                    p[2][nb] = zone
                else:
                    # append to first base
                    bases = Internal.getBases(t)
                    bases[0][2].append(zone)
    return None

#==============================================================================
# update le graph si zoneName sur proc est a envoyer a popp
#==============================================================================
def updateGraph__(graph, proc, popp, zoneName):
    if popp != proc:
        if proc not in graph: graph[proc] = {popp:[zoneName]}
        else:
            if popp not in graph[proc]: graph[proc][popp] = [zoneName]
            else:
                if zoneName not in graph[proc][popp]:
                    graph[proc][popp].append(zoneName)
    return graph

#==============================================================================
# merge deux graph
#==============================================================================
def mergeGraph(graph1, graph2):
    import copy
    graph = copy.deepcopy(graph1)
    for proc in graph2.keys():
        if proc in graph.keys():
            for popp in graph2[proc]:
                if popp in graph[proc]:
                    for zname in graph2[proc][popp]:
                        if zname not in graph[proc][popp]:  graph[proc][popp].append(zname)
                else:
                    graph[proc][popp] =  graph2[proc][popp]
        else:
            graph[proc] =  graph2[proc]

    return graph
#==============================================================================

def updateGraphSet__(graph, proc, popp, zoneName):
    if popp != proc:
        if proc not in graph: graph[proc] = {popp:set()}
        graph[proc][popp].add(zoneName)
    return graph

#==============================================================================
# Retourne le proc de z si z est local ou dans le procDict si il existe
#==============================================================================
def getProcLocal__(z, procDict=None):
    if procDict is not None: return procDict[z[0]]
    else:
        proc = Internal.getNodeFromName2(z, 'proc')
        if proc is not None: proc = Internal.getValue(proc)
        else: proc = -1
        return proc

#==============================================================================
# Retourne le proc de z si z est global ou dans le procDict si il existe
#==============================================================================
def getProcGlobal__(zoneName, t, procDict=None):
    if procDict is not None: return procDict[zoneName]
    else:
        z = Internal.getNodeFromName2(t, zoneName)
        proc = Internal.getNodeFromName2(z, 'proc')
        if proc is not None: proc = Internal.getValue(proc)
        else: proc = -1
        return proc

#==============================================================================
# Calcule le graph
# IN: t: arbre contenant des noeuds 'proc' entierement charge
# IN: type: type de graph voulu
# type='bbox' si intersection de bbox des zones (full/bbox)
# type='bbox2' si intersection de bbox et pas sur la meme base (full/bbox)
# type='bbox3' si intersection de bbox avec t2 (t:full/bbox et t2:full/bbox)
# type='match' si match entre zones (full/skel/load skel/partial[+procDict])
# type='nmatch' si nearmatch entre zones (full/skel/load skel/partial[+procDict])
# type='ID' si donnees d'interpolation entre zones (full/skel/load skel/partial+procDict)
# type='IBCD' si donnees IBC entres zones (full/skel/load skel/partial+procDict)
# type='ALLD' si toutes donnees (Interp+IBC) (full/skel/load skel/partial+procDict)
# type='proc' si la zone a un noeud proc different du proc courant (full/skel/load skel/partial)
# type='POST': t donor tree and t2 receptor tree. To interpolate data from donor to receptor where interpolation
#              data is stored in t in ID* subregions nodes. Requires procDict and procDict2
# intersectionsDict: dictionnaire Python contenant la liste de zones intersectantes,
# comme produit par "X.getIntersectingDomains". Attention, si type='bbox3',
# l'utilisateur doit fournir l'arbre t2 et intersectionDict doit decrire les
# intersections entre t et t2, comme produit par "X.getIntersectingDomains(t,t2)".
# exploc: True si explicite local
# OUT: graph: dictionnaire contenant des informations d'envoie
# des zones entre processeurs
# graph est construit de telle sorte que:
# graph[P0][P1] renvoie la liste des zones du processeur P0
# "connectees" avec au moins une zone du processeur P1
# Ex: addXZones: envoie les zones graph[rank][opp] au proc opp,
# attend ensuite les zones graph[opp][rank] pour tout opp.
#==============================================================================
def computeGraph(t, type='bbox', t2=None, procDict=None, rank=0,
                 intersectionsDict=None, exploc=False, procDict2=None,
                 it=None, reduction=True):
    """Return the communication graph for different block relation types."""
    zones = Internal.getZones(t)
    graph = {}
    if type == 'bbox': # zone de t de P0 intersectent une zone de t de P1
        if not intersectionsDict:
            try: import Connector.PyTree as X
            except: raise ImportError("computeGraph: requires Connector module.")
            intersectionsDict = X.getIntersectingDomains(t)
        for z in zones:
            proc = getProcLocal__(z, procDict)
            for z2 in zones:
                if z2[0] in intersectionsDict and z[0] in intersectionsDict[z2[0]]:
                    popp = getProcGlobal__(z2[0], t, procDict)
                    updateGraph__(graph, proc, popp, z[0])
        #import Connector.PyTree as X
        #for z in zones:
        #    proc = getProcLocal__(z, procDict)
        #    doms = X.getBBIntersectingDomains(t, z, tol=1.e-12)
        #    for d in doms:
        #        popp = getProcGlobal__(d[0], t, procDict)
        #        updateGraph__(graph, proc, popp, z[0])

    elif type == 'bbox2': # zone de t de P0 qui intersecte une zone de t de P1 mais qui n'est pas sur la meme base
        if not intersectionsDict:
            try: import Connector.PyTree as X
            except: raise ImportError("computeGraph: requires Connector module.")
            intersectionsDict = X.getIntersectingDomains(t)
        for z in zones:
            (p, c) = Internal.getParentOfNode(t, z)
            base = p[0]
            proc = getProcLocal__(z, procDict)
            for z2 in zones:
                if z2[0] in intersectionsDict and z[0] in intersectionsDict[z2[0]]:
                    (p, c) = Internal.getParentOfNode(t, z2)
                    based = p[0]
                    popp = getProcGlobal__(z2[0], t, procDict)
                    if popp != proc and base != based:
                        if proc not in graph: graph[proc] = {popp:[z[0]]}
                        else:
                            if popp not in graph[proc]: graph[proc][popp] = [z[0]]
                            else:
                                if z[0] not in graph[proc][popp]:
                                    graph[proc][popp].append(z[0])
        # import Connector.PyTree as X
        # for z in zones:
        #     (p, c) = Internal.getParentOfNode(t, z)
        #     base = p[0]
        #     proc = getProcLocal__(z, procDict)
        #     doms = X.getBBIntersectingDomains(t, z, tol=1.e-12)
        #     for d in doms:
        #         (p, c) = Internal.getParentOfNode(t, d)
        #         based = p[0]
        #         popp = getProcGlobal__(d[0], t, procDict)
        #         if (popp != proc and base != based):
        #             if proc not in graph: graph[proc] = {popp:[z[0]]}
        #             else:
        #                 if popp not in graph[proc]: graph[proc][popp] = [z[0]]
        #                 else:
        #                     if z[0] not in graph[proc][popp]:
        #                         graph[proc][popp].append(z[0])

    elif type == 'bbox3': # zone de t sur P0 qui intersecte une zone de t2 sur P1
        #if not t2: raise ValueError("computeGraph: type bbox3 requires a t2.")
        zones2 = Internal.getZones(t2)
        if not intersectionsDict:
            try: import Connector.PyTree as X
            except: raise ImportError("computeGraph: requires Connector module.")
            intersectionsDict = X.getIntersectingDomains(t, t2)
        for z in zones:
            proc = getProcLocal__(z, procDict)
            for z2 in zones2:
                if z[0] in intersectionsDict and z2[0] in intersectionsDict[z[0]]:
                    popp = getProcGlobal__(z2[0], t2, procDict2)
                    updateGraph__(graph, proc, popp, z[0])
                    #updateGraphSet__(graph, proc, popp, z[0])
        #for p in graph:
        #  for opp in graph[p]: graph[p][opp] = list(graph[p][opp])
        #import Connector.PyTree as X
        #for z in zones:
        #    proc = getProcLocal__(z, procDict)
        #    doms = X.getBBIntersectingDomains(t2, z, tol=1.e-12)
        #    for d in doms:
        #        popp = getProcGlobal__(d[0], t2, None)
        #        updateGraph__(graph, proc, popp, z[0])

    elif type == 'POST':
        if procDict is None:
            raise ValueError("computeGraph: type=POST requires procDict for donor tree")
        if procDict2 is None:
            raise ValueError("computeGraph: type=POST requires procDict for receptor tree")

        for z in zones:
            proc = getProcLocal__(z,procDict)
            subRegions = []
            for zsr in Internal.getNodesFromType1(z,'ZoneSubRegion_t'):
                sname = Internal.getName(zsr)[0:2]
                if sname=='ID':
                    rcvname = Internal.getValue(zsr)
                    popp = getProcGlobal__(rcvname, t2, procDict2)
                    updateGraph__(graph, proc, popp, z[0])

    elif type == 'ID_Unsteady': # base sur les interpolations data
        graph_steady={}; graph_unsteady={}
        for z in zones:
            proc = getProcLocal__(z, procDict)
            subRegions = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
            for s in subRegions:
                sname = s[0][0:2]
                if sname=='ID':
                    if '#' in s[0]:
                        numero_iter = int( s[0].split('#')[1].split('_')[0] )
                        donor = Internal.getValue(s)
                        idn = Internal.getNodesFromName1(s, 'InterpolantsDonor')
                        if idn != []: # la subRegion decrit des interpolations
                            popp = getProcGlobal__(donor, t, procDict)
                            if numero_iter not in graph_unsteady: graph_unsteady[numero_iter] = {}
                            #if(numero_iter==55):  print("subregion", s[0], z[0], donor, proc, popp)
                            updateGraph__(graph_unsteady[numero_iter], proc, popp, z[0])
                    else:
                        donor = Internal.getValue(s)
                        idn = Internal.getNodesFromName1(s, 'InterpolantsDonor')
                        if idn != []: # la subRegion decrit des interpolations
                            popp = getProcGlobal__(donor, t, procDict)
                            updateGraph__(graph_steady, proc, popp, z[0])

        return graph_steady, graph_unsteady

    elif type == 'ID': # base sur les interpolations data
        if not exploc:
            for z in zones:
                proc = getProcLocal__(z, procDict)
                subRegions = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
                for s in subRegions:
                    sname = s[0][0:2]
                    if sname=='ID':
                        if '#' in s[0] and it is not None:
                            numero_iter = int( s[0].split('#')[1].split('_')[0] )
                            if numero_iter == it:
                                donor = Internal.getValue(s)
                                idn = Internal.getNodesFromName1(s, 'InterpolantsDonor')
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph, proc, popp, z[0])
                        else:
                            donor = Internal.getValue(s)
                            idn = Internal.getNodesFromName1(s, 'InterpolantsDonor')
                            if idn != []: # la subRegion decrit des interpolations
                                popp = getProcGlobal__(donor, t, procDict)
                                updateGraph__(graph, proc, popp, z[0])

        else:
            maxlevel=1
            for z in zones:
                subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
                for s in subRegions2:
                    levrcv_ = Internal.getNodesFromName1(s,'LevelZRcv')
                    levrcv  = int(levrcv_[0][1][0])
                    levdnr_ = Internal.getNodesFromName1(s,'LevelZDnr')
                    levdnr  = int(levdnr_[0][1][0])
                    maximum = max(levrcv,levdnr)
                    if maximum > maxlevel:maxlevel=maximum
            nssiter = 4*maxlevel

            list_graph_=[]
            for ssiter in range(3,4):#range(1,2*nssiter+1):
                graph_ = {}
                for z in zones:
                    proc = getProcLocal__(z, procDict)
                    subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
                    subRegions = []
                    for s in subRegions2:
                        sname = s[0][0:2]
                        if sname == 'ID': subRegions.append(s)
                    for s in subRegions:
                        donor = Internal.getValue(s)
                        levdnr_ = Internal.getNodesFromName1(s,'LevelZDnr')
                        levdnr  = int(levdnr_[0][1][0])
                        levrcv_ = Internal.getNodesFromName1(s,'LevelZRcv')
                        levrcv  = int(levrcv_[0][1][0])
                        idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                        cycl = nssiter//levdnr

                        if levdnr > levrcv and ssiter <= nssiter:
                            if ssiter%cycl==cycl-1 or ssiter%cycl==cycl//2 and (ssiter//cycl)%2==1:
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr < levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==1 or ssiter%cycl==cycl//4 or ssiter%cycl==cycl//2-1 or ssiter%cycl==cycl//2+1 or ssiter%cycl==cycl//2+cycl//4 or ssiter%cycl==cycl-1):
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==cycl//2-1 or (ssiter%cycl==cycl//2 and (ssiter//cycl)%2==0) or ssiter%cycl==cycl-1):
                                #print(donor, ssiter)
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter > nssiter:
                            ssiter_ = ssiter - nssiter
                            if ssiter_%cycl==cycl//2 and (ssiter_//cycl)%2==1:
                                #if ssiter_%2==0 and ssiter_%cycl==cycl/2 and (ssiter_//cycl)%2==1:
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                list_graph_.append(graph_)

    elif type == 'IBCD': # base sur les IBC data
        if not exploc:
            for z in zones:
                proc = getProcLocal__(z, procDict)
                subRegions2 = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
                subRegions = []
                for s in subRegions2:
                    sname = s[0][0:2]
                    if sname == 'IB': subRegions.append(s)
                for s in subRegions:
                    donor = Internal.getValue(s)
                    idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                    if idn != []: # la subRegion decrit des IBC
                        popp = getProcGlobal__(donor, t, procDict)
                        updateGraph__(graph, proc, popp, z[0])

        else:
            maxlevel=1
            for z in zones:
                subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
                for s in subRegions2:
                    levrcv_ = Internal.getNodesFromName1(s,'LevelZRcv')
                    levrcv  = int(levrcv_[0][1][0])
                    levdnr_ = Internal.getNodesFromName1(s,'LevelZDnr')
                    levdnr  = int(levdnr_[0][1][0])
                    maximum = max(levrcv,levdnr)
                    if maximum > maxlevel:maxlevel=maximum
            nssiter = 4*maxlevel

            list_graph_=[]
            for ssiter in range(1,2*nssiter+1):
                graph_={}
                for z in zones:
                    proc = getProcLocal__(z, procDict)
                    subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
                    subRegions = []
                    for s in subRegions2:
                        sname = s[0][0:2]
                        if sname == 'IB': subRegions.append(s)

                    for s in subRegions:
                        donor = Internal.getValue(s)
                        levdnr_ = Internal.getNodesFromName1(s,'LevelZDnr')
                        levdnr  = int(levdnr_[0][1][0])
                        levrcv_ = Internal.getNodesFromName1(s,'LevelZRcv')
                        levrcv  = int(levrcv_[0][1][0])
                        idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                        cycl = nssiter//levdnr
                        if levdnr > levrcv and ssiter <= nssiter:
                            if ssiter%cycl==cycl-1 or ssiter%cycl==cycl//2 and (ssiter//cycl)%2==1:
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr < levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==1 or ssiter%cycl==cycl//4 or ssiter%cycl==cycl//2-1 or ssiter%cycl==cycl//2+1 or ssiter%cycl==cycl//2+cycl//4 or ssiter%cycl==cycl-1):
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==cycl//2-1 or (ssiter%cycl==cycl//2 and (ssiter//cycl)%2==0) or ssiter%cycl==cycl-1):
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter > nssiter:
                            #if (ssiter%8==6):
                            ssiter_ = ssiter - nssiter
                            if ssiter_%2==0 and ssiter_%cycl==cycl//2 and (ssiter_//cycl)%2==1:
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                list_graph_.append(graph_)


    elif type == 'INV_IBCD': # base sur les IBC data inverse (from receiver)
        if not exploc:
            for z in zones:
                proc = getProcLocal__(z, procDict)
                subRegions2 = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
                subRegions = []
                for s in subRegions2:
                    sname = s[0][0:2]
                    if sname=='IB': subRegions.append(s)
                for s in subRegions:
                    donor = Internal.getValue(s)
                    idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                    if idn != []: # la subRegion decrit des IBC
                        popp = getProcGlobal__(donor, t, procDict)
                        updateGraph__(graph, popp, proc, donor)


    elif type == '2_IBCD': # base sur les IBC data
        if not exploc:
            for z in zones:
                proc = getProcLocal__(z, procDict)
                subRegions2 = Internal.getNodesFromType1(z, 'ZoneSubRegion_t')
                subRegions = []
                for s in subRegions2:
                    sname = s[0][0:2]
                    if sname=='2_': subRegions.append(s)
                for s in subRegions:
                    donor = Internal.getValue(s)
                    idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                    if idn != []: # la subRegion decrit des IBC
                        popp = getProcGlobal__(donor, t, procDict)
                        updateGraph__(graph, proc, popp, z[0])


        else:
            maxlevel=1
            for z in zones:
                subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
                for s in subRegions2:
                    levrcv_ = Internal.getNodesFromName1(s,'LevelZRcv')
                    levrcv  = int(levrcv_[0][1][0])
                    levdnr_ = Internal.getNodesFromName1(s,'LevelZDnr')
                    levdnr  = int(levdnr_[0][1][0])
                    maximum = max(levrcv,levdnr)
                    if maximum > maxlevel:maxlevel=maximum
            nssiter = 4*maxlevel

            list_graph_=[]
            for ssiter in range(1,2*nssiter+1):
                graph_={}
                for z in zones:
                    proc = getProcLocal__(z, procDict)
                    subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
                    subRegions = []
                    for s in subRegions2:
                        sname = s[0][0:2]
                        if sname=='2_': subRegions.append(s)

                    for s in subRegions:
                        donor = Internal.getValue(s)
                        levdnr_ = Internal.getNodesFromName1(s,'LevelZDnr')
                        levdnr  = int(levdnr_[0][1][0])
                        levrcv_ = Internal.getNodesFromName1(s,'LevelZRcv')
                        levrcv  = int(levrcv_[0][1][0])
                        idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                        cycl = nssiter//levdnr
                        if levdnr > levrcv and ssiter <= nssiter:
                            if ssiter%cycl==cycl-1 or ssiter%cycl==cycl//2 and (ssiter//cycl)%2==1:
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr < levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==1 or ssiter%cycl==cycl//4 or ssiter%cycl==cycl//2-1 or ssiter%cycl==cycl//2+1 or ssiter%cycl==cycl//2+cycl//4 or ssiter%cycl==cycl-1):
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==cycl//2-1 or (ssiter%cycl==cycl//2 and (ssiter//cycl)%2==0) or ssiter%cycl==cycl-1):
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter > nssiter:
                            #if (ssiter%8==6):
                            ssiter_ = ssiter - nssiter
                            if ssiter_%2==0 and ssiter_%cycl==cycl//2 and (ssiter_//cycl)%2==1:
                                if idn != []: # la subRegion decrit des interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                list_graph_.append(graph_)


    elif type == 'ALLD': # base sur les Interpolations+IBC data
        for z in zones:
            proc = getProcLocal__(z, procDict)
            subRegions = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
            for s in subRegions:
                donor = Internal.getValue(s)
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations/IBC
                    popp = getProcGlobal__(donor, t, procDict)
                    updateGraph__(graph, proc, popp, z[0])

    elif type == 'match': # base sur les raccords matchs
        for z in zones:
            proc = getProcLocal__(z, procDict)
            GC = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
            GC2 = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for c in GC2:
                r = Internal.getNodeFromType1(c, 'GridConnectivityType_t')
                if r is not None and Internal.getValue(r) == 'Abutting1to1':
                    GC.append(c)
            for c in GC:
                donor = Internal.getValue(c)
                popp = getProcGlobal__(donor, t, procDict)
                updateGraph__(graph, proc, popp, z[0])

    elif type == 'nmatch': # base sur les nearmatchs
        for z in zones:
            proc = getProcLocal__(z, procDict)
            GC = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            for c in GC:
                gctype = Internal.getNodeFromType(c, 'GridConnectivityType_t')
                if Internal.getValue(gctype)=='Abutting':
                    donor = Internal.getValue(c)
                    popp = getProcGlobal__(donor, t, procDict)
                    updateGraph__(graph, proc, popp, z[0])
    elif type == 'proc':
        for z in zones:
            if not isZoneSkeleton__(z):
                popp = getProcLocal__(z, procDict)
                proc = rank
                updateGraph__(graph, proc, popp, z[0])

    if not exploc: return graph
    else: return list_graph_

# Split graph in one zone graphs
def splitGraph(graph):
    graphs = []
    c = 0; goon = True
    while goon:
        g = {}
        goon = False
        for p1 in graph:
            g[p1] = {}
            for p2 in graph[p1]:
                l = graph[p1][p2]
                if len(l) > c: g[p1][p2] = [l[c]]; goon = True
                else: g[p1][p2] = []
        if goon: graphs.append(g)
        c += 1
    return graphs

#==============================================================================
# Retourne le dictionnaire proc['blocName']
# a partir d'un arbre distribue contenant des noeuds proc
#==============================================================================
def getProcDict(t):
    """Return the dictionary proc['zoneName']."""
    proc = {}
    zones = Internal.getZones(t)
    for z in zones:
        nproc = getProc(z)
        proc[z[0]] = nproc
    return proc

#==============================================================================
# getProc in zone (if exists), otherwise return -1
# IN: a: zone node
#==============================================================================
def getProc(a):
    """Return the proc of zone."""
    nproc = Internal.getNodeFromName2(a, 'proc')
    if nproc is not None: nproc = Internal.getValue(nproc)
    else: nproc = -1
    return nproc

#==============================================================================
# setProc in zone
# IN: z: zone node
# IN: nproc: proc to set (int)
#==============================================================================
def setProc(t, nproc):
    """Set the proc to a zone or a set of zones."""
    tp = Internal.copyRef(t)
    _setProc(tp, nproc)
    return tp

def _setProc(t, nproc):
    zones = Internal.getZones(t)
    for z in zones:
        node = Internal.getNodeFromName1(z, '.Solver#Param')
        if node is not None: param = node
        else:
            param = ['.Solver#Param', None, [], 'UserDefinedData_t']
            z[2].append(param)
        v = numpy.zeros((1,1), Internal.E_NpyInt); v[0,0] = nproc
        node = Internal.getNodeFromName1(param, 'proc')
        if node is not None: node[1] = v
        else:
            a = ['proc', v, [], 'DataArray_t']
            param[2].append(a)
    return None

#==============================================================================
# Retourne le dictionnaire prop['blocName']
# a partir d'un arbre distribue contenant des noeuds 'propname'
#==============================================================================
def getPropertyDict(t, propname):
    """Return the dictionary prop['zoneName']."""
    prop = {}
    zones = Internal.getZones(t)
    for z in zones:
        nprop = getProperty(z, propname)
        prop[z[0]] = nprop
    return prop

#==============================================================================
# get property in zone (if exists), otherwise return -1
# IN: a: zone node
#==============================================================================
def getProperty(a, propname):
    """Return the property of zone."""
    nprop = Internal.getNodeFromName2(a, propname)
    if nprop is not None: nprop = Internal.getValue(nprop)
    else: nprop = -1
    return nprop
