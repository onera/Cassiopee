# This module allows manipulating SkeletonTree, PartialTree and
# BBoxTree
# - Skeleton tree: complete tree with numpy floats of size > 100 not
# loaded and equal to None
# - Partial tree: partial tree containing only loaded zones
# - BBox tree: tree identical to a pyTree but where zones are the
# BBox of the pyTree zones.

from . import Converter
from . import Internal
from . import PyTree
import numpy

#==============================================================================
# Reads a skeleton tree
# Warning: currently limited to hdf and adf
#==============================================================================
def convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5,
                             maxDepth=-1, dataShape=None, links=None):
    """Read a file and return a skeleton tree."""
    return PyTree.convertFile2PyTree(
        fileName, format, skeletonData=[maxFloatSize,maxDepth],
        dataShape=dataShape, links=links)

#==============================================================================
# Reads only one node of the tree or its children (depending on maxDepth)
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
# Reads a node of the tree or its children (depending on maxDepth)
# and completes t
#==============================================================================
def readPyTreeFromPaths(t, fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1, setOnlyValue=True,
                        dataShape=None, skipTypes=None, com=None, readIntMode=0):
    """Read nodes from file given their path and complete t."""
    tp = Internal.copyRef(t)
    _readPyTreeFromPaths(tp, fileName, paths, format, maxFloatSize, maxDepth, setOnlyValue, dataShape, skipTypes, com, readIntMode)
    return tp

def _readPyTreeFromPaths(t, fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1, setOnlyValue=True,
                         dataShape=None, skipTypes=None, com=None, readIntMode=0):
    """Read nodes from file given their path and complete t."""
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
# Writes only one node of the tree or its children (depending on maxDepth)
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
# Fix paths for writePaths
# 1. Must not contain CGNSTree first
# 2. Must not contain nodeName last
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
# Determine if a zone is a skeleton zone
# if ntype=0, based on coords + FlowSolution
# if ntype=1, based on coords only
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
# Converts a tree to a skeleton tree
# i.e. removes DataArray_t nodes and replaces them with None
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
# Converts a loaded skeleton tree to a partial tree
# rank=-1: removes skeleton zones (coord=None)
# rank>=0: removes zones with proc != rank
# Also removes skeleton fields
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
            # Completely removes the zone
            (p, c) = Internal.getParentOfNode(t, z)
            if Internal.isStdNode(t) == 0 and id(p) == id(t): del p[c]
            else: del p[2][c]
        else: # removes skeleton fields
            cx = Internal.getNodesFromType1(z, 'FlowSolution_t')
            for x in cx:
                t1 = Internal.getNodesFromType1(x, 'DataArray_t')
                for d in t1:
                    if d[1] is None:
                        p = Internal.getPath(x, d)
                        Internal._rmNodeByPath(z, p)
    return None

#==============================================================================
# Reads zones and replaces them in the tree
# IN: tree: input tree
# IN: rank: zones to load
# IN: or zoneNames: a list of zone names to load (['Base/Zone'])
# Warning: limited to adf and hdf
# Note: if a Proc node exists in the replaced zone, it is copied
# into the new zone
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
        # Paths of zones to be replaced
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
            if rank is not None: # copies solver data
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
# Writes zones to an already created file
# if zoneNames != None: writes specified zones
# if proc == rank: writes zones whose proc node matches
# if proc == -1: writes all zones of t
# Warning: limited to adf and hdf
#==============================================================================
def writeZones(t, fileName, format=None, proc=None, zoneNames=None, links=None, isize=8, rsize=8):
    """Write some zones in an existing file (adf or hdf)."""
    if zoneNames is None and proc is None: return None
    tp, ntype = Internal.node2PyTree(t)
    bases = Internal.getBases(tp)
    if proc is not None: # write zones by proc
        # Paths of nodes to be replaced (Zone and IntegralData)
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
# zones is a list of zones, but sorted by procs
# Ex: [ [zonesDeProc0], [zonesDeProc1], ...]
# IN: list of zones to replace in t
# IN: size: number of processors
# Puts zones in the tree by identifying zone names
# If zone name is base/zone, tries to find the base
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
# Updates the graph if zoneName on proc is to be sent to popp
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
# Merges two graphs
#==============================================================================
def mergeGraph(graph1, graph2):
    """Merge two graphs."""
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
# Returns the proc of z if z is local or in procDict if it exists
#==============================================================================
def getProcLocal__(z, procDict=None):
    if procDict is not None: return procDict[z[0]]
    else:
        proc = Internal.getNodeFromName2(z, 'proc')
        if proc is not None: proc = Internal.getValue(proc)
        else: proc = -1
        return proc

#==============================================================================
# Returns the proc of z if z is global or in procDict if it exists
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
# Computes the graph
# IN: t: tree containing fully loaded 'proc' nodes
# IN: type: type of graph wanted
# type='bbox' if bbox intersection of zones (full/bbox)
# type='bbox2' if bbox intersection and not on same base (full/bbox)
# type='bbox3' if bbox intersection with t2 (t:full/bbox and t2:full/bbox)
# type='match' if match between zones (full/skel/load skel/partial[+procDict])
# type='nmatch' if nearmatch between zones (full/skel/load skel/partial[+procDict])
# type='ID' if interpolation data between zones (full/skel/load skel/partial+procDict)
# type='IBCD' if IBC data between zones (full/skel/load skel/partial+procDict)
# type='ALLD' if all data (Interp+IBC) (full/skel/load skel/partial+procDict)
# type='proc' if zone has a proc node different from current proc (full/skel/load skel/partial)
# type='POST': t donor tree and t2 receptor tree. To interpolate data from donor to receptor where interpolation
#              data is stored in t in ID* subregions nodes. Requires procDict and procDict2
# intersectionsDict: Python dictionary containing the list of intersecting zones,
# as produced by "X.getIntersectingDomains". Note: if type='bbox3',
# the user must provide tree t2 and intersectionDict must describe the
# intersections between t and t2, as produced by "X.getIntersectingDomains(t,t2)".
# exploc: True if explicit local
# OUT: graph: dictionary containing send information
# of zones between processors
# graph is constructed such that:
# graph[P0][P1] returns the list of zones of processor P0
# "connected" with at least one zone of processor P1
# Ex: addXZones: sends zones graph[rank][opp] to proc opp,
# then waits for zones graph[opp][rank] for all opp.
#==============================================================================
def computeGraph(t, type='bbox', t2=None, procDict=None, rank=0,
                 intersectionsDict=None, exploc=False, procDict2=None,
                 it=None, reduction=True):
    """Return the communication graph for different block relation types."""
    zones = Internal.getZones(t)
    graph = {}
    if type == 'bbox': # zone of t from P0 intersect a zone of t from P1
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

    elif type == 'bbox2': # zone of t from P0 that intersects a zone of t from P1 but not on same base
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

    elif type == 'bbox3': # zone of t on P0 that intersects a zone of t2 on P1
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

    elif type == 'ID_Unsteady': # based on interpolation data
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
                        if idn != []: # the subRegion describes interpolations
                            popp = getProcGlobal__(donor, t, procDict)
                            if numero_iter not in graph_unsteady: graph_unsteady[numero_iter] = {}
                            #if(numero_iter==55):  print("subregion", s[0], z[0], donor, proc, popp)
                            updateGraph__(graph_unsteady[numero_iter], proc, popp, z[0])
                    else:
                        donor = Internal.getValue(s)
                        idn = Internal.getNodesFromName1(s, 'InterpolantsDonor')
                        if idn != []: # the subRegion describes interpolations
                            popp = getProcGlobal__(donor, t, procDict)
                            updateGraph__(graph_steady, proc, popp, z[0])

        return graph_steady, graph_unsteady

    elif type == 'ID': # based on interpolation data
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
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph, proc, popp, z[0])
                        else:
                            donor = Internal.getValue(s)
                            idn = Internal.getNodesFromName1(s, 'InterpolantsDonor')
                            if idn != []: # the subRegion describes interpolations
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
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr < levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==1 or ssiter%cycl==cycl//4 or ssiter%cycl==cycl//2-1 or ssiter%cycl==cycl//2+1 or ssiter%cycl==cycl//2+cycl//4 or ssiter%cycl==cycl-1):
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==cycl//2-1 or (ssiter%cycl==cycl//2 and (ssiter//cycl)%2==0) or ssiter%cycl==cycl-1):
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter > nssiter:
                            ssiter_ = ssiter - nssiter
                            if ssiter_%cycl==cycl//2 and (ssiter_//cycl)%2==1:
                                #if ssiter_%2==0 and ssiter_%cycl==cycl/2 and (ssiter_//cycl)%2==1:
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                list_graph_.append(graph_)

    elif type == 'IBCD': # based on IBC data
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
                    if idn != []: # the subRegion describes IBC
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
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr < levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==1 or ssiter%cycl==cycl//4 or ssiter%cycl==cycl//2-1 or ssiter%cycl==cycl//2+1 or ssiter%cycl==cycl//2+cycl//4 or ssiter%cycl==cycl-1):
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==cycl//2-1 or (ssiter%cycl==cycl//2 and (ssiter//cycl)%2==0) or ssiter%cycl==cycl-1):
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter > nssiter:
                            #if (ssiter%8==6):
                            ssiter_ = ssiter - nssiter
                            if ssiter_%2==0 and ssiter_%cycl==cycl//2 and (ssiter_//cycl)%2==1:
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                list_graph_.append(graph_)


    elif type == 'INV_IBCD': # based on inverse IBC data (from receiver)
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
                    if idn != []: # the subRegion describes IBC
                        popp = getProcGlobal__(donor, t, procDict)
                        updateGraph__(graph, popp, proc, donor)


    elif type == '2_IBCD': # based on IBC data
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
                    if idn != []: # the subRegion describes IBC
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
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr < levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==1 or ssiter%cycl==cycl//4 or ssiter%cycl==cycl//2-1 or ssiter%cycl==cycl//2+1 or ssiter%cycl==cycl//2+cycl//4 or ssiter%cycl==cycl-1):
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter <= nssiter:
                            if (ssiter%cycl==cycl//2-1 or (ssiter%cycl==cycl//2 and (ssiter//cycl)%2==0) or ssiter%cycl==cycl-1):
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                        if levdnr == levrcv and ssiter > nssiter:
                            #if (ssiter%8==6):
                            ssiter_ = ssiter - nssiter
                            if ssiter_%2==0 and ssiter_%cycl==cycl//2 and (ssiter_//cycl)%2==1:
                                if idn != []: # the subRegion describes interpolations
                                    popp = getProcGlobal__(donor, t, procDict)
                                    updateGraph__(graph_, proc, popp, z[0])
                list_graph_.append(graph_)


    elif type == 'ALLD': # based on Interpolations+IBC data
        for z in zones:
            proc = getProcLocal__(z, procDict)
            subRegions = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
            for s in subRegions:
                donor = Internal.getValue(s)
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # the subRegion describes interpolations/IBC
                    popp = getProcGlobal__(donor, t, procDict)
                    updateGraph__(graph, proc, popp, z[0])

    elif type == 'match': # based on match connections
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

    elif type == 'nmatch': # based on nearmatches
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
    """Split graph in one zone graph."""
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
# Returns the dictionary proc['blocName']
# from a distributed tree containing proc nodes
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
    """Set the proc to a zone or a set of zones."""
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
# Returns the dictionary prop['blocName']
# from a distributed tree containing 'propname' nodes
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
