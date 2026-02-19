# Part loader
from . import PyTree as C
from . import Internal
from . import Filter as Filter
from . import Mpi as Cmpi
import XCore.xcore
import numpy
from functools import partial

#============================================================
# add size to tree
#============================================================
def _addSizes2Zone(zone, zonePath, sizeData):
    """Create XX#Size node from sizeData."""
    for elmt in Internal.getNodesFromType1(zone, 'Elements_t'):
        elmtPath = zonePath+"/"+elmt[0]
        ecPath = elmtPath+"/ElementConnectivity"
        Internal.newIndexArray('ElementConnectivity#Size', value=sizeData[ecPath][2], parent=elmt)

    # zoneBC + BCDataSet
    zoneBCs = Internal.getNodesFromType1(zone, 'ZoneBC_t')
    for zoneBC in zoneBCs:
        zoneBCPath = zonePath+"/"+zoneBC[0]
        bcs = Internal.getNodesFromType1(zoneBC, 'BC_t')
        for bc in bcs:
            bcPath = zoneBCPath+"/"+bc[0]
            if Internal.getNodeFromName1(bc, 'PointList') is not None:
                plPath = bcPath+"/PointList"
                Internal.newIndexArray('PointList#Size', value=sizeData[plPath][2], parent=bc)
            for bcds in Internal.getNodesFromType1(bc, 'BCDataSet_t'):
                if Internal.getNodeFromName1(bcds, 'PointList') is not None:
                    plPath = bcPath+"/"+bcds[0]+"/PointList"
                    Internal.newIndexArray('PointList#Size', value=sizeData[plPath][2], parent=bcds)

    # ZoneGC
    zoneGCs = Internal.getNodesFromType1(zone, 'ZoneGridConnectivity_t')
    for zoneGC in zoneGCs:
        zoneGCPath = zonePath+"/"+zoneGC[0]
        gcs = Internal.getNodesFromType1(zoneGC, 'GridConnectivity_t') + \
            Internal.getNodesFromType1(zoneGC, 'GridConnectivity1to1_t')
        for gc in gcs:
            gcPath = zoneGCPath+"/"+gc[0]
            if Internal.getNodeFromName1(gc, 'PointList') is not None:
                plPath = gcPath+"/PointList"
                Internal.newIndexArray('PointList#Size', value=sizeData[plPath][2], parent=gc)
            if Internal.getNodeFromName1(gc, 'PointListDonor') is not None:
                pldPath = gcPath+"/PointListDonor"
                Internal.newIndexArray('PointListDonor#Size', value=sizeData[pldPath][2], parent=gc)

    # ZoneSubRegion
    zoneSubRegions = Internal.getNodesFromType1(zone, 'ZoneSubRegion_t')
    for zoneSubregion in zoneSubRegions:
        zoneSubregionPath = zonePath+"/"+zoneSubregion[0]
        if Internal.getNodeFromName1(zoneSubregion, 'PointList') is not None:
            plPath = zoneSubregionPath+"/PointList"
            Internal.newIndexArray('PointList#Size', value=sizeData[plPath][2], parent=zoneSubregion)

    # FlowSolutions
    FS = Internal.getNodesFromType1(zone, 'FlowSolution_t')
    for flowSol in FS:
        solPath = zonePath + "/" + Internal.getName(flowSol)
        if Internal.getNodeFromName1(flowSol, 'PointList') is not None:
            plPath = solPath+"/PointList"
            Internal.newIndexArray('PointList#Size', value=sizeData[plPath][2], parent=flowSol)

def _addSizes2Tree(sizeTree, sizeData):
    """Add size to tree."""
    bases = Internal.getNodesFromType1(sizeTree, 'CGNSBase_t')
    for b in bases:
        basePath = '/'+b[0]
        for z in Internal.getZones(b):
            zonePath = basePath+"/"+z[0]
            _addSizes2Zone(z, zonePath, sizeData)
    return None

def fixPointRanges(sizeTree):
    """Swap start and end in PointRange and PointRangeDonor"""
    bases = Internal.getBases(sizeTree)
    for base in bases:
        baseName = Internal.getName(base)
        zones = Internal.getZones(base)
        for zone in zones:
            zoneName = Internal.getName(zone)
            for zoneGC in Internal.getNodesFromType1(zone, 'ZoneGridConnectivity_t'):
                for gc in Internal.getNodesFromType1(zoneGC, 'GridConnectivity1to1_t'):
                    gcPath = baseName + '/' + zoneName
                    gcOppPath = Internal.getValue(gc)
                    if not '/' in gcOppPath:
                        gcOppPath = baseName + '/' + gcOppPath
                    if Internal.getNodeFromType1(gc, 'IndexRange_t') is not None:
                        transform = Internal.getValue(Internal.getNodeFromName1(gc, 'Transform'))
                        pointRange = Internal.getValue(Internal.getNodeFromName1(gc, 'PointRange'))
                        pointRange_d = Internal.getValue(Internal.getNodeFromName1(gc, 'PointRangeDonor'))

                        donor_dir = abs(transform) - 1
                        nb_points = pointRange[:,1] - pointRange[:,0]
                        nb_points_d = numpy.sign(transform)*(pointRange_d[donor_dir,1] - pointRange_d[donor_dir,0])
                        dir_to_swap = (nb_points != nb_points_d)

                        if gcPath < gcOppPath:
                            dir_to_swap = dir_to_swap[donor_dir]
                            pointRange_d[dir_to_swap, 0], pointRange_d[dir_to_swap, 1] = \
                                pointRange_d[dir_to_swap, 1], pointRange_d[dir_to_swap, 0]
                        elif gcPath > gcOppPath:
                            pointRange[dir_to_swap, 0], pointRange[dir_to_swap, 1] = \
                                pointRange[dir_to_swap, 1], pointRange[dir_to_swap, 0]
                        # If same base/zone, transform should be 1, 2, 3
                        else:
                            assert (dir_to_swap == False).all()

#=====================================================================
# load un squelette sur le rank 0 et bcast
# ajoute les noeuds size
#=====================================================================
def loadCollectiveSizeTree(filename):
    """Load skeleton, add size nodes and bcast."""
    if Cmpi.rank == 0:
        sizeData = {}
        t = C.convertFile2PyTree(filename,
                                 skeletonData=[3, 7],
                                 dataShape=sizeData,
                                 format='bin_hdf')
        _addSizes2Tree(t, sizeData)
        fixPointRanges(t)
    else: t = None
    t = Cmpi.bcast(t, root=0)
    return t

#============================================================
# distribution
#============================================================
def newDistribution(distributions=dict(), parent=None):
    """Create and return a CGNSNode to be used to store distribution data."""
    distriNode = Internal.newUserDefinedData(':CGNS#Distribution', None, parent)
    for name, value in distributions.items():
        Internal.newDataArray(name, value, parent=distriNode)
    return distriNode

def _createDistributionNodeFromDistrib(name, parent_node, distrib):
    distribUd = newDistribution(parent=parent_node)
    Internal.newDataArray(name, value=distrib, parent=distribUd)
    return None

def getDistribution(node, distriName=None):
    if distriName is not None:
        return Internal.getNodeFromPath(node, '/'.join([':CGNS#Distribution', distriName]))
    else:
        return Internal.getNodeFromName1(node, ':CGNS#Distribution')

def uniformDistributionAt(nelt, i, ninterval):
    step = nelt // ninterval
    remainder = nelt %  ninterval
    if i < remainder:
        inf = i * (step + 1)
        sup = inf + step + 1
    else:
        inf = i * step + remainder
        sup = inf + step
    return inf, sup

def uniformDistribution(nelt, comm):
    int_type = type(nelt)
    irank = int_type(Cmpi.rank)
    nrank = int_type(Cmpi.size)
    udist = uniformDistributionAt(nelt, irank, nrank)
    procIndices = numpy.empty(3, dtype=type(nelt))
    procIndices[0] = udist[0]
    procIndices[1] = udist[1]
    procIndices[2] = nelt
    return procIndices

def createDistributionNodeFromDistrib(name, parentNode, distrib):
    distribUd = newDistribution(parent=parentNode)
    Internal.newDataArray(name, value=distrib, parent=distribUd)

def computeElementsDistribution(zone, comm):
    if Internal.getZoneType(zone) == 1: pass
    else:
        elts = Internal.getNodesFromType1(zone, 'Elements_t')
        for elt in elts:
            er = Internal.getNodeFromName1(elt, 'ElementRange')
            n_totElmt = er[1][1] - er[1][0] + 1
            createDistributionNode(n_totElmt, comm, 'Element', elt)

def createDistributionNode(n_elt, comm, name, parent_node):
    distrib = uniformDistribution(n_elt, comm)
    createDistributionNodeFromDistrib(name, parent_node, distrib)

def computePLPRDistribution(node, comm):
    """Compute the distribution from PointList and PointRange."""
    PR = Internal.getNodeFromName1(node, 'PointRange')
    PR2 = Internal.getNodeFromName1(node, 'ElementRange')
    PL = Internal.getNodeFromName1(node, 'PointList')

    if PR:
        PR = PR[1].ravel('k')
        PRLength = PR[1] - PR[0] + 1
        createDistributionNode(PRLength, comm, 'Index', node)

    if PR2:
        PR2 = PR2[1].ravel('k')
        PRLength = PR2[1] - PR2[0] + 1
        createDistributionNode(PRLength, comm, 'Index', node)

    if PL:
        PLS = Internal.getNodeFromName1(node, 'PointList#Size')
        PLSize = Internal.getValue(PLS).prod()
        createDistributionNode(PLSize, comm, 'Index', node)

def _computeZoneDistribution(zone, comm):
    nVertex = C.getNPts(zone)
    nCell = C.getNCells(zone)
    distribVertex  = createDistributionNode(nVertex , comm, 'Vertex', zone)
    distribCell = createDistributionNode(nCell, comm, 'Cell'  , zone)
    computeElementsDistribution(zone, comm)
    #for zoneSubregion in Internal.getNodesFromType1(zone, 'ZoneSubRegion_t'):
    #  computePList_or_prange_distribution(zoneSubregion, comm)

    for flowSol in Internal.getNodesFromType1(zone, 'FlowSolution_t'):
        computePLPRDistribution(flowSol, comm)

    for zoneBC in Internal.getNodesFromType1(zone, 'ZoneBC_t'):
        for bc in Internal.getNodesFromType1(zoneBC, 'BC_t'):
            computePLPRDistribution(bc, comm)
            for bcds in Internal.getNodesFromType1(bc, 'BCDataSet_t'):
                computePLPRDistribution(bcds, comm)

    for zoneGC in Internal.getNodesFromType1(zone, 'ZoneGridConnectivity_t'):
        gcs = Internal.getNodesFromType1(zoneGC, 'GridConnectivity_t') + Internal.getNodesFromType1(zoneGC, 'GridConnectivity1to1_t')
        for gc in gcs:
            computePLPRDistribution(gc, comm)

def _addDistributionInfo(t, comm=Cmpi.KCOMM):
    for z in Internal.getZones(t):
        _computeZoneDistribution(z, comm)
    return None

#=============================================================
# Filters
#=============================================================
def PLPRSize(node):
    """Return size from PointList and Pointrange."""
    # CBX
    index = Internal.getNodesFromType(node, 'IndexArray_t')
    for ind in index:
        ps = Internal.getNodeFromName1(ind, 'PointList#Size')
        if ps is not None: return ps
        # No PL#Size, try to get info from :CGNS#Distribution and suppose size = 1,N
        distri = Internal.getVal(getDistribution(node, 'Index'))
        return numpy.array([1, distri[2]])
    index = Internal.getNodesFromType(node, 'IndexRange_t')
    for ind in index:
        # pas necessairement a plat?
        return ind[1].ravel('k')

def applyDataspaceToArrays(node, nodePath, dataSpace, hdfFilter):
    """Fill hdfFilter with dataSpace for DataArray."""
    for data in Internal.getNodesFromType1(node, 'DataArray_t'):
        path = nodePath+"/"+data[0]
        hdfFilter[path] = dataSpace

def applyDataspaceToPointlist(node, nodePath, dataSpace, hdfFilter):
    """Fill hdfFilter with dataSpace for PointList."""
    if Internal.getNodeFromName1(node, 'PointList') is not None:
        hdfFilter[nodePath + "/PointList"] = dataSpace
    if Internal.getNodeFromName1(node, 'PointListDonor') is not None:
        hdfFilter[nodePath + "/PointListDonor"] = dataSpace

def cell2Indexes(iCell, planSize, lineSize):
    """ Compute the (i,j,k) indices of a cell or a node from its global index."""
    k = iCell // planSize
    j = (iCell - k*planSize) // lineSize
    i = iCell - k*planSize - j*lineSize
    return i,j,k

def computeSlabs(arrayShape, gnumInterval):
    """ Compute HDF HyperSlabs to be used in order to contiguously load a part
    of a structured tridimensionnal array.

    arrayShape: Number of elements in x,y,z directions
    gnumInterval: semi open interval of elements to load in global numbering
    returns: list of hyperslabs

    Each slab is a list [[istart, iend], [jstart, jend], [kstart, kend]] of
    semi open intervals, starting at zero. The flattening order for the 3d
    array is increasing i, j, k.
    """
    hslabList = []
    nx, ny, nz = arrayShape
    lineSize = nx
    planSize = nx*ny

    ncellToLoad = gnumInterval[1] - gnumInterval[0]
    # print("{0} : cellInterval is [{1}:{2}[\n".format(iRank, gnumInterval[0], gnumInterval[1]))
    imin, jmin, kmin = cell2Indexes(gnumInterval[0],   planSize, lineSize)
    imax, jmax, kmax = cell2Indexes(gnumInterval[1]-1, planSize, lineSize)

    # print('toLoad : {0}  -- {1} {2} {3}  -- {4} {5} {6} \n'.format(
    #  ncellToLoad, imin, jmin, kmin, imax, jmax, kmax))

    istart = imin
    jstart = jmin
    kstart = kmin

    #If the line is full, merged it with next plan
    thisLineSize = min(nx, istart+ncellToLoad) - istart
    if thisLineSize != nx:
        jstart += 1
        if thisLineSize > 0:
            startLine  = [[istart, min(nx, istart+ncellToLoad)], [jmin, jmin+1], [kmin, kmin+1]]
            ncellToLoad -= thisLineSize
            hslabList.append(startLine)
            # print('startLine {0}, loaded {1} elmts\n'.format(
            # startLine, startLine[0][1] - startLine[0][0]))

    #If the plan is full, merged it with the block
    thisPlanSize = min(ny, jstart+(ncellToLoad // nx)) - jstart
    if thisPlanSize != ny:
        kstart += 1
        if thisPlanSize > 0:
            startPlane = [[0, nx], [jstart, min(ny, jstart+(ncellToLoad // nx))], [kmin, kmin+1]]
            ncellToLoad -= nx*thisPlanSize
            hslabList.append(startPlane)
            # print('startPlane {0}, loaded {1} lines ({2} elmts)\n'.format(
            # startPlane, startPlane[1][1] - startPlane[1][0], nx*(startPlane[1][1] - startPlane[1][0])))

    thisBlockSize = min(nz, kstart+(ncellToLoad // planSize)) - kstart
    if thisBlockSize > 0:
        centralBlock = [[0, nx], [0, ny], [kstart, min(nz, kstart+(ncellToLoad // planSize))]]
        ncellToLoad -= planSize*thisBlockSize
        hslabList.append(centralBlock)
        # print('centralBlock {0}, loaded {1} planes ({2} elmts)\n'.format(
        # centralBlock, centralBlock[2][1] - centralBlock[2][0], planSize*(centralBlock[2][1] - centralBlock[2][0])))

    if ncellToLoad >= nx:
        endPlane = [[0, nx], [0, (ncellToLoad // nx)], [kmax, kmax+1]]
        ncellToLoad -= nx*(endPlane[1][1] - endPlane[1][0])
        hslabList.append(endPlane)
        # print('endPlane {0}, loaded {1} lines ({2} elmts)\n'.format(
        # endPlane, endPlane[1][1] - endPlane[1][0], nx*(endPlane[1][1] - endPlane[1][0])))
    if ncellToLoad > 0:
        endLine = [[0, ncellToLoad], [jmax, jmax+1], [kmax, kmax+1]]
        ncellToLoad -= (endLine[0][1] - endLine[0][0])
        hslabList.append(endLine)
        # print('endLine {0}, loaded {1} elmts\n'.format(
        # endLine, endLine[0][1] - endLine[0][0]))
    assert(ncellToLoad == 0)

    return hslabList

def createCombinedDataspace(dataShape, distrib):
    """
    Create a dataspace from a flat distribution, but for arrays having a 3d (resp. 2d) stucture
    ie (Nx, Ny, Nz) (resp. (Nx, Ny)) numpy arrays.
    First, the 1d distribution is converted into slabs to load with the function compute_slabs.
    Those slabs are then combined to create the dataspace :
     for DSFile, we are expecting a list including all the slabs looking like
     [[startI_1, startJ_1, startK_1], [1,1,1], [nbI_1, nbJ_1, nbK_1], [1,1,1],
      [startI_2, startJ_2, startK_2], [1,1,1], [nbI_2, nbJ_2, nbK_2], [1,1,1], ...
      [startI_N, startJ_N, startK_N], [1,1,1], [nbI_N, nbJ_N, nbK_N], [1,1,1]]
     DSGlob me be the list of the tree dimensions sizes
     DSMmry and DSFrom have the same structure than flat / 1d dataspaces

    Mostly usefull for structured blocks.
    """
    slabList = computeSlabs(dataShape, distrib[0:2])
    dn_da = distrib[1] - distrib[0]
    DSFILEDA = []
    for slab in slabList:
        iS,iE, jS,jE, kS,kE = [item for bounds in slab for item in bounds]
        DSFILEDA.extend([[iS,jS,kS], [1,1,1], [iE-iS, jE-jS, kE-kS], [1,1,1]])
    DSMMRYDA = [[0]    , [1]    , [dn_da], [1]]
    DSFILEDA = list([list(DSFILEDA)])
    DSGLOBDA = [list(dataShape)]
    DSFORMDA = [[0]]
    return DSMMRYDA + DSFILEDA + DSGLOBDA + DSFORMDA

def createFlatDataspace(distrib):
    """Create the most basic dataspace (1d / flat) for a given distribution."""
    dn_da    = distrib[1] - distrib[0]
    DSMMRYDA = [[0         ], [1], [dn_da], [1]]
    DSFILEDA = [[distrib[0]], [1], [dn_da], [1]]
    DSGLOBDA = [[distrib[2]]]
    DSFORMDA = [[0]]
    return DSMMRYDA + DSFILEDA + DSGLOBDA + DSFORMDA

def createPeDataspace(distrib):
    """Create a dataspace from a flat distribution, of elements,
    but adapted to "ParentElements" arrays ie (N,2) numpy arrays.
    """
    dn_pe    = distrib[1] - distrib[0]
    DSMMRYPE = [[0              , 0], [1, 1], [dn_pe, 2], [1, 1]]
    DSFILEPE = [[distrib[0], 0], [1, 1], [dn_pe, 2], [1, 1]]
    DSGLOBPE = [[distrib[2], 2]]
    DSFORMPE = [[1]]
    return DSMMRYPE + DSFILEPE + DSGLOBPE + DSFORMPE

def createPointlistDataspace(distrib):
    """
    Create a dataspace from a flat distribution, but adapted to "fake 2d" arrays
    ie (1,N) numpy arrays.
    Mostly usefull for PointList arrays and DataArray of the related BCDataSets.
    """
    dnPL    = distrib[1] - distrib[0]
    DSMMRYPL = [[0,0          ], [1, 1], [1, dnPL], [1, 1]]
    DSFILEPL = [[0, distrib[0]], [1, 1], [1, dnPL], [1, 1]]
    DSGLOBPL = [[1, distrib[2]]]
    DSFORMPL = [[0]]
    return DSMMRYPL + DSFILEPL + DSGLOBPL + DSFORMPL

def createDataArrayFilter(distrib, dataShape=None):
    """
    Create an hdf dataspace for the given distribution. The kind of
    dataspace depends of the dataShape optional argument, representing
    the size of the array for which the dataspace is created in each dimension:
    - If dataShape is None or a single value, dataspace is 1d/flat
    - If dataShape is a 2d list [1, N], a dataspace adapted to pointlist is created
    - In other cases (which should correspond to true 2d array or 3d array), the
      dataspace is created from combine method (flat in memory, block in file).
    """
    if dataShape is None or len(dataShape) == 1: #Unstructured
        hdfDataSpace = createFlatDataspace(distrib)
    elif len(dataShape) == 2 and dataShape[0] == 1:
        hdfDataSpace = createPointlistDataspace(distrib)
    elif len(dataShape) == 2:
        hdfDataSpace = createFlatDataspace(distrib) # for pointrange (CB)
    else: # Structured
        hdfDataSpace = createCombinedDataspace(dataShape, distrib)

    return hdfDataSpace

def genElemts(zone):
    elmtsIni = Internal.getNodesFromType1(zone, 'Elements_t')
    for elmt in elmtsIni:
        yield elmt

def getSubregionExtent(subRegionNode, zone):
    """
    Return the path of the node (starting from zone node) related to subRegion_node
    node (BC, GC or itself)
    """
    #if Internal.getNodeFromName1(subRegion_node, "BCRegionName") is not None:
    #  for zbc, bc in iterNodesWithParentsByMatching(zone, "ZoneBC_t/BC_t"):
    #    if Internal.getName(bc) == Internal.getValue(Internal.getNodeFromName1(subRegion_node, "BCRegionName")):
    #      return Internal.getName(zbc) + '/' + Internal.getName(bc)
    #elif Internal.getNodeFromName1(subRegion_node, "GridConnectivityRegionName") is not None:
    #  gcPathes = ["ZoneGridConnectivity_t/GridConnectivity_t", "ZoneGridConnectivity_t/GridConnectivity1to1_t"]
    #  for gcPath in gcPathes:
    #    for zgc, gc in iterNodesWithParentsByMatching(zone, gcPath):
    #      if Internal.getName(gc) == Internal.getValue(Internal.getNodeFromName1(subRegion_node, "GridConnectivityRegionName")):
    #        return Internal.getName(zgc) + '/' + Internal.getName(gc)
    #else:
    return Internal.getName(subRegionNode)

def createZoneEsoElementsFilter(elmt, zonePath, hdfFilter, mode):
    distribElmt = Internal.getVal(getDistribution(elmt, 'Element'))
    dnElmt = distribElmt[1] - distribElmt[0]

    # For NGon only
    pe = Internal.getNodeFromName1(elmt, 'ParentElements')
    if pe:
        dataSpace = createPeDataspace(distribElmt)
        #hdfFilter[f"{zonePath}/{Internal.getName(elmt)}/ParentElements"] = dataSpace
        hdfFilter["%s/%s/ParentElements"%(zonePath,Internal.getName(elmt))] = dataSpace
        if Internal.getNodeFromName1(elmt, 'ParentElementsPosition'):
            #hdfFilter[f"{zonePath}/{Internal.getName(elmt)}/ParentElementsPosition"] = dataSpace
            hdfFilter["%s/%s/ParentElementsPosition"%(zonePath,Internal.getName(elmt))] = dataSpace
    eso = Internal.getNodeFromName1(elmt, 'ElementStartOffset')
    esoPath = None
    if eso:
        # Distribution for NGon/NFace -> ElementStartOffset is the same than DistrbutionFace, except
        # that the last proc have one more element
        nElmt = distribElmt[2]
        if mode == 'read':
            dnElmtIndex = dnElmt + 1 # + int(distribElmt[1] == nElmt)
        elif mode == 'write':
            dnElmtIndex = dnElmt + int((distribElmt[1] == nElmt) and (distribElmt[0] != distribElmt[1]))
        DSMMRYESO = [[0              ], [1], [dnElmtIndex], [1]]
        DSFILEESO = [[distribElmt[0]], [1], [dnElmtIndex], [1]]
        DSGLOBESO = [[nElmt+1]]
        DSFORMESO = [[0]]

        esoPath = zonePath+"/"+elmt[0]+"/ElementStartOffset"
        hdfFilter[esoPath] = DSMMRYESO + DSFILEESO + DSGLOBESO + DSFORMESO

    ec = Internal.getNodeFromName1(elmt, 'ElementConnectivity')
    if ec:
        if esoPath is None:
            raise RuntimeError("In order to load ElementConnectivity, the ElementStartOffset is mandatory")
        ecPath = zonePath+"/"+elmt[0]+"/ElementConnectivity"
        hdfFilter[ecPath] = partial(loadElementConnectivityFromEso, elmt, zonePath)

def loadElementConnectivityFromEso(elmt, zonePath, hdfFilter):
    distribUd = getDistribution(elmt)
    distribElmt = Internal.getNodeFromName1(distribUd, 'Element')[1]
    #dnElmt = distribElmt[1] - distribElmt[0]

    esoNode = Internal.getNodeFromName1(elmt, 'ElementStartOffset') # Maintenant il est chargÃ©
    if esoNode[1] is None: raise RuntimeError
    eso = esoNode[1]

    begFaceVertex = eso[0]
    endFaceVertex = eso[eso.shape[0]-1]
    dnFaceVertex = endFaceVertex - begFaceVertex

    # print("begFaceVertex::", begFaceVertex)
    # print("endFaceVertex::", endFaceVertex)
    distrib_n = None
    ecSize_n  = Internal.getNodeFromName1(elmt, 'ElementConnectivity#Size')
    if ecSize_n is not None:
        nFaceVertex = numpy.prod(ecSize_n[1])
    else:
        distribUd = getDistribution(elmt)
        distrib_n = Internal.getNodeFromName1(distribUd, "ElementConnectivity")
        nFaceVertex = distrib_n[1][2]

    # print("nFaceVertex::", nFaceVertex)

    nFace = distribElmt[2]
    #dnFaceIndex = dnElmt + int(distribElmt[1] == nFace)
    DSMMRYEC = [[0           ], [1], [dnFaceVertex], [1]]
    DSFILEEC = [[begFaceVertex], [1], [dnFaceVertex], [1]]
    DSGLOBEC = [[nFaceVertex ]]
    DSFORMEC = [[0]]

    ecPath = zonePath+"/"+elmt[0]+"/ElementConnectivity"
    hdfFilter[ecPath] = DSMMRYEC + DSFILEEC + DSGLOBEC + DSFORMEC

    if distrib_n is None:
        distrib = numpy.empty(3, dtype=eso.dtype)
        distrib[0] = begFaceVertex
        distrib[1] = endFaceVertex
        distrib[2] = nFaceVertex
        Internal.newDataArray("ElementConnectivity", value=distrib, parent=distribUd)

def createZoneStdElementsFilter(elmt, zonePath, hdfFilter):
    distribElmt = Internal.getVal(getDistribution(elmt, 'Element'))
    dnElmt = distribElmt[1] - distribElmt[0]

    _, elmt_npe = Internal.eltNo2EltName(elmt[1][0]) # nbre de vertex de l'element

    DSMMRYElmt = [[0                       ], [1], [dnElmt*elmt_npe], [1]]
    DSFILEElmt = [[distribElmt[0]*elmt_npe], [1], [dnElmt*elmt_npe], [1]]
    DSGLOBElmt = [[distribElmt[2]*elmt_npe]]
    DSFORMElmt = [[0]]

    path = zonePath+"/"+elmt[0]+"/ElementConnectivity"
    hdfFilter[path] = DSMMRYElmt + DSFILEElmt + DSGLOBElmt + DSFORMElmt

    pe = Internal.getNodeFromName1(elmt, 'ParentElements')
    if pe:
        dataSpace = createPeDataspace(distribElmt)
        #hdfFilter[f"{zonePath}/{Internal.getName(elmt)}/ParentElements"] = dataSpace
        hdfFilter["%s/%s/ParentElements"%(zonePath,Internal.getName(elmt))] = dataSpace

        if Internal.getNodeFromName1(elmt, 'ParentElementsPosition'):
            #hdfFilter[f"{zonePath}/{Internal.getName(elmt)}/ParentElementsPosition"] = dataSpace
            hdfFilter["%s/%s/ParentElementsPosition"%(zonePath,Internal.getName(elmt))] = dataSpace

def createZoneElementsFilter(zone, zonePath, hdfFilter, mode):
    """Prepare the hdfFilter for all the Element_t nodes found in the zone."""
    zoneElmts = genElemts(zone)
    for elmt in zoneElmts:
        if elmt[1][0] == 22 or elmt[1][0] == 23:
            createZoneEsoElementsFilter(elmt, zonePath, hdfFilter, mode)
        elif elmt[1][0] == 20:
            raise ValueError('MIXED elements not implemented.')
        else:
            createZoneStdElementsFilter(elmt, zonePath, hdfFilter)

def createZoneBCFilter(zone, zonePath, hdfFilter):
    """
    Fill up the hdf filter for the BC_t nodes present in
    the zone.
    Filter is created for the following nodes:
     - PointList (if present = unstruct. only)
     - All arrays founds in BCDataSets. Those arrays are supposed
       to be shaped as the PointList array. If a BCDataSet contains
       no PointList/PointRange node, the data is assumed to be consistent
       with the PointList/PointRange of the BC. Otherwise, the PointList/
       PointRange node of the BCDataSet is used to set the size of the BCData
       arrays. In this case, the PointList (if any) of the BCDataSet is
       written in the filter as well.
    """
    for zoneBC in Internal.getNodesFromType1(zone, 'ZoneBC_t'):
        zoneBCPath = zonePath+"/"+zoneBC[0]
        for bc in Internal.getNodesFromType1(zoneBC, 'BC_t'):
            bcPath = zoneBCPath+"/"+bc[0]

            distribBC = Internal.getVal(getDistribution(bc, 'Index'))
            bcShape = PLPRSize(bc)
            dataSpace = createDataArrayFilter(distribBC, bcShape)
            applyDataspaceToPointlist(bc, bcPath, dataSpace, hdfFilter)

            for bcds in Internal.getNodesFromType1(bc, "BCDataSet_t"):
                bcdsPath = bcPath + "/" + bcds[0]
                distrib_bcds_n = getDistribution(bcds)

                if distrib_bcds_n is None: #BCDS uses BC distribution
                    distribData = distribBC
                    dataShape = bcShape
                else: #BCDS has its own distribution
                    distribData = Internal.getNodeFromName1(distrib_bcds_n, 'Index')[1]
                    dataShape = PLPRSize(bcds)

                dataSpacePL = createDataArrayFilter(distribData, dataShape)
                #BCDataSet always use flat data array
                dataSpace_array = createDataArrayFilter(distribData, [dataShape.prod()])
                applyDataspaceToPointlist(bcds, bcdsPath, dataSpacePL, hdfFilter)
                for bcdata in Internal.getNodesFromType1(bcds, 'BCData_t'):
                    bcdataPath = bcdsPath + "/" + bcdata[0]
                    applyDataspaceToArrays(bcdata, bcdataPath, dataSpace_array, hdfFilter)


def createZoneGridConnectivityFilter(zone, zonePath, hdfFilter):
    """
    Fill up the hdf filter for the GC_t nodes present in the zone.
    For unstructured GC (GridConnectivity_t), the filter is set up for
    the PointList and PointListDonor arrays.
    Structured GC (GridConnectivity1to1_t) are skipped since there is
    no data to load for these nodes.
    """
    for zoneGC in Internal.getNodesFromType1(zone, 'ZoneGridConnectivity_t'):
        zoneGCPath = zonePath+"/"+zoneGC[0]
        for gc in Internal.getNodesFromType1(zoneGC, 'GridConnectivity_t'):
            gcPath = zoneGCPath+"/"+gc[0]
            distrib_ia = Internal.getVal(getDistribution(gc, 'Index'))
            gcShape = PLPRSize(gc)
            dataSpace = createDataArrayFilter(distrib_ia, gcShape)
            applyDataspaceToPointlist(gc, gcPath, dataSpace, hdfFilter)

def createFlowSolutionFilter(zone, zonePath, hdfFilter):
    """
    Fill up the hdf filter for the FlowSolution_t nodes present in the
    zone. The size of the dataspace are computed from the pointList node
    if present, or using allCells / allVertex if no pointList is present.
    Filter is created for the arrays and for the PointList if present
    """
    distribVertex = Internal.getVal(getDistribution(zone, 'Vertex'))
    distribCell = Internal.getVal(getDistribution(zone, 'Cell'))
    for flowSolution in Internal.getNodesFromType1(zone, 'FlowSolution_t'):
        flowSolutionPath = zonePath + "/" + Internal.getName(flowSolution)
        gridLocation = Internal.getNodeFromType1(flowSolution, 'GridLocation_t')
        if gridLocation is None: gridLocation = 'Vertex'
        else: gridLocation = Internal.getValue(gridLocation)
        distrib_ud_n = getDistribution(flowSolution)
        if distrib_ud_n:
            distribData = Internal.getNodeFromName1(distrib_ud_n, 'Index')[1]
            dataShape = PLPRSize(flowSolution)
            dataSpacePL = createDataArrayFilter(distribData, dataShape)
            dataSpace = createDataArrayFilter(distribData, [dataShape.prod()])
            applyDataspaceToPointlist(flowSolution, flowSolutionPath, dataSpacePL, hdfFilter)
        elif gridLocation == 'CellCenter':
            dataSpace = createDataArrayFilter(distribCell, zone[1][:,1])
        elif gridLocation == 'Vertex':
            dataSpace = createDataArrayFilter(distribVertex, zone[1][:,0])
        else:
            raise RuntimeError("GridLocation %s is not allowed without PL"%gridLocation)
        applyDataspaceToArrays(flowSolution, flowSolutionPath, dataSpace, hdfFilter)

def createZoneSubregionFilter(zone, zonePath, hdfFilter):
    """
    Fill up the hdf filter for the ZoneSubRegion_t nodes present in
    the zone.
    The size of the dataspace are computed from
    - the corresponding BC or GC if the subregion is related to one of them
      (this information is given by a BCRegionName or GridConnectivityRegionName
      node)
    - the PointList / PointSize node of the subregion otherwise.
    Filter is created for the following nodes:
     - All arrays present in ZoneSubRegion;
     - PointList array if the zone is unstructured and if the subregion
       is not related to a BC/GC.
    """
    zoneSubRegions = Internal.getNodesFromType1(zone, 'ZoneSubRegion_t')
    for zoneSubregion in zoneSubRegions:
        zoneSubregionPath = zonePath+"/"+zoneSubregion[0]

        # Search matching region
        matchingRegionPath = getSubregionExtent(zoneSubregion, zone)
        matchingRegion = Internal.getNodeFromPath(zone, matchingRegionPath)

        distrib_ud_n = getDistribution(matchingRegion)
        if not distrib_ud_n:
            raise RuntimeError("ZoneSubRegion {0} is not well defined".format(zoneSubregion[0]))
        distribData = Internal.getNodeFromName1(distrib_ud_n, 'Index')[1]

        dataShape = PLPRSize(matchingRegion)
        dataSpacePL = createDataArrayFilter(distribData, dataShape)
        dataSpaceAr = createDataArrayFilter(distribData, [dataShape.prod()])

        applyDataspaceToPointlist(zoneSubregion, zoneSubregionPath, dataSpacePL, hdfFilter)
        applyDataspaceToArrays(zoneSubregion, zoneSubregionPath, dataSpaceAr, hdfFilter)

def createZoneFilter(zone, zonePath, hdfFilter, mode):
    """
    Fill up the hdf filter for the following elements of the zone:
    Coordinates, Elements (NGon / NFace, Standards), FlowSolution
    (vertex & cells only), ZoneSubRegion, ZoneBC (including BCDataSet)
    and ZoneGridConnectivity.

    The bounds of the filter are determined by the :CGNS#Distribution
    node and, for the structured zones, by the size of the blocks.
    """
    # Coords
    distribVertex = Internal.getVal(getDistribution(zone, 'Vertex'))
    allVertexDataspace = createDataArrayFilter(distribVertex, zone[1][:,0])
    for GC in Internal.getNodesFromType1(zone, 'GridCoordinates_t'):
        GCPath = zonePath + "/" + Internal.getName(GC)
        applyDataspaceToArrays(GC, GCPath, allVertexDataspace, hdfFilter)
    createZoneElementsFilter(zone, zonePath, hdfFilter, mode)
    createZoneBCFilter(zone, zonePath, hdfFilter)
    createZoneGridConnectivityFilter(zone, zonePath, hdfFilter)
    createFlowSolutionFilter(zone, zonePath, hdfFilter)
    createZoneSubregionFilter(zone, zonePath, hdfFilter)

def createTreeHdfFilter(distTree, hdfFilter, mode='read'):
    bases = Internal.getBases(distTree)
    for b in bases:
        zones = Internal.getZones(b)
        for z in zones:
            zonePath = "/"+Internal.getName(b)+"/"+Internal.getName(z)
            createZoneFilter(z, zonePath, hdfFilter, mode)

#==========================================================
# load
#==========================================================
def updateTreeWithPartialLoadDict(distTree, partialDictLoad):
    for path, data in partialDictLoad.items():
        node = Internal.getNodeFromPath(distTree, path)
        node[1] = data

def loadTreeFromFilter(filename, distTree, comm, hdfFilter):
    hdfFilterWithDim = {key: value for (key, value) in hdfFilter.items() \
                        if isinstance(value, (list, tuple))}

    partialDictLoad = C.convertFile2PartialPyTreeFromPath(filename, hdfFilterWithDim, comm)
    updateTreeWithPartialLoadDict(distTree, partialDictLoad)

    # > Match with callable
    hdfFilterWithFunc = {key: value for (key, value) in hdfFilter.items() \
                         if not isinstance(value, (list, tuple))}
    unlockedOnce = True
    while (len(hdfFilterWithFunc) > 0 and unlockedOnce):
        nextHdfFilter = dict()
        unlockedOnce = False
        for key, f in hdfFilterWithFunc.items():
            try:
                f(nextHdfFilter)
                unlockedOnce = True
            except RuntimeError: # Not ready yet
                pass
        partialDictLoad = C.convertFile2PartialPyTreeFromPath(filename, nextHdfFilter, comm)
        updateTreeWithPartialLoadDict(distTree, partialDictLoad)
        hdfFilterWithFunc = {key: value for (key, value) in nextHdfFilter.items() \
                             if not isinstance(value, (list, tuple))}

    if unlockedOnce is False:
        raise RuntimeError("Something strange in the loading process")

def cleanDistributionInfo(distTree):
    bases = Internal.getBases(distTree)
    for b in bases:
        for zone in Internal.getNodesFromType1(b, 'Zone_t'):
            Internal._rmNodesByName1(zone, ':CGNS#Distribution')
            for elmt in Internal.getNodesFromType1(zone, 'Elements_t'):
                Internal._rmNodesByName1(elmt, ':CGNS#Distribution')
                Internal._rmNodesByName1(elmt, 'ElementConnectivity#Size')
            for zoneBC in Internal.getNodesFromType1(zone, 'ZoneBC_t'):
                for bc in Internal.getNodesFromType1(zoneBC, 'BC_t'):
                    Internal._rmNodesByName2(bc, ':CGNS#Distribution')
                    Internal._rmNodesByName2(bc, 'PointList#Size')
            for zoneGC in Internal.getNodesFromType1(zone, 'ZoneGridConnectivity_t'):
                for gc in Internal.getNodesFromType1(zoneGC, 'GridConnectivity_t') + \
                        Internal.getNodesFromType1(zoneGC, 'GridConnectivity1to1_t'):
                    Internal._rmNodesByName1(gc, ':CGNS#Distribution')
                    Internal._rmNodesByName1(gc, 'PointList#Size')
            for zoneSubregion in Internal.getNodesFromType1(zone, 'ZoneSubRegion_t'):
                Internal._rmNodesByName1(zoneSubregion, ':CGNS#Distribution')
                Internal._rmNodesByName1(zoneSubregion, 'PointList#Size')
            for zoneSol in Internal.getNodesFromType1(zone, 'FlowSolution_t'):
                Internal._rmNodesByName1(zoneSol, ':CGNS#Distribution')
                Internal._rmNodesByName1(zoneSol, 'PointList#Size')

def saveTreeFromFilter(fileName, distTree, comm, hdfFilter):
    hdfFilterWithDim = {key: value for (key, value) in hdfFilter.items() if isinstance(value, list)}
    hdfFilterWithFunc = {key: value for (key, value) in hdfFilter.items() if not isinstance(value, list)}
    #next_hdfFilter = dict()
    for key, f in hdfFilterWithFunc.items():
        f(hdfFilterWithDim)

    # Dont save distribution info, but work on a copy to keep it for further use
    copyDistTree = Internal.copyRef(distTree)
    cleanDistributionInfo(copyDistTree)
    C.convertPyTree2FilePartial(copyDistTree, fileName, comm, hdfFilterWithDim, ParallelHDF=True)

#========================================================
# resume
#========================================================
def loadAsChunks(fileName):
    distTree = loadCollectiveSizeTree(fileName)
    _addDistributionInfo(distTree)
    hdfFilter = {}
    createTreeHdfFilter(distTree, hdfFilter)
    loadTreeFromFilter(fileName, distTree, Cmpi.KCOMM, hdfFilter)
    Cmpi._setProc(distTree, Cmpi.rank)
    return distTree

#========================================================
def chunk2part(dt):
    BCType_l = set(Internal.KNOWNBCS)
    arrays = []
    zones = Internal.getZones(dt)

    z = zones[0]
    GC = Internal.getNodeFromName1(z, Internal.__GridCoordinates__)
    cx = Internal.getNodeFromName1(GC, 'CoordinateX')[1]
    cy = Internal.getNodeFromName1(GC, 'CoordinateY')[1]
    cz = Internal.getNodeFromName1(GC, 'CoordinateZ')[1]

    ngon = Internal.getNGonNode(z)
    ngonc = Internal.getNodeFromName1(ngon, 'ElementConnectivity')[1]
    ngonso = Internal.getNodeFromName1(ngon, 'ElementStartOffset')[1]

    nface = Internal.getNFaceNode(z)
    nfacec = Internal.getNodeFromName1(nface, 'ElementConnectivity')[1]
    nfaceso = Internal.getNodeFromName1(nface, 'ElementStartOffset')[1]

    fsolc = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    solc = []; solcNames = []
    if fsolc is not None:
        for f in fsolc[2]:
            if f[3] == 'DataArray_t':
                solc.append(f[1]); solcNames.append(f[0])

    fsol = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
    soln = []; solNames = []
    if fsol is not None:
        for f in fsol[2]:
            if f[3] == 'DataArray_t':
                soln.append(f[1]); solNames.append(f[0])

    zonebc = Internal.getNodeFromType1(z, 'ZoneBC_t')
    bcs = []
    bcNames = []; bcTypes = {}
    familyNames = {}
    if zonebc is not None:
        BCs = Internal.getNodesFromType1(zonebc, 'BC_t')
        for bc in BCs:
            bcname = bc[0]; bctype = Internal.getValue(bc)

            if bctype == 'FamilySpecified':
                fname = Internal.getNodeFromType1(bc, 'FamilyName_t')
                fn = Internal.getValue(fname)
                bcTypes[bcname] = fn
            else:
                bcTypes[bcname] = bctype
            bcNames.append(bcname)
            plist = Internal.getNodeFromName1(bc, 'PointList')
            bcs.append(plist[1][0])

    arrays.append([cx,cy,cz,ngonc,ngonso,nfacec,nfaceso,solc,soln,bcs])

    RES = XCore.xcore.chunk2partNGon(arrays)
    (mesh, commData, solc, sol, bcs, cells, faces, points) = RES
    Cmpi.barrier()

    # create zone
    zo = Internal.createZoneNode('Zone_' + '%d'%Cmpi.rank, mesh)

    # add solutions
    for n, name in enumerate(solNames):
        cont = Internal.createUniqueChild(zo, Internal.__FlowSolutionNodes__, 'FlowSolution_t')
        Internal.newDataArray(name, value=sol[n], parent=cont)

    for n, name in enumerate(solcNames):
        cont = Internal.createUniqueChild(zo, Internal.__FlowSolutionCenters__, 'FlowSolution_t')
        Internal._createUniqueChild(cont, 'GridLocation', 'GridLocation_t', value='CellCenter', )
        Internal.newDataArray(name, value=solc[n], parent=cont)

    for i in range(len(bcs)):
        if len(bcs[i]) != 0:
            cont = Internal.createUniqueChild(zo, 'ZoneBC', 'ZoneBC_t')
            val = bcTypes[bcNames[i]]
            if val not in BCType_l:
                Internal.newBC(name=bcNames[i], pointList=bcs[i], family=val, parent=cont)
            else:
                Internal.newBC(name=bcNames[i], pointList=bcs[i], btype=val, parent=cont)

    t = C.newPyTree(['Base', zo])
    Cmpi._setProc(t, Cmpi.rank)
    Internal._correctPyTree(t, level=7)

    return t, RES

#========================================================
def loadAndSplit(fileName):
    """Load and split a file containing one NGON."""
    dt = loadAsChunks(fileName)
    t, RES = chunk2part(dt)
    return t, RES

#========================================================
def mergeAndSave(distTree, fileName):
    """Save a distributed tree in one NGON zone."""
    hdfFilter = {}
    createTreeHdfFilter(distTree, hdfFilter)
    saveTreeFromFilter(fileName, distTree, Cmpi.KCOMM, hdfFilter)
    return None