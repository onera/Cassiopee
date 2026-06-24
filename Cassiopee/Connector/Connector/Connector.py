"""Find connectivity in grids.
"""
__version__ = '4.2'
__author__ = "Stephanie Peron, Christophe Benoit, Gaelle Jeanfaivre, Pascal Raud, Luis Bernardos"

from . import connector

__all__ = ['blankCells', '_blankCells', 'blankCellsTetra', 'blankCellsTri', 'blankIntersectingCells', 'chimeraTransfer', 'connectMatch',
           'getIntersectingDomainsAABB', 'maximizeBlankedCells', 'optimizeOverlap', 'setDoublyDefinedBC', 'setHoleInterpolatedPoints',
           'setIBCTransfers', 'setIBCTransfersD', 'setInterpTransfers', 'setInterpTransfersD', 'writeCoefs','maskXRay__',
           '_applyBCOverlapsStruct__', 'applyBCOverlapsStruct__', 'applyBCOverlapsNG__',
           'getInterpolatedPoints__', 'getEXPoints__', '_modCellN1', '_modCellN2', 'changeWall__']

#===============================================================================
def connectMatch(a1, a2, sameZone=0, tol=1.e-6, dim=3):
    """Find matching boundaries.
    Usage: connectMatch(a1, a2, sameZone, tol, dim)"""
    try: import Converter as C; import Transform as T
    except: raise ImportError("connectMatch requires Converter and Transform modules.")

    res = []
    if len(a1) != 5 or len(a2) != 5: print('Warning: connectMatch is valid only for structured grids.'); return res
    if dim == 2: nwins = 4
    elif dim == 3: nwins = 6
    else: raise ImportError("connectMatch: dim must be 2 or 3.")

    allWins=[]; dimsI=[]; dimsJ=[]; dimsK=[]; typeOfWins=[]; indirBlkOfWins=[]
    dimsI.append(a1[2]); dimsJ.append(a1[3]); dimsK.append(a1[4])
    dimsI.append(a2[2]); dimsJ.append(a2[3]); dimsK.append(a2[4])

    for win1 in range(1, nwins):
        imin1 = 1; jmin1 = 1; kmin1 = 1; imax1 = a1[2]; jmax1 = a1[3]; kmax1 = a1[4]
        if win1 == 1: imax1 = 1
        elif win1 == 2: imin1 = imax1
        elif win1 == 3: jmax1 = 1
        elif win1 == 4: jmin1 = jmax1
        elif win1 == 5: kmax1 = 1
        elif win1 == 6: kmin1 = kmax1
        win = T.subzone(a1,(imin1,jmin1,kmin1),(imax1,jmax1,kmax1))
        if win1 == 1 or win1 == 2: win = T.reorder(win,(3,1,2))
        elif win1 == 3 or win1 == 4: win = T.reorder(win,(1,3,2))
        allWins.append(win)
        indirBlkOfWins.append(0)
        typeOfWins.append(win1)
    if sameZone == 0:
        for win2 in range(1, nwins):
            imin2 = 1; jmin2 = 1; kmin2 = 1; imax2 = a2[2]; jmax2 = a2[3]; kmax2 = a2[4]
            if win2 == 1: imax2 = 1
            elif win2 == 2: imin2 = imax2
            elif win2 == 3: jmax2 = 1
            elif win2 == 4: jmin2 = jmax2
            elif win2 == 5: kmax2 = 1
            elif win2 == 6: kmin2 = kmax2
            win = T.subzone(a2,(imin2,jmin2,kmin2),(imax2,jmax2,kmax2))
            if win2 == 1 or win2 == 2: win = T.reorder(win,(3,1,2))
            elif win2 == 3 or win2 == 4: win = T.reorder(win,(1,3,2))
            allWins.append(win)
            indirBlkOfWins.append(1)
            typeOfWins.append(win2)

    allWins = C.extractVars(allWins, ['x','y','z'])
    allTags = C.node2Center(allWins)
    allTags = C.initVars(allTags, 'tag1', -1.) # defines the opposite window
    allTags = C.initVars(allTags, 'tag2', -2.) # defines the opposite index in opposite window

    allTags = identifyMatching(allTags, tol)
    allTags = C.extractVars(allTags, ['tag1','tag2'])

    # Gather matching cells into structured patches [ [[noz1,noz2],[imin1,imax1,...],[imin2,imax2,...],trirac] ]
    infos = gatherMatching(allWins, allTags, typeOfWins, indirBlkOfWins, dimsI, dimsJ, dimsK, dim, tol)
    for info in infos:
        range1 = info[1]
        range2 = info[2]
        topp0 = info[3]
        noz1 = info[0][0]; noz2 = info[0][1]
        if sameZone == 0 and noz1 != noz2 : res.append([range1,range2,topp0])
        elif sameZone == 1 and noz1 == noz2 and range1 != range2:  res.append([range1,range2,topp0])
    return res

#===============================================================================
def identifyMatching(A, tol=1.e-6):
    return connector.identifyMatching(A, tol)

#===============================================================================
def identifyMatchingP(A, AP, tol=1.e-6):
    return connector.identifyMatchingP(A, AP, tol)

#===============================================================================
def identifyMatchingNM(AR, AD, tol=1.e-6):
    return connector.identifyMatchingNM(AR, AD, tol)

#===============================================================================
def identifyDegenerated(A, tol=1.e-6):
    return connector.identifyDegenerated(A, tol)

#===============================================================================
# Returns [ [[noz1,noz2],[imin1,imax1,...],[imin2,imax2,...],trirac] ]
def gatherMatching(listOfWins, listOfTags, typeOfWins, blkOfWins,
                   allNI, allNJ, allNK, dim=3, tol=1.e-6):
    return connector.gatherMatching(listOfWins, listOfTags, typeOfWins,
                                    blkOfWins, allNI, allNJ, allNK, dim, tol)

def gatherMatchingNM(listOfWins, listOfTags, typeOfWins, blkOfWins,
                     allNI, allNJ, allNK, dim=3, tol=1.e-6):
    return connector.gatherMatchingNM(listOfWins, listOfTags, typeOfWins,
                                      blkOfWins, allNI, allNJ, allNK, dim, tol)

def gatherMatchingNGon__(tagsF, allExtIndices):
    return connector.gatherMatchingNGon(tagsF, allExtIndices)

def gatherDegenerated(listOfTags, typeOfWins, blkOfWins,
                      allNI, allNJ, allNK, dim=3):
    return connector.gatherDegenerated(listOfTags, typeOfWins, blkOfWins,
                                       allNI, allNJ, allNK, dim)
#------------------------------------------------------------------------------
# Check if the two matching windows have opposite normals
#------------------------------------------------------------------------------
def checkMatchWindows__(a, win, aopp, winopp, tol):
    """Check if matching windows are opposed in terms of normals."""
    return connector.checkMatchWindows(a, win, aopp, winopp, tol)

#------------------------------------------------------------------------------
def optimizeOverlap(nodes1, centers1, nodes2, centers2, prio1=0, prio2=0, isDW=0):
    """Optimize the overlap of grids defined by nodes1 and nodes2
    centers1 and centers2 define the coordinates of cell centers, modified by
    the double wall algorithm + cellN variable. 
    Usage: optimizeOverlap(nodes1, centers1, nodes2, centers2,prio1=0, prio2=0,isDW=0)"""
    import KCore.kcore as KCore
    import Converter as C
    posv1 = KCore.isNamePresent(centers1,'vol')
    posv2 = KCore.isNamePresent(centers2,'vol')
    if posv1 == -1:
        try: import Generator as G
        except: raise ImportError("optimizeOverlap requires Converter and Generator modules.")
        vol1 = G.getVolumeMap(nodes1); centers1 = C.addVars([centers1,vol1])
    if posv2 == -1:
        try: import Generator as G
        except: raise ImportError("optimizeOverlap requires Converter and Generator modules.")
        vol2 = G.getVolumeMap(nodes2); centers2 = C.addVars([centers2,vol2])

    extCenters1 = C.node2ExtCenter(nodes1)
    extCenters2 = C.node2ExtCenter(nodes2)
    hook1 = C.createHook([extCenters1],'extractMesh')
    hook2 = C.createHook([extCenters2],'extractMesh')
    res = optimizeOverlap__(extCenters1, centers1, extCenters2, centers2,
                            prio1, prio2, isDW, hook1, hook2)
    if posv1 == -1: res[0] = C.rmVars(res[0],'vol')
    if posv2 == -1: res[1] = C.rmVars(res[1],'vol')
    return res

def optimizeOverlap__(extCenters1, centers1, extCenters2, centers2,
                      prio1=0, prio2=0, isDW=0, hook1=None, hook2=None):
    return connector.optimizeOverlap(extCenters1, centers1, extCenters2,
                                     centers2, prio1, prio2, isDW, hook1, hook2)

#------------------------------------------------------------------------------
def maximizeBlankedCells(a, depth=2, dir=1, cellNName='cellN'):
    """Maximize the blanked region."""
    if isinstance(a[0], list):
        b = []
        for i in a:
            b.append(connector.maximizeBlankedCells(i, depth, dir, cellNName))
        return b
    else:
        return connector.maximizeBlankedCells(a, depth, dir, cellNName)

#-----------------------------------------------------------------------------
def _setHoleInterpolatedPoints(cellN, depth=2, dir=0, cellNName='cellN'):
    """Set interpolated points cellN=2 around cellN=0 points."""
    if depth == 0: return None
    if depth < 0:
        try: import Converter
        except: raise ImportError("_setHoleInterpolatedPoints: requires Converter module.")
    loc = 'nodes'
    if len(cellN) == 4:
        if cellN[3][-1]=='*': loc = 'centers'

    if loc == 'nodes':
        if depth < 0:
            Converter._initVars(cellN,'{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
            _getOversetHolesInterpNodes__(cellN, -depth, dir, cellNName)
            Converter._initVars(cellN,'{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
        else:
            _getOversetHolesInterpNodes__(cellN, depth, dir, cellNName)
    else: # unstructured with celln field at centers
        if depth < 0:
            Converter._initVars(cellN,'{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
            _getOversetHolesInterpCellCenters__(cellN, -depth, dir, cellNName)
            Converter._initVars(cellN,'{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
        else:
            _getOversetHolesInterpCellCenters__(cellN, depth, dir, cellNName)
    return None

def setHoleInterpolatedPoints(celln, depth=2, dir=0, cellNName='cellN', indices=None, BCField=None):
    """Set interpolated points cellN=2 around cellN=0 points."""
    if depth == 0: return celln
    try: import Converter
    except: raise ImportError("setHoleInterpolatedPoints: requires Converter module.")
    loc = 'nodes'
    if len(celln) == 4:
        if celln[3][-1]=='*': loc = 'centers'
    if loc == 'nodes':
        if depth < 0:
            celln = Converter.initVars(celln,
                                       '{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
            celln = getOversetHolesInterpNodes__(celln, -depth, dir, cellNName)
            celln = Converter.initVars(celln,
                                       '{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
        else: celln = getOversetHolesInterpNodes__(celln, depth, dir, cellNName)
    else: # unstructured with celln field at centers
        if depth < 0:
            celln = Converter.initVars(celln,
                                       '{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
            celln = getOversetHolesInterpCellCenters__(celln, -depth, dir, cellNName)
            celln = Converter.initVars(celln,
                                       '{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
        else: celln = getOversetHolesInterpCellCenters__(celln, depth, dir, cellNName, indices, BCField)
    return celln

#------------------------------------------------------------------------------
def blankCells(coords, cellnfields, body, blankingType=1, \
               delta=1.e-10, dim=3, masknot=0, tol=1.e-8, \
               XRaydim1=1000, XRaydim2=1000, cellNName='cellN'):
    """Blank cells in coords by a X-Ray mask defined by the body,
    within a distance delta.
    Usage: blankCells(coords, cellnfields, body, blankingType, delta, dim, maskNot, tol)"""
    try: import Converter as C
    except: raise ImportError("blankCells: requires Converter module.")
    # convert body to extended centers if structured and not node_in
    if blankingType != 0:
        # verify that all bodies are structured
        struct = 1
        for z in body:
            if len(z) != 5: struct = 0; break
        if struct == 1: body = C.node2ExtCenter(body)
    bodyt = []
    for z in body:
        z = C.convertArray2Tetra(z)
        bodyt.append(z)
    if blankingType == 2: # center_in: simply a node_in on centers
        coords = C.node2Center(coords); blankingType = 0
    return connector.blankCells(coords, cellnfields, bodyt, blankingType, \
                                delta, dim, masknot, tol, XRaydim1, XRaydim2, cellNName)
#------------------------------------------------------------------------------
# in place version: modifies cellnfields
def _blankCells(coords, cellnfields, body, blankingType=1, \
                delta=1.e-10, dim=3, masknot=0, tol=1.e-8, \
                XRaydim1=1000, XRaydim2=1000, cellNName='cellN'):
    """Blank cells in coords by a X-Ray mask defined by the body,
    within a distance delta. cellnfields is modified in place without copy.
    Usage: blankCells(coords, cellnfields, body, blankingType, delta, dim, maskNot, tol)"""
    try: import Converter as C
    except: raise ImportError("blankCells: requires Converter module.")
    # convert body to extended centers if structured and not node_in
    if blankingType != 0:
        # verify that all bodies are structured
        struct = 1
        for z in body:
            if len(z) != 5: struct = 0; break
        if struct == 1: body = C.node2ExtCenter(body)
    bodyt = []
    for z in body:
        z = C.convertArray2Tetra(z)
        bodyt.append(z)
    if blankingType == 2: # center_in: simply a node_in on centers
        coords = C.node2Center(coords); blankingType = 0

    return connector._blankCells(coords, cellnfields, bodyt, blankingType, \
                                 delta, dim, masknot, tol, XRaydim1, XRaydim2, cellNName)
#==============================================================================
# blankIntersectingCells
# IN: a: 3D structured mesh with wall orthogonal to k direction
# IN: tol: tolerance for intersecting cells
# OUT: returns the cellN, with cellN=0 for cells intersecting
# and the cells which are above.
#==============================================================================
def blankIntersectingCells(a, cellN, tol=1.e-12):
    """Blank intersecting cells in a zone.
    Usage: blankIntersectingCells(a, cellN, tol)"""
    res = connector.blankIntersectingCells(a, cellN, tol)
    return res

#==============================================================================
# blankCellsTetra
# IN: coords: 3D structured or unstructured mesh
# IN: meshT4: tetrahedral mesh of the mask
# IN: blankingType: 0(node_in) 1(cell_intersect) 2(cell_in)
# IN: tol: geometric tolerance
# OUT: returns the cellnfields, 0 for cells intersecting or inside the tet mesh
#==============================================================================
def blankCellsTetra(coords, cellnfields, meshT4, blankingType=1, tol=1.e-12, cellnval=0, overwrite=0, cellNName='cellN'):
    """Blank cells in coords (by setting the cellN to cellnval) falling inside a Tetra Mesh mask defined by meshT4.
    If overwrite is enabled (1), cells detected outside have a celln reset to 1.
    Usage: blankCellsTetra(coords, cellnfields, meshT4, connectT4, blankingType, tol, cellnval, overwrite)"""
    try:
        import Converter as C
        import Transform as T
        import Post as P
    except: raise ImportError("blankCellsTetra: requires Converter, Transform and Post  module.")
    cellnt = []
    maskSkin = P.exteriorFaces(meshT4)
    maskSkin = T.reorderAll(maskSkin, 1) #orient outward
    maskSkin = T.join(maskSkin)

    mask = connector.createTetraMask(meshT4, maskSkin, tol)

    for i in range(len(coords)):
        bt = blankingType
        if blankingType == 2: # center_in: simply a node_in on centers
            coords[i] = C.node2Center(coords[i])
            bt = 0
        cellnt.append(connector.blankCellsTetra(coords[i], cellnfields[i], mask, bt, cellnval, overwrite, cellNName))
    connector.deleteTetraMask(mask)
    return cellnt

#==============================================================================
# blankCellsTri
# IN: coords: 3D structured or unstructured mesh
# IN: meshT3: triangular mesh of the mask
# IN: blankingType: 0(node_in) 1(cell_intersect) 2(cell_in)
# IN: tol: geometric tolerance
# OUT: returns the cellnfields, 0 for cells intersecting or inside the tet mesh
#==============================================================================
def blankCellsTri(coords, cellnfields, meshT3, blankingType=1, tol=1.e-12,
                  cellnval=0, overwrite=0, cellNName='cellN'):
    """Blank cells in coords (by setting the cellN to cellnval) falling inside a Triangular surface mesh mask defined by meshT3.
    If overwrite is enabled (1), cells detected outside have a celln reset to 1.
    Usage: blankCellsTri(coords, cellnfields, meshT3, connectT4, blankingType, tol, cellnval, overwrite)"""
    try: import Converter as C; import Transform as T
    except: raise ImportError("blankCellsTetra: requires Converter module.")

    cellnt = []
    meshT3 = T.reorderAll(meshT3, 1) # orient outward
    meshT3 = T.join(meshT3)

    mask = connector.createTriMask(meshT3, tol)

    for i in range(len(coords)):
        bt = blankingType
        if blankingType == 2: # center_in: simply a node_in on centers
            coords[i] = C.node2Center(coords[i])
            bt = 0
        cellnt.append(connector.blankCellsTetra(coords[i], cellnfields[i], mask, bt, cellnval, overwrite, cellNName))
    connector.deleteTriMask(mask)
    return cellnt

def getIntersectingDomainsAABB(arrays, tol=1.e-10):
    """Return the intersection list of a list of bounding boxes."""
    return connector.getIntersectingDomainsAABB(arrays, tol)

#=============================================================================
# set cellN to 2 for nodes/cells in a neighborhood of depth nodes/cells
# for NGON arrays
#=============================================================================
def applyBCOverlapsNG__(a, faceList, depth, loc, val=2, cellNName='cellN'):
    if loc == 'nodes': locI = 0
    elif loc == 'centers': locI = 1
    else: raise ValueError("applyBCOverlapsUnstr: invalid location.")
    return connector.applyBCOverlapsNG(a, faceList, depth, locI, val, cellNName)

#=============================================================================
# set cellN to 2 for nodes/cells in a neighborhood of depth nodes/cells
# for structured arrays.
# in place
#=============================================================================
def _applyBCOverlapsStruct__(a, minIndex, maxIndex, depth, loc, val=2, cellNName='cellN'):
    if loc == 'nodes': locI = 0
    elif loc == 'centers': locI = 1
    else: raise ValueError("applyBCOverlapsStruct: invalid location.")
    return connector.applyBCOverlapStruct(a, minIndex, maxIndex, depth, locI, val, cellNName)

def applyBCOverlapsStruct__(a, minIndex, maxIndex, depth, loc, val=2, cellNName='cellN'):
    import Converter as C
    b = C.copy(a)
    _applyBCOverlapsStruct__(b, minIndex, maxIndex, depth, loc, val=val, cellNName=cellNName)
    return b

#=============================================================================
# Apply doubly defined boundary conditions for zone z
# cellN: cell nature field at centers: 2 if interpolated point
# listOfInterpZones: list of interpolation domains of z
# range: [imin,imax,jmin,jmax,kmin,kmax] of the BC of z
# depth: number of rows of interpolated points
# returns the modified cellN of the zone: if a BC point is not
# interpolable from listOfInterpZones, cellN is reset to 1
#=============================================================================
def setDoublyDefinedBC(z, cellN, listOfInterpZones, listOfCelln, winrange,
                       depth=1):
    """Set cellN to 2 to interpolated points of z near border of indices range defining a doubly defined BC.
    Usage: setDoublyDefinedBC(z, cellN, listOfInterpZones, listOfCelln, range, depth)"""
    try:
        import Converter as C
        listOfCelln = C.initVars(listOfCelln,'{cellN}=minimum(1.,{cellN})')
    except: pass
    return connector.setDoublyDefinedBC(z, cellN, listOfInterpZones,
                                        listOfCelln, winrange, depth)

#===============================================================================
# Compute donors and interpolation coefficients for a receiver zone
# from a list of donor zones. If several candidates are possible
# for interpolation from several donor zones, selects the
# interpolation cell with smallest volume
#
# IN: interpPts: points to interpolate (as list of arrays with coordinates modified for double wall)
# IN: zonesD: list of donor zones, containing cellN variable
# IN: order: order of interpolations
# IN: penalty=1: penalizes a donor cell in terms of volume if it is at the edge
# IN: nature=0: no vertex of the interpolation cell must have cellN=0
#     nature=1: all vertices of the interpolation cell must have cellN=1
# IN: extrap=0: no extrapolated points calculation
#     extrap=1: calculation and storage of possible extrapolated points
# IN: hook: hook on the ADT (not reconstructed in setInterpData), order must follow that of zonesD
# IN : interpDataType : 1 for ADT, 0 if donor are cartesian (optimized)
# OUT: res = [[rcvInd1],[donorInd1D],[donorType],[coefs],extrapInd1D, orphanInd1D]
#      res[0] list of 1D indices of interpolated/extrapolated points per donor zone (numbering from interpPts)
#      res[1] list of 1D indices of donor cells, whose storage is defined by the type
#      res[2] list of types for each interpolated point by a donor block
#      res[3] list per donor block of interpolation coefficients, as a 1D array
#      res[4] list per donor block of extrapolated points
#      res[5] orphanInd1D: 1D index of orphan points (number in interpPts)
#
#===============================================================================
def setInterpData__(interpPts, zonesD, order=2, penalty=1, extrap=1, nature=0, method='lagrangian', interpDataType=1, hook=None, dim=3):
    if method == 'lagrangian':
        if isinstance(interpPts[0], list): # list of arrays
            if interpPts[0][1].shape[1] >0: return connector.setInterpDataDW(interpPts, zonesD, order, nature, penalty, hook)
            else: return None
        else: # no list of arrays = no double wall
            if interpPts[1].shape[1]>0:
                if not isinstance(interpDataType,list): interpDataTypeL=[interpDataType]*len(zonesD)
                else: interpDataTypeL = interpDataType
                return connector.setInterpData(interpPts, zonesD, order, nature, penalty, extrap, interpDataTypeL, hook)
            else: return None

    elif method == 'leastsquares':
        if isinstance(interpPts[0], list): # list of arrays
            print('Warning: setInterpData: only 1st zone in 1st arg taken into account.')
            if interpPts[0][1].shape[1]>0:
                return connector.setInterpDataLS(interpPts[0], zonesD, order, nature, penalty, hook, dim)
            else: return None
        else:
            if interpPts[1].shape[1]>0:
                return connector.setInterpDataLS(interpPts, zonesD, order, nature, penalty, hook, dim)
            else: return None
    else:
        raise ValueError("setInterpData__: %s: not a valid interpolation method."%method)

#-----------------------------------------------------------------------------
# Compute and store interpolation coefficients for given points
# as 'NODE' form
# IN: rcvzonename: name of the interpolated zone containing interpPts
# IN: nir, njr dimensions (I,J) of the interpolated zone
# IN: interpPts: x,y,z at least of interpolated points
# IN: listOfInterpolationZones: list of interpolation domains - at nodes
# IN: listOfInterpolationCelln: list of associated celln fields - at centers
# IN: isEX=0 if processing at centers, isEX=1 if interpPts = EX points
# IN: zoneId: global number (in parallel) of the zone, starts at 1
# IN: cfMax, max authorized value (for sum of absolute coefficient values) for extrapolation
# OUT: interpolation coefficients
# OUT: indices for interpolated cells
# OUT: indices for interpolation cells
# OUT: indices of extended centers of interpolation cells
# OUT: volumes of interpolation cells
# OUT: interpTypes: type of interpolation applied to reconstruct the formula afterwards
# OUT: indices for orphan cells
# OUT: if isEX=1, direction for EX points
#-----------------------------------------------------------------------------
def setInterpolations__(rcvzonename,nir, njr, interpPts,
                        listOfInterpolationZones=[],
                        listOfInterpolationCelln=[], isEX=0, cfMax=30., zoneId=-1,
                        check=True):
    """Compute and store interpolation coefficients."""
    resInterp = connector.setInterpolations(nir, njr, interpPts,
                                            listOfInterpolationZones,
                                            listOfInterpolationCelln, isEX, zoneId, cfMax)
    # Summary:
    nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
    nbinterpolated = interpPts[0][1].shape[1]
    if len(resInterp[5])>0:
        for nozd in range(len(resInterp[6])): nbextrapolated += resInterp[6][nozd].size
        nborphan = resInterp[7][0].size
        nbinterpolated = nbinterpolated-nbextrapolated-nborphan
        if check: # sequential context
            if isEX == 0:
                print('Zone %s: interpolated=%d ; extrapolated=%d ; orphan=%d'%(rcvzonename, nbinterpolated, nbextrapolated, nborphan))
                if  nborphan>0: print('Warning: zone %s has %d orphan points !'%(rcvzonename, nborphan))
            else:
                print('Zone %s: EX interpolated=%d ; EX extrapolated=%d ; EX orphan=%d'%(rcvzonename, nbinterpolated, nbextrapolated, nborphan))
                if  nborphan>0: print('Warning: zone %s has %d EX orphan points !'%(rcvzonename, nborphan))
    return resInterp

#-----------------------------------------------------------------------------
# Write interpolation coefficients to files for elsA
# IN: ntotZones: total number of zones
# IN: listRcvId: list of receiver zone IDs per donor block
# IN: listCellIndicesRcv: list of arrays of interpolated cell indices, per donor block
# IN: listOfDirectionEX: list of arrays of EX point indirection, per donor block
# IN: listCellIndicesDonor: list of arrays of donor cell indices, per donor block
# IN: listInterpolantsDonor: list of arrays of interpolation coefficients, per donor block
# IN: listInterpTypes: list of interpolation types (100 classic, 102 and 103 for azimuthal periodicity)
# IN: listCellN: list of cellNatureField fields, per donor block
# IN: listDonorDim: list of donor block dimensions
# IN: nbInterpCellsForDonor: list of number of interpolated cells in a receiver block, per donor block
# IN: PrefixFile: output filename for the solver
# IN: isEX=0 if processing at centers, isEX=1 if interpPts = EX points
# IN: solver: solver name for which the file is generated. 1: elsA, 2: Cassiopee
# IN: nGhostCells: number of ghost cells for elsA data writing
#-----------------------------------------------------------------------------
def writeCoefs(ntotZones,listRcvId,listCellIndicesRcv,listOfDirectionEX,listCellIndicesDonor,listInterpolantsDonor,
               listInterpTypes, listCellN, listDonorDim, nbInterpCellsForDonor,
               PrefixFile, isEX=0, solver=1, nGhostCells=2):
    """Write interpolation coefficients in files for elsA."""
    connector.writeCoefs(ntotZones,listRcvId,listCellIndicesRcv,listOfDirectionEX,listCellIndicesDonor,listInterpolantsDonor, listInterpTypes, listCellN, listDonorDim, nbInterpCellsForDonor, PrefixFile, isEX, solver, nGhostCells)

#-----------------------------------------------------------------------------
# Returns the modified celln of interpolated cells in the vicinity of masked points
# IN: zc: contains at least celln (celln=0 for interpolated points)
# IN: depth: number of rows of interpolated cells
#-----------------------------------------------------------------------------
def getOversetHolesInterpCellCenters__(zc, depth=2, dir=0, cellNName='cellN',
                                       indices=None, BCField=None):
    """Set cellN=2 for the fringe of interpolated cells around cells of celln
    equal to 0."""
    return connector.getOversetHolesInterpCellCenters(zc, depth, dir, cellNName, indices, BCField)

# version getFromArray2: works only for structured
def _getOversetHolesInterpCellCenters__(zc, depth=2, dir=0, cellNName='cellN'):
    """Set cellN=2 for the fringe of interpolated cells around cells of celln
    equal to 0."""
    return connector._getOversetHolesInterpCellCenters(zc, depth, dir, cellNName)

#-----------------------------------------------------------------------------
# Returns the modified celln of interpolated nodes in the vicinity of masked points
# IN: zc: contains at least celln (celln=0 for interpolated points)
# IN: depth: number of rows of interpolated points
#-----------------------------------------------------------------------------
def getOversetHolesInterpNodes__(z, depth=2, dir=0, cellNName='cellN'):
    """Set cellN=2 for the fringe of interpolated nodes around nodes of celln
    equal to 0."""
    return connector.getOversetHolesInterpNodes(z, depth, dir, cellNName)

# version getFromArray2: works only for structured
def _getOversetHolesInterpNodes__(z, depth=2, dir=0, cellNName='cellN'):
    """Set cellN=2 for the fringe of interpolated nodes around nodes of celln
    equal to 0."""
    return connector._getOversetHolesInterpNodes(z, depth, dir, cellNName)

#------------------------------------------------------------------------------
# Returns the coordinates of EX points and indices in the initial mesh
# at centers of adjacent interpolated cells that constitute the fringe
# of EX points around masked cells, as well as those near boundaries
# IN: z: mesh at nodes
# IN: celln: celln at centers with 0 and 2 for masked and interpolated points
# OUT: infoEX: x,y,z,indcell1,indcell2, nodemin,dirEX of EX points: 'NODE' elements
# indcell1 corresponds to the interpolated center from which the EX point originates
# indcell2 corresponds to the second center to which the EX point interface is attached
# nodemin: index of the node with smallest index of the interface whose center is EX
# dirEX: direction of the interface whose center is EX (1, 2, or 3)
#------------------------------------------------------------------------------
def getEXPoints__(z, celln):
    """Return the coordinates of the EX points."""
    return connector.getEXPoints(z, celln)

#------------------------------------------------------------------------------
def maskXRay__(body, delta=0., dim=3, isNot=0, tol=1.e-8):
    """Create the pierce points of a X-Ray mask defined by body."""
    try: import Converter as C; import Transform as T
    except: raise ImportError("maskXRay: requires Converter and Transform modules.")
    body = C.convertArray2Tetra(body)
    body = T.join(body)
    return connector.maskXRay([body], delta, dim, isNot, tol)

# cellN modification
def _modCellN1(a, cellNName='cellN'):
    """Change cellN: 0->-1, 2->1"""
    if isinstance(a[0], list):
        for i in a:
            connector._modCellN1(i, cellNName)
    else:
        return connector._modCellN1(a, cellNName)
    return None

# cellN modification
def _modCellN2(a, cellNName='cellN'):
    """Change cellN: -1->0"""
    if isinstance(a[0], list):
        for i in a:
            connector._modCellN2(i, cellNName)
    else:
        return connector._modCellN2(a, cellNName)
    return None

#------------------------------------------------------------------------------
# IN: x,y,z, cellN located at the same place
# OUT: 'NODE' array with 'x,y,z,indcell'
#      with indcell the index of the corresponding node of z such that cellN(indcell)=2
#-----------------------------------------------------------------------------
def getInterpolatedPoints__(a):
    if isinstance(a[0], list):
        b = []
        for i in a: b.append(connector.getInterpolatedPoints(i))
        return b
    else:
        return connector.getInterpolatedPoints(a)

#------------------------------------------------------------------------------
# New changeWall algorithm, without double wall tolerance
# IN: z: mesh to project
# IN: firstWallPoints: indices of the first points of the mesh near its wall
#      with dir1, dir2, dir3, hmax of the first centers near the walls
# IN: projectionSurfaces: list of TRI surfaces corresponding to the projection surfaces
# OUT: z modified for points with cellN=2
#------------------------------------------------------------------------------
def changeWall__(z, firstWallPoints, projectionSurfaces, planarTol=0.):
    if projectionSurfaces == [] or firstWallPoints == []: return z
    else: return connector.changeWall(z, firstWallPoints, projectionSurfaces, planarTol)

#-----------------------------------------------------------------------------
# New changeWall algorithm, without double wall tolerance
# IN: EXPts: (x,y,z,indcell1,indcell2,nodemin,EXdir) for EX points
# IN: zc: mesh at centers of the zone containing EX points to project
# IN: zn: mesh at nodes of the zone containing EX points to project
# IN: firstWallCenters: indices at centers, dir1, dir2, dir3, hmax of the first centers near the walls
# IN: projectionSurfaces: list of TRI surfaces corresponding to the projection surfaces
# OUT: zc modified for points with cellN=2
#-----------------------------------------------------------------------------
def changeWallEX__(EXPts, zc, zn, firstWallCenters, projectionSurfaces, planarTol=0.):
    if projectionSurfaces == [] or firstWallCenters == []: return zc
    else: return connector.changeWallEX(EXPts, zn, zc, firstWallCenters, projectionSurfaces, planarTol)

#----------------------------------------------------------------------------------------------
# To determine the double wall projection boundaries, in the case where boundaries
# are split into sub-windows, we extend the sub-window to extended centers (input
# array) which we modify locally by taking the midpoint with the interior point
# IN: array: 2D structured boundary at extended centers
# IN: iminL, imaxL, jminL, jmaxL: indices of the sub-window in the topological window
# OUT: array modified at the edges (if iminL > 1 then we take the midpoint with iminL+1 for example)
#----------------------------------------------------------------------------------------------
def modifyBorders__(array,iminL,imaxL,jminL,jmaxL):
    return connector.modifyBorders(array, iminL,imaxL,jminL,jmaxL)
#-----------------------------------------------------------------------------
# Compute Chimera transfers
# IN: cellRcv: numpy array of interpolated cell indices
# IN: cellDonor: numpy array of donor cell indices
# IN: interpDonor: numpy array of interpolation coefficients
# IN: interpType: numpy array of interpolation types (102 hardcoded for elsA)
# IN: lRcvArrays: list of fields on which the transfer is applied (array)
# IN: lDonorArrays: list of fields of the interpolation domain (array)
#-----------------------------------------------------------------------------
def chimeraTransfer(cellRcv, cellDonor, interpType,
                    interpDonor, lRcvArrays, lDonorArrays):
    return connector.chimeraTransfer(cellRcv, cellDonor, interpType,
                                     interpDonor,lRcvArrays,lDonorArrays)

#------------------------------------------------------------------------------
# Compute interpolations from fields given in donorFields
# rcvFields: fields to interpolate as Converter.array
# donorFields: donor fields as arrays in the Converter sense
# indicesRcv: indices of points to interpolate in zonesR, as numpy array1D,
# indicesDnr: indices of donors, to retrieve according to type, as numpy array1D
# donorType: interpolation type as numpy array 1D
# coefs: interpolation coefficients as numpy array 1D
# variables: if None: transfer of all variables except cellN
#            if ['a','b']: transfer of variables ['a','b'] only
#------------------------------------------------------------------------------
def setInterpTransfers(rcvFields, donorFields, indicesRcv, indicesDnr,
                       donorType, coefs, variables=[]):
    return connector.setInterpTransfers(rcvFields, donorFields, variables, indicesRcv,
                                        indicesDnr, donorType, coefs)

#------------------------------------------------------------------------------
# Compute IBC interpolations from fields given in donorFields
# rcvFields: fields to interpolate as Converter.array
# donorFields: donor fields as arrays in the Converter sense
# indicesRcv: indices of points to interpolate in zonesR, as numpy array1D,
# indicesDnr: indices of donors, to retrieve according to type, as numpy array1D
# donorType: interpolation type as numpy array 1D
# coefs: interpolation coefficients as numpy array 1D
# bcType: 'slip','noslip'
# varType: 1=[ro,rou,rov,row,roE], 2=[ro,u,v,w,t], 3=[ro,u,v,w,p]
#------------------------------------------------------------------------------
def setIBCTransfers(rcvFields, donorFields, indicesRcv, indicesDnr, donorType,
                    coefs, xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI,
                    bcType='slip',
                    variables=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                    varType=1):
    return connector.setIBCTransfers(rcvFields, donorFields, variables,
                                     indicesRcv, indicesDnr, donorType, coefs,
                                     xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI,
                                     bcType, varType)

#------------------------------------------------------------------------------
# Compute interpolations from fields given in donorFields
# donorFields: donor fields as arrays in the Converter sense
# indicesDnr: indices of donors, to retrieve according to type, as numpy array1D
# donorType: interpolation type as numpy array 1D
# coefs: interpolation coefficients as numpy array 1D
# Returns an array with interpolated fields
#------------------------------------------------------------------------------
def setInterpTransfersD(donorFields, indicesDnr, donorType, coefs):
    return connector.setInterpTransfersD(donorFields, indicesDnr, donorType, coefs)

#------------------------------------------------------------------------------
# Compute IBC interpolations from fields given in donorFields
# donorFields: donor fields as arrays in the Converter sense
#              these fields are ordered according to varType
# indicesDnr: indices of donors, to retrieve according to type, as numpy array1D
# donorType: interpolation type as numpy array 1D
# coefs: interpolation coefficients as numpy array 1D
# bcType: BC type (slip, noslip)
# varType: 1=[ro,rou,rov,row,roE], 2=[ro,u,v,w,t], 3=[ro,u,v,w,p]
#------------------------------------------------------------------------------
def setIBCTransfersD(donorFields, indicesDnr, donorType, coefs,
                     xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI,
                     bcType='slip', varType=1):
    return connector.setIBCTransfersD(donorFields, indicesDnr, donorType, coefs,
                                      xPC,yPC,zPC,xPW,yPW,zPW,xPI,yPI,zPI,
                                      bcType, varType)
