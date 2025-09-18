"""Find connectivity in grids.
"""
__version__ = '4.1'
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
# Retourne [ [[noz1,noz2],[imin1,imax1,...],[imin2,imax2,...],trirac] ]
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
# Verifie si les deux windows en raccord match ont des normales opposees
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
    else: # non structure avec champ celln en centres
        if depth < 0:
            Converter._initVars(cellN,'{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
            _getOversetHolesInterpCellCenters__(cellN, -depth, dir, cellNName)
            Converter._initVars(cellN,'{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
        else:
            _getOversetHolesInterpCellCenters__(cellN, depth, dir, cellNName)
    return None

def setHoleInterpolatedPoints(celln, depth=2, dir=0, cellNName='cellN'):
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
    else: # non structure avec champ celln en centres
        if depth < 0:
            celln = Converter.initVars(celln,
                                       '{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
            celln = getOversetHolesInterpCellCenters__(celln, -depth, dir, cellNName)
            celln = Converter.initVars(celln,
                                       '{%s} = 1-{%s}+({%s}>1.5)*3'%(cellNName, cellNName, cellNName))
        else: celln = getOversetHolesInterpCellCenters__(celln, depth, dir, cellNName)
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
    # passe body en centres etendus si structure et pas node_in
    if blankingType != 0:
        # verif que ts les body sont structures
        struct = 1
        for z in body:
            if len(z) != 5: struct = 0; break
        if struct == 1: body = C.node2ExtCenter(body)
    bodyt = []
    for z in body:
        z = C.convertArray2Tetra(z)
        bodyt.append(z)
    if blankingType == 2: # center_in: simplement un node_in sur les centres
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
    # passe body en centres etendus si structure et pas node_in
    if blankingType != 0:
        # verif que ts les body sont structures
        struct = 1
        for z in body:
            if len(z) != 5: struct = 0; break
        if struct == 1: body = C.node2ExtCenter(body)
    bodyt = []
    for z in body:
        z = C.convertArray2Tetra(z)
        bodyt.append(z)
    if blankingType == 2: # center_in: simplement un node_in sur les centres
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
        #print 'coords : %d / %d' %(i+1, len(coords))
        bt = blankingType
        if blankingType == 2: # center_in: simplement un node_in sur les centres
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
        #print('coords : %d / %d' %(i+1, len(coords)))
        bt = blankingType
        if blankingType == 2: # center_in: simplement un node_in sur les centres
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
# Application des conditions aux limites doublement definies pour une zone z
# cellN: cellnaturefield aux centres: 2 si point interpole
# listOfInterpZones: liste des domaines d interpolation de z
# range: [imin,imax,jmin,jmax,kmin,kmax] de la CL de z
# depth: nb de rangees de points interpoles
# retourne le cellN modifie de la zone: si un point de la CL n est pas
# interpolable a partir de la listOfInterpZones, le cellN est remis a 1
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
# Calcul les donneurs et coefs d interpolation pour une zone receveuse
# a partir d'une liste de zones donneuses. Si plusieurs candidats possibles
# pour l'interpolation issus de plusieurs zones donneuses, selectionne la
# cellule d interpolation de plus petit volume
#
# IN: interpPts: pts a interpoler (sous forme de liste d arrays avec coordonnees modifiees pour le double wall)
# IN: zonesD: liste de zones donneuses, contenant la variable cellN
# IN: order: ordre des interpolations
# IN: penalty=1: penalise une cellule donneuse en terme de volume si elle est au bord
# IN: nature=0: aucun sommet de la cellule d'interpolation ne doit etre avec un cellN=0
#     nature=1: toutes les sommets de la cellule d'interpolation doivent etre de cellN=1
# IN: extrap=0: pas de calcul des point extrapoles
#     extrap=1: calcul et stockage des eventuels pts extrapoles
# IN: hook: hook sur l'adt (pas reconstruit dans setInterpData), l'ordre doit suivre celui de zonesD
# IN : interpDataType : 1 for ADT, 0 if donor are cartesian (optimized)
# OUT: res = [[rcvInd1],[donorInd1D],[donorType],[coefs],extrapInd1D, orphanInd1D]
#      res[0] liste des indices 1D des pts interpoles/extrapoles par zone donneuse (numerotation de interpPts)
#      res[1] liste des indices 1D des molecules donneuses, dont le stockage est defini par le type
#      res[2] liste des types pour chq pt interpole par un bloc donneur
#      res[3] liste par bloc donneur des coefs d interpolation, sous forme d un tableau 1D
#      res[4] liste par bloc donneur des pts extrapoles
#      res[5] orphanInd1D: indice 1D des points orphelins (numero ds interpPts)
#
#===============================================================================
def setInterpData__(interpPts, zonesD, order=2, penalty=1, extrap=1, nature=0, method='lagrangian', interpDataType=1, hook=None, dim=3):
    if method == 'lagrangian':
        if isinstance(interpPts[0], list): # liste d'arrays
            if interpPts[0][1].shape[1] >0: return connector.setInterpDataDW(interpPts, zonesD, order, nature, penalty, hook)
            else: return None
        else: # pas de liste d'arrays = pas de double wall
            if interpPts[1].shape[1]>0:
                if not isinstance(interpDataType,list): interpDataTypeL=[interpDataType]*len(zonesD)
                else: interpDataTypeL = interpDataType
                return connector.setInterpData(interpPts, zonesD, order, nature, penalty, extrap, interpDataTypeL, hook)
            else: return None

    elif method == 'leastsquares':
        if isinstance(interpPts[0], list): # liste d'arrays
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
# Calcul et stockage des coefficients d'interpolation pour les pts donnes
# sous forme de 'NODE'
# IN: rcvzonename: nom de la zone interpolee contenant les pts interpPts
# IN: nir, njr dimensions (I,J) de la zone interpolee
# IN: interpPts: x,y,z au moins des pts interpoles
# IN: listOfInterpolationZones: liste des domaines d interpolation - en noeuds
# IN: listOfInterpolationCelln: liste des champs celln associes - en centres
# IN: isEX=0 si traitement aux centres, isEX=1 si interpPts = pts EX
# IN: zoneId: numero global (en parallele) de la zone, commence a 1
# IN: cfMax, valeur max. authorisee (pour la somme des valeurs absolues des coefs) pour l extrapolation
# OUT: coefficients d'interpolation
# OUT: indices pour les cellules interpolees
# OUT: indices pour les cellules d interpolation
# OUT: indices des centres etendus des cellules d interpolation
# OUT: volumes des cellules d interpolation
# OUT: interpTypes: type d'interpolation appliquee pour reconstruire la formule ensuite
# OUT: indices pour les cellules orphelines
# OUT: si isEX=1, direction pour les points EX
#-----------------------------------------------------------------------------
def setInterpolations__(rcvzonename,nir, njr, interpPts,
                        listOfInterpolationZones=[],
                        listOfInterpolationCelln=[], isEX=0, cfMax=30., zoneId=-1,
                        check=True):
    """Compute and store interpolation coefficients."""
    resInterp = connector.setInterpolations(nir, njr, interpPts,
                                            listOfInterpolationZones,
                                            listOfInterpolationCelln, isEX, zoneId, cfMax)
    # Bilan :
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
# Ecriture des coefficients d'interpolation dans des fichiers pour elsA
# IN: ntotZones: nombre total de zones
# IN: listRcvId: liste des Id des zones interpolees par bloc donneur
# IN: listCellIndicesRcv: liste des tableaux d'indices des cellules interpolees, par bloc donneur
# IN: listOfDirectionEX: liste des tableaux d'indirection des points EX, par bloc donneur
# IN: listCellIndicesDonor: liste des tableaux d'indices des cellules donneuses, par bloc donneur
# IN: listInterpolantsDonor: liste des tableaux de coefficients d'interpolation, par bloc donneur
# IN: listInterpTypes: liste des types d interpolation (100 classique, 102 et 103 pour le periodique en azimut)
# IN: listCellN: liste des champs cellNatureField, par bloc donneur
# IN: listDonorDim: liste des dimensions des blocs donneurs
# IN: nbInterpCellsForDonor: liste du nombre de cellules interpolees dans un bloc interpolee, par bloc donneur
# IN: PrefixFile : nom du fichier de sortie pour le solveur.
# IN: isEX=0 si traitement aux centres, isEX=1 si interpPts = pts EX
# IN: solver : nom du solveur pour lequel le fichier est genere. 1: elsA, 2: Cassiopee
# IN: nGhostCells: nb de cellules fictives pour l ecriture des donnees pour elsA
#-----------------------------------------------------------------------------
def writeCoefs(ntotZones,listRcvId,listCellIndicesRcv,listOfDirectionEX,listCellIndicesDonor,listInterpolantsDonor,
               listInterpTypes, listCellN, listDonorDim, nbInterpCellsForDonor,
               PrefixFile, isEX=0, solver=1, nGhostCells=2):
    """Write interpolation coefficients in files for elsA."""
    connector.writeCoefs(ntotZones,listRcvId,listCellIndicesRcv,listOfDirectionEX,listCellIndicesDonor,listInterpolantsDonor, listInterpTypes, listCellN, listDonorDim, nbInterpCellsForDonor, PrefixFile, isEX, solver, nGhostCells)

#-----------------------------------------------------------------------------
# Retourne le celln modifie des cellules interpolees au voisinage des pts masques
# IN: zc: contient au moins le celln (celln=0 pour les pts interpoles)
# IN: depth: nb de rangees de cellules interpolees
#-----------------------------------------------------------------------------
def getOversetHolesInterpCellCenters__(zc, depth=2, dir=0, cellNName='cellN'):
    """Set cellN=2 for the fringe of interpolated cells around cells of celln
    equal to 0."""
    return connector.getOversetHolesInterpCellCenters(zc, depth, dir, cellNName)

# version getFromArray2: ne marche qu en structure
def _getOversetHolesInterpCellCenters__(zc, depth=2, dir=0, cellNName='cellN'):
    """Set cellN=2 for the fringe of interpolated cells around cells of celln
    equal to 0."""
    return connector._getOversetHolesInterpCellCenters(zc, depth, dir, cellNName)
#-----------------------------------------------------------------------------
# Retourne le celln modifie des noeuds interpoles au voisinage des pts masques
# IN: zc: contient au moins le celln (celln=0 pour les pts interpoles)
# IN: depth: nb de rangees de pts interpoles
#-----------------------------------------------------------------------------
def getOversetHolesInterpNodes__(z, depth=2, dir=0, cellNName='cellN'):
    """Set cellN=2 for the fringe of interpolated nodes around nodes of celln
    equal to 0."""
    return connector.getOversetHolesInterpNodes(z, depth, dir, cellNName)

# version getFromArray2: ne marche qu'en structure
def _getOversetHolesInterpNodes__(z, depth=2, dir=0, cellNName='cellN'):
    """Set cellN=2 for the fringe of interpolated nodes around nodes of celln
    equal to 0."""
    return connector._getOversetHolesInterpNodes(z, depth, dir, cellNName)

#------------------------------------------------------------------------------
# Retourne les coordonnees des pts EX et les indices dans le maillage initial
# en centres des cellules interpolees adjacentes qui constituent la frange
# de pts EX autour des cellules masquees, ainsi que ceux pres des frontieres
# IN: z: maillage en noeuds
# IN: celln: celln en centres avec 0 et 2 pour les pts masques et interpoles
# OUT: infoEX: x,y,z,indcell1,indcell2, nodemin,dirEX des pts EX : elts 'NODE'
# indcell1 correspond au centre interpole d ou provient le point EX
# indcell2 correspond au deuxieme centre auquel est rattachee l interface du pt EX
# nodemin: indice du noeud d'indice le plus petit de l'interface dont EX est le centre
# dirEX: direction de l'interface dont EX est le centre (1,2 ou 3)
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
# IN: x,y,z, cellN localises au meme endroit
# OUT: array 'NODE' avec 'x,y,z,indcell'
#      avec indcell l'indice du noeud correspondant de z tq cellN(indcell)=2
#-----------------------------------------------------------------------------
def getInterpolatedPoints__(a):
    if isinstance(a[0], list):
        b = []
        for i in a: b.append(connector.getInterpolatedPoints(i))
        return b
    else:
        return connector.getInterpolatedPoints(a)

#------------------------------------------------------------------------------
# Nouvel algo de changeWall, sans tolerance double wall
# IN: z: maillage a projeter
# IN: firstWallPoints: indices des premiers pts du maillage pres de sa paroi
#      avec dir1, dir2, dir3, hmax des premiers centres pres des parois
# IN: projectionSurfaces: liste des surfaces TRI correspondant aux surfaces de projection
# OUT: z modifie pour les pts de cellN=2
#------------------------------------------------------------------------------
def changeWall__(z, firstWallPoints, projectionSurfaces, planarTol=0.):
    if projectionSurfaces == [] or firstWallPoints == []: return z
    else: return connector.changeWall(z, firstWallPoints, projectionSurfaces, planarTol)

#-----------------------------------------------------------------------------
# Nouvel algo de changeWall, sans tolerance double wall
# IN: EXPts: (x,y,z,indcell1,indcell2,nodemin,EXdir) pour les pts EX
# IN: zc: maillage en centres de la zone contenant les pts EX a projeter
# IN: zn: maillage en noeuds de la zone contenant les pts EX a projeter
# IN: firstWallCenters: indices en centres, dir1, dir2, dir3, hmax des premiers centres pres des parois
# IN: projectionSurfaces: liste des surfaces TRI correspondant aux surfaces de projection
# OUT: zc modifie pour les pts de cellN=2
#-----------------------------------------------------------------------------
def changeWallEX__(EXPts, zc, zn, firstWallCenters, projectionSurfaces, planarTol=0.):
    if projectionSurfaces == [] or firstWallCenters == []: return zc
    else: return connector.changeWallEX(EXPts, zn, zc, firstWallCenters, projectionSurfaces, planarTol)

#----------------------------------------------------------------------------------------------
# Pour determiner les frontieres de projection double wall, dans le cas ou les frontieres
# sont decoupees en sous fenetres, on etend la sous-fenetre en centres etendus (input
# array) qu on modifie localement en prenant le pt milieu avec le point interieur
# IN: array : frontiere structuree 2D en centres etendus
# IN: iminL, imaxL, jminL, jmaxL : indices de la sous fenetre dans la fenetre topologique
# OUT : array modifie aux bords (si iminL > 1 alors on prend le pt milieu avec iminL+1 par ex)
#----------------------------------------------------------------------------------------------
def modifyBorders__(array,iminL,imaxL,jminL,jmaxL):
    return connector.modifyBorders(array, iminL,imaxL,jminL,jmaxL)
#-----------------------------------------------------------------------------
# Calcul des transferts Chimere
# IN: cellRcv     : tableau numpy des indices des cellules interpolees
# IN: cellDonor   : tableau numpy des indices des cellules donneuses
# IN: interpDonor : tableau numpy des coefficients d'interpolation
# IN: interpType  : tableau numpy du type des interpolations 102 cablees pour elsA
# IN: lRcvArrays  : liste des champs sur lesquels le transfert est applique (array)
# IN: lDonorArrays: liste des champs du domaine d'interpolation (array)
#-----------------------------------------------------------------------------
def chimeraTransfer(cellRcv, cellDonor, interpType,
                    interpDonor, lRcvArrays, lDonorArrays):
    return connector.chimeraTransfer(cellRcv, cellDonor, interpType,
                                     interpDonor,lRcvArrays,lDonorArrays)

#------------------------------------------------------------------------------
# Calcul des interpolations a partir de champs donnes dans donorFields
# rcvFields: champs a interpoler sous forme de Converter.array
# donorFields: champs donneurs sous forme d'arrays au sens de Converter
# indicesRcv: indices des pts a interpoler dans zonesR, sous forme de numpy array1D,
# indicesDnr: indices des donneurs, a retrouver selon le type, sous forme de numpy array1D
# donorType: type d'interpolation sous forme de numpy array 1D
# coefs: coefficients d'interpolation sous forme de numpy array 1D
# variables: si None: transfert de toutes les variables sauf le cellN
#            si ['a','b']: transfert des variables ['a','b'] uniquement
#------------------------------------------------------------------------------
def setInterpTransfers(rcvFields, donorFields, indicesRcv, indicesDnr,
                       donorType, coefs, variables=[]):
    return connector.setInterpTransfers(rcvFields, donorFields, variables, indicesRcv,
                                        indicesDnr, donorType, coefs)

#------------------------------------------------------------------------------
# Calcul des interpolations IBC a partir de champs donnes dans donorFields
# rcvFields: champs a interpoler sous forme de Converter.array
# donorFields: champs donneurs sous forme d'arrays au sens de Converter
# indicesRcv: indices des pts a interpoler dans zonesR, sous forme de numpy array1D,
# indicesDnr: indices des donneurs, a retrouver selon le type, sous forme de numpy array1D
# donorType: type d'interpolation sous forme de numpy array 1D
# coefs: coefficients d'interpolation sous forme de numpy array 1D
# bcType : 'slip','noslip'
# varType : 1=[ro,rou,rov,row,roE], 2=[ro,u,v,w,t], 3=[ro,u,v,w,p]
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
# Calcul des interpolations a partir de champs donnes dans donorFields
# donorFields: champs donneurs sous forme d'arrays au sens de Converter
# indicesDnr: indices des donneurs, a retrouver selon le type, sous forme de numpy array1D
# donorType: type d'interpolation sous forme de numpy array 1D
# coefs: coefficients d'interpolation sous forme de numpy array 1D
# Retourne un array avec les champs interpoles
#------------------------------------------------------------------------------
def setInterpTransfersD(donorFields, indicesDnr, donorType, coefs):
    return connector.setInterpTransfersD(donorFields, indicesDnr, donorType, coefs)

#------------------------------------------------------------------------------
# Calcul des interpolations IBC a partir des champs donnes dans donorFields
# donorFields: champs donneurs sous forme d'arrays au sens de Converter
#              ces champs sont ordonnes dans le sens de varType
# indicesDnr: indices des donneurs, a retrouver selon le type, sous forme de numpy array1D
# donorType: type d'interpolation sous forme de numpy array 1D
# coefs: coefficients d'interpolation sous forme de numpy array 1D
# bcType : type de CL (slip, noslip)
# varType : 1=[ro,rou,rov,row,roE], 2=[ro,u,v,w,t], 3=[ro,u,v,w,p]
#------------------------------------------------------------------------------
def setIBCTransfersD(donorFields, indicesDnr, donorType, coefs,
                     xPC, yPC, zPC, xPW, yPW, zPW, xPI, yPI, zPI,
                     bcType='slip', varType=1):
    return connector.setIBCTransfersD(donorFields, indicesDnr, donorType, coefs,
                                      xPC,yPC,zPC,xPW,yPW,zPW,xPI,yPI,zPI,
                                      bcType, varType)
