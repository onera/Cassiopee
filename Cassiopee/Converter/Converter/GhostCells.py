# - ghostcells -
# Un module d'ajout/suppression de cellules fictives dans les arbres python
try: range = xrange
except: pass

from . import Internal
from . import Converter
from . import converter
import numpy
from . import PyTree
import math
# global definition: coordinates name
coordinate = ['CoordinateX', 'CoordinateY', 'CoordinateZ']

#===============================================================================
# Add ghost cells in a pyTree or a zone
# Returns a new zone with ghost cells
# IN: t: top tree
# IN: b: tree/basis/zone to modify (defined in t)
# IN: d: number of ghost cells to add
# IN: modified: container names or the list of variables to be modified
# If modified != all: only the provided field is filled with ghost cells 
#                     -> flow field and zone dimensions are not consistent 
#                      
# adaptBCs = 0: BCs are not modified
#          = 1: Physical BCs/BCMatch/BCOverlap ranges are modified 
#               BCMatch ranges are not relevant geometrically (graphically)
# fillCorner  : method to fill edges and corners  
#          = 1: edges and corners are filled (grid coordinates+flow solution) according to 
#               the grid connectivity -> geometrically, the corners and edges can be wrong
#          = 0: neighbouring vectors are extrapolated to build edge cells, no filling with flow field
# fillCorner = 1 : non implemented for nearmatch
#===============================================================================
def addGhostCells(t, b, d, adaptBCs=0, modified=[], fillCorner=1):
    bp = Internal.copyRef(b)
    _addGhostCells(t, bp, d, adaptBCs, modified, fillCorner)
    return bp

def _addGhostCells(t, b, d, adaptBCs=0, modified=[], fillCorner=1):
    if modified == []: # Default value
        modified = [Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__]

    # Multi-base: if t contains several zones with the same name, only the first zone found is kept
    tp = Internal.copyRef(t)
    nodesRef = getFirstReferencedZones__(tp) # instead of tp
    stdNode = Internal.isStdNode(b)
    zones = Internal.getZones(b)

    # =========================================================================
    # Algorithm:
    # 3 loops are done
    # 1st loop on zones:
    # ---------
    #    - initialization with 0-order extrapolation
    #    - fill ghost cells (except corners cells) for join border
    # 2nd loop (only on joins which need it):
    # ---------
    #    - fill edge ghost cells for join border (using opposite ghost grid)
    # 3rd loop (only on joins which need it):
    # ---------
    #    - fill corner ghost cells for join border (using opposite ghost grid)
    # =========================================================================
    # In first loop, a dictionary dictjoins is built: 
    # it contains the valid zones and the name of the related donor zones
    # the corresponding list of valid connectivities is defined by 'validjoins' 
    # This dictionary will be used in second loop
    validjoins = []; dictjoins = {}; indexjoin = 0

    c = 0
    isghost = [True]*len(zones) # pour le cas NGON avec modified=uniquement sur les centres
    for zp in zones:
        dimZone = Internal.getZoneDim(zp)
        if dimZone[0] == 'Unstructured':
            typeElt = dimZone[3]
            if typeElt == 'NGON':
                # CB
                #isghost[c] = _addGhostCellsNGON__(zp, b, d, isghost[c], stdNode, modified)
                pass
            else: 
                #print('Warning: addGhostCells: not yet implemented for BE.')
                pass # other elt types
        else : # Structured
            validjoins,dictjoins,indexjoin = _addGhostCellsStruct__(zp, b, d, modified, nodesRef, validjoins, dictjoins, indexjoin)
        c += 1

    #---------------------------------------------------------------------------
    # remplissage des coins et aretes
    # geometriquement, les coins et aretes peuvent ne pas etre bien representes
    # en presence de points triples notamment
    # les champs sont remplis
    #---------------------------------------------------------------------------
    if fillCorner == 1: 
        # Second loop only for joins (deals with ghost cells corner)
        _fillGhostCellCorners__(b, modified, d, validjoins, dictjoins, loop=2)
        # Third loop only for joins (deals with last ghost cells corner)
        _fillGhostCellCorners__(b, modified, d, validjoins, dictjoins, loop=3)

    # -----------------------------------------------------------------------------
    # Modify BCMatch/BCOverlap and ZoneBCs
    #---------------------------------------------------------------------------
    if adaptBCs == 1:
        for zp in zones:
            dimZone = Internal.getZoneDim(zp)
            _adaptBCsForZone(zp, dimZone, d, nodesRef)

    #-------------------------------
    # put ghost dimensions for zone
    #-------------------------------   
    for zp in zones: _updateZoneDim__(zp, d)

    #------------------------------------------------------------
    # add Rind for zone
    #------------------------------------------------------------
    c = 0
    for zp in zones:
        if isghost[c]:
            rindplanes = [d,d,d,d,d,d]
            Internal.createChild(zp, 'ZoneRind', 'Rind_t', 
                                 value=rindplanes, children=[])
        c += 1
    #------------------------------------------------------------
    # Corners: geometrical extrapolation from neighbouring cells
    #------------------------------------------------------------
    if fillCorner == 0:
        for zp in zones:
            dimZone = Internal.getZoneDim(zp)
            if dimZone[0] == 'Structured':
                _fillCornerGhostCellsStruct2__(zp,d)
    return None

#-------------------------------------------------------------------------------
# addGhostCells for a structured zone (first loop)
# IN: zp: zone to add ghost cells
# IN: bp: tree/basis/list of zones
# IN: d: nb of ghost cells in all the directions
# IN: modified: list of fields to be modified (GridCoordinates, FlowSolution or some fields only)
# IN: nodesRef: dictionary that identifies the first zone of a given name 
# (if several zones with the same name in a multi-base tree.
# IN/OUT: validjoins: list of nodes of type GridConnectivity1to1_t f
# IN/OUT: dictjoins: dictionary of [zone, oppositezonename]
# IN/OUT: indexjoin: increment that relates the dictjoins to validjoins 
#-------------------------------------------------------------------------------
def _addGhostCellsStruct__(zp, bp, d, modified, nodesRef, validjoins, dictjoins, indexjoin):
    # In first loop, a dictionnary named 'validjoins' will defined the joins
    # with ghost cells implemented, and their parent base.
    # This dictionnary will be used in second loop
    dimZone = Internal.getZoneDim(zp)

    # 1. init with 0th order extrapolation
    _initWithExtrapStruct__(zp, dimZone, modified, d)

    # 2. correct with bcdataset
    _setBCDataInGhostCellsStruct__(zp, dimZone, modified, d)
    #Internal._rmNodesByType(zp,'BCDataSet_t')

    # 3. fill ghost cells for BCMatch and BCNearMatch
    return _fillGhostCellsForStructBCMatch__(zp, d, dimZone, modified, nodesRef, validjoins, dictjoins, indexjoin)

#-------------------------------------------------------------------------------
# addGhostCells for a NGON zone
# IN: zp: zone to add ghost cells
# IN: bp: tree/basis/list of zones
# IN: stdNode: 1 if bp is a standard node (not a list of zones)
# IN: isghost: True if the fields are modified (to add a ZoneRind_t node )
# IN: d: nb of ghost cells in all the directions
# IN: modified: list of fields to be modified (GridCoordinates, FlowSolution or some fields only)
#-------------------------------------------------------------------------------
def _addGhostCellsNGON__(zp, bp, d, isghost, stdNode, modified):
    (parentN, dN) = Internal.getParentOfNode(bp, zp)
    fieldsn, fieldsc = getFieldsInContainer__(zp, modified)
    #
    # fieldsn and fieldsc contain coords, solution at nodes and centers
    #
    if fieldsc == [] and fieldsn == []: 
        print('Warning: addGhostCells: variables not found. No ghost cells added.')
    elif fieldsc == [] and fieldsn != []: # nodes only 
        out = Converter.converter.addGhostCellsNGonNodes(fieldsn,d)
        #PyTree.setFields([out], zp, 'nodes', writeDim=True)
    elif fieldsc != [] and fieldsn == []:# centers only
        isghost = False
        # out = Converter.converter.addGhostCellsNGonCenters(fieldsc,d)
        # PyTree.setFields([out], zp, 'centers', writeDim=False)
    else: # both
        out = Converter.converter.addGhostCellsNGonBoth(fieldsn,fieldsc,d)
        #PyTree.setFields([out[0]], zp, 'nodes', writeDim=True)
        #PyTree.setFields([out[1]], zp, 'centers', writeDim=False)
    # zp is now modified
    if parentN is not None: 
        if stdNode == 0: parentN[dN] = zp
        else: parentN[2][dN] = zp
    else: bp = zp # liste de zones   
    return isghost

#==============================================================================
# Remove ghost Cells on GridConnectivity and all fields of FlowSolution
# and FlowSolution#Centers
# IN: z: les zones doivent avoir des dimensions en noeuds avec ghost cells
# IN: modified: si [], modifie GridCoordinates+FlowSolution+FlowSolution#Centers
# adaptBCs=0: BC non modifiees par addGC
#         =1: ranges des BC phys + BCMatch + BCOverlap modifies par addGC
# Return a new zone
#==============================================================================
def rmGhostCells(t, b, d, adaptBCs=0, modified=[]):
    bp = Internal.copyRef(b)
    return _rmGhostCells(t, bp, d, adaptBCs, modified)

def _rmGhostCells(t, b, d, adaptBCs=0, modified=[]):
    if modified == []: # Default value
        modified = [Internal.__GridCoordinates__,Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__]

    # Multi-base: if t contains several zones with the same name, only the first zone found is kept
    nodesRef = getFirstReferencedZones__(t)
    stdNode = Internal.isStdNode(b)
    nodes = Internal.getZones(b)

    for zp in nodes:
        dim = Internal.getZoneDim(zp)
        if dim[0] == 'Unstructured':
            typeElt = dim[3]
            if typeElt == 'NGON':
                pass
                #_rmGhostCellsNGON__(zp, b, d, stdNode, modified)
            else: 
                #print('Warning: rmGhostCells: not yet implemented for BE zones.')
                pass # other elt types 

        else: # Structured:
            # zone dimension with ghost cells
            img = dim[1]; jmg = dim[2]; kmg = dim[3]
            # Traitement particulier 2D
            if kmg == 2:
                try:
                    import Transform.PyTree as T               
                    zpp = T.subzone(zp, (1,1,1), (-1,-1,1)); kmg = 1
                    zp[2] = zpp[2] # force in place
                    zp[1] = zpp[1]
                    dim = Internal.getZoneDim(zp)
                    #print('Warning: rmGhostCells: matching boundaries will be lost.')
                except: pass
            # zone dimension without ghost cells
            im = img-2*d; jm = jmg-2*d; km = kmg-2*d
            # check dimensions
            dimZone = dim[4]
            if dimZone == 3:
                if im <= 0 or jm <= 0 or km <= 0:
                    print("Warning: rmGhostCells: cannot remove %d ghost cells. Try less."%d)
                    return b
            elif dimZone == 2:
                if im <= 0 or jm <= 0:
                    print("Warning: rmGhostCells: cannot to remove %d ghost cells. Try less."%d)
                    return b
            else:
                if im <= 0:
                    print("Warning: rmGhostCells: cannot remove %d ghost cells. Try less."%d)
                    return b

            #-----------------------
            # general data
            #-----------------------
            _rmGhostCellsStruct__(zp, dim, modified, d)

            #-------------------------------
            # Adapts physical BC
            #-------------------------------
            if adaptBCs == 1: _adaptBCStruct__(zp, dim, -d, nodesRef)               

            #---------------------------------
            # Remove ghost dimensions for zone
            #---------------------------------
            _updateZoneDim__(zp, -d)

        #------------------------------------
        # remove Rind_t node for zone
        #------------------------------------ 
        Internal._rmNodesByName(zp, 'ZoneRind')
    return b

#===============================================================================
# IN: zp: zone for which fields are extracted
# IN: modified: list of names of the container 
#               __GridCoordinates__, __FlowSolution__ and list of variables
# IN: coords=True: return the coordinates array if modified=Internal.__FlowSolutionCenters__
# OUT: return the list of arrays corresponding to fields at nodes and at centers
#===============================================================================
def getFieldsInContainer__(zp, modified, coords=True):
    fieldsn = []; fieldsc = []
    for name in modified:
        if name == Internal.__GridCoordinates__ or name == Internal.__FlowSolutionNodes__: 
            if fieldsn == []: fieldsn = PyTree.getFields(name,zp)[0]
            else: 
                sol = PyTree.getFields(name,zp)[0]
                if sol != []: fieldsn = Converter.addVars([fieldsn,sol])
        elif name == Internal.__FlowSolutionCenters__: 
            fieldsc = PyTree.getFields(name,zp)[0]
            if coords:
                fieldsn = PyTree.getFields(Internal.__GridCoordinates__,zp)[0]
        else:
            if not isinstance(name, list): name = [name]
            for vname in name:
                if vname not in PyTree.getVarNames(zp)[0]: continue
                vars0 = vname.split(':',1)
                if len(vars0) == 1: # nodes
                    if fieldsn == []: fieldsn = PyTree.getField(vname,zp)[0]
                    else:
                        fn = PyTree.getField(vname,zp)[0]
                        fieldsn = Converter.addVars([fieldsn,fn])
                else: 
                    if vars0[0] == 'centers': 
                        if fieldsc == []: fieldsc = PyTree.getField(vname,zp)[0]
                        else: 
                            fc = PyTree.getField(vname,zp)[0]
                            fieldsc = Converter.addVars([fieldsc,fc])
                    elif vars0[0] == 'nodes': # nodes
                        if fieldsn == []: fieldsn = PyTree.getField(vname,zp)[0]
                        else: 
                            fn = PyTree.getField(vname,zp)[0]
                            fieldsn = Converter.addVars([fieldsn,fn])
                    else:
                        if fieldsn == []: fieldsn = PyTree.getField(vname,zp)[0]
                        else:
                            fn = PyTree.getField(vname,zp)[0]
                            fieldsn = Converter.addVars([fieldsn,fn])

    return fieldsn, fieldsc

def _rmGhostCellsNGON__(zp, bp, d, stdNode, modified):
    (parentN, dN) = Internal.getParentOfNode(bp, zp)
    fieldsn, fieldsc = getFieldsInContainer__(zp, modified)
    #
    # fieldsn and fieldsc contain coords, solution at nodes and centers
    #
    if fieldsc == [] and fieldsn == []:
        print('Warning: rmGhostCells: variables not found. No ghost cells removed.')
    elif fieldsc == [] and fieldsn != []:
        out = Converter.converter.rmGhostCellsNGonNodes(fieldsn, d)
        #PyTree.setFields([out], zp, 'nodes', writeDim=True)
    elif fieldsc != [] and fieldsn == []:
        pass
        # out = Converter.converter.rmGhostCellsNGonCenters(fieldsc,d)
        # PyTree.setFields([out], zp, 'centers', writeDim=False)
    else:
        out = Converter.converter.rmGhostCellsNGonBoth(fieldsn,fieldsc,d)
        #PyTree.setFields([out[0]], zp, 'nodes', writeDim=True)
        #PyTree.setFields([out[1]], zp, 'centers', writeDim=False)
    # zp is now modified
    if parentN is not None: 
        if stdNode == 0: parentN[dN] = zp
        else: parentN[2][dN] = zp
    else: bp = zp # liste de zones
    return None

#==============================================================================
# Traitement de raccord : construction des cellules fictives au niveau des 
# aretes entre deux fenetres topologiques et des coins du maillage GC
# au depart : cellules degenerees
# cas structure seulement
#==============================================================================
def _fillCornerGhostCellsStruct2__(zp, d):
    return Converter.converter.fillCornerGhostCells2(zp, d, Internal.__GridCoordinates__,
                                                     Internal.__FlowSolutionNodes__,
                                                     Internal.__FlowSolutionCenters__)

#==============================================================================
# Traitement de raccord, remplissage des cellules fictives en 3 etapes:
# 1- treatment=0: cellules fictives correspondant aux cellules reelles du
# bloc donneur (bloc de raccord oppose),
# 2- treatment=1: cellules fictives correspondant aux aretes de la grille
# fictive,
# 3- treatment=2: cellules fictives correspondant aux coins de la grille
# fictive.
#------------------------------------------------------------------------------
# IN: d: depth (number of ghost cells in a mesh direction)
# IN/OUT: zp: current zone
# IN/OUT: join: current join where treatment is applied
#==============================================================================
def fillJoinGhostCellsStruct__(zp, join, modified, joininfo,
                               d, treatment=0):
    zdonor = joininfo[0]; zdonor = Internal.copyRef(zdonor)
    prange = joininfo[1]
    prangedonor = joininfo[2]
    trirac = joininfo[3]
    # nmratio=[1,2,1] means that the receptor zone is twice as coarse as the donor zone in the j direction
    nmratio = joininfo[4]
    translVect = joininfo[5]
    rotationData = joininfo[6]
    typegc = 0
    if len(nmratio) == 3: typegc = 1

    # Periodicity by rotation: vectors that must be modified
    vectorsR = []
    for v in ['Velocity','Momentum']:
        vectname = [v+'X',v+'Y',v+'Z']
        vectorsR.append(vectname)

    # current zone dimensions
    dim = Internal.getZoneDim(zp)
    dimZone = dim[4]

    # donor zone dimensions
    dimdonor = Internal.getZoneDim(zdonor)
    if translVect != [] and rotationData == []: # translation seulement
        try: import Transform
        except: raise ImportError("addGhostCells: periodicity requires Transform module.")
        tvx = -translVect[0]; tvy = -translVect[1]; tvz = -translVect[2]
        if treatment > 0: zdonor = updateZoneDim__(zdonor,d) 
        isCoordPresent = 0
        if Internal.__GridCoordinates__ in modified: isCoordPresent = 3
        else:
            for coord in coordinate:
                if coord in modified: isCoordPresent += 1
        if isCoordPresent == 3:
            coords = PyTree.getFields(Internal.__GridCoordinates__,zdonor)
            coords = Transform.translate(coords, (tvx,tvy,tvz))
            PyTree.setFields(coords, zdonor, 'nodes', writeDim=False)

    elif translVect == [] and rotationData != []: # rotation seulement
        try: import Transform
        except: raise ImportError("addGhostCells: periodicity requires Transform module.")
        xc = rotationData[0]; yc = rotationData[1]; zc = rotationData[2]
        vx = rotationData[3]; vy = rotationData[4]; vz = rotationData[5]
        angle = rotationData[6]
        if vy > 0.: anglep = angle
        else: anglep = -angle
        if treatment > 0: zdonor = updateZoneDim__(zdonor,d) 
        fieldsn, fieldsc = getFieldsInContainer__(zdonor, modified, coords=False)
        if fieldsn != []:
            fieldsn = Transform.rotate(fieldsn, (xc,yc,zc), (vx,vy,vz), anglep, vectors=vectorsR) 
            PyTree.setFields([fieldsn], zdonor, 'nodes', writeDim=False)
        if fieldsc != []:
            fieldsc = Transform.rotate(fieldsc, (xc,yc,zc), (vx,vy,vz), anglep, vectors=vectorsR) 
            PyTree.setFields([fieldsc], zdonor, 'centers', writeDim=False)
    elif  translVect != [] and rotationData != []:
        print('Periodicity by rotation and translation cannot be applied at the same time.')

    PR = prange[0][1]; PRD = prangedonor[0][1]
    borderinfoN = getInfoForFillJoinsStruct__(PR, PRD, trirac, dim, dimdonor, d, 'Vertex', dimZone, nmratio=nmratio,
                                              corner=treatment)
    borderinfoC = getInfoForFillJoinsStruct__(PR, PRD, trirac, dim, dimdonor, d, 'CellCenter', dimZone, nmratio=nmratio,
                                              corner=treatment)
    for name in modified:
        containers, containersD, locations = getContainers__(name, zp, zdonor)
        noc = 0

        for cont in containers:
            loc = locations[noc]
            if loc == 'CellCenter': borderinfo = borderinfoC
            else: borderinfo = borderinfoN
            if treatment == 0:
                cont[1] = fillJoinBorderStruct__(PR, PRD, trirac, dim, dimdonor, cont, containersD[noc],
                                                 d, loc, dimZone, borderinfo, typegc)
            else:             
                cont[1] = fillJoinCornerStruct__(PR, PRD, trirac, dim, dimdonor, cont, containersD[noc],
                                                 d, loc, dimZone, borderinfo, treatment, typegc)
            noc += 1
    return None

#==============================================================================
# Method called by addGhostCells.
# Fill ghost cells (except in corners) for a field at a matching join
# Return new field - for structured zones only
#==============================================================================
def fillJoinBorderStruct__(prange, prangedonor, trirac, dim, dimdonor, frecv,
                           fdonor, d, loc, dimZone, borderinfo, typegc):
    # copy of recv field (this field has already ghostcells)
    a = frecv[1]; f = fdonor[1]
    [arrayborder,dim1,dim2,listdonor,direction,dirdonor,incrrecv,incrdonor,im,imdonor,shiftDir1,shiftDir2,isFine]=borderinfo    
    # fill ghost values (except corner values) for join
    if typegc == 0:
        Converter.converter.fillJoin(a, f, arrayborder, listdonor, loc, dim1, dim2, direction, dirdonor,
                                     incrrecv, incrdonor, dimZone, d, im, imdonor)
    else:
        if loc != 'CellCenter':
            Converter.converter.fillJoinNMNodes(a, f, arrayborder, listdonor, dim1, dim2, direction, dirdonor,
                                                incrrecv, incrdonor, dimZone, d, im, imdonor,
                                                shiftDir1, shiftDir2, isFine)
        else: 
            Converter.converter.fillJoinNMCenters(a, f, arrayborder, listdonor, 
                                                  dim1, dim2, direction, dirdonor,
                                                  incrrecv, incrdonor, dimZone, d, im, imdonor,
                                                  shiftDir1, shiftDir2, isFine)
    return a
#===============================================================================
#
#===============================================================================
def getInfoForFillJoinsStruct__(prange, prangedonor, trirac, dim, dimdonor, 
                                d, loc, dimZone, nmratio=[], corner=0):
    if corner == 0: dloc = 0
    else: dloc = d
    # recv zone dimension
    im = dim[1]; jm = dim[2]; km = dim[3]

    # donor zone dimension
    imdonor = dimdonor[1]
    if dimZone != 1: jmdonor = dimdonor[2]
    if dimZone == 3: kmdonor = dimdonor[3]
    if loc == 'CellCenter':
        imdonor = imdonor-1
        im = im-1 ; jm = jm-1 ; km = km-1        
        if dimZone != 1: jmdonor = jmdonor-1
        if dimZone == 3: kmdonor = kmdonor-1        

    # get directions of recv and donor borders
    direction = getDirBorderStruct__(prange, dimZone)
    dirdonor = getDirBorderStruct__(prangedonor, dimZone)
    # get list of border nodes and direction of border
    if dimZone != 1:
        [arrayborder, dim1border, dim2border] = getBorderIndicesStruct__(prange,dim,direction,d,loc,dimZone)
        listdonor = getJoinDonorIndicesStruct__(prange,prangedonor,dimdonor,dirdonor,trirac,dloc,loc,dimZone)

    incrrecv = 1 ; incrdonor = 1
    absdirection = abs(direction)
    absdirdonor = abs(dirdonor)
    shiftDir1 = 1; shiftDir2 = 1
    isFine = 0
    if nmratio != []:
        nmratioFact = nmratio[0]*nmratio[1]*nmratio[2]      
        if nmratioFact < 1: 
            isFine = 1
            iratio = int(math.ceil(1./nmratio[0]))
            jratio = int(math.ceil(1./nmratio[1]))
            kratio = int(math.ceil(1./nmratio[2]))
        else: 
            isFine = -1
            iratio = int(nmratio[0])
            jratio = int(nmratio[1])
            kratio = int(nmratio[2])
        if absdirection == 1:
            shiftDir1 = jratio 
            shiftDir2 = kratio    
        elif absdirection == 2:
            shiftDir1 = iratio 
            shiftDir2 = kratio
        else:
            shiftDir1 = iratio
            shiftDir2 = jratio

    if corner == 0: # for fillJoinBorder
        # liste des directions de la zone receveur
        dir_rcv = [0]*dimZone
        dir_rcv[absdirection-1] = 1
        # liste des directions de la zone donneuse
        dir_donor = [0]*dimZone
        dir_donor[absdirdonor-1] = 1

        if dimZone == 3:
            incrrecv = increment__(dir_rcv[0],dir_rcv[1],dir_rcv[2],im,jm,km,d)
            incrdonor = increment__(dir_donor[0],dir_donor[1],dir_donor[2],imdonor,jmdonor,kmdonor,0)
        elif dimZone == 2:
            incrrecv = increment__(dir_rcv[0],dir_rcv[1],im,jm,d)
            incrdonor = increment__(dir_donor[0],dir_donor[1],imdonor,jmdonor,0)

        return [arrayborder, dim1border, dim2border,listdonor, direction, dirdonor, incrrecv, incrdonor, im, imdonor, shiftDir1, shiftDir2, isFine]

    else: # for fillJoinCorner
        incrrecvI = 0; incrrecvJ = 0
        incrdonorI = 0; incrdonorJ = 0
        # compute increments:
        #    For current (receiver) join:
        #    incrrecv: index between two cells in the direction of the frontier
        #    incrrecvI: index between two cells following the local axe I of the frontier
        #    incrrecvJ: index between two cells following the local axe J of the frontier
        #
        #    For donor join:
        #    incrdonor: index between two cells in the direction of the frontier
        #    incrdonorI: index between two cells following the local axe I of the frontier
        #    incrdonorJ: index between two cells following the local axe J of the frontier

        if dimZone == 3:
            if absdirection == 1:
                incrrecv = increment__(1,0,0,im,jm,km,d)
                incrrecvI = increment__(0,1,0,im,jm,km,d)
                incrrecvJ = increment__(0,0,1,im,jm,km,d)
                DonorI = trirac[1]
                DonorJ = trirac[2]
            elif absdirection == 2:
                incrrecv = increment__(0,1,0,im,jm,km,d)
                incrrecvI = increment__(1,0,0,im,jm,km,d)
                incrrecvJ = increment__(0,0,1,im,jm,km,d)
                DonorI = trirac[0]
                DonorJ = trirac[2]
            elif absdirection == 3:
                incrrecv = increment__(0,0,1,im,jm,km,d)
                incrrecvI = increment__(1,0,0,im,jm,km,d)
                incrrecvJ = increment__(0,1,0,im,jm,km,d)
                DonorI = trirac[0]
                DonorJ = trirac[1]

            if absdirdonor == 1:
                incrdonor = increment__(1,0,0,imdonor,jmdonor,kmdonor,d)
            elif absdirdonor == 2:
                incrdonor = increment__(0,1,0,imdonor,jmdonor,kmdonor,d)
            elif absdirdonor == 3:
                incrdonor = increment__(0,0,1,imdonor,jmdonor,kmdonor,d)

            # computation of incrdonorI and incrdonorJ in the same local (I,J) frame as for recv join
            absI = abs(DonorI); signI = absI//DonorI
            absJ = abs(DonorJ); signJ = absJ//DonorJ
            incrdonorI = signI*(absI-2)*(absI-3)//2 -  signI*(absI-1)*(absI-3)*(imdonor+2*d) + signI*(absI-1)*(absI-2)//2*(imdonor+2*d)*(jmdonor+2*d)
            incrdonorJ = signJ*(absJ-2)*(absJ-3)//2 -  signJ*(absJ-1)*(absJ-3)*(imdonor+2*d) + signJ*(absJ-1)*(absJ-2)//2*(imdonor+2*d)*(jmdonor+2*d)
        elif dimZone == 2:
            if abs(direction) == 1:
                incrrecv = increment__(1,0,im,jm,d)
                incrrecvI = increment__(0,1,im,jm,d)
                DonorI = trirac[1]
            elif abs(direction) == 2:
                incrrecv = increment__(0,1,im,jm,d)
                incrrecvI = increment__(1,0,im,jm,d)
                DonorI = trirac[0]
            if abs(dirdonor) == 1:
                incrdonor = increment__(1,0,imdonor,jmdonor,d)
            elif abs(dirdonor) == 2:
                incrdonor = increment__(0,1,imdonor,jmdonor,d)
            # computation of incrdonorI and incrdonorJ in the same local (I,J) frame as for recv join
            absI = abs(DonorI); signI = absI//DonorI
            incrdonorI = signI*(2-absI) + signI*(absI-1)*(imdonor+2*d)
            incrdonorJ = 0

        return [arrayborder, dim1border, dim2border, listdonor, incrrecvI, incrrecvJ, incrdonorI, incrdonorJ, direction, dirdonor, incrrecv, incrdonor, shiftDir1, shiftDir2, isFine]

#==============================================================================
# Method called by addGhostCells.
# Fill corner ghost cells for a field at a matching join
# Return new field
# -----
# passage: 1 corresponds to the first passage for corner ghost cells
#          2 corresponds to the second (and last) passage for unfilled ghost cells 
#==============================================================================
def fillJoinCornerStruct__(prange, prangedonor, trirac, dim, dimdonor, frecv,
                           fdonor, d, loc, dimZone, borderinfo, passage, typegc):
    # copy of recv field (this field has already ghostcells)
    a = frecv[1]
    f = fdonor[1]

    [arrayborder, dim1, dim2, listdonor,incrrecvI,incrrecvJ,incrdonorI,incrdonorJ,direction,dirdonor,incrrecv,incrdonor, shiftDir1, shiftDir2, isFine] = borderinfo

    # fill ghost corner values for join
    if typegc == 0:
        Converter.converter.fillCornerGhostCells(a, f, arrayborder, listdonor, loc, dim1, dim2, incrrecvI, incrrecvJ,
                                                 incrdonorI, incrdonorJ, direction, dirdonor, incrrecv, incrdonor, 
                                                 dimZone, d, passage)
    else: # non implemented for nearmatch
        pass
    return a    

#==============================================================================
# Return a list with donor indices for a matching join between structured zones
#==============================================================================
def getJoinDonorIndicesStruct__(prange,prangedonor,dimdonor,dirdonor,trirac,d,loc,dimZone, shift=0):
    im = dimdonor[1]
    jm = dimdonor[2]
    # trirac
    t1 = trirac[0]; t2 = trirac[1]; t3 = 0
    if dimZone == 3: km = dimdonor[3]; t3 = trirac[2]
    if loc == 'CellCenter':
        im = im-1 ; jm = jm-1
        if dimZone == 3: km = km-1

    # Donor border window dimension
    [wimin,wimax,wjmin,wjmax,wkmin,wkmax] = Internal.range2Window(prangedonor)
    if loc == 'CellCenter':
        wimax = wimax-1
        wjmax = wjmax-1
        if dimZone == 3: wkmax = wkmax-1
    # direction of recv border
    direction = getDirBorderStruct__(prange,dimZone)
    dirdonor = getDirBorderStruct__(prangedonor,dimZone)
    # get list of border nodes and direction of border
    [arrayborder, dim1border, dim2border] = getBorderIndicesStruct__(prangedonor,dimdonor,dirdonor,d,loc,dimZone, shift)
    if dimZone == 3:
        dim1 = 3
        array = numpy.zeros((dim1,2), dtype=Internal.E_NpyInt, order='F')
        array[0,1] = wimax-wimin+1
        array[1,1] = wjmax-wjmin+1
        array[2,1] = wkmax-wkmin+1
    else:
        dim1 = 2
        array = numpy.zeros((dim1,2), dtype=Internal.E_NpyInt, order='F')
        array[0,1] = wimax-wimin+1
        array[1,1] = wjmax-wjmin+1
    # reorder indices of opposite (donor) matching join wrt the trirac
    #dimb1 = arrayborder.shape[0]
    dimb1 = dim1border
    return Converter.converter.getJoinDonorIndices(array,arrayborder,t1,t2,t3,direction,dimZone,dim1,dimb1)

#==============================================================================
# Return an array with adjacent nodes/cells indices of a border
#==============================================================================
def getBorderIndicesStruct__(prange, dim, direction, d, loc, dimZone, shift=0):
    # dimensions of zone
    im = dim[1] ; jm = dim[2] ; km = dim[3]
    if loc == 'CellCenter':
        im = im-1 ; jm = jm-1 ; km = km-1        

    # border window dimension
    [wimin,wimax,wjmin,wjmax,wkmin,wkmax] = Internal.range2Window(prange)
    if loc == 'CellCenter':
        wimax = wimax-1
        if dimZone != 1: wjmax = wjmax-1
        if dimZone == 3: wkmax = wkmax-1

    dim1 = wjmax-wjmin+1
    dim2 = wkmax-wkmin+1
    if abs(direction) == 2:
        dim1 = wimax-wimin+1
    elif abs(direction) == 3:
        dim1 = wimax-wimin+1
        dim2 = wjmax-wjmin+1

    # 3D Treatment
    if dimZone == 3:
        arrayIndices = numpy.empty((dim1*dim2), dtype=Internal.E_NpyInt, order='F')
        Converter.converter.getJoinBorderIndices(arrayIndices, dim1, im, jm, km,
                                                 wimin, wimax, wjmin, wjmax, wkmin, wkmax,
                                                 direction, dimZone, d, shift)
    # 2D Treatment
    else:
        arrayIndices = numpy.empty((dim1), dtype=Internal.E_NpyInt, order='F')
        Converter.converter.getJoinBorderIndices(arrayIndices, dim1, im, jm, km,
                                                 wimin, wimax, wjmin, wjmax, wkmin, wkmax,
                                                 direction, dimZone, d, shift)
    return [ arrayIndices, dim1, dim2 ]

#==============================================================================
# Get direction of a border:
# -1: Imin, 1: Imax
# -2: Jmin, 2: Jmax
# -3: Kmin, 3: Kmax
# Same routine for 2D or 3D
#==============================================================================
def getDirBorderStruct__(prange, dimZone):
    # border window dimension
    [wimin,wimax,wjmin,wjmax,wkmin,wkmax] = Internal.range2Window(prange)
    # direction of border
    if wimin == wimax:
        if wimin == 1: direction = -1    # direction I, border min
        else: direction = 1    # direction I, border max
    if dimZone != 1:
        if wjmin == wjmax:
            if wjmin == 1: direction = -2    # direction J, border min
            else: direction = 2    # direction J, border max
    if dimZone == 3:
        if wkmin == wkmax:
            if wkmin == 1: direction = -3    # direction K, border min
            else: direction = 3    # direction K, border max
    return direction

#==============================================================================
# return increment depending on direction
#==============================================================================
def increment__(*thetuple):
    if len(thetuple) == 7:
        (incri,incrj,incrk,im,jm,km,d) = thetuple
        incr = incri + incrj*(im+2*d) + incrk*(im+2*d)*(jm+2*d)
    elif len(thetuple) == 5:
        (incri,incrj,im,jm,d) = thetuple
        incr = incri + incrj*(im+2*d)
    return incr

#==============================================================================
# Get direction of bc/connectivity window with range 'prange'
# Validity only for window with PointRange. Not implemented with PointList
# Return direction as follows : 0: imin, 1: imax, 2:jmin, 3:jmax, 4:kmin; 5: kmax
# IN: dim: zone dimension
# IN: prange: range of bc/connectivity window
#==============================================================================
def getDirection__(dim, prange):
    [wimin,wimax,wjmin,wjmax,wkmin,wkmax] = Internal.range2Window(prange[0][1])
    di = wimax-wimin
    dj = wjmax-wjmin
    dk = wkmax-wkmin
    # set direction of boundary window: 0: imin, 1: imax, 2:jmin, 3:jmax, 4:kmin; 5: kmax
    if dim == 3:
        if di == 0 and wimax == 1: direction = 0
        elif di == 0 and wimax != 1: direction = 1
        elif dj == 0 and wjmax == 1: direction = 2
        elif dj == 0 and wjmax != 1: direction = 3
        elif dk == 0 and wkmax == 1: direction = 4
        else: direction = 5
    elif dim == 2:
        if di == 0 and wimax == 1: direction = 0
        elif di == 0 and wimax != 1: direction = 1
        elif dj == 0 and wjmax == 1: direction = 2
        else: direction = 3
    else:
        if wimax == 1: direction = 0
        else: direction = 1
    return direction

#==============================================================================
# Change PointRange taking into account a change of depth=d in grid
# Implemented only for structured grids
# IN: dimZone: dimension of the zone (1, 2 or 3)
# IN: direction: direction of the border
# IN: d: depth change of the grid
# IN: corners: 0: range not defined on corners, 1: range defined on corners
# IN/OUT: prange: PointRange node
#==============================================================================
def changePointRange__(prange, dimZone, direction, d, corners=0, extend=0):
    pr = numpy.copy(prange[0][1])
    ijk = int(direction/2)
    minmax = direction%2

    for ind in range(dimZone[4]):
        if ind != ijk:
            #pr[ind][0] = max(pr[ind][0]+d,1)
            #pr[ind][1] = pr[ind][1]+(1+corners)*d
            if ind == 0: N = dimZone[1]
            elif ind == 1: N = dimZone[2]
            else: N = dimZone[3]
            indm = pr[ind][0]+d
            # if (indm <= d+1): indm = 1 # corner
            if indm <= d+1: 
                if extend == 0: indm = d+1 #  nouvelle version : ne pas etendre les BC physiques
                else: indm = 1 # corner

            pr[ind][0] = indm
            indm = pr[ind][1]+d
            # if (indm >= N+d): indm = N+2*d # corner
            if indm >= N+d: 
                if extend == 0: indm = N+d #  nouvelle version : ne pas etendre les BC physiques
                else: indm = N+2*d # corner

            pr[ind][1] = indm

    if minmax == 1: # decalage BC max
        pr[ijk][0] = pr[ijk][0]+2*d
        pr[ijk][1] = pr[ijk][1]+2*d
    prange[0][1] = pr
    return None

#-----------------------------------------------------------------------------
# return the containers defined by name in z 
#-----------------------------------------------------------------------------
def getContainers__(name, z, zdonor=None):
    containers = []; locations = []; containersD = []
    #-----------------------
    # deals with coordinates
    #-----------------------
    if name == Internal.__GridCoordinates__:
        for name in coordinate:
            coord = Internal.getNodesFromName2(z, name)
            if coord != []: 
                if zdonor is None: containers.append(coord[0]); locations.append('Vertex')
                else: 
                    coorddonor = Internal.getNodesFromName2(zdonor,name)
                    containersD.append(coorddonor[0])
                    containers.append(coord[0]); locations.append('Vertex')
    #------------------------
    # deals with FlowSolution
    #------------------------
    elif (name  == Internal.__FlowSolutionNodes__) or (name == Internal.__FlowSolutionCenters__):
        loc = 'Vertex'
        f = Internal.getNodesFromName2(z, name)
        if f != []:
            l = Internal.getNodesFromType1(f[0], 'GridLocation_t')
            if l != []: loc = Internal.getValue(l[0])
            else:
                if name == Internal.__FlowSolutionNodes__: loc = 'Vertex'
                else: loc = 'CellCenter'
            fdonor = []
            if zdonor is not None: fdonor = Internal.getNodesFromName2(zdonor, name)

            for i in f[0][2]:
                if i[3] == 'DataArray_t': 
                    if zdonor is None: containers.append(i); locations.append(loc)
                    else:
                        if fdonor != []: # zdonor != None
                            idonor = Internal.getNodesFromName1(fdonor[0], i[0])
                            if idonor != []: 
                                containers.append(i); locations.append(loc); containersD.append(idonor[0])

    #-----------------------------------------
    # deals with variable or list of variables
    #-----------------------------------------
    else:
        if not isinstance(name, list): name = [name]
        for v in name:
            vars = v.split(':',1)
            loc = 'Vertex'; container = Internal.__FlowSolutionNodes__
            if len(vars) == 2: # center value
                if vars[0] == 'centers':
                    loc = 'CellCenter'; container = Internal.__FlowSolutionCenters__
                    vars = vars[1]
                elif vars[0] == 'nodes': vars = vars[1]
                else: vars = v
            else: vars = v # node value
            # search vars in nodes/centers container
            containerNode = Internal.getNodesFromName1(z, container)
            if containerNode != []:
                ploc = 'Vertex'
                gloc  = Internal.getNodesFromType1(containerNode[0], 'GridLocation_t')
                if gloc != []: ploc = Internal.getValue(gloc[0])
                if zdonor is not None:
                    containerNodeDonor = Internal.getNodesFromName1(zdonor,container)
                    plocdonor = 'Vertex'
                    glocdonor = Internal.getNodesFromType1(containerNodeDonor[0], 'GridLocation_t')
                    if glocdonor != []: plocdonor = Internal.getValue(glocdonor[0])

                if ploc == loc:
                    flist = Internal.getNodesFromName1(containerNode[0], vars)
                    for f in flist:
                        if f != []:
                            if zdonor is None:
                                if f[3] == 'DataArray_t': containers.append(f); locations.append(loc)
                            else:
                                if plocdonor == loc:
                                    fdonor = Internal.getNodesFromName1(containerNodeDonor[0], f[0])
                                    if fdonor != []:
                                        containers.append(f); locations.append(loc); containersD.append(fdonor[0])

    if zdonor is None: return containers,locations
    else: return containers, containersD, locations

#-----------------------------------------------------------------------------
# init all ghost cells for zone z with extrapolation of zeroth order.
#-----------------------------------------------------------------------------
def _initWithExtrapStruct__(z, dim, modified, d):
    dimZone = dim[4]
    imz = dim[1]; jmz = dim[2]; kmz = dim[3]

    for name in modified:
        containers, locations = getContainers__(name, z)
        noc = 0
        for cont in containers:
            loc = locations[noc]
            b = cont[1]
            if loc == 'CellCenter': im = imz-1 ; jm = jmz-1 ; km = kmz-1
            else: im = imz; jm = jmz; km = kmz
            img = im+2*d ; jmg = jm+2*d ; kmg = km+2*d
            # copy real values and fill ghost values
            if dimZone == 3:
                a = numpy.empty((img,jmg,kmg), b.dtype, order='F')
                Converter.converter.cpyReal2Ghost(a, b, d, im, jm, km)
            elif dimZone == 2:
                a = numpy.empty((img,jmg), b.dtype, order='F')
                Converter.converter.cpyReal2Ghost(a, b, d, im, jm, 0)
            else:
                a = numpy.empty((img), b.dtype, order='F')
                Converter.converter.cpyReal2Ghost(a, b, d, im, 0, 0)
            cont[1] = a
            noc += 1
    return None

# Warning: z dimension is not already modified - containers of flow solutions are not consistent yet with zone dim
def _setBCDataInGhostCellsStruct__(z, dim, modified, d):
    if dim[0] == 'Unstructured': 
        print('addGhostCells: _setBCDataInGhostCellsStruct__: not valid for unstructured zones.')
        return None
    dimZone = dim[4]
    imz = dim[1] ; jmz = dim[2] ; kmz = dim[3]
    for name in modified:
        containers = Internal.getBCDataSetContainers(name, z) 
        if containers is not None:
            for cont in containers:
                varname = cont[0]
                loc = cont[1]
                locI = 0
                if loc == 'CellCenter': locI = 1
                dataSetInfo = cont[2]
                bcranges = []; dataBC = []
                for e in dataSetInfo:
                    if e[1] is not None: # firewall pour les mauvaises BCDataSet
                        bcranges.append(e[0])
                        dataBC.append(e[1])
                Converter.converter._setBCDataInGhostCellsStruct(z, bcranges, dataBC, imz, jmz, kmz, d, locI, varname, 
                                                                 Internal.__GridCoordinates__, Internal.__FlowSolutionNodes__,Internal.__FlowSolutionCenters__)
    return None


#-----------------------------------------------------------------------------
# Remove ghost cells in numpy arrays for zone z
#-----------------------------------------------------------------------------
def _rmGhostCellsStruct__(z, dim, modified, d):
    dimZone = dim[4]
    imz = dim[1] ; jmz = dim[2] ; kmz = dim[3]

    for name in modified:
        containers, locations = getContainers__(name,z)
        noc = 0
        for cont in containers:
            loc = locations[noc]
            b = cont[1]

            if loc == 'CellCenter': img = imz-1 ; jmg = jmz-1 ; kmg = kmz-1
            else: img = imz; jmg = jmz; kmg = kmz
            im = img-2*d ; jm = jmg-2*d ; km = kmg-2*d

            if dimZone == 3:
                a = numpy.empty((im,jm,km), b.dtype, order='F')
                Converter.converter.cpyGhost2Real(a,b,d,im,jm,km)
            elif dimZone == 2:
                a = numpy.empty((im,jm), b.dtype, order='F')
                Converter.converter.cpyGhost2Real(a,b,d,im,jm,0)
            else:
                a = numpy.empty((im), b.dtype, order='F')
                Converter.converter.cpyGhost2Real(a,b,d,im,0,0)
            cont[1] = a
            noc += 1

    return None

#==============================================================================
# IN: t: arbre contenant des containers "rind" (avec ghost cells)
# IN: d: nbre de ghost cells dans le rind
# IN: modified: la liste du nom des containers a modifier
# Modifie un ou plusieurs containers avec Rind pour le transformer en 
# container standard (centres ou noeuds)
#==============================================================================
def rmRindCells(t, d=-1, modified=[]):
    tp = Internal.copyRef(t)
    _rmRindCells(tp, d, modified)
    return tp

def _rmRindCells0(t):
    try:
        import Transform.PyTree as T
    except:
        print("Warning: rmRindCells requires Transform module.")
        return None

    for zp in Internal.getZones(t):
        dim = Internal.getZoneDim(zp)
        dimZone = dim[4]
        if dim[0] == 'Structured':
            imz = dim[1] ; jmz = dim[2] ; kmz = dim[3]
            rindnode = Internal.getNodeFromType(zp,'Rind_t')
            if rindnode is not None:
                rindnode = Internal.getValue(rindnode)
                rindimin = rindnode[0]; rindimax = rindnode[1]
                rindjmin = rindnode[2]; rindjmax = rindnode[3]
                rindkmin = rindnode[4]; rindkmax = rindnode[5]
                zpp = T.subzone(zp, (1+rindimin,1+rindjmin, 1+rindkmin),
                                (imz-rindimax,jmz-rindjmax,kmz-rindkmax))
                zp[2]=zpp[2]; zp[1]=zpp[1] # force in place    
    return None

def _rmRindCells(t, d=-1, modified=[]):
    if d == -1: return _rmRindCells0(t)

    stdNode = Internal.isStdNode(t)
    nodes = Internal.getZones(t)
    for zp in nodes:
        dim = Internal.getZoneDim(zp)
        dimZone = dim[4]
        if dim[0] == 'Structured':
            imz = dim[1] ; jmz = dim[2] ; kmz = dim[3]
            for name in modified:  # recherche les containers
                containers, locations = getContainers__(name, zp)
                noc = 0
                for cont in containers:
                    loc = locations[noc]
                    b = cont[1]
                    if loc == 'CellCenter': # output sera un champ en centres
                        im = imz-1; jm = jmz-1; km = kmz-1
                    else: # output sera un champ en noeuds
                        im = imz; jm = jmz; km = kmz

                    if dimZone == 3:
                        a = numpy.empty((im,jm,km), b.dtype, order='F')
                        Converter.converter.cpyGhost2Real(a,b,d,im,jm,km)
                    elif dimZone == 2:
                        a = numpy.empty((im,jm), b.dtype, order='F')
                        Converter.converter.cpyGhost2Real(a,b,d,im,jm,0)
                    else:
                        a = numpy.empty((im), b.dtype, order='F')
                        Converter.converter.cpyGhost2Real(a,b,d,im,0,0)
                    cont[1] = a
                    noc += 1               
    return None

#-----------------------------------------------------------------------------
# fill all fields for ghost cell corners for structured zones 
# only for BCMatch currently
# loop = 2: ghost cells corners
# loop = 3: only for 3D : deals with the very last ghost cell corners 
#-----------------------------------------------------------------------------
def _fillGhostCellCorners__(bp, modified, d, validjoins, dictjoins, loop=2):
    treatment = 0
    if loop == 2: treatment = 1
    elif loop == 3: treatment = 2

    indexjoin = 0
    for j in validjoins:
        zp = dictjoins[indexjoin][0]          # zrcv: current zone
        dim = Internal.getZoneDim(zp)
        dimZone = dim[4]
        joininfo = dictjoins[indexjoin][1]
        if dim[0] == 'Structured':
            if (loop == 2) or (loop==3 and dimZone==3): 
                zdonorname = joininfo[0][0]
                zdonor = Internal.getNodesFromName2(bp, zdonorname)        
                # Check if donor exists.
                # If not, zdonor does not have ghost cells, thus this second pass is useless.
                zd = []
                for itz in zdonor:
                    if itz[3] == 'Zone_t':
                        zd = itz
                        break
                if zd == []: continue
                zdonor = zd
                indexjoin += 1
                joininfo[0] = zdonor
                fillJoinGhostCellsStruct__(zp, j, modified, joininfo, d, treatment=treatment)
        else: pass
    return None

#-----------------------------------------------------------------------------
# Fill the ghost cells using the BCMatch connectivity between structured zones
#-----------------------------------------------------------------------------
def _fillGhostCellsForStructBCMatch__(zp, d, dimZone, modified, nodesRef,
                                      validjoins, dictjoins, indexjoin):
    dimZone = dimZone[4]
    # BCMatch
    joinlist = Internal.getNodesFromType2(zp, 'GridConnectivity1to1_t')
    joindirs = []
    for i in joinlist:
        joininfo = getJoinInfo__(i, 'BCMatch', dimZone, nodesRef, d)
        if joininfo != []:
            # store valid joins (and its parent base) for second loop treatment
            validjoins.append(i)
            #zdonor = joininfo[0]
            fillJoinGhostCellsStruct__(zp, i, modified, joininfo, d, treatment=0)
            dictjoins[indexjoin] = (zp,joininfo)
            indexjoin += 1

    if dimZone == 1: return validjoins, dictjoins, indexjoin 

    # BCNearMatch
    joinlist = Internal.getNodesFromType2(zp, 'GridConnectivity_t')
    joindirs = []
    for i in joinlist:
        typegc = Internal.getNodeFromName1(i, 'GridConnectivityType')
        if typegc is not None:
            val = Internal.getValue(typegc)
            if val == 'Abutting':
                joininfo = getJoinInfo__(i, 'BCNearMatch', dimZone, nodesRef, d)
                if joininfo != []:
                    # store valid joins (and its parent base) for second loop treatment
                    validjoins.append(i)
                    #zdonor = joininfo[0]
                    fillJoinGhostCellsStruct__(zp, i, modified, joininfo, d, treatment=0)
                    dictjoins[indexjoin] = (zp,joininfo)
                    indexjoin += 1
    return validjoins, dictjoins, indexjoin 

#-----------------------------------------------------------------------------
# adaptBCStruct: modifies ranges for structured zones:
# BCMatch, BCOverlap and Physical BCs defined in BC_t
#-----------------------------------------------------------------------------
def _adaptBCsForZone(zp, dimZone, d, nodesRef):
    if dimZone[0] == 'Structured': return _adaptBCStruct__(zp, dimZone, d, nodesRef)
    else: pass
    return None

def _adaptBCStruct__(zp, dimZone, d, nodesRef):
    dimPb = dimZone[4]

    # BCMatch
    joinlist = Internal.getNodesFromType2(zp, 'GridConnectivity1to1_t')
    for i in joinlist:
        prange = Internal.getNodesFromName1(i, 'PointRange')
        if prange != []:
            direction = getDirection__(dimZone[4], prange)
            changePointRange__(prange, dimZone, direction, d, extend=0)
            prangedonor = Internal.getNodesFromName1(i, 'PointRangeDonor')
            if prangedonor != []:
                zdonorname = Internal.getValue(i)
                if zdonorname in nodesRef: zdonor = nodesRef[zdonorname]
                else: zdonor = []
                # multi-bases approach
                # if several zdonor found corresponding to name v, zdonor=first zone found
                if zdonor == []: continue
                zdonor = zdonor[0]
                dimZoneDonor = Internal.getZoneDim(zdonor)
                #print('h2',zdonor[0], dimZoneDonor)
                direction = getDirection__(dimZoneDonor[4], prangedonor)
                changePointRange__(prangedonor, dimZoneDonor, direction, d, extend=0)
    # BCNearMatch
    joinlist = Internal.getNodesFromType2(zp, 'GridConnectivity_t')
    for i in joinlist:
        r = Internal.getNodeFromType2(i, 'GridConnectivityType_t')
        if r != []:
            val = Internal.getValue(r)
            if val== 'Abutting':     
                prange = Internal.getNodesFromName1(i, 'PointRange')
                if prange != []:
                    udd =  Internal.getNodesFromType1(i, 'UserDefinedData_t')
                    if udd == []: continue
                    direction = getDirection__(dimZone[4], prange)
                    changePointRange__(prange, dimZone, direction, d, extend=0)
                    prangedonor = Internal.getNodesFromName1(udd, 'PointRangeDonor')
                    if prangedonor != []:
                        zdonorname = Internal.getValue(i)
                        if zdonorname in nodesRef: zdonor = nodesRef[zdonorname]
                        else: zdonor = []
                        # multi-bases approach
                        # if several zdonor found corresponding to name v, zdonor=first zone found
                        if zdonor == []: continue
                        zdonor = zdonor[0]
                        dimZoneDonor = Internal.getZoneDim(zdonor)
                        #print('h3',zdonor[0], dimZoneDonor)
                        direction = getDirection__(dimZoneDonor[4], prangedonor)
                        changePointRange__(prangedonor, dimZoneDonor, direction, d, extend=0)
    # 'Overset'
    gcnnodes = Internal.getNodesFromType2(zp, 'GridConnectivity_t')
    for i in gcnnodes:
        r = Internal.getNodesFromType(i, 'GridConnectivityType_t')
        if r != []:
            val = Internal.getValue(r[0])
            if val == 'Overset':                           
                prange = Internal.getNodesFromName1(i, 'PointRange')
                if prange != []:
                    direction = getDirection__(dimPb, prange)
                    # change PointRange for extended mesh
                    changePointRange__(prange, dimZone, direction, d, extend=1)

    # Physical BCs
    bclist = Internal.getNodesFromType2(zp, 'BC_t')
    for bc in bclist:
        prange = Internal.getNodesFromName1(bc, 'PointRange')
        if prange != []:
            direction = getDirection__(dimPb, prange)
            # change PointRange for extended mesh
            changePointRange__(prange, dimZone, direction, d, extend=0)
    return None

#-----------------------------------------------------------------------------
# For multi-base pyTrees, if t contains several zones with the same name
# the first zone found is kept
# Returns the dictionary corresponding to the first zones found 
#-----------------------------------------------------------------------------
def getFirstReferencedZones__(t):
    nodesRef = {}
    bases = Internal.getBases(t)
    for ba in bases:
        nodes = Internal.getNodesFromType1(ba, 'Zone_t')
        for n in nodes:
            if n[0] not in nodesRef: nodesRef[n[0]] = [n]
    if len(bases) == 0: # t is a zone or a list of zones
        nodes = Internal.getNodesFromType1(t, 'Zone_t')
        for n in nodes:
            if n[0] not in nodesRef: nodesRef[n[0]] = [n] 
    return nodesRef

#---------------------------------
# modifies the dimensions of zone
# if d > 0: add ghost cells
# if d < 0: rm ghost cells
#---------------------------------
def _updateZoneDim__(zp, d):
    dimZone = Internal.getZoneDim(zp)
    if dimZone[0] != 'Unstructured':
        zp[1] = numpy.copy(zp[1])
        for j in range(zp[1].shape[1]-1):
            for i in range(zp[1].shape[0]):
                zp[1][i,j] = zp[1][i,j] + 2*d

def updateZoneDim__(z,d):
    zp = Internal.copyRef(z)
    _updateZoneDim__(zp,d)
    return zp

#=============================================================================
# Returns join info for BCMatch and BCNearMatch joins as a list:
# [zdonor,prange,prangedonor,trirac,nmratio,translvect,rotationData] 
# nmratio=[] for BCMatch and [2,1,1] for BCNearMatch
# translvect!=[] if translation-periodic
# rotationData = [] if no periodicity by rotation
# rotationData = [xc,yc,zc, axisX, axisY, axisZ, angle] if rotation-per 
#=============================================================================
def getJoinInfo__(join, jointype, dimZone, nodesRef, d):
    prange = Internal.getNodesFromName1(join, 'PointRange')
    nmr = []; translVect = []; rotationData = []
    if jointype == 'BCMatch':
        prangedonor = Internal.getNodesFromName1(join, 'PointRangeDonor')
        if prangedonor == []: return []
        transfo = Internal.getNodesFromName1(join, 'Transform')
    elif jointype == 'BCNearMatch':
        udd = Internal.getNodesFromType1(join,'UserDefinedData_t')
        if udd == []: return []
        prangedonor = Internal.getNodesFromName1(udd, 'PointRangeDonor')
        if prangedonor == []: return []
        transfo = Internal.getNodesFromName1(udd, 'Transform')
        nmrnode = Internal.getNodesFromName1(udd, 'NMRatio')
        if nmrnode == []: return []
        for nm in nmrnode[0][1]: nmr.append(nm)
    v = Internal.getValue(join)
    if v in nodesRef: zdonor = nodesRef[v]
    else: zdonor = []

    # if several zdonor found corresponding to name v, zdonor=first zone found
    if zdonor == []: return []
    zdonor = zdonor[0]
    dimdonor = Internal.getZoneDim(zdonor)
    #print('h1',zdonor[0], dimdonor)

    # get periodic info (if any)
    rotationData, translVect = Internal.getPeriodicInfo__(join)

    if dimdonor[0] != 'Unstructured':
        # get trirac
        # -----------------------
        if dimZone == 3: trirac = [1,2,3]
        elif dimZone == 2: trirac = [1,2]                
        else: trirac = [1]
        if transfo != []:
            trirac[0] = transfo[0][1][0]
            if dimZone != 1: trirac[1] = transfo[0][1][1]
            if dimZone == 3: trirac[2] = transfo[0][1][2]
        # get direction of donor border
        dirdonor = getDirBorderStruct__(prangedonor[0][1], dimZone)
        # check donor zone dimension
        # check if dimensions of donor border match with dimension of current zone
        check = True
        [wimin,wimax,wjmin,wjmax,wkmin,wkmax] = Internal.range2Window(prange[0][1])        
        [wimindonor,wimaxdonor,wjmindonor,wjmaxdonor,wkmindonor,wkmaxdonor]=Internal.range2Window(prangedonor[0][1])
        if dimZone == 3:
            delta = [wimax-wimin,wjmax-wjmin,wkmax-wkmin]
            deltadonor = [wimaxdonor-wimindonor,wjmaxdonor-wjmindonor,wkmaxdonor-wkmindonor]
            if jointype == 'BCMatch':
                if ((delta[1] != deltadonor[abs(trirac[1])-1]) or (delta[2] != deltadonor[abs(trirac[2])-1])): 
                    check = False
        elif dimZone == 2:
            delta = [wimax-wimin,wjmax-wjmin]
            deltadonor = [wimaxdonor-wimindonor,wjmaxdonor-wjmindonor]
            if jointype == 'BCMatch':
                if delta[1] != deltadonor[abs(trirac[1])-1]: check = False
        else:
            delta = [wimax-wimin]
            deltadonor = [wimaxdonor-wimindonor]

        if jointype == 'BCMatch':
            if delta[0] != deltadonor[abs(trirac[0])-1]: check = False
        if check:
            # check if donor zone is deep enough for "d" ghost cells
            imdonor = dimdonor[1]
            if dimZone != 1: jmdonor = dimdonor[2]
            else: jmdonor = 0 ### useless value for 1d: test is also on dirborder.
            if dimZone == 3: kmdonor = dimdonor[3]
            else: kmdonor = 0 ### useless value for 2d: test is also on dirborder.
            if (((abs(dirdonor) == 1) and (imdonor < d+1)) or 
                ((abs(dirdonor) == 2) and (abs(jmdonor) < d+1)) or 
                ((abs(dirdonor) == 3) and (kmdonor < d+1))):
                print("Warning: addGhostCells: dimension of matching zone is too small for", d, "ghost cells. Skip join treatment for ghost extension.")
                return []
            else:
                return [zdonor, prange, prangedonor, trirac, nmr, translVect, rotationData]

    return []

#==============================================================================
# Rend une subzone de a a partir d'une facelist
#==============================================================================
##def getLayer(a, faceList=None, nlayers=2):
#    import Transform.PyTree as T
#    Internal._adaptNFace2PE(a, remove=False) 
#    PE = Internal.getNodeFromPath(a, 'NGonElements/ParentElements')[1]    
#    n = faceList.size
#    ptFace = faceList.flat
#    elts = set()
#    faces = set()
#    for i in range(n):
#        ind = ptFace[i]-1
#        e1 = PE[ind,0]; e2 = PE[ind,1]
#        if e1 != 0: elts.add(e1-1)
#        if e2 != 0: elts.add(e2-1)
#    b = T.subzone(a, list(elts), type='elements')
#    return b

#==============================================================================
#def addGhostCellsP(t):
#    import Transform.PyTree as T
#    import Generator.PyTree as G
#    tp = Internal.copyRef(t)
#    tpp, ntype = Internal.node2PyTree(tp)
#    bases = Internal.getBases(tpp)
#    for b in bases:
#        c = 0
#        for z in b[2]:
#            if z[3] == 'Zone_t':
#                gc = Internal.getNodesFromType2(z, 'GridConnectivity_t')
#                zp = Internal.copyRef(z)
#                for g in gc:
#                    faceList = Internal.getNodeFromName1(g, 'PointListDonor')[1]
#                    donor = Internal.getValue(g)
#                    zdonor = Internal.getNodeFromName2(t, donor)
#                    layer = getLayer(zdonor, faceList)
#                    zp = T.join(zp, layer)
#                    zp = G.close(zp)
#                b[2][c] = zp
#            c += 1
#    zones = Internal.getZones(tpp)
#    #Internal._adapt2FastP(tp)
#    tp = Internal.pyTree2Node(tpp, ntype)
#    return tp
#
#==============================================================================
# Rend une subzone de a a partir d'une facelist sur 1 ou 2 voisinage
#==============================================================================
def getLayer(zD, zR, elts_old, mask, xyz0, no_layer, faceListD=None, faceListR=None):
    import Transform.PyTree as T

    zbc_R = Internal.getNodeFromType1(zR  , 'ZoneBC_t')

    NG_PE_D = Internal.getNodeFromPath(zD, 'NGonElements/ParentElements')[1]    
    NG_EC_D = Internal.getNodeFromPath(zD, 'NGonElements/ElementConnectivity')[1]    
    NG_IDX_D= Internal.getNodeFromPath(zD, 'NGonElements/FaceIndex')[1]    
    FA_IDX_D= Internal.getNodeFromPath(zD, 'NFaceElements/ElementIndex')[1]   
    FA_EC_D = Internal.getNodeFromPath(zD, 'NFaceElements/ElementConnectivity')[1]    
    coordxD = Internal.getNodeFromPath(zD, 'GridCoordinates/CoordinateX')[1]
    coordyD = Internal.getNodeFromPath(zD, 'GridCoordinates/CoordinateY')[1]
    coordzD = Internal.getNodeFromPath(zD, 'GridCoordinates/CoordinateZ')[1]

    if no_layer > 1:
        NG_PE_R = Internal.getNodeFromPath(zR, 'NGonElements/ParentElements')[1]    
        NG_EC_R = Internal.getNodeFromPath(zR, 'NGonElements/ElementConnectivity')[1]    
        NG_IDX_R= Internal.getNodeFromPath(zR, 'NGonElements/FaceIndex')[1]    
        FA_EC_R = Internal.getNodeFromPath(zR, 'NFaceElements/ElementConnectivity')[1]    
        FA_IDX_R= Internal.getNodeFromPath(zR, 'NFaceElements/ElementIndex')[1]    
        coordxR = Internal.getNodeFromPath(zR, 'GridCoordinates/CoordinateX')[1]
        coordyR = Internal.getNodeFromPath(zR, 'GridCoordinates/CoordinateY')[1]
        coordzR = Internal.getNodeFromPath(zR, 'GridCoordinates/CoordinateZ')[1]

        zbc_D = Internal.getNodeFromType1(zD   , 'ZoneBC_t')
        bc_D  = Internal.getNodesFromType1(zbc_D, 'BC_t')
        bc_R  = Internal.getNodesFromType1(zbc_R, 'BC_t')

    n = faceListD.size
    ptFaceD = faceListD.flat
    ptFaceR = faceListR.flat

    one = Internal.E_NpyInt(1)

    elts = set()
    for i in range(n):
        faceD = ptFaceD[i]-1
        if no_layer == 1:
            if NG_PE_D[ faceD, 0] == 0: e1D = NG_PE_D[faceD,1]
            else: e1D = NG_PE_D[faceD,0]
            elts.add(e1D-one)

            # calcul centre maille
            x0 = 0.; y0 = 0.; z0 = 0.
            eltD = FA_IDX_D[e1D-1]
            for face in range( 1, FA_EC_D[eltD ]+1):
                f       = FA_EC_D[ eltD  + face] -1
                ind_f   = NG_IDX_D[ f] 
                NvtxD   = NG_EC_D[ ind_f]
                for vtx in range( NvtxD ):
                    vtxD = NG_EC_D[ ind_f +vtx +1 ] -1
                    x0 = x0 + coordxD[vtxD]
                    y0 = y0 + coordyD[vtxD]
                    z0 = z0 + coordzD[vtxD]

            code = x0 + 1000000.*y0 +   1.e12*z0
            xyz0.add( code )

            #print('elts', elts)
        else:
            if (NG_PE_D[ faceD, 0] > NG_PE_D[ faceD, 1]):  e1D = NG_PE_D[faceD,1]
            else:                                          e1D = NG_PE_D[faceD,0] 
            #recuperation adresse element dans ElementConnectivity pour recuperer faces de l'element
            ind_e1D = FA_IDX_D[e1D-1]

            #print 'e1',e1D, faceD+1, mask
            #print 'GD',NG_PE_D[faceD,1],NG_PE_D[faceD,0]

            for fD in range( 1, FA_EC_D[ ind_e1D]+1):
                no_faceD = FA_EC_D[ ind_e1D + fD] -1
                if no_faceD != faceD:
                    e3D = NG_PE_D[no_faceD,0]
                    e4D = NG_PE_D[no_faceD,1]
                    #on test si la cellule est un raccord ou un BC
                    if e3D*e4D != 0:       
                        #print 'e34',e3D, e4D, NG_PE_D[faceD,0], NG_PE_D[faceD,1]
                        #print 'face',no_faceD+1, faceD+1, face
                        if   (e3D != e1D and ( e3D < mask[0] or e3D > mask[1]) ): 
                            #calcul centre maille
                            x0 = 0.
                            y0 = 0.
                            z0 = 0.
                            eltD =  FA_IDX_D[e3D-1]
                            for face in range(1, FA_EC_D[eltD ]+1):
                                f       = FA_EC_D[ eltD  + face] -1
                                ind_f   = NG_IDX_D[ f] 
                                NvtxD   = NG_EC_D[ ind_f]
                                for vtx in range(NvtxD):
                                    vtxD = NG_EC_D[ ind_f +vtx +1] -1
                                    x0=x0 + coordxD[vtxD]
                                    y0=y0 + coordyD[vtxD]
                                    z0=z0 + coordzD[vtxD]

                            size1_xyz = len( xyz0)
                            code = x0 + 1000000.*y0 + 1.e12*z0
                            xyz0.add(code)
                            size2_xyz = len( xyz0)
                            if size2_xyz > size1_xyz: elts.add(e3D-one)

                        elif e4D < mask[0] or e4D > mask[1]: 
                            #calcul centre maille
                            x0 = 0.; y0 = 0.; z0 = 0.
                            eltD =  FA_IDX_D[e4D-1]
                            for face in range(1, FA_EC_D[eltD ]+1):
                                f       = FA_EC_D[ eltD  + face] -1
                                ind_f   = NG_IDX_D[ f]
                                NvtxD   = NG_EC_D[ ind_f]
                                for vtx in range(NvtxD):
                                    vtxD = NG_EC_D[ ind_f +vtx +1 ] -1
                                    x0 = x0 + coordxD[vtxD]
                                    y0 = y0 + coordyD[vtxD]
                                    z0 = z0 + coordzD[vtxD]

                            size1_xyz = len( xyz0)
                            code = x0 + 1000000.*y0 + 1.e12*z0
                            xyz0.add(code)
                            size2_xyz = len(xyz0)
                            if size2_xyz > size1_xyz: elts.add(e4D-one)
                    else:
                        #on cherche le matching entre la face BC dans donneur et receveur
                        faceR= ptFaceR[i]-1
                        e1R  = NG_PE_R[faceR,0]
                        e2R  = NG_PE_R[faceR,1]
                        #element couche externe a forcemment un no + grand
                        if e1R > e2R: elR = e1R
                        else          : elR = e2R  
                        eltR =  FA_IDX_R[elR-1]
                        face_bingo = []
                        faceSearch = True
                        faceR      = 1
                        while faceSearch:
                        #for faceR in range( 1, FA_EC_R[eltR]+1):
                            no_faceR = FA_EC_R[ eltR + faceR] -1
                            e3R      = NG_PE_R[no_faceR,0]
                            e4R      = NG_PE_R[no_faceR,1]
                            ind_fD   = NG_IDX_D[no_faceD]
                            ind_fR   = NG_IDX_R[no_faceR]
                            NvtxD    = NG_EC_D[ ind_fD]
                            NvtxR    = NG_EC_R[ ind_fR]
                            #on elimine les faces non externes, celle du raccord ou celle possedant un nbr de vextex different
                            if no_faceR != faceR and e3R*e4R == 0 and NvtxD== NvtxR:

                                searchD    = True
                                bingo      = 0
                                cd =1
                                while searchD:
                                    vtxD = NG_EC_D[ ind_fD+cd]-1
                                    searchR= True
                                    cr =1
                                    while searchR:
                                        vtxR = NG_EC_R[ ind_fR+cr]-1
                                        #if(elR==26 and no_faceR ==104):
                                        #  print 'xD=',  coordxD[vtxD], coordyD[vtxD],coordzD[vtxD]
                                        #  print 'xr=',  coordxR[vtxR], coordyR[vtxR],coordzR[vtxR]
                                        #  print 'cr=', cr
                                        if coordxD[vtxD] == coordxR[vtxR]:
                                            if coordyD[vtxD] == coordyR[vtxR]:
                                                if coordzD[vtxD] == coordzR[vtxR]:  
                                                    bingo  +=1 
                                                    searchR = False
                                        cr +=1
                                        if cr == (NvtxR+1): searchR = False

                                    cd +=1
                                    if cd == (NvtxD+1): searchD = False
                                    if bingo == NvtxR: 
                                        face_bingo.append(no_faceR)
                                        faceSearch = False

                            faceR += 1
                            if (faceR == (FA_EC_R[eltR]+1)) : faceSearch = False

                        #on modifie le Nombre de face BC si la face n 'est pas un raccord
                        for fbingo in  face_bingo:     
                            for gD in bc_D:
                                ptlist = Internal.getNodeFromName1(gD, 'PointList')[1]
                                if no_faceD+1 in ptlist: 
                                    ldone    = False
                                    c        = 0
                                    clD      = Internal.getValue(gD)
                                    searchBC = True
                                    #Recherche existence CL sur la zone receuveuse, on ajoute la face en cas d'absence
                                    while searchBC:
                                        clR    = Internal.getValue(bc_R[c])
                                        if clD == clR and ldone == False:

                                            ptlistR = Internal.getNodeFromName1(bc_R[c], 'PointList')[1]
                                            if fbingo+1 not in ptlistR: 

                                                sizeR   = numpy.size(ptlistR)
                                                datap   = numpy.empty(sizeR+1, dtype=Internal.E_NpyInt)
                                                datap[0:sizeR] = ptlistR[0:sizeR]
                                                datap[sizeR]   = fbingo +1
                                                Internal.getNodeFromName1(bc_R[c], 'PointList')[1] = datap

                                            searchBC = False
                                            ldone    = True
                                        c += 1
                                        if c == len(bc_R): searchBC = False

                                    #creation du noeud BC si type inconnu sur zone receveuse
                                    if not ldone:
                                        print('type BC inconnu: creation from scratch',gD)
                                        #ptlistR = Internal.getNodeFromName1(tmp, 'PointList')[1]
                                        datap = numpy.empty(1, dtype=Internal.E_NpyInt)
                                        datap[0] = fbingo +1
                                        tmp = [ gD[0],  gD[1], [ gD[2][0] ], gD[3]]
                                        tmp[2].append( [ 'PointList',  datap, [], 'IndexArray_t'] )
                                        print('tmp',tmp)
                                        zbc_R[2]+=tmp

    #print('zoneD=', zD[0], 'elts new ', elts,'elts old ', elts_old)
    elts = elts-elts_old

    #print('layer=', no_layer, elts)
    b = T.subzone(zD, list(elts), type='elements')

    return b, elts, zbc_R

#==============================================================================
def addGhostCellsP(t, dims_woghost, list_elts, mask_elts, xyz0, no_layer):
    import Transform.PyTree as T
    import Generator.PyTree as G

    # -- Copy un arbre en gardant des references sur les numpy
    tp = Internal.copyRef(t)
    # -- converti un node en arbre. dans le cas present, ne fait rien, car node = tree: 
    #  tpp=tp
    tpp, ntype = Internal.node2PyTree(tp)

    # Dictionnaire pour ancien et nouveau nom de zone
    znameD = {}

    bases = Internal.getBases(tpp)
    zonenames =[]
    for b in bases:
        for z in b[2]:
            if z[3] == 'Zone_t': zonenames.append(z[0])

    for b in bases:
        c = 0
        for z in b[2]:
            if z[3] == 'Zone_t':
                zgc= Internal.getNodesFromType1(z  , 'ZoneGridConnectivity_t')
                coordxD = Internal.getNodeFromPath(z, 'GridCoordinates/CoordinateX')[1]
                coordyD = Internal.getNodeFromPath(z, 'GridCoordinates/CoordinateY')[1]
                coordzD = Internal.getNodeFromPath(z, 'GridCoordinates/CoordinateZ')[1]
                NG_EC_D = Internal.getNodeFromPath(z, 'NGonElements/ElementConnectivity')[1]

                gc = Internal.getNodesFromType1(zgc, 'GridConnectivity_t')
                zp = Internal.copyRef(z)
                #zp1= Internal.copyRef(z)
                c1            = 0
                nbElts_inf = dims_woghost[c][1]+1
                #print('Traitement zone',z[0],c)
                for g in gc:
                    faceListD= Internal.getNodeFromName1(g, 'PointListDonor')[1]
                    faceListR= Internal.getNodeFromName1(g, 'PointList')[1]
                    donor    = Internal.getValue(g)
                    zdonor   = Internal.getNodeFromName2(t, donor)

                    zgcD= Internal.getNodesFromType1(zdonor  , 'ZoneGridConnectivity_t')
                    gcD = Internal.getNodesFromType1(zgcD, 'GridConnectivity_t')
                    rac = 0
                    gcSearch = True
                    while  gcSearch:
                        if Internal.getValue( gcD[rac]) == z[0]:
                            gcSearch =  False
                        else: rac +=1
                    maskD  =  zonenames.index(zdonor[0])

                    if no_layer ==1 :
                        Internal._adaptNFace2Index(zdonor) 
                        Internal._adaptNGon2Index(zdonor) 

                    #if c==2:
                    # print'maskD=',maskD, 'rac=', rac, 'plage=',  mask_elts[ maskD ][ rac ], 'nolayer=',no_layer

                    layer, elts_new, zbc  = getLayer(zdonor,z,  list_elts[c][c1], mask_elts[ maskD ][ rac ],xyz0[c], no_layer, faceListD, faceListR)

                    list_elts[c][c1] = elts_new
                    #determine la plage des elements ajoute par raccord pour eviter d'avoir des cellules en double en layer2

                    if no_layer ==1 :
                        nbElts_add = len(elts_new)-1
                        mask_elts[c][c1] = [ nbElts_inf, nbElts_inf + nbElts_add ]
                        nbElts_inf      = nbElts_inf + nbElts_add +1

                    #if no_layer==1:
                    # print'mask: zone= ',c, 'rac=', c1, 'plage=',  mask_elts[ c ][ c1 ]

                    zp = T.join(zp, layer)
                    zp = G.close(zp)
                    zp[2]+= zgc
                    zp[2]+=[zbc]
                    c1   += 1

                    Internal._adaptNFace2PE(zp, remove=False)
                    NG_PE = Internal.getNodeFromPath(zp, 'NGonElements/ParentElements')[1]

                znameD[ z[0] ] = zp[0]
                Internal._adaptNFace2PE(zp, remove=False)
                Internal._adaptNFace2Index(zp) 
                Internal._adaptNGon2Index(zp) 
                b[2][c] = zp
            c += 1

    # on renome les nom de zones dans les connectivite
    zones = Internal.getZones(tpp)
    for z in zones:
        gc = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for g in gc:
            old_name = Internal.getValue(g)
            #print('oldname ',  old_name, 'et newname: ',znameD[ old_name ])
            Internal.setValue(g, znameD[ old_name ] )

    tp = Internal.pyTree2Node(tpp, ntype)
    return tp

#==============================================================================
def adapt2FastP(t, nlayers=2):

    dims_woghost=[]
    list_elts = []
    mask_elts = []
    xyz0 = []
    for layer in range(1,nlayers+1):

        zones = Internal.getZones(t)
        #sauvegarde taille de zone sans ghost et de la list des elememnt ajoutes
        for z in zones:
        #Nombre de face totale
            NGON_range = Internal.getNodeFromPath(z, 'NGonElements/ElementRange')[1]
            nface_tot  = NGON_range[1]-NGON_range[0] +1
            #Nombre de face raccord
            gc         = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            nface_rac = 0
            elts =[]
            mask =[]
            for g in gc:
                elts.append( set() )
                mask.append( None )
                nface_rac = nface_rac + numpy.size( Internal.getNodeFromName1(g, 'PointListDonor')[1] )
            #Nombre de face BC      
            gc         = Internal.getNodesFromType2(z, 'BC_t')
            nface_bc  = 0
            for g in gc:
                nface_bc = nface_bc + numpy.size( Internal.getNodeFromName1(g, 'PointList')[1] )

            #print('zone=', z[0],'Nface total=', nface_tot, 'Nface_rac=', nface_rac, 'Nface_bc=', nface_bc, 'layer=', layer -1)

            list_elts.append( elts )
            mask_elts.append( mask )
            xyz0.append( set() )  

            Nvtx  = z[1][0][0]
            Nelts = z[1][0][1]
            dims_woghost.append( [ Nvtx , Nelts,  nface_tot, nface_rac, nface_bc ] )

        #creation nouvel arbre avec ghost
        t = addGhostCellsP(t, dims_woghost, list_elts, mask_elts, xyz0, layer )
        #print 'LAYER=',layer
        #for list_el in list_elts:
        #   print 'list-el', list_el
        #print 'dim_wo apr', dims_woghost

    zones = Internal.getZones(t)
    c = 0
    nbzone = len(zones)
    for z in zones: 
            #print 'verif',z[0]
        NGON_range = Internal.getNodeFromPath(z, 'NGonElements/ElementRange')[1]
        nface_tot  = NGON_range[1]-NGON_range[0] +1
        #Nombre de face BC      
        gc         = Internal.getNodesFromType2(z, 'BC_t')
        nface_bc  = 0
        for g in gc:
            nface_bc = nface_bc + numpy.size( Internal.getNodeFromName1(g, 'PointList')[1] )

        #print('zone=', z[0],'Nface total=', nface_tot, 'Nface_rac=', 0 , 'Nface_bc=', nface_bc, 'layer=', 2)

        Nvtx  = z[1][0][0]
        Nelts = z[1][0][1]
        dims_woghost.append( [ Nvtx , Nelts,  nface_tot, 0, nface_bc ] )

        data_ng = numpy.zeros(6, dtype=Internal.E_NpyInt)
        data_ng[0] = dims_woghost[c][2] - dims_woghost[c][3] - dims_woghost[c][4] # couche 0: faceInterne= nbface tot(0) -nface_rac(0)-nfacebc(0)
        data_ng[1] = dims_woghost[c][3]                                           # couche 0: nbface_rac(0)
        data_ng[2] = dims_woghost[c][4]                                           # couche 0: nbfacebc(0)
        l1 = c + nbzone
        l2 = c + nbzone*2
        if (nlayers ==2): 
            data_ng[4] = dims_woghost[l2][4] - dims_woghost[c ][4]                    # couche 1: nbfacebc(1)
            data_ng[3] = dims_woghost[l1][2] - dims_woghost[c ][2] - data_ng[4]       # couche 1: faceInterne(1)= nbfacetot -nfacebc(1) 
            data_ng[5] = dims_woghost[l2][2] - dims_woghost[l1][2]                    # couche 2: nbintern(2)

        node =  Internal.getNodeFromName1(z, 'NGonElements')
        Internal.createUniqueChild(node, 'IntExt', 'DataArray_t', data_ng)

        data_nf = numpy.zeros(3, dtype=Internal.E_NpyInt)
        data_nf[0] = dims_woghost[c][1]                                           # couche 0: nombre element
        l1 = c + nbzone
        l2 = c + nbzone*2
        data_nf[1] = dims_woghost[l1][1] -dims_woghost[c][1]                      # couche 1: nombre element
        if (nlayers ==2): 
            data_nf[2] = dims_woghost[l2][1] -dims_woghost[l1][1]                     # couche 2: nombre element
        node =  Internal.getNodeFromName1(z, 'NFaceElements')
        Internal.createUniqueChild(node, 'IntExt', 'DataArray_t', data_nf)

        #print('zone=', z[0], 'Elts0=', data_nf[0], 'Elts1=', data_nf[1],'Elts2=', data_nf[2])
        c +=1

    #print 'dim_wo final ', dims_woghost

    c = 0
    for z in zones:

        dim = Internal.getZoneDim(z)

        NG_PE = Internal.getNodeFromPath(z, 'NGonElements/ParentElements')[1]    
        NG_EC = Internal.getNodeFromPath(z, 'NGonElements/ElementConnectivity')[1]    
        NG_IDX= Internal.getNodeFromPath(z, 'NGonElements/FaceIndex')[1]    
        FA_EC = Internal.getNodeFromPath(z, 'NFaceElements/ElementConnectivity')[1]    
        FA_IDX= Internal.getNodeFromPath(z, 'NFaceElements/ElementIndex')[1]   

        zbc  = Internal.getNodeFromType1(z  , 'ZoneBC_t')
        bcs  = Internal.getNodesFromType1(zbc, 'BC_t')
        ptlist_bc=[]
        for bc in bcs: 
            ptlist_bc.append( Internal.getNodeFromName1( bc , 'PointList')[1] ) 
            #print 'zone=',z[0], ' bc=', bc[0]

        zgc  = Internal.getNodeFromType1(   z, 'ZoneGridConnectivity_t')
        gcs=[]
        if zgc is not None: gcs  = Internal.getNodesFromType1(zgc, 'GridConnectivity_t')
        ptlist_rac =[]
        ptlist_racD=[]
        for gc in gcs: 
            ptlist_rac.append( Internal.getNodeFromName1( gc , 'PointList')[1] )
            #print 'zone=',z[0], ' rac =', gc[0]
            donor      = Internal.getValue(gc)
            zdonor     = Internal.getNodesFromName1(zones, donor)[0]
            zgc_D  = Internal.getNodeFromType1( zdonor, 'ZoneGridConnectivity_t')
            gcs_D  = Internal.getNodesFromType1( zgc_D, 'GridConnectivity_t')
            for gc_D in gcs_D: 
                receiver      = Internal.getValue(gc_D)
                if (receiver == z[0]):
                    ptlist_racD.append( Internal.getNodeFromName1( gc_D , 'PointListDonor')[1] )
                    #print 'zone=',zdonor[0], ' racD=', gc[0]



        FA_intext=  Internal.getNodeFromPath(z,'NFaceElements/IntExt')[1]
        FA_RG    = Internal.getNodeFromPath(z, 'NFaceElements/ElementRange')[1]    
        NG_intext=  Internal.getNodeFromPath(z,'NGonElements/IntExt')[1]
        NF       = Internal.getNodeFromPath(z, 'NFaceElements/ElementConnectivity')
        nf       = converter.adapt2FastP(NG_EC, FA_EC, NG_PE,  NG_intext, FA_intext, ptlist_bc, ptlist_rac, ptlist_racD, dim[2])
        NF[1]    = nf
        #Mise a jour du nombre d'elemnt
        print('FA_RG[1]', FA_RG[1], 'NG_intext[2]',NG_intext[2],'NG_intext[4]', NG_intext[4])
        print('voir christophe pour probleme visu')
        FA_RG[1] = FA_RG[1] + NG_intext[2] + NG_intext[4]
        z[1][0,1]= z[1][0,1]+ NG_intext[2] + NG_intext[4]  
        Internal._adaptNFace2Index(z) 

    return t

#==============================================================================
def adapt2FastP2(t, nlayers=2):

    dims_woghost=[]

    for layer in range(1,nlayers+1):

        zones = Internal.getZones(t)
        #sauvegarde taille de zone sans ghost et de la list des elements ajoutes
        for z in zones:
        #Nombre de face totale
            NGON_range = Internal.getNodeFromPath(z, 'NGonElements/ElementRange')[1]
            nface_tot  = NGON_range[1]-NGON_range[0] +1
            #Nombre de face raccord
            gc         = Internal.getNodesFromType2(z, 'GridConnectivity_t')
            nface_rac = 0

            for g in gc:
                nface_rac = nface_rac + numpy.size( Internal.getNodeFromName1(g, 'PointListDonor')[1] )
            #Nombre de face BC      
            gc         = Internal.getNodesFromType2(z, 'BC_t')
            nface_bc  = 0
            for g in gc:
                nface_bc = nface_bc + numpy.size( Internal.getNodeFromName1(g, 'PointList')[1] )


            Nvtx  = z[1][0][0]
            Nelts = z[1][0][1]
            dims_woghost.append( [ Nvtx , Nelts,  nface_tot, nface_rac, nface_bc ] )

            print('zone=', z[0],'Nface total=', nface_tot, 'Nface_rac=', nface_rac, 'Nface_bc=', nface_bc, 'layer=', layer -1, Nvtx, Nelts)


    # creation nouvel arbre avec ghost
    t = addGhostCellsNG(t, nlayers)

    # 
    zones = Internal.getZones(t)
    c = 0
    nbzone = len(zones)
    for z in zones: 
            #print 'verif',z[0]
        NGON_range = Internal.getNodeFromPath(z, 'NGonElements/ElementRange')[1]
        nface_tot  = NGON_range[1]-NGON_range[0] +1
        #Nombre de face BC      
        gc        = Internal.getNodesFromType2(z, 'BC_t')
        nface_bc  = 0
        for g in gc:
            nface_bc = nface_bc + numpy.size( Internal.getNodeFromName1(g, 'PointList')[1] )

        '''
        #print 'zone=', z[0],'Nface total=', nface_tot, 'Nface_rac=', 0 , 'Nface_bc=', nface_bc, 'layer=', 2

        Nvtx  = z[1][0][0]
        Nelts = z[1][0][1]
        dims_woghost.append( [ Nvtx , Nelts,  nface_tot, 0, nface_bc ] )

        data_ng = numpy.zeros(6, Internal.E_NpyInt)
        data_ng[0] = dims_woghost[c][2] - dims_woghost[c][3] - dims_woghost[c][4] # couche 0: faceInterne= nbface tot(0) -nface_rac(0)-nfacebc(0)
        data_ng[1] = dims_woghost[c][3]                                           # couche 0: nbface_rac(0)
        data_ng[2] = dims_woghost[c][4]                                           # couche 0: nbfacebc(0)
        l1 = c + nbzone
        l2 = c + nbzone*2
        if nlayers == 2: 
           data_ng[4] = dims_woghost[l2][4] - dims_woghost[c ][4]                    # couche 1: nbfacebc(1)
           data_ng[3] = dims_woghost[l1][2] - dims_woghost[c ][2] - data_ng[4]       # couche 1: faceInterne(1)= nbfacetot -nfacebc(1) 
           data_ng[5] = dims_woghost[l2][2] - dims_woghost[l1][2]                    # couche 2: nbintern(2)
   
        node =  Internal.getNodeFromName1(z, 'NGonElements')
        Internal.createUniqueChild(node, 'IntExt', 'DataArray_t', data_ng)

        data_nf = numpy.zeros(3, Internal.E_NpyInt)
        data_nf[0] = dims_woghost[c][1]                                           # couche 0: nombre element
        l1 = c + nbzone
        l2 = c + nbzone*2
        data_nf[1] = dims_woghost[l1][1] -dims_woghost[c][1]                      # couche 1: nombre element
        if nlayers == 2: 
           data_nf[2] = dims_woghost[l2][1] -dims_woghost[l1][1]                     # couche 2: nombre element
        node =  Internal.getNodeFromName1(z, 'NFaceElements')
        Internal.createUniqueChild(node, 'IntExt', 'DataArray_t', data_nf)

        print('zone=', z[0], 'Elts0=', data_nf[0], 'Elts1=', data_nf[1],'Elts2=', data_nf[2])
        '''
        c += 1

    return t

#===============================================================================
# Add ghost cells in a pyTree
# Returns a pyTree with it zones extended with ghost cells
# IN: t: top tree
#===============================================================================
def addGhostCellsNG(t, nlayers=2):

    bases = Internal.getBases(t) 
    zones = Internal.getZones(t)
    nbz = len(zones)

    # merge des BC de meme type
    for z in zones:
        bc_type={}
        bcs   = Internal.getNodesFromType2(z, 'BC_t')
        for bc in bcs:
            tmp = Internal.getValue(bc)
            if tmp in bc_type:
                bc_type[tmp].append(bc)
            else:
                bc_type[tmp]=[bc]
        for key in bc_type:
            #print(key, len(bc_type[key]))
            if len(bc_type[key]) != 1:
                size_fen = 0
                min_face= 10000000
                max_face=-10000000
                for bc in bc_type[key]:
                    ptlist  = Internal.getNodeFromName1(bc, 'PointList')[1]
                    min_face= min(min_face, numpy.amin(ptlist)) 
                    max_face= max(max_face, numpy.amax(ptlist)) 
                    size_fen+= numpy.size(ptlist)
                    #print("minmax", min_face, max_face, size_fen)
                #bc contigu
                if min_face + size_fen - 1 ==  max_face:
                    ptlistNew = numpy.arange(min_face, min_face + size_fen, dtype=Internal.E_NpyInt).reshape((1,size_fen))
                    #on supprime les noeuds merger et on gonfle le numpy du premier noeud
                    c=0
                    for bc in bc_type[key]:
                        if c==0:
                            Internal.getNodeFromName1(bc, 'PointList')[1] = ptlistNew
                        else:
                            Internal._rmNodesByName(z, bc[0])
                        c+=1
                else:
                    print(key,": cette CL n'est pas contigue et doit etre mergee")
                    ptlistNew = numpy.empty((1,size_fen), dtype=Internal.E_NpyInt)
                    c = 0
                    ipt_ptlist = 0
                    for bc in bc_type[key]:
                        if c==0:
                            ptlistOld = Internal.getNodeFromName1(bc, 'PointList')[1]
                            ptlistNew[:,ipt_ptlist:ipt_ptlist+numpy.size(ptlistOld)]= ptlistOld[:]
                            Internal.getNodeFromName1(bc, 'PointList')[1] = ptlistNew
                            ipt_ptlist += numpy.size(ptlistOld)
                        else:
                            ptlistOld = Internal.getNodeFromName1(bc, 'PointList')[1]
                            ptlistNew[:,ipt_ptlist:ipt_ptlist+numpy.size(ptlistOld)]= ptlistOld[:]
                            ipt_ptlist += numpy.size(ptlistOld)
                            Internal._rmNodesByName(z, bc[0])
                        c+=1

    ozones = []
    znames = []
    # Overall BC names and type
    BCNames = []
    BCTypes = []
    JNames = []
    JTypes = []

    # zone name to id
    name2id = dict()
    i = 0
    for z in zones :
        name2id[z[0]] = i
        i += 1

    zone_ngons = []
    F2Es = []
    PointLists = []
    donIds = []
    ptL_sizes = []
    PointListsD = []

    bc_PointLists = []
    bc_ptL_sizes = []

    nbz_unstruct=0
    for z in zones:
        zname = z[0]
        znames.append(zname)

        if Internal.getZoneType(z)==1: continue

        m = PyTree.getFields(Internal.__GridCoordinates__, z)[0]
        zone_ngons.append(m)

        F2Esep = Internal.getNodeFromName2(z, 'ParentElements')
        if F2Esep is not None: F2Esep=F2Esep[1]

        F2Es.append(F2Esep)

        raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        nb_racs = len(raccords)

        bcs = Internal.getNodesFromType2(z, 'BC_t')
        nb_bcs = len(bcs)

        z_PointLists    = []
        z_PointListsD   = []
        z_bc_PointLists = []
        z_donIds        = numpy.empty(nb_racs, Internal.E_NpyInt)
        z_ptL_sizes     = numpy.empty(nb_racs, Internal.E_NpyInt)
        z_bc_ptL_sizes  = numpy.empty(nb_bcs,  Internal.E_NpyInt)

        j=0
        for rac in raccords:

            rt = Internal.getNodeFromType1(rac, 'GridConnectivityType_t')
            jn = "".join(Internal.getValue(rt))
            JNames.append(rac[0])
            JTypes.append(jn)

            donnorName = "".join(Internal.getValue(rac))

            id = name2id[donnorName]
            z_donIds[j] = id

            ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]
            sz  = len(ptList)

            z_ptL_sizes[j] = sz
            z_PointLists.append(ptList)

            ptListD = Internal.getNodeFromName1(rac, 'PointListDonor')[1][0]
            z_PointListsD.append(ptListD)

            j = j+1

        b = 0
        for bc in bcs:
            BCNames.append(bc[0])
            BCTypes.append(Internal.getValue(bc))

            #o = Internal.getNodeFromName1(bc, 'PointList')[1]
            #print(f"Internal get node : {o[0]}")
            ptList = Internal.getNodeFromName1(bc, 'PointList')[1][0]
            z_bc_PointLists.append(ptList)
            z_bc_ptL_sizes[b] = len(ptList)
            b = b+1

        PointLists.append(z_PointLists)
        donIds.append(z_donIds)
        ptL_sizes.append(z_ptL_sizes)
        PointListsD.append(z_PointListsD)

        bc_PointLists.append(z_bc_PointLists)
        bc_ptL_sizes.append(z_bc_ptL_sizes)

        nbz_unstruct+=1
    #print "nb of bc in python side as input : %d"%len(bc_ptL_sizes)

    nodes = converter.addGhostCellsNG(zone_ngons, F2Es, ptL_sizes, PointLists, PointListsD, donIds, bc_ptL_sizes, bc_PointLists, nlayers)

    ozones = []

    # Create zones with ghost upon exit

    for i in range(nbz_unstruct):

        #print "zone %d"%i

        zwgh = PyTree.convertArrays2ZoneNode(znames[i], [nodes[i][0]])
        #zwgh[0] = znames[i] # oblige de forcer, car la fonction getZoneName est appelee par convertArrays2ZoneNode

        ######################################################################################
        # FIXME : can be retrieved from K_CONVERTER::addGhostCellsNG instead of regenerating it
        Internal._adaptNFace2PE(zwgh, remove=False) 
        Internal._adaptNFace2Index(zwgh) 
        Internal._adaptNGon2Index(zwgh) 
        #######################################################################################

        # MAJ BCs
        nb_bcs = len(nodes[i][2])
        for j in range(nb_bcs):

            (oid, ids) = nodes[i][2][j]
            # print oid
            # print BCNames[oid], BCTypes[oid]
            PyTree._addBC2Zone(zwgh, BCNames[oid], BCTypes[oid], faceList=ids) #BCTypes[oid]

        ozones.append(zwgh)

    # MAJ JOINS

    # ozones must have been created to set the joins between them

    tp = PyTree.newPyTree([bases[0][0]])
    # Attache les zones

    for i in range(nbz_unstruct):

        zwgh = ozones[i]

        nb_joins = len(nodes[i][1])
        for j in range(nb_joins):

            (oid, zid, ids) = nodes[i][1][j]

            # get the corresponding join in corresponding joined ghosted zone
            for k in range(len(nodes[zid][1])):
                (donnor_oid, donnor_zid, donnor_ids) = nodes[zid][1][k]
                if (donnor_zid == i): break

            PyTree._addBC2Zone(zwgh, JNames[oid], JTypes[oid],
                  zoneDonor=ozones[zid],
                  faceList=ids ,           #elementList=elementList, elementRange=elementRange, data=data, subzone=subzone,
                  faceListDonor=donnor_ids #, elementListDonor=elementListDonor, elementRangeDonor=elementRangeDonor
                  )
        node =  Internal.getNodeFromName1(zwgh, 'NFaceElements')
        Internal.createUniqueChild(node, 'IntExt', 'DataArray_t',  nodes[i][3])
        node =  Internal.getNodeFromName1(zwgh, 'NGonElements')
        Internal.createUniqueChild(node, 'IntExt', 'DataArray_t',  nodes[i][4])

    tp[2][1][2] += ozones

    #return ozones
    return tp




#===============================================================================================================================================
#===============================================================================================================================================

def _adaptBCStruct3__(join, dimZoneR,dimZoneD, d):

    #BCMatch
    prange = Internal.getNodesFromName1(join, 'PointRange')
    direction = getDirection__(dimZoneR[4], prange)
    changePointRange__(prange, dimZoneR, direction, d, extend=0)

    prangedonor = Internal.getNodesFromName1(join, 'PointRangeDonor')

    direction = getDirection__(dimZoneD[4], prangedonor)
    changePointRange__(prangedonor, dimZoneD, direction, d, extend=0)


    return None



