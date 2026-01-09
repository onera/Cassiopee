"""Grid intersection module.
"""
#
# Python Interface to create PyTrees defining meshes
#
from . import Intersector as XOR
from . import intersector
import time

__version__ = XOR.__version__

import numpy

try:
    import Converter.PyTree as C
    import Converter.Distributed as CD
    import Converter.Internal as Internal
    import Converter
    import Transform.PyTree as T
    import Generator.PyTree as G
    import Post.PyTree as P
except:
    raise ImportError("Intersector.PyTree: requires Converter.PyTree module.")

#=============================================================================
# Ajout d'un champ 'var' pour une zone qui ne le contient pas.
#=============================================================================
def addVar__(t, var, loc='centers'):
    tp = Internal.copyRef(t)
    _addVar__(tp, var, loc)
    return tp

def _addVar__(t, var, loc='centers'):
    if loc == 'centers': var = loc + ':'+ var
    C._addVars(t, var)
    return None

#=============================================================================
# Retourne le nombre de cellules d'un maillage (voir C.getNCells)
#=============================================================================
def nb_cells(a):
    ncellsTot = 0
    zones = Internal.getNodesFromType2(a, 'Zone_t')
    for z in zones:
        dim = Internal.getZoneDim(z)
        ncells = dim[2]
        ncellsTot += ncells
    return ncellsTot

#==============================================================================
# nb_faces: Returns the number of faces in t
# IN: t: 3D NGON mesh
# OUT: returns the adapted feature
#==============================================================================
def nb_faces(t):

    ncfacesTot = 0
    zones = Internal.getNodesFromType2(t, 'Zone_t')
    for z in zones:
        GEl = Internal.getElementNodes(z)
        NGON = 0; found = False
        for c in GEl:
            if c[1][0] == 22: found = True; break
            NGON += 1
        if found:
            node = GEl[NGON]
            er = Internal.getNodeFromName1(node, 'ElementRange')
            ncfacesTot += er[1][1]-er[1][0]+1
    return ncfacesTot

#==============================================================================
# updateNugaData: update a field by applying an oid
# IN: t: 3D NGON mesh
# OUT: returns the adapted feature
#==============================================================================
def updateNugaData(field, oid):
    new_fld = numpy.empty(len(oid), Internal.E_NpyInt)
    i = -1
    for oi in oid:
        i += 1
        new_fld[i] = field[oi]
    return new_fld

#==============================================================================
# getTreeDim: XXX
# IN: t: 3D NGON mesh
# OUT: XXX
#==============================================================================
def getTreeDim(t):
    zs = Internal.getZones(t)
    d = 0
    nb_elts_per_zone = 0
    for z in zs:
        dims = Internal.getZoneDim(z)
        if d == 0: d = dims[4]
        if dims[4] != d: return 0 # mixed type : not handled
        if dims[3] == 'NGON': nb_elts_per_zone += nb_cells(z)

    if d == 3 and nb_elts_per_zone == len(zs): d = 21 # NGON with one cell per zone => NUGA NGON => 2D

    return d

#==============================================================================
# getZoneNSTypeAndDim: XXX
# IN: t: 3D NGON mesh
# OUT: XXX
#==============================================================================
def getZoneNSTypeAndDim(z):
    dims = Internal.getZoneDim(z)
    d = dims[4]
    if d == 2:
        if dims[3] == 'NGON': return ('NGON_CASSIOPEE', 2)
        else: return ('BASIC', 2)
    else: # 3D
        #print(nb_cells(z))
        if dims[3] == 'NGON' and nb_cells(z) == 1: return ('NGON_NUGA',2)
        elif dims[3] == 'NGON': return ('NGON', 3)
        else: return ('BASIC', 3)

#==============================================================================
# InputType: XXX
# IN: t: 3D NGON mesh
# OUT: XXX
#==============================================================================
# return values : 0(single numpy), 1(list of numpies), 2(single zone), 3 (PyTree), 4(list of zones), 5(list of list of numpies)
def InputType(t): # fixme : based on first block only
    if isinstance(t, list):
        if isinstance(t[0], numpy.ndarray): return 1
        # note (Imad): on check que t[0] est une liste et que le premier élément est un numpy
        if isinstance(t[0], list):
            if isinstance(t[0][0], numpy.ndarray): return 5
        isnod = Internal.isStdNode(t)
        #print(isnod)
        if isnod == -1 or isnod == 0:
            #print('is std node')
            if XOR.isSingleZone(t):
                #print('is isSingleZone node')
                return 2
            if Internal.isTopTree(t):
                #print('is isSingleZone node')
                return 3
            if XOR.isSingleZone(t[0]):
                #print('is isSingleZone t0')
                return 4
        else: return -1
    else:
        if isinstance(t, numpy.ndarray): return 0
        else: return -1 # error

#==============================================================================
# NGONBlock

# GOAL: converts into a single connex NGON zone
# IN: t             : 3D NGON mesh
# IN: TOL           : tolerance
# IN: nb_comps      : nb of connex parts
# OUT: returns the adapted feature
#==============================================================================
def NGONBlock(t, nb_comps, mixed_type=False, keep_BC=False, tol=0.):

    if keep_BC == False:
        C._deleteGridConnectivity__(t)
        C._deleteZoneBC__(t)

    if mixed_type:
        t = C.breakConnectivity(t)
        for zone in Internal.getZones(t):
            if Internal.getZoneDim(zone)[3] in ['TRI','QUAD']: Internal.rmNode(t,zone)

    t = C.convertArray2NGon(t, recoverBC=keep_BC)

    #C.convertPyTree2File(t, 'tNG.cgns')
    # compute relevant tolerance : 1% of the minimum edge length
    TOL = tol
    if TOL == 0.:
        TOL = 0.1 * edgeLengthExtrema(t)

    t = concatenate(t, tol=TOL)

    valid = True
    if nb_comps > 0:
        valid = isConformalNConnex(t, nb_comps)
    if not valid:
        #C.convertPyTree2File(t, 'bad_oper.cgns')
        raise ValueError('NGONBlock: concatenation failed to produce a single connex and conformal block.')

    #Internal._correctPyTree(t) # to add families if required and missing
    zs = Internal.getZones(t)
    z = zs[0]
    return z

#==============================================================================
# isConformalNConnex
# GOAL: check if a component is conformal and connex
# IN: t         : 3D NGON mesh
# IN: nb_comps  : nb of connex parts
# OUT: returns the adapted feature
#==============================================================================
def isConformalNConnex(t, nb_comps):

    F = P.exteriorFaces(t)
    F = T.splitConnexity(F)
    zs = Internal.getZones(F)
    if len(zs) != nb_comps:
        print('nb of contours: '+str(len(zs)))
        print ('expected to be: '+str(nb_comps))
        return False
    F = P.exteriorFaces(F)
    zs = Internal.getZones(F)
    #print('nb of edge zones: ' + str(len(zs)))
    if len(zs) > 0:
        if nb_faces(zs[0]) != 0: return False
    return True


#==============================================================================
# computeMeshAssemblyParameters

# GOAL: Computes the tolerance (based on min edge length) and thr configuration span to nomalize it
# IN: t: 3D NGON mesh
# OUT: returns the tolerance and the max bounding box size
#==============================================================================
def computeMeshAssemblyParameters(t):
    t1 = C.convertArray2NGon(t)
    L = edgeLengthExtrema(t1) # min edge length
    tol = 0.05*L
    box = G.bbox(t)
    dx = abs(box[3] - box[0])
    dy = abs(box[4] - box[1])
    dz = abs(box[5] - box[2])
    H = max(dx, max(dy,dz))
    print('Computed TOL (on initial config. size) : ' +str(tol))
    print('Configuration span : '+str(H))
    return (tol, H)

#==============================================================================
# computeVmin

# GOAL : Computes the minimum cell volume in a liste of zones

# IN: t: 3D NGON mesh

# OUT: returns the tolerance and the max bounding box size
#==============================================================================
def computeVmin(zones):
    import sys;
    vmin=sys.float_info.max

    for z in zones:
        tmp = Internal.createElsaHybrid(z, method=1, methodPE=1)
        (vm, cellid, zoneid) = checkCellsVolume(tmp)
        vmin = min(vmin, vm)

    return vmin

def computeMetrics(NGzones):
    import sys
    MEL = GRmin = Vmin = sys.float_info.max

    for z in NGzones:

        zi = Internal.createElsaHybrid(z, method=1, methodPE=1)
        (Vmini, ivmin, vzoneid, GRmini, igrmin, grzoneid) = checkCellsVolumeAndGrowthRatio(zi)

        MELi = edgeLengthExtrema(zi)

        MEL = min(MEL, MELi)
        Vmin = min(Vmin, Vmini)
        GRmin = min(GRmin, GRmini)

    return (0.01*MEL, Vmin, GRmin)

#==============================================================================
# checkAssemblyForlSolver

# GOAL : Checks that the solver will be happy with the mesh

# IN: t             : 3D NGON mesh

# OUT: messages
#==============================================================================
def checkAssemblyForlSolver(t, fullcheck=False, nb_comps=1):

    MAXFLUXVAL=1.e-9
    VOLMIN_SOLVER = 1.e-20

    # VERIFICATION 1 : CONFORMITE
    conformal = isConformalNConnex(t, nb_comps)
    if not conformal:
        raise ValueError('checkAssemblyForSolver: non conformal mesh.')
    else:
        print('OK: conformal.')

    # VERIFICATION 2 : CELLULES FERMEES
    print("Check cells closure ...")
    err = checkCellsClosure(t)
    if err == 1:
        print('Boolean ERROR : open cells')
        import sys; sys.exit()
    # else :
    #   print ('OK : all cells are closed')

    # VERIFICATION 3 : ORIENTATION/QUALITE
    print("Check orientation/quality ...")
    (maxflux, cellid, zoneid) = checkCellsFlux(t)
    print('max flux : '+str(maxflux))
    if maxflux > MAXFLUXVAL:
        print('Boolean WARNING : quality might be not good enough for solver')
        #import sys; sys.exit()

    if fullcheck == False: return

    # VERIFICATION 4 : VOL MIN
    print("Check min cell volume ...")
    (vmin, cellid, zoneid) = checkCellsVolume(t)
    print('vol min : ', vmin)
    if vmin < VOLMIN_SOLVER:
        print('Boolean ERROR : too small cells detected : under solver threshold')

    # print("Check cell volume extrema...")
    # res = statsSize(t, 1)
    # dMax = res[0]
    # smin = res[1]
    # smax = res[2]
    # vmin = res[3]
    # vmax = res[4]
    # print ('vmin')
    # print (vmin)
    # print ('vmax')
    # print (vmax)

    #VERIFICATION 5 : BIG CELLS
    print("Check big cells (over 40) ...")
    N=40
    big = checkForBigCells(t, N)
    ofile = 'big_'+str(N)+'.cgns'
    C.convertPyTree2File(big, ofile)
    ofile = 'biggest.cgns'
    z = extractBiggestCell(t)
    C.convertPyTree2File(z, ofile)

    #VERIDICATION 6 : check for degen
    print("Check for degen cells...")
    checkForDegenCells(t)

    #VERIF 7
    print("Check for over connected cells..")
    detectOverConnectedFaces(t)

    #VERIF 8
    print("extract pathos..")
    z=extractPathologicalCells(t)
    ofile = 'pathos.cgns'
    #C.convertPyTree2File(z, ofile)

    #VERIF 9
    print("Check for identical cells..")
    detectIdenticalCells(t, TOL=1.e-10, clean=0)

#=============================================================================
# Concatenation des PointList d'un type de BC donne dans une liste de zones
#=============================================================================
def concatenateBC(bctype, zones, wallpgs, cur_shift):
    i = 0
    for z in zones:
        c = C.getFields(Internal.__GridCoordinates__, z)
        if c == []: continue
        if len(c[0]) == 5: continue # structure (skip)

        #print(' -- zone : %d / %d' %(i+1, len(zones)))
        i = i+1
        bnds = Internal.getNodesFromType(z, 'BC_t')
        #print(" -- this zone has %d boundaries"%(len(bnds)))
        #print(' -- cur shift %d' %(cur_shift))

        # GET THE WALL PGS FROM THE POINTLISTS
        for bb in bnds:
            if Internal.isValue(bb, bctype) == False: continue

            wpgs = bb[2][1][1][0] # POINTLIST NUMPY fixme : can be somewhere else in the array
            #print(wpgs)
            # SYNC THE POINTLIST BEFORE APPENDING  : SHIFT WITH THE CURRENT NB OF STORED POLYGONS
            id2 = numpy.empty(len(wpgs), Internal.E_NpyInt)
            id2[:] = wpgs[:] + cur_shift
            wallpgs.append(id2)

        c = c[0]
        #z_nb_pts= len(c[1][0])
        z_nb_pgs = c[2][0][0]
        #print(z_nb_pts)
        #print(z_nb_pgs)
        cur_shift += z_nb_pgs
    return (wallpgs, cur_shift)

# update BC and JOINS point lists given an indirection "new id to old id"
def updatePointLists(z, zones, oids):
    if len(oids) == 0:
        C._deleteZoneBC__(z)
        C._deleteGridConnectivity__(z)
        return
    updateBCPointLists1(z, oids)
    updateJoinsPointLists1(z, zones, oids)

def updateBCPointLists1(z, oids):
    bnds = Internal.getNodesFromType(z, 'BC_t')
    zname = z[0]

    ptLists = []
    for bb in bnds:
        ptl = Internal.getNodesFromType(bb, 'IndexArray_t')
        ptLists.append(ptl[0][1][0])

    if ptLists == []: return

    # recalcul des pointlist
    ptLists = XOR.updatePointLists(oids, ptLists)

    i=0
    bc_to_remove=[]
    #print('update the BC pointlists')
    for bb in bnds:
        torem = False
        if isinstance(ptLists[i], numpy.ndarray):
            if len(ptLists[i]) == 0: torem = True
        elif ptLists[i] is None: torem = True

        if torem:
            bc_to_remove.append(Internal.getPath(z, bb))
            continue

        ptl = Internal.getNodesFromType(bb, 'IndexArray_t')
        ptLists[i] = ptLists[i].reshape(1, ptLists[i].size)
        ptl[0][1] = ptLists[i]
        #print('shape dans updateBCPointLists1 : '+str(numpy.shape(ptl[0][1])))
        #print(ptl[0][1])
        i=i+1

    for p in bc_to_remove: Internal._rmNodeFromPath(z, p)

def updateJoinsPointLists1(z, zones, oids):

    joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
    zname=z[0]

    ptLists = []
    for j in joins:
        ptl = Internal.getNodeFromName1(j, 'PointList')
        ptLists.append(ptl[1])#ptl[0][1][0]

    if ptLists == []: return

    # recalcul des pointlist
    ptLists = XOR.updatePointLists(oids, ptLists)

    i=0
    # update the Join pointlist and synchronize with other zones (their pointListDonnor)
    for j in joins:
        donnorName = "".join(Internal.getValue(j))

        # Match has disappeared > remove from tree
        if not isinstance(ptLists[i], numpy.ndarray):
            Internal._rmNode(z, j)
            i=i+1
            continue

        ptl = Internal.getNodeFromName1(j, 'PointList')
        # print("donnorName:", donnorName)
        dz = Internal.getNodeFromName(zones, donnorName)

        if dz is not None:
            joinsD = Internal.getNodesFromType(dz, 'GridConnectivity_t')
            for jd in joinsD:
                #dname = "".join(jd[1])
                dname = "".join(Internal.getValue(jd))
                # print("dname / zname : ", dname, zname)
                if (dname != zname) : continue
                ptlD = Internal.getNodeFromName1(jd, 'PointListDonor')

                PG0 = ptl[1][0][0] # first polygon in the point list
                PG0D = ptlD[1][0][0] # first polygon in the point list
                if (PG0 != PG0D) : continue # not the right join (in case of multiple joins for 2 zones) : the first PG must be the same (assume one PG only in one join)

                ptLists[i] = numpy.reshape(ptLists[i], (1,len(ptLists[i])))

                ptl[1]  = ptLists[i]
                ptlD[1] = ptLists[i]

                break
            i=i+1

def getJoinsPtList(z, zname2id):
	raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')
	nb_racs = len(raccords)

	jzid=[]
	jptlist=[]

	j=0
	for rac in raccords:
		rt = Internal.getNodeFromType1(rac, 'GridConnectivityType_t')

		donnorName = "".join(Internal.getValue(rac))
		#print(donnorName)
		zid = zname2id[donnorName]
		#print('id is ' + str(id))
		ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]
		#print (ptList)
		sz  = len(ptList)

		jzid.append(zid)
		jptlist.append(ptList)

		j = j+1

	return (jzid, jptlist)

#==============================================================================
# wrapper to catch errors
#==============================================================================
def getProperty(z, pname):
    val = CD.getProperty(z, pname)
    #if val == -1 :
    #  print('ERROR : ' + str(pname) + ' is not set. You must call _setZonesAndJoinsUId on the input tree')
    return val

#==============================================================================
# getBCsPtLists : XXX
# IN: t : 3D NGON PyTree
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical zone hook
#==============================================================================
def getBCsPtLists(t):

    zone_to_bcptlists = {}
    zones = Internal.getZones(t)
    #
    for z in zones:
        zid = getProperty(z, 'zid')
        # BC and Joins point list (owned side)
        bcptlists = getBCPtList(z)
        zone_to_bcptlists[zid]=bcptlists

    return zone_to_bcptlists


def getBCPtList(z):
    #
    bnds = Internal.getNodesFromType(z, 'BC_t')

    ptList = []
    for bb in bnds :
        ptlnod = Internal.getNodesFromType(bb, 'IndexArray_t')
        if ptlnod == [] : continue
        x = ptlnod[0][1]

        if type(x[0]) is numpy.ndarray: # ptlnod[0][1] is a list with one ptlist : [ [...] ]
            ptList.append(x[0])
        else: # # ptlnod[0][1] is directly the ptlist : [...]
            ptList.append(x)

    return ptList

def getBCPtListOfType(z, typesList, families=[]):
    #
    nodes = []
    for btyp in typesList:
        nodes += Internal.getNodesFromValue(z, btyp)
        if families != []:nodes += C.getFamilyBCs(z, families)

    #print(nodes)
    ptList = []
    for n in nodes:
        ptlnod = Internal.getNodesFromType(n, 'IndexArray_t')
        if ptlnod == []: continue
        x = ptlnod[0][1]
        #print(x)
        if type(x[0]) is numpy.ndarray: # ptlnod[0][1] is a list with one ptlist : [ [...] ]
            ptList.append(x[0])
        else: # # ptlnod[0][1] is directly the ptlist : [...]
            ptList.append(x)

    if ptList != []: ptList = numpy.concatenate(ptList).ravel() # create a single list

    return ptList

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
def updateJoinsPointLists2(z, zones, jzids, ptLists):
    # zone name to id
    name2id = dict()
    i = 0
    for zz in zones:
        name2id[zz[0]] = i
        #print('zone ' + z[0] + ' has ' + str(i))
        i += 1

    #print('processed zone  : ' + z[0])

    joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
    zname = z[0]

    # update the Join pointlist and synchronize with other zones (their pointListDonnor)
    for j in joins:
        #donnorName = b''.join(j[1]).decode() WORK, try below
        donnorName = "".join(Internal.getValue(j))
        #print(donnorName)
        ptl = Internal.getNodeFromName1(j, 'PointList')
        dz = Internal.getNodeFromName(zones, donnorName)
        joinsD = Internal.getNodesFromType(dz, 'GridConnectivity_t')
        for jd in joinsD:
            dname = "".join(Internal.getValue(jd))
            if dname != zname: continue
            ptlD = Internal.getNodeFromName1(jd, 'PointListDonor')

            PG0 = ptl[1][0][0] # first polygon in the poitn list
            PG0D = ptlD[1][0][0] # first polygon in the poitn list
            if PG0 != PG0D: continue # not the right join (in case of multiple joins for 2 zones) : the first PG must be the same (assume one PG only in one join)

            id = name2id[donnorName]

            # print('donnor is ' + donnorName + ' with id  : ' + str(id))
            # print(ptl[1])
            # print(numpy.shape(ptl[1]))

            # find rank for this id and set list
            i=-1
            for k in jzids:
                i+=1
                if k != id: continue

                #print('new j')

                ptLists[i] = numpy.reshape(ptLists[i], (1,len(ptLists[i])))
                # print(ptLists[i])

                ptl[1]= ptLists[i]
                ptlD[1] = ptLists[i]

                break
            break

def updateBCPointLists2(z, ptLists):
    bnds = Internal.getNodesFromType(z, 'BC_t')

    i=0
    # update the BC pointlists
    for bb in bnds :
        ptl = Internal.getNodesFromType(bb, 'IndexArray_t')
        ptl[0][1] = ptLists[i]
        #print(ptl[0][1])
        i=i+1

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
def updateJoinsPointLists3(z, zidDict, rid_to_ptlist, ptlname): # 'PointList', 'PointListDonor'

    joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
    zid = getProperty(z, 'zid')

    # update the Join pointlist and synchronize with other zones (their PointListDonor)

    processed_rid = set()

    for j in joins:

        ptl    = Internal.getNodeFromName1(j, ptlname)

        rid    = getProperty(j, 'rid')
        donnorName = "".join(Internal.getValue(j))
        jzid   = zidDict[donnorName]

        if rid not in rid_to_ptlist : continue

        L1     = rid_to_ptlist[rid]

        if jzid==zid:
            half   = int(len(rid_to_ptlist[rid]) / 2)

            if ptlname == 'PointList':        # appel apres conformisation
                if rid in processed_rid:        ## deuxieme passe (i.e deuxime demi-raccord)
                    L1 = L1[half:]
                else:                           ## premiere passe
                    processed_rid.add(rid)
                    L1 = L1[:half]
            elif ptlname == 'PointListDonor': # appel apres echange de pointlist
                if rid in processed_rid:        ## deuxieme passe (i.e deuxime demi-raccord)
                    L1 = L1[:half]
                else:                           ## premiere passe
                    processed_rid.add(rid)
                    L1 = L1[half:]

        L1     = numpy.reshape(L1, (1,len(L1)))
        ptl[1] = L1

#------------------------------------------------------------------------------
#
#------------------------------------------------------------------------------
def get_transfo_to_ptlist_dico(z):

    # -----

    dico_rotation_to_ptList = {}

    key1 = [0.]*6
    key1 = tuple(key1)

    raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')

    for rac in raccords:
        rac_match_perio = Internal.getNodesFromType2(rac, 'GridConnectivityProperty_t')

        ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]

        if rac_match_perio: #perio BC
            center_rotat  = Internal.getNodesFromName3(rac, 'RotationCenter')[0][1]
            angle_rotat   = Internal.getNodesFromName3(rac, 'RotationAngle' )[0][1]

            # Condition to set negative angle as a translation like dict

            key2 = numpy.concatenate((center_rotat, angle_rotat))
            key2 = tuple(key2)

            if numpy.sum(angle_rotat) < 0: #negative rotation
                if key1 in dico_rotation_to_ptList:
                    dico_rotation_to_ptList[key1] = numpy.concatenate((dico_rotation_to_ptList[key1], ptList))
                else:
                    dico_rotation_to_ptList[key1] = ptList
            else: #translation or positive rotation
                if key2 in dico_rotation_to_ptList:
                    dico_rotation_to_ptList[key2] = numpy.concatenate((dico_rotation_to_ptList[key2], ptList))
                else:
                    dico_rotation_to_ptList[key2] = ptList
        else: #match BC
            if key1 in dico_rotation_to_ptList: #following match BC
                dico_rotation_to_ptList[key1] = numpy.concatenate((dico_rotation_to_ptList[key1], ptList))
            else: #first match BC
                dico_rotation_to_ptList[key1] = ptList

    return dico_rotation_to_ptList




#------------------------------------------------------------------------------
# Conformisation d'une soupe de TRI ou de BAR
#------------------------------------------------------------------------------
def conformUnstr(surface1, surface2=None, tol=0., left_or_right=0, itermax=10):
    """Conformizes a TRI or BAR soup (surface1) with taking into account surface2 if it's provided.
    Usage: conformUnstr(s1, s2, tol, left_or_right, itermax)"""
    s1 = C.getFields(Internal.__GridCoordinates__, surface1)[0]
    if surface2 is not None:
        s2 = C.getFields(Internal.__GridCoordinates__, surface2)[0]
    else: s2 = None
    s = XOR.conformUnstr(s1, s2, tol, left_or_right, itermax)
    return C.convertArrays2ZoneNode('conformized', [s])

def intersection(surface1, surface2, tol=0., itermax=10):
    """Computes the intersection between two closed-surface or two volume meshes.
    Usage: intersection(s1, s2, tol)"""
    s1 = C.getFields(Internal.__GridCoordinates__, surface1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, surface2)[0]
    s = XOR.intersection(s1, s2, tol, itermax)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the intersection between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanIntersection(a1, a2, tol)
    Usage for volumes: booleanIntersection(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanIntersection(s1, s2, tol, preserve_right, solid_right, agg_mode, improve_qual, itermax)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanUnion(a1, a2, tol=0., jtol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, multi_zone=False, simplify_pgs=True, hard_mode=0, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the union between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanUnion(a1, a2, tol)
    Usage for volumes: booleanUnion(a1, a2, tol, preserve_right, solid_right)"""

    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]

    if multi_zone:
        typzone1 = s1[3]
        typzone2 = s2[3]
        if typzone1 == 'NGON' and typzone2 == 'NGON': # only for Volume/Volume
            # compute the join tolerance
            if jtol == 0.:
                L1 = edgeLengthExtrema(a1)
                L2 = edgeLengthExtrema(a2)
                jtol = 0.01*min(L1,L2)
                #print('jtol is : '+str(jtol))

            return booleanUnionMZ(a1, a2, tol, jtol, agg_mode, improve_qual, simplify_pgs, hard_mode)

    #multi_zone option is ignored from here

    cur_shift=0
    extrudepgs=[]
    if (solid_right == 1) :
        zones = Internal.getZones(a2)
        (extrudepgs, cur_shift) = concatenateBC('UserDefined', zones, extrudepgs, cur_shift)
    if (extrudepgs != []) : extrudepgs = numpy.concatenate(extrudepgs) # create a single list
    #print("nb of pgs to pass : %s" %(len(extrudepgs)))

    res = XOR.booleanUnionWithHisto(s1, s2, tol, preserve_right, solid_right, agg_mode, improve_qual, extrudepgs, simplify_pgs, hard_mode, itermax)

    is_zone_list  = 0
    if (len(res) != 5) : is_zone_list = 1
    elif (res[0][3] != 'NGON' and res[0][3] != 'TRI' and res[0][3] != 'BAR') : is_zone_list = 1
    if (len(res) == 1) : is_zone_list = 2 # Not NGON and not DEBUG

    if is_zone_list  == 0:

        phnids1 = res[1]
        phnids2 = res[2]
        pgoids1 = res[3] # warning pgoids (not pgnids)
        pgoids2 = res[4] # warning pgoids (not pgnids)

        # Restore BCs
        # ===========
        newz    = C.convertArrays2ZoneNode('union', [res[0]])

        BC1s    = Internal.getNodesFromType(a1, "BC_t")
        BC2s    = Internal.getNodesFromType(a2, "BC_t")

        if (BC1s != []) or (BC2s != []):
            zoneBC  = Internal.createUniqueChild(newz, "ZoneBC", 'ZoneBC_t')
            zs      = Internal.getZones(newz)
            zs      = zs[0]

            # Check name to prevent collision
            name1s  = []
            for bc1 in BC1s:
                name1s.append(bc1[0])

            name2s  = []
            for bc2 in BC2s:
                name2s.append(bc2[0])

            for name1 in name1s:
                if name1 in name2s:
                    BC1s = Internal.renameNode(BC1s, name1, name1+'b')

            # Update pointList
            for bc1 in BC1s:
                # copy to avoid modif input tree
                bc  = Internal.copyTree(bc1)
                Internal._addChild(zoneBC, bc)

            updatePointLists(zs, newz, pgoids1)

            newBC1s = Internal.getNodesFromType(newz, "BC_t")

            Internal.rmNode(newz, zoneBC)

            zoneBC  = Internal.createUniqueChild(newz, "ZoneBC", 'ZoneBC_t')
            for bc2 in BC2s:
                # copy to avoid modif input tree
                bc = Internal.copyTree(bc2)
                Internal._addChild(zoneBC, bc)

            updatePointLists(zs, newz, pgoids2)

            Internal._addChild(zoneBC, newBC1s)

            # Search for BC families
            fam1_nodes = Internal.getNodesFromType(a1, 'Family_t')
            fam2_nodes = Internal.getNodesFromType(a2, 'Family_t')

            for node in fam1_nodes:
                Internal._addChild(newz, node)

            for node in fam2_nodes:
                Internal._addChild(newz, node)


        # Restore fields
        # ==============
        # a. Build varnames list
        varnames = []

        cont1 = Internal.getNodeFromName(a1, Internal.__FlowSolutionCenters__)
        cont2 = Internal.getNodeFromName(a2, Internal.__FlowSolutionCenters__)
        cont  = [cont1,cont2]

        flds = None

        if cont1 is not None and cont2 is not None:
            flds = Internal.getNodesFromType1(cont, 'DataArray_t')
        elif cont1 is not None:
            flds = Internal.getNodesFromType1(cont1, 'DataArray_t')
        elif cont2 is not None:
            flds = Internal.getNodesFromType1(cont2, 'DataArray_t')

        if flds is not None:
            for fld in flds:
                if fld[0] not in varnames:
                    varnames.append(fld[0])

        # b. Build new var arrays
        ncells  = C.getNCells(newz)
        nsize1  = len(phnids1)
        nsize2  = len(phnids2)

        for varname in varnames:
            varnew  = numpy.zeros(ncells, numpy.float64)
            # Operande 1
            # ----------
            if cont1 is not None:
                varnode = Internal.getNodeFromName(cont1, varname)
                if varnode is not None:
                    for k in range(nsize1):
                        ik = phnids1[k]
                        if ik != -1:
                            varnew[ik] = varnode[1][k]
            # Operande 2
            # ----------
            if cont2 is not None:
                varnode = Internal.getNodeFromName(cont2, varname)
                if varnode is not None:
                    for k in range(nsize2):
                        ik = phnids2[k]
                        if ik != -1:
                            varnew[ik] = varnode[1][k]

            # Create new var
            # --------------
            C._initVars(newz, 'centers:'+varname, 0)
            contnew    = Internal.getNodeFromName(newz, Internal.__FlowSolutionCenters__)
            varnode    = Internal.getNodeFromName(contnew, varname)
            varnode[1] = varnew

        return newz

    if is_zone_list  == 2:
        return( C.convertArrays2ZoneNode('union', [res[0]]) )

    # debug : mutli zones
    ozones = []

    for i in range(len(res)):
        if len(res[i][0][1]) != 0: ozones.append(C.convertArrays2ZoneNode(res[i][1], [res[i][0]])) #(zname, array)

    return ozones

def booleanUnionMZ(t1, t2, xtol=0., jtol=0., agg_mode=1, improve_qual=False, simplify_pgs=True, hard_mode=0): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    tp1 = Internal.copyRef(t1)
    tp2 = Internal.copyRef(t2)
    return _booleanUnionMZ(tp1, tp2, xtol, jtol, agg_mode, improve_qual, simplify_pgs, hard_mode)


def _booleanUnionMZ(t1, t2, xtol=0., jtol=0., agg_mode=1, improve_qual=False, simplify_pgs=True, hard_mode=0): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the union between two closed volume meshes.
    Usage for volumes: booleanUnion2(a1, a2, tol, agg_mode)"""
    t1_is_tree = Internal.isTopTree(t1)
    t2_is_tree = Internal.isTopTree(t2)

    m1s = []
    z1s = Internal.getZones(t1)
    for z in z1s:
        m1s.append(C.getFields(Internal.__GridCoordinates__, z)[0])
    m2s = []
    z2s = Internal.getZones(t2)
    for z in z2s:
        m2s.append(C.getFields(Internal.__GridCoordinates__, z)[0])

    res = XOR.booleanUnionMZ(m1s, m2s, xtol, jtol, agg_mode, improve_qual, simplify_pgs, hard_mode)

    #if t1 and t2 does not overlap, merge
    if res == -7:
        if t1_is_tree == True and t2_is_tree == True:
            return Internal.merge([t1, t2])
        else:
            return z1s+z2s


    i=0
    paths = []
    zs = []

    iz=-1
    # New name dictionary and zone name updating
    newname = {}
    for z in z1s:
        iz += 1
        newname[z[0]] = 'dom1_z_'+str(iz)
        z[0]          = newname[z[0]]

    # Matches updating (for zone names)
    for z in z1s:
        # Update joins names
        joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
        for j in joins:
            oldname = Internal.getValue(j)
            if oldname in newname:
                Internal.setValue(j,newname[oldname])

    iz=-1
    for z in z1s:
        iz +=1

        mesh = res[i]
        pg_oids=res[i+1]
        ph_oids=res[i+2]

        # print(" ")
        # print("OP1 - pgoids: ", pg_oids)

        if mesh != []:

            # MAJ du maillage de la zone
            C.setFields([mesh], z, 'nodes')

            # MAJ BCs
            updatePointLists(z, z1s, pg_oids)

            # MAJ CHAMP CENTRE
            cont = Internal.getNodesFromName2(z, Internal.__FlowSolutionCenters__)
            fields = Internal.getNodesFromType1(cont, 'DataArray_t')

            for f in fields:
                pt = f[1].ravel('k')
                f[1] = numpy.empty( (ph_oids.size), numpy.float64)
                f[1][:] = pt[ph_oids[:]]

            zs.append(z)

        else: # stocke le chemin des zones a supprimer
            paths.append(Internal.getPath(t1, z))

        i += 3

    for p in paths: Internal._rmNodeFromPath(t1, p)
    paths = []

    # Search for BC family nodes
    fam_nodes = Internal.getNodesFromType(t2, 'Family_t')

    for node in fam_nodes:
        zs.append(node)

    iz=-1
    # New name dictionary and zone name updating
    newname2 = {}
    for z in z2s:
        iz += 1
        newname2[z[0]] = 'dom2_z_'+str(iz)
        z[0]           = newname2[z[0]]

    # Matches updating (for zone names)
    for z in z2s:
        # Update joins names
        joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
        for j in joins:
            oldname = Internal.getValue(j)
            if oldname in newname2:
                Internal.setValue(j,newname2[oldname])

    iz = -1
    for z in z2s:
        iz +=1

        mesh = res[i]
        pg_oids=res[i+1]
        ph_oids=res[i+2]

        # print(" ")
        # print("OP2 - pgoids: ", pg_oids)
        if mesh != []:

            # MAJ du maillage de la zone
            C.setFields([mesh], z, 'nodes')

            # MAJ BCs
            updatePointLists(z, z2s, pg_oids)

            # # MAJ CHAMP CENTRE
            cont = Internal.getNodesFromName2(z, Internal.__FlowSolutionCenters__)
            fields = Internal.getNodesFromType1(cont, 'DataArray_t')

            for f in fields:
                pt = f[1].ravel('k')
                f[1] = numpy.empty( (ph_oids.size), numpy.float64)
                f[1][:] = pt[ph_oids[:]]

            zs.append(z)

        else: # stocke le chemin des zones a supprimer
            paths.append(Internal.getPath(t2, z))

        i += 3

    for p in paths: Internal._rmNodeFromPath(t2, p)

    # Search for BC family nodes
    fam_nodes = Internal.getNodesFromType(t1, 'Family_t')

    for node in fam_nodes:
        zs.append(node)


    # Create new matches
    dict1_ptl = res[i]
    dict2_ptl = res[i+1]

    joins = Internal.getNodesFromType(zs, 'GridConnectivity_t')
    zones = Internal.getZones(zs)

    if joins == []:
        glob = 0

        for zoneR in dict1_ptl.keys():
            dict1_zR = dict1_ptl[zoneR]
            for zoneD in dict1_zR:
                faceListR = dict1_zR[zoneD]
                dict2_zD  = dict2_ptl[zoneD]
                faceListD = dict2_zD[zoneR]

                z1OppName = zones[zoneD][0]
                z2OppName = zones[zoneR][0]

                name1     = 'match%d_%d'%(zoneR+1,glob); glob += 1
                name2     = 'match%d_%d'%(zoneD+1,glob); glob += 1

                C._addBC2Zone(zones[zoneR],name1,'BCMatch',faceList=faceListR+1,\
                              zoneDonor=z1OppName, faceListDonor=faceListD+1)

                C._addBC2Zone(zones[zoneD],name2,'BCMatch',faceList=faceListD+1,\
                              zoneDonor=z2OppName, faceListDonor=faceListR+1)
    else:
        matchName = []
        for j in joins:
            matchName.append(j[0])

        glob = 0

        for zoneR in dict1_ptl.keys():
            dict1_zR = dict1_ptl[zoneR]
            for zoneD in dict1_zR:
                faceListR = dict1_zR[zoneD]
                dict2_zD  = dict2_ptl[zoneD]
                faceListD = dict2_zD[zoneR]

                z1OppName = zones[zoneD][0]
                z2OppName = zones[zoneR][0]

                name1     = 'match%d_%d'%(zoneR+1,glob); glob += 1
                name2     = 'match%d_%d'%(zoneD+1,glob); glob += 1

                # Assure unique match name
                while (name1 in matchName):
                    name1 = 'match%d_%d'%(zoneR+1,glob); glob += 1

                while (name2 in matchName):
                    name2 = 'match%d_%d'%(zoneD+1,glob); glob += 1

                C._addBC2Zone(zones[zoneR],name1,'BCMatch',faceList=faceListR+1,\
                              zoneDonor=z1OppName, faceListDonor=faceListD+1)

                C._addBC2Zone(zones[zoneD],name2,'BCMatch',faceList=faceListD+1,\
                              zoneDonor=z2OppName, faceListDonor=faceListR+1)


    if t1_is_tree == True and t2_is_tree == True:
        return Internal.merge([t1, t2])
    else:
        return zs


def booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanMinus(a1, a2, tol)
    Usage for volumes: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanMinus(s1, s2, tol, preserve_right, solid_right, agg_mode, improve_qual, itermax)
    return C.convertArrays2ZoneNode('minus', [s])

def diffSurf(a1, a2, tol=0., preserve_right=1, agg_mode=1, improve_qual=False, outward_surf=True, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between a volume mesh and a surface mesh."""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.diffSurf(s1, s2, tol, preserve_right, agg_mode,improve_qual, outward_surf, itermax)
    return C.convertArrays2ZoneNode('VmS', [s])

def booleanModifiedSolid(solid, a2, tol=0., preserve_solid=1, agg_mode=1, improve_qual=False, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the transformed input solid after solving the intersection of its skin with a2.
    Usage: booleanModifiedSolid(a1, a2, tol, preserve_right, solid_right)"""
    sld = C.getFields(Internal.__GridCoordinates__, solid)[0]
    operand = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanModifiedSolid(operand, sld, tol, preserve_solid, agg_mode, improve_qual, itermax)
    return C.convertArrays2ZoneNode('modified_solid', [s])

#==============================================================================
# XcellN
# IN: t: 3D NGON mesh
# IN: priorities: one-to-one priorities between components
# IN: output_type: 0: binary mask; 1: continuous mask (xcelln); 2: clipped surface.
# IN: rtol: relative tolerance
# OUT: returns a 3D NGON surface mesh with the xcelln field (if output_type=0/1, the clipped surface with solution if output_type=2)
#==============================================================================
def XcellN(t, priorities, output_type=0, rtol=0.05):
    """Computes the weight coefficients of visibility for overset grid configurations as a field called xcelln, for both surface and volume mesh of any kind.
    Usage: XcellN(t, priorities [, output_type, rtol])"""
    tp = Internal.copyRef(t)
    _XcellN(tp, priorities, output_type, rtol)
    return tp

#==============================================================================
# _XcellN (in-place version)
# IN: t: 3D NGON SURFACE mesh
# IN: priorities: one-to-one priorities between components
# IN: output_type: 0: binary mask; 1: continuous mask (xcelln); 2: clipped surface.
# OUT: returns a 3D NGON surface mesh with the xcelln field (if output_type=0/1, the clipped surface with solution if output_type=2)
#==============================================================================
def _XcellN(t, priorities, output_type=0, rtol=0.05):
    """Computes the weight coefficients of visibility for overset grid configurations as a field called xcelln, for any kind of surface mesh.
    Usage : _XcellN(t, priorities [, rtol])"""

    DIM = getTreeDim(t)
    #print ('DIM ? ' +str(DIM))
    if DIM != 2 and DIM != 3: # and DIM != 21: # 21 NUGA SURFACE
        raise ValueError('XcellN: the input file has an unsupported format or contain mixed 2D/3D zones.')

    _XcellN_(t, priorities, output_type, rtol)

#==============================================================================
# _XcellN_ (in-place version)
# IN: t: 3D NGON SURFACE mesh
# IN: priorities: one-to-one priorities between components
# IN: output_type: 0 : binary mask; 1 : continuous mask (xcelln) ; 2 : clipped surface.
# OUT: returns a 3D NGON surface mesh with the xcelln field (if output_type=0/1, the clipped surface with solution if output_type=2)
#==============================================================================
def _XcellN_(t, priorities, output_type=0, rtol=0.05):
    """Computes the weight coefficients of visibility for overset grid configurations as a field called xcelln, for any kind of surface mesh.
    Usage : _XcellNSurf(t, priorities [, rtol])"""

    DBG = False
    TIMER = False

    try: import Transform.PyTree as T
    except: raise ImportError("XcellN: requires Transform module.")
    try: import Post.PyTree as P
    except: raise ImportError("XcellN: requires Post module.")
    try: import Generator.PyTree as G
    except: raise ImportError("XcellN: requires Generator module.")

    WALLBCS = ['BCWall', 'BCWallInviscid','BCWallViscous', 'BCWallViscousIsothermal', 'BCSymmetryPlane']

    DIM = getTreeDim(t)

    bases = Internal.getBases(t)
    #(BCs,BCNames,BCTypes) = C.getBCs(t)

    if len(bases) == 1:
        raise ValueError('XcellN: Only one base in the file. Each component must be separated in a given Base. No check between zones of the same component.')

    min_compid = min(min(priorities, key=min))
    max_compid = max(max(priorities, key=max))
    if max_compid < 0:
        raise ValueError('XcellN: Negativle values passes as priorities. mus be component id (0-based).')
    if max_compid >= len(bases):
        raise ValueError('XcellN: Greatest component specified in priorities exceeds nb of components.')

    # 1 PREPARE INPUTS : NUGA NGON zones + oriented BAR boundaries with wall ids
    if TIMER:
        xcelln_time = time.time()
        xcelln_time2 = time.time()
    # 1.1 convert to NUGA NGON
    tNG = convertTree2NUGANGON(t, True) # keepBC
    #C.convertPyTree2File(tNG, 'tNG.cgns')

    # 1.2 reorient
    if DIM == 3: _reorient(tNG)

    if TIMER:
        print ('XCellN: Preparing Inputs ::: NGON convert & reorientation ::: CPU time : ',time.time()-xcelln_time2,'s')
        xcelln_time2 = time.time()

    #C.convertPyTree2File(tNG, 'tNGo.cgns')

    # 1.3 get boundaries and wall ids
    wall_ids = []
    boundaries = []
    wallfamilies = []

    for z in Internal.getZones(tNG):
        for btyp in WALLBCS:
            wallfamilies += C.getFamilyBCNamesOfType(t, btyp)

    wallfamilies = list(set(wallfamilies)) # make unique occurences
    #print(wallfamilies)

    basesNG = Internal.getBases(tNG)
    bid=-1
    for bNG in basesNG:
        bid += 1
        zj = concatenate(bNG, tol=1.e-10) # discard inner joins
        #if DBG: C.convertPyTree2File(zj, 'zj_'+str(bid)+'.cgns')
        b_bounds = externalFaces(zj) # keeping orientation
        if DBG: C.convertPyTree2File(b_bounds, 'bound_b_'+str(bid)+'.cgns')

        # empty set (e.g. sphere)
        nbf = nb_cells(b_bounds)
        if nbf == 0:
            boundaries.append(None)
            wall_ids.append(None)
            continue

        m_bounds = C.getFields(Internal.__GridCoordinates__, b_bounds)[0]
        boundaries.append(m_bounds)

        # extract wall ids from input tree t and map them to tNG
        bS = bases[bid]
        dims = Internal.getZoneDim(Internal.getZones(bS)[0])
        walls = C.extractBCOfType(bS, 'BCWall') # includes all wall types and families
        wallf = None
        if len(walls) != 0:
            if dims[0] == 'Structured': walls = C.convertArray2Hexa(walls)
            walls = T.join(walls)
            if DIM == 2:
                hook = C.createHook(b_bounds, function='elementCenters') # identifying edges
            elif DIM == 3:
                hook = C.createHook(b_bounds, function='faceCenters')
            wallf = C.identifyElements(hook, walls) # wallf are ids in boundaries
            wallf = wallf[wallf >= 1]
            if wallf.size > 0:
                wallf -= 1 # make it 0 based
            else: wallf = None

        wall_ids.append(wallf)

    if TIMER:
        print ('XCellN: Preparing Inputs ::: BC and Walls ::: CPU time : ',time.time()-xcelln_time2,'s')
        xcelln_time2 = time.time()

    # 1.4 get the zones in a single list with parent base id
    ngons = []
    basenum = []
    zwall_ids = [] # for double wall mgt
    base_id = -1
    for b in basesNG:
        base_id += 1
        zones = Internal.getZones(b)
        for z in zones:
            c = C.getFields(Internal.__GridCoordinates__, z)[0]
            ngons.append(c)
            basenum.append(base_id)
            zwallf = getBCPtListOfType(z, WALLBCS, wallfamilies)
            zwallf = numpy.asarray(zwallf, dtype=Internal.E_NpyInt)
            if zwallf.size > 0: zwallf -= 1 # make it 0 based
            zwall_ids.append(zwallf)

    if TIMER:
        print ('XCellN: Preparing Inputs: CPU time: ',time.time()-xcelln_time,'s')
        xcelln_time = time.time()

    # 2. COMPUTE THE COEFFS PER ZONE (PARALLEL OMP PER ZONE)
    #print(wall_ids)
    #print(boundaries)
    xcellns = XOR.XcellN(ngons, zwall_ids, basenum, boundaries, wall_ids, priorities, output_type, rtol)

    if TIMER:
        print ('XCellN: Computing: CPU time: ',time.time()-xcelln_time,'s')
        xcelln_time = time.time()

    # 3. OUTPUT
    if output_type == 2: # output as clipped NGON
        i = 0
        # delete all data that are no longer valid
        C._deleteFlowSolutions__(t, 'nodes') #no histo for nodes currently
        C._deleteGridConnectivity__(t)
        C._deleteZoneBC__(t)

        bases = Internal.getBases(t)
        paths = []
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                # grab solution ptrs
                cont = Internal.getNodesFromName2(z, Internal.__FlowSolutionCenters__)
                fields = Internal.getNodesFromType1(cont, 'DataArray_t')

                # updating mesh and solution at centers
                mesh = xcellns[i]
                pg_oids = xcellns[i+1]

                if mesh[1].size > 0:
                    C.setFields([mesh], z, 'nodes')

                    for f in fields:
                        pt = f[1].ravel('k')
                        f[1] = numpy.empty( (pg_oids.size), numpy.float64)
                        f[1][:] = pt[pg_oids[:]]
                else: # stocke le chemin des zones a supprimer
                    paths.append(Internal.getPath(t, z))
                i += 2
        for p in paths: Internal._rmNodeFromPath(t, p)
        G._close(t)

        if TIMER:
            print ('XCellN: Writing output: CPU time: ',time.time()-xcelln_time,'s')
            xcelln_time = time.time()

        return None

    #SET XCELLN FIELD tO ZONES
    # fixme : conversion here because tNG NUGA does not work with node2center
    #         use convertArray2NGon to ensure that t is converted whatever its type
    #         e.g convert2Hexa does not work if t is NGON
    tNG = C.convertArray2NGon(t)
    basesNG = Internal.getBases(tNG)

    zid = -1; bid = -1
    has_structured_bases = False
    for b in bases:
        bid +=1
        zones = Internal.getZones(b)
        dims = Internal.getZoneDim(zones[0])
        if dims[0] == 'Structured': # set field on tNG
            has_structured_bases = True
            zonesNG = Internal.getZones(basesNG[bid])
            for z in zonesNG:
                zid += 1
                #print('set field on tNG')
                C.setFields([xcellns[zid]], z, 'centers', False)
            mc = C.node2Center(basesNG[bid])
            hookC = C.createGlobalHook([mc], 'nodes')
            hookN = C.createGlobalHook([basesNG[bid]], 'nodes')

            C._identifySolutions(b, basesNG[bid], hookN, hookC, tol=1000.)
            C.freeHook(hookC)
            C.freeHook(hookN)
        else: # set field on t
            for z in zones:
                zid += 1
                C.setFields([xcellns[zid]], z, 'centers', False)

    if TIMER:
        print('XCellN: Writing output: CPU time: ',time.time()-xcelln_time,'s')

    return None

#==============================================================================
# P1ConservativeInterpolation
# IN : XXX
# IN : XXX
# OUT: XXX
#==============================================================================
def P1ConservativeInterpolation(tR, tD):

    tp = Internal.copyRef(tR)
    _P1ConservativeInterpolation(tp, tD)
    return tp

#==============================================================================
# _P1ConservativeInterpolation
# IN : XXX
# IN : XXX
# OUT: XXX
#==============================================================================
def _P1ConservativeInterpolation(tR, tD):

    zR = Internal.getZones(tR)
    zD = Internal.getZones(tD)

    for zr in zR:
        mR = C.getFields(Internal.__GridCoordinates__, zr)[0]

        for zd in zD:
            mD = C.getFields(Internal.__GridCoordinates__, zd)[0]
            fldD = C.getFields(Internal.__FlowSolutionCenters__, zd)[0]

            fldR = intersector.P1ConservativeInterpolation(mR, mD, fldD)
            C.setFields([fldR], zr, 'centers', False)

#==============================================================================
# superMesh
# IN: surfz: 3D NGON surface mesh to clip
# IN: sclip: 3D NGON surface mesh (clipper)
# IN: tol: tolerance (abolute if positive, relative otherwise)
# IN: proj_on_first: if True(False), each sclip(surfz) face is projected on surfz(sclip).
# OUT: returns the polyclipping of surfz by sclip
#==============================================================================
def superMesh(surfz, sclip, tol=-1.e-4, proj_on_first=True):
    """Polyclips surfz surface with sclip surface.
    Usage: superMesh(surfz, sclip, priorFirst, tol)"""
    m1 = C.getFields(Internal.__GridCoordinates__, surfz)[0]
    m2 = C.getFields(Internal.__GridCoordinates__, sclip)[0]

    res = XOR.superMesh(m1, m2, tol, proj_on_first)
    if res == []: return None # empty result

    mesh = res[0]
    anc  = res[1]

    xmatch = C.convertArrays2ZoneNode('xmatch', [mesh])
    nnuga = Internal.getNodeFromName(surfz, 'NUGA')
    if nnuga != None:
        nnug = Internal.addChild(xmatch, nnuga)
        nuga_fields = Internal.getChildren(nnug)
        for nf in nuga_fields:
            nf0 = nb_faces(surfz)
            fld_sz = len(nf[1])
            if nf0 == fld_sz : #valid face field stored => convert it
                nf[1] = updateNugaData(nf[1], anc)

    return xmatch


#==============================================================================
# superMesh2:
# IN: surfz: 3D NGON surface mesh to clip
# IN: sclip: 3D NGON surface mesh (clipper)
# IN: tol: tolerance (abolute if positive, relative otherwise)
# IN: proj_on_first: if True(False), each sclip(surfz) face is projected on surfz(sclip).
# OUT: returns face indices and face surface values
#==============================================================================
def superMesh2(surfz, sclip, tol=-1.e-4, proj_on_first=True):
    """Polyclips surfz surface with sclip surface.
    Usage: superMesh(surfz, sclip, priorFirst, tol)"""
    m1 = C.getFields(Internal.__GridCoordinates__, surfz)[0]
    m2 = C.getFields(Internal.__GridCoordinates__, sclip)[0]

    res = intersector.superMeshCompSurf(m1, m2, tol, proj_on_first)

    if res == []: return (None, None, None) # empty result

    mesh    = res[0]
    anc     = res[1]
    anc2    = res[2]
    surf    = res[3]
    isMatch = res[4]

    xmatch = C.convertArrays2ZoneNode('xmatch', [mesh])
    nnuga = Internal.getNodeFromName(surfz, 'NUGA')
    if nnuga != None:
        nnug = Internal.addChild(xmatch, nnuga)
        nuga_fields = Internal.getChildren(nnug)
        for nf in nuga_fields:
            nf0 = nb_faces(surfz)
            fld_sz = len(nf[1])
            if nf0 == fld_sz : #valid face field stored => convert it
                nf[1] = updateNugaData(nf[1], anc)

    # C.convertPyTree2File(xmatch,'xmatch.cgns')

    return (anc, anc2, surf, isMatch)
#==============================================================================
def extractBCMatchTNC(ancA,ancB,weight,fields, iminA, jminA, kminA,
                      imaxA, jmaxA, kmaxA):

    fld = intersector.computeTNCFields(ancA,ancB,weight,fields,
                                       (iminA, jminA, kminA, imaxA, jmaxA, kmaxA))

    return fld


#==============================================================================
# replaceFaces
# IN : z            :  3D NGON zone
# IN : face_soup    : surface mesh with oid indir referring to z faces' ids
# OUT: returns z with replaced faces
#==============================================================================
def replaceFaces(z, face_soup):
    """Replaces z faces with faces in face_soup having an indir for mapping.
    Usage: replaceFaces(z, face_soup)"""
    zp = Internal.copyRef(z)
    _replaceFaces(zp, face_soup)
    return zp

#==============================================================================
# _replaceFaces
# IN : z            :  3D NGON zone
# IN : face_soup    : surface mesh with oid indir referring to z faces' ids
# OUT: returns z with replaced faces
#==============================================================================
def _replaceFaces(z, face_soup):
    """Replaces z faces with faces in face_soup having an indir for mapping (in-place)
    Usage: _replaceFaces(z, face_soup)"""
    nnuga = Internal.getNodeFromName(face_soup, 'NUGA')
    if nnuga == None:
        print('_replaceFaces : error : NUGA node is missing')
        return

    vfoid = Internal.getNodeFromName(nnuga, 'vfoid')
    if vfoid == None:
        print('_replaceFaces : error : vfoid node is missing')
        return

    vfoid = vfoid[1]

    m    = C.getFields(Internal.__GridCoordinates__, z)[0]
    soup = C.getFields(Internal.__GridCoordinates__, face_soup)[0]

    res = intersector.replaceFaces(m, soup, vfoid)

    mesh = res[0]
    oids = res[1]

    # MAJ du maillage de la zone
    C.setFields([mesh], z, 'nodes')

    # MAJ POINT LISTS #
    updateBCPointLists1(z, oids)


#==============================================================================
# triangulateExteriorFaces
# IN : t            :  3D NGON mesh
# IN : in_or_out    : 0 means "ONLY INTERNALS", 1 means "ONLY EXTERNALS", any other value means "BOTH"
# IN : improve_qual :  to improve mesh quality
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def triangulateExteriorFaces(t, in_or_out=2, improve_qual=0):
    """Triangulates exterior polygons of a volume mesh.
    Usage: triangulateExteriorFaces(t)"""
    return C.TZA1(t, 'nodes', 'nodes', True, XOR.triangulateExteriorFaces, in_or_out, improve_qual)

#==============================================================================
# _triangulateExteriorFaces
# IN : t            :  3D NGON mesh
# IN : in_or_out    : 0 means "ONLY INTERNALS", 1 means "ONLY EXTERNALS", any other value means "BOTH"
# IN : improve_qual :  to improve mesh quality
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def _triangulateExteriorFaces(t, in_or_out=2, improve_qual=0):
    """Triangulates exterior polygons of a volume mesh (in-place).
    Usage: triangulateExteriorFaces(t)"""
    return C._TZA1(t, 'nodes', 'nodes', True, XOR.triangulateExteriorFaces, in_or_out, improve_qual)


#==============================================================================
# triangulateSpecifiedFaces
# IN : t            : 3D NGON mesh
# IN : pgs          : list of polygons
# IN : improve_qual : to improve mesh quality
# OUT: returns a 3D NGON Mesh
#==============================================================================
def triangulateSpecifiedFaces(t, pgs, improve_qual=1):
    """Triangulates polygons specified by pgs.
    Usage: triangulateSpecifiedFaces(t, pgs, improve_qual)"""
    tp = Internal.copyRef(t)
    _triangulateSpecifiedFaces(tp,pgs, improve_qual)
    return tp

#==============================================================================
# _triangulateSpecifiedFaces
# IN: a: 3D NGON mesh
# IN: pgs : list of polygons
# OUT: returns a 3D NGON Mesh
#==============================================================================
def _triangulateSpecifiedFaces(t, pgs, improve_qual=1):
    """Triangulates polygons specified by pgs (in-place).
    Usage: triangulateSpecifiedFaces(t, pgs, improve_qual)"""
    zones = Internal.getZones(t)
    if (len(pgs) != len(zones)) :
        print('triangulateSpecifiedFaces: input error: nb of polygons packs differ from nb of zones.')
        return None

    i = 0
    for z in zones:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue
        m = Converter.convertArray2NGon(m)
        m = XOR.triangulateSpecifiedFaces(m, pgs[i], improve_qual)
        mesh = m[0]
        pg_oids= m[1]

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        # MAJ POINT LISTS #
        updatePointLists(z, zones, pg_oids)
        i += 1

#==============================================================================
# triangulateNonBasicFaces
# IN: mesh: 3D NGON mesh
# IN : quality improvement flag
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def triangulateNFaces(t, improve_qual=1, min_nvertices=5, discard_joins=True):
    """Triangulates nob basic polygons of a volume mesh.
    Usage: triangulateNonBasicFaces(t)"""
    tp = Internal.copyRef(t)
    _triangulateNFaces(tp,improve_qual, min_nvertices, discard_joins)
    return tp

def _triangulateNFaces(t, improve_qual=1, min_nvertices=5, discard_joins=True):

    zones = Internal.getZones(t)

    for z in zones:
        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue

        #coords = Converter.convertArray2NGon(coords)

        ptLists=[]
        if discard_joins:
            joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
            for j in joins:
                ptl = Internal.getNodeFromName1(j, 'PointList')
                ptLists.append(ptl[1])

        if ptLists != []:
            #print ptLists
            ptLists = numpy.concatenate(ptLists) # create a single list
            ptLists = ptLists -1 # 0-based
            ptLists = numpy.concatenate(ptLists) # create a single list
            #print ptLists
        else:
            ptLists=None

        res = XOR.triangulateNFaces(coords, improve_qual, min_nvertices, ptLists)

        mesh = res[0]
        pg_oids=res[1]

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        # MAJ POINT LISTS #
        updatePointLists(z, zones, pg_oids)


#==============================================================================
# triangulateBC
# IN: t: 3D NGON mesh
# IN : btype : boundary type to mesh
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def triangulateBC(t, bctype):

    tp = Internal.copyRef(t)
    _triangulateBC(tp,bctype)
    return tp
#==============================================================================
# _triangulateBC
# IN: t: 3D NGON mesh
# IN : btype : boundary type to mesh
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def _triangulateBC(t, bctype):

    zones = Internal.getZones(t)

    for z in zones:

        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue

        coords = Converter.convertArray2NGon(coords)

        bnds = Internal.getNodesFromType(z, 'BC_t')

        bcpgs = []
        for bb in bnds :
            if (Internal.isValue(bb, bctype) == False) : continue
            bcpgs.append(bb[2][1][1][0]) # POINTLIST NUMPY

        if bcpgs == []: continue
        bcpgs = numpy.concatenate(bcpgs) # create a single list
        bcpgs = bcpgs -1

        res = XOR.triangulateSpecifiedFaces(coords, bcpgs)

        mesh = res[0]
        pg_oids=res[1]

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        # MAJ POINT LISTS #
        updatePointLists(z, zones, pg_oids)

    return t


#==============================================================================
# convexifyFaces
# IN: mesh: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all the faces made convex
#==============================================================================
def convexifyFaces(t, convexity_TOL=1.e-8):
    """Convexifies any non-convex polygon in a mesh.
    Usage: convexifyFaces(t)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.convexifyFaces(m, convexity_TOL)
    return C.convertArrays2ZoneNode('allPGconvex', [m])

#==============================================================================
# externalFaces : Returns erternal faces for CASSIOPEE NGON types and NUGA NGON
#==============================================================================
def externalFaces(t, discard_joins=False, geo_dim=-1):
    """Returns erternal faces for CASSIOPEE NGON types and NUGA NGON.
    Usage: externalFaces(t)"""
    zs = Internal.getZones(t)
    efs = []

    i=-1
    for z in zs:
        i+=1
        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue
        ids=None
        if discard_joins == True:
            joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
            if joins != []:
                ids=[]
                for j in joins :
                    ptl = Internal.getNodeFromName1(j, 'PointList')
                    #print(ptl[1][0])
                    ids.append(ptl[1][0])

                if (ids != []):
                    #print(ids)
                    ids = numpy.concatenate(ids) # create a single list
                    ids = ids -1 # 0-based
                    #ids = numpy.concatenate(ids) # create a single list
                else:
                    ids=None
        res = XOR.externalFaces(coords, ids, geo_dim)
        #print(res)
        ef = res[0]
        vfoid=res[1] #ids in volume mesh

        zid = getProperty(z, 'zid')
        #print('zid : '+str(zid))
        zid_defined = True
        if zid == -1 :
            zid=i
            zid_defined=False

        zonnod = C.convertArrays2ZoneNode('ef_z'+str(zid), [ef])

        if zid_defined == True : # add it to the zone
            l = Internal.newIntegralData(name='zid', parent=zonnod) #ceate a zid for each block
            l[1] = zid

        # node : NUGA -> vfoid
        Internal.createChild(zonnod, 'NUGA', 'UserDefinedData_t', value=None, children=[])
        nnuga = Internal.getNodeFromName(zonnod, 'NUGA')
        Internal.newDataArray('vfoid', value=vfoid, parent=nnuga)

        efs.append(zonnod)
    return efs

#==============================================================================
# reorient : reorients outward the external polygons of a mesh
#==============================================================================
def reorient(t, dir=1):
    """Reorients outward the external polygons of a mesh.
    Usage: reorient(t)"""
    return C.TZA1(t, 'nodes', 'nodes', True, XOR.reorient, dir)

def _reorient(t, dir=1):
    return C._TZA1(t, 'nodes', 'nodes', True, XOR.reorient, dir)

#==============================================================================
# triangulateBC
# IN: t: 3D NGON mesh
# IN : btype : boundary type to reorient
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def reorientBC(t, bctype, dir):

    tp = Internal.copyRef(t)
    _reorientBC(tp,bctype, dir)
    return tp
#==============================================================================
# _reorientBC
# IN: t: 3D NGON mesh
# IN : btype : boundary type to reorient
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def _reorientBC(t, bctype, dir):

    zones = Internal.getZones(t)

    for z in zones:

        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue

        coords = Converter.convertArray2NGon(coords)

        bnds = Internal.getNodesFromType(z, 'BC_t')

        bcpgs = []
        for bb in bnds :
            if (Internal.isValue(bb, bctype) == False) : continue
            bcpgs.append(bb[2][1][1][0]) # POINTLIST NUMPY

        if bcpgs == []: continue
        bcpgs = numpy.concatenate(bcpgs) # create a single list
        bcpgs = bcpgs -1

        res = XOR.reorientSpecifiedFaces(coords, bcpgs, dir)

        mesh = res[0]
        pg_oids=res[1]

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

    return t

#==============================================================================
# reorientSpecifiedFaces
# IN: a: 3D NGON mesh
# IN: pgs : list of polygons
# OUT: returns a 3D NGON Mesh with consistent orientation for specified polygons
#==============================================================================
def reorientSpecifiedFaces(t, pgs, dir):

    tp = Internal.copyRef(t)
    _reorientSpecifiedFaces(tp,pgs, dir)
    return tp

def _reorientSpecifiedFaces(t, pgs, dir):

    zones = Internal.getZones(t)
    if (len(pgs) != len(zones)) :
        print('reorientSpecifiedFaces: input error: nb of polygons packs differ from nb of zones.')
        return None

    i=0
    for z in zones:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        m = XOR.reorientSpecifiedFaces(m, pgs[i], dir)
        C.setFields([m], z, 'nodes') # replace the mesh in the zone
        i = i+1

#==============================================================================
# reorientSurf
# IN: a: 3D NGON mesh
# IN: pgs : list of polygons
# OUT: returns a 3D NGON Mesh with consistent orientation for specified polygons
#==============================================================================
def reorientSurf(t, dir):

    tp = Internal.copyRef(t)
    _reorientSurf(tp, dir)
    return tp

def _reorientSurf(t, dir):

    zones = Internal.getZones(t)

    for z in zones:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        m = XOR.reorientSurf(m, dir)
        C.setFields([m], z, 'nodes') # replace the mesh in the zone

#==============================================================================
# prepareCellsSplit
# IN : t            : 3D NGON mesh
# IN : PH_set       : PH to process. 0 for concave cells or 1 for non-centroid-star_shaped cells
# IN : split_policy : 0 : convexify concave pgs. 1 : starify concave pgs. 2 : starify any pgs at concave-chains ends.
# OUT: returns a 3D NGON Mesh with some face polygon splits :
#      split (convexify, starify) some targeted polygons on targeted cells
#      (typically bad polyhedra -concaves, non-centroid-star-shaped-)
#      to prepare the split of those bad cells.
#==============================================================================
def prepareCellsSplit(t, PH_set=1, split_policy=0, PH_conc_threshold=1./3., PH_cvx_threshold=0.05, PG_cvx_threshold=1.e-8):
    """Splits some prescribed polygons following a prescribed splitting policy.
    Usage: prepareCellsSplit(t, PH_set, split_policy, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.prepareCellsSplit(m, PH_set, split_policy, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)
    return C.convertArrays2ZoneNode('cvxSplitReady', [m])

#==============================================================================
# splitNonStarCells
# IN : t                 : 3D NGON mesh
# IN : PH_conc_threshold : concavity dihedral angle threshold for cells
# IN : PH_cvx_threshold  : convexity dihedral angle threshold for cells
# IN : PG_cvx_threshold  : convexity angle threshold for polygons
# OUT: returns a 3D NGON Mesh with none (or at least less) non-centroid-star_shaped cells.
#==============================================================================
def splitNonStarCells(t, PH_conc_threshold=1./3., PH_cvx_threshold=0.05, PG_cvx_threshold=1.e-8):
    """Splits some non-centroid-star_shaped cells.
    Usage: splitNonStarCells(t, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.splitNonStarCells(m, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)
    return C.convertArrays2ZoneNode('splitNonStarCells', [m])

#==============================================================================
# simplifyCells : agglomerate superfluous polygons that overdefine cells
# IN: mesh: 3D NGON mesh
# IN: angular_threshold : should be as small as possible to avoid introducing degeneracies
# IN: discard_joins : when set to True, faces at joins are not modified
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def _simplifyCells(t, treat_externals, angular_threshold=1.e-12, discard_joins=True):
    """Simplifies over-defined polyhedral cells (agglomerate some elligible polygons).
    Usage: simplifyCells(t, treat_externals, angular_threshold, discard_joins)"""

    zones = Internal.getZones(t)
    ozones = []

    for z in zones:

        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        ids=None

        if treat_externals == 1 and discard_joins == True:
            joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
            if joins != []:
                ids=[]
                for j in joins :
                    ptl = Internal.getNodeFromName1(j, 'PointList')
                    ids.append(ptl[1])

                if (ids != []):
                    #print ids
                    ids = numpy.concatenate(ids) # create a single list
                    ids = ids -1 # 0-based
                    ids = numpy.concatenate(ids) # create a single list
                    #print ids
                else:
                    ids=None

        m = XOR.simplifyCells(m, treat_externals, angular_threshold, discarded_ids=ids)
        C.setFields([m], z, 'nodes') # replace the mesh in the zone


def simplifyCells(t, treat_externals, angular_threshold=1.e-12, discard_joins=True):
    """Simplifies over-defined polyhedral cells (agglomerate some elligible polygons).
    Usage: simplifyCells(t, treat_externals, angular_threshold, discard_joins)"""
    tp = Internal.copyRef(t)
    _simplifyCells(tp, treat_externals, angular_threshold, discard_joins)
    return tp
#==============================================================================
# simplifyFaces : remove superfluous nodes
# IN: mesh: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def simplifyFaces(t):
    """simplify over-defined polygons
    Usage: simplifyFaces(t)"""

    zones = Internal.getZones(t)
    ozones = []

    for z in zones:

        m = C.getFields(Internal.__GridCoordinates__, t)[0]

        m = XOR.simplifyFaces(m)

        ozones.append(C.convertArrays2ZoneNode(z[0], [m]))

    return ozones

#==============================================================================
# simplifySurf : agglomerate superfluous polygons that overdefine the surface
# IN: mesh: 3D NGON mesh
# IN: angular_threshold : should be as small as possible to avoid introducing degeneracies
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def simplifySurf(t, angular_threshold=1.e-12):
    """Simplifies over-defined surfaces (agglomerate some elligible polygons).
    Usage: simplifySurf(t, angular_threshold)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.simplifySurf(m, angular_threshold)
    return C.convertArrays2ZoneNode('simplifiedSurf', [m])

#==============================================================================
# agglomerateSmallCells : agglomerate prescribed cells
# IN: t: 3D NGON mesh
# IN: vmin : volume threshold
# IN: vratio : aspect ratio threshold
# IN: angular_threshold : for simplying cells by agglomerating adjacent polygons
# IN: method = 0 (XXX)
# OUT: returns a 3D NGON Mesh with less cells and with a smoother aspect ratio
#==============================================================================
def agglomerateSmallCells(t, vmin=0., vratio=0.01, angular_threshold=1.e-12, method=0):
    """Agglomerates prescribed cells.
    Usage: agglomerateSmallCells(t, vmin, vratio)"""

    nb_phs0  = nb_cells(t)

    m = C.getFields(Internal.__GridCoordinates__, t)[0]

    res = XOR.agglomerateSmallCells(m, vmin, vratio, angular_threshold, method)
    #print("NB ZONES %d"%(len(res)))

    z = C.convertArrays2ZoneNode('agglomeratedCells', [res[0]])

    nb_phs1  = nb_cells(z)
    if nb_phs1 < nb_phs0:
        nbc = int(nb_phs0 - nb_phs1)
        print('agglomerateSmallCells : Nb of small cells agglomerated : ' + str(nbc))
    else:
        print('agglomerateSmallCells : None agglomeration.')

    debug = False

    if (debug == False) : return z

    zones = []

    nb_zones = len(res)-1
    if (nb_zones == 0) : return z

    for i in range(nb_zones):
        zones.append(C.convertArrays2ZoneNode('agg', [res[i+1]]))

    C.convertPyTree2File(zones, 'agglo.cgns')

    return z

#==============================================================================
# shellAgglomerateSmallCells : eradicate small cells by agglomerating all surrounding cells
# IN: t: 3D NGON mesh
# IN: vmin : volume threshold
# IN: vratio : aspect ratio threshold
# OUT: returns a 3D NGON Mesh with less cells and with a smoother aspect ratio
#==============================================================================
def shellAgglomerateSmallCells(t, vmin=0., vratio=1000.):
    """Agglomerates prescribed cells.
    Usage: shellAgglomerateSmallCells(t, vmin, vratio)"""

    nb_phs0  = nb_cells(t)

    m = C.getFields(Internal.__GridCoordinates__, t)[0]

    res = XOR.shellAgglomerateSmallCells(m, vmin, vratio)

    #print("NB ZONES %d"%(len(res)))

    z = C.convertArrays2ZoneNode('agglomeratedCells', [res[0]])

    nb_phs1  = nb_cells(z)
    if nb_phs1 < nb_phs0:
        nbc = int(nb_phs0 - nb_phs1)
        print('shellAgglomerateSmallCells : Nb of small cells agglomerated : ' + str(nbc))
    else:
        print('shellAgglomerateSmallCells : None agglomeration.')

    return z

#==============================================================================
# agglomerateNonStarCells : Agglomerates non-centroid-star-shaped cells
# IN: t: 3D NGON mesh
# IN: vmin : volume threshold
# IN: vratio : aspect ratio threshold
# OUT: returns a 3D NGON Mesh with less cells and with a smoother aspect ratio
#==============================================================================
def agglomerateNonStarCells(t, angular_threshold=1.e-12):
    """Agglomerates non-centroid-star-shaped cells.
    Usage: agglomerateNonStarCells(t)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]

    res = XOR.agglomerateNonStarCells(m, angular_threshold)
    #print("NB ZONES %d"%(len(res)))

    z = C.convertArrays2ZoneNode('agglomeratedCells', [res[0]])

    debug = False

    if (debug == False) : return z

    zones = []

    nb_zones = len(res)-1
    if (nb_zones == 0) : return z

    for i in range(nb_zones):
        zones.append(C.convertArrays2ZoneNode('agg', [res[i+1]]))

    C.convertPyTree2File(zones, 'agglo.cgns')

    return z

#==============================================================================
# agglomerateCellsWithSpecifiedFaces : Agglomerates cells sharing specified polygons
# IN: a: 3D NGON mesh
# IN: pgs : list of polygons
# OUT: returns a 3D NGON Mesh
#==============================================================================
def agglomerateCellsWithSpecifiedFaces(t, pgs, simplify=1, amax=1.e-12, treat_externals=1): # 0 : dno not simplify, 1 : simplify only internals, 2 : simlplify evrywhere

    tp = Internal.copyRef(t)
    _agglomerateCellsWithSpecifiedFaces(tp,pgs, simplify, amax, treat_externals)
    return tp

def _agglomerateCellsWithSpecifiedFaces(t, pgs, simplify=1, amax=1.e-12, treat_externals=1):

    zones = Internal.getZones(t)
    if len(pgs) != len(zones):
        print('agglomerateCellsWithSpecifiedFaces: input error: nb of polygons packs differ from nb of zones : %s versus %s.'%(len(pgs), len(zones)))
        return None

    if simplify < 0 : simplify = 0
    if simplify > 1 : simplify = 1

    if amax < 0. : amax = 3.15 # value greater than Pi, to have no restrictions

    i=-1
    for z in zones:

        i+=1
        jids = []

        if simplify == 1 and treat_externals == 1: # keep track of join PGs to freeze them when simplifying
            joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
            if joins != []:
                for j in joins :
                    ptl = Internal.getNodeFromName1(j, 'PointList')
                    jids.append(ptl[1][0])

                if (jids != []):
                    #print(ids)
                    jids = numpy.concatenate(jids) # create a single list
                    jids = jids -1 # 0-based
                    #print(jids)
                else:
                    jids=None

        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res = XOR.agglomerateCellsWithSpecifiedFaces(m, pgs[i]) # FIXME : need oids to use discarded_ids

        m    = res[0]
        nids = res[1] # 0-based

        # update join ids
        if jids != [] and len(nids) != 0:
            for k in range(len(jids)):
                oj = jids[k]
                nj = nids[oj]
                jids[k]=nj
            jids = jids[jids > -1]
        else:
            jids = None

        if simplify == 1:
            m = XOR.simplifyCells(m, treat_externals, angular_threshold=amax, discarded_ids=jids)

        C.setFields([m], z, 'nodes') # replace the mesh in the zone

#==============================================================================
# agglomerateUncomputableCells : XXX
#==============================================================================
# def agglomerateUncomputableCells(mesh):
#     m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
#     m = XOR.agglomerateUncomputableCells(m)
#     return C.convertArrays2ZoneNode('agglomeratedCells', [m])

#==============================================================================
# collapseUncomputableFaces : XXX
#==============================================================================
def collapseUncomputableFaces(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.collapseUncomputableFaces(m)
    return C.convertArrays2ZoneNode('agglomeratedCells', [m])

#==============================================================================
# collapseSmallCells : XXX
#==============================================================================
def collapseSmallCells(t, vmin=0., grmin=-1.):
    """XXX"""
    tp = Internal.copyRef(t)
    _collapseSmallCells(tp, vmin, grmin)
    return tp

def _collapseSmallCells(t, vmin=0., grmin=-1.):
    zones = Internal.getZones(t)
    for z in zones:
        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue
        collapsed = XOR.collapseSmallCells(coords, vmin, grmin)
        C.setFields([collapsed], z, 'nodes')

#==============================================================================
# removeNonManifoldExternalCells : removes any outer cell that has a non manifold edge
#==============================================================================
def removeNonManifoldExternalCells(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.removeNonManifoldExternalCells(m)
    return C.convertArrays2ZoneNode('manifoldOuter', [m])

#==============================================================================
# immerseNodes : Regularize the mesh to be cut (by boolean operations) by immersing its nodes lying on cutter surface
# IN: t : 3D NGON mesh
# IN: s : unstructured surface mesh
# IN: TOL : tolerance (negative value means relative)
# OUT: returns a 3D NGON Mesh with moved vertices
#==============================================================================
def immerseNodes(t, s, TOL):
    """Regularize the mesh to be cut (by boolean operations) by immersing its nodes lying on the cutter surface.
    Usage: immerseNodes(t,s, TOL)"""
    tp = Internal.copyRef(t)
    _immerseNodes(tp, s, TOL)
    return tp

def _immerseNodes(t, s, TOL):
    """Regularize the mesh to be cut (by boolean operations) by immersing its nodes lying on the cutter surface.
    Usage: immerseNodes(t,s, TOL)"""

    zones = Internal.getZones(t)

    DIM = getTreeDim(t)
    if DIM != 3:
        print ('immerseNodes : t is not a 3D mesh')
        return

    szones = Internal.getZones(s)
    sdims = Internal.getZoneDim(szones[0])
    #print(sdims[0])
    # go unstructured
    if sdims[0] == 'Structured':
        s = C.convertArray2Hexa(s)
    # single surface mesh
    s = T.join(s)
    # convert to NUGA surface
    if sdims[3] == 'NGON':
        _convertNGON2DToNGON3D(s)
    else:
        _convertBasic2NGONFaces(s)

    # now we have a NUGA unstructured surface
    surf = C.getFields(Internal.__GridCoordinates__, s)[0]

    for z in zones:

        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue

        dims = Internal.getZoneDim(z)
        if dims[3] != 'NGON' :
            coords = Converter.convertArray2NGon(coords)

        moved = intersector.immerseNodes(coords, surf, TOL)

        # MAJ du maillage de la zone
        C.setFields([moved], z, 'nodes')

    return t

#==============================================================================
# closeCells : Closes any polyhedral cell in a mesh (processes hanging nodes on edges)
# IN: t: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all cells closed
#==============================================================================
def closeCells(t):
    """Closes any polyhedral cell in a mesh (processes hanging nodes on edges).
    Usage: closeCells(t)"""
    tp = Internal.copyRef(t)
    _closeCells(tp)
    return tp

def _closeCells(t):
    """Closes any polyhedral cell in a mesh (processes hanging nodes on edges).
    Usage: _closeCells(t)"""

    zs  = Internal.getZones(t)

    meshes = []
    zids   = []

    for z in zs:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m ==[] : continue

        meshes.append(m)

        zid=getProperty(z, 'zid')
        zids.append(zid)

    zidDict={}
    zs = Internal.getZones(t)
    for z in zs:
        zidDict[z[0]]=getProperty(z, 'zid')

    #
    zid_to_rid_to_list_owned = getJoinsPtLists(t, zidDict)

    rid_to_zones = getRidToZones(t, zidDict)

    # SEQ or open-MP processing across joins
    meshes = intersector.closeCells(meshes, zids, zid_to_rid_to_list_owned, rid_to_zones)

    # associate zid and closed mesh
    zid_to_m = {}
    for i in range(len(zids)) :
        zid = zids[i]
        m = meshes[i]
        zid_to_m[zid]=m

    for z in zs :
        zid = getProperty(z, 'zid')
        m = zid_to_m[zid]
        # MAJ du maillage de la zone
        C.setFields([m], z, 'nodes')

    return t

#==============================================================================
# getJoinsPtList : XXX
# IN: t : 3D NGON PyTree
# IN: XXX
# OUT: XXX
#==============================================================================
def getJoinsPtLists(t, zidDict):

    zone_to_rid_to_list_owned = {}

    zones = Internal.getZones(t)

    for z in zones:
        zid = getProperty(z, 'zid')

        rid_to_list_owned = {}

        raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for rac in raccords:

            rid = getProperty(rac, 'rid')

            donnorName = "".join(Internal.getValue(rac))
            ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]

            jzid = zidDict[donnorName]
            if jzid==zid: # auto-match case
                if rid in rid_to_list_owned:
                    rid_to_list_owned[rid] = numpy.concatenate((rid_to_list_owned[rid], ptList))
                else:
                    rid_to_list_owned[rid] = ptList
            else:         # other cases
                rid_to_list_owned[rid] = ptList

        zone_to_rid_to_list_owned[zid] = rid_to_list_owned

    return zone_to_rid_to_list_owned


#==============================================================================
# getRidToZones : Computes a dictionary rid <-> (zid1,zid2)
# IN: t : 3D NGON PyTree
# IN: zidDict : dictionary zname <-> zid
# OUT: Returns the dictionary rid <-> (zid1,zid2)
#==============================================================================
def getRidToZones(t, zidDict):
    """ Function returning ridDict (zones concerned 
    by rac rid, given tree t as well as zidDict.
    """

    ridDict = {}

    zones = Internal.getZones(t)
    for z in zones:
        zid = getProperty(z, 'zid')

        raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t')
        for rac in raccords:
            rid = getProperty(rac, 'rid')

            z1 = getProperty(z, 'zid')

            donnorName = "".join(Internal.getValue(rac))
            z2 = zidDict[donnorName]

            ridDict[rid] = (z1,z2)

    return ridDict


#==============================================================================
# adaptCells : Adapts an unstructured mesh a with respect to a sensor
# IN: t: 3D NGON mesh
# IN: sensdata: sensor data (a bunch of vertices or a mesh for a geom sensor, a mesh for a xsensor, punctual values for a nodal or cell sensor)
# IN: sensor_type: geom_sensor (0) , xsensor (1), nodal_sensor (2), cell_sensor(3), xsensor(4)
# IN smoothing_type: First-neighborhood (0) Shell-neighborhood(1)
# IN itermax: max number of level in the hierarchy
# IN: subdiv_type: isotropic currently
# IN: sensor_metric_policy (specific for xsensor): which reference cell size (edge length) to use ? min (0), mean (1), max(2) or min_or_max(3)
# IN: hmesh: hierarchical mesh hook
# IN: sensor: sensor hook
# IN: conformize: dump a conformal polyhedral mesh
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCells(t, sensdata=None, sensor_type=0, smoothing_type=0, itermax=-1, subdiv_type=0, sensor_metric_policy=0, hmesh=None, sensor=None, conformize=1):
    """Adapts an unstructured mesh a with respect to a sensor.
    Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, sensor_metric_policy=0, hmesh=None, sensor=None)"""
    tp = Internal.copyRef(t)
    _adaptCells(tp, sensdata, sensor_type, smoothing_type, itermax, subdiv_type, sensor_metric_policy, hmesh, sensor, conformize)
    return tp

def _adaptCells(t, sensdata=None, sensor_type=0, smoothing_type=0, itermax=-1, subdiv_type=0, sensor_metric_policy=0, hmesh=None, sensor=None, conformize=1):
    """Adapts an unstructured mesh a with respect to a sensor.
    Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, sensor_metric_policy=0, hmesh=None, sensor=None)"""

    if sensor_type == 1: sensor_type=4 # initial xsensor does not exist anymore

    if sensor_type != 5 and sensdata is None and sensor is None:
        print('INPUT ERROR : no source data to initialize a sensor')
        return

    if hmesh is None and sensor is not None:
        print ('INPUT ERROR : you must also give as an argument the hmesh targeted by the sensor')
        return

    zs = Internal.getZones(t)
    NBZ = len(zs)

    owesHmesh=0

    if hmesh is None :
        #print("create hm : ", NBZ)
        hmesh = createHMesh(t, subdiv_type)
        if hmesh == [None] : # no basic elts in t (assumed because hmesh creation as done [None]
            return
        owesHmesh=1

    owesSensor=0
    if sensor is None :
        #print("create sensor")
        sensor = createSensor(hmesh, sensor_type, smoothing_type, itermax, sensor_metric_policy)
        owesSensor=1

    err=0

    if sensdata is None and sensor_type == 5 : #extract metric fields
        sensdata = [] # going to fill it
        for z in zs :
            mxx = C.getField("mxx", z)
            if mxx == [] : err = 1; break
            mxy = C.getField("mxy", z)
            if mxy == [] : err = 1; break
            mxz = C.getField("mxz", z)
            if mxz == [] : err = 1; break
            myy = C.getField("myy", z)
            if myy == [] : err = 1; break
            myz = C.getField("myz", z)
            if myz == [] : err = 1; break
            mzz = C.getField("mzz", z)
            if mzz == [] : err = 1; break
            sensdata.append([mxx[0][1], mxy[0][1], mxz[0][1], myy[0][1], myz[0][1], mzz[0][1]])

        if err == 1 :
            print('INPUT ERROR : metric fields are missing or incomplete in the inout tree')
            return

    #print(sensdata)

    if sensdata is not None:
        #print("assignData2Sensor")
        if sensor_type == 1 or sensor_type == 4:
            sensdata = C.convertArray2NGon(sensdata)
        err = assignData2Sensor(sensor, sensdata)
        if err == 1:
            print('INPUT ERROR : sensor data list must be sized as nb of sensors')
            return

    zidDict={}
    zs = Internal.getZones(t)
    for z in zs:
        zidDict[z[0]]=getProperty(z, 'zid')

    #
    zone_to_rid_to_list_owned = getJoinsPtLists(t, zidDict)

    rid_to_zones = getRidToZones(t, zidDict)

    #print('adaptCells..')
    intersector.adaptCells(hmesh, sensor, zone_to_rid_to_list_owned, rid_to_zones) #fxme agglo

    if owesHmesh == 1: #and owesSensor == 1 :
        #print("_conformizeHMesh")
        _conformizeHMesh(t, hmesh, conformize)

    if owesHmesh == 1:
        #   #print('delete owned hmesh')
        deleteHMesh(hmesh)
    if owesSensor == 1:
        #print('delete owned sensor')
        deleteSensor(sensor)

#==============================================================================
# adaptCellsNodal (deprecated): Adapts a polyhedral mesh a1 with repsect to the nodal subdivision values.
# IN: t: 3D NGON mesh
# IN: nodal_vals: nb of subdivision required expressed at mesh nodes
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCellsNodal(t, sensdata=None, smoothing_type=0, subdiv_type=0, hmesh=None, sensor=None):
    """Adapts an unstructured mesh a with respect to a sensor (DEPRECATED : use adaptCells with a sensor of type 2 instead.)
    Usage: adaptCellsNodal(t, sensdata=None, smoothing_type = 0, subdiv_type=0, hmesh=None, sensor=None)"""
    tp = Internal.copyRef(t)
    _adaptCells(tp, sensdata, 2, smoothing_type, -1, 0, subdiv_type, hmesh, sensor)
    return tp

#==============================================================================
# setZonesAndJoinsUId : Assigns a zid and  the dictionary zid <-> rank
# IN: zidDict : dictionary zname <-> zid
# IN: procDict : dictionary zname <-> rank
# OUT: Returns the dictionary zid <-> rank
#==============================================================================
def setZonesAndJoinsUId(t):
    """Adapts an unstructured mesh a with respect to a sensor (DEPRECATED : use adaptCells with a sensor of type 2 instead.)
    Usage: adaptCellsNodal(t, sensdata=None, smoothing_type = 0, subdiv_type=0, hmesh=None, sensor=None)"""
    tp = Internal.copyRef(t)
    _setZonesAndJoinsUId(tp)
    return tp

def _setZonesAndJoinsUId(t):

    zs = Internal.getZones(t)
    iz = -1
    ir = -1

    # ---------------------------------------------------------------------
    # dict_VD = {z1 : {z2 : {0 : {[idx_min_1, idx_min_2], [ir_1, ir_2]}, 1 : {[idx_min_1, idx_min_2], [ir_3, ir_4]}}}
    #                 {z4 : {0 : {[idx_min], [ir_5]}, 1 : {[idx_min], [ir_6]}, 1 : {[idx_min], [ir_7]}, 1 : {[idx_min], [ir_8]}}}
    #           {z2 : {z3 : {0 : {[idx_min, ...], [ir, ...]}}}
    #           {z3 : {z4 : {0 : {[idx_min, ...], [ir, ...]}}}}
    # ---------------------------------------------------------------------

    dict_VD = {}
    zidDict = {}

    for z in zs: # loop on blocks
        iz += 1

        zidDict[z[0]] = iz #build zidDict

        l = Internal.newIntegralData(name='zid', parent=z) #ceate a zid for each block
        l[1] = iz

    for z in zs: # loop on blocks

        raccords = Internal.getNodesFromType2(z, 'GridConnectivity_t') # check if rac exists for block z
        nb_racs  = len(raccords)

        for rac in raccords: # loop on racs

            # GET ZONES INDICES
            z1 = getProperty(z, 'zid')

            donnorName = "".join(Internal.getValue(rac))
            z2 = zidDict[donnorName]

            # re-arrange
            if z1 > z2:
                c  = z2
                z2 = z1
                z1 = c

            ptList = Internal.getNodeFromName1(rac, 'PointList')[1][0]
            ptList_donor = Internal.getNodeFromName1(rac, 'PointListDonor')[1][0]

            is_perio = Internal.getNodeFromName1(rac, 'GridConnectivityProperty')
            idx_min  = numpy.minimum(numpy.min(ptList), numpy.min(ptList_donor))

            # 0==match connec; 1==perio connec
            if is_perio:
                idx_perio = 1
            else:
                idx_perio = 0

            # -----

            if z1 in dict_VD:
                if z2 in dict_VD[z1]:
                    # if 'idx_min' in dict_VD[z1][z2]:
                    if idx_perio in dict_VD[z1][z2]:
                        if idx_min in dict_VD[z1][z2][idx_perio]:
                            r    = Internal.newIntegralData(name='rid', parent=rac)
                            r[1] = dict_VD[z1][z2][idx_perio][idx_min]
                        else:
                            ir += 1
                            dict_VD[z1][z2][idx_perio][idx_min]=ir

                            r    = Internal.newIntegralData(name='rid', parent=rac)
                            r[1] = ir
                    else:
                        ir += 1
                        dict_VD[z1][z2][idx_perio] = {}
                        dict_VD[z1][z2][idx_perio][idx_min]=ir

                        r    = Internal.newIntegralData(name='rid', parent=rac)
                        r[1] = ir
                else:
                    ir += 1
                    dict_VD[z1][z2] = {}
                    dict_VD[z1][z2][idx_perio] = {}
                    dict_VD[z1][z2][idx_perio][idx_min]=ir

                    r    = Internal.newIntegralData(name='rid', parent=rac)
                    r[1] = ir
            else:
                ir += 1
                dict_VD[z1] = {}
                dict_VD[z1][z2] = {}
                dict_VD[z1][z2][idx_perio] = {}
                dict_VD[z1][z2][idx_perio][idx_min]=ir

                r    = Internal.newIntegralData(name='rid', parent=rac)
                r[1] = ir

#==============================================================================
# createHMesh : Returns a hierarchical mesh hook per zone
# IN: t : 3D NGON PyTree
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical zone hook
#==============================================================================
def createHMesh(t, subdiv_type=0):
    """Returns a hierarchical mesh hook per zone.
    Usage: createHMesh(t, subdiv_type= 0)"""

    zones = Internal.getZones(t)

    hmeshs = []
    for z in zones:

        m = C.getFields(Internal.__GridCoordinates__, z)[0]

        dico_rotation_to_ptList = get_transfo_to_ptlist_dico(z)

        m = intersector.initForAdaptCells(m, dico_rotation_to_ptList)

        zid = getProperty(z, 'zid')
        hm = intersector.createHMesh(m, subdiv_type, zid)
        hmeshs.append(hm)

    return hmeshs

#==============================================================================
# deleteHMesh : Releases a hierachical zone hook
# IN: hook : Python hook
# OUT: Nothing
#==============================================================================
def deleteHMesh(hooks):
    """Releases a hierachical zone hook.
    Usage: deleteHMesh(hooks)"""
    for h in hooks:
        if h == None : continue
        intersector.deleteHMesh(h)

#==============================================================================
# createSensor : Returns a sensor hook
# IN: t : 3D NGON PyTree
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical PyTree hook
#==============================================================================
def createSensor(hmeshs, sensor_type=0, smoothing_type=0 , itermax=-1, sensor_metric_policy=0):
    """Returns a sensor hook.
    Usage: createSensor(hmeshs, sensor_type = 0, smoothing_type=0 , itermax = -1)"""
    sensors = []
    for hmesh in hmeshs:
        if hmesh == None :
            sensors.append(None)
        else :
            sensors.append(intersector.createSensor(hmesh,sensor_type,smoothing_type,itermax,sensor_metric_policy))
    return sensors

#==============================================================================
# deleteSensor : Releases a sensor hook
# IN: hook : Python hook
# OUT: Nothing
#==============================================================================
def deleteSensor(hooks):
    """Releases a sensor hook.
    Usage: deleteSensor(hooks)"""
    for h in hooks:
        if h == None : continue
        intersector.deleteSensor(h)

#==============================================================================
# assignData2Sensor : Assigns data to a sensor
# IN: hook : Python hook
# OUT: Nothing
#==============================================================================
def assignData2Sensor(hooks, sensdata):
    """Assigns data to a sensor.
      Usage: assignData2Sensor(hooks, sensdata)"""
    sens_data_typ = InputType(sensdata)
    #print(sens_data_typ)
    if sens_data_typ == -1 :
        print('assignData2Sensor (geom) ERROR : wrong input data for sensor')
        return 1
    elif sens_data_typ == 2 : # single mesh/coords
        data = C.getFields(Internal.__GridCoordinates__, sensdata)[0]
        #print (data)
        for h in hooks:
            intersector.assignData2Sensor(h, data) # assign the single source point input cloud to all sensors

    elif sens_data_typ == 3: # top tree => join
        data = T.join(sensdata)
        data = C.getFields(Internal.__GridCoordinates__, data)[0]
        for h in hooks:
            intersector.assignData2Sensor(h, data) # assign the single source point input cloud to all sensors

    elif sens_data_typ == 4: # list of zones
        if len(sensdata) != len(hooks):
            data = T.join(sensdata)
            data = C.getFields(Internal.__GridCoordinates__, data)[0]
            for h in hooks:
                intersector.assignData2Sensor(h, data) # assign the single source point input cloud to all sensors
        else:
            i=-1
            for h in hooks:
                i+=1
                data = C.getFields(Internal.__GridCoordinates__, sensdata[i])[0]
                intersector.assignData2Sensor(h, data)

    elif sens_data_typ == 0: # single numpy
        if len(hooks) != 1:
            print('assignData2Sensor (nodal or centered) ERROR : data list must be sized as number of sensors')
            return 1
        intersector.assignData2Sensor(hooks[0], sensdata)

    elif sens_data_typ == 1 or sens_data_typ == 5: # list of numpies or list of list of numpies
        if len(sensdata) != len(hooks):
            print('assignData2Sensor (nodal or centered) ERROR : data list must be sized as number of sensors')
            return 1
        i=-1
        for h in hooks:
            i+=1
            intersector.assignData2Sensor(h, sensdata[i])
    return 0

def createCom(t, subdiv_type=0, nbz=1):
    return intersector.createCom(t, subdiv_type, nbz)

def deleteCom(hook):
    intersector.deleteCom(hook)

#==============================================================================
# conformizeHMesh : Converts the basic element leaves of a hierarchical mesh (hooks is a list of hooks to hiearchical zones) to a conformal polyhedral mesh.
#                   Each hiearchcial zone is referring to a zone in the original mesh t. So the mesh is replaced in the returned tree and the BCs/Joins/Fields are transferred.
# IN: t : PyTree before adaptation
# IN: hook : list of hooks to hiearchical zones (same size as nb of zones in t).
# OUT: Nothing
#==============================================================================
def conformizeHMesh(t, hooks, conformize=1):
    """Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
    Usage: conformizeHMesh(t, hooks)"""
    tp = Internal.copyRef(t)
    _conformizeHMesh(tp, hooks, conformize)
    return tp

#==============================================================================
# _conformizeHMesh: converts the basic element leaves of a hierarchical mesh (hooks is a list of hooks to hiearchical zones) to a conformal polyhedral mesh.
#                   Each hiearchcial zone is referring to a zone in the original mesh t. So the mesh is replaced in the tree and the BCs/Joins/Fields are transferred.
# IN: t: PyTree before adaptation
# IN: hook: list of hooks to hiearchical zones (same size as nb of zones in t).
# OUT: Nothing
#==============================================================================
def _conformizeHMesh(t, hooks, conformize=1):
    """Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
    Usage: _conformizeHMesh(t, hooks)"""
    nb_hooks = len(hooks)
    zones = Internal.getZones(t)
    nb_zones = len(zones)

    if nb_zones != nb_hooks:
        print('must give one hook per zone')
        return

    zidDict = {}
    for z in zones: # loop on blocks
        zidDict[z[0]] = getProperty(z, 'zid')

    ### 1. UPDATE BC & JOIN POINTLISTS WITH CURRENT ENABLING STATUS AND MPI EXCHANGES
    zone_to_rid_to_list_owned = getJoinsPtLists(t, zidDict)
    zone_to_bcptlists = getBCsPtLists(t)
    rid_to_zones = getRidToZones(t, zidDict)

    i=0
    for z in zones:

        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue
        if hooks[i] == None : continue

        zid = getProperty(z, 'zid')

        # BC and Joins point list (owned side)
        bcptlists = zone_to_bcptlists[zid]
        rid_to_ptlist = zone_to_rid_to_list_owned[zid]

        fieldsC = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
        if fieldsC == [] : fieldsC = None

        fieldsN = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        if fieldsN == [] : fieldsN = None

        try:
            pfieldsF = Internal.getNodeFromName(z, 'CADData')
            fieldsF = Internal.getChildFromName(pfieldsF, 'fcadid')
            fieldsF = [fieldsF[1]]
            if fieldsF == [] : fieldsF = None # Unnecessary ?

        except TypeError:
            fieldsF = None

        res = intersector.conformizeHMesh(hooks[i], bcptlists, rid_to_ptlist, rid_to_zones, fieldsC, fieldsN, fieldsF, conformize)

        # res[0] : mesh
        # res[1] : List : updated bc pointlist
        # res[2] : Dict : updated 'jzone to PtList'
        # res[3] : List : center fields
        # res[4] : List : node fields
        # res[5] : List : face fields

        mesh = res[0]

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        if len(res) < 2 : continue

        # MAJ BCs
        bcptlists = res[1]

        if bcptlists != [] :
            updateBCPointLists2(z, bcptlists)
        else:
            C._deleteZoneBC__(z)

        # MAJ Joins
        rid_to_ptlist = res[2]

        zone_to_rid_to_list_owned[zid] = rid_to_ptlist #dico z_to_rid_to_ptlist_owned est mis à jour avant de le donner à exchangePointList

        if rid_to_ptlist != {}:
            updateJoinsPointLists3(z, zidDict, rid_to_ptlist, 'PointList')

        # MAJ FIELDS

        C._deleteFlowSolutions__(z)

        ## center fields
        fieldz = [res[3]]
        #print (fieldz)
        if fieldz != [None]:
            for f in fieldz:
                C.setFields([f], z, 'centers', False)

        # # ## node fields
        fieldz = [res[4]]
        #print (fieldz)
        if fieldz != [None]:
            for f in fieldz:
                C.setFields([f], z, 'nodes', False)

        ## face fields
        fieldz = [res[5]]
        #print(fieldz)
        if fieldz is not [None]:
            Internal.newDataArray('fcadid', value=fieldz[0], parent=Internal.getNodeFromName(z, 'CADData'))

        i +=1
    _transposePointLists(t, hooks, zidDict, rid_to_zones, zone_to_bcptlists, zone_to_rid_to_list_owned)

#==============================================================================
# _transposePointLists : XXX
#==============================================================================
def _transposePointLists(t, hooks, zidDict, rid_to_zones=None, zone_to_bcptlists=None, zone_to_rid_to_list_owned=None):

    nb_hooks = len(hooks)
    zones = Internal.getZones(t)
    nb_zones = len(zones)

    if nb_zones != nb_hooks:
        print('must give one hook per zone')
        return

    if zone_to_rid_to_list_owned == None :
        zone_to_rid_to_list_owned = getJoinsPtLists(t, zidDict)

    if rid_to_zones == None:
        rid_to_zones = getRidToZones(t, zidDict)

    #
    zone_to_rid_to_list_opp = intersector.transposePointLists(rid_to_zones, zone_to_rid_to_list_owned)

    if zone_to_rid_to_list_opp == {} : return # single block
    #
    for z in zones:
        zid = getProperty(z, 'zid')
        if zid not in zone_to_rid_to_list_opp : continue

        rid_to_list_opp = zone_to_rid_to_list_opp[zid]
        if rid_to_list_opp == {} : continue

        updateJoinsPointLists3(z, zidDict, rid_to_list_opp, 'PointListDonor')

#==============================================================================
# interpolateHMeshNodalField : XXX
#==============================================================================
def interpolateHMeshNodalField(t, hooks, fname):
    """XXX"""
    tp = Internal.copyRef(t)
    _interpolateHMeshNodalField(tp, hooks, fname)
    return tp

#==============================================================================
# _interpolateHMeshNodalField : XXX
#==============================================================================
def _interpolateHMeshNodalField(t, hooks, fname):
    """XXX"""
    nb_hooks = len(hooks)
    zones = Internal.getZones(t)
    nb_zones = len(zones)

    if nb_zones != nb_hooks:
        print('must give one hook per zone')
        return
    i=0
    for z in zones:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue
        if hooks[i] == None : continue

        fieldsN = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        if fieldsN == [] : fieldsN = None

        #print(fieldsN)

        bcptlists = getBCPtList(z)
        #print(bcptlists)

        #recuperation du field ayant pour nom fname : u ou v
        for j in range(len(fname)):
            #print(fname[j])
            fieldsN[j] = C.getFields(Internal.__FlowSolutionNodes__, C.extractVars(z, fname[j]))

            fieldN = fieldsN[j]
            ofield = intersector.interpolateHMeshNodalField(hooks[i], fieldsN[j][0], bcptlists)

            ## MAJ node fields
            C.setFields([ofield], z, 'nodes', False)

        i=i+1

#==============================================================================
# adaptBox : Adapts a bounding box to a cloud of interior points
#==============================================================================
def adaptBox(t, box_ratio=10., smoothing_type=0, itermax=-1):
    """Adapts a bounding box to a cloud of interior points.
    Usage: adaptBox(t, box_ratio)"""

    z = Internal.getZones(t)[0]

    m = C.getFields(Internal.__GridCoordinates__, z)[0]

    m = intersector.adaptBox(m, box_ratio, smoothing_type, itermax)
    return C.convertArrays2ZoneNode('adapted', [m])

#==============================================================================
# extractUncomputables :XXX
#==============================================================================
def extractUncomputables(t,neigh_level=2):
    """Extracts any entity that will probably cause trouble to a CFD solver.
    Usage: extractUncomputables(t, neigh_level)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    res = XOR.extractUncomputables(m, neigh_level)

    zones = []
    nb_zones = len(res)

    if (nb_zones == 1) :
        zones.append(C.convertArrays2ZoneNode('remaining', [res[0]]))
    else:
        zones.append(C.convertArrays2ZoneNode('upgs', [res[0]]))
        zones.append(C.convertArrays2ZoneNode('uphs', [res[1]]))
        zones.append(C.convertArrays2ZoneNode('uphs_wv1', [res[2]]))
        zones.append(C.convertArrays2ZoneNode('remaining', [res[3]]))

    return zones

#==============================================================================
# extractPathologicalCells : Extracts cells that will cause a CFD solver run failure
# IN: t          : 3D NGON mesh
# IN: neigh_level: number of neighbor layers (surounding pathologies) to extract as well
# OUT: returns a 3D NGON Mesh with seprate zones for pathologies
#==============================================================================
def extractPathologicalCells(t, neigh_level=0):
    """Extracts all cells that will probably cause trouble to a CFD solver.
    Usage: extractPathologicalCells(t, neigh_level)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    res = XOR.extractPathologicalCells(m, neigh_level)

    zones = []

    nb_zones = len(res)
    if (nb_zones == 1) :
        zones.append(C.convertArrays2ZoneNode('okstar', [res[0]]))
    else:

        if (len(res[0][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('open_cells', [res[0]]))
        if (len(res[1][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('with_degen_PGs', [res[1]]))
        if (len(res[2][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('delaunay_failure', [res[2]]))
        if (len(res[3][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('non_star', [res[3]]))
        #if (len(res[4][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('neighbors', [res[4]]))
        #if (len(res[5][1][0]) != 0) : zones.append(C.convertArrays2ZoneNode('good', [res[5]]))

    return zones

#==============================================================================
# extractOuterLayers     : Extracts prescribed outer cell layers into a single zone
# IN : t:                : 3D NGON mesh
# IN : N:                : Number of layers to extract
# IN : discard_external  : for volume mesh with holes (e.g. external flow), set it to 1 to extract only layers around bodies.
# IN : output_remaining  : to add a zone in the ouptut tree with remaining elts
# OUT: r
#==============================================================================
def extractOuterLayers(t, N, discard_external=0, output_remaining=False):
    """ Extracts prescribed outer cell layers.
    Usage: extractOuterLayers(t, N, discard_external)"""
    zs = Internal.getZones(t)

    zones = []

    i=-1
    for z in zs:
        i +=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res = XOR.extractOuterLayers(m, N, discard_external)
        nb_zones = len(res)
        if (nb_zones == 1 and output_remaining == True) :
            zones.append(C.convertArrays2ZoneNode('remaining_'+str(i), [res[0]]))
        else:
            zones.append(C.convertArrays2ZoneNode('outers_'+str(i), [res[0]]))
            if output_remaining == True:
                zones.append(C.convertArrays2ZoneNode('remaining_'+str(i), [res[1]]))

    return zones

#==============================================================================
# extractNthCell : XXX
#==============================================================================
def extractNthCell(t, nth):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.extractNthCell(m, nth)

    zones = []
    nb_zones = len(m)-1

    zones.append(C.convertArrays2ZoneNode('cell_%d'%(nth), [m[0]]))

    if (nb_zones == 0) : return zones

    # here it has neighbors
    for i in range(nb_zones):
        zones.append(C.convertArrays2ZoneNode('neigh', [m[i+1]]))

    return zones

#==============================================================================
# extractNthFace : XXX
#==============================================================================
def extractNthFace(t, nth):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.extractNthFace(m, nth)

    zones = []
    nb_zones = len(m)-1

    zones.append(C.convertArrays2ZoneNode('face_%d'%(nth), [m[0]]))

    if (nb_zones == 0) : return zones

    # here it has parent elements
    for i in range(nb_zones):
        zones.append(C.convertArrays2ZoneNode('ph', [m[i+1]]))

    return zones

#==============================================================================
# removeNthCell : XXX
#==============================================================================
def removeNthCell(t, nth):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.removeNthCell(m, nth)
    return C.convertArrays2ZoneNode('mes_wo_%d'%(nth), [m])

#==============================================================================
# removeNthFace : XXX
#==============================================================================
def removeNthFace(t, nth):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.removeNthFace(m, nth)
    return C.convertArrays2ZoneNode('mes_wo_%d'%(nth), [m])

#==============================================================================
# extractBiggestCell : XXX
#==============================================================================
def extractBiggestCell(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.extractBiggestCell(m)

    zones = []
    nb_zones = len(m)-1

    zones.append(C.convertArrays2ZoneNode('cell', [m[0]]))

    if (nb_zones == 0) : return zones

    # here it has neighbors
    for i in range(nb_zones):
        zones.append(C.convertArrays2ZoneNode('neigh', [m[i+1]]))

    return zones

#==============================================================================
# extractBadVolCells : extract cells with bad volume / growth ratio
# IN: a          : 3D NGON mesh
# OUT: returns a single cell NGON mesh
#==============================================================================
def extractBadVolCells(t, ar=0.125, vmin=0., nneighs=0):
    """ Extracts bad cells based on gowth ratio
    Usage: extractBadVolCells(a)"""
    import sys;
    zones = Internal.getZones(t)
    om = []
    i=-1

    for z in zones:
        i+=1
        GEl = Internal.getElementNodes(z)
        NGON = 0; found = False
        for c in GEl:
            if c[1][0] == 22: found = True; break
            NGON += 1
        PE = None
        if found:
            node = GEl[NGON]
            PE = Internal.getNodeFromName1(node, 'ParentElements')
            if PE is None:
                print ('skipping zone %d as it does not have ParentElement'%i)
                continue
        else:
            print ('skipping zone %d as it does not have ParentElement'%i)
            continue

        print('extracting bad cells for zone %d'%i)
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res=XOR.extractBadVolCells(m, PE[1], ar, vmin, nneighs)
        res = C.convertArrays2ZoneNode('badcells_z_'+str(i), [res])
        om.append(res)

    return om

#==============================================================================
# extractOverConnectedCells : XXX
#==============================================================================
def extractOverConnectedCells(t, nneighs=0):
    """ XXX"""
    import sys;
    zones = Internal.getZones(t)
    om = []
    i=-1

    for z in zones:
        i+=1

        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res=XOR.extractOverConnectedCells(m, nneighs)
        res = C.convertArrays2ZoneNode('overcon_cells_z_'+str(i), [res])
        if nb_cells(res) > 0:
            print('extracting fully-over-connected cells for zone %d'%i)
            om.append(res)

    return om

#==============================================================================
# detectIdentitcalCells : detects (and optionally removes) geometrically identical cells
#======================================================================
def detectIdenticalCells(t, TOL=1.e-15, clean=0):
    """Detects (and optionannly removes) geometrically identical cells.
    Usage: detectIdenticalCells(t)"""
    return C.TZA1(t, 'nodes', 'nodes', True, XOR.detectIdenticalCells, TOL, clean)

def _detectIdenticalCells(t, TOL=1.e-15, clean=0):
    return C._TZA1(t, 'nodes', 'nodes', True, XOR.detectIdenticalCells, TOL, clean)

#==============================================================================
# detectOverConnectedFaces : detects Faces that belong to more than 2 cells in a mesh.
#======================================================================
def detectOverConnectedFaces(t):
    """Detects Faces that belong to more than 2 cells in a mesh."""
    return C.TZA1(t, 'nodes', 'nodes', True, XOR.detectOverConnectedFaces)

def _detectOverConnectedFaces(t, TOL=1.e-15, clean=0):
    return C._TZA1(t, 'nodes', 'nodes', True, XOR.detectOverConnectedFaces, TOL, clean)

#==============================================================================
# collapseMicroRegions : XXX
#======================================================================
def collapseSmallEdges(t, eratio, lmax=-1.):
    """Collapse small edges."""
    return C.TZA1(t, 'nodes', 'nodes', True, XOR.collapseSmallEdges, eratio, lmax)


def _collapseSmallEdges(t, eratio, lmax=-1.):
    return C._TZA1(t, 'nodes', 'nodes', True, XOR.collapseSmallEdges, eratio, lmax)

#==============================================================================
# getOverlappingFaces   : returns the list of polygons in a1 and a2 that are overlapping.
# IN : t1:              : NGON mesh (surface or volume).
# IN : t2:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# IN: amax              : maximal angular value (in rad) between the normals of each pair of colliding polygons.
#                         In ragnge [0,PI]. A value of 0. means pure parallelism. A value of PI means any collision.
# IN: dir2              : if specified, direction vector used for all a2's polygons instead of their own normals.
# OUT: 2 lists of overlapping polygons, the first one for a1, the seoncd one for a2.
#==============================================================================
def getOverlappingFaces(t1, t2, RTOL=0.1, amax=0.1, dir2=(0.,0.,0.)):
    """ Returns the list of polygons in a1 and a2 that are overlapping.
    Usage: getOverlappingFaces(t1, t2, RTOL, amax, dir2)"""

    try: import Transform as T
    except: raise ImportError("getOverlappingFaces: requires Transform module.")

    zones2 = Internal.getZones(t2)
    m2 = concatenate(zones2); m2 = G.close(m2)
    m2 = C.getFields(Internal.__GridCoordinates__, m2)[0]

    zones1 = Internal.getZones(t1)
    pgids = []

    i=-1
    for z in zones1:
        i+=1
        m1 = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m1 == []: continue

        pgids.append(XOR.getOverlappingFaces(m1,m2, RTOL, amax, dir2))

    return pgids

#==============================================================================
# getCollidingTopFaces  : Returns the list of top faces (HEXA and PRISM only) in t1 colliding t2.
# IN : t1:              : NGON mesh (surface or volume).
# IN : t2:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# OUT: list of t1 involved faces
#==============================================================================
def getCollidingTopFaces(t1, t2, RTOL=0.1):
    """ Returns the list of top faces (HEXA and PRISM only) in t1 colliding t2.
    Usage: getCollidingTopFaces(t1, t2, RTOL)"""

    try: import Transform as T
    except: raise ImportError("getCollidingTopFaces: requires Transform module.")

    zones2 = Internal.getZones(t2)
    m2 = concatenate(zones2); m2 = G.close(m2)
    m2 = C.getFields(Internal.__GridCoordinates__, m2)[0]

    zones1 = Internal.getZones(t1)
    pgids = []

    i=-1
    for z in zones1:
        i+=1
        m1 = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m1 == []: continue

        pgids.append(XOR.getCollidingTopFaces(m1,m2, RTOL))

    return pgids

#==============================================================================
# getCollidingCells     : returns the list of cells in a1 and a2 that are colliding.
# IN : t1:              : NGON mesh (surface or volume).
# IN : t2:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# OUT: 2 lists of colliding cells, the first one for t1, the seoncd one for t2.
#==============================================================================
def getCollidingCells(t1, t2, RTOL=1.e-12, only_externals=False):
    """ Returns the list of cells in a1 and a2 that are colliding.
    Usage: getOverlappingFaces(t1, t2, RTOL)"""

    try: import Transform as T
    except: raise ImportError("getCollidingCells: requires Transform module.")

    zones2 = Internal.getZones(t2)
    m2 = concatenate(zones2); m2 = G.close(m2)
    m2 = C.getFields(Internal.__GridCoordinates__, m2)[0]

    zones1 = Internal.getZones(t1)
    pgids = []

    i=-1
    for z in zones1:
        i+=1
        m1 = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m1 == []: continue

        pgids.append(XOR.getCollidingCells(m1,m2, RTOL, only_externals))

    return pgids

#==============================================================================
# getAttachedCells     : returns the cells in t1 attached to polygons in s2.
# IN : t1:              : NGON mesh (volume).
# IN : t2:              : NGON mesh (surface).
# OUT: attached cells.
#==============================================================================
def getAttachedCells(t1, s2):
    """ Returns the cells in t1 attached to polygons in s2.
    Usage: getAttachedCells(t1, s2)"""
    centsS = centroids(s2)

    pgids = getFaceIdsWithCentroids(t1,s2)

    return getCells(t1, pgids)

#==============================================================================
# getFacesWithCentroids     : returns the cells in t1 attached to polygons in s2.
# IN : t1:              : NGON mesh (volume).
# IN : t2:              : NGON mesh (surface).
# OUT: attached cells.
#==============================================================================
def getFaceIdsWithCentroids(t1, cents):
    """ Returns the faces in t1 having their centroids in cents.
    Usage: getFacesWithCentroids(t1, cents)"""
    c = C.getFields(Internal.__GridCoordinates__, cents)[0]

    zones1 = Internal.getZones(t1)
    pgids=[]
    for z in zones1:
        m1 = C.getFields(Internal.__GridCoordinates__, z)[0]
        pgids.append(XOR.getFaceIdsWithCentroids(m1,c))

    return pgids

#==============================================================================
# getFaceIdsCollidingVertex     : returns the Faces in t1 colliding a cloud of vertices
# IN : t1:              : NGON mesh (volume).
# IN : t2:              : NGON mesh (surface).
# OUT: attached cells.
#==============================================================================
def getFaceIdsCollidingVertex(t1, vtx):
    """ Returns the faces in t1 having their centroids in cents.
    Usage: getFaceIdsCollidingVertex(t1, [x,y,z])"""

    zones1 = Internal.getZones(t1)
    pgids=[]
    for z in zones1:
        m1 = C.getFields(Internal.__GridCoordinates__, z)[0]
        pgids.append(XOR.getFaceIdsCollidingVertex(m1, vtx))

    return pgids

#==============================================================================
# getCells              : Returns the cells in t1 having faces or cell ids
# IN : t1:              : NGON mesh (volume).
# IN : ids:             : face or cell ids.
# IN : are_face_ids:    : boolean to tell if ids are referring to faces or cells.
# OUT: attached cells.
#==============================================================================
def getCells(t1, ids, are_face_ids=True):
    """ Returns the cells in t1 having faces or cell ids.
    Usage: getCells(t1, pgids, are_face_ids = True)"""

    zones1 = Internal.getZones(t1)
    cells = []

    if len(ids) != len(zones1):
        print('getCells Input error : specified ids list must be sized as nb of zones')
        return None

    i=-1
    for z in zones1:
        i+=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res = XOR.getCells(m, ids[i], are_face_ids)
        narr = res[0] # mesh
        oids = res[1] # old cell ids

        if len(narr[2][0]) == 4 : continue #no connectivity

        #print(narr[2])
        #print(oids)
        newz = C.convertArrays2ZoneNode('zone%s'%i,[narr])
        if oids[1] != [] : C.setFields([oids], newz, 'centers', False)

        cells .append(newz)

    return cells

#==============================================================================
# getFaces              : returns the faces in t1 with ids in pgids.
# IN : t1:              : NGON mesh (volume).
# IN : pgids:           : polygon ids.
# OUT: attached cells.
#==============================================================================
def getFaces(t1, pgids):
    """ Returns the faces in t1 with ids in pgids.
    Usage: getFaces(t1, pgids)"""

    zones1 = Internal.getZones(t1)
    cells = []

    i=-1
    for z in zones1:
        i+=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        cells .append(C.convertArrays2ZoneNode('face%s'%i, [XOR.getFaces(m, pgids[i])]))

    return cells


#==============================================================================
# getNthNeighborhood     : returns the list of cells in the N-thneighborhood of t cells given in ids
# IN : t :               : NGON mesh.
# IN : N :               : number of neighborhood required
# IN : ids :             : input cells ids
# OUT: Returns the list of cells in the N-th neighborhood.
#==============================================================================
def getNthNeighborhood(t, N, ids):
    """ Returns the list of cells in the N-th neighborhood of cells given in ids.
    Usage: getNthNeighborhood(t, N, ids)"""

    zones = Internal.getZones(t)
    if len(ids) != len(zones) :
        print ('getNthNeighborhood : input ERROR : ids list and nb of zones are different')
        return

    idsNeigh = []

    i=-1
    for z in zones:
        i+=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue

        idsNeigh.append(XOR.getNthNeighborhood(m, N, ids[i]))

    return idsNeigh

#==============================================================================
# estimateAdapReq     : estimates an cell-specified adaptation requirement from
#                       on a istotropic metric field based on donnor connectivity
# IN : t :               : NGON mesh.
# IN : donnor :          : donnor mesh (Volume; Surface or BAR)
# IN : metric_policy :   : 0 : ISO_MIN (based on min edge length at nodes), ISO_MEAN : 1 (average edge length), ISO_MAX : 2 (based on max edge length)
# IN : rtol :            : relative tolerance for collisions detection
# IN : minv :            : lower threshold for nb of subdivisions
# IN : maxv :            : upper threshold for nb of subdivisions
# OUT: Returns a list of integers sized as the nb of cells in t giving the nb of subdivisions
#      per cell in the range [minv, maxv].
#==============================================================================
def estimateAdapReq(t, donnor, metric_policy=2, rtol=1.e-12, minv=0, maxv=5):
    """ estimates an cell-specified adaptation requirement from on a istotropic metric field based on donnor connectivity.
     Usage : estimateAdapReq(t, donnor [, metric_policy, rtol, minv, maxv])"""
    zones2 = Internal.getZones(donnor)
    m2 = concatenate(zones2); m2 = G.close(m2)
    m2 = C.getFields(Internal.__GridCoordinates__, m2)[0]

    zones = Internal.getZones(t)
    cell_vals = []
    i=-1
    for z in zones:
        i+=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue

        cell_vals.append(XOR.estimateAdapReq(m, m2, metric_policy, rtol, minv, maxv))

    return cell_vals

#==============================================================================
# getAnisoInnerFaces   : returns the list of polygons in a1 that are connecting 2 aniso elements.
# IN : t1:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# IN: dir2              : if specified, direction vector used for all a2's polygons instead of their own normals.
# OUT: 2 lists of overlapping polygons, the first one for a1, the seoncd one for a2.
#==============================================================================
def getAnisoInnerFaces(t1, aniso_ratio=0.05):
    """ Returns the list of polygons in a1 that are connecting 2 aniso elements.
    Usage: getAnisoInnerFaces(t1, aniso_ratio)"""

    try: import Transform as T
    except: raise ImportError("getAnisoInnerFaces: requires Transform module.")

    zones = Internal.getZones(t1)

    pgids = []

    for z in zones:

        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue

        pgids.append(XOR.getAnisoInnerFaces(m, aniso_ratio))

    return pgids

#==============================================================================
# diffMesh : Returns the diff between 2 meshes as 2 zones.
# IN: t1, t2:               : two 3D NGON mesh to check
# OUT: Returns 2 zones : one with the a1 cells that are not in a2, the other is the reciprocal.
#==============================================================================
def diffMesh(t1, t2):
    """ Returns the difference between 2 meshes as 2 zones.
    Usage: diffMesh(t1, t2)"""

    zones = []

    z1s = Internal.getZones(t1)
    z2s = Internal.getZones(t2)

    if len(z1s) != len(z2s):
        print("inputs have different nb of zones") ; return zones

    zid=-1
    for z1 in z1s:
        zid +=1
        z2 = z2s[zid]
        m1 = C.getFields(Internal.__GridCoordinates__, z1)[0]
        m2 = C.getFields(Internal.__GridCoordinates__, z2)[0]

        res = XOR.diffMesh(m1, m2)

        nb_zones = len(res)

        if (nb_zones == 0) : # fixme : never happen
            continue

        if res[0][1].size == 0 and res[0][1].size == 0 :
            continue

        zones.append(C.convertArrays2ZoneNode('z1', [res[0]]))
        zones.append(C.convertArrays2ZoneNode('z2', [res[1]]))

    if len(zones) == 0:
        print("No differences")

    return zones

#==============================================================================
# selfX : Checks self-intersections in a mesh
# IN: t:               : 3D NGON mesh
# OUT: Returns the first two cell ids that collide.
#==============================================================================
def selfX(t):
    """ Checks  self-intersections in a mesh.
    Usage: selfX(t)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.selfX(m)
    return C.convertArrays2ZoneNode('selfX', [m])

#==============================================================================
# checkCellsClosure : Returns the first cell id that is non-closed.
# IN: t:               : 3D NGON mesh
# OUT: Returns the first cell id that is non-closed
#==============================================================================
def checkCellsClosure(t):
    """ Returns the first cell id that is open
    Usage: checkCellsClosure(a)"""
    zones = Internal.getZones(t)
    iz=-1
    for z in zones:
        iz +=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        err = XOR.checkCellsClosure(m)
        if err == 1 :
            print('checkCellsClosure found an open cell in ' + str(iz) + '-th zone')
            return 1
    return 0

#==============================================================================
# checkCellsFlux : Computes the cell fluxes using the ParentElement node
#==============================================================================
def checkCellsFlux(t):
    """ Returns the cell id for which the Gauss flux is the greatest
    Usage: checkCellsFlux(a, PE)"""
    import sys;
    zones = Internal.getZones(t)
    maxflux=-sys.float_info.max
    cellid = -1
    zoneid=-1
    i=0
    for z in zones:
        GEl = Internal.getElementNodes(z)
        NGON = 0; found = False
        for c in GEl:
            if c[1][0] == 22: found = True; break
            NGON += 1
        PE = None
        if found:
            node = GEl[NGON]
            PE = Internal.getNodeFromName1(node, 'ParentElements')
            if PE is None:
                print ('skipping zone %d as it does not have ParentElement'%i)
                continue
        else:
            print ('skipping zone %d as it does not have ParentElement'%i)
            continue

        print('checking nullity of fluxes for zone %d'%i)
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res=XOR.checkCellsFlux(m, PE[1])
        if res[1] > maxflux:
            maxflux=res[1]
            cellid=res[0]
            zoneid=i
        i+=1
    return (maxflux, cellid, zoneid)


#==============================================================================
# checkAngularExtrema : Computes the min/max dihedral angles in a mesh
#==============================================================================
def checkAngularExtrema(t):
    """ Returns the cell id for which the Gauss flux is the greatest
    Usage: checkCellsFlux(a, PE)"""
    import sys;
    zones = Internal.getZones(t)
    maxA=-sys.float_info.max
    minA = sys.float_info.max
    mincellid = -1
    minzoneid=-1
    maxzoneid=-1
    maxcellid=-1
    i=0
    for z in zones:
        GEl = Internal.getElementNodes(z)
        NGON = 0; found = False
        for c in GEl:
            if c[1][0] == 22: found = True; break
            NGON += 1
        PE = None
        if found:
            node = GEl[NGON]
            PE = Internal.getNodeFromName1(node, 'ParentElements')
            if PE is None:
                print ('skipping zone %d as it does not have ParentElement'%i)
                continue
        else:
            print ('skipping zone %d as it does not have ParentElement'%i)
            continue

        print('checking nullity of fluxes for zone %d'%i)
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res=XOR.checkAngularExtrema(m, PE[1])
        if res[1] > maxA:
            maxA=res[1]
            maxcellid=res[0]
            maxzoneid=i
        if res[3] < minA:
            minA=res[3]
            mincellid=res[2]
            minzoneid=i
        i+=1
    return (minA, mincellid, minzoneid, maxA, maxcellid, maxzoneid)

#==============================================================================
# checkCellsVolume : Computes the cell fluxes using the ParentElement node
#==============================================================================
def checkCellsVolume(t):
    """ XXX"""
    import sys;
    zones = Internal.getZones(t)
    vmin=sys.float_info.max
    cellid = -1
    zoneid=-1
    i=-1
    for z in zones:
        i+=1
        GEl = Internal.getElementNodes(z)
        NGON = 0; found = False
        for c in GEl:
            if c[1][0] == 22: found = True; break
            NGON += 1
        PE = None
        if found:
            node = GEl[NGON]
            PE = Internal.getNodeFromName1(node, 'ParentElements')
            if PE is None:
                print ('skipping zone %d as it does not have ParentElement'%i)
                continue
        else:
            print ('skipping zone %d as it does not have ParentElement'%i)
            continue

        print('checking vol min for zone %d'%i)
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res=XOR.checkCellsVolume(m, PE[1])
        if res[1] < vmin:
            vmin=res[1]
            cellid=res[0]
            zoneid=i

    return (vmin, cellid, zoneid)

#==============================================================================
# checkCellsVolumeAndGrowthRatio : Computes the cell fluxes using the ParentElement node
#==============================================================================
def checkCellsVolumeAndGrowthRatio(t):
    """ XXX"""
    import sys;
    zones = Internal.getZones(t)
    vmin=sys.float_info.max
    ivmin = -1
    vzoneid=-1
    grmin=sys.float_info.max
    igrmin = -1
    grzoneid=-1
    i=-1
    for z in zones:
        i+=1
        GEl = Internal.getElementNodes(z)
        NGON = 0; found = False
        for c in GEl:
            if c[1][0] == 22: found = True; break
            NGON += 1
        PE = None
        if found:
            node = GEl[NGON]
            PE = Internal.getNodeFromName1(node, 'ParentElements')
            if PE is None:
                print ('skipping zone %d as it does not have ParentElement'%i)
                continue
        else:
            print ('skipping zone %d as it does not have ParentElement'%i)
            continue

        print('checking min vol and growth ratio for zone %d'%i)
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res=XOR.checkCellsVolumeAndGrowthRatio(m, PE[1])
        #print('vmin for this zone : '+str(res[1]))
        if res[1] < vmin:
            vmin=res[1]
            ivmin=res[0]
            vzoneid=i
        #print('grmin for this zone : '+str(res[3]))
        if res[3] < grmin:
            grmin=res[3]
            igrmin=res[2]
            grzoneid=i

    return (vmin, ivmin, vzoneid, grmin, igrmin, grzoneid)
#==============================================================================
# checkForDegenCells : XXX
#==============================================================================
def checkForDegenCells(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    return XOR.checkForDegenCells(m)

#==============================================================================
# checkForBigCells : XXX
#==============================================================================
def checkForBigCells(t, n):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.checkForBigCells(m, n)
    return C.convertArrays2ZoneNode('big', [m])

#==============================================================================
# oneph : XXX
#==============================================================================
def oneph(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    return XOR.oneph(m)

#==============================================================================
# edgeLengthExtrema : XXX
#==============================================================================
def edgeLengthExtrema(t):
    zones = Internal.getZones(t)
    Lmin = 10000000
    for z in zones:
        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        L = XOR.edgeLengthExtrema(coords)
        Lmin = min(L, Lmin)
    #print('min over zones is ', Lmin)
    return Lmin

#==============================================================================
# computeGrowthRatio : Returns a field of aspect ratio
# IN: t    : 3D NGON mesh
# IN: vmim : volume threshold
# OUT: Returns the first cell id that is non-closed
#==============================================================================
def computeGrowthRatio(t, vmin=0.):
    tp = Internal.copyRef(t)
    _computeGrowthRatio(tp, vmin)
    return tp

def _computeGrowthRatio(t, vmin=0.):
    """ Returns a field of aspect ratio.
    Usage: computeGrowthRatio(t, vmin)"""
    zones = Internal.getZones(t)
    i=-1
    for z in zones:
        i+=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        ar = XOR.computeGrowthRatio(m, vmin)
        C.setFields([ar], z, 'centers', False)

#==============================================================================
# extrudeUserDefinedBC : XXX
#==============================================================================
def extrudeBC(t, height=0.25, mean_or_min=1, create_ghost=1, bndType='UserDefined'):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    cur_shift=0
    extrudepgs=[]
    zones = Internal.getZones(t)
    #print("nb of zones %d"%(len(zones)))
    (extrudepgs, cur_shift) = concatenateBC(bndType, [zones], extrudepgs, cur_shift)
    if (extrudepgs != []) : extrudepgs = numpy.concatenate(extrudepgs) # create a single list
    else : return t
    #print("nb of pgs to pass : %s" %(len(extrudepgs)))


    mo = XOR.extrudeBC(m, extrudepgs, height, mean_or_min, create_ghost)

    return C.convertArrays2ZoneNode('extruded', [mo])

#==============================================================================
# extrudeSurf : XXX
#==============================================================================
def extrudeSurf(t, layer_height, nlayers=1, strategy=1):
    """XXX"""
    return C.TZA1(t, 'nodes', 'nodes', True, XOR.extrudeSurf, layer_height, nlayers, strategy)

def _extrudeSurf(t, layer_height, nlayers=1, strategy=1):
    return C._TZA1(t, 'nodes', 'nodes', True, XOR.extrudeSurf, layer_height, nlayers, strategy)

#==============================================================================
# extrudeRevolSurf : XXX
#==============================================================================
def extrudeRevolSurf(t, ax_pt, ax_dir, nlayers=1):
    """XXX"""
    return C.TZA1(t, 'nodes', 'nodes', True, XOR.extrudeRevolSurf, ax_pt, ax_dir, nlayers)

def _extrudeRevolSurf(t, ax_pt, ax_dir, nlayers=1):
    return C._TZA1(t, 'nodes', 'nodes', True, XOR.extrudeRevolSurf, ax_pt, ax_dir, nlayers)

#==============================================================================
# statsUncomputableFaces : XXX
#==============================================================================
def statsUncomputableFaces(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    return XOR.statsUncomputableFaces(m)

#==============================================================================
# statsSize : XXX
#==============================================================================
def statsSize(t, compute_metrics=1):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    return XOR.statsSize(m, compute_metrics)

#==============================================================================
# removeBaffles : XXX
#==============================================================================
def removeBaffles(t):
    return C.TZA1(t, 'nodes', 'nodes', True, XOR.removeBaffles)

#==============================================================================
# removeBaffles : XXX
#==============================================================================
def _removeBaffles(t):
    return C._TZA1(t, 'nodes', 'nodes', True, XOR.removeBaffles)

#==============================================================================
# convert2Polyhedron : XXX
#==============================================================================
def convert2Polyhedron(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.convert2Polyhedron(m)
    return C.convertArrays2ZoneNode('ph', [m])

#==============================================================================
# oneZonePerCell : output a PyTree with a zone per cell
#==============================================================================
def oneZonePerCell(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.oneZonePerCell(m)
    zones = []
    nb_zones = len(m)

    print(nb_zones)

    if nb_zones == 0: return zones

    # here it has parent elements
    for i in range(nb_zones):
        zones.append(C.convertArrays2ZoneNode('cell%s'%i, [m[i]]))

    return zones

#==============================================================================
# oneZonePerFace : output a PyTree with a zone per face
#==============================================================================
def oneZonePerFace(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.oneZonePerFace(m)
    zones = []
    nb_zones = len(m)

    print(nb_zones)

    if nb_zones == 0: return zones

    # here it has parent elements
    for i in range(nb_zones):
        zones.append(C.convertArrays2ZoneNode('face%s'%i, [m[i]]))

    return zones

#==============================================================================
# convertNGON2DToNGON3D : Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
# IN: t    : 3D NGON mesh
# OUT: Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
#==============================================================================
def convertNGON2DToNGON3D(t):
    tp = Internal.copyRef(t)
    _convertNGON2DToNGON3D(tp)
    return tp

def _convertNGON2DToNGON3D(t):
    """ Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
    Usage: convertNGON2DToNGON3D(t)"""
    C._deleteFlowSolutions__(t) # fixme
    zones = Internal.getZones(t)
    zo = []
    for z in zones:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        m = XOR.convertNGON2DToNGON3D(m)
        C.setFields([m], z, 'nodes')
    return None
def convertBasic2NGONFaces(t):
    tp = Internal.copyRef(t)
    _convertBasic2NGONFaces(tp)
    return tp

def _convertBasic2NGONFaces(t):
    """ Converts a Basic type format for faces (QUAD or TRI) to nuga Face/Node Format.
    Usage: _convertBasic2NGONFaces(t)"""
    zones = Internal.getZones(t)
    zo = []
    for z in zones:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        m = XOR.convertBasic2NGONFaces(m)
        C.setFields([m], z, 'nodes')
    return None
def convertTree2NUGANGON(t, keep_BC=False):
    tp = Internal.copyRef(t)
    _convertTree2NUGANGON(tp, keep_BC)
    return tp

def _convertTree2NUGANGON(t, keep_BC=False):

    if keep_BC == True:
        (BCs,BCNames,BCTypes) = C.getBCs(t)

    zones = Internal.getZones(t)
    for z in zones:
        dims = Internal.getZoneDim(z)
        firstzonetype = dims[0]
        #print(firstzonetype)
        # Go Unstructured
        if firstzonetype == 'Structured' :
            #print('Structured -> Basic')
            C._convertArray2Hexa(z)

        (typ, d) = getZoneNSTypeAndDim(z)

        #print((typ,d))
        if typ == 'BASIC':
            #print('BASic -> NGON cassiopee')
            C._convertArray2NGon(z)
            if d == 2 : typ = 'NGON_CASSIOPEE'
        # now it s a NGON (cassiopee or NUGA)
        if typ == 'NGON_CASSIOPEE':
            #print('NGON cassiopee -> NUGA')
            _convertNGON2DToNGON3D(z)

    if keep_BC == True:
        C._recoverBCs(t,(BCs,BCNames,BCTypes))

def centroids(t):
    tj = T.join(t)
    m = C.getFields(Internal.__GridCoordinates__, tj)[0]
    c = XOR.centroids(m)
    return C.convertArrays2ZoneNode('centroids', [c])

def volumes(t, algo=1, all_pgs_convex=False):
    tp = Internal.copyRef(t)
    _volumes(tp, algo, all_pgs_convex)
    return tp

def _volumes(t, algo=1, all_pgs_convex=False):
    zones = Internal.getZones(t)
    for z in zones:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        vols = XOR.volumes(m, algo, all_pgs_convex)
        C.setFields([vols], z, 'centers', False)

def merge(tz, sz, tol=1.e-15): #target zone, list source zones
    m = C.getFields(Internal.__GridCoordinates__, tz)[0]
    #print m
    s = C.getFields(Internal.__GridCoordinates__, sz)[0]
    #print s
    return XOR.merge(m, s, tol)

def concatenate(t, tol=1.e-15):
    zones = Internal.getZones(t)


    ms = []
    for z in zones:
        m  = C.getFields(Internal.__GridCoordinates__, z)[0]
        ms.append(m)

    res = XOR.concatenate(ms, tol)
    m   = res[0]

    nzone    = len(zones)
    z_pgnids = []
    z_phnids = []

    for i in range(nzone):
        z_pgnids.append(res[i+1])

    for i in range(nzone):
        z_phnids.append(res[i+nzone+1])

    newz = C.convertArrays2ZoneNode('assembly', [m])

    # Restore BCs
    zoneBC = Internal.createUniqueChild(newz, "ZoneBC", 'ZoneBC_t')
    i = 0
    for z in zones:
        BCs  = Internal.getNodesFromType(z, "BC_t")
        for bc in BCs:
            pointList = Internal.getNodeFromName(bc, Internal.__FACELIST__)
            pointList = pointList[1][0]
            for k in range(len(pointList)):
                iface = pointList[k]-1
                nids  = z_pgnids[i][iface]
                pointList[k] = nids+1

            # Restore BCDATASET
            BCDataSets = Internal.getNodesFromType(bc, "BCDataSet_t")

            for BCDataSet in BCDataSets:
                # Check size
                for BCData in Internal.getNodesFromType(BCDataSet, "BCData_t"):
                    for data in Internal.getNodesFromType(BCData, "DataArray_t"):
                        if ( len(data[1]) != len(pointList)):
                            raise ValueError('Concatenate. BC size differs for BCDataSet.')

                Internal._addChild(bc,BCDataSet)

            Internal._addChild(zoneBC, bc)

        i = i+1

    famBC = Internal.getNodesFromType(t, "Family_t")
    if famBC is not None:
        Internal._addChild(newz, famBC)

    # Restore fields
    # --------------
    # 1/ Compute new elt size and varname list
    i        = 0
    new_size = 0
    varnames = []
    for z in zones:
        size_z   = numpy.count_nonzero(z_phnids[i] != -1) # new connectivity size
        new_size = new_size + size_z

        cont     = Internal.getNodeFromName(z, Internal.__FlowSolutionCenters__)

        if cont is not None:
            flds     = Internal.getNodesFromType1(cont, 'DataArray_t')

            for fld in flds:
                if fld[0] not in varnames:
                    varnames.append(fld[0])

        i        = i+1

    for varname in varnames:
        # a. concatenate all old data for varname variable
        varglob = None
        for z in zones:
            ncells  = C.getNCells(z)
            cont    = Internal.getNodeFromName(z, Internal.__FlowSolutionCenters__)

            if cont is not None:

                varnode = Internal.getNodeFromName(cont, varname)
                # Put zero if variable does not exist
                if varnode is None:
                    varnode = numpy.zeros(ncells, numpy.float64)
                else:
                    varnode = varnode[1]

            else: # no FlowSolution in z
                varnode = numpy.zeros(ncells, numpy.float64)

            if varglob is None:
                varglob = varnode.ravel('k')
            else:
                varglob = numpy.concatenate((varglob,varnode))



        # b. Create new variable
        C._initVars(newz, 'centers:'+varname, 0)

        # c. Extract and replace new data array
        varnew = numpy.zeros(new_size, numpy.float64)
        i      = 0
        ik     = 0
        for z in zones:
            nelt = len(z_phnids[i])
            for k in range(nelt):
                if (z_phnids[i][k] != -1):
                    varnew[ik] = varglob[z_phnids[i][k]]
                    ik         = ik + 1
            i = i+1

        contnew    = Internal.getNodeFromName(newz, Internal.__FlowSolutionCenters__)
        varnode    = Internal.getNodeFromName(contnew, varname)
        varnode[1] = varnew

    return newz

def drawOrientation(t):
    zones = Internal.getZones(t)
    zo = []
    i=-1
    for z in zones:
        i+=1
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        d=XOR.drawOrientation(m)
        zo.append(C.convertArrays2ZoneNode('z_'+str(i), [d]))
    return zo

#
def volume(t, fieldname=None):

    try: import Generator.PyTree as G
    except: raise ImportError("volume: requires Generator module.")

    fldname = None
    if fieldname is not None: fldname = 'centers:'+fieldname

    zones = Internal.getZones(t)
    v = 0.
    for z in zones:
        xcelln = None
        if fldname is not None: xcelln = C.getField(fldname, z)
        if xcelln is not None: xcelln = xcelln[0]

        z = C.convertArray2NGon(z)
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        v += XOR.volume(m, xcelln)
    return v


#==============================================================================
# syncMacthPeriodicFaces : force periodicity for faces that are supposed to be periodic
# IN: a                  : 3D NGON mesh
# IN: rotationCenter: coordinates of the center of rotation for the periodicity
# IN: rotationAngle: rotation axis for the periodicity (its norm gives the angle of rotation)
# IN: translation: translation vector for a translation periodicity
# IN: TOL: tolerance. A negative value give a relative tolerance base on min edge length
# OUT: returns a 3D NGON Mesh with synchronised faces
#==============================================================================
def syncMacthPeriodicFaces(t, rotationCenter=[0.,0.,0.],
                           rotationAngle=[0.,0.,0.],
                           translation=[0.,0.,0.], tol=-0.01,
                           unitAngle=None, reorient=True):
    tp = Internal.copyRef(t)
    _syncMacthPeriodicFaces(tp, rotationCenter, rotationAngle, translation, tol, unitAngle, reorient)
    return tp

def _syncMacthPeriodicFaces(t, rotationCenter=[0.,0.,0.],
                            rotationAngle=[0.,0.,0.],
                            translation=[0.,0.,0.], tol=-0.01,
                            unitAngle=None, reorient=True):

    # WARNING : currently implemented by applying it individually by zone

    if unitAngle in ['Radian', None]: rotationAngleR=rotationAngle
    elif unitAngle == 'Degree': rotationAngleR=[v*Internal.__DEG2RAD__ for v in rotationAngle]
    else: raise ValueError('syncMacthPeriodicFaces: value for unitAngle is not valid.')

    zs = Internal.getZones(t)
    for z in zs :
        if reorient == True : _reorient(z)
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        res = XOR.syncMacthPeriodicFaces(m, rotationCenter, rotationAngleR,
                                         translation, tol)
        msh = res[0]
        #print(msh)
        C.setFields([msh], z, 'nodes')



#~ def conservativeTransfer(a1, a2, tol=0., reconstruction_type=0):
        #~
        #~ s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
        #~ flowsol1 = C.getFields(Internal.__FlowSolutionCenters__, a1)[0]
        #~ s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
#~
        #~ flowsol2 = XOR.conservativeTransfer(s1, flowsol1, s2, tol, reconstruction_type)
        #~
        #~ C._deleteFlowSolutions__(a2)
        #~ return C.setFields([flowsol2], a2, 'centers')
        #~
#~ def totalMass(a1):
#~
        #~ s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
        #~ flowsol1 = C.getFields(Internal.__FlowSolutionCenters__, a1)[0]
        #~
        #~ return XOR.totalMass(s1, flowsol1)
        #~
#~
#~ def normL1(a1,a2):
        #~
        #~ s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
        #~ flowsol1 = C.getFields(Internal.__FlowSolutionCenters__, a1)[0]
        #~
        #~ s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
        #~ flowsol2= C.getFields(Internal.__FlowSolutionCenters__, a2)[0]
#~
        #~ return generator.normL1(s1, flowsol1, s2, flowsol2)
#~
def testmain(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    c = XOR.testmain(m)
