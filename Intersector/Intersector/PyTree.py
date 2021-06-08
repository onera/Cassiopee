"""Grid generation module.
"""
# 
# Python Interface to create PyTrees defining meshes
#
from . import Intersector as XOR
from . import intersector
import time

__version__ = XOR.__version__

import numpy

try: range = xrange
except: pass

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter
    import Transform.PyTree as T
    import Generator.PyTree as G
    import Post.PyTree as P
    import Connector.PyTree as X
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
# Retourne le nombre de cellules d'un maillage
#=============================================================================
def nb_cells(a):
  ncellsTot = 0
  zones = Internal.getNodesFromType2(a, 'Zone_t')
  for z in zones:
      dim = Internal.getZoneDim(z)
      np = dim[1]
      ncells = dim[2]
      ncellsTot += ncells
  return ncellsTot

#==============================================================================
# nb_faces : Returns teh numbe rof faces in t

# IN: t             : 3D NGON mesh

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
# getTreeDim : XXX

# IN: t             : 3D NGON mesh

# OUT:  XXX
#==============================================================================
def getTreeDim(t):
  zs = Internal.getZones(t)
  d = 0
  nb_elts_per_zone = 0
  for z in zs:
    dims = Internal.getZoneDim(z)
    if d == 0: d = dims[4]
    if dims[4] != d:return 0; # mixed type : not handled
    if dims[3] == 'NGON' : nb_elts_per_zone += nb_cells(z)

  if d == 3 and nb_elts_per_zone == len(zs) : d = 21 # NGON with one cell per zone => NUGA NGON => 2D

  return d

#==============================================================================
# getZoneNSTypeAndDim : XXX

# IN: t             : 3D NGON mesh

# OUT:  XXX
#==============================================================================
def getZoneNSTypeAndDim(z):
  dims = Internal.getZoneDim(z)
  d = dims[4]
  if d == 2:
    if dims[3] == 'NGON' : return ('NGON_CASSIOPEE', 2)
    else : return ('BASIC', 2)
  else: # 3D
    #print(nb_cells(z))
    if dims[3] == 'NGON' and nb_cells(z) == 1 : return ('NGON_NUGA',2)
    elif dims[3] == 'NGON' : return ('NGON', 3)
    else : return ('BASIC', 3)


#==============================================================================
# InputType : XXX

# IN: t             : 3D NGON mesh

# OUT:  XXX
#==============================================================================
# return values : 0(single numpy), 1(list of numpies), 2(single zone), 3 (PyTree), 4(list of zones)
def InputType(t): # fixme : based on first block only
  if isinstance(t, list):
    if isinstance(t[0], numpy.ndarray): return 1
    isnod = Internal.isStdNode(t)
    #print(isnod)
    if isnod == -1 or isnod == 0:
      #print('is std node')
      if XOR.isSingleZone(t):
        #print('is isSingleZone node')
        return 2
      if Internal.isTopTree(t) :
        #print('is isSingleZone node')
        return 3
      if XOR.isSingleZone(t[0]):
        #print('is isSingleZone t0')
        return 4
    else : return -1
  else:
    if isinstance(t, numpy.ndarray): return 0
    else: return -1 # error

#==============================================================================
# NGONBlock

# GOAL : converts into a single connex NGON zone

# IN: t             : 3D NGON mesh
# IN: TOL           : tolerance
# IN: nb_comps      : nb of connex parts

# OUT: returns the adapted feature
#==============================================================================
def NGONBlock(t, TOL, nb_comps, keep_BC=False):

    if keep_BC == True:
      (BCs,BCNames,BCTypes) = C.getBCs(t)

    t = C.convertArray2NGon(t); #t = G.close(t)
    t = T.join(t, tol=TOL)
    t = G.close(t, tol=TOL)

    # FIXME : commented because
    valid = isConformalNConnex(t, nb_comps)
    if valid == False:
      C.convertPyTree2File(t, 'bad_oper.cgns')
      print('Invalid operand. Might increase the tolerance')
      import sys; sys.exit()

    if keep_BC == True:
      C._recoverBCs(t,(BCs,BCNames,BCTypes))

    #C.convertPyTree2File(t, 't.cgns')
    zs = Internal.getZones(t)
    z = zs[0]
    return z


#==============================================================================
# isConformalNConnex

# GOAL : check if a component is conformal and connex

# IN: t             : 3D NGON mesh
# IN: nb_comps      : nb of connex parts

# OUT: returns the adapted feature
#==============================================================================
def isConformalNConnex(t, nb_comps):

    F = P.exteriorFaces(t)
    F = T.splitConnexity(F)
    zs = Internal.getZones(F)
    #print('nb of surfs zones: ' + str(len(zs)))
    if len(zs) != nb_comps : return False
    F = P.exteriorFaces(F)
    zs = Internal.getZones(F)
    #print('nb of edge zones: ' + str(len(zs)))
    if len(zs) > 0:
        if nb_faces(zs[0]) != 0 : return False
    return True


#==============================================================================
# computeMeshAssemblyParameters

# GOAL : Computes the tolerance (based on min edge length) and thr configuration span to nomalize it 

# IN: t             : 3D NGON mesh

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

# IN: t             : 3D NGON mesh

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
  if conformal == False:
    print('ERROR : non conformal mesh')
    import sys; sys.exit()
  else:
    print('OK : conformal')

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

  # VERIFICATION 4 : VOL MIN
  print("Check min cell volume ...")
  (vmin, cellid, zoneid) = checkCellsVolume(t)
  print('vol min : ', vmin)
  if vmin < VOLMIN_SOLVER:
      print('Boolean ERROR : too small cells detected : under solver threshold')

  if fullcheck == False:
    return

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
# Concatenation des PointList d un type de BC donne dans une liste de zones
#=============================================================================
def concatenateBC(bctype, zones, wallpgs, cur_shift):
    i=0
    for z in zones:
      c = C.getFields(Internal.__GridCoordinates__, z)

      if c == []: continue

      #print(' -- zone : %d / %d' %(i+1, len(zones)))
      i=i+1
      bnds = Internal.getNodesFromType(z, 'BC_t')
      #print(" -- this zone has %d boundaries"%(len(bnds)))
      #print(' -- cur shift %d' %(cur_shift))

      # GET THE WALL PGS FROM THE POINTLISTS
      for bb in bnds:
        if Internal.isValue(bb, bctype) == False: continue
          
        wpgs = bb[2][1][1][0] # POINTLIST NUMPY fixme : can be somewhere else in the array
        #print(wpgs)
        # SYNC THE POINTLIST BEFORE APPENDING  : SHIFT WITH THE CURRENT NB OF STORED POLYGONS
        id2 = numpy.empty(len(wpgs), numpy.int32)
        id2[:] = wpgs[:] + cur_shift
        wallpgs.append(id2)

      c = c[0]
      #print(c)
      #z_nb_pts= len(c[1][0])
      z_nb_pgs= c[2][0][0]
      #print(z_nb_pts)
      #print(z_nb_pgs)
      cur_shift += z_nb_pgs
    return (wallpgs, cur_shift)

# update BC and JOINS point lists given an indirection "new id to old id"
def updatePointLists(z, zones, oids):
    #print('updateBCPointLists1')
    updateBCPointLists1(z, oids)
    #print('updateJoinsPointLists')
    updateJoinsPointLists1(z, zones, oids)

def updateBCPointLists1(z, oids):
    bnds = Internal.getNodesFromType(z, 'BC_t')
    zname=z[0]

    ptLists = []
    for bb in bnds :
      ptl = Internal.getNodesFromType(bb, 'IndexArray_t')
      ptLists.append(ptl[0][1][0])

    if ptLists == []: return

    # recalcul des pointlist
    ptLists = XOR.updatePointLists(oids, ptLists)

    i=0
    # update the BC pointlists 
    for bb in bnds :
      ptl = Internal.getNodesFromType(bb, 'IndexArray_t')
      ptl[0][1] = ptLists[i]
      #print(ptl[0][1])
      i=i+1

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
      #donnorName = "".join(j[1])
      donnorName = "".join(Internal.getValue(j))
      ptl = Internal.getNodeFromName1(j, 'PointList')
      #print(donnorName)
      dz = Internal.getNodeFromName(zones, donnorName)
      joinsD = Internal.getNodesFromType(dz, 'GridConnectivity_t')
      for jd in joinsD:
        #dname = "".join(jd[1])
        dname = "".join(Internal.getValue(jd))
        #print(dname)
        if (dname != zname) : continue
        ptlD = Internal.getNodeFromName1(jd, 'PointListDonor')
        
        PG0 = ptl[1][0][0] # first polygon in the point list 
        PG0D = ptlD[1][0][0] # first polygon in the point list
        if (PG0 != PG0D) : continue # not the right join (in case of multiple joins for 2 zones) : the first PG must be the same (assume one PG only in one join)
        
        ptLists[i] = numpy.reshape(ptLists[i], (1,len(ptLists[i]))) #checkme : seems to be useless 

        ptl[1]= ptLists[i]
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

def getBCPtList(z):
  #
  bnds = Internal.getNodesFromType(z, 'BC_t')
  zname=z[0]

  ptLists = []
  for bb in bnds :
    ptList = Internal.getNodesFromType(bb, 'IndexArray_t')
    ptLists.append(ptList[0][1][0])
    #ptLists.append(ptList)

  #print (ptLists)
  return ptLists

def getBCPtListOfType(z, typesList, families = []):
    #
    nodes = []
    for btyp in typesList:
        nodes += Internal.getNodesFromValue(z, btyp)
        if families != []:nodes += C.getFamilyBCs(z, families)

    #print(nodes)
    ptList = []
    for n in nodes:
        ptlnod = Internal.getNodesFromType(n, 'IndexArray_t')
        #print (ptlnod)
        ptList.append(ptlnod[0][1][0])

    if (ptList != []) : ptList = numpy.concatenate(ptList).ravel() # create a single list

    return ptList

#------------------------------------------------------------------------------
# 
#------------------------------------------------------------------------------
def updateJoinsPointLists2(z, zones, jzids, ptLists):
  # zone name to id
  name2id = dict()
  i = 0
  for zz in zones :
    name2id[zz[0]] = i
    #print('zone ' + z[0] + ' has ' + str(i))
    i += 1

  #print('processed zone  : ' + z[0])

  joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
  zname=z[0]

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
      if (dname != zname) : continue
      ptlD = Internal.getNodeFromName1(jd, 'PointListDonor')
      
      PG0 = ptl[1][0][0] # first polygon in the poitn list 
      PG0D = ptlD[1][0][0] # first polygon in the poitn list
      if (PG0 != PG0D) : continue # not the right join (in case of multiple joins for 2 zones) : the first PG must be the same (assume one PG only in one join)
      
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

def intersection(surface1, surface2, tol=0.):
    """Computes the intersection between two closed-surface or two volume meshes.
    Usage: intersection(s1, s2, tol)"""
    s1 = C.getFields(Internal.__GridCoordinates__, surface1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, surface2)[0]
    s = XOR.intersection(s1, s2, tol)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the intersection between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanIntersection(a1, a2, tol)
    Usage for volumes: booleanIntersection(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanIntersection(s1, s2, tol, preserve_right, solid_right, agg_mode, improve_qual)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanUnion(a1, a2, tol=0., jtol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, multi_zone=False, simplify_pgs=True): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
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

          return booleanUnionMZ(a1, a2, tol, jtol, agg_mode, improve_qual, simplify_pgs)

    #multi_zone option is ignored from here

    cur_shift=0
    extrudepgs=[]
    if (solid_right == 1) :
        zones = Internal.getZones(a2)
        (extrudepgs, cur_shift) = concatenateBC('UserDefined', zones, extrudepgs, cur_shift)
    if (extrudepgs != []) : extrudepgs = numpy.concatenate(extrudepgs) # create a single list
    #print("nb of pgs to pass : %s" %(len(extrudepgs)))

    res = XOR.booleanUnion(s1, s2, tol, preserve_right, solid_right, agg_mode, improve_qual, extrudepgs, simplify_pgs)
    
    is_zone_list  = 0
    if (len(res) != 4) : is_zone_list = 1
    elif (res[3] != 'NGON' and res[3] != 'TRI' and res[3] != 'BAR') : is_zone_list = 1

    if (is_zone_list == 0) : return C.convertArrays2ZoneNode('union', [res])

    # debug : mutli zones
    ozones = []

    for i in range(len(res)):
      if len(res[i][0][1]) != 0: ozones.append(C.convertArrays2ZoneNode(res[i][1], [res[i][0]])) #(zname, array)
    
    return ozones

def booleanUnionMZ(t1, t2, xtol=0., jtol=0., agg_mode=1, improve_qual = False, simplify_pgs = True): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    tp1 = Internal.copyRef(t1)
    tp2 = Internal.copyRef(t2)
    return _booleanUnionMZ(tp1, tp2, xtol, jtol, agg_mode, improve_qual, simplify_pgs)


def _booleanUnionMZ(t1, t2, xtol=0., jtol=0., agg_mode=1, improve_qual = False, simplify_pgs = True): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the union between two closed volume meshes.
    Usage for volumes: booleanUnion2(a1, a2, tol, agg_mode)"""
    m1s = []
    z1s = Internal.getZones(t1)
    for z in z1s:
      m1s.append(C.getFields(Internal.__GridCoordinates__, z)[0])
    m2s = []
    z2s = Internal.getZones(t2)
    for z in z2s:
      m2s.append(C.getFields(Internal.__GridCoordinates__, z)[0])

    res = XOR.booleanUnionMZ(m1s, m2s, xtol, jtol, agg_mode, improve_qual, simplify_pgs)

    i=0
    zs = []
    for z in z1s:
        mesh = res[i]
        pg_oids=res[i+1]

        #print mesh

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        # MAJ POINT LISTS #
        #updatePointLists(z, z1s, pg_oids)
        i=i+2
        zs.append(z)

    for z in z2s:
        mesh = res[i]
        pg_oids=res[i+1]

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        # MAJ POINT LISTS #
        #updatePointLists(z, z2s, pg_oids)
        i=i+2
        zs.append(z)

    t1_is_tree = Internal.isTopTree(t1)
    t2_is_tree = Internal.isTopTree(t2)

    if t1_is_tree == True and t2_is_tree == True:
      return Internal.merge([t1, t2])
    else:
      return zs
      

def booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanMinus(a1, a2, tol)
    Usage for volumes: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanMinus(s1, s2, tol, preserve_right, solid_right, agg_mode)
    return C.convertArrays2ZoneNode('minus', [s])

def diffSurf(a1, a2, tol=0., preserve_right=1, agg_mode=1, improve_qual=False, outward_surf=True): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between a volume mesh and a surface mesh."""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.diffSurf(s1, s2, tol, preserve_right, agg_mode,improve_qual, outward_surf)
    return C.convertArrays2ZoneNode('VmS', [s])
    
def booleanModifiedSolid(solid, a2, tol=0., preserve_solid=1, agg_mode=1, improve_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the transformed input solid after solving the intersection of its skin with a2.
    Usage: booleanModifiedSolid(a1, a2, tol, preserve_right, solid_right)"""
    sld = C.getFields(Internal.__GridCoordinates__, solid)[0]
    operand = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanModifiedSolid(operand, sld, tol, preserve_solid, agg_mode, improve_qual)
    return C.convertArrays2ZoneNode('modified_solid', [s])

#==============================================================================
# XcellN
# IN : t: 3D NGON mesh
# IN : priorities : one-to-one priorities between components
# IN : output_type : 0 : binary mask; 1 : continuous mask (xcelln) ; 2 : clipped surface.
# IN : rtol : relative tolerance
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
# IN : priorities : one-to-one priorities between components
# IN : output_type : 0 : binary mask; 1 : continuous mask (xcelln) ; 2 : clipped surface. 
# OUT: returns a 3D NGON surface mesh with the xcelln field (if output_type=0/1, the clipped surface with solution if output_type=2)
#==============================================================================
def _XcellN(t, priorities, output_type=0, rtol=0.05):
  """Computes the weight coefficients of visibility for overset grid configurations as a field called xcelln, for any kind of surface mesh.
  Usage : _XcellN(t, priorities [, rtol])"""

  DIM = getTreeDim(t)
  #print ('DIM ? ' +str(DIM))
  if DIM != 2 and DIM != 3 : # and DIM != 21: # 21 NUGA SURFACE
    print ('XcellN : Input error : the input file has an unsupported format or contain mixed 2D/3D zones.')
    return

  _XcellN_(t, priorities, output_type, rtol)
  
#==============================================================================
# _XcellN_ (in-place version)
# IN: t: 3D NGON SURFACE mesh
# IN : priorities : one-to-one priorities between components
# IN : output_type : 0 : binary mask; 1 : continuous mask (xcelln) ; 2 : clipped surface. 
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

  if len(bases) == 1 :
    print('XcellN : Input error : Only one base in the file. Each component must be separated in a given Base. No check between zones of the same component.')
    return

  min_compid = min(min(priorities, key=min))
  max_compid = max(max(priorities, key=max))
  if max_compid < 0 :
    print('XcellN : Input error : Negativle values passes as priorities. mus be component id (0-based).')
    return
  if max_compid >= len(bases):
    print('XcellN : Input error : Greatest component specified in priorities exceeds nb of components.')
    return

  # 1 PREPARE INPUTS : NUGA NGON zones + oriented BAR boundaries with wall ids
  if TIMER == True:
    xcelln_time = time.time()
    xcelln_time2 = time.time()
  # 1.1 convert to NUGA NGON
  tNG = convertTree2NUGANGON(t, True) # keepBC
  #C.convertPyTree2File(tNG, 'tNG.cgns')
  
  # 1.2 reorient
  if DIM == 3 :
    _reorient(tNG)

  if TIMER == True:
    print ('XCellN : Preparing Inputs ::: NGON convert & reorientation ::: CPU time : ',time.time()-xcelln_time2,'s')
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
    bid +=1
    zj = concatenate(bNG, tol = 1.e-10) # discard inner joins
    #if DBG == True : C.convertPyTree2File(zj, 'zj_'+str(bid)+'.cgns')
    b_bounds = externalFaces(zj) # keeping orientation
    if DBG == True : C.convertPyTree2File(b_bounds, 'bound_b_'+str(bid)+'.cgns')
    
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
      if DIM == 2 :
        hook = C.createHook(b_bounds, function='elementCenters') # identifying edges
      elif DIM == 3 :
        hook = C.createHook(b_bounds, function='faceCenters')
      wallf = C.identifyElements(hook, walls) # wallf are ids in boundaries
      #print(wallf)
      wallf = wallf[wallf >= 1]
      if wallf != [] :
        wallf -= 1 # make it 0 based
      else : wallf = None

    wall_ids.append(wallf)

  if TIMER == True:
    print ('XCellN : Preparing Inputs ::: BC and Walls ::: CPU time : ',time.time()-xcelln_time2,'s')
    xcelln_time2 = time.time()

  # 1.4 get the zones in a single list with parent base id
  ngons = []
  basenum = []
  zwall_ids = [] # for double wall mgt
  base_id=-1
  for b in basesNG:
    base_id += 1
    zones = Internal.getZones(b)
    for z in zones:
        c = C.getFields(Internal.__GridCoordinates__, z)[0]
        ngons.append(c)
        basenum.append(base_id)
        zwallf = getBCPtListOfType(z, WALLBCS, wallfamilies)
        if (zwallf != []) : zwallf -= 1 # make it 0 based
        zwall_ids.append(zwallf)

  #print(zwall_ids)
  #import sys; sys.exit()

  if TIMER == True:
    print ('XCellN : Preparing Inputs : CPU time : ',time.time()-xcelln_time,'s')
    xcelln_time = time.time()

  # 2. COMPUTE THE COEFFS PER ZONE (PARALLEL OMP PER ZONE)
  xcellns = XOR.XcellN(ngons, zwall_ids, basenum, boundaries, wall_ids, priorities, output_type, rtol)

  if TIMER == True:
    print ('XCellN : Computing : CPU time : ',time.time()-xcelln_time,'s')
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

    if TIMER == True:
      print ('XCellN : Writing output : CPU time : ',time.time()-xcelln_time,'s')
      xcelln_time = time.time()

    return None

  #SET XCELLN FIELD tO ZONES
  # fixme : conversion here because tNG NUGA does not work with node2center
  #         use convertArray2NGon to ensure that t is converted whatever its type
  #         e.g convert2Hexa does not work if t is NGON
  tNG = C.convertArray2NGon(t)
  basesNG = Internal.getBases(tNG)

  zid = -1
  bid = -1
  has_structured_bases = False
  for b in bases:
    bid +=1
    zones = Internal.getZones(b)
    dims = Internal.getZoneDim(zones[0])
    if dims[0] == 'Structured' : # set field on tNG
      has_structured_bases = True
      zonesNG = Internal.getZones(basesNG[bid])
      for z in zonesNG:
        zid +=1
        #print('set field on tNG')
        C.setFields([xcellns[zid]], z, 'centers', False)
      mc = C.node2Center(basesNG[bid])
      hookC = C.createGlobalHook([mc], 'nodes')
      hookN = C.createGlobalHook([basesNG[bid]], 'nodes')

      C._identifySolutions(b, basesNG[bid], hookN, hookC, tol=1000.)
      C.freeHook(hookC)
      C.freeHook(hookN)
    else : # set field on t
      for z in zones:
        zid +=1
        #print('set field on t')
        C.setFields([xcellns[zid]], z, 'centers', False)

  if TIMER == True:
      print ('XCellN : Writing output : CPU time : ',time.time()-xcelln_time,'s')
    
  return None

def P1ConservativeInterpolation(tR, tD):
     
    tp = Internal.copyRef(tR)
    _P1ConservativeInterpolation(tp, tD)
    return tp

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
# triangulateExteriorFaces
# IN: mesh: 3D NGON mesh
# IN : in_or_out : 0 means "ONLY INTERNALS", 1 means "ONLY EXTERNALS", any other value means "BOTH"
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def triangulateExteriorFaces(t, in_or_out=2, improve_qual=0):
    """Triangulates exterior polygons of a volume mesh.
    Usage: triangulateExteriorFaces(t)"""
    return C.TZA(t, 'nodes', 'nodes', XOR.triangulateExteriorFaces, t, in_or_out, improve_qual)

def _triangulateExteriorFaces(t, in_or_out=2):
    return C._TZA(t, 'nodes', 'nodes', XOR.triangulateExteriorFaces, t, in_or_out, improve_qual)


#==============================================================================
# triangulateSpecifiedFaces 
# IN: a: 3D NGON mesh
# IN: pgs : list of polygons
# OUT: returns a 3D NGON Mesh
#==============================================================================
def triangulateSpecifiedFaces(t, pgs, improve_qual=1):
     
    tp = Internal.copyRef(t)
    _triangulateSpecifiedFaces(tp,pgs, improve_qual)
    return tp

def _triangulateSpecifiedFaces(t, pgs, improve_qual=1):

    zones = Internal.getZones(t)
    if (len(pgs) != len(zones)) :
        print('triangulateSpecifiedFaces: input error: nb of polygons packs differ from nb of zones.')
        return None

    i=0
    for z in zones:
      m = C.getFields(Internal.__GridCoordinates__, z)[0]
      if m == []: continue
      m = Converter.convertArray2NGon(coords)
      m = XOR.triangulateSpecifiedFaces(m, pgs[i], improve_qual)
      mesh = res[0]
      pg_oids=res[1]

      # MAJ du maillage de la zone
      C.setFields([mesh], z, 'nodes') 

      # MAJ POINT LISTS #
      updatePointLists(z, zones, pg_oids)
      i = i+1

#==============================================================================
# triangulateNonBasicFaces
# IN: mesh: 3D NGON mesh
# IN : quality improvement flag
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def triangulateNFaces(t, improve_qual=1, min_nvertices=5, discard_joins= True):
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
        if (discard_joins == True):
            joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
            for j in joins :
                ptl = Internal.getNodeFromName1(j, 'PointList')
                ptLists.append(ptl[1])

        if (ptLists != []):
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
    m = XOR.convexifyFaces(m)
    return C.convertArrays2ZoneNode('allPGconvex', [m])

#==============================================================================
# externalFaces : Returns erternal faces for CASSIOPEE NGON types and NUGA NGON
#==============================================================================
def externalFaces(t):
    """Returns erternal faces for CASSIOPEE NGON types and NUGA NGON.
    Usage: externalFaces(t)"""
    zs = Internal.getZones(t)
    efs = []

    i=-1
    for z in zs:
      i+=1
      coords = C.getFields(Internal.__GridCoordinates__, z)[0]
      if coords == []: continue
      ef = XOR.externalFaces(coords)
      efs.append(C.convertArrays2ZoneNode('ef_z'+str(i), [ef]))
    return efs

#==============================================================================
# reorient : reorients outward the external polygons of a mesh
#==============================================================================
def reorient(t, dir=1):
    """Reorients outward the external polygons of a mesh.
    Usage: reorient(t)"""
    return C.TZA(t, 'nodes', 'nodes', XOR.reorient, t, dir)

def _reorient(t, dir=1):
    return C._TZA(t, 'nodes', 'nodes', XOR.reorient, t, dir)

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
def prepareCellsSplit(t, PH_set=1, split_policy=0, PH_conc_threshold=1./3., PH_cvx_threshold=0.05, PG_cvx_threshold = 1.e-8):
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
def splitNonStarCells(t, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8):
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
def simplifyCells(t, treat_externals, angular_threshold = 1.e-12, discard_joins=True):
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

      ozones.append(C.convertArrays2ZoneNode(z[0], [m]))

    return ozones

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
def simplifySurf(t, angular_threshold = 1.e-12):
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
# OUT: returns a 3D NGON Mesh with less cells and with a smoother aspect ratio
#==============================================================================
def agglomerateSmallCells(t, vmin=0., vratio=1000., angular_threshold=1.e-12, force=False):
    """Agglomerates prescribed cells.
    Usage: agglomerateSmallCells(t, vmin, vratio)"""

    nb_phs0  = nb_cells(t)

    m = C.getFields(Internal.__GridCoordinates__, t)[0]

    res = XOR.agglomerateSmallCells(m, vmin, vratio, angular_threshold, force)
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

# def agglomerateSmallCells(mesh, vmin=0., vratio=1000.):
#     m = C.getFields(Internal.__GridCoordinates__, mesh)[0]
#     print("one")
#     res = XOR.agglomerateSmallCells(m, vmin, vratio)
#     print("NB ZONES %d"%len(res))

#     z = C.convertArrays2ZoneNode('agglomeratedCells', [res[0]])

#     debug = True

#     if (debug == False) : return z

#     nb_aggs = res[1]
#     nb_cels = nb_cells(z);
#     nb_points = len(res[0][1][0])
    
#     #print "NB AGG OVER NB CELLS : %d / %d "%(nb_aggs, nb_cels)

#     C._initVars(z, 'centers:cellN', 0)

#     cellN = Internal.getNodesFromName(z, 'cellN')[0][1]
#     cellN[0:nb_aggs] = numpy.arange(1,nb_aggs+1)

#     C.convertPyTree2File(z, 'agglo.cgns')

#     return z

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
def agglomerateCellsWithSpecifiedFaces(t, pgs, simplify=1, amax = 1.e-12, treat_externals=1): # 0 : dno not simplify, 1 : simplify only internals, 2 : simlplify evrywhere
     
    tp = Internal.copyRef(t)
    _agglomerateCellsWithSpecifiedFaces(tp,pgs, simplify, amax, treat_externals)
    return tp

def _agglomerateCellsWithSpecifiedFaces(t, pgs, simplify=1, amax = 1.e-12, treat_externals=1):

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
      if jids != []:
        for k in range(len(jids)):
          oj = jids[k]
          nj = nids[oj]
          jids[k]=nj
        jids = jids[jids > -1]
      else:
        jids = None

      if simplify == 1:
        m = XOR.simplifyCells(m, treat_externals, angular_threshold=amax, discarded_ids = jids)

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
def collapseSmallCells(t, vmin, ar):
    """Regularize the mesh to be cut (by boolean operations) by immersing its nodes lying on the cutter surface.
    Usage: immerseNodes(t,s, TOL)"""
    tp = Internal.copyRef(t)
    _collapseSmallCells(tp, vmin, ar)
    return tp

def _collapseSmallCells(t, vmin, ar):
    zones = Internal.getZones(t)
    for z in zones:
      coords = C.getFields(Internal.__GridCoordinates__, z)[0]
      if coords == []: continue
      collapsed = intersector.collapseSmallCells(coords, vmin, ar)
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
    if sdims[3] == 'NGON' :
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

    zones = Internal.getZones(t)

    for z in zones:
        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue

        coords = Converter.convertArray2NGon(coords)
        mesh = intersector.closeCells(coords)

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

    return t

#==============================================================================
# adaptCells : Adapts an unstructured mesh a with respect to a sensor
# IN: t : 3D NGON mesh
# IN: sensdata : sensor data (a bunch of vertices or a mesh for a geom sensor, a mesh for a xsensor, punctual values for a nodal or cell sensor)
# IN: sensor_type : geom_sensor (0) , xsensor (1), nodal_sensor (2), cell_sensor(3), xsensor2(4)
# IN smoothing_type : First-neighborhood (0) Shell-neighborhood(1)
# IN itermax : max number of level in the hierarchy
# IN: subdiv_type : isotropic currently
# IN: hmesh : hierarchical mesh hook
# IN: sensor : sensor hook
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None):
  """Adapts an unstructured mesh a with respect to a sensor.
  Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None)"""
  tp = Internal.copyRef(t)
  _adaptCells(tp, sensdata, sensor_type, smoothing_type, itermax, subdiv_type, hmesh, sensor)
  return tp

def _adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None):
    """Adapts an unstructured mesh a with respect to a sensor.
    Usage: adaptCells(t, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None)"""

    if sensdata is None and sensor is None:
      print('INPUT ERROR : no source data to initialize a sensor')
      return

    if hmesh is None and sensor is not None:
      print ('INPUT ERROR : you must also give as an argument the hmesh targeted by the sensor')
      return

    NBZ = len(Internal.getZones(t))

    owesHmesh=0
    com = None
    if hmesh == [None]: # no basic elts in t
      return 
    if hmesh is None :
      #print("create hm : ", NBZ)
      if NBZ == 1:
        hmesh = createHMesh(t, subdiv_type)
        if hmesh == [None] : # no basic elts in t (assumed because hmesh creation as done [None]
          return
      else :
      	(hmesh, com) = createHZones(t, subdiv_type)
      owesHmesh=1

    owesSensor=0
    if sensor is None : 
      #print("create sensor")
      sensor = createSensor(hmesh, sensor_type, smoothing_type, itermax)
      owesSensor=1

    err=0
    if sensdata is not None:
      #print("assignData2Sensor")
      if sensor_type == 4:
        sensdata = C.convertArray2NGon(sensdata)
      err = assignData2Sensor(sensor, sensdata)
      if err == 1:
        print('INPUT ERROR : sensor data list must be sized as nb of sensors')
        return

    #print('adaptCells..')
    intersector.adaptCells(hmesh, sensor)

    if owesHmesh == 1 : #and owesSensor == 1 :
    	#print("_conformizeHMesh")
    	_conformizeHMesh(t, hmesh)

    if owesHmesh == 1 :
    	#print('delete owned hmesh')
    	deleteHMesh(hmesh)
    if owesSensor == 1 : 
    	#print('delete owned sensor')
    	deleteSensor(sensor)
    	#print('delete done')
    if com is not None:
    	#print('delete com')
    	deleteCom(com)

#==============================================================================
# adaptCellsNodal (deprecated) : Adapts a polyhedral mesh a1 with repsect to the nodal subdivision values.
# IN: t : 3D NGON mesh
# IN: nodal_vals : nb of subdivision required expressed at mesh nodes
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCellsNodal(t, sensdata=None, smoothing_type = 0, subdiv_type=0, hmesh=None, sensor=None):
    """Adapts an unstructured mesh a with respect to a sensor (DEPRECATED : use adaptCells with a sensor of type 2 instead.)
    Usage: adaptCellsNodal(t, sensdata=None, smoothing_type = 0, subdiv_type=0, hmesh=None, sensor=None)"""
    tp = Internal.copyRef(t)
    _adaptCells(tp, sensdata, 2, smoothing_type, -1, subdiv_type, hmesh, sensor)
    return tp

#==============================================================================
# createHMesh : Returns a hierarchical zone hook 
# IN: z : 3D NGON zone
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical zone hook 
#==============================================================================
def createHMesh(z, subdiv_type= 0):
    """Returns a hierarchical zone hook .
    Usage: createHMesh(z, subdiv_type= 0)"""
    m = C.getFields(Internal.__GridCoordinates__, z)[0]
    bcptlists = getBCPtList(z)

    hmeshs = []
    hmeshs.append(intersector.createHMesh(m, subdiv_type, bcptlists, 0, None, None, None))
    return hmeshs

#==============================================================================
# createHZones : Returns a hierarchical PyTree hook 
# IN: t : 3D NGON PyTree
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical PyTree hook 
#==============================================================================
def createHZones(t, subdiv_type= 0):
    """Returns a hierarchical PyTree hook.
    Usage: createHZones(t, subdiv_type= 0)"""
    zones = Internal.getZones(t)
    hmeshs = []

    nbz = len(zones)

    if nbz == 1:
      hmeshs = createHMesh(zones[0], subdiv_type)
      return (hmeshs, None)

    # zone name to id
    name2id = dict()
    i = 0
    for z in zones :
      name2id[z[0]] = i
      #print('zone ' + z[0] + ' has ' + str(i))
      i += 1

    # create COM
    m = C.getFields(Internal.__GridCoordinates__, zones[0])[0] # fixme : first zone tells for all
    com = createCom(m, subdiv_type, nbz)

    # create HZones
    zid=0
    for z in zones:
      (jzids, jptlists) = getJoinsPtList(z, name2id)
      bcptlists = getBCPtList(z)

      m = C.getFields(Internal.__GridCoordinates__, z)[0]
      hmeshs.append(intersector.createHMesh(m, subdiv_type, bcptlists, zid, jzids, jptlists, com))
      zid+=1
    
    return (hmeshs, com)

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
def createSensor(hmeshs, sensor_type = 0, smoothing_type=0 , itermax = -1):
    """Returns a sensor hook.
    Usage: createSensor(hmeshs, sensor_type = 0, smoothing_type=0 , itermax = -1)"""
    sensors = []
    for hmesh in hmeshs:
      if hmesh == None : 
        sensors.append(None)
      else :
        sensors.append(intersector.createSensor(hmesh,sensor_type,smoothing_type,itermax))
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

  elif sens_data_typ == 1: # list of numpies
    if len(sensdata) != len(hooks):
      print('assignData2Sensor (nodal or centered) ERROR : data list must be sized as number of sensors')
      return 1
    i=-1
    for h in hooks:
      i+=1
      intersector.assignData2Sensor(h, sensdata[i])
  return 0

def createCom(t, subdiv_type = 0, nbz = 1):
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
def conformizeHMesh(t, hooks):
    """Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
    Usage: conformizeHMesh(t, hooks)"""
    tp = Internal.copyRef(t)
    _conformizeHMesh(tp, hooks)
    return tp

#==============================================================================
# _conformizeHMesh : Converts the basic element leaves of a hierarchical mesh (hooks is a list of hooks to hiearchical zones) to a conformal polyhedral mesh.
#                   Each hiearchcial zone is referring to a zone in the original mesh t. So the mesh is replaced in the tree and the BCs/Joins/Fields are transferred.
# IN: t : PyTree before adaptation
# IN: hook : list of hooks to hiearchical zones (same size as nb of zones in t).
# OUT: Nothing 
#==============================================================================
def _conformizeHMesh(t, hooks):
    """Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
    Usage: _conformizeHMesh(t, hooks)"""
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

        fieldsC = C.getFields(Internal.__FlowSolutionCenters__, z)[0]
        if fieldsC == [] : fieldsC = None

        fieldsN = C.getFields(Internal.__FlowSolutionNodes__, z)[0]
        if fieldsN == [] : fieldsN = None

        fieldsF = None
        # todo : get fields from BCDataSets

        res = intersector.conformizeHMesh(hooks[i], fieldsC, fieldsN, fieldsF)

        # res[0] : mesh
        # res[1] : ranges for what is in res from 3 to end in res
        #          res[1][0] -> res[1][1]-1 : joins
        #          res[1][1] -> res[1][2]-1 : bcs
        #          res[1][2] -> res[1][3]-1 : fields 
        # res[2] : joined zones ids
        # res[3] -> end :  joins, bcs, fields

        mesh = res[0]

        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        if len(res) < 2 : continue

        ranges = res[1]
        jzids = res[2]
        
        # MAJ Joins
        jptlists = res[ranges[0] : ranges[1]]
        if jptlists != []:
          updateJoinsPointLists2(z, zones, jzids, jptlists)

        # MAJ BCs
        bcptlists = res[ranges[1] : ranges[2]]
        #print(bcptlists)
        if bcptlists != [] :
          updateBCPointLists2(z, bcptlists)
        else:
          C._deleteZoneBC__(z)

        C._deleteFlowSolutions__(z)

        ## MAJ center fields
        fieldz = res[ranges[2]: ranges[3]]
        for f in fieldz:
          C.setFields([f], z, 'centers', False)

        ## MAJ node fields
        fieldz = res[ranges[3]:ranges[4]]
        #print (fieldz)
        for f in fieldz:
          C.setFields([f], z, 'nodes', False)

        ## MAJ face fields 
        fieldz = res[ranges[4]:]
        # todo (in BCDataSets)
        
        i=i+1

#==============================================================================
# adaptBox : Adapts a bounding box to a cloud of interior points
#==============================================================================
def adaptBox(t, box_ratio = 10., smoothing_type=0, itermax=-1):
    """Adapts a bounding box to a cloud of interior points.
    Usage: adaptBox(t, box_ratio)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
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
# detectIdentitcalCells : detects (and optionally removes) geometrically identical cells 
#======================================================================
def detectIdenticalCells(t, TOL=1.e-15, clean=0):
    """Detects (and optionannly removes) geometrically identical cells.
    Usage: detectIdenticalCells(t)"""
    return C.TZA(t, 'nodes', 'nodes', XOR.detectIdenticalCells, t, TOL, clean)

def _detectIdenticalCells(t, TOL=1.e-15, clean=0):
    return C._TZA(t, 'nodes', 'nodes', XOR.detectIdenticalCells, t, TOL, clean)

#==============================================================================
# detectOverConnectedFaces : detects Faces that belong to more than 2 cells in a mesh.
#======================================================================
def detectOverConnectedFaces(t):
    """Detects Faces that belong to more than 2 cells in a mesh."""
    return C.TZA(t, 'nodes', 'nodes', XOR.detectOverConnectedFaces, t)

def _detectOverConnectedFaces(t, TOL=1.e-15, clean=0):
    return C._TZA(t, 'nodes', 'nodes', XOR.detectOverConnectedFaces, t, TOL, clean)

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
def getOverlappingFaces(t1, t2, RTOL = 0.1, amax = 0.1, dir2=(0.,0.,0.)):
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
# getCollidingCells     : returns the list of cells in a1 and a2 that are colliding.
# IN : t1:              : NGON mesh (surface or volume).
# IN : t2:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# OUT: 2 lists of colliding cells, the first one for t1, the seoncd one for t2.
#==============================================================================
def getCollidingCells(t1, t2, RTOL = 1.e-12, only_externals = False):
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
def getCells(t1, ids, are_face_ids = True):
   """ Returns the cells in t1 having faces or cell ids.
   Usage: getCells(t1, pgids, are_face_ids = True)"""

   zones1 = Internal.getZones(t1)
   cells = []

   i=-1
   for z in zones1:
    i+=1
    m = C.getFields(Internal.__GridCoordinates__, z)[0]
    cells .append(C.convertArrays2ZoneNode('cell%s'%i, [XOR.getCells(m, ids[i], are_face_ids)]))

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
def estimateAdapReq(t, donnor, metric_policy=2, rtol= 1.e-12, minv=0, maxv=5):
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
def getAnisoInnerFaces(t1, aniso_ratio = 0.05):
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
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    return XOR.checkCellsClosure(m)

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
        print('vmin for this zone : '+str(res[1]))
        if res[1] < vmin:
          vmin=res[1]
          ivmin=res[0]
          vzoneid=i
        print('grmin for this zone : '+str(res[3]))
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
def extrudeBC(t, height = 0.25, mean_or_min = 1, create_ghost = 1, bndType = 'UserDefined'):
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
def extrudeSurf(t, layer_height, nlayers = 1, strategy = 1):
    """XXX"""
    return C.TZA(t, 'nodes', 'nodes', XOR.extrudeSurf, t, layer_height, nlayers, strategy)

def _extrudeSurf(t, layer_height, nlayers = 1, strategy = 1):
    return C._TZA(t, 'nodes', 'nodes', XOR.extrudeSurf, t, layer_height, nlayers, strategy)

#==============================================================================
# extrudeRevolSurf : XXX
#==============================================================================
def extrudeRevolSurf(t, ax_pt, ax_dir, nlayers = 1):
    """XXX"""
    return C.TZA(t, 'nodes', 'nodes', XOR.extrudeRevolSurf, t, ax_pt, ax_dir, nlayers)

def _extrudeRevolSurf(t, ax_pt, ax_dir, nlayers = 1):
    return C._TZA(t, 'nodes', 'nodes', XOR.extrudeRevolSurf, t, ax_pt, ax_dir, nlayers)

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
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.removeBaffles(m)
    return C.convertArrays2ZoneNode('wobaffle', [m])

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

def merge(tz, sz, tol = 1.e-15): #target zone, list source zones
    m = C.getFields(Internal.__GridCoordinates__, tz)[0]
    #print m
    s = C.getFields(Internal.__GridCoordinates__, sz)[0]
    #print s
    return XOR.merge(m, s, tol)

def concatenate(t, tol = 1.e-15):
  zones = Internal.getZones(t)
  ms = []
  for z in zones:
    m = C.getFields(Internal.__GridCoordinates__, z)[0]
    ms.append(m)
  m = XOR.concatenate(ms, tol)
  return C.convertArrays2ZoneNode('assembly', [m])

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

    try : import Generator.PyTree as G
    except: raise ImportError("volume: requires Generator module.")

    fldname = None
    if fieldname != None: fldname = 'centers:'+fieldname

    zones = Internal.getZones(t)
    v = 0.
    for z in zones:
        xcelln = None
        if fldname != None : xcelln = C.getField(fldname, z)
        if xcelln != None : xcelln = xcelln[0]
        
        z = C.convertArray2NGon(z); z = G.close(z)
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        v += XOR.volume(m, xcelln)
    return v
  
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
