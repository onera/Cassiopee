"""Grid generation module.
"""
# 
# Python Interface to create PyTrees defining meshes
#
from . import Intersector as XOR
from . import intersector

__version__ = XOR.__version__

import numpy

try: range = xrange
except: pass

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter
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
  import Converter.Internal as I
  ncellsTot = 0
  zones = I.getNodesFromType2(a, 'Zone_t')
  for z in zones:
      dim = I.getZoneDim(z)
      np = dim[1]
      ncells = dim[2]
      ncellsTot += ncells
  return ncellsTot

def getTreeDim(t):
  zs = Internal.getZones(t)
  d = 0
  for z in zs:
    dims = Internal.getZoneDim(z)
    if d == 0: d = dims[4]
    if dims[4] != d:return 0; # mixed type : not handled
  return d

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

# update BC and JOINS point lists givzn an indirection "new id to old id"
def updatePointLists(z, zones, oids):
    bnds = Internal.getNodesFromType(z, 'BC_t')
    joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
    zname=z[0]

    ptLists = []
    for bb in bnds :
      ptl = Internal.getNodesFromType(bb, 'IndexArray_t')
      ptLists.append(ptl[0][1][0])
    for j in joins:
      ptl = Internal.getNodeFromName1(j, 'PointList')
      ptLists.append(ptl[1])#ptl[0][1][0]

    if (ptLists == []) : return

    # recalcul des pointlist
    ptLists = XOR.updatePointLists(oids, ptLists)

    i=0
    # update the BC pointlists 
    for bb in bnds :
      ptl = Internal.getNodesFromType(bb, 'IndexArray_t')
      ptl[0][1] = ptLists[i]
      #print(ptl[0][1])
      i=i+1

    # update the Join pointlist and synchronize with other zones (their pointListDonnor)
    for j in joins:
      donnorName = "".join(j[1])
      ptl = Internal.getNodeFromName1(j, 'PointList')
      #print(donnorName)
      dz = Internal.getNodeFromName(zones, donnorName)
      joinsD = Internal.getNodesFromType(dz, 'GridConnectivity_t')
      for jd in joinsD:
        dname = "".join(jd[1])
        if (dname != zname) : continue
        ptlD = Internal.getNodeFromName1(jd, 'PointListDonor')
        
        PG0 = ptl[1][0][0] # first polygon in the poitn list 
        PG0D = ptlD[1][0][0] # first polygon in the poitn list
        if (PG0 != PG0D) : continue # not the right join (in case of multiple joins for 2 zones) : the first PG must be the same (assume one PG only in one join)
        ptl[1]= ptLists[i]
        ptlD[1] = ptLists[i]
        break
      i=i+1

def transferFields(z, oids):
  nparrs = []
  fields = C.getFields(Internal.__FlowSolutionCenters__, z)

  nf = len(fields)
  ofields = []

  for f in fields:
    nparrs.append(f[1][0])

  C._deleteFlowSolutions__(z, 'centers') # destroy and rebuild one by one ??

  fi=0
  for f in fields:
    C._initVars(z, 'centers:'+f[0], 1.)
    of = C.getField('centers:'+f[0], z)
    
    n = len(oids) # nb of element has decreased : len <= len(fields)
    print(len(nparrs[fi]))
    #onumparr = numpy.empty((n,), dtype=numpy.int32)
    for i in range(n):
      of[0][1][0][i] = nparrs[fi][oids[i]]
    fi += 1

    #ofields.append(f)

  #print(ofields)
  #return ofields

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

def booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, multi_zone=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the union between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanUnion(a1, a2, tol)
    Usage for volumes: booleanUnion(a1, a2, tol, preserve_right, solid_right)"""

    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]

    if multi_zone == True :
        typzone1 = s1[3]
        typzone2 = s2[3]
        if typzone1 == 'NGON' and typzone2 == 'NGON': # only for Volume/Volume
          # compute the join tolerance
          L1 = edgeLengthExtrema(a1)
          L2 = edgeLengthExtrema(a2)
          jtol = 0.1*min(L1,L2)

          return booleanUnionMZ(a1, a2, tol, jtol, agg_mode, improve_qual)

    #multi_zone option is ignored from here

    cur_shift=0
    extrudepgs=[]
    if (solid_right == 1) :
        zones = Internal.getZones(a2)
        (extrudepgs, cur_shift) = concatenateBC('UserDefined', zones, extrudepgs, cur_shift)
    if (extrudepgs != []) : extrudepgs = numpy.concatenate(extrudepgs) # create a single list
    #print("nb of pgs to pass : %s" %(len(extrudepgs)))

    res = XOR.booleanUnion(s1, s2, tol, preserve_right, solid_right, agg_mode, improve_qual, extrudepgs)
    
    is_zone_list  = 0
    if (len(res) != 4) : is_zone_list = 1
    elif (res[3] != 'NGON' and res[3] != 'TRI' and res[3] != 'BAR') : is_zone_list = 1

    if (is_zone_list == 0) : return C.convertArrays2ZoneNode('union', [res])

    # debug : mutli zones
    ozones = []

    for i in range(len(res)):
      if len(res[i][0][1]) != 0: ozones.append(C.convertArrays2ZoneNode(res[i][1], [res[i][0]])) #(zname, array)
    
    return ozones

def booleanUnionMZ(t1, t2, xtol=0., jtol=0., agg_mode=1, improve_qual = False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    tp1 = Internal.copyRef(t1)
    tp2 = Internal.copyRef(t2)
    return _booleanUnionMZ(tp1, tp2, xtol, jtol, agg_mode, improve_qual)


def _booleanUnionMZ(t1, t2, xtol=0., jtol=0., agg_mode=1, improve_qual = False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
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

    res = XOR.booleanUnionMZ(m1s, m2s, xtol, jtol, agg_mode, improve_qual)

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

def diffSurf(a1, a2, tol=0., preserve_right=1, agg_mode=1, improve_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between a volume mesh and a surface mesh."""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.diffSurf(s1, s2, tol, preserve_right, agg_mode,improve_qual)
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
    Usage : XcellN(t, priorities [, output_type, rtol])"""
    tp = Internal.copyRef(t)
    _XcellN(tp, priorities, output_type, rtol)
    return tp

#==============================================================================
# _XcellN (in-place version)
# IN: t: 3D NGON mesh
# IN : priorities : one-to-one priorities between components
# IN : output_type : 0 : binary mask; 1 : continuous mask (xcelln) ; 2 : clipped surface.
# OUT: returns a 3D NGON surface mesh with the xcelln field (if output_type=0/1, the clipped surface with solution if output_type=2)
#==============================================================================
def _XcellN(t, priorities, output_type=0, rtol=0.05):
    """Computes the weight coefficients of visibility for overset grid configurations as a field called xcelln, for both surface and volume mesh of any kind.
    Usage : _XcellN(t, priorities [, output_type, rtol])"""
    d = getTreeDim(t)
    if d == 2:
      _XcellNSurf(t, priorities, output_type, rtol)
    elif d == 3:
      print ('XcellN : not implemented yet for 3D')
    else :
      print ('XcellN : the input file contain mixed 2D/3D zones : not handled currently')

#==============================================================================
# XcellNSurf
# IN: t: 3D NGON SURFACE mesh
# IN : priorities : one-to-one priorities between components
# IN : output_type : 0 : binary mask; 1 : continuous mask (xcelln) ; 2 : clipped surface.
# OUT: returns a 3D NGON surface mesh with the xcelln field (if output_type=0/1, the clipped surface with solution if output_type=2)
#==============================================================================
def XcellNSurf(t, priorities, output_type=0, rtol=0.05):
    """Computes the weight coefficients of visibility for overset grid configurations as a field called xcelln, for any kind of surface mesh.
    Usage : XcellNSurf(t, priorities [, output_type, rtol])"""
    tp = Internal.copyRef(t)
    _XcellNSurf(tp, priorities, output_type, rtol)
    return tp

#==============================================================================
# _XcellNSurf (in-place version)
# IN: t: 3D NGON SURFACE mesh
# IN : priorities : one-to-one priorities between components
# IN : output_type : 0 : binary mask; 1 : continuous mask (xcelln) ; 2 : clipped surface. 
# OUT: returns a 3D NGON surface mesh with the xcelln field (if output_type=0/1, the clipped surface with solution if output_type=2)
#==============================================================================
def _XcellNSurf(t, priorities, output_type=0, rtol=0.05):
  """Computes the weight coefficients of visibility for overset grid configurations as a field called xcelln, for any kind of surface mesh.
  Usage : _XcellNSurf(t, priorities [, rtol])"""

  try: import Transform.PyTree as T
  except: raise ImportError("XcellN: requires Transform module.")
  try: import Post.PyTree as P
  except: raise ImportError("XcellN: requires Post module.")
  try: import Generator.PyTree as G
  except: raise ImportError("XcellN: requires Generator module.")


  allbcs = Internal.KNOWNBCS
  allbcs.remove('BCOverlap')
  allbcs.remove('BCMatch')
  allbcs.remove('BCWallViscous')
  allbcs.remove('BCWallInviscid')
  allbcs.remove('BCWallViscousIsothermal')

  #t = T.reorderAll(t, dir=1)
  #C.convertPyTree2File(t, 'reorederedt.cgns')
  bases = Internal.getBases(t)
   #(BCs,BCNames,BCTypes) = C.getBCs(t)

  if len(bases) == 1 :
    print('Only one base in the file. Each component must be separated in a given Base. No check between zones of the same component.')
    return

  boundaries = []
  wall_ids = []
  
  # PREPARE INPUTS FOR ALL ZONES : boundaries and wall ids
  for b in bases:
    bj = T.join(b)

    b_bounds = P.exteriorFaces(bj)
    b_bounds = C.convertArray2Tetra(b_bounds) # to have BARs
    m_bounds = C.getFields(Internal.__GridCoordinates__, b_bounds)[0]

    walls = []
    for btype in allbcs:
      walls += C.extractBCOfType(b, btype)
    wallf = None
    if len(walls) is not 0: 
        walls = T.join(walls)
        hook = C.createHook(b_bounds, function='elementCenters') # identifying edges
        wallf = C.identifyElements(hook, walls) # wallf are ids in boundaries
        wallf -= 1 # make it 0 based
    boundaries.append(m_bounds)
    wall_ids.append(wallf)

  # get the zones in a single list with parent base id
  ngons = []
  basenum = []

  structured_tree=True
  ngon_tree=True
  zs = Internal.getZones(t)
  for z in zs:
    dims = Internal.getZoneDim(z)
    if dims[0] == 'Unstructured' : 
      structured_tree = False
    if dims[3] != 'NGON' : 
      ngon_tree = False
  #print('ngon_tree?',ngon_tree, 'structured_tree?', structured_tree )

  if ngon_tree == False:
    tNG = C.convertArray2NGon(t)
    #tNG = G.close(tNG, tol=1.e-15)
    bases = Internal.getBases(tNG)

  base_id=-1
  for b in bases:

    base_id += 1
    zones = Internal.getZones(b)
    
    for z in zones:
        z = convertNGON2DToNGON3D(z)
        c = C.getFields(Internal.__GridCoordinates__, z)[0]
        ngons.append(c)
        basenum.append(base_id)

  #print(wall_ids)

  # COMPUTE THE COEFFS PER ZONE (PARALLEL OMP PER ZONE)
  xcellns = XOR.XcellNSurf(ngons, basenum, boundaries, wall_ids, priorities, output_type, rtol)

  if output_type == 2: # output as ckipped NGON
    i=0
    C._deleteFlowSolutions__(t)#, 'nodes') #no histo for nodes currently
    bases = Internal.getBases(t)
    for b in bases:
      zones = Internal.getZones(b)
      for z in zones:
        # udpating mesh and solution
        mesh = xcellns[i]
        pg_oids = xcellns[i+1]
        
        #transferFields(z, pg_oids)
        C.setFields([mesh], z, 'nodes')

        i +=2
    return

  #APPLY IT TO ZONES
  if structured_tree == False:
    bases = Internal.getBases(t)

  i=0
  for b in bases:
    zones = Internal.getZones(b)
    for z in zones:
      C.setFields([xcellns[i]], z, 'centers', False)
      i = i+1

  if structured_tree == True :
    # back to STRUCT
    print('XcellN : back to original mesh type (STRUCT, Basic...)')
    mc = C.node2Center(tNG)
    hookC = C.createGlobalHook([mc], 'nodes')
    hookN = C.createGlobalHook([tNG], 'nodes')

    C._identifySolutions(t, tNG, hookN, hookC, tol=1000.)
    C.freeHook(hookC)
    C.freeHook(hookN)

  #C.convertPyTree2File(tNG, 'tNG.cgns')
  #C.convertPyTree2File(t, 't.cgns')

#==============================================================================
# unify
# IN: t: background Mesh (NGON 3D)
# IN: priorities: XXX
# OUT: XXX
#==============================================================================
def unify(t, priorities):
     
    tp = Internal.copyRef(t)
    _unify(tp, priorities)
    return tp

def _unify(t, priorities):

    try: import Transform.PyTree as T
    except: raise ImportError("unify: requires Transform module.")
    try: import Post.PyTree as P
    except: raise ImportError("unify: requires Post module.")
    try: import Converter.Internal as I
    except: raise ImportError("unify: requires Post module.")
    try: import Converter.PyTree as C
    except: raise ImportError("unify: requires Converter module.")
    try: import Generator.PyTree as G
    except: raise ImportError("unify: requires Generator module.")

    # create the masks and flag the walls
    print('unify : prepare masks...')

    # get walls in STRUCT from STRUCt tree (because convertArray2Hexa does not preserve walls)
    walls = []
    bs = I.getBases(t)
    for b in bs:
        # walls (BEFORE NGON CONVERSION)
        ws = extractBaseWalls(b) # QUAD
        #print ws
        if (ws != []):
          ws = C.convertArray2NGon(ws)
          ws = T.join(ws)
          ws = G.close(ws) #ficme : necessary ?
          ws = convertNGON2DToNGON3D(ws)

        walls.append(ws)

    # build mask and associate walls
    masks = []
    basewallfaces = []
    thx8 = C.convertArray2Hexa(t) #FIXME : enable if structured
    C._deleteFlowSolutions__(thx8)
    bs = I.getBases(thx8)
    ib=0
    #ms = []
    for b in bs:
     
        # mask ORIENTED OUTWARD
        b = T.join(b); b = G.close(b)
        b = C.convertArray2NGon(b)
        b = reorientExternalFaces(b)
        b = P.exteriorFaces(b)
        m = convertNGON2DToNGON3D(b)
        #m = simplifyCells(m, 1, 5.e-4)
        #C.convertPyTree2File(m, "m.cgns")
                
        c = C.getFields(Internal.__GridCoordinates__, m)[0]
        masks.append(c)
        #ms.append(m)

        wallf = []
        if (walls[ib] != []):
          hook = C.createHook(m, function='faceCenters')
          wallf = C.identifyFaces(hook, walls[ib]) # wallf are ids in m
          #print(wallf)
        basewallfaces.append(wallf)
        ib = ib+1

    #C.convertPyTree2File(ms, "masks.cgns")
    
    # compute
    print('unify : compute celln...')
    C._initVars(t, 'centers:xcelln', 1.)

    tNG = C.convertArray2NGon(t)
    tNG = G.close(tNG, tol=1.e-15)

    bases = I.getBases(tNG)
    i=0;
    for b in bases:
 
        print("unify :    for component %d over %d ..."%(i+1, len(bases)))

        zones = I.getZones(b)
        ngons = []
        basenum = []
        for z in zones:
            c = C.getFields(Internal.__GridCoordinates__, z)[0]
            ngons.append(c)
            basenum.append(i)

        xcellns = intersector.unify(ngons, basenum, masks, priorities, basewallfaces)
        C.setFields(xcellns, b, 'centers', False)
        i = i+1
    
    
    # back to STRUCT
    print('unify : back to STRUCT...')
    mc = C.node2Center(tNG)
    hookC = C.createGlobalHook([mc], 'nodes')
    hookN = C.createGlobalHook([tNG], 'nodes')

    C._identifySolutions(t, tNG, hookN, hookC, tol=1000.)
    C.freeHook(hookC)
    C.freeHook(hookN)

    #C.convertPyTree2File(tNG, 'tNG.cgns')
    #C.convertPyTree2File(t, 't.cgns')

def extractBaseWalls(b):
    try: import Transform.PyTree as T
    except: raise ImportError("unify: requires Transform module.")
    try: import Converter.Internal as I
    except: raise ImportError("unify: requires Post module.")
    try: import Converter.PyTree as C
    except: raise ImportError("unify: requires Converter module.")
    walls = []
    zbc = C.extractBCOfType(b, 'BCWall')
    for zb in zbc:
      if (zb == []) : continue 
      walls.append(zb)

    zbc = C.extractBCOfType(b, 'BCWallInviscid')
    for zb in zbc:
      if (zb == []) : continue 
      walls.append(zb)

    zbc = C.extractBCOfType(b, 'BCWallViscous')
    for zb in zbc:
      if (zb == []) : continue 
      walls.append(zb)

    return walls

# def getWallsCentroids(t):
#     try: import Transform.PyTree as T
#     except: raise ImportError("unify: requires Transform module.")
#     try: import Converter.Internal as I
#     except: raise ImportError("unify: requires Post module.")
#     try: import Converter.PyTree as C
#     except: raise ImportError("unify: requires Converter module.")
#     cloud = []
#     bs = I.getBases(t)
#     for b in bs:
#         zbc = C.extractBCOfType(b, 'BCWall')
#         if (zbc == []) : continue 
#         for zb in zbc:
#             if (zb == []) : continue 
#             zb = C.node2Center(zb)
#             cloud.append(zb)
#         zbc = C.extractBCOfType(b, 'BCWallInviscid')
#         if (zbc == []) : continue 
#         for zb in zbc:
#             if (zb == []) : continue 
#             zb = C.node2Center(zb)
#             cloud.append(zb)
#         zbc = C.extractBCOfType(b, 'BCWallViscous')
#         if (zbc == []) : continue 
#         for zb in zbc:
#             if (zb == []) : continue 
#             zb = C.node2Center(zb)
#             cloud.append(zb)
#     #cloud = C.convertArray2Node(cloud)
#     #cloud = T.join(cloud)
#     return cloud
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
# reorientExternalFaces : reorients outward the external polygons of a mesh
#==============================================================================
def reorientExternalFaces(t):
    """Reorients outward the external polygons of a mesh.
    Usage: reorientExternalFaces(t)"""
    return C.TZA(t, 'nodes', 'nodes', XOR.reorientExternalFaces, t)

def _reorientExternalFaces(t):
    return C._TZA(t, 'nodes', 'nodes', XOR.reorientExternalFaces, t)

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
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def simplifyCells(t, treat_externals, angular_threshold = 1.e-12):
    """Simplifies over-defined polyhedral cells (agglomerate some elligible polygons).
    Usage: simplifyCells(t, treat_externals, angular_threshold)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.simplifyCells(m, treat_externals, angular_threshold)
    return C.convertArrays2ZoneNode('simplifiedCells', [m])

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
def agglomerateSmallCells(t, vmin=0., vratio=1000.):
    """Agglomerates prescribed cells.
    Usage: agglomerateSmallCells(t, vmin, vratio)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]

    res = XOR.agglomerateSmallCells(m, vmin, vratio)
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
def agglomerateNonStarCells(t):
    """Agglomerates non-centroid-star-shaped cells.
    Usage: agglomerateNonStarCells(t)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]

    res = XOR.agglomerateNonStarCells(m)
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
def agglomerateCellsWithSpecifiedFaces(t, pgs, simplify=2, amax = 1.e-12): # 0 : dno not simplify, 1 : simplify only internals, 2 : simlplify evrywhere
     
    tp = Internal.copyRef(t)
    _agglomerateCellsWithSpecifiedFaces(tp,pgs, simplify, amax)
    return tp

def _agglomerateCellsWithSpecifiedFaces(t, pgs, simplify=2, amax = 1.e-12):

    zones = Internal.getZones(t)
    if len(pgs) != len(zones):
    	print('agglomerateCellsWithSpecifiedFaces: input error: nb of polygons packs differ from nb of zones : %s versus %s.'%(len(pgs), len(zones)))
    	return None

    if simplify < 0 : simplify = 0
    if simplify > 2 : simplify = 2
    i=0
    for z in zones:
      m = C.getFields(Internal.__GridCoordinates__, z)[0]
      m = XOR.agglomerateCellsWithSpecifiedFaces(m, pgs[i])
      
      simp = simplify
      if (simp != 0) : 
        simp -= 1
        m = XOR.simplifyCells(m, treat_externals=simp, angular_threshold=amax)# treat externals iff simplify==1
      C.setFields([m], z, 'nodes') # replace the mesh in the zone
      i = i+1
 
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
# removeNonManifoldExternalCells : removes any outer cell that has a non manifold edge
#==============================================================================
def removeNonManifoldExternalCells(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.removeNonManifoldExternalCells(m)
    return C.convertArrays2ZoneNode('manifoldOuter', [m])

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
# adaptCells : Adapts a polyhedral mesh a1 with repsect to a2 points
# IN: t1 : 3D NGON mesh
# IN: t2 : source points (any kind of mesh)
# IN: sensor_type : basic (0) or xsensor (1)
# IN smoothing_type : First-neighborhood (0) Shell-neighborhood(1)
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCells(t, t2, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None):
     
    tp = Internal.copyRef(t)
    _adaptCells(tp, t2, sensor_type, smoothing_type, itermax, subdiv_type, hmesh)
    return tp

def _adaptCells(t, t2, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None):
    """Adapts a polyhedral mesh t1 with repsect to t2 points.
    Usage: adaptCells(t1, t2, sensor_type, smoothing_type, itermax, subdiv_type, hmesh)"""

    source = C.getFields(Internal.__GridCoordinates__, t2)[0]

    zones = Internal.getZones(t)

    i=-1
    for z in zones:
        i+=1
        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue

        coords = Converter.convertArray2NGon(coords)

        if hmesh is not None:
            res = intersector.adaptCells(coords, source, sensor_type, smoothing_type, itermax, subdiv_type, hmesh[i])
        else:
            res = intersector.adaptCells(coords, source, sensor_type, smoothing_type, itermax, subdiv_type, None)

        mesh = res[0]
        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')
        if (len(res) > 1):
            pg_oids=res[1]
            # MAJ POINT LISTS #
            updatePointLists(z, zones, pg_oids)

    return t


#==============================================================================
# Dynamic cells adaptation 
#==============================================================================
def adaptCellsDyn(t,ts,hmeshs,hsensors):
    tp = Internal.copyRef(t)
    _adaptCellsDyn(tp, ts, hmeshs, hsensors)
    return tp

def _adaptCellsDyn(t,ts,hmeshs,hsensors):
    source = C.getFields(Internal.__GridCoordinates__, ts)[0]

    zones   = Internal.getZones(t)
    i       = 0
    hmesh   = hmeshs[i]
    hsensor = hsensors[i]
    for z in zones:        
        res = intersector.adaptCellsDyn(source,hmesh,hsensor)
        i   = i+1
    

#==============================================================================
# adaptCellsNodal : Adapts a polyhedral mesh a1 with repsect to the nodal subdivision values.
# IN: t : 3D NGON mesh
# IN: nodal_vals : nb of subdivision required expressed at mesh nodes
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCellsNodal(t, nodal_vals, hmesh=None):
    tp = Internal.copyRef(t)
    _adaptCellsNodal(tp, nodal_vals, hmesh)
    return tp

def _adaptCellsNodal(t, nodal_vals, hmesh=None):
    """Adapts a polyhedral mesh a1 with repsect to the nodal subdivision values.
    Usage: _adaptCellsNodal(t, nodal_vals)"""

    zones = Internal.getZones(t)
    nb_zones = len(zones)
    nb_nodals = len(nodal_vals)

    if nb_zones != nb_nodals:
        print('must give one nodal list (sized as cooridnates) per zone')
        return

    i=-1
    for z in zones:
        i+=1
        coords = C.getFields(Internal.__GridCoordinates__, z)[0]
        if coords == []: continue

        nval = nodal_vals[i] 
        #print(nval)

        if hmesh is not None:
            res = intersector.adaptCellsNodal(coords, nval, hmesh[i])
        else:
            res = intersector.adaptCellsNodal(coords, nval, None)

        mesh = res[0]
        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        if (len(res) > 1):
            pg_oids=res[1]
            # MAJ POINT LISTS #
            updatePointLists(z, zones, pg_oids)

    return t

def createHMesh(t, subdiv_type = 0):
    zones = Internal.getZones(t)
    i=0
    hmeshs = []
    for z in zones:
        m = C.getFields(Internal.__GridCoordinates__, z)[0]
        if m == []: continue
        hmeshs.append(intersector.createHMesh(m, subdiv_type))
        i=i+1
    return hmeshs

def createGeomSensor(hmeshs, sensor_type = 0, itermax = -1):

    sensors = []

    for hmesh in hmeshs:
        sensors.append(intersector.createSensor(hmesh,sensor_type,itermax))
  
    return sensors

def deleteHMesh(hooks):
    nb_hooks = len(hooks)
    i=0
    for h in range(nb_hooks):
        intersector.deleteHMesh(hooks[i])
        i=i+1

def conformizeHMesh(t, hooks):
    tp = Internal.copyRef(t)
    _conformizeHMesh(tp, hooks)
    return tp

def _conformizeHMesh(t, hooks):
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
        res = intersector.conformizeHMesh(hooks[i])
        mesh = res[0]
        # MAJ du maillage de la zone
        C.setFields([mesh], z, 'nodes')

        if (len(res) > 1):
            pg_oids=res[1]
            # MAJ POINT LISTS #
            updatePointLists(z, zones, pg_oids)
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
# extractOuterLayers : Extracts prescribed outer cell layers
# IN: t:               : 3D NGON mesh
# IN: N:               : Number of layers to extract
# IN: discard_external : for volume mesh with holes (e.g. external flow), set it to 1 to extract only layers around bodies.
# OUT: r
#==============================================================================
def extractOuterLayers(t, N, discard_external=0):
    """ Extracts prescribed outer cell layers.
    Usage: extractOuterLayers(t, N, discard_external)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    res = XOR.extractOuterLayers(m, N, discard_external)
    
    zones = []
        
    nb_zones = len(res)
    if (nb_zones == 1) :
      zones.append(C.convertArrays2ZoneNode('remaining', [res[0]]))
    else:
      zones.append(C.convertArrays2ZoneNode('outers', [res[0]]))
      zones.append(C.convertArrays2ZoneNode('remaining', [res[1]]))

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
   m2 = concatenate(zones2)
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
    m1 = C.getFields(Internal.__GridCoordinates__, t1)[0]
    m2 = C.getFields(Internal.__GridCoordinates__, t2)[0]
    
    res = XOR.diffMesh(m1, m2)
    
    zones = []
    nb_zones = len(res)

    if (nb_zones == 0) : 
        print("No difference.") ; return zones
    
    zones.append(C.convertArrays2ZoneNode('z1', [res[0]]))
    zones.append(C.convertArrays2ZoneNode('z2', [res[1]]))
    
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
# checkCellsClosure : XXX
#==============================================================================
def checkCellsClosure(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    return XOR.checkCellsClosure(m)

#==============================================================================
# checkForDegenCells : XXX
#==============================================================================
def checkForDegenCells(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    return XOR.checkForDegenCells(m)

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
# computeAspectRatio : Returns a field of aspect ratio
# IN: t    : 3D NGON mesh
# IN: vmim : volume threshold
# OUT: Returns the first cell id that is non-closed
#==============================================================================
def computeAspectRatio(t, vmin=0.):
    """ Returns a field of aspect ratio.
    Usage: computeAspectRatio(t, vmin)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    ar = XOR.computeAspectRatio(m, vmin)
    z = C.convertArrays2ZoneNode('w_aspect', [m])
    C.setFields([ar], z, 'centers', False)
    return z

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
# oneZonePerCell : XXX
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

def centroids(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
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

def concatenate(zones, tol = 1.e-15):
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