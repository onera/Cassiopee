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

#=============================================================================
# Concatenation des PointList d un type de BC donne dans une liste de zones
#=============================================================================
def concatenateBC(bctype, zones, wallpgs, cur_shift):
    i=0
    for z in zones:
      c = C.getFields(Internal.__GridCoordinates__, z)

      if (c == []): continue

      #print(' -- zone : %d / %d' %(i+1, len(zones)))
      i=i+1
      bnds = Internal.getNodesFromType(z, 'BC_t')
      #print(" -- this zone has %d boundaries"%(len(bnds)))
      #print(' -- cur shift %d' %(cur_shift))

      # GET THE WALL PGS FROM THE POINTLISTS
      for bb in bnds:
        if Internal.isValue(bb, bctype) == False: continue
          
        wpgs = bb[2][1][1][0] # POINTLIST NUMPY
        #print wpgs
        # SYNC THE POINTLIST BEFORE APPENDING  : SHIFT WITH THE CURRENT NB OF STORED POLYGONS
        id2 = numpy.empty(len(wpgs), numpy.int32)
        id2[:] = wpgs[:] + cur_shift
        wallpgs.append(id2)

      c = c[0]
      #print c
      #z_nb_pts= len(c[1][0])
      z_nb_pgs= c[2][0][0]
      #print z_nb_pts
      #print z_nb_pgs
      cur_shift += z_nb_pgs
    return (wallpgs, cur_shift)

# update BC and JOINS point lists givzn an indirection "new id to old id"
def updatePointLists(z, zones, oids):
    bnds = Internal.getNodesFromType(z, 'BC_t')
    joins = Internal.getNodesFromType(z, 'GridConnectivity_t')
    zname=z[0]

    ptLists = []
    for bb in bnds :
      ptLists.append(bb[2][1][1][0])
    for j in joins:
      ptLists.append(j[2][1][1][0])

    if (ptLists == []) : return

    # recalcul des pointlist
    ptLists = XOR.updatePointLists(oids, ptLists)

    i=0
    # update the BC pointlists 
    for bb in bnds :
      bb[2][1][1] = ptLists[i]
      #print bb[2][1][1]
      i=i+1

    # update the Join pointlist and synchronize with other zones (their pointListDonnor)
    for j in joins:
      donnorName = "".join(j[1])
      #print donnorName
      dz = Internal.getNodeFromName(zones, donnorName)
      joinsD = Internal.getNodesFromType(dz, 'GridConnectivity_t')
      for jd in joinsD:
        dname = "".join(jd[1])
        if (dname != zname) : continue
        
        PG0 = j[2][1][1][0][0] # first polygon in the poitn list 
        PG0D = jd[2][2][1][0][0] # first polygon in the poitn list
        if (PG0 != PG0D) : continue # not the right join (in case of multiple joins for 2 zones) : the first PG must be the same (assume one PG only in one join)
        j[2][1][1]= ptLists[i]
        jd[2][2][1] = ptLists[i]
        break
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

def booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the intersection between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanIntersection(a1, a2, tol)
    Usage for volumes: booleanIntersection(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanIntersection(s1, s2, tol, preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual)
    return C.convertArrays2ZoneNode('inter', [s])

def booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the union between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanUnion(a1, a2, tol)
    Usage for volumes: booleanUnion(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]

    cur_shift=0
    extrudepgs=[]
    if (solid_right == 1) :
        zones = Internal.getZones(a2)
        (extrudepgs, cur_shift) = concatenateBC('UserDefined', zones, extrudepgs, cur_shift)
    if (extrudepgs != []) : extrudepgs = numpy.concatenate(extrudepgs) # create a single list
    #print "nb of pgs to pass : %s" %(len(extrudepgs))

    res = XOR.booleanUnion(s1, s2, tol, preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, extrudepgs)
    
    is_zone_list  = 0
    if (len(res) != 4) : is_zone_list = 1
    elif (res[3] != 'NGON' and res[3] != 'TRI' and res[3] != 'BAR') : is_zone_list = 1

    if (is_zone_list == 0) : return C.convertArrays2ZoneNode('union', [res])

    # debug : mutli zones
    ozones = []

    for i in range(len(res)):
      if len(res[i][0][1]) != 0: ozones.append(C.convertArrays2ZoneNode(res[i][1], [res[i][0]])) #(zname, array)
    
    return ozones

def booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanMinus(a1, a2, tol)
    Usage for volumes: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanMinus(s1, s2, tol, preserve_right, solid_right, agg_mode)
    return C.convertArrays2ZoneNode('minus', [s])

def diffSurf(a1, a2, tol=0., preserve_right=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between a volume mesh and a surface mesh."""
    s1 = C.getFields(Internal.__GridCoordinates__, a1)[0]
    s2 = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.diffSurf(s1, s2, tol, preserve_right, agg_mode,improve_conformal_cloud_qual)
    return C.convertArrays2ZoneNode('VmS', [s])
    
def booleanModifiedSolid(solid, a2, tol=0., preserve_solid=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the transformed input solid after solving the intersection of its skin with a2.
    Usage: booleanModifiedSolid(a1, a2, tol, preserve_right, solid_right)"""
    sld = C.getFields(Internal.__GridCoordinates__, solid)[0]
    operand = C.getFields(Internal.__GridCoordinates__, a2)[0]
    s = XOR.booleanModifiedSolid(operand, sld, tol, preserve_solid, agg_mode, improve_conformal_cloud_qual)
    return C.convertArrays2ZoneNode('modified_solid', [s])

#==============================================================================
# XcellN
# IN: t: background Mesh (NGON 3D)
# IN: prioritaryMesh: hiding Mesh (NGON 3D)
# IN: blankingMatrix
# OUT: returns the cellnfields, between 0 (fully hidden) and 1 (fully visible)
#==============================================================================
def XcellN(t, prioritaryMesh, blankingMatrix=[]):
    try: import Transform as T
    except: raise ImportError("XcellN: requires Transform module.")

    nb = -1
    a = Internal.copyRef(t)
    # ajout du celln aux centres si n'existe pas pour une zone
    loc = 'centers'
    a = addVar__(a, var='cellN', loc=loc)
    bases = Internal.getBases(a)
    if blankingMatrix == []: blankingMatrix = numpy.ones((len(bases), len(prioritaryMesh)), numpy.int32)
    for b in bases:
        nb += 1
        #print 'bgm base : %d / %d' %(nb+1, len(bases))
        coords = C.getFields(Internal.__GridCoordinates__, b)
        if coords == []: continue

        coords = Converter.convertArray2NGon(coords)

        if loc == 'centers': cellN = C.getField('centers:cellN', b)
        else: cellN = C.getField('cellN', b)

        bc = []
        wallpgs = [] # LIST OF BODY WALLS IDS USED TO IGNORE BGM CELLS INSIDE THEM 
        ghostpgs = [] # LIST OF BOUNDARIES TO EXTRUDE TO PREVENT UNECESSARY X COMPUTATIONS
        cur_shift=0
        for nb2 in range(len(prioritaryMesh)):
            blanking = blankingMatrix[nb, nb2]
            #if (prioritaryMesh[nb2] == []): print('empty')
            if (prioritaryMesh[nb2] == []): continue
            
            #print('hiding base : %d / %d' %(nb2+1, len(prioritaryMesh)))
            zones = Internal.getZones(prioritaryMesh[nb2])
            i=0
            for z in zones:
                c = C.getFields(Internal.__GridCoordinates__, z)

                if c == []: continue

                #print(' -- hiding base %d zone : %d / %d' %(nb2+1, i+1, len(zones)))

                c = c[0]
                bc.append(c)

            (wallpgs, cur_shift_new) = concatenateBC('BCWall', zones, wallpgs, cur_shift)

            (ghostpgs, cur_shift_new) = concatenateBC('UserDefined', zones, ghostpgs, cur_shift)
            cur_shift=cur_shift_new

        if (wallpgs != []) : wallpgs = numpy.concatenate(wallpgs) # create a single list
        #print("nb of wall pgs %s"%(len(wallpgs)))
        if (ghostpgs != []) : ghostpgs = numpy.concatenate(ghostpgs) # create a single list
        #print("nb of ghost pgs %s"%(len(ghostpgs)))
                
        if bc == []:
            #print('Warning : no xcelln to compute for base %d'%(nb))
            continue

        bc = Converter.convertArray2NGon(bc); bc = T.join(bc);

        cellN = XOR.XcellN(coords, cellN, bc, wallpgs, ghostpgs)
        bc = None
        coords = None
        C.setFields(cellN, b, loc, False)
    return a

#==============================================================================
# triangulateExteriorFaces
# IN: mesh: 3D NGON mesh
# IN : in_or_out : 0 means "ONLY INTERNALS", 1 means "ONLY EXTERNALS", any other value means "BOTH"
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def triangulateExteriorFaces(t, in_or_out=2):
    """Triangulates exterior polygons of a volume mesh.
    Usage: triangulateExteriorFaces(t)"""
    return C.TZA(t, 'nodes', 'nodes', XOR.triangulateExteriorFaces, t, in_or_out)

def _triangulateExteriorFaces(t, in_or_out=2):
    return C._TZA(t, 'nodes', 'nodes', XOR.triangulateExteriorFaces, t, in_or_out)


#==============================================================================
# triangulateSpecifiedFaces 
# IN: a: 3D NGON mesh
# IN: pgs : list of polygons
# OUT: returns a 3D NGON Mesh
#==============================================================================
def triangulateSpecifiedFaces(t, pgs):
     
    tp = Internal.copyRef(t)
    _triangulateSpecifiedFaces(tp,pgs)
    return tp

def _triangulateSpecifiedFaces(t, pgs):

    zones = Internal.getZones(t)
    if (len(pgs) != len(zones)) :
        print('triangulateSpecifiedFaces: input error: nb of polygons packs differ from nb of zones.')
        return None

    i=0
    for z in zones:
      m = C.getFields(Internal.__GridCoordinates__, z)[0]
      m = XOR.triangulateSpecifiedFaces(m, pgs[i])
      C.setFields([m], z, 'nodes') # replace the mesh in the zone
      i = i+1

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
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.reorientExternalFaces(m)
    return C.convertArrays2ZoneNode('oriented', [m])
    
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
#     print "one"
#     res = XOR.agglomerateSmallCells(m, vmin, vratio)
#     print "NB ZONES %d"%(len(res))

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
def agglomerateCellsWithSpecifiedFaces(t, pgs, simplify=2): # 0 : dno not simplify, 1 : simplify only internals, 2 : simlplify evrywhere
     
    tp = Internal.copyRef(t)
    _agglomerateCellsWithSpecifiedFaces(tp,pgs, simplify)
    return tp

def _agglomerateCellsWithSpecifiedFaces(t, pgs, simplify=2):

    zones = Internal.getZones(t)
    if len(pgs) != len(zones):
    	print('agglomerateCellsWithSpecifiedFaces: input error: nb of polygons packs differ from nb of zones.')
    	return None

    i=0
    for z in zones:
      m = C.getFields(Internal.__GridCoordinates__, z)[0]
      m = XOR.agglomerateCellsWithSpecifiedFaces(m, pgs[i])
      if (simplify != 0) : 
        simplify = simplify -1
        m = XOR.simplifyCells(m, simplify)# treat externals iff simplify==1
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
# closeOctalCells : Closes any polyhedral cell in an octree
# IN: t: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all cells closed
#==============================================================================
def closeOctalCells(t):
    """Closes any polyhedral cell in an octree.
    Usage: closeOctalCells(t)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.closeOctalCells(m)
    return C.convertArrays2ZoneNode('closed', [m])

#==============================================================================
# adaptCells : Adapts a polyhedral mesh a1 with repsect to a2 points
# IN: t1 : 3D NGON mesh
# IN: t2 : source points (any kind of mesh)
# IN: sensor_type : basic (0) or xsensor (1)
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCells(t1, t2, sensor_type = 0, itermax=-1, force_basic=0):
    """Adapts a polyhedral mesh t1 with repsect to t2 points.
    Usage: adaptCells(t1, t2, sensor_type)"""
    m1 = C.getFields(Internal.__GridCoordinates__, t1)[0]
    m2 = C.getFields(Internal.__GridCoordinates__, t2)[0]
    m = intersector.adaptCells(m1, m2, sensor_type,itermax, force_basic)
    return C.convertArrays2ZoneNode('adapted', [m])

#==============================================================================
# adaptBox : Adapts a bounding box to a cloud of interior points
#==============================================================================
def adaptBox(t, box_ratio = 10., itermax=-1):
    """Adapts a bounding box to a cloud of interior points.
    Usage: adaptBox(t, box_ratio)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = intersector.adaptBox(m, box_ratio, itermax)
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
# getOverlappingFaces   : returns the list of polygons in a1 and a2 that are overlapping.
# IN : t1:              : NGON mesh (surface or volume).
# IN : t2:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# IN: ps_min            : minimal value for the dot product of the normals of each pair of colliding polygons. A value of 1. means pure parallelism.
# IN: dir2              : if specified, direction vector used for all a2's polygons instead of their own normals.
# OUT: 2 lists of overlapping polygons, the first one for a1, the seoncd one for a2.
#==============================================================================
def getOverlappingFaces(t1, t2, RTOL = 0.1, ps_min = 0.95, dir2=(0.,0.,0.)):
   """ Returns the list of polygons in a1 and a2 that are overlapping.
   Usage: getOverlappingFaces(t1, t2, RTOL, ps_min, dir2)"""

   try: import Transform as T
   except: raise ImportError("getOverlappingFaces: requires Transform module.")

   zones1 = Internal.getZones(t1)
   zones2 = Internal.getZones(t2)

   pgids = []
   t2j = T.join(zones2)
   m2 = C.getFields(Internal.__GridCoordinates__, t2j)[0]

   for z in zones1:
        
     m1 = C.getFields(Internal.__GridCoordinates__, z)[0]
     if m1 == []: continue

     pgids.append(XOR.getOverlappingFaces(m1,m2, RTOL, ps_min, dir2))

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
# edgeLengthExtrema : XXX
#==============================================================================
def edgeLengthExtrema(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    return XOR.edgeLengthExtrema(m)

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
def extrudeUserDefinedBC(t, height = 0.25, mean_or_min = 1, create_ghost = 1):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    cur_shift=0
    extrudepgs=[]
    zones = Internal.getZones(t)
    #print("nb of zones %d"%(len(zones)))
    (extrudepgs, cur_shift) = concatenateBC('UserDefined', [zones], extrudepgs, cur_shift)
    if (extrudepgs != []) : extrudepgs = numpy.concatenate(extrudepgs) # create a single list
    #print("nb of pgs to pass : %s" %(len(extrudepgs)))

    mo = XOR.extrudeUserDefinedBC(m, extrudepgs, height, mean_or_min, create_ghost)

    return C.convertArrays2ZoneNode('union', [mo])

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

    if (nb_zones == 0) : return zones

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
    """ Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
    Usage: convertNGON2DToNGON3D(t)"""
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    m = XOR.convertNGON2DToNGON3D(m)
    return C.convertArrays2ZoneNode('nuga', [m])

def centroids(t):
    m = C.getFields(Internal.__GridCoordinates__, t)[0]
    c = XOR.centroids(m)
    return C.convertArrays2ZoneNode('centroids', [c])

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
