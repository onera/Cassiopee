"""Grid generation module.
"""
__version__ = '2.9'
__author__ = "Sam Landier, Christophe Benoit, Stephanie Peron, Luis Bernardos"
# 
# Python Interface to create arrays defining meshes
#
from . import intersector
import Generator as G
import Converter as C

try: range = xrange
except: pass


def updatePointLists(oids, pointLists):
    return intersector.updatePointLists(oids, pointLists)


def conformUnstr(a1, a2=None, tol=0., left_or_right=0, itermax=10):
    """Conformizes a1 (optionally with a2).
    If a2 is specified, the third argument "left_or_right" tells wheter the ouput contains only a1 modified (0), a2 modified (1) or both (2).
    Usage: conformUnstr(a1, a2, tol, left_or_right)"""
    c = intersector.conformUnstr(a1, a2, tol, left_or_right, itermax)
    return G.close(c)

def intersection(a1, a2, tol=0.):
    """Computes the intersection trace (a polyline) between two input closed surfaces.
    Usage: intersection(a1, a2, tol)"""
    try:
        import Converter
        a1 = Converter.convertArray2Tetra(a1)
        a2 = Converter.convertArray2Tetra(a2)
        a1 = G.close(a1); a2 = G.close(a2)
    except: pass
    c = intersector.booleanIntersectionBorder(a1, a2, tol, 1, 1, 0, False) #last 4 args are dummy for now
    return G.close(c)

def booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the intersection between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanIntersection(a1, a2, tol)
    Usage for volumes: booleanIntersection(a1, a2, tol, preserve_right, solid_right, agg_mode)"""
    if a1[3] != 'NGON' and a2[3] != 'NGON':
      try:
          import Converter
          a1 = Converter.convertArray2Tetra(a1)
          a2 = Converter.convertArray2Tetra(a2)
          a1 = G.close(a1); a2 = G.close(a2)
      except: pass
    c = intersector.booleanIntersection(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual)
    return G.close(c)

def booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_conformal_cloud_qual=False, extrude_pgs=[]): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the union between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanUnion(a1, a2, tol)
    Usage for volumes: booleanUnion(a1, a2, tol, preserve_right, solid_right)"""
    if a1[3] != 'NGON' and a2[3] != 'NGON':
      try:
        import Converter
        a1 = Converter.convertArray2Tetra(a1)
        a2 = Converter.convertArray2Tetra(a2)
        a1 = G.close(a1); a2 = G.close(a2)
      except: pass
      c = intersector.booleanUnion(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, extrude_pgs)
      return G.close(c)
    else: 
      c = intersector.booleanUnion(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual, extrude_pgs)
      return c #close is done inside

def booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between two closed-surface or two volume meshes.
    Usage for surfaces or bars: booleanMinus(a1, a2, tol)
    Usage for volumes: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    if a1[3] != 'NGON' and a2[3] != 'NGON':
      try:
        import Converter
        a1 = Converter.convertArray2Tetra(a1)
        a2 = Converter.convertArray2Tetra(a2)
        a1 = G.close(a1); a2 = G.close(a2)
      except: pass
    c = intersector.booleanMinus(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_conformal_cloud_qual)
    return G.close(c)

def booleanModifiedSolid(solid, a2, tol=0., preserve_solid=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the transformed input solid after solving the intersection of its skin with a2.
    Usage: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    c = intersector.booleanModifiedSolid(a2, solid, tol, 1, preserve_solid, agg_mode, improve_conformal_cloud_qual)
    return G.close(c)
    
def diffSurf(a1, a2, tol=0., preserve_right=1, agg_mode=1, improve_conformal_cloud_qual=False): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between a volume mesh and a surface mesh."""
    solid_right=1
    c = intersector.DiffSurf(a1, a2, tol, solid_right, preserve_right, agg_mode, improve_conformal_cloud_qual)
    return G.close(c)
    
#==============================================================================
# XcellN
# IN: coords: 3D structured or unstructured mesh
# OUT: returns the cellnfields, 0 for fully inside, 1 for fully outside, in between when intersecting
#==============================================================================
def XcellN(coords, cellnfields, maskingMesh, wall_pgl=[], ghost_pgl=[]):
    """Computes the acurate cellN based on overlapped volumes.
    Usage: XcellN(coords, cellnfields, maskingMesh)"""
    cellnt = []
    #C.convertArrays2File([maskingMesh], "mask.plt")
    #print(pgl)
    for i in range(len(coords)):
      #print 'coords : %d / %d' %(i+1, len(coords))
      #C.convertArrays2File([coords[i]], "bloc%d.plt"%(i))
      #print pgl
      #print coords[i]
      cn = intersector.XcellN(coords[i], cellnfields[i], maskingMesh, wall_pgl, ghost_pgl)
      #C.convertArrays2File([cn], "walls%d.plt"%(i))
      cellnt.append(cn)
    return cellnt
    
#==============================================================================
# P1ConservativeChimeraCoeffs
# IN: aR: receiver mesh
# IN: cellnR: receiver cellN (only cells with value equal to 2 will be considered)
# IN: aD: donor (source) mesh
# OUT: (indices, coeffs, delimiter, receiver original ids)
#==============================================================================
def P1ConservativeChimeraCoeffs(aR, cellnR, aD):
    return intersector.P1ConservativeChimeraCoeffs(aR, cellnR, aD)

#==============================================================================
# triangulateExteriorFaces
# IN: coords: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all the external faces triangulated
#==============================================================================
def triangulateExteriorFaces(a, in_or_out=2):
    """Triangulates exterior polygons of a volume mesh.
    Usage: triangulateExteriorFaces(a, in_or_out)"""
    return intersector.triangulateExteriorFaces(a, in_or_out)

#==============================================================================
# triangulateSpecifiedFaces
# IN: coords: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with BC polygons triangulated
#==============================================================================
def triangulateSpecifiedFaces(a, pgs):
    """Triangulates specified polygons of a volume mesh.
    Usage: triangulateExteriorFaces(a, in_or_out)"""
    return intersector.triangulateSpecifiedFaces(a, pgs)
#synonym
def triangulateBC(a, pgs):
    """Triangulates specified polygons of a volume mesh.
    Usage: triangulateExteriorFaces(a, in_or_out)"""
    return intersector.triangulateSpecifiedFaces(a, pgs)

#==============================================================================
# reorientExternalFaces
# IN: 
# OUT: 
#==============================================================================
def reorientExternalFaces(a):
    """Reorients outward the external polygons of a mesh.
    Usage: reorientExternalFaces(a)"""
    return intersector.reorientExternalFaces(a)
    
#==============================================================================
# convexifyFaces
# IN: coords: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all the faces made convex
#==============================================================================
def convexifyFaces(a, convexity_TOL = 1.e-8):
    """Convexifies any non-convex polygon in a mesh.
    Usage: convexifyFaces(a)"""
    return intersector.convexifyFaces(a, convexity_TOL)

#==============================================================================
# prepareCellsSplit
# IN : mesh         : 3D NGON mesh
# IN : PH_set       : PH to process. 0 for concave cells or 1 for non-centroid-star_shaped cells
# IN : split_policy : 0 : convexify concave pgs. 1 : starify concave pgs. 2 : starify any pgs at concave-chains ends.
# OUT: returns a 3D NGON Mesh with some face polygon splits : 
#      split (convexify, starify) some targeted polygons on targeted cells 
#      (typically bad polyhedra -concaves, non-centroid-star-shaped-)
#      to prepare the split of those bad cells.
#==============================================================================
def prepareCellsSplit(a, PH_set = 1, split_policy = 0, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8):
    """Splits some prescribed polygons following a prescribed splitting policy.
    Usage: prepareCellsSplit(a, PH_set, split_policy, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)"""
    return intersector.prepareCellsSplit(a, PH_set, split_policy, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)
   
#==============================================================================
# splitNonStarCells
# IN : a                 : 3D NGON mesh
# IN : PH_conc_threshold : concavity dihedral angle threshold for cells
# IN : PH_cvx_threshold  : convexity dihedral angle threshold for cells
# IN : PG_cvx_threshold  : convexity angle threshold for polygons
# OUT: returns a 3D NGON Mesh with none (or at least less) non-centroid-star_shaped cells.
#==============================================================================
def splitNonStarCells(a, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8):
    """Splits some non-centroid-star_shaped cells.
    Usage: splitNonStarCells(a, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)"""
    return intersector.splitNonStarCells(a, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)

#==============================================================================
# simplifyCells : agglomerate superfluous polygons that overdefine cells
# IN: a                 : 3D NGON mesh
# IN: angular_threshold : should be as small as possible to avoid introducing degeneracies
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def simplifyCells(a, treat_externals, angular_threshold = 1.e-12):
    """Simplifies over-defined polyhedral cells (agglomerate some elligible polygons).
    Usage: simplifyCells(a, treat_externals, angular_threshold)"""
    return intersector.simplifyCells(a, treat_externals, angular_threshold)

#==============================================================================
# agglomerateSmallCells : agglomerate prescribed cells
# IN: a: 3D NGON mesh
# IN: vmin : volume threshold
# IN: vratio : aspect ratio threshold
# OUT: returns a 3D NGON Mesh with less cells and with a smoother aspect ratio
#==============================================================================
def agglomerateSmallCells(a, vmin=0., vratio=1000.):
    """Agglomerates prescribed cells.
    Usage: agglomerateSmallCells(a, vmin, vratio)"""
    return intersector.agglomerateSmallCells(a, vmin, vratio)

#==============================================================================
# agglomerateCellsWithSpecifiedFaces : Agglomerates cells sharing specified polygons
# IN: a: 3D NGON mesh
# IN: pgs : list of polygons
# OUT: returns a 3D NGON Mesh
#==============================================================================
def agglomerateCellsWithSpecifiedFaces(a, pgs):
    """Agglomerates cells to make disappear specified polygons
    Usage: agglomerateCellsWithSpecifiedFaces(a)"""
    return intersector.agglomerateCellsWithSpecifiedFaces(a, pgs)

#==============================================================================
# agglomerateNonStarCells : Agglomerates non-centroid-star-shaped cells
# IN: a: 3D NGON mesh
# IN: vmin : volume threshold
# IN: vratio : aspect ratio threshold
# OUT: returns a 3D NGON Mesh with less cells and with a smoother aspect ratio
#==============================================================================
def agglomerateNonStarCells(a):
    """Agglomerates non-centroid-star-shaped cells.
    Usage: agglomerateNonStarCells(a)"""
    return intersector.agglomerateNonStarCells(a)

#==============================================================================
# collapseUncomputableFaces : XXX
#==============================================================================
def collapseUncomputableFaces(a):
    return intersector.collapseUncomputableFaces(a)

#==============================================================================
# removeNonManifoldExternalCells : removes any outer cell that has a non manifold edge
#==============================================================================
def removeNonManifoldExternalCells(a):
    return intersector.removeNonManifoldExternalCells(a)

#==============================================================================
# closeOctalCells : Closes any polyhedral cell in an octree
# IN: a : 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all cells closed
#==============================================================================
def closeOctalCells(a):
    """Closes any polyhedral cell in an octree.
    Usage: closeOctalCells(a)"""
    return intersector.closeOctalCells(a)

#==============================================================================
# adaptCells : Adapts a polyhedral mesh a1 with repsect to a2 points
# IN: a1 : 3D NGON mesh
# IN: a2 : source points (any kind of mesh)
# OUT: returns a 3D NGON Mesh with adapted cells
#==============================================================================
def adaptCells(a1, a2, sensor_type=0, itermax=-1, force_basic=0):
    """Adapts a polyhedral mesh a1 with repsect to a2 points.
    Usage: adaptCells(a1, a2, sensor_type)"""
    return intersector.adaptCells(a1, a2, sensor_type, itermax, force_basic)

#==============================================================================
# adaptBox : Adapts a bounding box to a cloud of interior points.
#==============================================================================
def adaptBox(a, box_ratio=10., itermax=-1):
    """Adapts a bounding box to a cloud of interior points"""
    return intersector.adaptBox(a, box_ratio, itermax)

    
#==============================================================================
# extractUncomputables : Extracts any entity that will probably cause trouble to a CFD solver
# IN: t:         : 3D NGON mesh
# IN: neigh_level: number of neighbor layers (surounding pathologies) to extract as well
# OUT: returns a 3D NGON Mesh with seprate zones for pathologies
#==============================================================================
def extractUncomputables(a, neigh_level=2):
    """Extracts any entity that will probably cause trouble to a CFD solver.
    Usage: extractUncomputables(a, neigh_level)"""
    return intersector.extractUncomputables(a, neigh_level)

#==============================================================================
# extractPathologicalCells : Extracts cells that will cause a CFD solver run failure
# IN: a          : 3D NGON mesh
# IN: neigh_level: number of neighbor layers (surounding pathologies) to extract as well
# OUT: returns a 3D NGON Mesh with seprate zones for pathologies
#==============================================================================
def extractPathologicalCells(a, neigh_level=0):
    """Extracts all cells that will probably cause trouble to a CFD solver.
    Usage: extractPathologicalCells(t, neigh_level)"""
    return intersector.extractPathologicalCells(a, neigh_level)

#==============================================================================
# extractOuterLayers : Extracts prescribed outer cell layers
# IN: a:               : 3D NGON mesh
# IN: N:               : Number of layers to extract
# IN: discard_external : for volume mesh with holes (e.g. external flow), set it to 1 to extract only layers around bodies.
# OUT: returns a 3D NGON Mesh with seprate zones for pathologies
#==============================================================================
def extractOuterLayers(a, N, discard_external=0):
    """ Extracts prescribed outer cell layers.
    Usage: extractOuterLayers(a, N, discard_external)"""
    return intersector.extractOuterLayers(a, N, discard_external)

#==============================================================================
# XXX : XXX
#==============================================================================
def extractNthCell(a, nth):
    return intersector.extractNthCell(a, nth)

#==============================================================================
# XXX : XXX
#==============================================================================
def extractNthFace(a, nth):
    return intersector.extractNthFace(a, nth)

#==============================================================================
# XXX : XXX
#==============================================================================
def removeNthCell(a, nth):
    return intersector.removeNthCell(a, nth)

#==============================================================================
# getOverlappingFaces   : returns the list of polygons in a1 and a2 that are overlapping.
# IN : a1:              : NGON mesh (surface or volume).
# IN : a2:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# IN: ps_min            : minimal value for the dot product of the normals of each pair of colliding polygons. A value of 1. means pure parallelism.
# IN: dir2              : if specified, direction vector used for all a2's polygons instead of their own normals.
# OUT: 2 lists of overlapping polygons, the first one for a1, the seoncd one for a2.
#==============================================================================
def getOverlappingFaces(a1, a2, RTOL = 0.1, ps_min = 0.95, dir2=(0.,0.,0.)):
    """ Returns the list of polygons in a1 and a2 that are overlapping.
    Usage: getOverlappingFaces(a1, a2, RTOL, ps_min, dir2)"""
    return intersector.getOverlappingFaces(a1,a2, RTOL, ps_min, dir2)

#==============================================================================
# selfX : Checks self-intersections in a mesh
# IN: a:               : 3D NGON mesh
# OUT: Returns the first two cell ids that collide.
#==============================================================================
def selfX(a):
    """ Checks  self-intersections in a mesh.
    Usage: selfX(a)"""
    return intersector.selfX(a)

#==============================================================================
# diffMesh : Returns the diff between 2 meshes as 2 zones.
# IN: a1, a2:               : two 3D NGON mesh to check
# OUT: Returns 2 zones : one with the a1 cells that are not in a2, the other is the reciprocal.
#==============================================================================
def diffMesh(a1, a2):
    """ Returns the difference between 2 meshes as 2 zones.
    Usage: diffMesh(a1, a2)"""
    return intersector.diffMesh(a1, a2)

#==============================================================================
# checkCellsClosure : Returns the first cell id that is non-closed.
# IN: a:               : 3D NGON mesh
# OUT: Returns the first cell id that is non-closed
#==============================================================================
def checkCellsClosure(a):
    """ Returns the first cell id that is non-closed.
    Usage: checkCellsClosure(a)"""
    return intersector.checkCellsClosure(a)

#==============================================================================
# checkForDegenCells : check if there are any cell with less than 4 faces.
#==============================================================================
def checkForDegenCells(a):
    return intersector.checkForDegenCells(a)

#==============================================================================
# edgeLengthExtrema : returns min and max edge lengths
#==============================================================================
def edgeLengthExtrema(a):
    return intersector.edgeLengthExtrema(a)

#==============================================================================
# computeAspectRatio : Returns a field of aspect ratio
# IN: a    : 3D NGON mesh
# IN: vmim : volume threshold
# OUT: Returns the first cell id that is non-closed
#==============================================================================
def computeAspectRatio(a, vmin=0.):
    """ Returns a field of aspect ratio.
    Usage: computeAspectRatio(a, vmin)"""
    return intersector.computeAspectRatio(a, vmin)



#==============================================================================
# extrudeUserDefinedBC : XXX
#==============================================================================
def extrudeUserDefinedBC(a, extrude_pgs=[], height = 0.25, mean_or_min = 1, create_ghost=1):
    return intersector.extrudeUserDefinedBC(a, extrude_pgs, height, mean_or_min, create_ghost)

#==============================================================================
# statsUncomputableFaces : XXX
#==============================================================================
def statsUncomputableFaces(a):
	return intersector.statsUncomputableFaces(a)

#==============================================================================
# statsSize : XXX
#==============================================================================
def statsSize(a, compute_metrics = 1):
  return intersector.statsSize(a, compute_metrics)

#==============================================================================
# removeBaffles : XXX
#==============================================================================
def removeBaffles(a):
  return intersector.removeBaffles(a)

#==============================================================================
# convert2Polyhedron : XXX
#==============================================================================
def convert2Polyhedron(a):
  return intersector.convert2Polyhedron(a)

#==============================================================================
# oneZonePerCell : XXX
#==============================================================================
def oneZonePerCell(a):
  return intersector.oneZonePerCell(a)



#==============================================================================
# convertNGON2DToNGON3D : Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
# IN: a    : 3D NGON mesh
# OUT: Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
#==============================================================================
def convertNGON2DToNGON3D(a):
  """ Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
    Usage: convertNGON2DToNGON3D(a)"""
  return intersector.convertNGON2DToNGON3D(a)

#~ def conservativeTransfer(a1, flowsol, a2, tol=0., reconstruction_type=0):
    #~ c = intersector.conservative_transfer(a1, flowsol, a2, tol, reconstruction_type)
    #~ return c
    #~ 
#~ def totalMass(a1, flowsol):
    #~ intersector.total_mass(a1, flowsol)
    #~ return a1
    
def centroids(a):
    return intersector.centroids(a)
