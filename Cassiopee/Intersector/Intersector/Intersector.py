"""Intersections module.
"""
__version__ = '4.0'
__author__ = "Sam Landier"
#
# Python Interface to create arrays defining meshes
#
from . import intersector
import Generator as G
import Converter as C

try: range = xrange
except: pass

def isSingleZone(a):
    #print(a)
    if len(a) != 4: return False
    if a[3] != 'Zone_t': return False
    return True

def isCGNSTree(a):
    #print(a)
    if len(a) != 4: return False
    if a[3] != 'Zone_t': return False
    return True

def updatePointLists(oids, pointLists):
    return intersector.updatePointLists(oids, pointLists)

def conformUnstr(a1, a2=None, tol=0., left_or_right=0, itermax=10):
    """Conformizes a1 (optionally with a2).
    If a2 is specified, the third argument "left_or_right" tells wheter the ouput contains only a1 modified (0), a2 modified (1) or both (2).
    Usage: conformUnstr(a1, a2, tol, left_or_right)"""
    c = intersector.conformUnstr(a1, a2, tol, left_or_right, itermax)
    return G.close(c)

def intersection(a1, a2, tol=0., itermax=10):
    """Computes the intersection trace (a polyline) between two input closed surfaces.
    Usage: intersection(a1, a2, tol)"""
    try:
        import Converter
        a1 = Converter.convertArray2Tetra(a1)
        a2 = Converter.convertArray2Tetra(a2)
        a1 = G.close(a1); a2 = G.close(a2)
    except: pass
    c = intersector.booleanIntersectionBorder(a1, a2, tol, itermax)
    return G.close(c)

def booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
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
    c = intersector.booleanIntersection(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_qual, True, itermax)# last is dummy (outward)
    return G.close(c)

def booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, extrude_pgs=[], simplify_pgs=True, hard_mode=0, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
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
        c = intersector.booleanUnion(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_qual, extrude_pgs, 0, 0, itermax)
        return G.close(c[0])
    else:
        c = intersector.booleanUnion(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_qual, extrude_pgs, simplify_pgs, hard_mode, itermax)
        return c[0] #close is done inside

def booleanUnionWithHisto(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, extrude_pgs=[], simplify_pgs=True, hard_mode=0, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
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
        c = intersector.booleanUnion(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_qual, extrude_pgs, 0, 0, itermax)
        return G.close(c)
    else:
        c = intersector.booleanUnion(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_qual, extrude_pgs, simplify_pgs, hard_mode, itermax)
        return c #close is done inside






def booleanUnionMZ(a1, a2, xtol=0., jtol=0., agg_mode=1, improve_qual=False, simplify_pgs=True, hard_mode=0): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the union between two volume meshes.
    Usage for volumes: booleanUnionMZ(a1, a2, tol, agg_mode)"""
    c = intersector.booleanUnionMZ(a1, a2, xtol, jtol, agg_mode, improve_qual, simplify_pgs, hard_mode)
    return c #close is done inside

def booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, improve_qual=False, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
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
    c = intersector.booleanMinus(a1, a2, tol, preserve_right, solid_right, agg_mode, improve_qual, True, itermax)# last is dummy (outward)
    return G.close(c)

def booleanModifiedSolid(solid, a2, tol=0., preserve_solid=1, agg_mode=1, improve_qual=False, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the transformed input solid after solving the intersection of its skin with a2.
    Usage: booleanMinus(a1, a2, tol, preserve_right, solid_right)"""
    c = intersector.booleanModifiedSolid(a2, solid, tol, 1, preserve_solid, agg_mode, improve_qual, True, itermax)# last is dummy (outward)
    return G.close(c)

def diffSurf(a1, a2, tol=0., preserve_right=1, agg_mode=1, improve_qual=False, outward_surf=True, itermax=10): #agg_mode : 0(NONE), 1(CONVEX), 2(FULL)
    """Computes the difference between a volume mesh and a surface mesh."""
    solid_right=1
    c = intersector.DiffSurf(a1, a2, tol, solid_right, preserve_right, agg_mode, improve_qual, outward_surf, itermax)
    return G.close(c)

#==============================================================================
# XcellN
# IN: t: 3D NGON SURFACE mesh
# IN: priorities: one-to-one priorities between components
# IN: output_type: 0: binary mask; 1: continuous mask (xcelln); 2: clipped surface.
#      If set to False, the field has any value in [0,1] upon exit, the values in between are the surface ratio of the visible cells
# OUT: returns a 3D NGON surface mesh with the xcelln field (if output_type=0/1, the clipped surface with solution if output_type=2)
#==============================================================================
def XcellN(coord, zwall_ids, basenum, masks, wall_ids, priorities, output_type=0, rtol=0.05):
    """Computes the weight coefficients of visibility for overset grid configurations as a field called xcelln, for any kind of surface mesh.
    Usage : XcellN(t, priorities [, output_type, rtol])"""
    return intersector.XcellN(coord, zwall_ids, basenum, masks, wall_ids, priorities, output_type, rtol)

#==============================================================================
# P1ConservativeInterpolation
# IN: aR : receiver mesh
# IN : aD : donnor mesh
# IN : fldD : donnor center fields
# OUT: returns aR with transferred center fields from aD in a conservative manner
#==============================================================================
def P1ConservativeInterpolation(aR, aD, fldD):
    """Does conservative interpolations from aD center fields (fldD) to aR
    Usage : interpolate(aR, aD, fldD)"""
    return intersector.P1ConservativeInterpolation(aR, aD, fldD)

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
# superMesh
# IN: surfz: 3D NGON surface mesh to clip
# IN: sclip: 3D NGON surface mesh (clipper)
# IN: tol: tolerance (abolute if positive, relative otherwise)
# IN: proj_on_first: if True(False), each sclip(surfz) face is projected on surfz(sclip).
# OUT: returns the polyclipping of surfz by sclip
#==============================================================================
def superMesh(surfz, sclip, tol=-1.e-4, proj_on_first=True):
    """Polyclips surfz with sclip.
    Usage: superMesh(surfz, sclip, tol, proj_on_first)"""
    return intersector.superMesh(surfz, sclip, tol, proj_on_first)

#==============================================================================
# triangulateExteriorFaces
# IN: coords: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all the external faces triangulated
#==============================================================================
def triangulateExteriorFaces(a, in_or_out=2, improve_qual=0):
    """Triangulates exterior polygons of a volume mesh.
    Usage: triangulateExteriorFaces(a, in_or_out)"""
    return intersector.triangulateExteriorFaces(a, in_or_out, improve_qual)

#==============================================================================
# triangulateSpecifiedFaces
# IN: coords: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with BC polygons triangulated
#==============================================================================
def triangulateSpecifiedFaces(a, pgs, improve_qual=1):
    """Triangulates specified polygons of a volume mesh.
    Usage: triangulateSpecifiedFaces(a, pgs, improve_qual)"""
    return intersector.triangulateSpecifiedFaces(a, pgs, improve_qual)
#synonym
def triangulateBC(a, pgs, improve_qual=1):
    """Triangulates specified polygons of a volume mesh.
    Usage: triangulateBC(a, pgs, improve_qual=1)"""
    return intersector.triangulateSpecifiedFaces(a, pgs, improve_qual)

#==============================================================================
# triangulateNFaces
# IN: mesh: 3D NGON mesh
# IN : quality improvement flag
# OUT: returns a 3D NGON mesh with all the external faces triangulated
#==============================================================================
def triangulateNFaces(a, improve_qual=1, min_nvertices=5, discard_joins=True):
    """Triangulates nob basic polygons of a volume mesh.
    Usage: triangulateNFaces(t, improve_qual=1, min_nvertices=5, discard_joins=True)"""
    return intersector.triangulateNFaces(a, improve_qual, min_nvertices, discard_joins)

#==============================================================================
# externalFaces : Returns erternal faces for CASSIOPEE NGON types and NUGA NGON
#==============================================================================
def externalFaces(a, discarded_ids=None, geo_dim=-1):
    """Returns erternal faces for CASSIOPEE NGON types and NUGA NGON.
    Usage: externalFaces(t)"""
    return intersector.externalFaces(a, discarded_ids, geo_dim)

#==============================================================================
# reorient
# IN:
# OUT:
#==============================================================================
def reorient(a, dir=1):
    """Reorients outward the external polygons of a mesh.
    Usage: reorient(a)"""
    return intersector.reorient(a, dir)

#==============================================================================
# reorientSpecifiedFaces
# IN:
# OUT:
#==============================================================================
def reorientSpecifiedFaces(a, pgs, dir):
    """Reorients outward (dir = 1) or inward (dir=-1) the specified polygons of a mesh.
    Usage: reorientExternalFaces(a)"""
    return intersector.reorientSpecifiedFaces(a, pgs, dir)

#==============================================================================
# reorientSpecifiedFaces
# IN:
# OUT:
#==============================================================================
def reorientSurf(a, dir):
	return intersector.reorientSurf(a, dir)
#==============================================================================
# reorientBC
# IN:
# OUT:
#==============================================================================
def reorientBC(a, pgs, dir):
    """Reorients outward (dir = 1) or inward (dir=-1) the specified polygons of a mesh.
    Usage: reorientExternalFaces(a)"""
    return intersector.reorientSpecifiedFaces(a, pgs, dir)

#==============================================================================
# convexifyFaces
# IN: coords: 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all the faces made convex
#==============================================================================
def convexifyFaces(a, convexity_TOL=1.e-8):
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
def prepareCellsSplit(a, PH_set=1, split_policy=0, PH_conc_threshold=1./3., PH_cvx_threshold=0.05, PG_cvx_threshold=1.e-8):
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
def splitNonStarCells(a, PH_conc_threshold=1./3., PH_cvx_threshold=0.05, PG_cvx_threshold=1.e-8):
    """Splits some non-centroid-star_shaped cells.
    Usage: splitNonStarCells(a, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)"""
    return intersector.splitNonStarCells(a, PH_conc_threshold, PH_cvx_threshold, PG_cvx_threshold)

#==============================================================================
# simplifyCells : agglomerate superfluous polygons that overdefine cells
# IN: a                 : 3D NGON mesh
# IN: angular_threshold : should be as small as possible to avoid introducing degeneracies
# IN : discarded_ids : list of ids (0-based) to skip
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def simplifyCells(a, treat_externals, angular_threshold=1.e-12, discarded_ids=None):
    """Simplifies over-defined polyhedral cells (agglomerate some elligible polygons).
    Usage: simplifyCells(a, treat_externals, angular_threshold)"""
    return intersector.simplifyCells(a, treat_externals, angular_threshold, discarded_ids)

#==============================================================================
# simplifyFaces : remove superfluous nodes
# IN: a                 : 3D NGON mesh
# OUT: returns a 3D NGON Mesh with less nodes
#==============================================================================
def simplifyFaces(a):
    """Simplifies over-defined polygons.
    Usage: simplifyFaces(a)"""
    return intersector.simplifyFaces(a)

#==============================================================================
# simplifySurf : agglomerate superfluous polygons that overdefine cells
# IN: a                 : 3D NGON mesh
# IN: angular_threshold : should be as small as possible to avoid introducing degeneracies
# OUT: returns a 3D NGON Mesh with less polygons (but same shape)
#==============================================================================
def simplifySurf(a, angular_threshold=1.e-12):
    """Simplifies over-defined surfaces (agglomerate some elligible polygons).
    Usage: simplifySurf(a, angular_threshold)"""
    return intersector.simplifySurf(a, angular_threshold)

#==============================================================================
# agglomerateSmallCells : agglomerate prescribed cells
# IN: a: 3D NGON mesh
# IN: vmin : volume threshold
# IN: vratio : growth ratio threshold
# IN: angular_threshold : for simplying cells by agglomerating adjacent polygons
# IN: method = 0 (XXX)
# OUT: returns a 3D NGON Mesh with less cells and with a smoother growth ratio
#==============================================================================
def agglomerateSmallCells(a, vmin=0., vratio=0.01, angular_threshold=1.e-12, method=0):
    """Agglomerates prescribed cells.
    Usage: agglomerateSmallCells(a, vmin, vratio)"""
    return intersector.agglomerateSmallCells(a, vmin, vratio,angular_threshold, method)

#==============================================================================
# shellAgglomerateSmallCells : eradicate small cells by agglomerating all surrounding cells
# IN: a: 3D NGON mesh
# IN: vmin : volume threshold
# IN: vratio : growth ratio threshold
# OUT: returns a 3D NGON Mesh with less cells and with a smoother growth ratio
#==============================================================================
def shellAgglomerateSmallCells(a, vmin=0., vratio=1000.):
    """eradicate small cells by agglomerating all surrounding cells.
    Usage: shellAgglomerateSmallCells(a, vmin, vratio)"""
    return intersector.shellAgglomerateSmallCells(a, vmin, vratio)

#==============================================================================
# agglomerateCellsWithSpecifiedFaces : Agglomerates cells sharing specified polygons
# IN: a: 3D NGON mesh
# IN: pgs : list of polygons
# OUT: returns a 3D NGON Mesh
#==============================================================================
def agglomerateCellsWithSpecifiedFaces(a, pgs):
    """Agglomerates cells to make disappear specified polygons.
    Usage: agglomerateCellsWithSpecifiedFaces(a)"""
    return intersector.agglomerateCellsWithSpecifiedFaces(a, pgs)

#==============================================================================
# agglomerateNonStarCells : Agglomerates non-centroid-star-shaped cells
# IN: a: 3D NGON mesh
# IN: vmin : volume threshold
# IN: vratio : growth ratio threshold
# OUT: returns a 3D NGON Mesh with less cells and with a smoother growth ratio
#==============================================================================
def agglomerateNonStarCells(a, angular_threshold=1.e-12):
    """Agglomerates non-centroid-star-shaped cells.
    Usage: agglomerateNonStarCells(a)"""
    return intersector.agglomerateNonStarCells(a, angular_threshold)

#==============================================================================
# collapseUncomputableFaces : XXX
#==============================================================================
def collapseUncomputableFaces(a):
    return intersector.collapseUncomputableFaces(a)

#==============================================================================
# collapseSmallCells : XXX
#==============================================================================
def collapseSmallCells(a, vmin=0., grmin=-1.):
    return intersector.collapseSmallCells(a, vmin, grmin)

#==============================================================================
# removeNonManifoldExternalCells : removes any outer cell that has a non manifold edge
#==============================================================================
def removeNonManifoldExternalCells(a):
    return intersector.removeNonManifoldExternalCells(a)

#==============================================================================
# immerseNodes : Regularize the mesh to be cut (by boolean operations) by immersing its nodes lying on cutter surface.
# IN: a : 3D NGON mesh
# IN: s : unstructured surface mesh
# IN: TOL : tolerance (negative value means relative)
# OUT: returns a 3D NGON Mesh with moved vertices
#==============================================================================
def immerseNodes(a, s, TOL):
    """Regularize the mesh to be cut (by boolean operations) by immersing its nodes lying on the cutter surface.
    Usage: immerseNodes(a,s, TOL)"""
    return intersector.immerseNodes(a, s, TOL)

#==============================================================================
# closeCells : Closes any polyhedral cell in a mesh (processes hanging nodes on edges).
# IN: a : 3D NGON mesh
# OUT: returns a 3D NGON Mesh with all cells closed
#==============================================================================
def closeCells(a):
    """Closes any polyhedral cell in a mesh (processes hanging nodes on edges).
    Usage: closeCells(a)"""
    return intersector.closeCells([a], [0], {}, {})[0]

#==============================================================================
# adaptCells : Adapts an unstructured mesh a with respect to a sensor
# IN: a : 3D NGON unstructured mesh (with some basic elements)
# IN: sensdata : sensor data (a bunch of vertices or a mesh for a geom sensor, a mesh for a xsensor, punctual values for a nodal or cell sensor)
# IN: sensor_type : geom_sensor (0) , xsensor (1), nodal_sensor (2), cell_sensor(3)
# IN: smoothing_type : First-neighborhood (0) Shell-neighborhood(1)
# IN: itermax : max number of level in the hierarchy
# IN: subdiv_type : isotropic currently
# IN: hmesh : hierarchical mesh hook
# IN: sensor : sensor hook
# OUT: returns a 3D NGON conformal polyhedral mesh with adapted cells
#==============================================================================
def adaptCells(a, sensdata=None, sensor_type=0, smoothing_type=0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None):
    """Adapts an unstructured mesh a with repsect to a sensor.
    Usage: adaptCells(a1, [sensdata, sensor_type, smoothing_type, itermax, subdiv_type, hmesh, sensor])"""
    owesHMesh=0
    if hmesh is None:
        hmesh = createHMesh(a, subdiv_type)
        if hmesh == None : return [a] # no basic elt in a
        owesHMesh=1
    sensor = createSensor(hmesh, sensor_type, smoothing_type, itermax)

    if sensor_type == 1 or sensor_type == 4: #xsensor2 need an NGON
        sensdata = C.convertArray2NGon(sensdata)

    assignData2Sensor(sensor, sensdata)
    am = intersector.adaptCells(hmesh, sensor, {}, {})
    if owesHMesh == 1 :
        deleteHMesh(hmesh)
    deleteSensor(sensor)
    return am

#==============================================================================
# adaptBox : Adapts a bounding box to a cloud of interior points.
#==============================================================================
def adaptBox(a, box_ratio=10., smoothing_type=0, itermax=-1):
    """Adapts a bounding box to a cloud of interior points"""
    return intersector.adaptBox(a, box_ratio, smoothing_type, itermax)

#==============================================================================
# createHMesh : Returns a hierarchical zone hook
# IN: a : 3D NGON array
# IN: subdiv_type : isotropic currently
# OUT: Returns a hierarchical zone hook
#==============================================================================
def createHMesh(a, subdiv_type=0): # 0 : ISO, 1: ISO_HEX
    """Returns a hierarchical zone hook.
    Usage: createHMesh(a, subdiv_type= 0)"""
    a = intersector.initForAdaptCells(a, {})
    return intersector.createHMesh(a, subdiv_type, 0)

#==============================================================================
# deleteHMesh : Releases a hierachical zone hook
# IN: hmesh : hierarchcial mesh hook
# OUT: Nothing
#==============================================================================
def deleteHMesh(hmesh):
    """Releases a hierachical zone hook.
    Usage: deleteHMesh(hooks)"""
    return intersector.deleteHMesh(hmesh)

#==============================================================================
# conformizeHMesh : Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
#
# IN: hmesh : hierarchcial mesh hook
# OUT: Nothing
#==============================================================================
def conformizeHMesh(hmesh):
    """Converts the basic element leaves of a hierarchical mesh to a conformal polyhedral mesh.
    Usage: conformizeHMesh(hooks)"""
    return intersector.conformizeHMesh(hmesh)

def interpolateHMeshNodalField(hooks, fieldN):
    return intersector.interpolateHMeshNodalField(hooks, fieldN)

def createSensor(hmesh, sensor_type=0, smoothing_type=0 , itermax=-1, sensor_metric_policy=0):
    return intersector.createSensor(hmesh, sensor_type, smoothing_type, itermax, sensor_metric_policy)

def assignData2Sensor(hmesh, sensdata):
    return intersector.assignData2Sensor(hmesh, sensdata)

def deleteSensor(hmesh):
    return intersector.deleteSensor(hmesh)

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
# extractNthCell : Extracts nth cell in a
# IN: a          : 3D NGON mesh
# IN: nth        : cell number
# OUT: returns the cell and its neighbors
#==============================================================================
def extractNthCell(a, nth):
    """ Extracts nth cell in a.
    Usage: extractNthCell(a, nth)"""
    return intersector.extractNthCell(a, nth)

#==============================================================================
# extractBiggestCell : Extracts the biggest cell in a mesh
# IN: a          : 3D NGON mesh
# OUT: returns a single cell NGON mesh
#==============================================================================
def extractBiggestCell(a):
    """ Extracts the biggest cell in a.
    Usage: extractBiggestCell(a)"""
    return intersector.extractBiggestCell(a)

#==============================================================================
# extractBadVolCells : extract cells with bad volume / growth ratio
# IN: a          : 3D NGON mesh
# OUT: returns a single cell NGON mesh
#==============================================================================
def extractBadVolCells(a, PE, ar=0.125, vmin=0., nneighs=0):
    """ Extracts bad cells based on gowth ratio
    Usage: extractBadVolCells(a, PE, ar, vmin, nneighs)"""
    return intersector.extractBadVolCells(a, PE, ar, vmin, nneighs)

#==============================================================================
# extractOverConnectedCells : XXX
# IN: a          : 3D NGON mesh
# OUT: rXXX
#==============================================================================
def extractOverConnectedCells(a, nneighs=0):
    """ XXX
    Usage: extractOverConnectedCells(a, nneighs)"""
    return intersector.extractOverConnectedCells(a, nneighs)

#==============================================================================
# extractNthFace : Extracts the nth face in a NGON mesh
# IN: a          : 3D NGON mesh
# IN: nth        : face number
# OUT: returns the face and its left and right cells
#==============================================================================
def extractNthFace(a, nth):
    """ Extracts the nth face in a.
    Usage: extractNthFace(a, nth)"""
    return intersector.extractNthFace(a, nth)

#==============================================================================
# removeNthCell : Remove the nth cell in a mesh
# IN: a          : 3D NGON mesh
# IN: nth        : cell number
# OUT: returns the mesh without the prescribed cell
#==============================================================================
def removeNthCell(a, nth):
    """ Removes the nth cell in a.
    Usage: removeNthCell(a, nth)"""
    return intersector.removeNthCell(a, nth)

#==============================================================================
# removeNthFace : Remove the nth cell in a mesh
# IN: a          : 3D NGON mesh
# IN: nth        : cell number
# OUT: returns the mesh without the prescribed cell
#==============================================================================
def removeNthFace(a, nth):
    """ Removes the nth cell in a.
    Usage: removeNthCell(a, nth)"""
    return intersector.removeNthFace(a, nth)

#==============================================================================
# detectIdenticalCells : detects (and optionally removes) geometrically identical cells
# IN: a          : 3D NGON mesh
# IN: nth        : cell number
# OUT: returns the mesh without the prescribed cell
#==============================================================================
def detectIdenticalCells(a, TOL=1.e-15, clean=0):
    return intersector.detectIdenticalCells(a, TOL, clean)

#==============================================================================
# detectOverConnectedFaces : detects Faces that belong to more than 2 cells in a mesh.
#======================================================================
def detectOverConnectedFaces(a):
    """Detects Faces that belong to more than 2 cells in a mesh."""
    return intersector.detectOverConnectedFaces(a)

#==============================================================================
# collapseSmallEdges : XXX
#======================================================================
def collapseSmallEdges(a, eratio, lmax=-1):
    """XXX"""
    return intersector.collapseSmallEdges(a, eratio, lmax)


#==============================================================================
# getOverlappingFaces   : returns the list of polygons in a1 and a2 that are overlapping.
# IN : a1:              : NGON mesh (surface or volume).
# IN : a2:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# IN: amax              : maximal angular value (in rad) between the normals of each pair of colliding polygons.
#                         In ragnge [0,PI]. A value of 0. means pure parallelism. A value of PI means any collision.
# IN: dir2              : if specified, direction vector used for all a2's polygons instead of their own normals.
# OUT: 2 lists of overlapping polygons, the first one for a1, the seoncd one for a2.
#==============================================================================
def getOverlappingFaces(a1, a2, RTOL=0.1, amax=0.1, dir2=(0.,0.,0.)):
    """ Returns the list of polygons in a1 and a2 that are overlapping.
    Usage: getOverlappingFaces(a1, a2, RTOL, amax, dir2)"""
    return intersector.getOverlappingFaces(a1,a2, RTOL, amax, dir2)


#==============================================================================
# getCollidingTopFaces  : Returns the list of TRI/QUAD in a1 (HEXA and PRISM only) that collide a2.
# IN : a1:              : NGON mesh.
# IN : a2:              : NGON mesh.
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# OUT: list of t1 involved faces
#==============================================================================
def getCollidingTopFaces(a1, a2, RTOL=0.1):
    """ Returns the list of TRI/QUAD in a1 (HEXA and PRISM only) that collide a2.
    Usage: getCollidingTopFaces(a1, a2, RTOL)"""
    try :
        res  = intersector.getCollidingTopFaces(a1,a2, RTOL)
    except:
        res = []
    return res

#==============================================================================
# getCollidingCells     : returns the list of cells in a1 and a2 that are colliding.
# IN : t1:              : NGON mesh (surface or volume).
# IN : t2:              : NGON mesh (surface or volume).
# IN : RTOL:            : Relative tolerance (in ]0., 1.[).
# IN : only_externals:  : Boolean
# OUT: 2 lists of colliding cells, the first one for a1, the seoncd one for a2.
#==============================================================================
def getCollidingCells(a1, a2, RTOL=1.e-12, only_externals=False):
    """ Returns the list of cells in a1 and a2 that are colliding.
   Usage: getCollidingCells(a1, a2, RTOL = 1.e-12, only_externals = False)"""
    return intersector.getCollidingCells(a1, a2, RTOL, only_externals)

#==============================================================================
# getFaceIdsWithCentroids     : returns the cells in t1 attached to polygons in s2.
# IN : t1:              : NGON mesh (volume).
# IN : t2:              : NGON mesh (surface).
# OUT: attached cells.
#==============================================================================
def getFaceIdsWithCentroids(t1, cents):
    """ Returns the faces in t1 having their centroids in cents.
    Usage: getFacesWithCentroids(t1, cents)"""
    return intersector.getFaceIdsWithCentroids(t1,cents)

#==============================================================================
# getFaceIdsCollidingVertices     : returns the cells in t1 attached to polygons in s2.
# IN : t1:              : NGON mesh (volume).
# IN : cloud:              : vertex cloud.
# OUT: attached cells.
#==============================================================================
def getFaceIdsCollidingVertex(t1, vtx):
    """ Returns the faces in t1 having their centroids in cents.
    Usage: getFaceIdsCollidingVertex(t1, vtx)"""
    return intersector.getFaceIdsCollidingVertex(t1, vtx)

#==============================================================================
# getFaces              : returns the faces with ids in pgids.
# IN : t1:              : NGON mesh (volume).
# IN : pgids:           : polygon ids.
# OUT: attached cells.
#==============================================================================
def getFaces(t1, pgids):
    """ Returns the faces with ids in pgids.
    Usage: getFaces(t1, pgids)"""
    return intersector.getFaces(t1, pgids)

#==============================================================================
# getCells              : returns the cells in t1 having specified faces or cell ids.
# IN : t1:              : NGON mesh (volume).
# IN : ids:             : polygon ids.
# OUT: selected cells.
#==============================================================================
def getCells(t1, ids, are_face_ids=True):
    """ Returns the cells in t1 having specified faces or cell ids.
    Usage: getCells(t1, ids, are_face_ids)"""
    return intersector.getCells(t1, ids, are_face_ids)

#==============================================================================
# getNthNeighborhood     : returns the list of cells in the N-thneighborhood of t cells given in ids
# IN : a :               : NGON mesh.
# IN : N :               : number of neighborhood required
# IN : ids :             : input cells ids
# OUT: Returns the list of cells in the N-th neighborhood.
#==============================================================================
def getNthNeighborhood(a, N, ids):
    """ Returns the list of cells in the N-th neighborhood of cells given in ids.
    Usage: getNthNeighborhood(t, N, ids)"""
    return intersector.getNthNeighborhood(a, N, ids)

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
    """ estimates an cell-specified adaptation requirement based on a istotropic metric field based on donnor connectivity.
    Usage : estimateAdapReq(t, donnor [, metric_policy, rtol, minv, maxv])"""
    return intersector.estimateAdapReq(t, donnor, metric_policy, rtol, minv, maxv)

#==============================================================================
# getAnisoInnerFaces   : returns the list of polygons in a1 that are connecting 2 aniso elements.
# IN : a1:              : NGON mesh (surface or volume).
# IN : aniso_ratio:            : xxx

# OUT: 1 list of elected polygons
#==============================================================================
def getAnisoInnerFaces(a1, aniso_ratio=0.05):
    """ Returns the list of polygons in a1 that are connecting 2 aniso elements.
    Usage: getAnisoInnerFaces(a1, aniso_ratio)"""
    return intersector.getAnisoInnerFaces(a1,aniso_ratio)


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
    """ Returns the first cell id that is open.
    Usage: checkCellsClosure(a)"""
    return intersector.checkCellsClosure(a)

#==============================================================================
# checkCellsFlux : Computes the Gauss fluxes using the input orientation (ParentElement).
#                  Should be clsoe to zero machine for a closed and well oriented mesh.
# IN: a:               : 3D NGON mesh
# OUT: A message telling the cell id for which the Gauss flux is the greatest.
#==============================================================================
def checkCellsFlux(a, PE):
    """ Returns the cell id for which the Gauss flux is the greatest.
    Usage: checkCellsFlux(a, PE)"""
    return intersector.checkCellsFlux(a, PE)

#==============================================================================
# checkAngularExtrema : Returns the min/max dihedral angles vals and associted cell ids..
# IN: a:         : 3D NGON mesh
# OUT: A message telling the cell id for which the Gauss flux is the greatest.
#==============================================================================
def checkAngularExtrema(a, PE):
    """ Returns the min/max dihedral angles vals and associted cell ids.
    Usage: checkAngularExtrema(a, PE)"""
    return intersector.checkAngularExtrema(a, PE)
#==============================================================================
# checkCellsVolume : Computes the minimum volume using the input orientation (ParentElement).
#
# IN: a:               : 3D NGON mesh
# OUT: A message telling the cell id for which the volume is the smallest.
#==============================================================================
def checkCellsVolume(a, PE):
    """ Computes the minimum volume using the input orientation (ParentElement).
    Usage: checkCellsVolume(a, PE)"""
    return intersector.checkCellsVolume(a, PE)

#==============================================================================
# checkCellsVolumeAndGrowthRatio : Computes the minimum volume and growth ratio using the input orientation (ParentElement).
#
# IN: a:               : 3D NGON mesh
# OUT: A message telling the cell id for which the volume is the smallest.
#==============================================================================
def checkCellsVolumeAndGrowthRatio(a, PE):
    """ Computes the minimum volume and growth ratio using the input orientation (ParentElement).
    Usage: checkCellsVolumeAndGrowthRatio(a, PE)"""
    return intersector.checkCellsVolumeAndGrowthRatio(a, PE)

#==============================================================================
# checkForDegenCells : check if there are any cell with less than 4 faces.
#==============================================================================
def checkForDegenCells(a):
    """ Checks if there are any cell with less than 4 faces.
    Usage: checkForDegenCells(a)"""
    return intersector.checkForDegenCells(a)

#==============================================================================
# checkForBigCells : XXX
#==============================================================================
def checkForBigCells(a, n):
    return intersector.checkForBigCells(a,n)

#==============================================================================
# oneph : XXX
#==============================================================================
def oneph(a):
    return intersector.oneph(a)
#==============================================================================
# edgeLengthExtrema : returns min and max edge lengths
#==============================================================================
def edgeLengthExtrema(a):
    return intersector.edgeLengthExtrema(a)

def edgeLengthMax(a):
    return intersector.edgeLengthMax(a)

#==============================================================================
# computeGrowthRatio : Returns a field of growth ratio
# IN: a    : 3D NGON mesh
# IN: vmim : volume threshold
# OUT: Returns the first cell id that is non-closed
#==============================================================================
def computeGrowthRatio(a, vmin=0.):
    """ Returns a field of growth ratio.
    Usage: computeGrowthRatio(a, vmin)"""
    return intersector.computeGrowthRatio(a, vmin)



#==============================================================================
# extrudeBC : XXX
#==============================================================================
def extrudeBC(a, extrude_pgs=[], height=0.25, mean_or_min=1, create_ghost=1):
    return intersector.extrudeBC(a, extrude_pgs, height, mean_or_min, create_ghost)

#==============================================================================
# extrudeSurf : XXX
#==============================================================================
def extrudeSurf(a, layer_height, nlayers=1, strategy=1):
    return intersector.extrudeSurf(a, layer_height, nlayers, strategy)

#==============================================================================
# extrudeRevolSurf : XXX
#==============================================================================
def extrudeRevolSurf(a, ax_pt, ax_dir, nlayers=1):
    return intersector.extrudeRevolSurf(a, ax_pt, ax_dir, nlayers)

#==============================================================================
# statsUncomputableFaces : XXX
#==============================================================================
def statsUncomputableFaces(a):
	return intersector.statsUncomputableFaces(a)

#==============================================================================
# statsSize : XXX
#==============================================================================
def statsSize(a, compute_metrics=1):
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
# oneZonePerFace : XXX
#==============================================================================
def oneZonePerFace(a):
    return intersector.oneZonePerFace(a)

#==============================================================================
# convertNGON2DToNGON3D : Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
# IN: a    : 3D NGON mesh
# OUT: Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
#==============================================================================
def convertNGON2DToNGON3D(a):
    """ Converts a Cassiopee NGON Format for polygons (Face/Edge) to a Face/Node Format.
      Usage: convertNGON2DToNGON3D(a)"""
    return intersector.convertNGON2DToNGON3D(a)

#==============================================================================
# convertBasic2NGONFaces : Converts a Basic (TRI,QUAD) surface to a NGON nuga Format.
# IN: a    : Surface mesh
# OUT: Converts a Basic (TRI,QUAD) surface to a NGON nuga Format.
#==============================================================================
def convertBasic2NGONFaces(a):
    """ Converts a Basic (TRI,QUAD) surface to a NGON nuga Format.
      Usage: convertBasic2NGONFaces(a)"""
    return intersector.convertBasic2NGONFaces(a)

#==============================================================================
# centroids : Computes cells centroids in a
# IN: a    : a 3D polyhedral mesh
# OUT: a NODE aray
#==============================================================================
def centroids(a):
    """ Computes cells centroids in a.
    Usage: centroids(a)"""
    return intersector.centroids(a)

#==============================================================================
# volumes : Computes cells volumes in a
# IN: a    : a 3D polyhedral mesh
# OUT: a numpy 1D array of floats.
#==============================================================================
def volumes(a, algo=1, all_pgs_convex=False):
    """ Computes cells volumes in a.
    Usage: volumes(a, algo=1, all_pgs_convex=False)"""
    return intersector.volumes(a, algo, all_pgs_convex)

def merge(a, s, tol=1.e-15): #target arr, source arr
    return intersector.merge(a, s, tol)

def concatenate(la, tol=1.e-15):
    return intersector.concatenate(la, tol)

#==============================================================================
# volume : Computes the volume of a with an optional weighting
# IN: a  : 3D NGON mesh
# IN: xcelln : name of the weighting field (at centers)/
# OUT: volume computed
#==============================================================================
def volume(a, xcelln=None):
    """ Computes the volume of a with an optional weighting.
    Usage: volume(a, xcelln)"""
    return intersector.volume(a, xcelln)

#==============================================================================
# XXX
#==============================================================================
def drawOrientation(a):
    """XXX
    Usage: XXX"""
    return intersector.drawOrientation(a)

#==============================================================================
# syncMacthPeriodicFaces : force periodicity for faces that are supposed to be periodic
# IN: a                 : 3D NGON mesh
# IN: rotationCenter : coordinates of the center of rotation for the periodicity
# IN: rotationAngle : rotation axis for the periodicity (its norm gives the angle of rotation)
# IN : translation : translation vector for a translation periodicity
# IN : TOL : tolerance. A negative value give a relative tolerance base on min edge length
# OUT: returns a 3D NGON Mesh with synchronised faces
#==============================================================================
def syncMacthPeriodicFaces(a, rotationCenter=[0.,0.,0.],
                           rotationAngle=[0.,0.,0.],
                           translation=[0.,0.,0.], tol=-0.01):
    """ Force periodicity for faces that are supposed to be periodic.
      Usage: syncMacthPeriodicFaces(a, rotationCenter, rotationAngle, translation, TOL)"""
    return intersector.syncMacthPeriodicFaces(a, rotationCenter, rotationAngle,
                                              translation, tol)

#~ def conservativeTransfer(a1, flowsol, a2, tol=0., reconstruction_type=0):
    #~ c = intersector.conservative_transfer(a1, flowsol, a2, tol, reconstruction_type)
    #~ return c
    #~
#~ def totalMass(a1, flowsol):
    #~ intersector.total_mass(a1, flowsol)
    #~ return a1

def testmain(a):
    intersector.testmain(a)
