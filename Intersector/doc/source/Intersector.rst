.. Intersector documentation master file

Intersector: Mesh-Intersection-Based Services
=============================================

Preamble
########

This module provides pre and post processing services relying on mesh-intersection computations on arbitrary polyhedral meshes.

It also gives auxiliary functions that transform topologically and geometrically polyhedral meshes which are useful in the process of mesh generation by intersection.
    
A mesh can be stored as an array (as defined in the Converter documentation)
or in a zone node of a CGNS/python tree (pyTree).

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use with the array interface, you have to import Intersector module::

   import Intersector as XOR

For use with the pyTree interface::

    import Intersector.PyTree as XOR

.. py:module:: Intersector


List of functions
##################


**-- Main Functions**

.. autosummary::

   Intersector.conformUnstr 
   Intersector.booleanUnion
   Intersector.booleanIntersection  
   Intersector.booleanMinus
   Intersector.diffSurf
   Intersector.intersection
   Intersector.XcellN
.. Intersector.booleanModifiedSolid
.. Intersector.P1ConservativeChimeraCoeffs 

**-- Collision predicates**

.. autosummary::

   Intersector.selfX
   Intersector.getOverlappingFaces

**-- Transformation Functions**

.. autosummary::

   Intersector.triangulateBC
   Intersector.triangulateExteriorFaces
   Intersector.reorientExternalFaces
   Intersector.convexifyFaces  
   Intersector.prepareCellsSplit
   Intersector.splitNonStarCells
   Intersector.simplifyCells
   Intersector.agglomerateSmallCells
   Intersector.agglomerateNonStarCells
   Intersector.agglomerateCellsWithSpecifiedFaces
   Intersector.closeOctalCells
   Intersector.adaptCells
   Intersector.adaptBox

**-- Extraction Functions**

.. autosummary::

   Intersector.extractPathologicalCells
   Intersector.extractOuterLayers
..   Intersector.extractUncomputables
..   Intersector.extractNthCell
..   Intersector.extractNthFace
..   Intersector.removeNthCell

**-- Check Functions**

.. autosummary::

   
   Intersector.diffMesh
   Intersector.checkCellsClosure
   Intersector.computeAspectRatio
   
.. Intersector.statsUncomputableFaces
.. Intersector.statsSize
.. Intersector.edgeLengthExtrema
.. Intersector.checkForDegenCells
.. Intersector.extrudeUserDefinedBC
   
.. Intersector.removeBaffles
.. Intersector.convert2Polyhedron
.. Intersector.oneZonePerCell

**-- Conversion Functions**

.. autosummary::
   
   Intersector.convertNGON2DToNGON3D

Contents
#########

Main Functions
--------------------------


.. py:function:: Intersector.conformUnstr(a1, a2=None, tol=0., left_or_right=0, itermax=10)

    Makes conformal a TRI or a BAR soup (i.e. a set of elements not necessarily connected as a mesh) by detecting and solving all the collisions between elements. 
    
    Colliding elements are cut to get a conformal set. Mixing types BAR and TRI is not currently handled.

    :param a1:  First input mesh (BAR or TRI).
    :type  a1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2:   Second input mesh (BAR or TRI). If s2 is 'None' self-intersections are solved over s1.
    :type  a2:   [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol:   Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol:   float
    :param left_or_right: Tells the function what to ouput : the transformed s1 (0), s2(1) or both (2).
    :type  left_or_right: 0,1 or 2
    :param itermax: Number of intersection/merging iterations. 10 is the default value.
    :type  itermax: int

    *Tips and Notes:*

    * Set itermax to 1. to improve speed and the Delaunay kernel robustness. The result might have poorer quality triangles though.
    
    * Tolerance :

      - if tol > 0. : the value is used as an absolute overall tolerance
      - if tol = 0. : a value is computed as being 5% of the smallest edge length.
      - if tol < 0. : MIN(5%, -tol) is used as a ratio to apply to the smallest edge length to get the tolerance.
    
    *Example of use:*

    * `Makes conform a TRI or BAR soup (array) <Examples/Intersector/conformUnstr.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/conformUnstr.py

    * `Makes conform a TRI or BAR soup (pyTree) <Examples/Intersector/conformUnstrPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/conformUnstrPT.py


---------------------------------------


.. py:function:: Intersector.booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, extrude_pgs=[])

    Creates a conformal union between two components, either TRI surfaces or Polyhedral volumes. 

    :param a1:  First mesh operand.
    :type  a1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2:   Second mesh operand.
    :type  a2:   [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol:   Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol:   float
    :param preserve_right: Indicates the merging direction, either a1->a2 or a2->a1. If set to 1(0), it means a1->a2 (a2->a1), i.e. a2(a1)'s points are preserved.
    :type  preserve_right: 0 or 1
    :param solid_right: Indicates that the second operand is not penetrable, i.e. it is prioritized over the first operand a1.
    :type  solid_right: 0 or 1
    :param agg_mode: Option for agglomerating cut polygons : 0 to keep them as split triangles, 1 to get convex agglomerations and 2 to get a full agglomeration.
    :type  agg_mode: 0,1 or 2.
    :param extrude_pgs: Optional list of polygons to extrude.
    :type  extrude_pgs: list of int

    
    **Prerequisites :**

    * External polygons must be oriented consistently and outwardly (use Transform.reorderAll before)

    
    *Tips and Notes:*

    * For assembling meshes, set solid_right to 1 and pass the prioritized mesh as second operand.
    
    * extrude_pgs : required whenever a1 and a2 are in contact and a2 is prioritized : avoids to compute useless intersections by telling what are the indices of contact polygons in a2.


    *Example of use:*

    * `Union of two spherical surface meshes (array) <Examples/Intersector/booleanUnion.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/booleanUnion.py

    * `Union of two volume meshes (pyTree) <Examples/Intersector/booleanUnionNGPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/booleanUnionNGPT.py

---------------------------------------


.. py:function:: Intersector.booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1)

    Computes a conformal intersection between two components, either TRI surfaces or Polyhedral volumes. 

    :param a1:  First mesh operand.
    :type  a1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2:   Second mesh operand.
    :type  a2:   [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol:   Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol:   float
    :param preserve_right: Indicates the merging direction, either a1->a2 or a2->a1. If set to 1(0), it means a1->a2 (a2->a1), i.e. a2(a1)'s points are preserved.
    :type  preserve_right: 0 or 1
    :param solid_right: Indicates that the second operand is not penetrable, i.e. it is prioritized over the first operand a1.
    :type  solid_right: 0 or 1
    :param agg_mode: Option for agglomerating cut polygons : 0 to keep them as split triangles, 1 to get convex agglomerations and 2 to get a full agglomeration.
    :type  agg_mode: 0,1 or 2.
        
    **Prerequisites :**

    * External polygons must be oriented consistently and outwardly (use Transform.reorderAll before)


    *Example of use:*

    * `Intersection of two spherical surface meshes (array) <Examples/Intersector/booleanIntersection.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/booleanIntersection.py

    * `Intersection of two volume meshes (pyTree) <Examples/Intersector/booleanIntersectionNGPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/booleanIntersectionNGPT.py

---------------------------------------


.. py:function:: Intersector.booleanMinus(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1)

    Computes a conformal difference between two components, either TRI surfaces or Polyhedral volumes. 

    :param a1:  First mesh operand.
    :type  a1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2:   Second mesh operand.
    :type  a2:   [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol:   Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol:   float
    :param preserve_right: Indicates the merging direction, either a1->a2 or a2->a1. If set to 1(0), it means a1->a2 (a2->a1), i.e. a2(a1)'s points are preserved.
    :type  preserve_right: 0 or 1
    :param solid_right: Indicates that the second operand is not penetrable, i.e. it is prioritized over the first operand a1.
    :type  solid_right: 0 or 1
    :param agg_mode: Option for agglomerating cut polygons : 0 to keep them as split triangles, 1 to get convex agglomerations and 2 to get a full agglomeration.
    :type  agg_mode: 0,1 or 2.
        
    **Prerequisites :**

    * External polygons must be oriented consistently and outwardly (use Transform.reorderAll before)
    

    *Example of use:*

    * `Difference of two spherical surface meshes (array) <Examples/Intersector/booleanMinus.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/booleanMinus.py

    * `Difference of two volume meshes (pyTree) <Examples/Intersector/booleanMinusNGPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/booleanMinusNGPT.py


---------------------------------------


.. py:function:: Intersector.intersection(a1, a2, tol=0.)


    Returns the 'BAR' contour defining the intersection between two TRI-surfaces.

    :param a1:  First mesh operand.
    :type  a1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2:   Second mesh operand.
    :type  a2:   [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol:   Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol:   float

    *Example of use:*

    * `Circular trace of the intersection between two spheres (array) <Examples/Intersector/intersection.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/intersection.py

    * `Circular trace of the intersection between two spheres (pyTree) <Examples/Intersector/intersectionPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/intersectionPT.py

---------------------------------------


.. py:function:: Intersector.diffSurf(a1, a2, tol=0., preserve_right=1, agg_mode=1)

    Cut-cell function : Computes a conformal difference between a volume mesh and a surface mesh

    :param a1:  First mesh operand.
    :type  a1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2:   Second mesh operand.
    :type  a2:   [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol:   Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol:   float
    :param preserve_right: Indicates the merging direction, either a1->a2 or a2->a1. If set to 1(0), it means a1->a2 (a2->a1), i.e. a2(a1)'s points are preserved.
    :type  preserve_right: 0 or 1
    :param solid_right: Indicates that the second operand is not penetrable, i.e. it is prioritized over the first operand a1.
    :type  solid_right: 0 or 1
    :param agg_mode: Option for agglomerating cut polygons : 0 to keep them as split triangles, 1 to get convex agglomerations and 2 to get a full agglomeration.
    :type  agg_mode: 0,1 or 2.
        
    **Prerequisites :**

    * External polygons must be oriented consistently and outwardly (use Transform.reorderAll before)

    * The surface format must be an NGON Face/Node (apply before Intersector.convertNGON2DToNGON3D on the surface)
    

    *Example of use:*

    * `Cut-cell mesh with an octree and a sphere (array) <Examples/Intersector/diffSurf.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/diffSurf.py

    * `Cut-cell mesh with an octree and a sphere (pyTree) <Examples/Intersector/diffSurfPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/diffSurfPT.py

---------------------------------------


.. py:function:: Intersector.XcellN(a, cellnfields, maskingMesh, wall_pgl=[], ghost_pgl=[])

    Computes the cell nature field of a background mesh (a) in an overset configuration : similarly to the blankCells functions, the input maskingMesh are volume meshes that hide a.
    
    The computed celln is accurate, giving a floating value ranging from 0 (fully masked) to 1 (fully visible).

    The input grids (a and makingMesh) are defined by coordinates located at nodes as a list of arrays.

    :param           a:  Mesh where to compute XcellN
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param cellnfields:   celln array for a
    :type  cellnfields:   [array, list of arrays]
    :param           maskingMesh:  Prioritized mesh that hides a
    :type            maskingMesh:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param wall_pgl:   Optional list of polygons to treat as walls.
    :type  wall_pgl:   list of int
    :param ghost_pgl:   Optional list of polygons to extrude.
    :type  ghost_pgl:   list of int

    *Tips and Notes:*

    * Warning: location of celln must be located at centers.
    * Warning: In order to set the celln to 0. inside blanking bodies, you need to create BCWall type boundaries on the body faces.

    *Example of use:*

    * `Computing an accurate cellN (array) <Examples/Intersector/XcellN.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/XcellN.py

    * `Computing an accurate cellN (pyTree) <Examples/Intersector/XcellNPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/XcellNPT.py

---------------------------------------


Transformation Functions
--------------------------

.. py:function:: Intersector.triangulateBC(a, bctype)

    Triangulates the prescribed BC type polygons of a volume mesh.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param      bctype:  boundary type ('BCWall', ...).
    :type       bctype:  string

    *Example of use:*

    * `BC polygons triangulation (pyTree) <Examples/Intersector/triangulateBCPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/triangulateBCPT.py


.. py:function:: Intersector.triangulateExteriorFaces(a, in_or_out=2)

    Triangulates the prescribed external polygons of a volume mesh.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param   in_or_out:  In case of a non-connex mesh (i.e. whith holes like an external airflow mesh around bodies), set to 0 for processing only body walls, set to 1 for processing only the outer boundary, or 2 for processing all of them.
    :type    in_or_out:  0,1 or 2

    *Example of use:*

    * `External polygons triangulation (array) <Examples/Intersector/triangulateExteriorFaces.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/triangulateExteriorFaces.py

    * `External polygons triangulation (pyTree) <Examples/Intersector/triangulateExteriorFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/triangulateExteriorFacesPT.py

---------------------------------------


.. py:function:: Intersector.reorientExternalFaces(a)

    Reorients outward the external polygons of a mesh.

    *Example of use:*

    * `Reorientation (array) <Examples/Intersector/reorientExternalFaces.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/reorientExternalFaces.py

    * `Reorientation (pyTree) <Examples/Intersector/reorientExternalFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/reorientExternalFacesPT.py

---------------------------------------


.. py:function:: Intersector.convexifyFaces(a, convexity_TOL = 1.e-8)

    Makes a convex decomposition of any concave polygon in a mesh.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param   convexity_TOL:  convexity angle threshold
    :type    convexity_TOL:  float

    *Example of use:*

    * `Convexify polygons (array) <Examples/Intersector/convexifyFaces.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/convexifyFaces.py

    * `Convexify polygons (pyTree) <Examples/Intersector/convexifyFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/convexifyFacesPT.py

---------------------------------------


.. py:function:: Intersector.prepareCellsSplit(a, PH_set = 1, split_policy = 0, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8)

    Prepares the bad cells split (:any:`splitNonStarCells`) by splitting some of their polygons with a prescribed policy : convexification, starification.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param      PH_set:  PH to process. 0 for concave cells or 1 for non-centroid-star_shaped cells
    :type       PH_set:  0 or 1
    :param      split_policy:  0 : convexify concave pgs. 1 : starify concave pgs from worst vertex. 2 : starify concave pgs from concave-chains ends.
    :type       split_policy:  0,1 or 2
    :param      PH_conc_threshold:  Concavity dihedral angle threshold for cells
    :type       PH_conc_threshold:  float
    :param      PH_cvx_threshold:  Convexity dihedral angle threshold for cells
    :type       PH_cvx_threshold:  float
    :param      PG_cvx_threshold:  Convexity angle threshold for polygons
    :type       PG_cvx_threshold:  float

    *Example of use:*

    * `prepareCellsSplit (array) <Examples/Intersector/prepareCellsSplit.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/prepareCellsSplit.py

    * `prepareCellsSplit (pyTree) <Examples/Intersector/prepareCellsSplitPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/prepareCellsSplitPT.py



---------------------------------------


.. py:function:: Intersector.splitNonStarCells(a, PH_conc_threshold = 1./3., PH_cvx_threshold = 0.05, PG_cvx_threshold = 1.e-8)

    First strategy to eradicate bad cells : Splits non-centroid-star-shaped (NCSS) cells into two cells. These cells might be NCSS as well so this function should be called several times to get rid off the pathologies. Some call agglomerateSmallCells should be done afterwards to balance the aspect ratio.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param      PH_conc_threshold:  Concavity dihedral angle threshold for cells
    :type       PH_conc_threshold:  float
    :param      PH_cvx_threshold:  Convexity dihedral angle threshold for cells
    :type       PH_cvx_threshold:  float
    :param      PG_cvx_threshold:  Convexity angle threshold for polygons
    :type       PG_cvx_threshold:  float

    *Tips and Notes:*

    * Call :any:`prepareCellsSplit` before this function to ensure to process as much pathologies as possible.

    *Example of use:*

    * `splitNonStarCells (array) <Examples/Intersector/splitCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/splitCells.py

    * `splitNonStarCells (pyTree) <Examples/Intersector/splitCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/splitCellsPT.py

---------------------------------------


.. py:function:: Intersector.simplifyCells(a, treat_externals, angular_threshold = 1.e-12)

    Agglomerates superfluous polygons that over-defines cells. After agglomerating (e.g. after calling :any:`agglomerateSmallCells`) , we end up with cells that are multiply-connected, i.e. they share more than one polygon.
    If 2 cells share 2 polygons that are connected (sharing an edge) and their dihedral angle is below the angular_threshold, then the 2 polygon are agglomerated upon exit.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param      treat_externals:  Process outer polygons (1) or not (0).
    :type       treat_externals:  0 or 1
    :param      angular_threshold:  Largest angular deviation admitted between polygons in order to allow the agglomeration for them.
    :type       angular_threshold:  float

    *Example of use:*

    * `simplifyCells (array) <Examples/Intersector/simplifyCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/simplifyCells.py

    * `simplifyCells (pyTree) <Examples/Intersector/simplifyCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/simplifyCellsPT.py

---------------------------------------


.. py:function:: Intersector.agglomerateSmallCells(a, vmin=0., vratio=1000.)

    Agglomerates cells that are too small (below vmin) or having a poor aspect ratio with a neighbor (below vratio) with the best neighbor available. The agglomeration process does not create non-star-shaped agglomerates.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param      vmin:  volume threshold.
    :type       vmin:  float
    :param      vratio: aspect ratio threshold.
    :type       vratio:  float


    *Tips and Notes:*

    * See :any:`computeAspectRatio` to get the definition of the computed aspect ratio.

    *Example of use:*

    * `agglomerateSmallCells (array) <Examples/Intersector/agglomerateSmallCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/agglomerateSmallCells.py

    * `agglomerateSmallCells (pyTree) <Examples/Intersector/agglomerateSmallCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/agglomerateSmallCellsPT.py


---------------------------------------


.. py:function:: Intersector.agglomerateNonStarCells(a)

    Agglomerate cells that are non-centroid-star-shaped. The agglomeration process does not create non-star-shaped agglomerates.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

    *Example of use:*

    * `agglomerateNonStarCells (array) <Examples/Intersector/agglomerateNonStarCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/agglomerateNonStarCells.py

    * `agglomerateNonStarCells (pyTree) <Examples/Intersector/agglomerateNonStarCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/agglomerateNonStarCellsPT.py

---------------------------------------

.. py:function:: Intersector.agglomerateCellsWithSpecifiedFaces(a, pgs, simplify)

    Agglomerate cells that are non-centroid-star-shaped. The agglomeration process does not create non-star-shaped agglomerates.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           pgs:  list of polygons to remove
    :type            pgs:  list of integers


    *Example of use:*

    * `agglomerateCellsWithSpecifiedFaces (array) <Examples/Intersector/getOverlappingFaces.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/getOverlappingFaces.py

    * `agglomerateCellsWithSpecifiedFaces (pyTree) <Examples/Intersector/getOverlappingFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/getOverlappingFacesPT.py

    *Tips and Notes:*

    * When asembling 2 meshes m1 and m2 where m2 is priorized , to improve the assembly quality, do before calling the boolean union:

    1) getOverlappingFaces (m1, skin(m2)) where skin(m2) is the external polygonal skin of m2

    2) agglomerateCellsWithSpecifiedFaces on m1 with the above list of polygons
    

---------------------------------------


.. py:function:: Intersector.closeOctalCells(a)

    Closes any polyhedral cell in a 2:1 octree.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

    *Example of use:*

    * `closeOctalCells (array) <Examples/Intersector/closeOctalCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/closeOctalCells.py

    * `closeOctalCells (pyTree) <Examples/Intersector/closeOctalCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/closeOctalCellsPT.py

    *Tips and Notes:*

    * Do this transformation whenever you need to use a surface algorithm on the octree (e.g. :any:`reorientExternalFaces`)

---------------------------------------


.. py:function:: Intersector.adaptCells(a1, a2, sensor_type=0)

    Adapts a1 cells (only TETRA and HEXA cells currently) with respect to a2 points. Adaptation is a per-cell octal 2:1 decomposition.
    With a sensor_type equal to 0, a2 points are only considered : a1 will be refined such any a1 cell contains at most 1 a2's point.
    With a sensor_type equal to 1, a2's connectivity is also taken into account by adding refinement wherever a1 cells are crossed by a2 edges.

    :param           a1:  Input mesh (NGON format)
    :type            a1:  [array] or [ single zone pyTree (currently)]
    :param           a2:  Source points or source mesh
    :type            a2:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           sensor_type:  type of sensor. Using only the point cloud (0) or both points and connectivity via instersections (1)
    :type            sensor_type:  int

    *Example of use:*

    * `adaptCells (array) <Examples/Intersector/adaptCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCells.py

    * `adaptCells (pyTree) <Examples/Intersector/adaptCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCellsPT.py

    *Tips and Notes:*

    * Do this transformation before calling any Volume-Volume boolean operations in order to improve the mesh quality of the result.
    * When the input mesh has mixed-basic-type elements, only Tets and Hexas will be considered currently for adaptation. but the result wil be conformal, the pyramids and prisms will modified to respect the conformity. 

---------------------------------------


.. py:function:: Intersector.adaptBox(a, box_ratio)

    Adapts the bounding box of a cloud of points. Adaptation is an octal 2:1 decomposition.

    :param           a:  Input points cloud
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           box_ratio:  ratio to scale the box
    :type            a2:  float

    *Example of use:*

    * `adaptBox (array) <Examples/Intersector/adaptBox.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptBox.py

    * `adaptBox (pyTree) <Examples/Intersector/adaptBoxPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptBoxPT.py


---------------------------------------

Extraction Functions
--------------------------


.. py:function:: Intersector.extractPathologicalCells(a, neigh_level=0)

    Extracts cells that will potentially cause a failure when running a CFD solver. There are 4 zones upon exit, one for each pathology:

      - Non-centroid-star-shaped Cells
      - Cells having degenrated polygons for which the normal cannot be computed
      - Cells having degenrated polygons for a delaunay triangulation fails
      - Open Cells

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           neigh_level:  Number of neighbor layers (surounding pathologies) to extract as well
    :type            neigh_level:  int

    *Example of use:*

    * `extractPathologicalCells (array) <Examples/Intersector/extractPathologicals.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/extractPathologicals.py

    * `extractPathologicalCells (pyTree) <Examples/Intersector/extractPathologicalsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/extractPathologicalsPT.py

---------------------------------------


.. py:function:: Intersector.extractOuterLayers(a, N, discard_external=0)

    Extracts prescribed outer cell layers.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           N:  Number of layers to extract
    :type            N:  int
    :param           discard_external:  For volume mesh with holes (e.g. external flow), set it to 1 to extract only layers around bodies, or 0 to extract over all the outer polygons.
    :type            discard_external:  0 or 1

    *Example of use:*

    * `extractOuterLayers (array) <Examples/Intersector/extractOuters.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/extractOuters.py

    * `extractOuterLayers (pyTree) <Examples/Intersector/extractOutersPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/extractOutersPT.py

---------------------------------------

.. py:function:: Intersector.getOverlappingFaces(t1, t2, RTOL, ps_min, dir2)

    Detects all the overlapping polygons in t1 and t2. Returns the result as a list sized as the number of zones in t1. Each entry gives 2 lists : the first is the pg ids in t1 ith-zone, the second is the pg ids in t2 (joined). 

    :param           t1:  Input mesh
    :type            t1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           t2:  Input mesh
    :type            t2:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           RTOL:  relative tolerance
    :type            RTOL:  float
    :param           RTOL:  minimal dot product between normals of a pair of polygon to consider them as potentially overlapping.
    :type            RTOL:  float


    *Example of use:*

    * `getOverlappingFaces (array) <Examples/Intersector/getOverlappingFaces.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/getOverlappingFaces.py

    * `getOverlappingFaces (pyTree) <Examples/Intersector/getOverlappingFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/getOverlappingFacesPT.py

    *Tips and Notes:*

    * When asembling 2 meshes m1 and m2 where m2 is priorized , to improve the assembly quality, do before calling the boolean union:
    
    1) getOverlappingFaces (m1, skin(m2)) where skin(m2) is the external polygonal skin of m2

    2) agglomerateCellsWithSpecifiedFaces on m1 with the above list of polygons

---------------------------------------

Check Functions
--------------------------


.. py:function:: Intersector.selfX(a)

    Checks self-intersections in a mesh. Returns the first two cell indices that collide.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

    *Example of use:*

    * `selfX (array) <Examples/Intersector/selfX.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/selfX.py

    * `selfX (pyTree) <Examples/Intersector/selfXPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/selfXPT.py

---------------------------------------


.. py:function:: Intersector.diffMesh(a1, a2)

    Extracts the diff between 2 meshes. Returns 2 zones : one zone with the a1 cells that are not in a2, the second one is the reciprocal.

    :param           a1:  Input mesh
    :type            a1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           a2:  Input mesh
    :type            a2:  [array, list of arrays] or [pyTree, base, zone, list of zones]

    *Example of use:*

    * `diffMesh (array) <Examples/Intersector/diffMesh.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/diffMesh.py

    * `diffMesh (pyTree) <Examples/Intersector/diffMeshPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/diffMeshPT.py

---------------------------------------

.. py:function:: Intersector.checkCellsClosure(a)

    Checks that input mesh cells are closed, i.e. each cell' edge is shared by exactly two polygons.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

    *Example of use:*

    * `checkCellsClosure (array) <Examples/Intersector/checkCellsClosure.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/checkCellsClosure.py

    * `checkCellsClosure (pyTree) <Examples/Intersector/checkCellsClosurePT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/checkCellsClosurePT.py

---------------------------------------


.. py:function:: Intersector.computeAspectRatio(a, vmin=0.)

    For each cell, the aspect ratio with each of its neighbors is computed as the ratio of the biggest volume to the smallest one.

    The maximum over all the neighbors is chosen:

    Aspect Ratio for Cell i  =  MAX_k ( MAX(vi, vk) / MIN(vi, vk) ) where k is a neighbor. 

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param        vmin:  volume threshold.
    :type         vmin:  float

    *Example of use:*

    * `computeAspectRatio (array) <Examples/Intersector/computeAspectRatio.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/computeAspectRatio.py

    * `computeAspectRatio (pyTree) <Examples/Intersector/computeAspectRatioPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/computeAspectRatioPT.py

---------------------------------------

Conversion Functions
--------------------------

.. py:function:: Intersector.convertNGON2DToNGON3D(a)


Converts a polygon surface stored in the Cassiopee NGON format (Face/Edge) to a Face/Node format.


---------------------------------------


.. toctree::
   :maxdepth: 2   


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

