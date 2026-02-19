.. Intersector documentation master file

:tocdepth: 2

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
   :nosignatures:

   Intersector.conformUnstr 
   Intersector.booleanUnion
   Intersector.booleanIntersection  
   Intersector.booleanMinus
   Intersector.diffSurf
   Intersector.intersection
   Intersector.PyTree.XcellN
   Intersector.adaptCells
.. Intersector.booleanModifiedSolid
.. Intersector.P1ConservativeChimeraCoeffs

**-- Collision predicates**

.. autosummary::
   :nosignatures:

   Intersector.getOverlappingFaces
   Intersector.getCollidingCells
   Intersector.selfX

**-- Transformation Functions**

.. autosummary::
   :nosignatures:

   Intersector.triangulateBC
   Intersector.triangulateExteriorFaces
   Intersector.reorient
   Intersector.convexifyFaces
   Intersector.syncMacthPeriodicFaces
   Intersector.prepareCellsSplit
   Intersector.splitNonStarCells
   Intersector.agglomerateSmallCells
   Intersector.agglomerateNonStarCells
   Intersector.agglomerateCellsWithSpecifiedFaces
   Intersector.simplifyCells
   Intersector.closeCells

**-- Adaptation Specific Functions**

.. autosummary::
   :nosignatures:

   Intersector.adaptBox
   Intersector.createHMesh
   Intersector.deleteHMesh
   Intersector.conformizeHMesh
   Intersector.createSensor
   Intersector.assignData2Sensor
   Intersector.deleteSensor

**-- Metric Functions**

.. autosummary::
   :nosignatures:

   Intersector.edgeLengthExtrema
   Intersector.volumes
   Intersector.centroids
   Intersector.computeGrowthRatio
.. Intersector.statsSize

**-- Extraction Functions**

.. autosummary::
   :nosignatures:

   Intersector.extractPathologicalCells
   Intersector.extractOuterLayers
   Intersector.getCells
..   Intersector.extractUncomputables
..   Intersector.extractNthCell
..   Intersector.extractNthFace
..   Intersector.removeNthCell

**-- Check Functions**

.. autosummary::
   :nosignatures:

   Intersector.diffMesh
   Intersector.checkCellsClosure
   Intersector.checkCellsFlux
   Intersector.checkCellsVolume
   Intersector.checkForDegenCells
.. Intersector.statsUncomputableFaces
.. Intersector.extrudeUserDefinedBC   
.. Intersector.removeBaffles
.. Intersector.convert2Polyhedron
.. Intersector.oneZonePerCell

**-- Conversion Functions**

.. autosummary::
   :nosignatures:
   
   Intersector.convertNGON2DToNGON3D

Contents
#########

Main Functions
--------------------------


.. py:function:: Intersector.conformUnstr(a1, a2=None, tol=0., left_or_right=0, itermax=10)

    Makes conformal a TRI or a BAR soup (i.e. a set of elements not necessarily connected as a mesh) by detecting and solving all the collisions between elements. 
    
    Colliding elements are cut to get a conformal set. Mixing types BAR and TRI is not currently handled.

    :param a1: First input mesh (BAR or TRI).
    :type  a1: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2: Second input mesh (BAR or TRI). If s2 is 'None' self-intersections are solved over s1.
    :type  a2: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol: Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol: float
    :param left_or_right: Tells the function what to ouput : the transformed s1 (0), s2(1) or both (2).
    :type  left_or_right: 0, 1 or 2
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


.. py:function:: Intersector.booleanUnion(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1, extrude_pgs=[], multi_zone=False)

    Creates a conformal union between two components, either TRI surfaces or Polyhedral volumes. 

    :param a1:  First mesh operand.
    :type  a1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2:  Second mesh operand.
    :type  a2:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol: Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol: float
    :param preserve_right: Indicates the merging direction, either a1->a2 or a2->a1. If set to 1(0), it means a1->a2 (a2->a1), i.e. a2(a1)'s points are preserved.
    :type  preserve_right: 0 or 1
    :param solid_right: Indicates that the second operand is not penetrable, i.e. it is prioritized over the first operand a1.
    :type  solid_right: 0 or 1
    :param agg_mode: Option for agglomerating cut polygons: 0 to keep them as split triangles, 1 to get convex agglomerations and 2 to get a full agglomeration.
    :type  agg_mode: 0,1 or 2.
    :param extrude_pgs: Optional list of polygons to extrude.
    :type  extrude_pgs: list of int
    :param multi_zone: If set to True, preserve input zoning of a1 and a2 upon exit.
    :type  multi_zone: True or False

    
    **Prerequisites :**

    * External polygons must be oriented consistently and outwardly (use Transform.reorderAll before)

    
    *Tips and Notes:*

    * For assembling meshes, set solid_right to 1 and pass the prioritized mesh as second operand.
    
    * extrude_pgs: required whenever a1 and a2 are in contact and a2 is prioritized: avoids to compute useless intersections by telling what are the indices of contact polygons in a2.


    *Example of use:*

    * `Union of two spherical surface meshes (array) <Examples/Intersector/booleanUnion.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/booleanUnion.py

    * `Union of two volume meshes (pyTree) <Examples/Intersector/booleanUnionNGPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/booleanUnionNGPT.py

---------------------------------------


.. py:function:: Intersector.booleanIntersection(a1, a2, tol=0., preserve_right=1, solid_right=1, agg_mode=1)

    Computes a conformal intersection between two components, either TRI surfaces or Polyhedral volumes. 

    :param a1: First mesh operand.
    :type  a1: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2: Second mesh operand.
    :type  a2: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol: Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol: float
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

    :param a1: First mesh operand.
    :type  a1: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param a2: Second mesh operand.
    :type  a2: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol: Merging tolerance when points (existing or computed by intersections) are too close.
    :type  tol: float
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


.. py:function:: Intersector.PyTree.XcellN(t, priorities, output_type=0, rtol=0.05)

    Computes the visibility coefficient for each cell in an overset surface grid configuration t with one-to-one priorities. t can be structured or unstructured.

    Depending on the output_type argument, this function computes:

    * a ternary blanking information (0 (hidden), 1(visible) and 0.5(colliding), when output_type=0

    * a continuous blanking information (any val in [0,1] based on the ratio of the visible part of the cell), when output_type=1

    * a clipped polygonal surface (NGON format) where all the hidden surface parts have been removed, when output_type=2

    .. image:: _static/xcelln_conf.jpg
      :width: 21%
    .. image:: _static/xcelln_mode0.jpg
      :width: 21%
    .. image:: _static/xcelln_mode1.jpg
      :width: 21%
    .. image:: _static/xcelln_mode2.jpg
      :width: 21%

    *From left to right: Sphere made of 2 overset patches, results with output_type=0,1 and 2 displayed for the non-prioritized patch.*

    When output_type is 0 or 1, a 'xcelln' field is added to each zone of the PyTree.
    When output_type is 2, the fields defined at center are transferred in the output mesh. 

    :param t:  an overset surface mesh where components are separated in different bases.
    :type  t:  [pyTree, base, zone, list of zones].
    :param priorities:  list of one-to-one pairs of priorities between components.
    :type  priorities:  list of pairs of integers.
    :param output_type: ternary blanking field.
    :type  output_type: 0, 1 or 2.
    :param rtol: relative tolerance for detecting and computing intersections.
    :type  rtol: float.
    
        
    **Prerequisites :**

    * All the surface patches must be oriented consistently and outwardly (use Transform.reorderAll before)

    **Tips and Notes:**

    * each component must set in a separate base.
    * priorisation and computation is done between components only : not between zones of the same base.
    * wall boundaries are considered whatever the priorisation:

    .. image:: _static/xcelln_prior_conf.jpg
      :width: 28%
    .. image:: _static/xcelln_prior1.jpg
      :width: 28%
    .. image:: _static/xcelln_prior2.jpg
      :width: 28%

    *Example of wall treatment on the fuselage/collar zone of a CRM configuration when output_type is 2: if the collar is prioritized (middle) or not (right), the part of the fuselage falling inside the wing is clipped.*
    

    *Example of use:*

    * `xcelln field on structured configuration <Examples/Intersector/xcellnSphYinYangPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/xcellnSphYinYangPT.py

---------------------------------------


.. py:function:: Intersector.adaptCells(a, sensdata=None, sensor_type = 0, smoothing_type = 0, itermax=-1, subdiv_type=0, hmesh=None, sensor=None, conformize=1)

    Adapts a cells (any basic cells - TETRA, PENTA, PYRA, HEXA - currently) with respect to the sensor and its data (sensdata).
    
    Adaptation is a per-cell octal 2:1 decomposition.
    
    With a sensor_type equal to :

    * 0 (geometrical sensor), sensdata must contain vertices : a will be refined until any a cell contains at most 1 data vertex.
    
    * 1 (xsensor), sensdata must be a mesh *m*, ant its connectivity is also taken into account by adding refinement wherever a cells are crossed by m edges.
    
    * 2 (nodal sensor), sensdata are nodal values giving the number of required subdivision around that node.
    
    * 3 (cell sensor), sensdata are cell values giving the number of required subdivision per cell.

    :param           a:  Input mesh (NGON format)
    :type            a:  [array] or [ single zone pyTree (currently)]
    :param           sensdata:  Data for the sensor
    :type            sensdata:  [array, list of arrays] or [pyTree, base, zone, list of zones] for sensor type 0 and 1. Numpy of integers for other sensors.
    :param           sensor_type:  type of sensor. geometrical (0), xensor (1), nodal sensor (2), cell sensor (3)
    :type            sensor_type:  int
    :param           smoothing_type:  first-neighborhood (0), shell-neighborhood(1)
    :type            smoothing_type:  int
    :param           itermax:  maximum nb of generations
    :type            itermax:  int
    :param           subdiv_type:  type of adaptation, currently only isotropic (0).
    :type            subdiv_type:  int
    :param           hmesh:  structure that holds the hierarchical genealogy structure in case of successive adaptations on a mesh. Instantiated with Intersector.createHMesh
    :type            hmesh:  hook
    :param           sensor:  structure that holds the sensor and its data in case of successive adaptations on a mesh
    :type            sensor:  hook
    :param           conformize:  conformal output mesh enabled (1), disabled (0). Enabled by default.
    :type            conformize:  int

    *Example of use:*

    * `adaptCells (array) <Examples/Intersector/adaptCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCells.py

    * `adaptCells (pyTree) <Examples/Intersector/adaptCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCellsPT.py

    *Tips and Notes:*

    * Do this transformation before calling any Volume-Volume boolean operations in order to improve the mesh quality of the result.
    * When the input mesh has any kind of polyhedral elements, only basic elements will be considered currently for adaptation. but the result wil be conformal, the non-handled elements will modified to respect the conformity. 

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


.. py:function:: Intersector.reorient(a)

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

.. py:function:: Intersector.syncMacthPeriodicFaces(t, rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.], translation=[0.,0.,0.], tol=-0.01, unitAngle=None, reorient=True)

    Force periodicity when some faces should be periodic but are not due to connectivity inconsistencies (tipically overdefined or splitted faces).
    This function ensures that a following call to Connector.connectMatchPeriodic succeeds when the mesh has been procuced with the intersection strategy. 
    Periodicity can be defined either by rotation or translation.

    :param           t:  Input mesh
    :type            t:  [pyTree, base, zone, list of zones]
    :param           tol:  tolerance. Negative value (in [-1; 0]) specifies a relative value base on min edge length in the mesh.
    :type            a:  float

    *Example of use:*

    * `Synchronize faces(pyTree) <Examples/Intersector/syncPerioFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/syncPerioFacesPT.py

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

    First strategy to eradicate bad cells : Splits non-centroid-star-shaped (NCSS) cells into two cells. These cells might be NCSS as well so this function should be called several times to get rid off the pathologies. Some call agglomerateSmallCells should be done afterwards to balance the growth ratio.

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
    If 2 cells share 2 polygons that are connected (sharing an edge) and their dihedral angle is below the angular_threshold, then the polygons are agglomerated upon exit. The angular threshold (expressed in radian) is the maximum aboslute deviation around the planar position.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param      treat_externals:  Process outer polygons (1) or not (0).
    :type       treat_externals:  0 or 1
    :param      angular_threshold:  Largest angular deviation admitted between adjacent polygons in order to allow their agglomeration.
    :type       angular_threshold:  float

    *Example of use:*

    * `simplifyCells (array) <Examples/Intersector/simplifyCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/simplifyCells.py

    * `simplifyCells (pyTree) <Examples/Intersector/simplifyCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/simplifyCellsPT.py

---------------------------------------


.. py:function:: Intersector.agglomerateSmallCells(a, vmin=0., vratio=0.01, angular_threshold = 1.e-12)

    Agglomerates cells that are too small (below vmin) or having a poor growth ratio with a neighbor (below vratio) with the best neighbor available. The agglomeration process does not create non-star-shaped agglomerates.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param      vmin:  volume threshold.
    :type       vmin:  float
    :param      vratio: growth ratio threshold.
    :type       vratio:  float
    :param      angular_threshold:  Largest angular deviation admitted between adjacent polygons in order to allow their agglomeration.
    :type       angular_threshold:  float


    *Tips and Notes:*

    * See :any:`computeGrowthRatio` to get the definition of the computed growth ratio.

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

    * When assembling 2 meshes m1 and m2 where m2 is priorized, to improve the assembly quality, do before calling the boolean union:

    1) getOverlappingFaces (m1, skin(m2)) where skin(m2) is the external polygonal skin of m2

    2) agglomerateCellsWithSpecifiedFaces on m1 with the above list of polygons
    

---------------------------------------


.. py:function:: Intersector.closeCells(a)

    Closed any polyhedral cell in a mesh which is open because it has, and only has, hanging nodes on its edges.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

    *Example of use:*

    * `closeCells (array) <Examples/Intersector/closeCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/closeCells.py

    * `closeCells (pyTree) <Examples/Intersector/closeCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/closeCellsPT.py

    *Tips and Notes:*

    * Do this transformation whenever you need to use a surface algorithm on the octree (e.g. :any:`reorient`)


---------------------------------------

Adaptation Specific Functions
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


.. py:function:: Intersector.createHMesh(a, subdiv_type= 0)

    Builds a hierarchcial mesh structure for a and returns a hook 

    :param           a:  Input points cloud
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           subdiv_type:  type of adaptation, currently only isotropic (0).
    :type            subdiv_type:  int

    *Example of use:*

    * `createHMesh <Examples/Intersector/adaptCellsDynPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCellsDynPT.py

    
---------------------------------------

.. py:function:: Intersector.deleteHMesh(hmesh)

    Deletes a hierarchcial mesh

    :param           hmesh:  hmesh hook
    :type            hmesh:  hook

    *Example of use:*

    * `deleteHMesh <Examples/Intersector/adaptCellsDynPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCellsDynPT.py


---------------------------------------

.. py:function:: Intersector.conformizeHMesh(a, hooks)

    Converts the basic element leaves of a hierarchical mesh (`hooks` is a list of pointers to hiearchical zones) to a conformal polyhedral mesh.
    Each hierarchical zone is referring to a zone in the original Pytree `t`. So the mesh is replaced in the returned tree and the BCs/Joins/Fields are transferred.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           hooks:  list of pointers to hiearchical zones
    :type            a:  hooks

    *Example of use:*

    * `conformizeHMesh <Examples/Intersector/adaptCellsDynPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCellsDynPT.py

---------------------------------------

.. py:function:: Intersector.createSensor(hmesh, sensor_type = 0, smoothing_type=0 , itermax = -1)

    Creates a sensor and returns a hook on it.

    :param           hmesh:  hmesh hook
    :type            hmesh:  hook
    :param           sensor_type:  type of sensor. geometrical (0), xensor (1), nodal sensor (2), cell sensor (3)
    :type            sensor_type:  int
    :param           smoothing_type:  first-neighborhood (0), shell-neighborhood(1)
    :type            smoothing_type:  int
    :param           itermax:  maximum nb of generations
    :type            itermax:  int

    *Example of use:*

    * `createSensor <Examples/Intersector/adaptCellsDynPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCellsDynPT.py

---------------------------------------

.. py:function:: Intersector.assignData2Sensor(sensor, sensdata)

    Assigns data to a sensor.

    :param           sensor:  sensor hook
    :type            sensor:  hook
    :param           sensdata:  Data for the sensor
    :type            sensdata:  [array, list of arrays] or [pyTree, base, zone, list of zones] for sensor type 0 and 1. Numpy of integers for other sensors.
    
    *Example of use:*

    * `assignData2Sensor <Examples/Intersector/adaptCellsDynPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/adaptCellsDynPT.py

---------------------------------------

.. py:function:: Intersector.deleteSensor(sensor)

   Deletes a sensor

   :param           sensor:  sensor hook
   :type            sensor:  hook

   *Example of use:*

   * `deleteSensor <Examples/Intersector/adaptCellsDynPT.py>`_:

   .. literalinclude:: ../build/Examples/Intersector/adaptCellsDynPT.py


---------------------------------------

Metric Functions
---------------------------------------

.. py:function:: Intersector.edgeLengthExtrema(a)

    Returns the minimum edge length in a. 

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

---------------------------------------

.. py:function:: Intersector.volumes(a)

    Returns the cell volumes as a field (PyTree) or a numpy of floats. 

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

---------------------------------------

.. py:function:: Intersector.centroids(a)

    Returns the cell centroids as a points cloud. 

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

---------------------------------------

.. py:function:: Intersector.computeGrowthRatio(a, vmin=0.)

    For each cell, the growth ratio with each of its neighbors is computed as the ratio of the biggest volume to the smallest one.

    The maximum over all the neighbors is chosen:

    Growth Ratio for Cell i  =  MAX_k ( MAX(vi, vk) / MIN(vi, vk) ) where k is a neighbor. 

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param        vmin:  volume threshold.
    :type         vmin:  float

    *Example of use:*

    * `computeGrowthRatio (array) <Examples/Intersector/computeGrowthRatio.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/computeGrowthRatio.py

    * `computeGrowthRatio (pyTree) <Examples/Intersector/computeGrowthRatioPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/computeGrowthRatioPT.py

---------------------------------------


Extraction Functions
---------------------------------------


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


.. py:function:: Intersector.getCells(a, ids, are_face_ids = True)

    Returns the cells in t1 having specified faces or cell ids.

    :param     a:  Input mesh
    :type      a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param     ids:  face or cell ids
    :type      ids:  numpy of ints
    :param     are_face_ids:  Tells whether the ids are referring to faces or cells.
    :type      are_face_ids:  boolean

    *Example of use:*

    See :any:`getCollidingCells` for an example.



---------------------------------------

.. py:function:: Intersector.getOverlappingFaces(t1, t2, RTOL, ps_min, dir2)

    Detects all the overlapping polygons in t1 and t2. Returns the result as a list sized as the number of zones in t1. Each entry gives 2 lists : the first is the pg ids in t1 ith-zone, the second is the pg ids in t2 (joined). 

    :param           t1:  Input mesh
    :type            t1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           t2:  Input mesh
    :type            t2:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           RTOL:  relative tolerance
    :type            RTOL:  float
    :param           ps_min:  minimal dot product between normals of a pair of polygon to consider them as potentially overlapping.
    :type            ps_min:  float
    :param           dir2:  given direction to compare t1's faces with. If None, t2's normals are used.
    :type            dir2:  tuple


    *Example of use:*

    * `getOverlappingFaces (array) <Examples/Intersector/getOverlappingFaces.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/getOverlappingFaces.py

    * `getOverlappingFaces (pyTree) <Examples/Intersector/getOverlappingFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/getOverlappingFacesPT.py

    *Tips and Notes:*

    * When assembling 2 meshes m1 and m2 where m2 is priorized, to improve the assembly quality, do before calling the boolean union:
    
    1) getOverlappingFaces (m1, skin(m2)) where skin(m2) is the external polygonal skin of m2

    2) agglomerateCellsWithSpecifiedFaces on m1 with the above list of polygons

---------------------------------------

.. py:function:: Intersector.getCollidingCells(t1, t2, RTOL)

    Returns the list of cells in t1 and t2 that are colliding. Possible combinations of mesh types for (t1,t2) are (volume,volume), (volume,surface), (surface, polyline).

    :param           t1:  Input mesh
    :type            t1:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           t2:  Input mesh
    :type            t2:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param           RTOL:  relative tolerance
    :type            RTOL:  float
    

    *Example of use:*

    * `getCollidingCells (array) <Examples/Intersector/getCollidingCells.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/getCollidingCells.py

    * `getCollidingCells (pyTree) <Examples/Intersector/getCollidingCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/getCollidingCellsPT.py


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

.. py:function:: Intersector.checkCellsFlux(a)

    Computes the cell fluxes using the ParentElement elsA's node for orientation.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

    *Example of use:*

    * `checkCellsFlux (pyTree) <Examples/Intersector/checkCellsFluxPT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/checkCellsFluxPT.py

---------------------------------------

.. py:function:: Intersector.checkCellsVolume(a)

    Computes the minimum volume using the ParentElement elsA's node for orientation.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

    *Example of use:*

    * `checkCellsVolume (pyTree) <Examples/Intersector/checkCellsVolumePT.py>`_:

    .. literalinclude:: ../build/Examples/Intersector/checkCellsVolumePT.py

---------------------------------------

.. py:function:: Intersector.checkForDegenCells(a)

    Checks if there are any cell with less than 4 faces.

    :param           a:  Input mesh
    :type            a:  [array, list of arrays] or [pyTree, base, zone, list of zones]

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

