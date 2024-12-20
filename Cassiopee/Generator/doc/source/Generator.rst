.. Generator documentation master file

:tocdepth: 2

Generator: mesh generation module
==================================


Preamble
########

Generator module works on arrays (as defined in Converter) or on CGNS/python trees (pyTrees) 
containing grid information (coordinates must be defined).

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use with the array interface, you have to import Generator module::

   import Generator as G

For use with the pyTree interface::

    import Generator.PyTree as G


.. py:module:: Generator

List of functions
##################

**-- Basic grid generation**

.. autosummary::
   :nosignatures:

   Generator.cart
   Generator.cartr1
   Generator.cartr2
   Generator.cartHexa
   Generator.cartTetra
   Generator.cartPenta
   Generator.cartPyra
   Generator.cartNGon
   Generator.cylinder
   Generator.cylinder2
   Generator.cylinder3

**-- General purpose grid generators**

.. autosummary::
   :nosignatures:

   Generator.delaunay
   Generator.constrainedDelaunay
   Generator.checkDelaunay
   Generator.T3mesher2D
   Generator.tetraMesher
   Generator.TFI
   Generator.TFITri
   Generator.TFIO
   Generator.TFIHalfO
   Generator.TFIMono
   Generator.TFIStar
   Generator.hyper2D
   Generator.PolyLine.polyLineMesher
   Generator.PolyC1.polyC1Mesher
   Generator.pointedHat
   Generator.stitchedHat
   Generator.surfaceWalk
   Generator.collarMesh

**-- Cartesian grid generators**

.. autosummary::
   :nosignatures:

   Generator.gencartmb
   Generator.octree
   Generator.octree2Struct
   Generator.adaptOctree
   Generator.expandLayer
   Generator.PyTree.cartRx
   Generator.PyTree.cartRx3

**-- Operations on meshes**

.. autosummary::
   :nosignatures:

   Generator.close
   Generator.zip
   Generator.selectInsideElts
   Generator.map
   Generator.mapSplit
   Generator.refine
   Generator.mapCurvature
   Generator.densify
   Generator.grow
   Generator.stack
   Generator.addNormalLayers
   Generator.TTM
   Generator.snapFront
   Generator.snapSharpEdges

**-- Operations on surface meshes**

.. autosummary::
   :nosignatures:

   Generator.fittingPlaster
   Generator.gapfixer
   Generator.gapsmanager
   Generator.mmgs

**-- Information on generated meshes**

.. autosummary::
   :nosignatures:

   Generator.barycenter
   Generator.bbox
   Generator.bboxOfCells
   Generator.BB
   Generator.CEBBIntersection
   Generator.bboxIntersection
   Generator.checkPointInCEBB
   Generator.getVolumeMap
   Generator.getNormalMap
   Generator.getSmoothNormalMap
   Generator.getOrthogonalityMap
   Generator.getRegularityMap
   Generator.getAngleRegularityMap
   Generator.getTriQualityMap
   Generator.getCellPlanarity
   Generator.getCircumCircleMap
   Generator.getInCircleMap
   Generator.getEdgeRatio
   Generator.getMaxLength
   Generator.checkMesh

**-- Operations on distributions**

.. autosummary::
   :nosignatures:

   Generator.enforceX
   Generator.enforceMoinsX
   Generator.enforcePlusX
   Generator.enforceLine
   Generator.enforcePoint
   Generator.enforceCurvature
   Generator.addPointInDistribution



Contents
#########

Basic grid generation
--------------------------


.. py:function:: Generator.cart((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))

    Create a structured Cartesian mesh with ni x nj x nk points starting from point (xo,yo,zo) and of step (hi,hj,hk).

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param (hi,hj,hk):  values of advancing step in the three directions
    :type  (hi,hj,hk):  3-tuple of floats
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 1D, 2D or 3D structured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Cartesian mesh generation (array) <Examples/Generator/cart.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cart.py

    * `Cartesian mesh generation (pyTree) <Examples/Generator/cartPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartPT.py

---------------------------------------

.. py:function:: Generator.cartr1((xo,yo,zo), (hi,hj,hk), (ri,rj,rk), (ni,nj,nk))

    Create a structured Cartesian mesh with geometric distribution of factors r.
    (hi,hj,hk) are the steps of first cell.

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param (hi,hj,hk):  Step of first cell in the three directions
    :type  (hi,hj,hk):  3-tuple of floats
    :param (ri,rj,rk):  geometric factors in the three directions
    :type  (ri,rj,rk):  3-tuple of floats
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 1D, 2D or 3D structured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Geometric Cartesian mesh generation (array) <Examples/Generator/cartr1.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartr1.py

    * `Geometric Cartesian mesh generation (pyTree) <Examples/Generator/cartr1PT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartr1PT.py


---------------------------------------

.. py:function:: Generator.cartr2((xo,yo,zo), (hi,hj,hk), (ri,rj,rk), (xf,yf,zf))

    Create a structured Cartesian mesh with geometric distribution of factors r.
    (hi,hj,hk) are the steps of first cell.
    (xf,yf,zf) are coordinates of last point.
    (hi,hj,hk) or (ri,rj,rk) can be slightly modified to match (xf,yf,zf).

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param (hi,hj,hk):  Step of first cell in the three directions
    :type  (hi,hj,hk):  3-tuple of floats
    :param (ri,rj,rk):  geometric factors in the three directions
    :type  (ri,rj,rk):  3-tuple of floats
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 1D, 2D or 3D structured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Geometric Cartesian mesh with fixed last point (array) <Examples/Generator/cartr2.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartr2.py

    * `Geometric Cartesian mesh with fixed last point (pyTree) <Examples/Generator/cartr2PT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartr2PT.py

---------------------------------------

.. py:function:: Generator.cartHexa((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))

    Create an unstructured hexahedral mesh defined from a Cartesian grid of ni x nj x nk points starting from point (xo,yo,zo) and of step (hi,hj,hk). Type of elements are 'QUAD' for 2D arrays and 'HEXA' for 3D arrays.

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param (hi,hj,hk):  values of advancing step in the three directions
    :type  (hi,hj,hk):  3-tuple of floats
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 1D, 2D or 3D unstructured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Cartesian hexa mesh generation (array) <Examples/Generator/cartHexa.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartHexa.py

    * `Cartesian hexa mesh generation (pyTree) <Examples/Generator/cartHexaPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartHexaPT.py

---------------------------------------

.. py:function:: Generator.cartTetra((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))

    Create an unstructured tetrahedral mesh defined from a Cartesian grid of ni x nj x nk points starting from point (xo,yo,zo) and of step (hi,hj,hk). Type of elements are 'TRI' for 2D arrays and 'TETRA' for 3D arrays.

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param (hi,hj,hk):  values of advancing step in the three directions
    :type  (hi,hj,hk):  3-tuple of floats
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 1D, 2D or 3D unstructured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Cartesian tetra mesh generation (array) <Examples/Generator/cartTetra.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartTetra.py

    * `Cartesian tetra mesh generation (pyTree) <Examples/Generator/cartTetraPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartTetraPT.py

---------------------------------------

.. py:function:: Generator.cartPenta((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))

    Create an unstructured prismatic mesh defined from a regular Cartesian mesh. The initial Cartesian mesh is defined by ni x nj x nk points starting from point (xo,yo,zo) and of step (hi,hj,hk). Type of elements is 'PENTA'.

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param (hi,hj,hk):  values of advancing step in the three directions
    :type  (hi,hj,hk):  3-tuple of floats
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 1D, 2D or 3D unstructured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Cartesian penta mesh generation (array) <Examples/Generator/cartPenta.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartPenta.py

    * `Cartesian penta mesh generation (pyTree) <Examples/Generator/cartPentaPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartPentaPT.py

---------------------------------------

.. py:function:: Generator.cartPyra((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))

    Create an unstructured pyramidal mesh defined from a regular Cartesian mesh. The initial Cartesian mesh is defined by ni x nj x nk points starting from point (xo,yo,zo) and of step (hi,hj,hk). Type of elements is 'PYRA'.

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param (hi,hj,hk):  values of advancing step in the three directions
    :type  (hi,hj,hk):  3-tuple of floats
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 1D, 2D or 3D unstructured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Cartesian pyra mesh generation (array) <Examples/Generator/cartPyra.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartPyra.py

    * `Cartesian pyra mesh generation (pyTree) <Examples/Generator/cartPyraPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartPyraPT.py

---------------------------------------

.. py:function:: Generator.cartNGon((xo,yo,zo), (hi,hj,hk), (ni,nj,nk))

    Create a NGON mesh defined from a regular Cartesian mesh. The initial Cartesian mesh is defined by ni x nj x nk points starting from point (xo,yo,zo) and of step (hi,hj,hk). Type of elements is 'NGON'.

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param (hi,hj,hk):  values of advancing step in the three directions
    :type  (hi,hj,hk):  3-tuple of floats
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 1D, 2D or 3D unstructured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Cartesian ngon mesh generation (array) <Examples/Generator/cartNGon.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartNGon.py

    * `Cartesian ngon mesh generation (pyTree) <Examples/Generator/cartNGonPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartNGonPT.py

---------------------------------------

.. py:function:: Generator.cylinder((xo,yo,zo), R1, R2, tetas, tetae, H, (ni,nj,nk))

    Create a regular cylindrical grid (or a portion of cylinder between tetas and tetae) with ni x nj x nk points, of center-bottom point (xo,yo,zo), of inner radius R1, outer radius R2 and height H. For a direct mesh, use tetae < tetas.

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param R1:  value of inner radius
    :type  R1:  float
    :param R2:  value of outer radius
    :type  R2:  float
    :param tetas:  start angle (in degree)
    :type  tetas:  float
    :param tetae:  end angle (in degree)
    :type  tetae:  float
    :param (ni,nj,nk):  number of points in each direction
    :type  (ni,nj,nk):  3-tuple of integers
    :return: a 3D structured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Regular cylinder mesh generation (array) <Examples/Generator/cylinder.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cylinder.py

    * `Regular cylinder mesh generation (pyTree) <Examples/Generator/cylinderPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cylinderPT.py

---------------------------------------

.. py:function:: Generator.cylinder2((xo,yo,zo), R1, R2, tetas, tetae, H, arrayR, arrayTeta, arrayZ)

    Create an irregular cylindrical grid (or a portion of cylinder between tetas and tetae) with ni x nj x nk points, of center-bottom point (xo,yo,zo), of inner radius R1, outer radius R2, height H and with distributions in r, teta, z. Distributions are arrays defining 1D meshes (x and i varying) giving a distribution in [0,1]. Their number of points gives ni, nj, nk.

    :param (xo,yo,zo):  coordinates of the starting point
    :type  (xo,yo,zo):  3-tuple of floats
    :param R1:  value of inner radius
    :type  R1:  float
    :param R2:  value of outer radius
    :type  R2:  float
    :param tetas:  start angle (in degree)
    :type  tetas:  float
    :param tetae:  end angle (in degree)
    :type  tetae:  float
    :param H:  value of cylinder height
    :type  H:  float
    :param arrayR:  distribution along radius
    :type  arrayR:  array
    :param arrayTeta:  distribution along azimuth
    :type  arrayTeta:  array
    :param arrayZ:  distribution along height
    :type  arrayZ:  array
    :return: a 3D structured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Irregular cylinder mesh generation (array) <Examples/Generator/cylinder2.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cylinder2.py

    * `Irregular cylinder mesh generation (pyTree) <Examples/Generator/cylinder2PT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cylinder2PT.py

---------------------------------------

.. py:function:: Generator.cylinder3(a, tetas, tetae, arrayTeta)

    Create an irregular cylindrical grid (or a portion of cylinder between tetas and tetae) from a xz plane mesh defined by a and a teta distribution defined by arrayTeta.

    :param a:  definition of the xz plane mesh
    :type  a:  [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param tetas:  start angle (in degree)
    :type  tetas:  float
    :param tetae:  end angle (in degree)
    :type  tetae:  float
    :param arrayTeta:  distribution along azimuth
    :type  arrayTeta:  array
    :return: a 3D structured mesh
    :rtype: array or pyTree zone

    *Example of use:*

    * `Irregular cylinder mesh generation from a xz plane (array) <Examples/Generator/cylinder3.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cylinder3.py

    * `Irregular cylinder mesh generation from a xz plane (pyTree) <Examples/Generator/cylinder3PT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cylinder3PT.py

---------------------------------------------------------------------------

General purpose grid generator
---------------------------------------

.. py:function:: Generator.delaunay(a, tol=1.e-10, keepBB=0) 

    Create a 2D Delaunay type mesh from an array. The array can be a 2D structured array, or an unstructured array of type 'NODE', 'TRI' or 'QUAD'. Tol is a geometric tolerance. Points nearer than tol are merged. If keepBB is set to 1, the bounding box is kept in the final triangulation.

    :param a:  structured or unstructured 2D mesh
    :type  a:  [array] or [zone]
    :param tol:  geometric tolerance
    :type  tol:  float
    :param keepBB:  keep bounding box in result?
    :type  keepBB:  integer (0 or 1)
    :return: a 2D unstructured mesh
    :rtype: Identical to a

    *Example of use:*

    * `2D Delaunay mesh generation (array) <Examples/Generator/delaunay.py>`_:

    .. literalinclude:: ../build/Examples/Generator/delaunay.py

    * `2D Delaunay mesh generation (pyTree) <Examples/Generator/delaunayPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/delaunayPT.py

---------------------------------------

.. py:function:: Generator.constrainedDelaunay(c, tol=1.e-10, keepBB=0) 

    Create a constrained Delaunay triangulation of the convex hull of a contour c. Contour must be a BAR-array and must be in the plane (x,y). Tol is a geometric tolerance. Points nearer than tol are merged. If keepBB is set to 1, the bounding 
    box is kept in the final triangulation.

    :param c:  contour in xy plane
    :type  c:  BAR-array
    :param tol:  geometric tolerance
    :type  tol:  float
    :param keepBB:  keep bounding box in result?
    :type  keepBB:  integer (0 or 1)
    :return: a 2D unstructured mesh
    :rtype: Identical to a

    *Example of use:*

    * `2D Delaunay mesh generation from a contour (array) <Examples/Generator/constrainedDelaunay.py>`_:

    .. literalinclude:: ../build/Examples/Generator/constrainedDelaunay.py

    * `2D Delaunay mesh generation from a contour (pyTree) <Examples/Generator/constrainedDelaunayPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/constrainedDelaunayPT.py

---------------------------------------

.. py:function:: Generator.checkDelaunay(c, tri) 

    Check if the Delaunay triangulation defined in tri is inside the contour c.

    :param c:  contour in xy plane
    :type  c:  BAR-array
    :param tri:  2D Delaunay triangulation mesh
    :type  tri:  array or pyTree
    :return: contour
    :rtype: BAR-array

    *Example of use:*

    * `Check Delaunay triangulation wrt. contour (array) <Examples/Generator/checkDelaunay.py>`_:

    .. literalinclude:: ../build/Examples/Generator/checkDelaunay.py

    * `Check Delaunay triangulation wrt. contour (pyTree) <Examples/Generator/checkDelaunayPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/checkDelaunayPT.py

---------------------------------------

.. py:function:: Generator.T3mesher2D(a, triangulateOnly=0, grading=1.2, metricInterpType=0) 

    Creates a 2D Delaunay mesh given a BAR defined in a. If triangulateOnly=1 then only points of a are triangulated, if triangulateOnly=0, then interior points are inserted.
    
    The grading parameter allows to control the growth ratio of the mesh metric : a value greater(lesser) than 1. tends to produce a coarser (finer) mesh in the region far from the boundaries. A value equal to 1. provides a uniform mesh over the domain. This grading is related to the metric field, it is not the size ratio between two adjacent edges or triangles.
    
    The metricInterpType parameter controls the metrics interpolation type: either linear or geometric. A geometric metric interpolation tends to promote smaller sizes.

    :param c:  BAR-contour (soup of conformal edges defining an enclosed 2D-domain, can be non-manifold, i.e. having inner edges or subdomains)
    :type  c:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param triangulateOnly:  insertion or not of interior points for the triangulation
    :type  triangulateOnly:  integer (0 or 1)
    :param grading:  metric growth ratio
    :type  grading:  float (strictly positive value)
    :param metricInterpType:  metric interpolation type, linear(0) or geometric(1)
    :type  metricInterpType:  integer (0 or 1)
    :return: 2D mesh
    :rtype: Identical to input

    *Example of use:*

    * `2D Delaunay mesh generation from a BAR (array) <Examples/Generator/T3mesher2D.py>`_:

    .. literalinclude:: ../build/Examples/Generator/T3mesher2D.py

    * `2D Delaunay mesh generation from a BAR (pyTree) <Examples/Generator/T3mesher2DPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/T3mesher2DPT.py

---------------------------------------

.. py:function:: Generator.tetraMesher(a, maxh=-1., grading=0.4, algo=1, optionString="") 

    Create a 3D tetra mesh given a TRI surface defined in a. If the TRI surface has external normals, tetras are filled inside the surface. If algo=0, netgen is used, if algo=1, tetgen is used.

    :param c:  triangulated surface mesh
    :type  c:  [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param maxh: max cell size of generated mesh [tetgen]
    :type maxh: float 
    :param grading: max adjacent cell ratio [tetgen]
    :type grading: float 
    :param algo:  choice parameter between netgen and tetgen
    :type  algo:  integer (0 or 1)
    :param optionString: string of options identical to tetgen [tetgen]
    :type optionString: string
    :return: 3D mesh
    :rtype: [array] or [pyTree zone]

    *Example of use:*

    * `3D tetra mesh generation (array) <Examples/Generator/tetraMesher.py>`_:

    .. literalinclude:: ../build/Examples/Generator/tetraMesher.py

    * `3D tetra mesh generation (pyTree) <Examples/Generator/tetraMesherPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/tetraMesherPT.py

---------------------------------------

.. py:function:: Generator.TFI([imin, imax, jmin, jmax, kmin, kmax])

    Generate a mesh by transfinite interpolation (TFI). Generated mesh can be 2D or 3D structured, or unstructured TRI or PENTA mesh. Warning: the boundaries can be in a different order from the examples below, except for the PENTA TFI meshes.
    2D structured mesh is built from imin, imax, jmin, jmax boundaries.
    3D structured mesh is built from imin, imax, jmin, jmax, kmin, kmax boundaries.
    Dimensions must be equal for each pair (imin,imax), (jmin,jmax)...
    TRI mesh is built from imin, jmin, diag boundaries. Each boundary is a structured array with the same dimension. PENTA mesh is built from Tmin, Tmax triangles boundary and imin, imax, diag boundaries. Tmin, Tmax must be structured triangles of dimension nxn. imin, jmin, diag must be structured n*p arrays.

    :param imin:  I-direction minimum boundary
    :type  imin:  array
    :param imax:  I-direction maximum boundary
    :type  imax:  array
    :param jmin:  J-direction minimum boundary
    :type  jmin:  array
    :param jmax:  J-direction maximum boundary
    :type  jmax:  array
    :param kmin:  K-direction minimum boundary
    :type  kmin:  array
    :param kmax:  K-direction maximum boundary
    :type  kmax:  array
    :param diag:  third direction boundary for TRI or PENTA meshes
    :type  diag:  array
    :return: 2D or 3D mesh
    :rtype: array or pyTree

    *Example of use:*

    * `TFI mesh generation (array) <Examples/Generator/TFI.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFI.py

    * `TFI mesh generation (pyTree) <Examples/Generator/TFIPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIPT.py

---------------------------------------

.. py:function:: Generator.TFITri(a1, a2, a3)

    Generate three structured meshes by transfinite interpolation around three given curves a1, a2, a3. N3-N2+N1 must be odd.

    :param a1:  first curve
    :type  a1:  array
    :param a2:  second curve
    :type  a2:  array
    :param a3:  third curve
    :type  a3:  array
    :return: 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `TFI structured mesh generation between three curves (array) <Examples/Generator/TFITri.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFITri.py

    * `TFI structured mesh generation between three curves (pyTree) <Examples/Generator/TFITriPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFITriPT.py

---------------------------------------

.. py:function:: Generator.TFIO(a, weight=None)

    Generate five meshes by transfinite interpolation around one given curve a. The number of points of a must be odd.

    :param a:  curve
    :type  a:  array or zone
    :param weight: if provided, the weight sets the size of the inside square
    :type weight: None, float or list of floats
    :return: 2D structured mesh (butterfly O-H topology)
    :rtype: array or pyTree

    *Example of use:*

    * `Butterfly structured mesh generation by TFI (array) <Examples/Generator/TFIO.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIO.py

    * `Butterfly structured mesh generation by TFI (pyTree) <Examples/Generator/TFIOPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIOPT.py

---------------------------------------

.. py:function:: Generator.TFIHalfO(a1, a2)

    .. A1.O0.D0

    Generate four meshes by transfinite interpolation around two given curves a1 and a2 forming a half-O. 
    N1, the number of points of a1 and N2, the number of points of a2 must be odd.

    :param a1:  first curve
    :type  a1:  array or Zone
    :param a2:  second curve
    :type  a2:  array or Zone
    :return: 2D structured mesh (half butterfly C-H topology)
    :rtype: array or Zone

    *Example of use:*

    * `Half-Butterfly structured mesh generation by TFI (array) <Examples/Generator/TFIHalfO.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIHalfO.py

    * `Half-Butterfly structured mesh generation by TFI (pyTree) <Examples/Generator/TFIHalfOPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIHalfOPT.py

---------------------------------------

.. py:function:: Generator.TFIMono(a1, a2)

    Generate one mesh by transfinite interpolation around two given curves a1 and a2 forming a half-O. N1-N2 must be even.

    :param a1:  first curve
    :type  a1:  array
    :param a2:  second curve
    :type  a2:  array
    :return: 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `TFI structured mesh generation between two curves (array) <Examples/Generator/TFIMono.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIMono.py

    * `TFI structured mesh generation between two curves (pyTree) <Examples/Generator/TFIMonoPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIMonoPT.py

---------------------------------------

.. py:function:: Generator.TFIStar(a)

    Generate a set of meshes by transfinite interpolation a list of given curves.

    :param a:  input curves
    :type  a:  list of structured arrays or zones
    :return: list of 2D structured mesh
    :rtype: arrays or pyTree

    *Example of use:*

    * `TFI structured mesh generation of a set of curves (array) <Examples/Generator/TFIStar.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIStar.py

    * `TFI structured mesh generation of a set of curves (pyTree) <Examples/Generator/TFIStarPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TFIStarPT.py

---------------------------------------

.. py:function:: Generator.hyper2D(line, distrib, "C")

    Generate an hyperbolic mesh (2D) of "C" or "O" type from a from a line defined by line and from a distribution defined by distrib. The resulting mesh is nearly orthogonal.

    :param line:  starting line of the hyperbolic mesh
    :type  line:  array
    :param distrib:  distribution orthogonal to the line
    :type  distrib:  array
    :return: 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Hyperbolic structured mesh generation from a line (array) <Examples/Generator/hyper2d.py>`_:

    .. literalinclude:: ../build/Examples/Generator/hyper2d.py

    * `Hyperbolic structured mesh generation from a line (pyTree) <Examples/Generator/hyper2dPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/hyper2dPT.py

---------------------------------------

.. py:function:: Generator.PolyLine.polyLineMesher(a, h, hf, density)

    Generate a 2D mesh around a 2D polyline where a is the input polyline (BAR-array), h is the height of the mesh, hf is the height of the first cell and density is the number of points per unity of length.
    In the 'array' version, it returns a list where B[0] is the list of generated meshes, B[1] is the list of wall boundaries, B[2] is the list of overlap boundaries, B[3] is h, B[4] is density (eventually modified by the mesher).
    In the pyTree version, it returns a list [zones,hs,densities], where zones is a list of zones of a CGNS python tree, containing the blocks, wall boundaries, match and overlap boundaries; hs is the list of heights (modified if necessary), and densities the list of densities (also modified if necessary).

    :param a:  input polyline
    :type  a:  BAR-array
    :param h:  height of the mesh
    :type  h:  float
    :param hf:  first cell size
    :type  hf:  float
    :param density:  number of points per unity of length
    :type  density:  integer
    :return: 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Structured mesh generation from a polyline (array) <Examples/Generator/polyLineMesher.py>`_:

    .. literalinclude:: ../build/Examples/Generator/polyLineMesher.py

    * `Structured mesh generation from a polyline (pyTree) <Examples/Generator/polyLineMesherPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/polyLineMesherPT.py

---------------------------------------

.. py:function:: Generator.PolyC1.polyC1Mesher(A, h, hf, density, splitCrit=10.)

    Generate a 2D mesh around a 2D polyC1 curve where A is a list of i-arrays each representing a C1 curve. All i-arrays put together must represent a polyC1 curve. SplitCrit is a curvature radius triggering split. Other arguments are similar to polyLineMesher. The function return is also similar to polyLineMesher. 

    :param A:  list of 1D curves
    :type  A:  list of arrays
    :param h:  height of the mesh
    :type  h:  float
    :param hf:  first cell size
    :type  hf:  float
    :param density:  number of points per unity of length
    :type  density:  integer
    :param splitCrit:  threshold curvature radius below which the initial curve is split
    :type  splitCrit:  float
    :return: 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Structured mesh generation from a C1 line (array) <Examples/Generator/polyC1Mesher.py>`_:

    .. literalinclude:: ../build/Examples/Generator/polyC1Mesher.py

    * `Structured mesh generation from a C1 line (pyTree) <Examples/Generator/polyC1MesherPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/polyC1MesherPT.py

---------------------------------------

.. py:function:: Generator.pointedHat(a, (x,y,z))

    Create a structured mesh from a curve defined by a i-array and a point. For the pyTree version: if a contains a solution, it is not taken into account in b. 

    :param a:  closed 1D curve
    :type  a:  array
    :param (x,y,z):  coordinates of point
    :type  (x,y,z):  3-tuple of floats
    :return: 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `2D closing mesh generation from a closed curve and a point (array) <Examples/Generator/pointedHat.py>`_:

    .. literalinclude:: ../build/Examples/Generator/pointedHat.py

    * `2D closing mesh generation from a closed curve and a point (pyTree) <Examples/Generator/pointedHatPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/pointedHatPT.py

---------------------------------------

.. py:function:: Generator.stitchedHat(a, (offx,offy,offz), tol=1.e-6, tol2=1.e-5)

    Create a stitched mesh from a curve defined by a i-array. The surface is stitched in the middle. Tol is the accuracy of the search, tol2 is a merging tolerance and offx, offy, off z an optional offset. For the pyTree version: if a contains a solution, it is not taken into account in b. 

    :param a:  closed 1D curve
    :type  a:  array
    :param (offx,offy,offz):  coordinates of offset vector
    :type  (offx,offy,offz):  3-tuple of floats
    :param tol:  accuracy of search
    :type  tol:  float
    :param tol2:  merging tolerance
    :type  tol2:  float
    :return: 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `2D stitched mesh generation from a closed curve (array) <Examples/Generator/stitchedHat.py>`_:

    .. literalinclude:: ../build/Examples/Generator/stitchedHat.py

    * `2D stitched mesh generation from a closed curve (pyTree) <Examples/Generator/stitchedHatPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/stitchedHatPT.py

---------------------------------------

.. py:function:: Generator.surfaceWalk(surfaces, c, dj, constraints=[], niter=0, alphaRef=180., check=0, toldist=1.e-6)

    Surface extrusion starting from a curve, resulting into a surface mesh. dj is the distribution of points in the extrusion direction starting from c, niter the number of smoothing iterations. check=1 means that the extrusion stops at the layer before cells intersect alphaRef is the deviation angle wrt 180 degrees enabling to stop the extrusion before it crosses a sharp edge on the surface. toldist is a tolerance below which points are considered matching. Constraints can be set as 1D zones. 

    :param surfaces:  list of surfaces
    :type  surfaces:  list of arrays
    :param c:  starting curve for the extrusion
    :type  c:  array or pyTree
    :param dj:  distribution of points for the extrusion
    :type  dj:  1D-array
    :param constraints:  1D curves constraining the extrusion
    :type  constraints:  list of arrays
    :param niter:  number of smoothing iterations
    :type  niter:  integer
    :param alphaRef:  deviation angle (in degrees) stopping the extrusion
    :type  alphaRef:  float
    :param check:  activation key for stopping the extrusion (0 or 1)
    :type  check:  integer
    :param toldist:  merging points tolerance
    :type  toldist:  float
    :return: 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `2D mesh extrusion from a curve and walking on a surface (array) <Examples/Generator/surfaceWalk.py>`_:

    .. literalinclude:: ../build/Examples/Generator/surfaceWalk.py

    * `2D mesh extrusion from a curve and walking on a surface (pyTree) <Examples/Generator/surfaceWalkPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/surfaceWalkPT.py

---------------------------------------

.. py:function:: Generator.collarMesh(s1, s2, dj,dk, niterj=100, niterk=100, ext=5, alphaRef=30., type='union', contour=[], constraints1=[], constraints2=[], toldist=1.e-10, topology='overset')

    Create a collar mesh at junction(s) between two surfaces s1 and s2 in union or difference assembly, using a distribution along the surface dj and a distribution in the normal direction to the wall dk. niterj and niterk are the number of smoothing iterations for j and k directions. ext is the extension of the collar mesh for difference assembly. type is the assembly type, and can be 'union' or 'difference'. alphaRef is the deviation angle wrt 180 degrees above which the walk is stopped. contour is the starting contour to create the collar grids, constraints1 and constraints2 are 1D zones defining the curves the collar grid must follow on surfaces s1 and s2 respectively. toldist is the matching point tolerance. Parameter 'topology' can be 'overset' or 'extruded', only useful in case of difference. Topology set to 'overset' results in two overlapping collar grids, whereas it results in a collar grid extruded from the surface grid in the other case. 

    :param s1:  surface
    :type  s1:  array or pyTree
    :param s2:  surface
    :type  s2:  array or pyTree
    :param dj:  distribution of points along surfaces
    :type  dj:  1D-array
    :param dk:  distribution of points in the normal direction
    :type  dk:  1D-array
    :param niterj:  number of smoothing iterations in j direction
    :type  niterj:  integer
    :param niterk:  number of smoothing iterations in k direction
    :type  niterk:  integer
    :param ext:  extension of collar for difference assembly
    :type  ext:  integer
    :param alphaRef:  deviation angle (in degrees) stopping the extrusion
    :type  alphaRef:  float
    :param type:  type of the assembly (union or difference)
    :type  type:  string
    :param contour:  starting curve for the collar creation
    :type  contour:  list of arrays
    :param constraints1:  1D curves constraining the collar on s1 surface
    :type  constraints1:  list of arrays
    :param constraints2:  1D curves constraining the collar on s2 surface
    :type  constraints2:  list of arrays
    :param toldist:  merging points tolerance
    :type  toldist:  float
    :param topology:  choice of collar mesh topology (overset or extruded) in case of difference assembly 
    :type  topology:  string
    :return: 3D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `3D collar mesh between two surfaces (array) <Examples/Generator/collarMesh.py>`_:

    .. literalinclude:: ../build/Examples/Generator/collarMesh.py

    * `3D collar mesh between two surfaces (pyTree) <Examples/Generator/collarMeshPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/collarMeshPT.py


---------------------------------------------------------------------------

Cartesian grid generators
---------------------------------------

.. py:function:: Generator.gencartmb(A, h, Dfar, nlvl)

    Simple Cartesian generator. Create a set of Cartesian grids (B) around a list of body grids (A). Those grids are patched with a ratio of 2. The user controls the number of levels, and the number of points for each level of grid. h is the spatial step on the finest level. Dfar is the maximal distance to the body. nlvl is a list that provides the number of points per level (nlvl[0]: finest grid), except for the finest level. 

    :param A:  body grids
    :type  A:  array/list of arrays or pyTree/list of pyTrees
    :param h:  spatial step in the finest level
    :type  h:  float
    :param Dfar:  maximal distance to the body A
    :type  Dfar:  float
    :param nlvl:  list of number of points per level (except the finest one)
    :type  nlvl:  list of integers
    :return: 2D/3D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Generation of Cartesian mesh refined near body grids (array) <Examples/Generator/gencartmb.py>`_:

    .. literalinclude:: ../build/Examples/Generator/gencartmb.py

    * `Generation of Cartesian mesh refined near body grids (pyTree) <Examples/Generator/gencartmbPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/gencartmbPT.py

---------------------------------------

.. py:function:: Generator.octree(surfs, snearList=[], dfarList=[], dfar=-1., balancing=0, levelMax=1000, ratio=2, octant=None, mode=0)

    Create a QUAD quadtree mesh in 2D or an HEXA octree mesh in 3D starting from a list of bodies and snears. Each parameter snear is the required spatial step of the octree near the corresponding body; 
    the extension of the domaine can be provided by dfar, starting from the global bounding box of all surfaces defined by surfs.
    A list of extensions can be provided in dfarList, in order not to take into account a surface in the computation of the bounding box. 
    It must be set to -1 for the surface that must not be taken into account.
    Parameter balancing=1 means that the octree is balanced, i.e. adjacent elements are at worst twice as big/small; levelMax is the maximum number of levels required. If ratio=2, then a classical octree mesh is built. If ratio=3, a 27-tree mesh is built, in which case the spacing ratio is 3 (and not 2) between two adjacent elements. 
    Parameter balancing enables to balance the octree; balancing=0 means no balancing; balancing=1 means a classical balancing, whereas
    balancing=2 takes also into account elements sharing a common vertex. Paramater mode=1 expands the domain size to get exactly the minimum snear set as input, otherwise the real minimum snear might be slightly lower than expected to comply with the dfar parameter.

    :param surfs:  body grids
    :type  surfs:  list of arrays/pyTrees
    :param snears:  list of spatial step near the corresponding body
    :type  snears:  list of floats
    :param dfar:  maximal distance to the body grids
    :type  dfar:  float
    :param balancing:  activation key for balanced octree (0, 1 or 2)
    :type  balancing:  integer
    :param levelMax:  maximum number of levels
    :type  levelMax:  integer
    :param ratio:  spacing ratio between two adjacent elements
    :type  ratio:  integer
    :param mode:  activation key (0 or 1)
    :type  mode:  integer
    :return: 2D/3D unstructured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Generation of unstructured octree mesh refined near body grids (array) <Examples/Generator/octree.py>`_:

    .. literalinclude:: ../build/Examples/Generator/octree.py

    * `Generation of unstructured octree mesh refined near body grids (pyTree) <Examples/Generator/octreePT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/octreePT.py

---------------------------------------

.. py:function:: Generator.octree2Struct(octree, vmin=15, ext=0, optimized=1, merged=1, AMR=0, sizeMax=1000000)

    Convert an octree or a quadtree mesh into a set of Cartesian grids. Parameter ext is the extension of Cartesian grids in all the directions; vmin can be an integer defining the number of points in each Cartesian grid, or a list of integers, defining the number of points per refinement level. In that case, the first element of the list of vmin defines the finest level. Specifying all the levels is not mandatory. If optimized=1, the ext value is reduced by -1 at overlap borders for the coarsest grid for minimum overlapping. If merged=1, Cartesian grids are merged in order to reduce the number of created grids. If AMR=1, a set of AMR zones are generated. Parameter sizeMax can be used when merging is applied: in that case, the number of points per grid does not exceed sizeMax. Warning: to obtain multigrid blocks, optimized must be set to 0. 

    :param octree:  input unstructured octree grid
    :type  octree:  array or pyTree
    :param vmin:  number of points in all Cartesian grids or list of number of points for each octree level
    :type  vmin:  integer or list of integers
    :param ext:  extension of Cartesian grids (0 = no extension, N = extension of N cells in all direction)
    :type  ext:  integer
    :param optimized:  activation key for optimization of coarsest grid (0 or 1)
    :type  optimized:  integer
    :param merged:  activation key for automatic merging of Cartesian grids
    :type  merged:  integer
    :param AMR:  activation key for AMR generation (0 or 1)
    :type  AMR:  integer
    :param sizeMax:  maximum number of points in Cartesian grids after merging
    :type  sizeMax:  integer
    :return: 2D/3D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Generation of structured Cartesian mesh from an octree grid (array) <Examples/Generator/octree2Struct.py>`_:

    .. literalinclude:: ../build/Examples/Generator/octree2Struct.py

    * `Generation of structured Cartesian mesh from an octree grid (pyTree) <Examples/Generator/octree2StructPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/octree2StructPT.py

---------------------------------------

.. py:function:: Generator.adaptOctree(octree, indicator, balancing=1, ratio=2)

    Adapt an unstructured octree with respect to an indicator field located at element centers. If 'indicator' is strictly positive for an element, then the element must be refined as many times as required by the indicator number. If 'indicator' is strictly negative, the element is coarsened if possible as many times as required by the indicator number. If 'indicator' is 0., the element remains unchanged. balancing=1 means that the octree is balanced after adaptation. If ratio=2, then a classical octree mesh is built. If ratio=3, a 27-tree mesh is built, in which case the spacing ratio is 3 (and not 2) between two adjacent elements. For array interface indicator is an array, for pyTree version, indicator is the name of field stored as a solution located at centers. 
    Exists also as in place version (_adaptOctree) that modifies a and returns None. 


    :param octree:  input unstructured octree grid
    :type  octree:  array or pyTree
    :param indicator:  field of values to indicate where to refine, coarsen or maintain the octree grid
    :type  indicator:  array or variable name in the pyTree
    :param balancing:  activation key for balanced octree (0, 1 or 2, see the definition of octree function for the meaning)
    :type  balancing:  integer
    :param ratio:  spacing ratio between two adjacent elements
    :type  ratio:  integer
    :return: modified reference copy of t
    :rtype: same as input

    *Example of use:*

    * `Adaptation of an octree grid wrt. indicator (array) <Examples/Generator/adaptOctree.py>`_:

    .. literalinclude:: ../build/Examples/Generator/adaptOctree.py

    * `Adaptation of an octree grid wrt. indicator (pyTree) <Examples/Generator/adaptOctreePT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/adaptOctreePT.py

---------------------------------------

.. py:function:: Generator.expandLayer(octree, level=0, corners=0, balancing=0)

    Expand the layer of given level for an octree unstructured mesh. If corners=1, expand also in corners directions.
    Exists also as in place version (_expandLayer) that modifies a and returns None. 

    :param octree:  input unstructured octree grid
    :type  octree:  array or pyTree
    :param level:  level to be expanded (level=0 is the finest)
    :type  level:  integer
    :param corners:  activation key for expansion in corners (0 or 1)
    :type  corners:  integer
    :param balancing:  activation key for balanced octree (0, 1 or 2, see the definition of octree function for the meaning)
    :type  balancing:  integer
    :return: modified reference copy of t
    :rtype: same as input

    *Example of use:*

    * `Expansion of user-specified level in an octree grid (array) <Examples/Generator/expandLayer.py>`_:

    .. literalinclude:: ../build/Examples/Generator/expandLayer.py

    * `Expansion of user-specified level in an octree grid (pyTree) <Examples/Generator/expandLayerPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/expandLayerPT.py

---------------------------------------

.. py:function:: Generator.PyTree.cartRx(X0, H, N, Nb, depth=0, addCellN=False, addBCMatch=False, rank=None, size=None)

    Create a set of regular cartesian grids.
    If depth > 0, an overlap of depth cells is added between grids.
    When used in parallel, the set of grids is distributed. 

    :param X0: first point coordinates
    :type  X0: tuple of 3 floats
    :param H: steps of grids
    :type  H: tuple of 3 floats
    :param N: number of points in each direction for each block
    :type  N: tuple of 3 integers
    :param Nb: number of blocks in each direction
    :type  Nb: tuple of 3 floats
    :param depth: number of overlap cells between cartesian grids
    :type  depth: integer
    :param addCellN: if True and depth > 0, add cellN field
    :type  addCellN: boolean
    :param addBCMatch: if True and depth=0, create BCMatch
    :type  addBCMatch: boolean
    :param rank: current rank when running in parallel
    :type  rank: integer
    :param size: current number of procs when running in parallel
    :type  size: integer
    
    :return: regular cartesian set of grids
    :rtype: list of zones

    *Example of use:*

    * `Generation of a set of regular cartesian grids (pyTree) <Examples/Generator/cartRxPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartRxPT.py

---------------------------------------

.. py:function:: Generator.PyTree.cartRx3(XC0, XC1, HC, XF0, XF1, R, dim=3, rank=None, size=None)

    Create a set of regular and geometric cartesian grids with double steps.
    The mesh is made of a regular cartesian core and border grids are geometric cartesian grids.
    When used in parallel, the set of grids is distributed. 

    :param XC0: first point of cartesian core
    :type  XC0: tuple of 3 floats
    :param XC1: last of of cartesian core
    :type  XC1: tuple of 3 floats
    :param HC: core cartesian steps
    :type  HC: tuple of 3 floats
    :param XF0: first point of geometric region
    :type  XF0: tuple of 3 floats
    :param XF1: last point of geometric region
    :type  XF1: tuple of 3 floats
    :param R: geometric factor in each direction
    :type  R: tuple of 3 floats
    :param dim: space dimension
    :type  dim: 2 or 3
    :param rank: current rank when running in parallel
    :type  rank: integer
    :param size: current number of procs when running in parallel
    :type  size: integer
    
    :return: cartesian set of grids
    :rtype: list of zones

    *Example of use:*

    * `Generation of a set of regular and geometric cartesian grids (pyTree) <Examples/Generator/cartRx3PT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/cartRx3PT.py


---------------------------------------------------------------------------


Operations on meshes
-----------------------

.. py:function:: Generator.close(a, tol=1.e-12, rmOverlappingPts=True, rmOrphanPts=True, rmDuplicatedFaces=True, rmDuplicatedElts=True, rmDegeneratedFaces=True, rmDegeneratedElts=True, indices=None)

    Close a mesh defined by array a. Points of that mesh which are distant of tol maximum to one another are merged. For unstructured zones, connectivity is
    cleaned (removal of unreferenced and duplicated vertices, multiply-defined or totally degenerated faces and elements). If indices=[], the vertex indirection table following mesh cleaning is returned.
    
    Exists also as in place version (_close) that modifies a and returns None. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :param tol:  merging points tolerance
    :type  tol:  float
    :param rmOverlappingPts: whether to remove overlapping points within a radius of tol
    :type rmOverlappingPts: boolean
    :param rmOrphanPts: whether to remove orphan (unreferenced) points
    :type rmOrphanPts: boolean
    :param rmDuplicatedFaces: whether to remove duplicated faces
    :type rmDuplicatedFaces: boolean
    :param rmDuplicatedElts: whether to remove duplicated elements
    :type rmDuplicatedElts: boolean
    :param rmDegeneratedFaces: whether to remove degenerated faces
    :type rmDegeneratedFaces: boolean
    :param rmDegeneratedElts: whether to remove degenerated elements
    :type rmDegeneratedElts: boolean
    :param indices: vertex indirection table following mesh cleaning
    :type indices: [array, list of arrays]
    :return: modified reference copy of t
    :rtype: array or pyTree

    *Example of use:*

    * `Mesh closing (array) <Examples/Generator/close.py>`_:

    .. literalinclude:: ../build/Examples/Generator/close.py

    * `Mesh closing (pyTree) <Examples/Generator/closePT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/closePT.py
        
---------------------------------------

.. py:function:: Generator.zip(a, tol=1.e-12)

    Zip a set of zones defined by a at their borders, if they are at a distance less than tol (mean value of coordinates is set).
    
    Exists also as in place version (_zip) that modifies a and returns None. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :param tol:  merging points tolerance
    :type  tol:  float
    :return: modified reference copy of t
    :rtype: array or pyTree

    *Example of use:*

    * `Border closing (array) <Examples/Generator/zip.py>`_:

    .. literalinclude:: ../build/Examples/Generator/zip.py

    * `Border closing (pyTree) <Examples/Generator/zipPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/zipPT.py

---------------------------------------

.. py:function:: Generator.selectInsideElts(a, curves)

    Select elements of a TRI-array, whose centers are inside the given list of curves, defined by BAR-arrays.

    :param a:  input triangle 2D mesh
    :type  a:  array or pyTree
    :param curves:  list of curves
    :type  curves:  array or list of arrays
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Selection of TRI cells inside a specified curve (array) <Examples/Generator/selectInsideElts.py>`_:

    .. literalinclude:: ../build/Examples/Generator/selectInsideElts.py

    * `Selection of TRI cells inside a specified curve (pyTree) <Examples/Generator/selectInsideEltsPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/selectInsideEltsPT.py

---------------------------------------

.. py:function:: Generator.map(a, distrib, dir)

    Map a distribution on a curve or on a structured surface. Map a i-array distribution in a direction (dir=1,2,3) in a surface or volume mesh.

    :param a:  1D/2D/3D structured mesh
    :type  a:  array or pyTree
    :param distrib:  distribution of points
    :type  distrib:  array
    :param dir:  direction i/j/k for the distribution (dir=1,2,3)
    :type  dir:  integer
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Map distribution (array) <Examples/Generator/map.py>`_:

    .. literalinclude:: ../build/Examples/Generator/map.py

    * `Map distribution (pyTree) <Examples/Generator/mapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/mapPT.py

---------------------------------------

.. py:function:: Generator.mapSplit(a, distrib, splitCrit=100.)

    Split a i-array and map a distribution on the splitted i-array. SplitCrit is the curvature radius triggering split.

    :param a:  1D/2D/3D structured mesh
    :type  a:  array or pyTree
    :param distrib:  distribution of points
    :type  distrib:  array
    :param splitCrit:  curvature radius for array splitting
    :type  splitCrit:  float
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Split and map distribution (array) <Examples/Generator/mapSplit.py>`_:

    .. literalinclude:: ../build/Examples/Generator/mapSplit.py

    * `Split and map distribution (pyTree) <Examples/Generator/mapSplitPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/mapSplitPT.py

---------------------------------------

.. py:function:: Generator.refine(a, power, dir)

    Refine a structured array. The original distribution is kept but the number of points is multiplied by power. Dir is the direction of refinement (1, 2, 3). If dir=0, refine in all directions.
    
    Exists also as in place version (_refine) that modifies a and returns None. 

    :param a:  1D/2D/3D structured mesh
    :type  a:  [array] or [zone]
    :param power:  multiplication factor of number of points
    :type  power:  float
    :param dir:  direction i/j/k for the distribution (dir=0,1,2,3)
    :type  dir:  integer
    :rtype: Identical to a

    *Example of use:*

    * `Structured mesh refinement (array) <Examples/Generator/refine.py>`_:

    .. literalinclude:: ../build/Examples/Generator/refine.py

    * `Structured mesh refinement (pyTree) <Examples/Generator/refinePT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/refinePT.py

---------------------------------------

.. py:function:: Generator.mapCurvature(a, N, power, dir)

    Map a structured array following the curvature. N is the final number of points. Dir is the direction of remeshing (1, 2, 3). 

    :param a:  1D/2D structured mesh
    :type  a:  [array] or [zone]
    :param N:  number of points after new distribution
    :type  N:  integer
    :param power:  refinement factor
    :type  power:  float
    :param dir:  direction i/j/k for the distribution (dir=1,2,3)
    :type  dir:  integer
    :rtype: Identical to a

    *Example of use:*

    * `Map distribution wrt. curvature (array) <Examples/Generator/mapCurvature.py>`_:

    .. literalinclude:: ../build/Examples/Generator/mapCurvature.py

    * `Map distribution wrt. curvature (pyTree) <Examples/Generator/mapCurvaturePT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/mapCurvaturePT.py

---------------------------------------

.. py:function:: Generator.densify(a, h)

    Densify a i-array or a BAR-array with a new discretization step h. Discretization points from the original array are kept.
    
    Exists also as in place version (_densify) that modifies a and returns None. 

    :param a:  1D structured mesh
    :type  a:  [array] or [zone]
    :param h:  new cell size step for the points densification
    :type  h:  float
    :rtype: Identical to a

    *Example of use:*

    * `Curve densification (array) <Examples/Generator/densify.py>`_:

    .. literalinclude:: ../build/Examples/Generator/densify.py

    * `Curve densification (pyTree) <Examples/Generator/densifyPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/densifyPT.py

---------------------------------------

.. py:function:: Generator.grow(a, vector)

    Grow a surface array of one layer. Vector is the node displacement. For the array version, vector is defined by an array. For the PyTree version, vector = ['v1','v2','v3'] where variables 'v1', 'v2', 'v3' are defined as solutions in a, located at nodes. 

    :param a:  2D surface mesh
    :type  a:  array or Zone
    :param vector:  vector of node displacement
    :type  vector:  array or list of 3 variables contained in the solution
    :return: new 3D structured mesh
    :rtype: array or Zone

    *Example of use:*

    * `Extrusion of one layer from a surface mesh (array) <Examples/Generator/grow.py>`_:

    .. literalinclude:: ../build/Examples/Generator/grow.py

    * `Extrusion of one layer from a surface mesh (pyTree) <Examples/Generator/growPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/growPT.py

---------------------------------------

.. py:function:: Generator.stack(a, b=None)

    Can stack 1D meshes with same topology into 2D meshes.
    Can stack 2D meshes with same topology into 3D meshes.
    
    :param a:  a mesh or a list of topologically identical meshes
    :type  a:  array, zone or list of arrays, zones
    :param b:  a mesh
    :type  b:  array, zone    
    :return: a mesh
    :rtype: array or Zone

    *Example of use:*

    * `Mesh generation by stacking two meshes (array) <Examples/Generator/stack.py>`_:

    .. literalinclude:: ../build/Examples/Generator/stack.py

    * `Mesh generation by stacking two meshes (pyTree) <Examples/Generator/stackPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/stackPT.py

---------------------------------------

.. py:function:: Generator.addNormalLayers(a, d, check=0, niter=0)

    Normal extrusion from a surface mesh. d is a 1D distribution providing the height of each layer. If check=1, the extrusion stops before negative volume cells are created. Niter specifies the number of iterations for normals smoothing. 

    :param a:  2D surface mesh
    :type  a:  array or pyTree
    :param d:  distribution of normal extrusion
    :type  d:  1D array
    :param check:  activation key for negative volume criteria (0 or 1)
    :type  check:  integer
    :param niter:  number of iterations for normals smoothing
    :type  niter:  integer
    :return: new 3D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Normal extrusion from a surface mesh (array) <Examples/Generator/addNormalLayers.py>`_:

    .. literalinclude:: ../build/Examples/Generator/addNormalLayers.py

    * `Normal extrusion from a surface mesh (pyTree) <Examples/Generator/addNormalLayersPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/addNormalLayersPT.py

---------------------------------------

.. py:function:: Generator.TTM(a, niter=100)

    Smooth a mesh using elliptic generator. 

    :param a:  2D structured mesh
    :type  a:  array or pyTree
    :param niter:  number of smoothing iterations
    :type  niter:  integer
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `2D structured mesh smoothing (array) <Examples/Generator/TTM.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TTM.py

    * `2D structured mesh smoothing (pyTree) <Examples/Generator/TTMPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/TTMPT.py

---------------------------------------

.. py:function:: Generator.snapFront(a, S, optimized=1)

    Snap a mesh to a surface S. A front must be defined in a by a cellN field. Points of this front are snapped to the surface. If optimized=0, the exterior front cellN=1 is snapped, else if optimized=1 optimized front cellN=1 is snapped, else if optimized=2, front cellN=0 is snapped.
    
    Exists also as in place version (_snapFront) that modifies a and returns None. 


    :param a:  3D mesh
    :type  a:  array or pyTree
    :param S:  surface mesh
    :type  S:  list of zones
    :param optimized:  optimization key (0,1,2)
    :type  optimized:  integer
    :return: new unstructured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Mesh snapping to a surface (array) <Examples/Generator/snapFront.py>`_:

    .. literalinclude:: ../build/Examples/Generator/snapFront.py

    * `Mesh snapping to a surface (pyTree) <Examples/Generator/snapFrontPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/snapFrontPT.py

---------------------------------------

.. py:function:: Generator.snapSharpEdges(a, S, step=None, angle=30.)

    Snap a mesh to a surface S, constrained by sharp edges and corners. if step != None, sharp edges are refined with this step. Sharp Edges are calculated depending on angle. 
    
    Exists also as in place version (_snapSharpEdges) that modifies a and returns None. 

    :param a:  mesh to be deformed
    :type  a:  array or pyTree
    :param S:  surface mesh
    :type  S:  list of zones
    :param step:  step for sharp edges refinement
    :type  step:  float
    :param angle:  angle (in degrees) for sharp edges detection
    :type  angle:  float
    :return: new unstructured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Mesh snapping to sharp edges of a surface (array) <Examples/Generator/snapSharpEdges.py>`_:

    .. literalinclude:: ../build/Examples/Generator/snapSharpEdges.py

    * `Mesh snapping to sharp edges of a surface (pyTree) <Examples/Generator/snapSharpEdgesPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/snapSharpEdgesPT.py


---------------------------------------------------------------------------

Operation on surface meshes
---------------------------------------

.. py:function:: Generator.fittingPlaster(a, bumpFactor=0.)

    Create a structured surface patch that fits a curve a. BumpFactor controls the curvature of the patch. 

    :param a:  curve to fit
    :type  a:  array
    :param bumpFactor:  amplitude of the bump
    :type  bumpFactor:  float
    :return: created 2D structured mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Creation of patch that fits a curve (array) <Examples/Generator/fittingPlaster.py>`_:

    .. literalinclude:: ../build/Examples/Generator/fittingPlaster.py

    * `Creation of patch that fits a curve (pyTree) <Examples/Generator/fittingPlasterPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/fittingPlasterPT.py

---------------------------------------

.. py:function:: Generator.gapfixer(a, c, hardPoints=None, refine=1)

    Fill a surface defined by a BAR contour a drawn on a surface c. You can force the generated mesh to pass 
    through given hardPoints (NODES). If refine=0, no inside points are added. 

    :param a:  contour of the gap
    :type  a:  BAR array
    :param c:  surface to be filled
    :type  c:  array
    :param hardPoints:  mesh containing nodes to be enforced
    :type  hardPoints:  array or list of arrays
    :param refine:  activation key for including points in the gap mesh (0 or 1)
    :type  refine:  integer
    :return: TRI surface mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Surface gap filling from a contour (array) <Examples/Generator/gapfixer.py>`_:

    .. literalinclude:: ../build/Examples/Generator/gapfixer.py

    * `Surface gap filling from a contour (pyTree) <Examples/Generator/gapfixerPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/gapfixerPT.py

---------------------------------------

.. py:function:: Generator.gapsmanager(A, mode=0, coplanar=0)

    Fill multiple gaps in a set of surface components A. Also, eliminate overlap regions between components if any. Normals for all patches must be pointed outwards. Set mode=0 for nodal mesh, 1 for center mesh, and 2 otherwise. Set coplanar=1 if all components are lying on a same plane. 

    :param A:  surface mesh with gaps
    :type  A:  array or pyTree
    :param mode:  key for grid location (0 = nodes, 1 = centers, 2 = others)
    :type  mode:  integer
    :param coplanar:  activation key for coplanar components of A (0 or 1)
    :type  coplanar:  integer
    :return: new surface mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Surface gaps filling (array) <Examples/Generator/gapsmanager.py>`_:

    .. literalinclude:: ../build/Examples/Generator/gapsmanager.py

    * `Surface gaps filling (pyTree) <Examples/Generator/gapsmanagerPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/gapsmanagerPT.py


---------------------------------------

.. py:function:: Generator.mmgs(A, ridgeAngle=45., hmin=0., hmax=0., hausd=0.01, grow=1.1, anisotropy=0, optim=0, fixedConstraints=[], sizeConstraints=[])

    Refine a TRI surface mesh using MMGs.

    :param A:  surface TRI mesh
    :type  A:  [array, list of arrays] or [zone, pyTree]
    :param ridgeAngle: Angle between adjacent cells that MMGs consider to be a ridge (degrees) 
    :type  ridgeAngle: double
    :param hmin: the minimum mesh step in final mesh 
    :type  hmin: double
    :param hmax: the maximum mesh step in final mesh 
    :type  hmax: double
    :param hausd: the maximum chordal deviation in final mesh from initial mesh 
    :type  hausd: double
    :param grow: the maximum difference of steps between two adjacent cells 
    :type  grow: double
    :param optim: if 1, only optimize mesh keeping the same number of points 
    :type  optim: int
    :param fixedConstraints: curves or surface identifying points that must be in output mesh 
    :type  fixedConstraints: [list of arrays] or [list of zones]
    :param sizeConstraints: curves or surface defining sizemap 
    :type  sizeConstraints: [list of arrays] or [list of zones]
    :return: remeshed surface
    :rtype: identical to input

    *Example of use:*

    * `Surface mesh remeshing (array) <Examples/Generator/mmgs.py>`_:

    .. literalinclude:: ../build/Examples/Generator/mmgs.py

    * `Surface mesh remeshing (pyTree) <Examples/Generator/mmgsPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/mmgsPT.py


---------------------------------------------------------------------------

Information on generated meshes
---------------------------------

.. ---------------------------------------

.. .. py:function:: Generator.check(a)

    Check regularity, orthogonality for a mesh defined by an array. 

    :param a:  mesh to be checked
    :type  a:  array or pyTree
    :return: new surface mesh
    :rtype: array or pyTree

    *Example of use:*

    * `Mesh check (array) <Examples/Generator/check.py>`_:

    .. literalinclude:: ../build/Examples/Generator/check.py

    * `Mesh check (pyTree) <Examples/Generator/checkPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/checkPT.py

---------------------------------------

.. py:function:: Generator.barycenter(a, weight=None)

    Return the barycenter of a, with optional weight. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :param weight:  name of the weight variable in a
    :type  weight:  string
    :return: coordinates of the barycenter
    :rtype: 3-list of floats

    *Example of use:*

    * `Computation of the mesh barycenter (array) <Examples/Generator/barycenter.py>`_:

    .. literalinclude:: ../build/Examples/Generator/barycenter.py

    * `Computation of the mesh barycenter (pyTree) <Examples/Generator/barycenterPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/barycenterPT.py

---------------------------------------

.. py:function:: Generator.bbox(a)

    Return the bounding box [xmin, ymin, zmin, xmax, ymax, zmax] of a. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :return: coordinates of the bounding box
    :rtype: 6-list of floats

    *Example of use:*

    * `Computation of the mesh bounding box (array) <Examples/Generator/bbox.py>`_:

    .. literalinclude:: ../build/Examples/Generator/bbox.py

    * `Computation of the mesh bounding box (pyTree) <Examples/Generator/bboxPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/bboxPT.py

---------------------------------------

.. py:function:: Generator.bboxOfCells(a)

    Return the bounding box of each cell of a. The bounding box field is located at centers of cells. 
    
    Exists also as in place version (_bboxOfCells) that modifies a and returns None. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of the cell bounding boxes (array) <Examples/Generator/bboxOfCells.py>`_:

    .. literalinclude:: ../build/Examples/Generator/bboxOfCells.py

    * `Computation of the cell bounding boxes (pyTree) <Examples/Generator/bboxOfCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/bboxOfCellsPT.py

---------------------------------------

.. py:function:: Generator.BB(a, method='AABB', weighting=0, tol=0.)

    Return the bounding box of a as an array or a zone. If method is 'AABB', then it computes the Axis-Aligned Bounding-Box, if method is 'OBB' then it computes the Oriented Bounding-Box. The argument weighting may be 0, and the OBB is computed using a Cloud-Point approach, or 1, and it is computed using a Surface-Weighting approach. If weighting=1, then the provided array must be a surface composed of triangles. 
    
    Exists also as in place version (_BB) that modifies a and returns None. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :param method:  choice between axis-aligned or oriented bounding box
    :type  method:  string
    :param weighting:  activation key for surface weighting approach
    :type  weighting:  integer
    :param tol: extension of bounding box in all the directions
    :type tol: float
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Bounding box generation (array) <Examples/Generator/BB.py>`_:

    .. literalinclude:: ../build/Examples/Generator/BB.py

    * `Bounding box generation (pyTree) <Examples/Generator/BBPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/BBPT.py

---------------------------------------

.. py:function:: Generator.CEBBIntersection(a1, a2, tol=1.e-10)

    Test the Cartesian Elements Bounding Box (CEBB) intersection between a1 and a2. Tolerance is a float given by tol. Return 0 if no intersection, 1 otherwise. 

    :param a1:  input mesh
    :type  a1:  array or pyTree
    :param a2:  input mesh
    :type  a2:  array or pyTree
    :param tol:  tolerance of intersection
    :type  tol:  float
    :return: 0 if no intersection, 1 otherwise
    :rtype: integer

    *Example of use:*

    * `Intersection by Cartesian elements bounding box between two meshes (array) <Examples/Generator/CEBBIntersection.py>`_:

    .. literalinclude:: ../build/Examples/Generator/CEBBIntersection.py

    * `Intersection by Cartesian elements bounding box between two meshes (pyTree) <Examples/Generator/CEBBIntersectionPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/CEBBIntersectionPT.py

---------------------------------------

.. py:function:: Generator.bboxIntersection(a1, a2, tol=1.e-6, isBB=False, method='AABB')

    Test if a1 and a2 intersects. Three options are available: method='AABB' (intersection between two Axis-Aligned Bounding Boxes, by default); method='OBB' (intersection between two Oriented Bounding Boxes, the most general case); method='AABBOBB' (intersection between an AABB -a1- and an OBB -a2-).; If a1 and a2 are directly the corresponding bounding boxes, the user may switch isBB=True in order to avoid recalculating them. Return 0 if no intersection, 1 otherwise. 
    
    Exists also as in place version (_bboxIntersection) that modifies a1 and returns None. 

    :param a1:  input mesh
    :type  a1:  array or pyTree
    :param a2:  input mesh
    :type  a2:  array or pyTree
    :param tol:  tolerance of intersection
    :type  tol:  float
    :param isBB:  activation key if meshes already are bounding boxes
    :type  isBB:  boolean
    :param method:  intersection method
    :type  method:  string
    :return: 0 if no intersection, 1 otherwise
    :rtype: integer

    *Example of use:*

    * `Intersection by bounding box between two meshes (array) <Examples/Generator/bboxIntersection.py>`_:

    .. literalinclude:: ../build/Examples/Generator/bboxIntersection.py

    * `Intersection by bounding box between two meshes (pyTree) <Examples/Generator/bboxIntersectionPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/bboxIntersectionPT.py

---------------------------------------

.. py:function:: Generator.checkPointInCEBB(a, (x,y,z))

    Test if a given point is in the CEBB of a. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :param (x,y,z):  coordinates of point
    :type  (x,y,z):  3-tuple of floats
    :return: 0 if point is not in the CEBB of a, 1 otherwise
    :rtype: integer

    *Example of use:*

    * `Detection of point location in the bounding box of a mesh (array) <Examples/Generator/checkPointInCEBB.py>`_:

    .. literalinclude:: ../build/Examples/Generator/checkPointInCEBB.py

    * `Detection of point location in the bounding box of a mesh (pyTree) <Examples/Generator/checkPointInCEBBPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/checkPointInCEBBPT.py

---------------------------------------

.. py:function:: Generator.getVolumeMap(a, method=0, tol=1e-12)

    Return the volume field of an array. Volume is located at centers of cells. 
    
    Exists also as in place version (_getVolumeMap) that modifies a and returns None. 

    :param a:  input volume or surface mesh
    :type  a:  array or pyTree
    :param method: method of volumes computation (0 or 1). method = 1 is usually more robust on NGons
    :type method: int
    :param tol: tolerance used within internal routines when method = 1
    :type tol: float
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of cells volume (array) <Examples/Generator/getVolumeMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getVolumeMap.py

    * `Computation of cells volume (pyTree) <Examples/Generator/getVolumeMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getVolumeMapPT.py

---------------------------------------

.. py:function:: Generator.getNormalMap(a)

    Return the surface normals field of a surface array. It is located at centers of cells. 
    
    Exists also as in place version (_getNormalMap) that modifies a and returns None. 

    :param a:  input surface mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of surface normals (array) <Examples/Generator/getNormalMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getNormalMap.py

    * `Computation of surface normals (pyTree) <Examples/Generator/getNormalMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getNormalMapPT.py

---------------------------------------

.. py:function:: Generator.getSmoothNormalMap(a, niter=2,eps=0.4)

    Return the smoothed surface normals field of a surface array, located at nodes. niter is the number of smoothing operations, and eps is a smoothing weight. 
    
    Exists also as in place version (_getSmoothNormalMap) that modifies a and returns None. 

    :param a:  input surface mesh
    :type  a:  array or pyTree
    :param niter:  smoothing iterations number
    :type  niter:  integer
    :param eps:  smoothing weight
    :type  eps:  float
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of surface smoothed normals (array) <Examples/Generator/getSmoothNormalMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getSmoothNormalMap.py

    * `Computation of surface smoothednormals (pyTree) <Examples/Generator/getSmoothNormalMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getSmoothNormalMapPT.py

---------------------------------------

.. py:function:: Generator.getOrthogonalityMap(a)

    Return the orthogonality map of an array. The orthogonality map corresponds to the maximum deviation of all dihedral angles of an element. The orthogonality map is expressed in degree and located at centers. 
    
    Exists also as in place version (_getOrthogonalityMap) that modifies a and returns None. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of cells orthogonality (array) <Examples/Generator/getOrthogonalityMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getOrthogonalityMap.py

    * `Computation of cells orthogonality (pyTree) <Examples/Generator/getOrthogonalityMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getOrthogonalityMapPT.py

---------------------------------------

.. py:function:: Generator.getRegularityMap(a)

    Return the regularity map of an array. The regularity map corresponds to the maximum deviation of 
    the volume ratio of an element and all its neigbouring cells. The regularity map is located 
    at centers in a field "regularity". 
    
    Exists also as in place version (_getRegularityMap) that modifies a and returns None. 
    Exists also as parallel distributed version (G.Mpi.getRegularityMap).

    :param a: input mesh
    :type  a: array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of cells regularity with respect to volumes (array) <Examples/Generator/getRegularityMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getRegularityMap.py

    * `Computation of cells regularity with respect to volumes (pyTree) <Examples/Generator/getRegularityMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getRegularityMapPT.py

---------------------------------------

.. py:function:: Generator.getAngleRegularityMap(a)

    Return the angle regularity map of an array. 
    The angle regularity map corresponds to the maximum angle difference between an element and its 
    neigbouring cells. The angle regularity map is expressed in degree and located at centers in 
    a field "regularityAngle".
    
    Exists also as in place version (_getAngleRegularityMap) that modifies a and returns None. 
    Exists also as parallel distributed version (G.Mpi.getAngleRegularityMap).

    :param a: input mesh
    :type  a: array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of cells regularity with respect to angles (array) <Examples/Generator/getAngleRegularityMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getAngleRegularityMap.py

    * `Computation of cells regularity with respect to angles (pyTree) <Examples/Generator/getAngleRegularityMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getAngleRegularityMapPT.py

---------------------------------------

.. py:function:: Generator.getTriQualityMap(a)

    Return the quality map of a TRI array. The triangle quality is a value between 0. (degenerated triangle) and 1. (equilateral triangle). The quality map is located at centers. 
    
    Exists also as in place version (_getTriQualityMap) that modifies a and returns None. 

    :param a:  input surface TRI mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of triangles quality (array) <Examples/Generator/getTriQualityMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getTriQualityMap.py

    * `Computation of triangles quality (pyTree) <Examples/Generator/getTriQualityMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getTriQualityMapPT.py

---------------------------------------

.. py:function:: Generator.getCellPlanarity(a)

    Return a measure of cell planarity for each cell. It is located at centers of cells. 
    
    Exists also as in place version (_getCellPlanarity) that modifies a and returns None. 

    :param a:  input surface mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of cells planarity (array) <Examples/Generator/getCellPlanarity.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getCellPlanarity.py

    * `Computation of cells planarity (pyTree) <Examples/Generator/getCellPlanarityPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getCellPlanarityPT.py

---------------------------------------

.. py:function:: Generator.getCircumCircleMap(a)

    Return the map of circum circle radius of any cell of a 'TRI' array. 
    
    Exists also as in place version (_getCircumCircleMap) that modifies a and returns None. 

    :param a:  input surface mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of cells circumscribed circle radius (array) <Examples/Generator/getCircumCircleMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getCircumCircleMap.py

    * `Computation of cells circumscribed circle radius (pyTree) <Examples/Generator/getCircumCircleMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getCircumCircleMapPT.py

---------------------------------------

.. py:function:: Generator.getInCircleMap(a)

    Return the map of inscribed circle radius of any cell of a 'TRI' array. 
    
    Exists also as in place version (_getInCircleMap) that modifies a and returns None. 

    :param a:  input surface mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of cells inscribed circle radius (array) <Examples/Generator/getInCircleMap.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getInCircleMap.py

    * `Computation of cells inscribed circle radius (pyTree) <Examples/Generator/getInCircleMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getInCircleMapPT.py

---------------------------------------

.. py:function:: Generator.getEdgeRatio(a)

    Return the ratio between the longest and the smallest edges of a cell. 
    
    Exists also as in place version (_getEdgeRatio) that modifies a and returns None. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of maximum edge ratio of cells (array) <Examples/Generator/getEdgeRatio.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getEdgeRatio.py

    * `Computation of maximum edge ratio of cells (pyTree) <Examples/Generator/getEdgeRatioPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getEdgeRatioPT.py

---------------------------------------

.. py:function:: Generator.getMaxLength(a)

    Return the length of the longer edge of each cell. 
    
    Exists also as in place version (_getMaxLength) that modifies a and returns None. 

    :param a:  input mesh
    :type  a:  array or pyTree
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Computation of maximum edge length of cells (array) <Examples/Generator/getMaxLength.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getMaxLength.py

    * `Computation of maximum edge length of cells (pyTree) <Examples/Generator/getMaxLengthPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/getMaxLengthPT.py
    

---------------------------------------

.. py:function:: Generator.checkMesh(m, critVol=0., critOrtho=15., critReg=0.1, critAngReg=15., addGC=False, verbose=0)

    Return global information about a given mesh based on four unit tests. Two local tests are used (getVolumeMap and getOrthogonalityMap), as well as two global tests (getRegularityMap and getAngleRegularityMap) to assess the mesh quality.
    
    :param m:  input mesh
    :type  m:  array or pyTree
    :param critVol,...: critical values for each unit test. All critical values are used to check if a value is above a given threshold, except for the first unit test (getVolumeMap) for which we usually check for negative volume cells in the mesh
    :type  critVol,...:  floats
    :param addGC: if True, add one layer of ghost cells to the pyTree. Useful for global tests (getAngleRegularityMap and getRegularityMap) in the case of multiple zones
    :type  addGC:  Boolean
    :param verbose: if 2, print mesh info per sub-zone (local min, max, mean values and percentage of critical cells per test). If 1, only print mesh info when critical cells are encountered. If 0, do not print anything
    :type  verbose:  int
    :return: global min, max and mean values as well as the number of critical cells obtained for each test 
    :rtype: dictionary

    *Example of use:*

    * `Get global mesh info (array) <Examples/Generator/checkMesh.py>`_:

    .. literalinclude:: ../build/Examples/Generator/checkMesh.py

    * `Get global mesh info (pyTree) <Examples/Generator/checkMeshPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/checkMeshPT.py



---------------------------------------------------------------------------

Operations on distributions
---------------------------------------

.. py:function:: Generator.enforceX(a, x0, enforcedh, (supp,add))

    Enforce a region around a line x=x0. The size of the cell around the line is enforcedh. "supp" points are suppressed from the starting distribution on the left and right side. "add" points are added on the left and add points are added on the right. Add exactely add points.
    Adjust add in order to have a monotonic distribution with: Generator.enforceX(a, x0, enforcedh, supp, add).
    Exists also for Y and Z directions: Generator.enforceY, Generator.enforceZ. 

    :param a:  input structured mesh
    :type  a:  array or pyTree
    :param x0:  X-coordinate for refinement
    :type  x0:  float
    :param enforcedh:  cell size near refinement
    :type  enforcedh:  float
    :param supp:  number of suppressed points
    :type  supp:  integer
    :param add:  number of added points
    :type  add:  integer
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Structured mesh refinement around a line x=x0 (array) <Examples/Generator/enforceX.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforceX.py

    * `Structured mesh refinement around a line x=x0 (pyTree) <Examples/Generator/enforceXPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforceXPT.py

---------------------------------------

.. py:function:: Generator.enforceMoinsX(a, enforcedh, (supp,add))

    Same as before but with a one sided distribution (left). This can be usefull to create a boundary layer distribution in an Euler mesh.
    Adjust add in order to have a monotonic distribution with: Generator.enforceMoinsX(a, enforcedh, supp, add).
    Exists also for Y and Z directions: Generator.enforceMoinsY, Generator.enforceMoinsZ. 

    :param a:  input structured mesh
    :type  a:  array or pyTree
    :param enforcedh:  cell size near refinement
    :type  enforcedh:  float
    :param supp:  number of suppressed points
    :type  supp:  integer
    :param add:  number of added points
    :type  add:  integer
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Structured mesh refinement at left side (array) <Examples/Generator/enforceMoinsX.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforceMoinsX.py

    * `Structured mesh refinement at left side (pyTree) <Examples/Generator/enforceMoinsXPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforceMoinsXPT.py

---------------------------------------

.. py:function:: Generator.enforcePlusX(a, enforcedh, (supp,add))

    Same as before but with a one sided distribution (right).
    Adjust add in order to have a monotonic distribution with: Generator.enforcePlusX(a, x0, enforcedh, supp, add).
    Exists also for Y and Z directions: Generator.enforceMoinsY, Generator.enforcePlusZ. 

    :param a:  input structured mesh
    :type  a:  array or pyTree
    :param enforcedh:  cell size near refinement
    :type  enforcedh:  float
    :param supp:  number of suppressed points
    :type  supp:  integer
    :param add:  number of added points
    :type  add:  integer
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Structured mesh refinement at right side (array) <Examples/Generator/enforcePlusX.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforcePlusX.py

    * `Structured mesh refinement at right side (pyTree) <Examples/Generator/enforcePlusXPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforcePlusXPT.py

---------------------------------------

.. py:function:: Generator.enforceLine(a, line, enforcedh, (supp,add))

    Enforce a curvilinear line defined by the array line in a distribution defined by the array a. 

    :param a:  input 2D distribution
    :type  a:  array or pyTree
    :param line:  line
    :type  line:  array
    :param enforcedh:  cell size near refinement
    :type  enforcedh:  float
    :param supp:  number of suppressed points
    :type  supp:  integer
    :param add:  number of added points
    :type  add:  integer
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Distribution refinement around a line (array) <Examples/Generator/enforceLine.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforceLine.py

    * `Distribution refinement around a line (pyTree) <Examples/Generator/enforceLinePT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforceLinePT.py

---------------------------------------

.. py:function:: Generator.enforcePoint(a, x0)

    Enforce a point in the distribution. The index of enforced point is returned. 

    :param a:  input 2D distribution
    :type  a:  array or pyTree
    :param x0:  I-location of the refinement point
    :type  x0:  float
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Distribution refinement at a point (array) <Examples/Generator/enforcePoint.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforcePoint.py

    * `Distribution refinement at a point (pyTree) <Examples/Generator/enforcePointPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforcePointPT.py

---------------------------------------

.. py:function:: Generator.enforceCurvature(a, curve, power=0.5)

    Enforce the curvature of an i-curve in a distribution defined by a. Power reflects the power of stretching. 

    :param a:  input 2D distribution
    :type  a:  array or pyTree
    :param curve:  reference curve for curvature
    :type  curve:  array
    :param power:  stretching ratio
    :type  power:  float
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Distribution refinement wrt. a curve curvature (array) <Examples/Generator/enforceCurvature.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforceCurvature.py

    * `Distribution refinement wrt. a curve curvature (pyTree) <Examples/Generator/enforceCurvaturePT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/enforceCurvaturePT.py

---------------------------------------

.. py:function:: Generator.addPointInDistribution(a, ind)

    Add a point in a distribution at index ind. 

    :param a:  input distribution
    :type  a:  array or pyTree
    :param ind:  I-index of inserted point
    :type  ind:  integer
    :return: modified reference copy of a
    :rtype: array or pyTree

    *Example of use:*

    * `Point insertion in a distribution (array) <Examples/Generator/addPointInDistribution.py>`_:

    .. literalinclude:: ../build/Examples/Generator/addPointInDistribution.py

    * `Point insertion in a distribution (pyTree) <Examples/Generator/addPointInDistributionPT.py>`_:

    .. literalinclude:: ../build/Examples/Generator/addPointInDistributionPT.py


---------------------------------------------------------------------------


.. toctree::
   :maxdepth: 2   


Index
#######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

