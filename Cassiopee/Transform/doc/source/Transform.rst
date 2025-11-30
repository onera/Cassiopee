.. Transform documentation master file

:tocdepth: 2


Transform: mesh transformation module
=========================================

Preamble
########

Transform module performs simple transformations of meshes. 
It works on arrays (as defined in Converter documentation) or on CGNS/Python trees (pyTrees), if they provide grid coordinates. 
In the pyTree version, flow solution nodes and boundary conditions and grid connectivity are preserved if 
possible. In particular, splitting a mesh does not maintain the grid connectivity. 

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

To use the module with the Converter.array interface::

   import Transform as T

To use the module with the CGNS/Python interface::

    import Transform.PyTree as T


.. py:module:: Transform


List of functions
##################

**-- Basic operations**

.. autosummary::
   :nosignatures:

    Transform.oneovern
    Transform.reorder
    Transform.reorderAll 
    Transform.makeCartesianXYZ
    Transform.makeDirect
    Transform.addkplane
    Transform.collapse
    Transform.patch


**-- Mesh positioning**

.. autosummary::
   :nosignatures:
   
   Transform.rotate
   Transform.translate

**-- Mesh transformation**

.. autosummary::
   :nosignatures:
    
    Transform.cart2Cyl
    Transform.homothety
    Transform.contract
    Transform.scale
    Transform.symetrize
    Transform.perturbate
    Transform.smooth
    Transform.smoothField
    Transform.dual
    Transform.breakElements

**-- Mesh splitting and merging**

.. autosummary::
   :nosignatures:

    Transform.subzone
    Transform.join
    Transform.merge
    Transform.mergeCart
    Transform.splitNParts
    Transform.splitSize
    Transform.splitCurvatureAngle
    Transform.splitCurvatureRadius
    Transform.splitConnexity
    Transform.splitMultiplePts
    Transform.PyTree.splitFullMatch
    Transform.splitSharpEdges
    Transform.splitTBranches
    Transform.splitManifold
    Transform.splitBAR
    Transform.splitTRI

**-- Mesh deformation**

.. autosummary::
   :nosignatures:

    Transform.deform
    Transform.deformNormals
    Transform.deformPoint
    Transform.deformMesh

**-- Mesh projections**

.. autosummary::
   :nosignatures:
    
    Transform.projectAllDirs
    Transform.projectDir
    Transform.projectOrtho
    Transform.projectOrthoSmooth
    Transform.projectRay


Contents
#########


Basic operations
---------------------------------------


.. py:function:: Transform.oneovern(a, (Ni,Nj,Nk))

    .. A1.O0.D1

    Extract every Ni,Nj,Nk points in the three directions of a structured mesh a.
    
    Exists also as an in-place version (_oneovern) which modifies a and returns None.
   
    :param a: input data
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param (Ni,Nj,Nk): period of extraction in the three directions
    :type (Ni,Nj,Nk): 3-tuple of integers
    :return: a coarsened structured mesh
    :rtype: identical to input

    *Example of use:*

    * `Extract every two points in a structured mesh (array) <Examples/Transform/oneovern.py>`_:

    .. literalinclude:: ../build/Examples/Transform/oneovern.py

    * `Extract every two points in a structured mesh (pyTree) <Examples/Transform/oneovernPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/oneovernPT.py

---------------------------------------

.. py:function:: Transform.reorder(a, order) 

    .. A1.O0.D1
      
    For a structured grid, change the (i,j,k) ordering of a. 
    If you set order=(i2,j2,k2) for a (i,j,k) mesh, going along i2 direction of the resulting mesh
    will be equivalent to go along i direction of initial mesh.
   
    The transformation can be equivalently described by a matrix M, filled with
    a single non-zero value per line and column (equal to -1 or 1).
    Then, order=(orderi, orderj, orderk) means: 
    M[abs(orderi),1]=sign(orderi); M[abs(orderj),2]=sign(orderj); M[abs(orderk),3]=sign(orderk).

    For a 2D unstructured grid (TRI, QUAD, 2D NGON), order the element nodes such that all normals are oriented towards the 
    same direction. If order is set to (1,), all elements are oriented as element 0. If order is (-1,), all elements are oriented
    in the opposite sense of element 0.

    For a 3D unstructured grid (TETRA, PYRA, PENTA, HEXA), element normals must be pointing outward by definition such that the
    argument order can be left as None.

    3D NGON grids are not supported yet.
   
    Exists also as an in-place version (_reorder) which modifies a and returns None.
   
    :param a: initial mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param order: integers specifying transformation
    :type order: 3-tuple of signed integers or a tuple of a single 1 or -1 or None
    :return: a reoriented mesh
    :rtype: identical to input

    *Example of use:* 

    * `Reorder a mesh (array) <Examples/Transform/reorder.py>`_:

    .. literalinclude:: ../build/Examples/Transform/reorder.py

    * `Reorder a mesh (pyTree) <Examples/Transform/reorderPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/reorderPT.py

---------------------------------------

.. py:function:: Transform.reorderAll(a, dir=1)

    .. A1.O0.D1
    
    Order a set of surface grids a such that their normals points in the same direction. 
    All the grids in a must be of same nature (structured or unstructured). 
    Orientation of the first grid in the list is used to reorder the other grids. If dir=-1, the orientation 
    is the opposite direction of the normals of the first grid.
   
    In case of unstructured grids, reorientation is guaranteed to be outward (dir=1) if they 
    represent a closed volume.

    Exists also as an in-place version (_reorderAll) which modifies a and returns None.
   
    :param a: initial set of surface grids
    :type a: [list of arrays] or [list of zones]
    :param dir: 1 (default), -1 for a reversed orientation
    :type dir: signed integer
    :return: a reoriented mesh
    :rtype: identical to input

    *Example of use:*

    * `Reorder a set of grids (array) <Examples/Transform/reorderAll.py>`_:

    .. literalinclude:: ../build/Examples/Transform/reorderAll.py

    * `Reorder a set of grids (pyTree) <Examples/Transform/reorderAllPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/reorderAllPT.py

---------------------------------------

.. py:function:: Transform.makeCartesianXYZ(a)

    .. A1.O0.D1. Est il necessaire de documenter cette fonction?
    
    Reorder a structured Cartesian mesh in order to get i,j,k aligned with X,Y,Z respectively.
    Exists also as an in-place version (_makeCartesianXYZ) which modifies a and returns None.
   
    :param a: Cartesian mesh with (i,j,k) not ordered following (X,Y,Z)
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :return: a Cartesian mesh such that direction i is aligned with (0X), ...
    :rtype: identical to input

    *Example of use:*

    * `Align a Cartesian mesh with XYZ (array) <Examples/Transform/makeCartesianXYZ.py>`_:

    .. literalinclude:: ../build/Examples/Transform/makeCartesianXYZ.py

    * `Align a Cartesian mesh with XYZ (pyTree) <Examples/Transform/makeCartesianXYZPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/makeCartesianXYZPT.py

---------------------------------------

.. py:function:: Transform.makeDirect(a)

    .. A1.O0.D1
    
    Reorder an indirect structured mesh to get a direct mesh.
    Exists also as an in-place version (_makeDirect) which modifies a and returns None.
   
    :param a: structured mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :return: a direct structured mesh
    :rtype: identical to input

    *Example of use:*

    * `Make a mesh direct (array) <Examples/Transform/makeDirect.py>`_:

    .. literalinclude:: ../build/Examples/Transform/makeDirect.py

    * `Make a mesh direct (pyTree) <Examples/Transform/makeDirectPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/makeDirectPT.py

---------------------------------------

.. py:function:: Transform.addkplane(a, N=1)

    .. A1.O0.D1
    
    Add one or more planes at constant heights in z: z0+1, ...,z0+N.
    Exists also as an in-place version (_addkplane) which modifies a and returns None.
   
    :param a: any mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param N: number of layers in the k direction to be added
    :type N: integer
    :return: expanded mesh
    :rtype: identical to input

    *Example of use:*

    * `Add a k-plane in direction (Oz) (array) <Examples/Transform/addkplane.py>`_:

    .. literalinclude:: ../build/Examples/Transform/addkplane.py

    * `Add a k-plane in direction (Oz) (pyTree) <Examples/Transform/addkplanePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/addkplanePT.py

---------------------------------------

.. py:function:: Transform.collapse(a)

    .. A1.O0.D1
    
    Collapse the smallest edges of each element of a triangular mesh. Return
    each element as a BAR. 
    Exists also as an in-place version (_collapse) which modifies a and returns None.
   
    :param a: a TRI mesh
    :type a: [array list of arrays] or [zone, list of zones, base, pyTree] 
    :return: a BAR mesh
    :rtype: identical to input

    *Example of use:*

    * `Collapse smallest edge of a triangular mesh (array) <Examples/Transform/collapse.py>`_:

    .. literalinclude:: ../build/Examples/Transform/collapse.py

    * `Collapse smallest edge of a triangular mesh (pyTree) <Examples/Transform/collapsePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/collapsePT.py

---------------------------------------

.. py:function:: Transform.patch(a, b, position=None, nodes=None, order=None)

    .. A1.O0.D1. Cette fonction devrait etre optimum pour patcher des solutions.
    
    For a structured mesh a, patch (replace) a structured mesh b from starting point position=(i,j,k) of a. 
    Coordinates and fields are replaced.

    For an unstructured mesh a, patch an unstructured mesh (of same type) by replacing the nodes of indices nodes.
   
    Exists also as an in-place version (_patch) which modifies a and returns None.
   
    :param a: initial mesh
    :type a: array or zone
    :param b: patch mesh
    :type b: array or zone
    :param position: indices starting from 1 of the starting node to be replaced in a
    :type position: 3-tuple of integers 
    :param nodes: list of nodes of the unstructured mesh a to be replaced
    :type nodes: numpy array of integers (starting from 1)
    :param order: 3-tuple of integers indicating order of b relative to a (see reorder)
    :type order: None or 3-tuple of integers
    :return: a modified zone
    :rtype: an array or a zone

    *Example of use:*

    * `Patch a mesh in another one (array) <Examples/Transform/patch.py>`_:

    .. literalinclude:: ../build/Examples/Transform/patch.py

    * `Patch a mesh in another one (pyTree) <Examples/Transform/patchPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/patchPT.py


---------------------------------------------------------------------------

Mesh positioning
------------------------


.. py:function:: Transform.rotate(a, C, arg1, arg2=None, vectors=[])

    .. A2.O1.D1
    
    Rotate a mesh. Rotation can be also applied on some vector fields (e.g. velocity and momentum). 
    If the vector field is located at cell centers, then each vector component name must be 
    prefixed by 'centers:'. For instance:
    vectors = [['VelocityX', 'VelocityY', 'VelocityZ']] -or-
    vectors = [['centers:VelocityX', 'centers:VelocityY', 'centers:VelocityZ']] -or-
    vectors = [['centers:MomentumX', 'centers:MomentumY', 'centers:MomentumZ']]
    

    Exists also as an in-place version (_rotate) which modifies a and returns None.

    Rotation parameters can be specified either by:

    - a rotation axis (arg1) and a rotation angle in degrees (arg2)

    - two axes (arg1 and arg2): axis arg1 is rotated into axis arg2 

    - three Euler angles in degrees arg1=(alpha, beta, gamma). alpha is a rotation along X (Ox->Ox, Oy->Oy1, Oz->Oz1), beta is a rotation along Y (Ox1->Ox2, Oy1->Oy1, Oz1->Oz2), gamma is a rotation along Z (Ox2->Ox3, Oy2->Oy3, Oz2->Oz2):


    :param a: mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param C: center of rotation
    :type C: 3-tuple of floats
    :param arg1: rotation axis or original axis or rotation angles (in degrees)
    :type arg1: 3-tuple of floats or 3-tuple of 3-tuple of floats
    :param arg2: angle of rotation (in degrees) or destination axis or None 
    :type arg2: float or 3-tuple of floats or None
    :param vectors: for each vector, list of the names of the vector components
    :type vectors: [list of list of strings]
    :return: mesh after rotation
    :rtype: identical to input

    *Example of use:*

    * `Rotate a mesh (array) <Examples/Transform/rotate.py>`_:

    .. literalinclude:: ../build/Examples/Transform/rotate.py

    * `Rotate a mesh (pyTree) <Examples/Transform/rotatePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/rotatePT.py

---------------------------------------

.. py:function:: Transform.translate(a, T)

    .. A2.O1.D1
    
    Translate a  mesh of vector T=(tx,ty,tz).

    Exists also as an in-place version (_translate) which modifies a and returns None.

    :param a: mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param T: translation vector
    :type T: 3-tuple of floats
    :return: mesh after translation
    :rtype: identical to input

    *Example of use:*

    * `Translate a mesh (array) <Examples/Transform/translate.py>`_:

    .. literalinclude:: ../build/Examples/Transform/translate.py

    * `Translate a mesh (pyTree) <Examples/Transform/translatePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/translatePT.py


---------------------------------------


Mesh transformation
--------------------------


.. py:function:: Transform.cart2Cyl(a, C, axis, thetaShift=0.)

    .. A1.O0.D1
    
    Convert a mesh in Cartesian coordinates into a mesh in cylindrical coordinates. 
    One of the Cartesian axes, defined by parameter AXIS, must be the revolution axis 
    of the cylindrical frame.
    Axis can be one of (0,0,1), (1,0,0) or (0,1,0).
    ThetaShift is an optional rotation angle that can be procided to position the grid
    in the [-PI,PI] quadrant. 

    Exists also as an in-place version (_cart2Cyl) which modifies a and returns None.

    :param a: mesh with coordinates defined in the Cartesian frame
    :type a: [array, list of arrays] or [zone, list of zone, base, pyTree]
    :param C: center of revolution
    :type C: 3-tuple of floats
    :param axis: revolution axis
    :type axis: 3-tuple of floats
    :param thetaShift: angle (in degrees) for the grid to be lie in [-PI,PI] 
    :type thetaShift: float
    :return: mesh with coordinates in the cylindrical frame
    :rtype: identical to input

    *Example of use:*

    * `Cart2Cyl a mesh (array) <Examples/Transform/cart2Cyl.py>`_:

    .. literalinclude:: ../build/Examples/Transform/cart2Cyl.py

    * `Cart2Cyl a mesh (pyTree) <Examples/Transform/cart2CylPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/cart2CylPT.py


---------------------------------------


.. py:function:: Transform.homothety(a, C, alpha)

    .. A2.O1.D1
    
    Apply an homothety of center C and a factor alpha to a mesh a.

    Exists also as an in-place version (_homothety) which modifies a and returns None.

    :param a: mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param C: center of homothety
    :type C: 3-tuple of floats
    :param alpha: homothety factor
    :type alpha: float
    :return: mesh after homothety
    :rtype: identical to input

    *Example of use:*

    * `Homothety a mesh (array) <Examples/Transform/homothety.py>`_:

    .. literalinclude:: ../build/Examples/Transform/homothety.py

    * `Homothety a mesh (pyTree) <Examples/Transform/homothetyPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/homothetyPT.py


---------------------------------------


.. py:function:: Transform.contract(a, C, dir1, dir2, alpha)

    .. A2.O1.D1

    Make a contraction of factor alpha of a mesh with respect to a plane defined by a point C and 
    vectors dir1 and dir2.

    Exists also as an in-place version (_contract) which modifies a and returns None.

    :param a: mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param C: point of the contraction plane 
    :type C: 3-tuple of floats
    :param dir1: first vector defining the plane
    :type dir1: 3-tuple of floats
    :param dir2: second vector defining the plane
    :type dir2: 3-tuple of floats
    :param alpha: contraction factor
    :type alpha: float
    :return: mesh after contraction
    :rtype: identical to input

    *Example of use:*

    * `Contract a mesh (array) <Examples/Transform/contract.py>`_:

    .. literalinclude:: ../build/Examples/Transform/contract.py

    * `Contract a mesh (pyTree) <Examples/Transform/contractPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/contractPT.py


---------------------------------------


.. py:function:: Transform.scale(a, factor=1., X=None)

    .. A2.O1.D1
  
    Scale a  mesh of factor factor. If factor is a list of floats, scale
    with given factor for each canonical axis. If invariant reference point X is not
    given, it is set to the barycenter of a.

    Exists also as an in-place version (_scale) which modifies a and returns None.

    :param a: mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param factor: scaling factor 
    :type factor: float or list of 3 floats
    :param X: reference point
    :type X: None or tuple of 3 floats
    :return: mesh after scaling
    :rtype: identical to input

    *Example of use:*

    * `Scale a mesh (array) <Examples/Transform/scale.py>`_:

    .. literalinclude:: ../build/Examples/Transform/scale.py

    * `Scale a mesh (pyTree) <Examples/Transform/scalePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/scalePT.py


---------------------------------------


.. py:function:: Transform.symetrize(a, P, vector1, vector2)

    .. A2.O0.D1
    
    Symmetrize a mesh with respect to a plane defined by point P and vectors vector1 and vector2.

    Exists also as an in-place version (_symetrize) which modifies a and returns None.

    :param a: mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param C: point of the symmetry plane 
    :type C: 3-tuple of floats
    :param vector1: first vector of the symetry plane
    :type vector1: 3-tuple of floats
    :param vector2: second vector of the symetry plane
    :type vector2: 3-tuple of floats
    :return: mesh after symmetrization
    :rtype: identical to input

    *Example of use:*

    * `Symmetrize a mesh (array) <Examples/Transform/symetrize.py>`_:

    .. literalinclude:: ../build/Examples/Transform/symetrize.py

    * `Symmetrize a mesh (pyTree) <Examples/Transform/symetrizePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/symetrizePT.py


---------------------------------------


.. py:function:: Transform.perturbate(a, radius, dim=3)

    .. A1.O0.D1
    
    Perturbate randomly a mesh a with given radius. Mesh points are modified aleatoirely
    in all directions, with a distance less or equal to radius.
    If dim=2, Z coordinates are fixed.
    If dim=1, only the X coordinates are modified.

    Exists also as an in-place version (_perturbate) which modifies a and returns None.

    :param a: mesh
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param radius: radius of perturbation
    :type radius: float 
    :param dim: to select if 1, 2 or the 3 coordinates are modified.
    :type dim: integer
    :return: mesh after perturbation
    :rtype: identical to input

    *Example of use:*

    * `Perturbate a mesh (array) <Examples/Transform/perturbate.py>`_:

    .. literalinclude:: ../build/Examples/Transform/perturbate.py

    * `Perturbate a mesh (pyTree) <Examples/Transform/perturbatePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/perturbatePT.py


---------------------------------------


.. py:function:: Transform.smooth(a, eps=0.5, niter=4, type=0, fixedConstraints=[], projConstraints=[], delta=1., point=(0,0,0), radius=-1.)

    .. A1.O0.D0

    Perform a Laplacian smoothing on a set of structured grids or an unstructured mesh ('QUAD', 'TRI') with a weight eps, and niter smoothing iterations. 
    Type=0 means isotropic Laplacian, type=1 means scaled Laplacian, type=2 means taubin smoothing.
    Constraints can be defined in order to avoid smoothing of some points (for instance the exterior faces of a): 

    Exists also as an in-place version (_smooth) which modifies a and returns None.

    :param a: input mesh
    :type a: array or zone
    :param eps: smoother power
    :type eps: float
    :param niter: number of smoothing iterations
    :type niter: integer
    :param type: type of smoothing algorithm
    :type type: integer
    :param fixedConstraints: set of fixed regions
    :type fixedConstraints: [list of arrays] or [list of zones]
    :param projConstraints: smoothed mesh projected on them
    :type projConstraints: [list of arrays] or [list of zones]
    :param delta: strength of constraints
    :type delta: float
    :param point: center of the region to be smoothed in case of local smoothing
    :type point: 3-tuple of float
    :param radius: if local smoothing, radius of the region to be smoothed
    :type radius: float 
    :return: mesh after smoothing
    :rtype: array or zone

    *Example of use:*

    * `Smooth a mesh (array) <Examples/Transform/smooth.py>`_:

    .. literalinclude:: ../build/Examples/Transform/smooth.py

    * `Smooth a mesh (pyTree) <Examples/Transform/smoothPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/smoothPT.py


---------------------------------------

.. py:function:: Transform.smoothField(a, eps=0.1, niter=1, type=0, varNames=[])

    Perform a Laplacian smoothing on given fields.

    Exists also as an in-place version (_smoothField) which modifies a and returns None.

    :param a: input zone with fields
    :type a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param eps: smoother power
    :type eps: float
    :param niter: number of smoothing iterations
    :type niter: integer
    :param type: type of smoothing algorithm
    :type type: integer 0 (isotropic) or 1 (scale)
    :param varNames: variable names that must be smoothed
    :type varNames: list of strings 

    *Example of use:*

    * `Smooth field (array) <Examples/Transform/smoothField.py>`_:

    .. literalinclude:: ../build/Examples/Transform/smoothField.py

    * `Smooth field (pyTree) <Examples/Transform/smoothFieldPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/smoothFieldPT.py

---------------------------------------

.. py:function:: Transform.dual(a, extraPoints=1)

    .. A1.O0.D0
    
    Return the dual of a mesh a. If extraPoints=1, external face centers are added.

    Exists also as an in-place version (_dual) which modifies a and returns None.

    :param a: mesh
    :type a: array or zone
    :param extraPoints: 0/1 external face centers are added
    :type extraPoints: integer
    :return: dual mesh
    :rtype: array or zone

    *Example of use:*

    * `Dual a mesh (array) <Examples/Transform/dual.py>`_:

    .. literalinclude:: ../build/Examples/Transform/dual.py

    * `Dual a mesh (pyTree) <Examples/Transform/dualPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/dualPT.py

---------------------------------------


.. py:function:: Transform.breakElements(a)

    .. A1.O0.D0
    
    Break a NGON mesh into a set of grids, each of them being a basic element grid (with a single connectivity).

    :param a: NGON mesh
    :type a: array or zone
    :return: list of grids of basic elements 
    :rtype: [list of arrays or list of zones]

    *Example of use:*

    * `Break a NGON mesh into basic elements (array) <Examples/Transform/breakElements.py>`_:

    .. literalinclude:: ../build/Examples/Transform/breakElements.py

    * `Break a NGON mesh into basic elements (pyTree) <Examples/Transform/breakElementsPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/breakElementsPT.py


---------------------------------------


Mesh splitting and merging
--------------------------

.. py:function:: Transform.subzone(a, minIndex, maxIndex=None, type=None)

    .. A1.O0.D0
    
    Extract a subzone.

    Extract a subzone of a structured mesh a, where min and max ranges  must be specified. Negative indices can be used (as in Python): -1 means max index::

         b = T.subzone(a, (imin,jmin,kmin), (imax,jmax,kmax))

    Extract a subzone of an unstructured mesh a, where the vertex list of the subzone must be specified (indices start at 1)::

         b = T.subzone(a, [1,2,...])

    Extract a subzone of an unstructured mesh providing the indices of elements (index starts at 0)::

         b = T.subzone(a, [0,1,...], type='elements')

    Extract a subzone of an unstructured array providing the indices of faces (for unstructured zones with basic elements: indFace=indElt*numberOfFaces+noFace, for NGON zones: use the natural face indexing, starting from 1)::

         b = T.subzone(a, [1,2,...], type='faces')


    :param a: input data
    :type  a: array or zone
    :param minIndex: (imin,jmin,kmin) for a structured grid, list of indices otherwise
    :type  minIndex:  3-tuple of integers
    :param maxIndex: (imax,jmax,kmax) for a structured grid, None otherwise
    :type  maxIndex:  3-tuple of integers
    :param type: type of subzone to perform (None, 'elements', 'faces')
    :type  type: None or string
    :return: subzoned mesh
    :rtype: identical to a

    *Example of use:*

    * `Subzone (array) <Examples/Transform/subzone.py>`_:

    .. literalinclude:: ../build/Examples/Transform/subzone.py

    * `Subzone (pyTree) <Examples/Transform/subzonePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/subzonePT.py


---------------------------------------

.. py:function:: Transform.join(a, b=None, tol=1.e-10)

    .. A1.O0.D0
    
    Join two zones in one (if possible) or join a list of zones in one zone (if possible).
    For the pyTree version, boundary conditions are maintained for structured grids only.

    :param a: input data
    :type  a: [array, list of arrays] or [zone, list of zones, base, pyTree]
    :param tol: tolerance for abutting grids 
    :type tol: float
    :return: unique joined zone
    :rtype: array or zone

    *Example of use:*

    * `Join (array) <Examples/Transform/join.py>`_:

    .. literalinclude:: ../build/Examples/Transform/join.py

    * `Join (pyTree) <Examples/Transform/joinPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/joinPT.py

---------------------------------------


.. py:function:: Transform.merge(a, sizeMax=1000000000, dir=0, tol=1.e-10, alphaRef=180., mergeBCs=False)

    .. A1.O0.D0
    
    Join a set of zones such that a minimum number of zones is obtained at the end. 
    Parameter sizeMax defines the maximum size of merged grids.
    dir is the constraint direction along which the merging is prefered. 
    Default value is 0 (no prefered direction), 1 for i, 2 for j, 3 for k. 
    alphaRef can be used for surface grids and avoids merging adjacent zones
    sharing an angle deviating of alphaRef to 180.


    For the pyTree version, boundary conditions are maintained for structured grids only.


    :param a:  list of grids
    :type  a:  [list of arrays] or [list of zones, base, pyTree]
    :param sizeMax: maximum size of merged grids
    :type sizeMax: integer
    :param dir: direction of merging (structured grids only): 0:ijk; 1:i; 2:j; 3:k
    :type dir: integer
    :param tol: tolerance for abutting grids 
    :type tol: float
    :param alphaRef: angle max of deviation for abutting grids, above which grids are not merged (for surface grids only)
    :type alphaRef: float
    :param mergeBCs: if True, merge BCs and perform connectMatch
    :type mergeBCs: boolean
    :return: list of merged grids
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Merge (array) <Examples/Transform/merge.py>`_:

    .. literalinclude:: ../build/Examples/Transform/merge.py

    * `Merge (pyTree) <Examples/Transform/mergePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/mergePT.py

---------------------------------------


.. py:function:: Transform.mergeCart(a, sizeMax=1000000000, tol=1.e-10)

    .. A1.O0.D0
    
    Merge a set of Cartesian grids. This function is similar to the function Transform.merge but is optimized for Cartesian grids.

    :param a:  list of Cartesian grids
    :type  a:  [list of arrays] or [list of zones, base, pyTree]
    :param sizeMax: maximum size of merged grids
    :type sizeMax: integer
    :param tol: tolerance for abutting grids 
    :type tol: float
    :return: list of merged Cartesian grids
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Merge Cartesian grids (array) <Examples/Transform/mergeCart.py>`_:

    .. literalinclude:: ../build/Examples/Transform/mergeCart.py

    * `Merge Cartesian grids (pyTree) <Examples/Transform/mergeCartPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/mergeCartPT.py

---------------------------------------



.. py:function:: Transform.splitNParts(a, N, multigrid=0, dirs=[1,2,3], recoverBC=True, topTree=None)

    .. A1.O0.D0
    
    Split a set of M grids into N parts of same size roughly, provided M < N. 

    Argument multigrid enables to ensure the multigrid level by the splitting, provided the input grids are of that multigrid level. It can also be useful to split at nearmatch interfaces (multigrid=1 for 1:2 interfaces and multigrid 2 for 1:4 interfaces).

    For the pyTree version, boundary conditions and matching connectivity are split.

    Exists also as in place version (_splitNParts) that modifies a and returns None. In this case, a must be a pyTree.

    :param a:  list of grids
    :type  a:  [list of arrays] or [list of zones, base, pyTree]
    :param N: number of grids after splitting
    :type N: integer
    :param multigrid: for structured grids only. 0: no constraints; 1: grids are 2n+1 per direction ; 2: grids are 4n+1 per direction
    :type multigrid: integer
    :param dirs: directions where splitting is allowed (for structured grids only)
    :type dirs: list of integers (possible values:1,2,3 or a combination of them)
    :param recoverBC: BCs are recovered after split (True) or not (False)
    :type recoverBC: Boolean (True or False)
    :param topTree: if a is not the top tree, provides full tree for match updates
    :type topTree: CGNS Tree
    :return: list of splitted grids
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split a mesh in N parts (array) <Examples/Transform/splitNParts.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitNParts.py

    * `Split a mesh in N parts  (pyTree) <Examples/Transform/splitNPartsPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitNPartsPT.py

---------------------------------------

.. py:function:: Transform.splitSize(a, N, multigrid=0, dirs=[1,2,3], type=0, R=None, minPtsPerDir=5, topTree=None)

    .. A1.O0.D0
    
    - Split structured blocks if their number of points is greater than N. 

    - splitSize can also be used to split blocks in order to fit as better as possible on a number of R processors. 

    Argument multigrid enables to ensure the multigrid level by the splitting, provided the input grids are of that multigrid level.

    For the pyTree version, boundary conditions and matching connectivity are split.

    Exists also as in place version (_splitSize) that modifies a and returns None. In this case, a must be a pyTree.

    :param a:  list of grids
    :type  a:  [list of arrays] or [list of zones, base, pyTree]
    :param N: number of grids after splitting
    :type N: integer
    :param multigrid: for structured grids only. 0: no constraints; 1: grids are 2n+1 per direction ; 2: grids are 4n+1 per direction
    :type multigrid: integer
    :param dirs: directions where splitting is allowed (for structured grids only)
    :type dirs: list of integers (possible values:1,2,3 or a combination of them)
    :param type: only for split by size (not resources): 0: centered splitting; 1: upwind splitting when better
    :type type: integer
    :param R: number of resources (processors)
    :type R: integer
    :param minPtsPerDir: minimum number of points per direction
    :type minPtsPerDir: integer
    :param topTree: if a is not the top tree, provides full tree for match updates
    :type topTree: CGNS Tree
    :return: list of splitted grids
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split a mesh  by size (array) <Examples/Transform/splitSize.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitSize.py

    * `Split a mesh by size  (pyTree) <Examples/Transform/splitSizePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitSizePT.py

---------------------------------------


.. py:function:: Transform.splitCurvatureAngle(a, sensibility)

    .. A1.O0.D0
    
    Split a curve defined by a 1D structured grid with respect to the curvature angle. 
    If angle is lower than 180-sensibility in degrees or greater than 180+sensibility degrees, curve is split.

    :param a:  list of grids
    :type  a:  array or zone
    :param sensibility: sensibility angle (in degrees) to allow splitting
    :type sensibility: float
    :return: list of split curves
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split a curve wrt curvature angle (array) <Examples/Transform/splitCurvatureAngle.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitCurvatureAngle.py

    * `Split a curve wrt curvature angle (pyTree) <Examples/Transform/splitCurvatureAnglePT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitCurvatureAnglePT.py


---------------------------------------


.. py:function:: Transform.splitCurvatureRadius(a, Rs=100.)

    .. A1.O0.D0
    
    Split a curve defined by a 1D structured grid with respect to the curvature radius, using B-Spline approximations.
    The curve can be closed or not. 

    :param a: input mesh
    :type  a:  array or zone
    :param Rs: threshold curvature radius below which the initial curve is split
    :type Rs: float
    :return: list of split curves
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split a curve wrt curvature radius (array) <Examples/Transform/splitCurvatureRadius.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitCurvatureRadius.py

    * `Split a curve wrt curvature radius (pyTree) <Examples/Transform/splitCurvatureRadiusPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitCurvatureRadiusPT.py

---------------------------------------


.. py:function:: Transform.splitConnexity(a)

    .. A1.O0.D0 SLOW

    Split an unstructured mesh into connex parts.

    :param a: input unstructured mesh
    :type  a: array or zone
    :return: list of connex parts
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split an unstructured mesh into connex parts (array) <Examples/Transform/splitConnexity.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitConnexity.py

    * `Split an unstructured mesh into connex parts (pyTree) <Examples/Transform/splitConnexityPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitConnexityPT.py

---------------------------------------


.. py:function:: Transform.splitMultiplePts(a, dim=3)

    .. A1.O0.D0

    Split a structured mesh at external nodes connected to an even number of points, meaning that the geometrical point 
    connects an odd number of blocks.

    :param a: input set of structured grids
    :type  a: [list of arrays] or [list of zones]
    :return: set of structured grids after splitting
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split a mesh at odd connections (array) <Examples/Transform/splitMultiplePts.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitMultiplePts.py

    * `Split a mesh at odd connections (pyTree) <Examples/Transform/splitMultiplePtsPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitMultiplePtsPT.py

---------------------------------------


.. py:function:: Transform.PyTree.splitFullMatch(a)

    Split a structured mesh such that all match boundaries are full block faces.
    
    Exists also as in place version (_splitFullMatch) that modifies a and returns None. In this case, a must be a pyTree.

    :param a: input set of structured grids
    :type  a: [pyTree, base or list of zones]
    :return: split zones
    :rtype: identical to input

    *Example of use:*

    * `Split a mesh for matching on full faces (pyTree) <Examples/Transform/splitFullMatchPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitFullMatchPT.py

---------------------------------------


.. py:function:: Transform.splitSharpEdges(a, alphaRef=30.)

    .. A1.O0.D0

    Split a 1D or 2D mesh at edges sharper than alphaRef.
    If the input grid is structured, then it returns an unstructured grid (BAR or QUAD).

    :param a: input mesh
    :type  a: [array, list of arrays] or [zone, list of zones]
    :param alphaRef: angle (in degrees) below which the mesh must be split
    :type alphaRef: float
    :return: set of unstructured grids (with no sharp edges)
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split a surface at sharp edges (array) <Examples/Transform/splitSharpEdges.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitSharpEdges.py

    * `Split a  surface at sharp edges (pyTree) <Examples/Transform/splitSharpEdgesPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitSharpEdgesPT.py

---------------------------------------

.. py:function:: Transform.splitTBranches(a, tol=1.e-13)

    .. A1.O0.D0
    
    Split a curve defined by a 'BAR' if it has T-branches.

    :param a: input mesh
    :type  a: [array, list of arrays] or [zone, list of zones]
    :param tol: matching tolerance between points that define two branches
    :type tol: float 
    :return: set of BAR grids (with no T-branches) 
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split T-branches (array) <Examples/Transform/splitTBranches.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitTBranches.py

    * `Split T-branches (pyTree) <Examples/Transform/splitTBranchesPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitTBranchesPT.py

---------------------------------------

.. py:function:: Transform.splitManifold(a)

    .. A1.O0.D0

    Split a unstructured mesh (TRI or BAR only) into manifold pieces.

    :param a: input mesh (TRI or BAR)
    :type  a: [array, list of arrays] or [zone, list of zones, base, pyTree] 
    :return: set of TRI or BAR grids
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split into manifold parts (array) <Examples/Transform/splitManifold.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitManifold.py

    * `Split into manifold parts (pyTree) <Examples/Transform/splitManifoldPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitManifoldPT.py

---------------------------------------

.. py:function:: Transform.splitBAR(a, N, N2=-1)

    .. A1.O0.D1 

    Split a curve defined by a BAR at index N.
    If N2 is provided, split also at index N2.

    :param a: input mesh (BAR)
    :type  a: array or zone 
    :param N: index of split in a
    :type N: integer
    :param N2: optional second split index
    :type N2: integer
    :return: two BARS 
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split a BAR at given index (array) <Examples/Transform/splitBAR.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitBAR.py

    * `Split a BAR at given index (pyTree) <Examples/Transform/splitBARPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitBARPT.py

---------------------------------------


.. py:function:: Transform.splitTRI(a, idxList)

    .. A1.O0.D0

    Split a triangular mesh into several triangular grids delineated by the polyline of indices idxList in the original TRI mesh.

    :param a: input mesh (TRI)
    :type  a: array or zone 
    :param idxList: indices of split in a defining a polyline
    :type idxList: list of integers
    :return: a set of TRI grids
    :rtype: [list of arrays] or [list of zones]

    *Example of use:*

    * `Split a TRI mesh with indices (array) <Examples/Transform/splitTRI.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitTRI.py

    * `Split a TRI mesh with indices (pyTree) <Examples/Transform/splitTRIPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/splitTRIPT.py

---------------------------------------


Mesh deformation
--------------------------



.. py:function:: Transform.deform(a, vector=['dx','dy','dz'])

    .. A1.O0.D0
    
    Deform a surface by moving each point of a vector. The vector field must be defined in a, with the same location.

    Exists also as an in-place version (_deform) which modifies a and returns None.

    :param a: input surface mesh, containing the vector fields
    :type  a: [array, list of arrays] or [zone, list of zones] 
    :param vector: vector component names defined in a
    :type vector: list of 3 strings
    :return: deformed surface mesh
    :rtype: identical to input

    *Example of use:*

    * `Deform a surface (array) <Examples/Transform/deform.py>`_:

    .. literalinclude:: ../build/Examples/Transform/deform.py

    * `Deform a surface (pyTree) <Examples/Transform/deformPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/deformPT.py

---------------------------------------

.. py:function:: Transform.deformNormals(a, alpha, niter=1)

    .. A1.O0.D0
    
    Deform a surface mesh a by moving each point of the surface by a scalar field alpha times the surface normals in niter steps.

    Exists also as an in-place version (_deformNormals) which modifies a and returns None.

    :param a: input surface mesh, containing the vector fields
    :type  a: [array, list of arrays] or [zone, list of zones] 
    :param alpha: factor of growth wrt to normals
    :type alpha: float
    :param niter: number of steps (raise it to increase the smoothing of the resulting surface)
    :type niter: integer
    :return: deformed surface mesh
    :rtype: identical to input

    *Example of use:*

    * `Deform a surface following normals (array) <Examples/Transform/deformNormals.py>`_:

    .. literalinclude:: ../build/Examples/Transform/deformNormals.py

    * `Deform a surface following normals (pyTree) <Examples/Transform/deformNormalsPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/deformNormalsPT.py

---------------------------------------

.. py:function:: Transform.deformPoint(a, xyz, dxdydz, depth, width)

    .. A1.O0.D0
    
    Deform a surface mesh a by moving point P of vector V. Argument 'depth' controls the depth of deformation. 
    Argument 'width' controls the width of deformation.

    Exists also as an in-place version (_deformPoint) which modifies a and returns None.

    :param a: input surface mesh, containing the vector fields
    :type  a: [array, list of arrays] or [zone, list of zones] 
    :param P: point that is moved
    :type P: 3-tuple of floats
    :param V: vector of deformation
    :type V: 3-tuple of floats
    :param depth: to control the depth of deformation
    :type depth: float
    :param width: to control the width of deformation
    :type width: float
    :return: deformed surface mesh
    :rtype: identical to input

    *Example of use:*

    * `Deform a surface at a point P (array) <Examples/Transform/deformPoint.py>`_:

    .. literalinclude:: ../build/Examples/Transform/deformPoint.py

    * `Deform a surface at a point P (pyTree) <Examples/Transform/deformPointPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/deformPointPT.py

---------------------------------------

.. py:function:: Transform.deformMesh(a, surfDelta, beta=4., type='nearest')

    .. A1.O0.D0
    
    Deform a mesh defined by a given surface or a set of surfaces for which a deformation is defined at nodes as a vector field 'dx,dy,dz'.
    The surface surfDelta does not necessary match with a border of the meshes.
    Beta enables to extend the deformation region as multiplication factor of local deformation.

    Exists also as an in-place version (_deformMesh) which modifies a and returns None.

    :param a: input surface mesh, containing the vector fields
    :type  a: [array, list of arrays] or [zone, list of zones] 
    :param surfDelta: surface on which the deformation is defined
    :type surfDelta: [array,list of arrays] or [zone,list of zones] 
    :return: deformed mesh
    :rtype: identical to input

    *Example of use:*

    * `Deform a mesh (array) <Examples/Transform/deformMesh.py>`_:

    .. literalinclude:: ../build/Examples/Transform/deformMesh.py

    * `Deform a mesh (pyTree) <Examples/Transform/deformMeshPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/deformMeshPT.py

---------------------------------------


Mesh projections
--------------------------

.. py:function:: Transform.projectAllDirs(a, s, vect=['nx','ny','nz'], oriented=0)

    .. A1.O0.D0
    
    Project a surface mesh a onto a set of surfaces s according to a vector defined for each point of the mesh a.
    If oriented=0, both directions are used for projection, else the vector direction is used.

    Exists also as an in-place version (_projectAllDirs) which modifies a and returns None.

    :param a: input surface mesh, containing the vector fields
    :type  a: [array, list of arrays] or [zone, list of zones] 
    :param s: projection surface
    :type s: [array,list of arrays] or [zone,list of zones] 
    :param vect: vector component names 
    :type vect: list of 3 strings
    :param oriented: 0 for projection in the vector direction and also in its opposite direction
    :type oriented: integer (0 or 1)
    :return: projected mesh
    :rtype: identical to input

    *Example of use:*

    * `Project a mesh (array) <Examples/Transform/projectAllDirs.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectAllDirs.py

    * `Project a mesh (pyTree) <Examples/Transform/projectAllDirsPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectAllDirsPT.py

---------------------------------------


.. py:function:: Transform.projectDir(a, s, dir, smooth=0, oriented=0)

    .. A1.O0.D1
    
    Project a surface mesh a onto a set of surfaces s following a constant direction dir.
    If oriented=0, both directions are used for projection, else the vector direction is used.
    If smooth=1, points that cannot be projected are smoothed (available only for structured grids).

    Exists also as an in-place version (_projectDir) which modifies a and returns None.

    :param a: input surface mesh
    :type  a: [array list of arrays] or [zone, list of zones] 
    :param s: projection surface
    :type s: [array, list of arrays] or [zone, list of zones] 
    :param dir: constant vector that directs the projection 
    :type dir: 3-tuple of floats
    :param smooth: smoothing of unprojected points
    :type smooth: integer (0 or 1)
    :param oriented: 0 for projection in the vector direction and also in its opposite direction
    :type oriented: integer (0 or 1)
    :return: projected mesh
    :rtype: identical to input

    *Example of use:*

    * `Project a mesh following a direction (array) <Examples/Transform/projectDir.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectDir.py

    * `Project a mesh following a direction (pyTree) <Examples/Transform/projectDirPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectDirPT.py

---------------------------------------


.. py:function:: Transform.projectOrtho(a, s)

    .. A1.O0.D1
    
    Project a surface mesh a orthogonally onto a set of surfaces s.

    Exists also as an in-place version (_projectOrtho) which modifies a and returns None.

    :param a: input surface mesh, containing the vector fields
    :type  a: [array, list of arrays] or [zone, list of zones] 
    :param s: projection surface
    :type s: [array, list of arrays] or [zone, list of zones] 
    :return: projected mesh
    :rtype: identical to input

    *Example of use:*

    * `Project a mesh orthogonally (array) <Examples/Transform/projectOrtho.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectOrtho.py

    * `Project a mesh orthogonally (pyTree) <Examples/Transform/projectOrthoPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectOrthoPT.py

---------------------------------------


.. py:function:: Transform.projectOrthoSmooth(a, s, niter=1)

    .. A1.O0.D1
    
    Project a surface mesh a following smoothed normals onto a set of surfaces s.

    Exists also as an in-place version (_projectOrthoSmooth) which modifies a and returns None.

    :param a: input surface mesh
    :type  a: [array, list of arrays] or [zone, list of zones] 
    :param s: projection surface 
    :type s: [array, list of arrays] or [zone, list of zones] 
    :param niter: number of smoothing iterations
    :type niter: integer
    :return: projected mesh
    :rtype: identical to input

    *Example of use:*

    * `Project a mesh following smoothed normals (array) <Examples/Transform/projectOrthoSmooth.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectOrthoSmooth.py

    * `Project a mesh following smoothed normals (pyTree) <Examples/Transform/projectOrthoSmoothPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectOrthoSmoothPT.py

---------------------------------------


.. py:function:: Transform.projectRay(a, s, P)

    .. A1.O0.D1
    
    Project a surface mesh a onto a set of surfaces s following rays starting from a point P.

    Exists also as an in-place version (_projectRay) which modifies a and returns None.

    :param a: input surface mesh
    :type  a: [array, list of arrays] or [zone, list of zones] 
    :param s: projection surface 
    :type s: [array, list of arrays] or [zone, list of zones] 
    :param P: starting point of rays
    :type P: 3-tuple of floats
    :return: projected mesh
    :rtype: identical to input

    *Example of use:*

    * `Project a mesh following rays (array) <Examples/Transform/projectRay.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectRay.py

    * `Project a mesh following rays (pyTree) <Examples/Transform/projectRayPT.py>`_:

    .. literalinclude:: ../build/Examples/Transform/projectRayPT.py

---------------------------------------


.. toctree::
   :maxdepth: 2   


Indices and tables
###################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

