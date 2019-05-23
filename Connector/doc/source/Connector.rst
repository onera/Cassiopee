.. Connector documentation master file


Connector: Grid connectivity module
=============================================


Preamble
######## 

Connector module is used to compute connectivity between meshes.
It manipulates arrays (as defined in Converter documentation)
or CGNS/Python trees (pyTrees) as data structures.

This module is part of Cassiopee, a free open-source
pre- and post-processor for CFD simulations.


To use the Connector module with the array interface::

    import Connector as X 

With the pyTree interface::

    import Connector.PyTree as X 

.. py:module:: Connector


List of functions
##################


**-- Multiblock connectivity**

.. autosummary::

   Connector.connectMatch
   Connector.PyTree.connectMatchPeriodic   
   Connector.PyTree.connectNearMatch
   Connector.PyTree.setDegeneratedBC

**-- Overset grid connectivity**

.. autosummary::

   Connector.blankCells
   Connector.blankCellsTetra
   Connector.blankCellsTri
   Connector.blankIntersectingCells
   Connector.setHoleInterpolatedPoints
   Connector.optimizeOverlap
   Connector.maximizeBlankedCells

   Connector.PyTree.cellN2OversetHoles
   Connector.PyTree.setInterpData
   Connector.PyTree.getOversetInfo
   Connector.PyTree.extractChimeraInfo


**-- Overset grid connectivity for elsA solver**

.. autosummary::

    Connector.PyTree.setInterpolations
    Connector.PyTree.chimeraInfo

**-- Immersed boundary connectivity**

.. autosummary::

    


Contents
#########

Multiblock connectivity
-------------------------

.. py:function:: Connector.connectMatch

    Detect and set all matching windows, even partially.

    *Using the array interface:*

        ::

         res = X.connectMatch(a1, a2, sameZone=0, tol=1.e-6, dim=3)

        Detect and set all matching windows between two structured arrays a1 and a2.
        Return the subrange of indices of abutting windows and an index transformation from a1 to a2.
        If the CFD problem is 2D, then dim must be set to 2.
        Parameter sameZone must be set to 1 if a1 and a2 define the same zone.

    :param a1,a2:  Input data
    :type  a1,a2:  arrays
    

    *Using the PyTree interface:*

        ::

         t = X.connectMatch(t, tol=1.e-6, dim=3)
        
        Detect and set all matching windows in a zone node, a list of zone nodes or a complete pyTree.
        Set automatically the Transform node corresponding to the transformation from matching block 1 to block 2.
        If the CFD problem is 2D, then dim must be set to 2.

    :param t: input data
    :type  t: pyTree, base, zone, list of zones
    :rtype:  identical to input

    *Example of use:*

    * `Detect matching boundaries of a mesh (array) <Examples/Connector/connectMatch.py>`_:

    .. literalinclude:: ../build/Examples/Connector/connectMatch.py

    * `Add 1-to-1 abutting connectivity in a pyTree (pyTree) <Examples/Connector/connectMatchPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/connectMatchPT.py

---------------------------------------------------------------------------

.. py:function:: Connector.PyTree.connectMatchPeriodic(t, rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.], translation=[0.,0.,0.], tol=1.e-6, dim=3, unitAngle=None)

    Detect and set all periodic matching borders, even partially, in a zone node, a list of zone nodes, a base, or a full pyTree.
    Periodicity can be defined either by rotation or translation or by a composition of rotation and translation.


    :param t: input data
    :type  t: pyTree, base, zone, list of zones
    :rtype:  identical to input


    Set automatically the Transform node corresponding to the transformation from matching block 1 to block 2, and the 'GridConnectivityProperty/Periodic' for periodic matching BCs.

    If the CFD problem is 2D, then dim must be set to 2.

    For periodicity by rotation, the rotation angle units can be specified by argument unitAngle, which can be 'Degree','Radian',None.

    If unitAngle=None or 'Degree': parameter rotationAngle is assumed to be defined in degrees.

    If unitAngle='Radian': parameter rotationAngle is assumed in radians.


    .. note:: 

      - if the mesh is periodic in rotation and in translation separately (i.e. connecting with some blocks in rotation, and some other blocks in translation), the function must be applied twice.

      - Since *Cassiopee2.6*: 'RotationAngle' node in 'Periodic' node is always defined in Radians. A DimensionalUnits child node is also defined.
    

    *Example of use:*

    * `Add periodic  1-to-1 abutting grid connectivity in a pyTree (pyTree) <Examples/Connector/connectMatchPeriodicPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/connectMatchPeriodicPT.py

---------------------------------------------------------------------------

.. py:function:: Connector.PyTree.connectNearMatch(t, ratio=2, tol=1.e-6, dim=3)

    Detect and set all near-matching windows, even partially in a zone node, a list of zone nodes or a complete pyTree.
    A 'UserDefinedData' node is set, with the PointRangeDonor, the Transform and NMRatio nodes providing information for the opposite zone.
    .. warning:: connectMatch must be applied first if matching windows exist.

    Parameter ratio defines the ratio between near-matching windows and can be an integer (e.g. 2) or a list of 3 integers (e.g. [1,2,1]), specifying
    the nearmatching direction to test (less CPU-consuming).
    If the CFD problem is 2D, then dim must be set to 2.
    

    :param t: input data
    :type  t: pyTree, base, zone, list of zones
    :rtype:  identical to input

    *Example of use:*

    * `Add n-to-m abutting grid connectivity in a pyTree (pyTree) <Examples/Connector/connectNearMatchPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/connectNearMatchPT.py

---------------------------------------------------------------------------

.. py:function:: Connector.PyTree.setDegeneratedBC(t, dim=3, tol=1.e-10)

    Detect all degenerated lines in 3D zones and define a BC as a 'BCDegenerateLine' BC type. For 2D zones, 'BCDegeneratePoint'
    type is defined.
    If the problem is 2D according to (i,j), then parameter 'dim' must be set to 2.
    Parameter 'tol' defines a distance below which a window is assumed degenerated.
    
    *Example of use:*

    * `Add degenerated line as BCs in a pyTree (pyTree) <Examples/Connector/setDegeneratedBCPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setDegeneratedBCPT.py


---------------------------------------------------------------------------

Overset connectivity
-------------------------

.. py:function:: Connector.blankCells()

    Blank cells using X-Ray method.
    
    *Using the array interface:*

        ::

          cellns = X.blankCells(coords, cellns, body, blankingType=2, delta=1.e-10, dim=3, masknot=0, tol=1.e-8)

        Blank the cells of a list of grids defined by coords (located at nodes).
        The X-Ray mask is defined by bodies, which is a list of arrays.
        Cellnaturefield defined in cellns is modified (0: blanked points, 1: otherwise).
        Some parameters can be specified: blankingType, delta, masknot, tol. Their meanings are described in the table below:

        +---------------------------+----------------------------------------------------------------------------------+
        |  Parameter value          |  Meaning                                                                         | 
        +===========================+==================================================================================+
        |    blankingType=0         |  blank nodes inside bodies (node_in).                                            | 
        +---------------------------+----------------------------------------------------------------------------------+
        |    blankingType=2         |  blank cell centers inside bodies (center_in).                                   | 
        +---------------------------+----------------------------------------------------------------------------------+
        |    blankingType=1         |  blank cell centers intersecting with body (cell_intersect).                     | 
        +---------------------------+----------------------------------------------------------------------------------+
        |    blankingType=-2        |  blank cell centers using an optimized cell intersection (cell_intersect_opt)    |
        |                           |  and interpolation depth=2 (blanking region may be reduced where blanking point  | 
        |                           |  can be interpolated).                                                           | 
        +---------------------------+----------------------------------------------------------------------------------+
        |    blankingType=-1        |  blank cell centers using an optimized cell intersection (cell_intersect_opt)    | 
        |                           |  and interpolation depth=1.                                                      | 
        +---------------------------+----------------------------------------------------------------------------------+
        |    delta=0.               |  cells are blanked in the body                                                   | 
        +---------------------------+----------------------------------------------------------------------------------+
        | delta greater than 0.     |  the maximum distance to body, in which cells are blanked                        | 
        +---------------------------+----------------------------------------------------------------------------------+
        | masknot=0                 |  Classical blanking applied                                                      | 
        +---------------------------+----------------------------------------------------------------------------------+
        | masknot=1                 |  Inverted blanking applied: cells out of the body are blanked                    | 
        +---------------------------+----------------------------------------------------------------------------------+
        | dim=3                     |  body described by a surface and blanks 3D cells.                                | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  dim=2                    |  body blanks 2D or 1D zones.                                                     | 
        +---------------------------+----------------------------------------------------------------------------------+
        | tol=1.e-8 (default)       |  tolerance for the multiple definition of the body.                              | 
        +---------------------------+----------------------------------------------------------------------------------+


        .. note:: in case of blankingType=0, location of cellns and coords must be identical.
    
    *Using the pyTree interface:*
        ::

         B = X.blankCells(t, bodies, BM, depth=2, blankingType='cell_intersect', delta=1.e-10, dim=3, tol=1.e-8, XRaydim1=1000, XRaydim2=1000)


        **blankCells** function sets the cellN to 0 to blanked nodes or cell centers of both structured and unstructured grids.

        The location of the cellN field depends on the *blankingType* parameter: if 'node_in' is used, nodes are blanked, else centers are blanked.

        The mesh to be blanked is defined by a pyTree t, where each basis defines a Chimera component. The list of bodies blanking the grids is defined in bodies.

        Each element of the list bodies is a set of CGNS/Python zones defining a closed and watertight surface.

        The blanking matrix BM is a numpy array of size nbases x nbodies.

        BM(i,j)=1 means that ith basis is blanked by jth body.

        BM(i,j)=0 means no blanking, and BM(i,j)=-1 means that inverted hole-cutting is performed.

        blankingType can be 'cell_intersect', 'cell_intersect_opt', 'center_in' or 'node_in'. Parameter depth is only meaningfull for 'cell_intersect_opt'.

        XRaydim1 and XRaydim2 are the dimensions of the X-Ray hole-cutting in the x and y directions in 3D.

        If the variable 'cellN' does not exist in the input pyTree, it is initialized to 1, located at 'nodes' if 'node_in' is set, and at centers in other cases.

        .. warning:: 'cell_intersect_opt' can be CPU time-consuming when delta>0.

    *Example of use:*

    * `Blank cells (array) <Examples/Connector/blankCells.py>`_:

    .. literalinclude:: ../build/Examples/Connector/blankCells.py

    * `Blank cells (pyTree) <Examples/Connector/blankCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/blankCellsPT.py


---------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.blankCellsTetra()

    *Using the array interface:*

        ::

          cellns = X.blankCellsTetra(coords, cellns, body, blankingType=2, tol=1.e-12)
    

        Blanks the input grids nodes or cells (depending on the *blankingType* value) that fall inside a volume body mask.
        The blanking is achieved by setting the Cellnaturefield to *cellnval* (0 by default) in *cellns*.
        
        The input grids are defined by coords located at nodes as a list of arrays. The body mask is defined by sets of tetrahedra in any orientation, as a list of arrays.

        If the *blankingMode* is set to 1 (overwrite mode), Cellnaturefield is reset to 1 for any node/cell outside the body mask.
        Hence the value of 1 is forbidden for cellnval upon entry (it will be replaced by 0).

        The parameters meanings and values are described in the table below:

        +---------------------------+----------------------------------------------------------------------------------+
        |  Parameter value          |  Meaning                                                                         | 
        +===========================+==================================================================================+
        |  blankingType=0           |  blanks the nodes falling inside the body masks (node_in).                       | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  blankingType=2           |  blanks the cells having their center falling inside the body masks (center_in). | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  blankingType=1           |  blanks the cells that intersect or fall inside the body masks (cell_intersect). | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  tol=1.e-12               |  tolerance for detecting intersections (NOT USED CURRENTLY).                     | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  cellnval=0 (default)     |  value used for flagging as blanked.                                             | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  blankingMode=0 (default) |  Appending mode: cellns is only modified for nodes/cells falling inside the body |
        |                           |  mask by setting the value in cellns to cellnval.                                | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  blankingMode=1           |  Overwriting mode: cellns is modified for both nodes/cells falling inside        |
        |                           |  (set to cellnval) and outside (set to 1) the body mask.                         | 
        +---------------------------+----------------------------------------------------------------------------------+
    
        .. warning:: in case of blankingType=0, location of cellns and coords must be identical.
    

    *Using the pyTree interface:*

        ::

          B = X.blankCellsTetra(t, bodies, BM, blankingType='node_in', tol=1.e-12, cellnval=0, overwrite=0)
    
        Blanks the input grids nodes or cells (depending on the *blankingType* value) that fall inside a volume body mask.
        
        The blanking is achieved by setting the Cellnaturefield to *cellnval* (0 by default) in *cellns*.
        
        The mesh to be blanked is defined by a pyTree t, where each basis defines a Chimera component. The list of bodies blanking the grids is defined in bodies.

        Each element of the list bodies is a set of CGNS/Python zones defining a tetrahedra mesh.

        The blanking matrix BM is a numpy array of size nbases x nbodies.

        BM(i,j)=1 means that ith basis is blanked by jth body.

        BM(i,j)=0 means no blanking, and BM(i,j)=-1 means that inverted hole-cutting is performed.

        blankingType can be 'cell_intersect', 'center_in' or 'node_in'.

        If the variable 'cellN' does not exist in the input pyTree, it is initialized to 1, located at 'nodes' if 'node_in' is set, and at centers in other cases.

        If the *overwrite* is set to 1 (overwrite mode), Cellnaturefield is reset to 1 for any node/cell outside the body mask.

        Hence the value of 1 is forbidden for cellnval upon entry (it will be replaced by 0). 

    *Example of use:*

    * `Blank cells with a tetra mesh (array) <Examples/Connector/blankCellsTetra.py>`_:

    .. literalinclude:: ../build/Examples/Connector/blankCellsTetra.py

    * `Blank cells with a tetra mesh (pyTree) <Examples/Connector/blankCellsTetraPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/blankCellsTetraPT.py


-----------------------------------------------------------------------------------------------------------

.. py:function:: Connector.blankCellsTri()

    *Using the array interface:*

        ::

         cellns = X.blankCellsTri(coords, cellns, body, blankingType=2, tol=1.e-12, cellnval=0, blankingMode=0)
    
        Blanks the input grids nodes or cells (depending on the *blankingType* value) that fall inside a surfacic body mask.
        
        The blanking is achieved by setting the Cellnaturefield to *cellnval* (0 by default) in *cellns*.
        
        The input grids are defined by coords located at nodes as a list of arrays. The body mask is defined by triangular surface meshes in any orientation, as a list of arrays.

        If the *blankingMode* is set to 1 (overwrite mode), Cellnaturefield is reset to 1 for any node/cell outside the body mask.
        Hence the value of 1 is forbidden for cellnval upon entry (it will be replaced by 0).

        The parameters meanings and values are described in the table below:

        +---------------------------+----------------------------------------------------------------------------------+
        |  Parameter value          |  Meaning                                                                         | 
        +===========================+==================================================================================+
        |  blankingType=0           |  blanks the nodes falling inside the body masks (node_in).                       | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  blankingType=2           |  blanks the cells having their center falling inside the body masks (center_in). | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  blankingType=1           |  blanks the cells that intersect or fall inside the body masks (cell_intersect). | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  tol=1.e-12 (default)     |  tolerance for detecting intersections (NOT USED CURRENTLY).                     | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  cellnval=0 (default)     |  value used for flagging as blanked.                                             | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  blankingMode=0 (default) |  Appending mode: cellns is only modified for nodes/cells falling inside the body |
        |                           |  mask by setting the value in cellns to cellnval.                                | 
        +---------------------------+----------------------------------------------------------------------------------+
        |  blankingMode=1           |  Overwriting mode: cellns is modified for both nodes/cells falling inside        |
        |                           |  (set to cellnval) and outside (set to 1) the body mask.                         | 
        +---------------------------+----------------------------------------------------------------------------------+
    
        .. warning:: in case of blankingType=0, location of cellns and coords must be identical.
    

    *Using the pyTree interface:*

        ::

         B = X.blankCellsTri(t, bodies, BM, blankingType='node_in', tol=1.e-12, cellnval=0, overwrite=0)
    
        Blanks the input grids nodes or cells (depending on the *blankingType* value) that fall inside a volume body mask.

        The blanking is achieved by setting the Cellnaturefield to *cellnval* (0 by default) in *cellns*.
        
        The mesh to be blanked is defined by a pyTree t, where each basis defines a Chimera component. The list of bodies blanking the grids is defined in bodies.

        Each element of the list bodies is a set of CGNS/Python zones defining a triangular watertight closed surface. 
        
        The blanking matrix BM is a numpy array of size nbases x nbodies.

        BM(i,j)=1 means that ith basis is blanked by jth body.

        BM(i,j)=0 means no blanking, and BM(i,j)=-1 means that inverted hole-cutting is performed.

        blankingType can be 'cell_intersect', 'center_in' or 'node_in'.

        If the variable 'cellN' does not exist in the input pyTree, it is initialized to 1, located at 'nodes' if 'node_in' is set, and at centers in other cases.

        If the *overwrite* is set to 1 (overwrite mode), Cellnaturefield is reset to 1 for any node/cell outside the body mask.

        Hence the value of 1 is forbidden for cellnval upon entry (it will be replaced by 0).
    

    *Example of use:*

    * `Blank cells with a triangular surface mask (array) <Examples/Connector/blankCellsTri.py>`_:

    .. literalinclude:: ../build/Examples/Connector/blankCellsTri.py

    * `Blank cells with a triangular surface mask (pyTree) <Examples/Connector/blankCellsTriPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/blankCellsTriPT.py


------------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.setHoleInterpolatedPoints()

    *Using the array interface:*
    
        ::

         a = X.setHoleInterpolatedPoints(a, depth=2, dir=0, loc='centers', cellNName='cellN')
    
        Compute the fringe of interpolated points around a set of blanked points in a mesh a.
        Parameter depth is the number of layers of interpolated points to be set.
        If depth > 0 the fringe of interpolated points is set outside the blanked zones,
        whereas if depth < 0, the depth layers of blanked points are marked as to be interpolated.
        If dir=0, uses a directional stencil of depth points, if dir=1, uses a full depth x depth x depth stencil:
        Blanked points are identified by the variable 'cellN'; 'cellN' is set to 2 for the fringe of interpolated points.
        If cellN is located at cell centers, set loc parameter to 'centers', else loc='nodes'.
        

    *Using the pyTree interface:*

        ::

         t = X.setHoleInterpolatedPoints(t, depth=2, dir=0, loc='centers', cellNName='cellN')
        
        Compute the fringe of interpolated points around a set of blanked points in a pyTree t.
        Parameter depth is the number of layers of interpolated points that are built; if depth > 0 the fringe of interpolated points is outside the blanked zones, and if depth < 0,
        it is built towards the inside.
        Blanked points are identified by the variable 'cellN' located at mesh nodes or centers. 'cellN' is set to 2 for the fringe of interpolated points.
    

    *Example of use:*

    * `Set the fringe of interpolated points near blanked points (array) <Examples/Connector/setHoleInterpolatedPts.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setHoleInterpolatedPts.py

    * `Set the fringe of interpolated points near the blanked points (pyTree) <Examples/Connector/setHoleInterpolatedPtsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setHoleInterpolatedPtsPT.py


-------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.optimizeOverlap()

    *Using the array interface:*

        ::

         cellns = X.optimizeOverlap(nodes1, centers1, nodes2, centers2, prio1=0, prio2=0)

        Optimize the overlap between two zones defined by nodes1 and nodes2, centers1 and centers2 correspond to the mesh located at centers and the field 'cellN'.
        The field 'cellN' located at centers is set to 2 for interpolable points.
        Priorities can be defined for zones: prio1=0 means that the priority of zone 1 is high.
        If two zones have the same priority, then the cell volume criterion is used to set the cellN to 2 for one of the overlapping cells, the other not being modified.
        If the priorities are not specified, the cell volume criterion is applied also:

    *Using the pyTree interface:*

        ::

         t = X.optimizeOverlap(t, double_wall=0, priorities=[], planarTol=0.)

        Optimize the overlapping between all structured zones defined in a pyTree t.
        The 'cellN' variable located at cell centers is modified, such that cellN=2 for a cell interpolable from another zone. 
        Double wall projection technique is activated if 'double_wall'=1. Parameter planarTol can be useful for double wall cases, 
        in the case when double wall surfaces are planar but distant from planarTol to each other.  
        The overlapping is optimized between zones from separated bases, and is based on a priority to the cell of smallest size.
        One can impose a priority to a base over another base, using the list priorities.
        For instance, priorities = ['baseName1',0, 'baseName2',1] means that zones from base of name 'baseName1' are preferred over
        zones from base of name 'baseName2':
    
    *Example of use:*

    * `Optimize overlapping (array) <Examples/Connector/optimizeOverlap.py>`_:

    .. literalinclude:: ../build/Examples/Connector/optimizeOverlap.py

    * `Optimize overlapping (pyTree) <Examples/Connector/optmimizeOverlapPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/optimizeOverlapPT.py

-----------------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.maximizeBlankedCells(a, depth=2, dir=1)


    Change useless interpolated points status (2) to blanked points (0).
    If dir=0, uses a directional stencil of depth points, if dir=1,
    uses a full depth x depth x depth stencil.
    
    *Example of use:*

    * `Maximize blanked cells (array) <Examples/Connector/maximizeBlankedCells.py>`_:

    .. literalinclude:: ../build/Examples/Connector/maximizeBlankedCells.py

    * `Optimize overlapping (pyTree) <Examples/Connector/maximizeBlankedCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/maximizeBlankedCellsPT.py



-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.setDoublyDefinedBC(a, cellN, listOfInterpZones, listOfCelln, range, depth=2)

    
    When a border of zone z is defined by doubly defined BC in range=[i1,i2,j1,j2,k1,k2],
    one can determine whether a point is interpolated or defined by the physical BC. The array cellN defines the cell nature field at centers for zone z.
    If a cell is interpolable from a donor zone, then the cellN is set to 2 for this cell.
    The lists listOfInterpZones and listOfCelln are the list of arrays defining the interpolation domains, and corresponding cell nature fields. depth can be 1 or 2. If case of depth=2,
    if one point of the two layers is not interpolable, then celln is set to 1 for both points:
   
    *Example of use:*

    * `Set interpolated/BC points on doubly defined BCs (array) <Examples/Connector/setDoublyDefinedBC.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setDoublyDefinedBC.py

    * `Set interpolated/BC points on doubly defined BCs (pyTree) <Examples/Connector/setDoublyDefinedBCPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setDoublyDefinedBCPT.py


-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.blankIntersectingCells(a, cellN, tol=1.e-10)


    Blank intersecting cells of a 3D mesh. Only faces normal to k-planes for structured meshes and faces normal to triangular faces for prismatic meshes, and faces normal to 1234 and 5678 faces for hexahedral meshes are tested.
    The cellN is set to 0 for intersecting cells/elements. Input data are A the list of meshes, cellN the list of cellNatureField located at cell centers.
    Array version: the cellN must be an array located at centers, defined separately

    Blank intersecting cells of a 3D mesh. Only faces normal to k-planes for structured meshes and faces normal to triangular faces for prismatic meshes,
    and faces normal to 1234 and 5678 faces for hexahedral meshes are tested. Set the cellN to 0 for intersecting cells/elements. Input data are A the list of meshes, cellN the list of cellNatureField located at cell centers:
    The cellN variable is defined as a FlowSolution#Center node. The cellN is set to 0 for intersecting and negative volume cells:
    
    a = X.blankIntersectingCells(a, tol=1.e-10, depth=2)

    *Example of use:*

    * `Blank intersecting cells (array) <Examples/Connector/blankIntersectingCells.py>`_:

    .. literalinclude:: ../build/Examples/Connector/blankIntersectingCells.py

    * `Blank intersecting cells (pyTree) <Examples/Connector/blankIntersectingCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/blankIntersectingCellsPT.py

 
-----------------------------------------------------------------------------------------------------------------------------
    
.. py:function:: Connector.getIntersectingDomains(t, t2=None, method='AABB', taabb=None, tobb=None, taabb2=None, tobb2=None)
    
    Create a Python dictionary describing the intersecting zones. If t2 is not provided, then the computed dictionary states the self-intersecting zone names, otherwise, it
    computes the intersection between t and t2. Mode can be 'AABB', for Axis-Aligned Bounding Box method, 'OBB' for Oriented Bounding Box method, or 'hybrid', using a combination
    of AABB and OBB which gives the most accurate result. Depending on the selected mode, the user can provide the corresponding AABB and/or OBB PyTrees of t and/or t2, so that the
    algorithm will reuse those BB PyTrees instead of calculating them.
    
    *Example of use:*

    * `Create intersection dictionary (array) <Examples/Connector/getIntersectingDomains.py>`_:

    .. literalinclude:: ../build/Examples/Connector/getIntersectingDomainsAABB.py

    * `Create intersection dictionary (pyTree) <Examples/Connector/getIntersectingDomainsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/getIntersectingDomainsPT.py

-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.getCEBBIntersectingDomains(A, B, sameBase)

    Detect the domains defined in the list of bases B whose CEBB
    intersect domains defined in base A.
    Return the list of zone names for each basis.
    If sameBase=1, the intersecting domains are also searched in base:
 
    *Example of use:*

    * `detect CEBB intersection between bases (pyTree) <Examples/Connector/getCEBBIntersectingDomainsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/getCEBBIntersectingDomainsPT.py

-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: X.getCEBBTimeIntersectingDomains(base, func, bases, funcs, inititer=0, niter=1, dt=1, sameBase)

    in a Chimera pre-processing for bodies in relative motion,
    it can be useful to determine intersecting domains at any iteration.
    niter defines the number of iterations on which CEBB intersections are detected, starting from iteration inititer. 
    dt defines the timestep. 
    func defines a python function defining the motion of base, funcs is the list of python functions describing motions for bases.

    .. warning::  1. motions here are only relative motions. If all bases are translated with the same translation motion, it must not be defined in func.
             2. If no motion is defined on a basis, then the corresponding function must be []:

    *Example of use:*

    * `CEBB intersection between bases with motions (pyTree) <Examples/Connector/getCEBBTimeIntersectingDomainsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/getCEBBTimeIntersectingDomainsPT.py

-----------------------------------------------------------------------------------------------------------------------------

.. py:function::  X.applyBCOverlaps(t, depth=2, loc='centers')

    set the cellN to 2 for the fringe nodes or cells (depending on parameter 'loc'='nodes' or 'centers') near the overlap borders defined in the pyTree t.
    Parameter 'depth' defines the number of layers of interpolated points.

    *Example of use:* 

    * `set cellN to 2 near overlap BCs in a pyTree (pyTree) <Examples/Connector/applyBCOverlapsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/applyBCOverlapsPT.py

-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: X.setDoublyDefinedBC(t, depth=2)

    when a border is defined by doubly defined BC, one can determine whether a point is interpolated or defined by the physical BC. The cellN is set to 
    2 if cells near the doubly defined BC are interpolable from a specified donor zone:

    *Example of use:* 

    * `set interpolated/BC points on doubly defined BCs (pyTree) <Examples/Connector/setDoublyDefinedBCPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setDoublyDefinedBCPT.py

-----------------------------------------------------------------------------------------------------------------------------


.. py:function:: Connector.PyTree.cellN2OversetHoles(t)

    Compute the OversetHoles node into a pyTree from the cellN field, located at nodes or centers. For structured zones, defines it as a list of ijk indices, located at nodes or centers.
    For unstructured zones, defines the OversetHoles node as a list of indices ind, defining the cell vertices that are of cellN=0 if the cellN is located at nodes, and defining the cell centers that are of cellN=0 if the cellN is located at centers.

    The OversetHoles nodes can be then dumped to files, defined by the indices of blanked nodes or cells.

    *Example of use:*

    * `Create overset hole nodes (pyTree) <Examples/Connector/cellN2OversetHolesPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/cellN2OversetHolesPT.py

    * `Dump overset hole nodes to file (pyTree) <Examples/Connector/cellN2OversetHolesPT2.py>`_:

    .. literalinclude:: ../build/Examples/Connector/cellN2OversetHolesPT2.py
   
-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.PyTree.setInterpData(aR, aD, double_wall=0, order=2, penalty=1, nature=0, loc='nodes', storage='direct', topTreeRcv=None, topTreeDnr=None, sameName=0)

    Compute and store in a pyTree the interpolation information (donor and receptor points, interpolation type, interpolation coefficients) given receptors defined by aR,
    donor zones given by aD. If storage='direct', then aR with interpolation data stored in receptor zones are returned, and if storage='inverse', then aD with interpolation data stored in donor zones are returned.
    Donor zones can be structured or unstructured TETRA. receptor zones can be structured or unstructured.

    Interpolation order can be 2, 3 or 5 for structured donor zones, only order=2 for unstructured donor zones is performed.

    Parameter loc can 'nodes' or 'centers', meaning that receptor points are zone nodes or centers.

    penalty=1 means that a candidate donor cell located at a zone border is penalized against interior candidate cell.

    nature=0 means that a candidate donor cell containing a blanked point(cellN=0) is not valid. If nature=1 all the nodes of the candidate donor cell must be cellN=1 to be valid.

    double_wall=1 activates the double wall correction. If there are walls defined by families in aR or aD, the corresponding top trees topTreeRcv or/and topTreeDnr must be defined.

    If sameName=1, interpolation from donor zones with the same name as
    receptor zones are avoided.

    .. warning:: currently, no periodic Chimera taken into account by this function automatically.
    
    Interpolation data are stored as a ZoneSubRegion_t node, stored under the donor or receptor zone node depending of the storage.
    
    *Example of use:* 

    * `Compute interpolation connectivity (pyTree) <Examples/Connector/setInterpDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setInterpDataPT.py 

-----------------------------------------------------------------------------------------------------------------------------


.. py:function:: Connector.PyTree.getOversetInfo(tR, tD, type='interpolated')

    Set information on Chimera connectivity, i.e. interpolated, extrapolated or orphan cells, donor aspect ratio and ratio between volume of donor and receptor cells.
    This function is compliant with the storage as defined for setInterpData function.
    If type='interpolated', variable 'interpolated' is created and is equal to 1 for interpolated and extrapolated points, 0 otherwise.
    If type='extrapolated', variable 'extrapolated' is created and its value is the sum of the absolute values of coefficients, 0 otherwise.
    If type='orphan', variable 'orphan' is created and is equal to 1 for orphan points,  0 otherwise.
    If type='cellRatio', variable 'cellRatio' is created and is equal to max(volD/volR,volR/volD) for interpolated and extrapolated points (volR and volD are volume of receptors and donors).
    If type='donorAspect', variable 'donorAspect' is created and is equal to the ratio between the maximum and minimum length of donors, and 0 for points that are not interpolated.

    *Example of use:* 

    * `Get overset information (pyTree) <Examples/Connector/getOversetInfoPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/getOversetInfoPT.py


-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: Connector.PyTree.extractChimeraInfo(t, type='interpolated', loc='centers')

    Extract some Chimera information (interpolated, extrapolated, orphan points or extrapolated points with a sum of coefficients greater than a given value).
    Function chimeraInfo or oversetInfo must be applied first to compute the interpolated/extrapolated/orphan fields.
    Information is extracted as 'NODE'-type zones, whose names are suffixed by the original zone names.
    If no zone is extracted, returns []. 
    Location loc must be compliant with the location of interpolation data (i.e. location of Chimera receptor points).

    If type='interpolated', interpolated and extrapolated points are extracted.
    If type='extrapolated', extrapolated points are extracted.
    If type='orphan', orphan points are extracted. 
    If type='cf>value',  extrapolated points where the sum of absolute coefficients is greater than value are extracted.

    *Example of use:* 

    * `Extract chimera points (pyTree) <Examples/Connector/extractChimeraInfoPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/extractChimeraInfoPT.py

-----------------------------------------------------------------------------------------------------------------------------

Overset grid connectivity for elsA solver
------------------------------------------

.. py:function:: X.setInterpolations(t, loc='cell', double_wall=0, storage='inverse', prefixFile='', sameBase=0, solver='elsA', nGhostCells=2, parallelDatas=[], cfMax=30., planarTol=0., check=True)

    This function is specific to elsA solver.
    Set the Chimera connectivity (EX points and cell centers to be interpolated, index for donor interpolation cell and interpolation coefficients). 
    Double wall projection technique is activated if 'double_wall=1'.
    Parameter planarTol can be useful for double wall cases, in the case when double wall surfaces are planar but distant from planarTol to each other.  
    Parameter 'sameBase=1' means that donor zones can be found in the same base as the current zone. 

    parallelDatas=[graph,rank,listOfInterpCells] is defined only in a coupling context. It contains the graph of communication, the rank of the current processor and the list of interpolated cells/faces indices.  graph is a Python dictionary with the following structure : graph[proc1][proc2] gives the list of zones on processeur 1 which intersect with zones on processor 2. 
    loc='cell' or 'face' indicates the location of the interpolated points (cell center or face).
    Interpolations with location at 'face' correspond to interpolation for EX points according to elsA solver.

    storage can be 'direct' (interpolation data stored in receptor zones) or 'inverse' (data stored in donor zones). Value storage='direct' is only valid if Chimera connectivity is read by elsA solver as Chimera connectivity files.

    In a distributed mode, the storage must be necessarily 'inverse'. If the Chimera connectivity is read by elsA directly as CGNS ZoneSubRegion_t nodes, then the storage must be 'inverse' too. 
    prefixFile is the prefix for the name of the connectivity files generated for elsA solver (solver='elsA') or Cassiopee solver (solver='Cassiopee'). If prefixFile is not specified by the user, no Chimera connectivity file is dumped. 

    nGhostCells is the number of ghost cells that are required by elsA solver when writing connectivity files. Can be greater than 2 only if elsA reads the Chimera connectivity files.<br>
    cfMax is a threshold value for a valid extrapolation: if the sum of the extrapolated coefficients is greater than cfMax, then the cell is marked as orphan.<br>
    parallelDatas: used to perform setInterpolations in a distributed context : a list of communication information [graph, rank, interpPts], where interpPts is the list of interpolated cells/EX points. 

    Parameter check is a Boolean which displays the summary of interpolated, extrapolated and orphan points. 

    Exists also as an in-place function (X._setInterpolations):

    *Example of use:* 

    * `set interpolations  (pyTree) <Examples/Connector/setInterpolationsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setInterpolationsPT.py

-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: X.chimeraTransfer(t, storage='direct', variables=[], loc='cell', mesh='extended')

    compute Chimera transfers. This function is compliant with the storage as it is defined for setInterpolations function.

    Parameter storage can be 'direct' or 'inverse' and must be consistent with the storage computed by setInterpolations

    Parameter 'variables' specifies the variables for which the transfer is applied.

    Parameter loc can be 'cell' or 'face' to transfer the variables at cell centers or EX points.

    Parameter mesh can be 'standard' or 'extended'. For elsA simulations, it is mandatory to use mesh='extended' and storage='inverse'.

    Exists also as an in-place function (X._chimeraTransfer)

    *Example of use:* 

    * `compute Chimera transfers (pyTree) <Examples/Connector/chimeraTransferPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/chimeraTransferPT.py

-----------------------------------------------------------------------------------------------------------------------------

.. py:function:: X.chimeraInfo(t, type='interpolated')

    set information on Chimera connectivity, i.e. interpolated, extrapolated or orphan cells, donor aspect ratio and ratio between
    volume of donor and receptor cells.

    This function is compliant with the storage as it is defined for setInterpolations function.

    If type='interpolated', variable 'centers:interpolated' is created and is equal to 1 for interpolated and extrapolated cells, 0 otherwise.

    If type='extrapolated', variable 'centers:extrapolated' is created and its value is the sum of the absolute values of coefficients, 0 otherwise.

    If type='orphan', variable 'centers:orphan' is created and is equal to 1 for orphan cells,  0 otherwise.

    If type='cellRatio', variable 'centers:cellRatio' is created and is equal to max(volD/volR,volR/volD) for interpolated and extrapolated cells (volR and volD are volume of receptor and donor cells).

    If type='donorAspect', variable 'centers:donorAspect' is created and is equal to the ratio between the maximum and minimum length of donor cells, and 0 for cells that are not interpolated.

    Exists also as an in-place function (X._chimeraInfo)

    *Example of use:* 

    * `set Chimera information (pyTree) <Examples/Connector/chimeraInfoPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/chimeraInfoPT.py

-----------------------------------------------------------------------------------------------------------------------------

Immersed boundary (IBM) pre-processing
------------------------------------------

.. py:function:: X.setIBCData(aR, aD, order=2, penalty=0, nature=0, method='lagrangian', loc='nodes', storage='direct', he=0., hi=0., dim=3)

    Compute and store IBM information (donor and receptor points, interpolation type, interpolation coefficients, coordinates of corrected, wall and interpolated points) given receptors defined by aR,
    donor zones given by aD.

    * If storage='direct', then aR with interpolation data stored in receptor zones are returned, and if storage='inverse', then aD with interpolation data stored in donor zones are returned.


    * Donor zones can be structured or unstructured TETRA. receptor zones can be structured or unstructured ;

    * Interpolation order can be 2, 3 or 5 for structured donor zones, only order=2 for unstructured donor zones is performed ;

    * Parameter loc can 'nodes' or 'centers', meaning that receptor points are zone nodes or centers ;

    * penalty=1 means that a candidate donor cell located at a zone border is penalized against interior candidate cell ;

    * nature=0 means that a candidate donor cell containing a blanked point(cellN=0) is not valid. If nature=1 all the nodes of the candidate donor cell must be cellN=1 to be valid ;

    * Interpolation data are stored as a ZoneSubRegion_t node, stored under the donor or receptor zone node depending of the storage ;

    * aR must contain information about distances and normals to bodies, defined by 'TurbulentDistance','gradxTurbulentDistance','gradyTurbulentDistance' and 'gradzTurbulentDistance', located at nodes or cell centers ;

    * Corrected points are defined in aR as points with cellN=2, located at nodes or cell centers ;

    * Parameter he is a constant, meaning that the interpolated points are pushed away of a distance he from the IBC points if these are external to the bodies ;

    * Parameter hi is a constant. If hi=0., then the interpolated points are mirror points of IBC points. If hi>0., then these mirror points are then pushed away of hi from their initial position ;

    * hi and he can be defined as a field (located at nodes or centers) for any point.

     *Example of use:* 

    * `set IBM data in the pyTree (pyTree) <Examples/Connector/setIBCDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setIBCDataPT.py

------------------------------------------------------------------------------------------------

.. py:function:: X.ToolboxIBM.prepareIBMData(t, tb, DEPTH=2, loc='centers', frontType=1)

    Compute and store all the information required for IBM computations. 
    For Euler computations, corrected points are inside body, for viscous computations,
    corrected points are outside body.

    :param t:  pyTree defining the computational domain as a structured mesh
    :type  t:  pyTree 
    :param tb: pyTree defining the obstacles. Each obstacle must be defined in a CGNS basis as a surface closed mesh, whose normals must be oriented towards
               the fluid region
    :type  tb: pyTree 
    :param DEPTH: number of layers of IBM corrected points (usually 2 meaning 2 ghost cells are required)
    :type DEPTH: integer
    :param loc: location of IBM points: at nodes if loc='nodes' or at cell centers if loc='centers' ('centers'default)
    :type loc: string
    :param frontType: 0: constant distance front, 1: minimum distance front, 2: adapted distance front
    :type frontType: integer [0-2]
    :return: t, tc
    :rtype: pyTree

    The problem dimension (2 or 3D) and the equation model (Euler, Navier-Stokes or RANS) are required in tb.

    Output datas are:

    - the pyTree t with a 'cellN' field located at nodes or centers, marking as computed/updated/blanked points for the solver (cellN=1/2/0).

    - a donor pyTree tc storing all the information required to transfer then the solution at cellN=2 points:
   
    *Example of use:* 

    * `compute the IBM preprocessing (pyTree) <Examples/Connector/prepareIBMDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/prepareIBMDataPT.py

    .. warning:: requires Connector.ToolboxIBM.py module.

------------------------------------------------------------------------------------------------

.. py:function:: X.ToolboxIBM.extractIBMInfo(a)

    Extract the IBM particular points once the IBM data is computed and stored in a pyTree. These points are the IBM points that are marked
    as updated points for the IBM approach and the corresponding wall and interpolated (in fluid) points.

    If information is stored in the donor pyTree tc, then a=tc, else a must define the receptor pyTree t.

    *Example of use:* 

    * `extract all the IBM points (pyTree) <Examples/Connector/extractIBMInfoPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/extractIBMInfoPT.py

    .. warning:: requires Connector.ToolboxIBM.py module.

------------------------------------------------------------------------------------------------

.. py:function:: X.ToolboxIBM.extractIBMWallFields(a,tb=None)

    Extract the solution at walls. If IBM data is stored in donor pyTree tc, then a must be tc, else a is the pyTree t.

    If tb is None, then tw is the cloud of wall points. If tb is a triangular surface mesh, then the solution extracted at cloud points is interpolated
    on the vertices of the triangular mesh. a must contain the fields in the ZoneSubRegions of name prefixed by 'IBCD'. 

    *Example of use:* 

    * `Extract the solution on the wall  (pyTree) <Examples/Connector/extractIBMWallFieldsPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/extractIBMWallFieldsPT.py
 
    .. warning:: requires Connector.ToolboxIBM.py module.

-----------------------------------------------------------------------------------------------------------------------------

Overset and Immersed Boundary transfers with pyTrees
====================================================

    The following function enables to update the solution at some points, marked as interpolated for overset and IBM approaches.


 .. py:function:: X.setInterpTransfers(aR, topTreeD, variables=None, variablesIBC=['Density', 'MomentumX', 'MomentumY', 'MomentumZ', 'EnergyStagnationDensity'], bcType=0, varType=1, storage='unknown')

    General transfers from a set of donor zones defined by topTreeD to receptor zones defined in aR.

    Both Chimera and IBM transfers can be applied and are identified by the prefix name of the ZoneSubRegion node
    created when computing the overset or IBM interpolation data.

    Parameter variables is the list of variable names that are transfered by Chimera interpolation.

    Parameter variablesIBC defines the name of the 5 variables used for IBM transfers.

    Parameter bcType can be 0 or 1 (see table below for details).

    Parameter varType enables to define the meaning of variablesIBC, if their name is not standard (see table below for more details).

    Parameter storage enables to define how the information is stored (see table below).

    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |  Parameter value        |  Meaning                                                                                                  | 
    +=========================+===========================================================================================================+
    |    bcType=0             | IBM transfers model slip conditions                                                                       | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |    bcType=1             | IBM transfers model no-slip conditions                                                                    | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |    varType=1            | Density,MomentumX,MomentumY,MomentumZ,EnergyStagnationDensity                                             | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |    varType=2            | Density,VelocityX,VelocityY,VelocityZ,Temperature                                                         | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |    varType=3            | Density,VelocityX,VelocityY,VelocityZ,Pressure                                                            | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |    storage=0            | Interpolation data is stored in receptor zones aR                                                         | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |    storage=1            | Interpolation data is stored in donor zones topTreeD                                                      | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |    storage=-1           | Interpolation data storage is unknown or can be stored in donor and receptor zones.                       | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+
    |    extract=1            | Wall fields are stored in zone subregions (density and pressure, utau and yplus if wall law is applied).  | 
    +-------------------------+-----------------------------------------------------------------------------------------------------------+

    Exists also as an in-place function (X._setInterpTransfers):

     *Example of use:* 

    * `Transfers the solution from donor zones to receptor zones (pyTree) <Examples/Connector/setInterpTransfersPT.py>`_:

    .. literalinclude:: ../build/Examples/Connector/setInterpTransfersPT.py
    
-----------------------------------------------------------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2   

Index
###################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

