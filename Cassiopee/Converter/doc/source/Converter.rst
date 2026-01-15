.. Converter documentation master file

:tocdepth: 2

Converter: CFD data conversion module
=========================================

Preamble
########

This module provides functions for CFD data conversion (both file format
and grid topology).

This module is part of Cassiopee, a free open-source
pre- and post-processor for CFD simulations.

This module can manipulate two different data structures:
the first one is called an **array** and is very close to a numpy array,
the second one is called a **pyTree** and it implements the CGNS/python standard.

- An **array** is a simple definition of a mesh using classical numpy arrays. 

An array can be a structured array defined by a python list **[ 'x,y,z,...', an, ni, nj, nk ]**, where *ni, nj, nk* are the dimension of the grid and *an* is a (*nfld*, *nixnjxnk*) numpy array containing data (coordinates and fields). 

An array can also be an unstructured array 
defined by **[ 'x,y,z,...', an, cn, 'ELTTYPE' ]**, where *cn* is a numpy array 
storing the elements-to-nodes connectivity and *an* is a numpy array of data. 
If *an* stores fields on nodes, 'ELTTYPE' can be **'NODE', 'BAR', 'TRI', 'QUAD', 'TETRA', 'PYRA', 'PENTA', 'HEXA'**. 
If *an* stores field on elements, 'ELTTYPE' can be **'NODE*', 'BAR*', 'TRI*', 'QUAD*', 'TETRA*', 'PYRA*', 'PENTA*', 'HEXA*'**. 
Finally, an unstructured array can be 
of type **'NGON'**, describing meshes made of polyhedral elements. For those arrays, 
the connectivity *cn* is made of a faces-to-nodes connectivity (FN) and a 
elements-to-faces connectivity (EF). *cn* is then a flat numpy array 
[nfaces, sizeofFN, ..FN.., nelts, sizeofEF, ..EF..], where nfaces is the number 
of faces in mesh, nelts the number of elements in mesh. For each face, FN 
is [number of nodes, ..nodes indices..]. For each element, EF is 
[number of faces, ..face indices..].

To use the array interface::

    import Converter as C

- A **pyTree** is a CGNS/Python tree, that is a mapping of the CGNS standard in Python, using lists and numpy arrays. 

Each node of the tree is a Python list defined by 
**[ 'name', ar, [...], 'CGNSType_t' ]**, where *ar* is the value stored by this 
node (*ar* is a numpy array), and [...] 
designates a list of nodes that are the children of the current node.

**Important note**: Numpy arrays stored in pyTrees stores only one variable 
for each node. Numpy arrays for a structured zone can be accessed by *ar* [i,j,k] and by *ar* [ind] for a unstructured zone.

To use the pyTree interface::

    import Converter.PyTree as C

.. py:module:: Converter

Name standardisation
#####################

Some functions of Converter, Post and other modules perform specific treatments for 
given variables. For instance, the computeVariables function in the Post module can 
compute the pressure automatically if density and velocity are defined with their 
CGNS names. Recognised names are CGNS names, but some alternative names are also 
recognised. Naming convention is described in the following table. 

+--------------------------------------------------------------------+----------------------------------+-------------------+
| Description                                                        | CGNS                             | Altenative names  |
+====================================================================+==================================+===================+
| Coordinate in x direction                                          | CoordinateX                      | x, X              |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Coordinate in y direction                                          | CoordinateY                      | y, Y              |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Coordinate in z direction                                          | CoordinateZ                      | z, Z              |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Density                                                            | Density                          | ro                |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Momentum in x direction                                            | MomentumX                        | rou, rovx         |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Momentum in y direction                                            | MomentumY                        | rov, rovy         |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Momentum in z direction                                            | MomentumZ                        | row, rovz         |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Density times total energy                                         | EnergyStagnationDensity          | roE               |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Density times turbulence kinetic energy                            | TurbulentEnergyKineticDensity    | rok               |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Density times dissipation rate of turbulence kinetic energy        | TurbulentDissipationDensity      | roeps             |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Static pressure                                                    | Pressure                         | p, P              |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Dynamic pressure                                                   | PressureDynamic                  |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Enthalpy                                                           | Enthalpy                         |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Entropy                                                            | Entropy                          |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Stagnation pressure                                                | PressureStagnation               |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Stagnation temperature                                             | TemperatureStagnation            |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| x-component of the absolute velocity                               | VelocityX                        | vx, u             |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| y-component of the absolute velocity                               | VelocityY                        | vy, v             |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| z-component of the absolute velocity                               | VelocityZ                        | vz, w             |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Absolute velocity magnitude                                        | VelocityMagnitude                |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Absolute mach number                                               | Mach                             |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Molecular viscosity                                                | ViscosityMolecular               |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Cell nature field (0:blanked, 1:discretised, 2:interpolated)       |                                  | cellN, celln      |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Cell nature field (0:blanked, 1:discretised, -Id interp block)     |                                  | cellNF, cellnf    |
+--------------------------------------------------------------------+----------------------------------+-------------------+
| Cell nature field (-1:orphan, 0:blanked, 1:discretised,            |                                  | status            |
| 2: interpolated explicitely, 3: extrapolated, 4: interp implicit)  |                                  |                   |
+--------------------------------------------------------------------+----------------------------------+-------------------+



List of functions
##################


**-- Array creation and manipulations**

.. autosummary::
   :nosignatures:

    Converter.array
    Converter.getValue
    Converter.setValue
    Converter.addVars
    Converter.copy


**-- PyTree creation and manipulations**

.. autosummary::
   :nosignatures:

    Converter.PyTree.newPyTree
    Converter.PyTree.addBase2PyTree
    Converter.PyTree.getNobOfBase
    Converter.PyTree.getNobNozOfZone
    Converter.PyTree.breakConnectivity
    Converter.PyTree.mergeConnectivity
    Converter.PyTree.mergeByEltType
    Converter.PyTree.sliceNGonFaces
    Converter.PyTree.deleteEmptyZones
    Converter.PyTree.addState
    Converter.PyTree.addChimera2Base
    Converter.PyTree.addBC2Zone
    Converter.PyTree.fillEmptyBCWith
    Converter.PyTree.rmBCOfType
    Converter.PyTree.rmBCOfName
    Converter.PyTree.rmBCDataVars
    Converter.PyTree.extractBCOfType
    Converter.PyTree.extractBCOfName
    Converter.PyTree.getEmptyBC
    Converter.PyTree.getBCs
    Converter.PyTree.recoverBCs
    Converter.PyTree.extractBCFields
    Converter.PyTree.getConnectedZones
    Converter.PyTree.addFamily2Base
    Converter.PyTree.tagWithFamily
    Converter.PyTree.getFamilyZones
    Converter.PyTree.getFamilyBCs
    Converter.PyTree.getFamilyZoneNames
    Converter.PyTree.getFamilyBCNamesOfType
    Converter.PyTree.getFamilyBCNamesDict
    Converter.PyTree.getValue
    Converter.PyTree.setValue
    Converter.PyTree.setPartialFields
    Converter.PyTree.addVars
    Converter.PyTree.fillMissingVariables
    Converter.PyTree.cpVars


**-- Array / PyTree common manipulations**

.. autosummary::
   :nosignatures:

    Converter.getVarNames
    Converter.isNamePresent
    Converter.getNPts
    Converter.getNCells
    Converter.initVars
    Converter.extractVars
    Converter.rmVars
    Converter.convertArray2Tetra
    Converter.convertArray2Hexa
    Converter.convertArray2NGon
    Converter.convertArray2Node
    Converter.convertBAR2Struct
    Converter.convertTri2Quad
    Converter.convertHO2LO
    Converter.convertLO2HO
    Converter.conformizeNGon
    Converter.node2Center
    Converter.center2Node
    Converter.PyTree.addGhostCells
    Converter.PyTree.rmGhostCells
    Converter.PyTree.signNGonFaces

**-- Array / PyTree analysis**

.. autosummary::
   :nosignatures:
   
    Converter.diffArrays
    Converter.getMinValue
    Converter.getMaxValue
    Converter.getMeanValue
    Converter.getMeanRangeValue
    Converter.normL0
    Converter.normL2
    Converter.normalize
    Converter.magnitude
    Converter.randomizeVar
    Converter.isFinite

**-- Array / PyTree input/output**

.. autosummary::
   :nosignatures:

    Converter.convertFile2Arrays
    Converter.convertArrays2File
    Converter.PyTree.convertFile2PyTree
    Converter.PyTree.convertPyTree2File

**-- Preconditioning**

.. autosummary::
   :nosignatures:

    Converter.createHook
    Converter.createGlobalHook
    Converter.freeHook

**-- Geometrical/topological identification**

.. autosummary::
   :nosignatures:

    Converter.identifyNodes
    Converter.identifyFaces
    Converter.identifyElements
    Converter.identifySolutions
    Converter.nearestNodes
    Converter.nearestFaces
    Converter.nearestElements

    Converter.createGlobalIndex
    Converter.recoverGlobalIndex

**-- Client/server to exchange arrays/pyTrees**

.. autosummary::
   :nosignatures:

    Converter.createSockets
    Converter.listen
    Converter.send

**-- Converter arrays/3D arrays conversion**

.. autosummary::
   :nosignatures:
    
    Converter.Array3D.convertArrays2Arrays3D
    Converter.Array3D.convertArrays3D2Arrays


Contents
########

Array creation and manipulations
---------------------------------

.. py:function:: Converter.array(vars, ni, nj, nk)

    Create a structured array containing variables x,y,z on a nixnjxnk 
    structured grid.

    :param vars: variables stored in array
    :type vars: string
    :param ni,nj,nk: grid dimensions  
    :type  ni,nj,nk:  int 
    :rtype: structured array

.. py:function:: Converter.array(vars, np, ne, eltType)

    Create a unstructured array containing variables x,y,z on a 
    unstructured grid. The grid has np points, ne elements of
    type eltType. eltType can be 'NODE', 'BAR', 'TRI', 'QUAD', 
    'TETRA', 'PYRA', 'PENTA', 'HEXA', 'NGON'. 

    :param vars: variables stored in array
    :type vars: string
    :param np: number of points of grid  
    :type  np: int
    :param ne: number of elements of grid  
    :type  ne: int
    :param eltType: type of elements
    :type eltType: string
    :rtype: unstructured array

    *Example of use:*

    * `Array creation (array) <Examples/Converter/arrayK.py>`_:

    .. literalinclude:: ../build/Examples/Converter/arrayK.py

---------------------------------------------------------------------------

.. py:function:: Converter.getValue(array, ind)

    Return the list of values defined in array for point of index ind 
    (for both structured and unstructured arrays). For structured arrays, 
    you can specify (i,j,k) instead of ind. For unstructured arrays, 
    the index ind corresponds to the location type of point defining 
    array a: for instance, if array a describes a field at element vertices, 
    ind is a vertex index (ind starts at 0 and (i,j,k) start at 1).

    :param array: input array
    :type array: array
    :param ind: index 
    :type  ind: int or tuple of int
    :rtype: list of floats corresponding to field values

    *Example of use:*

    * `Get values for a given grid index (array) <Examples/Converter/getValue.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getValue.py

---------------------------------------------------------------------------

.. py:function:: Converter.setValue(array, ind, values)

    Set the values of one point of index ind in array. 
    values must be a list corresponding to the variables stored in array.

    :param array: input array
    :type array: array
    :param ind: index
    :type  ind: int or tuple of int
    :param values: values of field to set in this point
    :type values: list of floats

    *Example of use:*

    * `Set values at a given grid index (array) <Examples/Converter/setValue.py>`_:

    .. literalinclude:: ../build/Examples/Converter/setValue.py

---------------------------------------------------------------------------

.. py:function:: Converter.addVars(array, add='Density')

    Add variable(s) to an array. Variables argument can be a string 
    name ('ro') or a list of string names (['ro', 'rou']).
    Variables are set to 0.

    :param array: input array
    :type array: [array, list of arrays]
    :param add: variable to add 
    :type  add: string or list of strings
    :rtype: array with additional variables

.. py:function:: Converter.addVars([a,b,c])

    Concatenate array fields with the same dimensions. 
    Variables defined by a list of arrays are put in the same array. 

    :param arrays: input arrays
    :type arrays: list of arrays with same dimension
    :rtype: array with all variables concanated

    *Example of addVars(array, 'Density'):*

    * `Adding variables (array) <Examples/Converter/addVar.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addVar.py

    *Example of addVars([a,b,c]):*

    * `Adding many variables (array) <Examples/Converter/addVars.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addVars.py

---------------------------------------------------------------------------

.. py:function:: Converter.copy(array)

    Copy an array (return a new duplicated array).

    :param array: input array
    :type array: [array, list of arrays]
    :rtype: identical to input

    *Example of use:*

    * `Array copy (array) <Examples/Converter/copya.py>`_:

    .. literalinclude:: ../build/Examples/Converter/copya.py

---------------------------------------------------------------------------


pyTree creation and manipulation
----------------------------------

.. py:function:: Converter.PyTree.newPyTree(args)

    Create a new pyTree. You can specify base names, cell dimension in base, and 
    attached zones eventually. See below example for all possibilities of input.

    :param args: input
    :type args: [list of baseNames, list of baseNames and list of zones]
    :rtype: a new pyTree

    *Example of use:*

    * `Create pyTree (pyTree) <Examples/Converter/newPyTree.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newPyTree.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.addBase2PyTree(a, baseName, cellDim=3)

    Add a base named 'baseName' to a pyTree. Third argument specifies the cell 
    dimension (cellDim=3 for volume meshes, cellDim=2 for surface meshes).

    Exists also as in place version (_addBase2PyTree) that modifies a and 
    returns None.

    :param a: pyTree
    :type a: CGNS pyTree node
    :param baseName: name of created base
    :type baseName: string
    :param cellDim: cell dimension of zones in base
    :type cellDim: int
    :rtype: pyTree with new base added

    *Example of use:*

    * `Add base to pyTree (pyTree) <Examples/Converter/addBase2PyTree.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addBase2PyTree.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getNobOfBase(base, t)

    Get the number of a given base in tree base list, such that t[2][nob] = base.

    :param base: a base node
    :type base: CGNS base node
    :param t: tree containing base
    :type t: pyTree
    :rtype: the no of base in t children list

    *Example of use:*

    * `Get base number (pyTree) <Examples/Converter/getNobOfBasePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNobOfBasePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getNobNozOfZone(zone, t)

    Get the number (nob, noz) of a given zone a tree base and zone list , 
    such that t[2][nob][2][noz] = zone.

    :param zone: a zone node
    :type zone: CGNS zone node
    :param t: top tree containing zone
    :type t: pyTree
    :rtype: the no of base and zone in t children list

    *Example of use:*

    * `Get zone number (pyTree) <Examples/Converter/getNobNozOfZonePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNobNozOfZonePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.breakConnectivity(a)

    Break a multi-element zone (unstructured) into single type element zones.
    If a is a zone node, return a list of single type element zones.
    If a is a base, a tree or a list of zones return a base, a tree
    or a list of zones containing single type element zones.

    :param a: Input data
    :type a: [pyTree, base, zone, list of zones]
    :rtype: list of single element zones or identical to input

    *Example of use:*

    * `Break connectivity (pyTree) <Examples/Converter/breakConnectivityPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/breakConnectivityPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.mergeConnectivity(a, b, boundary=0)

    Merge two zones (unstructured) into a single zone with a multiple 
    connectivity. If boundary=1, b will be a BC connectivity in a 
    (b must be a subzone of a), if boundary=0, b will be a element connectivity. 

    :param a: first zone
    :type a: CGNS Zone node
    :param b: second zone
    :type b: CGNS Zone node
    :param boundary: 1 if b is boundary connectivity, 0 if b is an element connectivity
    :type boundary: 0 or 1
    :rtype: single zone with multiple connectivity

    *Example of use:*

    * `Merge connectivity (pyTree) <Examples/Converter/mergeConnectivityPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/mergeConnectivityPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.mergeByEltType(a)

    Merge an unstructured array by element type, thus ensuring each element type
    is listed only once. For example, if the input zone is TRI,TRI,QUAD,
    return a TRI,QUAD zone.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :rtype: Identical to input

    *Example of use:*

    * `Merge an unstructured array by element type (pyTree) <Examples/Converter/mergeByEltTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/mergeByEltType.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.sliceNGonFaces(z, indices=None)

    Slice an NGON connectivity using a list of face indices. Return two numpy
    arrays: an array of face vertices and an array of face offsets.

    :param a: input data
    :type a: CGNS Zone node
    :param indices: face indices
    :type indices: list of integers
    :rtype: zone

    *Example of use:*

    * `Slice NGON connectivity (pyTree) <Examples/Converter/sliceNGonFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/sliceNGonFacesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.deleteEmptyZones(a)

    Delete structured zones with a null ni, nj or nk, delete unstructured 
    zones with a null number of nodes or elements.

    Exists also as in place version (_deleteEmptyZones) that modifies a and 
    returns None.

    :param a: Input data 
    :type a: [pyTree, base, list of zones]
    :rtype: Identical to input

    *Example of use:*

    * `Delete empty zones (pyTree) <Examples/Converter/deleteEmptyZones.py>`_:

    .. literalinclude:: ../build/Examples/Converter/deleteEmptyZones.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.addState(a, state, value)

    Add a FlowEquation or a ReferenceState data.

    Exists also as in place version (_addState) that modifies a and 
    returns None.

    :param a: Input data
    :type a: [pyTree, base, zone, list of zones]
    :param state: the state to add or modify
    :type state: string
    :param value: the value of state
    :type value: int, float, string, numpy
    :rtype: Identical to input


.. py:function:: Converter.PyTree.addState(a, adim='adim1', MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8, UInf=None, TInf=None, PInf=None, RoInf=None, LInf=None, Mus=None, MutSMuInf=0.2, TurbLevelInf=1.e-4, EquationDimension=None, GoverningEquations=None)

    Add a full reference state built from Adim. See Initiator documentation.

    Exists also as in place version (_addState) that modifies a and 
    returns None.

    :param a: Input data
    :type a: [pyTree, base, zone, list of zones]
    :param adim: type of adimensioning
    :type adim: string
    :param MInf, alphaZ...: data for adimensioning
    :type MInf, alphaZ...: floats 
    :rtype: pyTree with new base added

    *Example of addState(a, state, value):*

    * `Add single state (pyTree) <Examples/Converter/addStatePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addStatePT.py

    *Example of addState(a, adim, ...):*

    * `Add reference state (pyTree) <Examples/Converter/addState2PT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addState2PT.py


---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.addChimera2Base(base, setting, value)

    Add a Chimera setting to a node in base. 
    Settings are added in a .Solver#Chimera user defined node.
    Chen using chimera, a CGNS base defines one component.
    They are used to define priority in grid assembly, setting of xray blanking,
    tolerance in double wall, and the kind of relationship for assembling
    components. '+' means union, '-' means difference, '0' means inactive,
    'N': means neutral.

    Exists also as in place version (_addChimera2Base) that modifies base and 
    returns None.


    :param base: input base node 
    :type base: CGNS base node
    :param setting: type of chimera setting
    :type setting: string in ['Priority', 'XRayTol', 'XRayDelta', 'DoubleWallTol', '+', '-', '0', 'N']
    :param value: value for setting
    :type value: int, float, string
    :rtype: reference copy of input

    *Example of use:*

    * `Add Chimera to base (pyTree) <Examples/Converter/addChimera2Base.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addChimera2Base.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.addBC2Zone(a, bndName, bndType, wrange=[], zoneDonor=[], rangeDonor=[], trirac=[1,2,3], rotationCenter=[], rotationAngle=[], translation=[], faceList=[], elementList=[], elementRange=[], data=None, subzone=None, faceListDonor=None, elementListDonor=None, elementRangeDonor=None, tol=1.e-12, unitAngle=None)

    Add a physical boundary condition (BC) or a grid connectivity (GC) to a structured/basic element/NGON zone of a PyTree. Parameter bndName is the 
    name of the BC or GC. Exists also as in place version (_addBC2Zone) modifying zone and returning None.

    :param a: zone in which the BC/GC is defined
    :type a: CGNS zone node
    :param bndName: name of the BC/GC
    :type bndName: string
    :param bndType: type of BC or GC defined either by a CGNS type or by a family of BCs of name myFamilyBC. Joins between stages must be defined by a familyBC prefixed by 'BCStage' 
    :type bndType: string as a CGNS BC type or ['BCMatch','BCNearMatch','BCOverlap','FamilySpecified:'+myFamilyBC]
    :param wrange: for structured grids only. Defines the range of the BC/GC
    :type wrange: a list of integers defining the window [imin,imax,jmin,jmax,kmin,kmax] or a string in ['imin','imax','jmin','jmax,'kmin','kmax']
    :param zoneDonor: donor zone(s)
    :type zoneDonor: zone node for abutting GC and list of [zone nodes, zone names, family of zones prefixed by 'FamilySpecified:]
    :param rangeDonor: range of donor zone for abutting GC, 'doubly_defined' for a doubly defined overlap GC.
    :type rangeDonor: a list of integers defining the window [imin,imax,jmin,jmax,kmin,kmax] or a string in ['imin','imax','jmin','jmax,'kmin','kmax'] or 'doubly_defined'
    :param trirac: for an abutting GC, defines the transformation from the window of zone to the donor window  
    :type trirac: list of three signed integers as a permutation of [1,2,3]
    :param rotationCenter: for GC with periodicity by rotation, coordinates of the rotation center
    :type rotationCenter: 3-tuple of floats
    :param rotationAngle: for GC with periodicity by rotation, angles of rotation in the three directions
    :type rotationAngle: 3-tuple of floats
    :param unitAngle:  defines the units of the rotationAngle (if None, the rotation angle is assumed in degrees)
    :type  unitAngle: 'Radian','Degree',None
    :param translation: for GC with periodicity by translation
    :type translation: 3-tuple of floats
    :param faceList: list of indices of faces  for unstructured BE/NGON
    :type faceList: list of integers (>0)
    :param elementList: list of indices of elements defining the BC (unstructured basic elements only) 
    :type elementList: list of integers
    :param elementRange: range of elements referencing an existing boundary connectivity (unstructured basic elements only)
    :type elementRange: list of two integers: [rangeMin, rangeMax]
    :param data: Dirichlet data set in a BCDataSet node with name 'State'
    :type data: numpy array of floats 
    :param subzone: zone corresponding to the window where the BC is defined (for unstructured zones only)
    :type subzone: pyTree zone 
    :param faceListDonor: list of donor faces (unstructured zones only)
    :type faceListDonor: list of integers
    :param elementListDonor: list of elements defining the donor window (unstructured basic elements only)
    :type elementListDonor: list of integers
    :param elementRangeDonor: range of elements defining an existing boundary connectivity corresponding to the donor window (unstructured basic elements only)
    :type elementRangeDonor: list of two integers [rangeMin, rangeMax]
    :param tol: tolerance for abutting GC
    :type tol: float
    :rtype: reference copy of input


    .. For structured grids, parameter 'wrange' specifies the window where the BC is applied. It can be defined by a list of indices [imin, imax, jmin, jmax, kmin, kmax] or by 'imin', 'jmin', ... if the BC is defined on the whole window ('imin' stands for i=1, 'imax' for i=imax). A physical BC for a structured grid is defined by::

    ..      b = C.addBC2Zone(a, bndName, bndType, wrange, data=None) 

    .. For a 1-to-1 abutting GC, donor zone and wrange and transformation must be specified::

    ..     b = C.addBC2Zone(a, bndName, 'BCMatch', wrange, zoneDonor=donorZone, rangeDonor=donorRange, trirac=[-1,2,3]) 

    .. For periodic 1-to-1 GC, you must specify rotationCenter, rotationAngle and translation::

    ..     b = C.addBC2Zone(a, bndName, 'BCMatch', wrange, zoneDonor=donorZone, rangeDonor=donorRange, trirac=[-1,2,3], rotationCenter=(0,0,0), rotationAngle(0,0,25.), translation=(0,0,0)) 

    .. For an overlap GC, donor zones can be provided as a list of zones/zone names/family of zones. 
    .. If the window range defines an overset GC and a physical BC, then 'doubly_defined' must defined as 'rangeDonor'::

    ..     b = C.addBC2Zone(a, bndName, 'BCOverlap', wrange, zoneDonor=donorZones, rangeDonor='doubly_defined',topTreeD=None) 

    .. For basic element unstructured zones, the location of the BC/GC can be specified either by a list of faces defined by 'faceList', either by 'elementList' or 'elementRange' referencing an existing boundary connectivity, or by a subzone::

    ..     b = C.addBC2Zone(a, bndName, bndType, faceList=[], data=None).or. b = C.addBC2Zone(a, bndName, bndType, elementList=[], data=None) .or. b = C.addBC2Zone(a, bndName, bndType, elementRange=[], data=None) .or. b = C.addBC2Zone(a, bndName, bndType, subzone=z, data=None) 

    .. For NGON zones, only faceList or subzone can be used::

    ..     b = C.addBC2Zone(a, bndName, bndType, faceList=[], data=None).or. b = C.addBC2Zone(a, bndName, bndType, subzone=z, data=None) 



    *Example of use:*

    * `Add boundary condition to zone (pyTree) <Examples/Converter/addBC2Zone.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addBC2Zone.py

    * `Add boundary condition to NGON zone (pyTree) <Examples/Converter/addBC2ZoneU.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addBC2ZoneU.py


---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.fillEmptyBCWith(a, bndName, bndType, dim=3)

    Fill empty boundary conditions of grids with the 
    given boundary condition. Parameter dim can be 2 or 3.

    Exists also as in place version (_fillEmptyBCWith) that modifies a and 
    returns None.

    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :param bndName: generic name of bnd
    :type bndName: string
    :param bndType: type of bnd
    :type bndType: string
    :param dim: dimension of problem
    :type dim: 2 or 3
    :rtype: reference copy of input

    *Example of use:*

    * `Fill empty BC (pyTree) <Examples/Converter/fillEmptyBCWith.py>`_:

    .. literalinclude:: ../build/Examples/Converter/fillEmptyBCWith.py    

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.rmBCOfType(a, bndType)

    Remove all boundaries of a given type. bndType accepts wildcards.
    bndType can also be a family 
    BC name. In this case, to remove a family named 'myFamily', you
    must set bndType to 'FamilySpecified:myFamily'.

    Exists also as in place version (_rmBCOfType) that modifies a and 
    returns None.

    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :param bndType: type of bnd to remove (accepts wildcards)
    :type bndType: string
    :rtype: reference copy of input

    *Example of use:*

    * `Remove BC of type (pyTree) <Examples/Converter/rmBCOfType.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmBCOfType.py    

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.rmBCOfName(a, bndName)

    Remove all boundaries of a given name. 
    bndName accepts wildcards.
    bndName can also be a family 
    BC name. In this case, to remove a family named 'myFamily', you
    must set bndName to 'FamilySpecified:myFamily'.

    Exists also as in place version (_rmBCOfName) that modifies a and 
    returns None.

    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :param bndName: name of bnd to remove (accepts wildcards)
    :type bndName: string
    :rtype: reference copy of input

    *Example of use:*

    * `Remove BC of name (pyTree) <Examples/Converter/rmBCOfNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmBCOfNamePT.py


---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.rmBCDataVars(a, varName)

    Remove variables given by varName in all BCDataSet. 
    a can be tree, zone or list of zones. 
    varName can be single variable name or a list of variable name.

    Exists also as in place version (_rmBCDataVars) that modifies a and 
    returns None.

    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :param varName: name of variable (or list of name) to remove
    :type varName: string or list of strings
    :rtype: reference copy of input

    *Example of use:*

    * `Remove BCDataSet variable (pyTree) <Examples/Converter/rmBCDataVarsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmBCDataVarsPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.extractBCOfType(a, bndType, topTree=None, reorder=True, shift=0)

    Extract all boundaries of a given type. Returns a list of zones.
    Each zone corresponds to one boundary condition. 
    If a BCDataSet exists in boundary condition,
    it is contained as zones flow solution. If flow solution
    exists and no BCDataSet exists, the flow solution is extracted.
    bndType accepts wildcards. 
    bndType can also be a family BC name. Int this case, to extract
    the BCs of 'myFamily', you must set
    bntType to 'FamilySpecified:myFamily'. 
    
    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :param bndType: type of BC to extract (accepts wildcards)
    :type bndType: string
    :param topTree: top tree if a is a zone and contains families of BCs.
    :type topTree: pyTree
    :param reorder: if True, extracted zones are reordered such that normals are oriented towards the interior of a.
    :type reorder: Boolean
    :param shift: if not 0, shift boundary of shift cells (only structured grids)
    :type shift: int
    :rtype: list of zones

    *Example of use:*

    * `Extract BCs of type (pyTree) <Examples/Converter/extractBCOfTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/extractBCOfTypePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.extractBCOfName(a, bndName, reorder=True, shift=0)

    Extract all boundaries of a given name. Returns a list of zones.
    Each zone corresponds to one boundary condition. 
    If a BCDataSet exists in boundary condition,
    it is contained as zones flow solution. If flow solution
    exists and no BCDataSet exists, the flow solution is extracted.
    bndName accepts wildcards. 
    bndName can also be a family BC name. Int this case, to extract
    the BCs of 'myFamily', you must set
    bntName to 'FamilySpecified:myFamily'. 
    
    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :param bndName: name of BC to extract (accepts wildcards)
    :type bndName: string
    :param reorder: if True, extracted zones are reordered such that normals are oriented towards the interior of a.
    :type reorder: Boolean
    :param shift: if not 0, shift boundary of shift cells (only structured grids)
    :type shift: int
    :rtype: list of zones

    *Example of use:*

    * `Extract BC(s) of name (pyTree) <Examples/Converter/extractBCOfNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/extractBCOfNamePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getEmptyBC(a, dim=3, splitFactor=180.)

    For each zone, undefined boundary conditions is a list of
    ranges [imin,imax,jmin,jmax,kmin,kmax] of undefined boundaries 
    for structured zones or is a list of face indices 
    for unstructured zones. 
    The complete return is a list of of undefined boundary condition
    for each zone.    

    Lists can be empty ([[],...,[]]) if all the boundary conditions 
    of a zone have been defined. 
    Parameter dim can be 2 or 3. 
    For unstructured grids, undefined boundaries can be split if the 
    angle between neighbouring elements exceeds splitFactor in degrees 
    (default no split). 
    
    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :param dim: dimension of problem
    :type dim: 2 or 3
    :rtype: list of ranges or face indices

    *Example of use:*

    * `Get empty BCs (pyTree) <Examples/Converter/getEmptyBC.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getEmptyBC.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getBCs(t, reorder=True)

    Return the BCs with their complete geometries, names and types.

    :param t: input data
    :type t: [pyTree, base, zone, list of zones]
    :param reorder: if True, extracted BCs are reordered such that normals are oriented towards the interior of a.
    :type reorder: Boolean
    :rtype: tuple (BCs, BCNames, BCTypes) where BCs is a list of BC nodes, BCNames a list of BC names and BCTypes a list of BC types.

    *Example of use:*

    * `Get bcs (pyTree) <Examples/Converter/getBCsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getBCsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.recoverBCs(t, (BCs, BCNames, BCTypes), tol=1.e-11)

    Recover given BCs onto a NGon tree. BCs are given by a tuple of geometries, names and types has obtained
    by getBCs.
    Exists also as in place version (_recoverBCs) that modifies t and returns None.

    :param t: input NGon data
    :type t: [pyTree, base, zone, list of zones]
    :param (BCs, BCNames, BCTypes): tuple (BCs, BCNames, BCTypes) where BCs is a list of BC nodes, BCNames a list of BC names and BCTypes a list of BC types.
    :type (BCs, BCNames, BCTypes): tuple (list,list,list)
    :rtype: reference copy of t

    *Example of use:*

    * `Recover boundary conditions (pyTree) <Examples/Converter/recoverBCsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/recoverBCsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.extractBCFields(a,varList=None)

    Extract fields defined at BCs of a zone *z*. If no BCDataSet is defined then a 0th-order extrapolation from interior cells is done. If a BCDataSet is defined, it has priority on the extrapolation. List of variables can be specified by the user. If not, the variables that are extracted are those defined in the FlowSolution node located at cell centers. Currently, this function works for structured and NGON zones. It returns the list of variables that could have been extracted and the indices of the face centers of the corresponding BCs.


   :param a: input data
   :type a: zone 
   :param varList: list of variables to be extracted at BCs 
   :type varList: list of strings defining variables or None
   :rtype: tuple (varList, fields, indicesBC) where varList is a list of extracted variables, fields the list of numpy arrays defining extracted fields at BCs, indicesBC the numpy array of indices of BC faces.                     

   *Example of use:*
   
   * `extractBCFields (pyTree) <Examples/Converter/extractBCFieldsPT.py>`_:

   .. literalinclude:: ../build/Examples/Converter/extractBCFieldsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getConnectedZones(a, topTree, type='all')

    Get zones connected to a given zone a by 'BCMatch' or 'BCNearMatch' or 'all' 
    (defined in zone GridConnectivity). 

    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :param topTree: the pyTree containing a
    :type topTree: pyTree 
    :rtype: list of zones (shared with a)

    *Example of use:*

    * `Get connected zones (pyTree) <Examples/Converter/getConnectedZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getConnectedZonesPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.addFamily2Base(a, familyName, bndType=None)

    Add a family node to a base node of a tree. 
    The family can designates a set of zone (family of zones) or
    a set of boundary conditions (family of BCs).
    If the family designates a family BC, then bndType must be defined
    with a CGNS BC type or with 'UserDefined'.
    This family name can then be referenced in zones or in boundary conditions.

    Exists also as in place version (_addFamily2Base) that modifies a and 
    returns None.

    :param a: input data 
    :type a: [pyTree, base]
    :rtype: reference copy of a

    *Example of use:*

    * `Add family to base (pyTree) <Examples/Converter/addFamily2Base.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addFamily2Base.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.tagWithFamily(a, familyName, add=False)

    Tag zones or a BC nodes with a family name.
    If a is a pyTree, base, zone or list of zones, family is supposed
    to be a family of zones.
    If a is a BC node or a list of BC nodes, family is supposed to
    be a familyBC. 
    If add=True and a family already exists, the family is
    added as a AdditionalFamilyName.

    Exists also as in place version (_tagWithFamily) that modifies a and 
    returns None.

    :param a: input data 
    :type a: [pyTree, base, zone, list of zones, BC node, list of BC nodes]
    :param add: if True, family is added otherwise it replaces an eventual existing family
    :type add: True or False
    :rtype: reference copy of a

    *Example of use:*

    * `Tag zones or BC with family (pyTree) <Examples/Converter/tagWithFamilyPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/tagWithFamilyPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getFamilyZones(a, familyName)

    Get all zones of given family (family of zones).

    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :rtype: list of zones (shared with a)

    *Example of use:*

    * `Get family zones (pyTree) <Examples/Converter/getFamilyZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getFamilyZonesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getFamilyBCs(a, familyName)

    Get all BC nodes corresponding to a given familyName (family of BCs).

    :param a: input data 
    :type a: [pyTree, base, zone, list of zones]
    :rtype: list of BC nodes (shared with a)

    *Example of use:*

    * `Get family BC nodes (pyTree) <Examples/Converter/getFamilyBCsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getFamilyBCsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getFamilyZoneNames(a)

    Return all family zone names defined in a. 
    
    :param a: input data 
    :type a: [pyTree, base]
    :rtype: list of familyZone names

    *Example of use:*

    * `Get familyZone names (pyTree) <Examples/Converter/getFamilyZoneNamesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getFamilyZoneNamesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getFamilyBCNamesOfType(a, bndType=None)

    Return all family BC names of a given type. If type is None, 
    return all family BC names.
    
    :param a: input data 
    :type a: [pyTree, base]
    :rtype: list of familyBC names

    *Example of use:*

    * `Get familyBC names of given type (pyTree) <Examples/Converter/getFamilyBCNamesOfTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getFamilyBCNamesOfTypePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getFamilyBCNamesDict(a)

    Return all family BC names contained in a as a dictionary 'familyName', 'BCType'.
    The dictionary is dict['familyName'] = 'BCType'.
     
    :param a: input data 
    :type a: [pyTree, base]
    :rtype: dictionary of familyBC names with their type

    *Example of use:*

    * `Get familyBC names as a dictionary (pyTree) <Examples/Converter/getFamilyBCNamesDictPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getFamilyBCNamesDictPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.getValue(a, var, ind)

    Return the field value(s) defined in a zone a for point of index ind 
    (for both structured and unstructured zones). 
    For structured zones, you can specify (i,j,k) instead of ind. 
    For unstructured zones, the index ind corresponds to the location 
    type of point defining zone a. For instance, if a describes a 
    field at element vertices, ind is a vertex index. 
    var is the name of the field variable or a list of field variables
    or a container name. Variable name can be preceded with 'centers:'
    or 'nodes:'. 
    This routine is slow and must not be used to access all points of a zone. In this case,
    it is better to access the field numpy with Internal.getNodeFromName
    for example.     
    
    :param a: input zone 
    :type a: zone
    :param var: field name
    :type var: string
    :param ind: index
    :type ind: int or tuple of ints
    :rtype: float or list of floats

    *Example of use:*

    * `Get field value (pyTree) <Examples/Converter/getValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getValuePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.setValue(a, var, ind, value)

    Set the values of one point of index ind in a zone a. 
    var is the name of the field variable
    or a container name. Variable name can be preceded with 'centers:'
    or 'nodes:'. 
    value can be a float or a list of floats corresponding to the values 
    of the variables to be modified. 
    This routine is slow and must not be used to access all points of a zone. In this case,
    it is better to use setPartialFields.
    
    :param a: input zone 
    :type a: zone
    :param var: field name
    :type var: string
    :param ind: index
    :type ind: int or tuple of ints

    *Example of use:*

    * `Set field value (pyTree) <Examples/Converter/setValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/setValuePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.setPartialFields(a, F, I, loc='nodes',startFrom=0)

    Set the values for a given list of indices. Field values are given as a 
    list of arrays in F (one array for each zone), indices are given as 
    a list of numpys in I (one numpy for each zone), 
    loc can be 'nodes' or 'centers'.

    Exists also as in place version (_setPartialFields) that modifies a 
    and returns None.
    
    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param F: list of arrays of fields values
    :type F: list of arrays
    :param I: list of indices
    :type I: list of numpys
    :param loc: location of fields in zone
    :type loc: 'centers' or 'nodes'
    :param startFrom: starting indice of I (e.g. 0 or 1)
    :type startFrom: integer 
    :rtype: reference copy of input


    *Example of use:*

    * `Set values for a list of indices (pyTree) <Examples/Converter/setPartialFieldsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/setPartialFieldsPT.py

---------------------------------------------------------------------------


.. py:function:: Converter.PyTree.addVars(a, vars)

    Add given variables. Variables are added to the flow container
    as described by Internal.__FlowSolutionNodes__ or
    Internal.__FlowSolutionCenters__. 
    Prefix the variable names with 'centers:' or 'nodes:' to specify variable 
    location.
    Exists also as in place 
    version (_addVars) that modifies a and returns None.

    :param a: input data
    :type a: [pyTree, base, list of zones]
    :param vars: list of variable names to add. 
    :type vars: list of strings
    :rtype: reference copy of intput

    *Example of use:*

    * `Add variables (pyTree) <Examples/Converter/addVarsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addVarsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.fillMissingVariables(a)

    Add missing variables and reorder variables for all zones, such that 
    all zones have the same variables at the end.

    Exists also as in place version (_fillMissingVariables) that modifies a 
    and returns None.

    :param a: input data
    :type a: [pyTree, base, list of zones]
    :rtype: reference copy of intput

    *Example of use:*

    * `Fill missing variables (pyTree) <Examples/Converter/fillMissingVariablesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/fillMissingVariablesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.cpVars(a1, var1, a2, var2)

    Copy a variable from zone a1, with name var1, to zone a2, with name var2. 
    The var location must be coherent. a1 and a2 can be identical.

    Exists also as in place version (_cpVars) that modifies a2 
    and returns None.

    :param a1: input zone 1
    :type a1: zone node
    :param var1: variable name (can be preceded of 'centers:' or 'nodes:')
    :type var1: string
    :param a2: receiver zone 2
    :type a2: zone node
    :param var2: variable name (can be preceded of 'centers:' or 'nodes:')
    :type var2: string
    :rtype: reference copy of a2

    *Example of use:*

    * `Variables copy (pyTree) <Examples/Converter/cpVars.py>`_:

    .. literalinclude:: ../build/Examples/Converter/cpVars.py

---------------------------------------------------------------------------


Array / PyTree common manipulations
------------------------------------

.. py:function:: Converter.getVarNames(a, excludeXYZ=False, loc='both')

    Return the list of variable names contained in a. 
    Localization of variables 
    can be specified ('nodes', 'centers', 'both'). Coordinates
    can be excluded. Only containers defined in Internal.__GridCoordinates__,
    Internal.__FlowSolutionNodes__ and Internal.__FlowSolutionCenters__ are
    scanned.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param exludeXYZ: if True, Coordinates are not scanned
    :type excludeXYZ: True or False
    :param loc: variable localisations
    :type loc: 'nodes', 'centers', 'both'
    :rtype: list of field names (one for each zone of a)

    *Example of use:*

    * `Get variable names (array) <Examples/Converter/getVarNames.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getVarNames.py

    * `Get variable names (pyTree) <Examples/Converter/getVarNamesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getVarNamesPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.isNamePresent(a, varName)

    Return -1 if a doesn't contain field varName, 0 if at 
    least one zone in a contains varName, 1 if all zones in a 
    contain varName. 

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varName: variable name (can be preceded by 'nodes:' or 'centers:')
    :type varName: string
    :rtype: -1, 0, 1

    *Example of use:*

    * `Is variable present (array) <Examples/Converter/isNamePresent.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isNamePresent.py

    * `Is variable present (pyTree) <Examples/Converter/isNamePresentPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isNamePresentPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.getNPts(a)

    Return the total number of points in a.

    Exists also as parallel distributed version (C.Mpi.getNPts).

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :rtype: int

    *Example of use:*

    * `Get number of points (array) <Examples/Converter/getNPts.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNPts.py

    * `Get number of points (pyTree) <Examples/Converter/getNPtsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNPtsPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.getNCells(a)

    Return the total number of cells in a. 
    
    Exists also as parallel distributed version (C.Mpi.getNCells).

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :rtype: int

    *Example of use:*

    * `Get number of cells (array) <Examples/Converter/getNCells.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNCells.py

    * `Get number of cells (pyTree) <Examples/Converter/getNCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNCellsPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.initVars(a, varNameString, value, isVectorized=False)

    Initialize one or several variables as given by varNameString.

    For initialisation by a formula string, only one variable can be set at a time.

    For initialisation by a function or by a constant, varNameString can be a string
    or a list of strings.

    If the function is vectorized (can be interpreted as a numpy formula), set isVectorized to True.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varNameString: string describing variable or formula
    :type varNameString: string or list of strings
    :param value: value in case of constant init or function.
    :type value: float or function and parameters
    :param isVectorized: when using functions, indicates that function is vectorized.
    :type isVectorized: boolean
    :rtype: identical to input

    *Example of use:*

    * `Init a variable to a constant value (array) <Examples/Converter/initVars.py>`_:

    .. literalinclude:: ../build/Examples/Converter/initVars.py

    * `Init a variable with a function (array) <Examples/Converter/initVar.py>`_:

    .. literalinclude:: ../build/Examples/Converter/initVar.py

    * `Init a variable with a formula (array) <Examples/Converter/initVarsByEq.py>`_:

    .. literalinclude:: ../build/Examples/Converter/initVarsByEq.py

    .. note:: When initializing variables using string formulas, function names must match the names of the numpy library.
    
    * `Init a variable to a constant value (pyTree) <Examples/Converter/initVarsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/initVarsPT.py

    * `Init a variable with a function (pyTree) <Examples/Converter/initVarPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/initVarPT.py

    * `Init a variable with a formula (pyTree) <Examples/Converter/initVarsByEqPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/initVarsByEqPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.extractVars(a, varNames)

    Extract variables defined in varNames from a (other variables are removed). 

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varNames: names of variable to extract (can starts with 'nodes:' or 'centers:'')
    :type varNames: string or list of strings.
    :rtype: Identical to input

    *Example of use:*

    * `Extract some variables from array (array) <Examples/Converter/extractVars.py>`_:

    .. literalinclude:: ../build/Examples/Converter/extractVars.py

    * `Extract some variables from zone (array) <Examples/Converter/extractVarsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/extractVarsPT.py


-----------------------------------------------------------------------------------

.. py:function:: Converter.rmVars(a, varNames)

    Remove variable(s) from a. varNames is a string name or a list of string names.

    Exists also as in place version (_rmVars) that modifies a and returns None.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varNames: names of variable to remove (can starts with 'nodes:' or 'centers:'')
    :type varNames: string or list of strings.
    :rtype: Identical to input

    *Example of use:*

    * `Remove some variables from array (array) <Examples/Converter/rmVars.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmVars.py

    * `Remove some variables from zone (array) <Examples/Converter/rmVarsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmVarsPT.py


-----------------------------------------------------------------------------------

.. py:function:: Converter.convertArray2Tetra(a, split='simple')

    Create tetra unstructured array from an any type of mesh. 
    2D elements are made triangular, else they are made tetrahedral. 
    If split='simple', conversion does not create new points. 
    If split='withBarycenters', barycenters of elements and faces are added.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param split: 'simple' or 'withBarycenters'
    :type split: string
    :rtype: Identical to input

    *Example of use:*

    * `Convert structured mesh to tetra (array) <Examples/Converter/convertStruct2Tetra.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertStruct2Tetra.py

    * `Convert hexa mesh to tetra (array) <Examples/Converter/convertHexa2Tetra.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertHexa2Tetra.py

    * `Convert structured mesh to tetra (pyTree) <Examples/Converter/convertArray2TetraPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArray2TetraPT.py


-----------------------------------------------------------------------------------

.. py:function:: Converter.convertArray2Hexa(a)

    Create hexa unstructured array from an any type of mesh. 
    2D elements are made quadrangular, else they are made hexahedral.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :rtype: Identical to input

    *Example of use:*

    * `Convert structured mesh to hexa (array) <Examples/Converter/convertStruct2Hexa.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertStruct2Hexa.py

    * `Convert structured mesh to hexa (pyTree) <Examples/Converter/convertArray2HexaPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArray2HexaPT.py


-----------------------------------------------------------------------------------

.. py:function:: Converter.convertArray2NGon(a, recoverBC=True, api=1)

    Create NGON array from an any type of mesh.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param recoverBC: BCs can be recovered or not on the NGON a (not valid for arrays).
    :type recoverBC: boolean
    :param api: CGNSv3 compact (=1), CGNSv3 (=2) or CGNSv4 (=3).
    :type api: integer (1, 2 or 3)
    :rtype: Identical to input

    *Example of use:*

    * `Convert structured mesh to NGON (array) <Examples/Converter/convertArray2NGon.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArray2NGon.py

    * `Convert structured mesh to NGON (pyTree) <Examples/Converter/convertArray2NGonPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArray2NGonPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.convertArray2Node(a)

    Create NODE array from an any type of mesh. A node array only contains
    node and no connectivity.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :rtype: Identical to input

    *Example of use:*

    * `Convert structured mesh to NODE (array) <Examples/Converter/convertArray2Node.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArray2Node.py

    * `Convert structured mesh to NODE (pyTree) <Examples/Converter/convertArray2NodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArray2NodePT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.convertBAR2Struct(a)

    Create a structured 1D array from a BAR array. The BAR array 
    must not contain branches. To split a branched BAR, you may consider
    T.splitTBranches.

    :param a: input data (BAR)
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :rtype: Identical to input

    *Example of use:*

    * `Convert BAR to i-array (array) <Examples/Converter/convertBAR2Struct.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertBAR2Struct.py

    * `Convert BAR to i-array (pyTree) <Examples/Converter/convertBAR2StructPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertBAR2StructPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.convertTri2Quad(a, alpha=30.)

    Convert a TRI-array to a QUAD-array. Neighbouring cells with an angle 
    lower than alpha can be merged. It returns the 
    QUAD-array b and the rest of not merged cells in a TRI-array c.

    :param a: input data (TRI)
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param alpha: angle for merging
    :type alpha: float
    :rtype: Identical to input

    *Example of use:*

    * `Convert TRI to QUAD-array (array) <Examples/Converter/convertTri2Quad.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertTri2Quad.py

    * `Convert TRI to QUAD-array (pyTree) <Examples/Converter/convertTri2QuadPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertTri2QuadPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.convertHO2LO(a, mode=0)

    Convert a high order element mesh into a low order (linear) mesh.
    If mode=1, only valid for BAR_3, TRI_6, QUAD_8, QUAD_9, TETRA_10, HEXA_20, HEXA_27, PENTA_18, PYRA_14.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param mode: if 0, coarse extraction, 1, tesselate all points
    :type mode: int
    :rtype: Identical to input

    *Example of use:*

    * `Convert High order mesh to Low order mesh (array) <Examples/Converter/convertHO2LO.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertHO2LO.py

    * `Convert High order to Low order mesh (pyTree) <Examples/Converter/convertHO2LOPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertHO2LOPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.convertLO2HO(a, mode=0, order=2)

    Convert a low order element mesh into a high order mesh. Points are added linearly on edges
    or faces.
    Order 2 can give: BAR_3, TRI_6, QUAD_8, QUAD_9, TETRA_10, HEXA_20, HEXA_27, PENTA_18, PYRA_14.
    Order 3 can give: BAR_4, TRI_9, ...
    
    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param mode: specify the type of generated high order elements
    :type mode: int
    :param order: specify the order of output elements
    :type order: int
    :rtype: Identical to input

    *Example of use:*

    * `Convert Low order mesh to High order mesh (array) <Examples/Converter/convertLO2HO.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertLO2HO.py

    * `Convert Low order to High order mesh (pyTree) <Examples/Converter/convertLO2HOPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertLO2HOPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.conformizeNGon(a, tol=1.e-6)

    Conformize the cell faces of a NGon, such that a face of a cell 
    corresponds to a unique face of another cell. Typically, a
    mesh with hanging nodes will be made conform.

    :param a: input data (NGON)
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tol: tolerance for face matching
    :type tol: float
    :rtype: Identical to input

    *Example of use:*

    * `NGON mesh conformization (array) <Examples/Converter/conformizeNGon.py>`_:

    .. literalinclude:: ../build/Examples/Converter/conformizeNGon.py

    * `NGON mesh conformization (pyTree) <Examples/Converter/conformizeNGonPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/conformizeNGonPT.py


-----------------------------------------------------------------------------------

.. py:function:: Converter.node2Center(a, var='')

    Change data location from nodes to centers. 
    If no variable is specified, the mesh coordinates are also put
    to centers, resulting in a "all in nodes" mesh.
    If a variable is specified, only this variable is passed to centers
    and stored in __FlowSolutionCenters__ container.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: modified variables
    :type var: string or list of strings or container name
    :rtype: Identical to input

    *Example of use:*

    * `Nodes to centers conversion (array) <Examples/Converter/node2Center.py>`_:

    .. literalinclude:: ../build/Examples/Converter/node2Center.py

    * `Nodes to centers conversion (pyTree) <Examples/Converter/node2CenterPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/node2CenterPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.center2Node(a, var='', cellNType=0)

    Change data location from centers to nodes. 
    If no variable is specified, the mesh coordinates are also put
    to nodes, resulting in a "all in nodes" mesh.
    If a variable is specified, only this variable is passed to nodes
    and stored in __FlowSolutionNodes__ container.
    cellNType indicates the treatement for blanked points when cellN 
    field is present. cellNType=0, means that, if a node receives 
    at least one cellN=0 value from a center, its cellN is set to 0. 
    cellNType=1 means that, only if all values of neighbouring centers 
    are cellN=0, its cellN is set to 0.

    Exists also as parallel distributed version (C.Mpi.center2Node).

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variables to modify
    :type var: string or list of strings or container name
    :param cellNType: describes the type of treatment for cellN variables.
    :type cellNType: int
    :rtype: Identical to input

    *Example of use:*

    * `Centers to nodes conversion (array) <Examples/Converter/center2Node.py>`_:

    .. literalinclude:: ../build/Examples/Converter/center2Node.py

    * `Centers to nodes conversion (pyTree) <Examples/Converter/center2NodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/center2NodePT.py


-----------------------------------------------------------------------------------

.. py:function:: Converter.PyTree.addGhostCells(t, b, d, adaptBCs=1, modified=[], fillCorner=1)

    Add ghost cells to structured grids.
    if modified is given, limit add ghost cells to given field containers. Otherwise, ghost cells
    are added to all containers.
    If adaptBCs=1, Zone BCs are adapted to fit grid with ghost cells.
    If fillCorner=1, edges and corners are filled 
    according to the grid connectivity (geometrically, the corners and edges can be wrong).
    If fillCorner=0, neighbouring vectors are extrapolated to build edge cells, no filling with flow field.
    Exists also as in place version (_addGhostCells) that modifies a and returns None.

    :param t: top tree
    :type t: pyTree
    :param b: zones to modify
    :type b: [pyTree, base, zone, list of zones]
    :param d: number of layers of ghost cells to add
    :type d: int
    :param adaptBCs: if 1, zone BCs are modified to fit ghost grids.
    :type adaptBCs: 0 or 1
    :param modified: a list of containers to modify. If [], all containers are modified.
    :type modified: list of container names
    :param fillCorner: method used to fill corners
    :type fillCorner: 0 or 1
    :rtype: Identical to input b

    *Example of use:*

    * `Add ghost cells (pyTree) <Examples/Converter/addGhostCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addGhostCellsPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.PyTree.rmGhostCells(t, b, d, adaptBCs=1, modified=[])

    Remove ghost cells to structured grids. See addGhostCells.

    Exists also as in place version (_rmGhostCells) that modifies a and returns None.

    :param t: top tree
    :type t: pyTree
    :param b: zones to modify
    :type b: [pyTree, base, zone, list of zones]
    :param d: number of layers of ghost cells to add
    :type d: int
    :param adaptBCs: if 1, zone BCs are modified to fit ghost grids.
    :type adaptBCs: 0 or 1
    :param modified: a list of containers to modify. If [], all containers are modified.
    :type modified: list of container names
    :rtype: Identical to input b

    *Example of use:*

    * `Remove ghost cells (pyTree) <Examples/Converter/rmGhostCellsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmGhostCellsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.signNGonFaces(t)

    For NGON zones, sign the NFACE connectivity with cell external normals.

    Exists also as in place version (_signNGonFaces) that modifies t and returns None.

    :param t: tree
    :type t: pyTree
    :rtype: t with signed NFACE

    *Example of use:*

    * `Sign NGon faces (pyTree) <Examples/Converter/signNGonFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/signNGonFacesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.PyTree.makeParentElements(t)

    For NGON zones, construct parent elements array from NFACE connectivity.
    Always checks the validity of input mesh. Raises error if mesh is invalid.

    Exists also as in place version (_makeParentElements) that modifies t and returns None.

    :param t: tree
    :type t: pyTree
    :rtype: t with parent elements array

    *Example of use:*

    * `Construct parent elements array from NFACE connectivity (pyTree) <Examples/Converter/makeParentElementsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/makeParentElementsPT.py

---------------------------------------------------------------------------


Array / PyTree analysis
------------------------------------

.. py:function:: Converter.diffArrays(a, b, removeCoordinates=True, atol=1.e-11, rtol=0.)

    Given a solution in a and a solution in b, both defined on the same mesh, 
    return the differences.

    :param a: input data 1
    :type a: [list of arrays] or [pyTree, base, zone, list of zones]
    :param b: input data 2
    :type b: [list of arrays] or [pyTree, base, zone, list of zones]
    :param removeCoordinates: if True, remove original coordinates (pyTree)
    :type removeCoordinates: boolean
    :param atol: absolute tolerance
    :type atol: float
    :param rtol: relative tolerance
    :type rtol: float
    :rtype: Identical to input 1

    *Example of use:*

    * `Difference between two solutions (array) <Examples/Converter/diffArrays.py>`_:

    .. literalinclude:: ../build/Examples/Converter/diffArrays.py

    * `Difference between two solutions (pyTree) <Examples/Converter/diffArraysPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/diffArraysPT.py

.. -----------------------------------------------------------------------------------

.. .. py:function:: Converter.getArgMin(a, var)

    Return the field values corresponding to the node
    where variable 'var' is minimum.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name
    :type var: string
    :rtype: list of field values (array), just var value (pyTree)

    *Example of use:*

    * `Get the field value where F is minimum (array) <Examples/Converter/getArgMin.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getArgMin.py

    * `Get the field value where x is minimum (pyTree) <Examples/Converter/getArgMinPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getArgMinPT.py

.. -----------------------------------------------------------------------------------

.. .. py:function:: Converter.getArgMax(a, var)

    Return the field values corresponding to the node
    where variable 'var' is maximum.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name
    :type var: string
    :rtype: list of field values (array), just var value (pyTree)

    *Example of use:*

    * `Get the field value where F is maximum (array) <Examples/Converter/getArgMax.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getArgMax.py

    * `Get the field value where x is maximum (pyTree) <Examples/Converter/getArgMaxPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getArgMaxPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.getMinValue(a, var)

    Return the minimum value of field 'var' on input.

    Exists also as parallel distributed version (C.Mpi.getMinValue).

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name
    :type var: string
    :rtype: minimum value

    *Example of use:*

    * `Get the minimum of F (array) <Examples/Converter/getMinValue.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getMinValue.py

    * `Get the minimum of F (pyTree) <Examples/Converter/getMinValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getMinValuePT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.getMaxValue(a, var)

    Return the maximum value of field 'var' on input.

    Exists also as parallel distributed version (C.Mpi.getMaxValue).

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name
    :type var: string
    :rtype: maximum value

    *Example of use:*

    * `Get the maximum of F (array) <Examples/Converter/getMaxValue.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getMaxValue.py

    * `Get the maximum of F (pyTree) <Examples/Converter/getMaxValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getMaxValuePT.py


-----------------------------------------------------------------------------------

.. py:function:: Converter.getMeanValue(a, var)

    Return the mean value of field 'var' on input.

    Exists also as parallel distributed version (C.Mpi.getMeanValue).

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name
    :type var: string
    :rtype: mean value

    *Example of use:*

    * `Get the mean of F (array) <Examples/Converter/getMeanValue.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getMeanValue.py

    * `Get the mean of F (pyTree) <Examples/Converter/getMeanValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getMeanValuePT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.getMeanRangeValue(a, var, rmin, rmax)

    Return the mean value of variable 'var' for a given range of value.
    The field 'var' is sorted. Then the mean value of 'var' on the given range is
    returned. For instance, getMeanRangeValue(a, 'F', 0., 0.3) return the mean
    value of F for the 30% lowest values.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name
    :type var: string
    :param rmin: min of range in [0,1]
    :type rmin: float
    :param rmax: max of range in [0,1]
    :type ramx:  float
    :return: mean value of var for the min-max %
    :rtype: float

    *Example of use:*

    * `Get the mean value of x in a sorted range of values (array) <Examples/Converter/getMeanRangeValue.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getMeanRangeValue.py

    * `Get the mean value of x in a sorted range of values (pyTree) <Examples/Converter/getMeanRangeValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getMeanRangeValuePT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.normL0(a, var)

    Return the infinite norm of field 'var' on input.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name
    :type var: string
    :rtype: float

    *Example of use:*

    * `Get the L0 norm of F (array) <Examples/Converter/normL0.py>`_:

    .. literalinclude:: ../build/Examples/Converter/normL0.py

    * `Get the L0 norm of F (pyTree) <Examples/Converter/normL0PT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/normL0PT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.normL2(a, var)

    Return the L2 norm of field 'var' on input.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name
    :type var: string
    :rtype: float

    *Example of use:*

    * `Get the L2 norm of F (array) <Examples/Converter/normL2.py>`_:

    .. literalinclude:: ../build/Examples/Converter/normL2.py

    * `Get the L2 norm of F (pyTree) <Examples/Converter/normL2PT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/normL2PT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.normalize(a, ['sx','sy',sz'])

    Normalize a vector defined by its 3 vector components. 
    The vector component values are modified such that the
    vector (a.sx,a.sy,a.sz) has a unit norm for each point.

    Exists also as in place version (_normalize) that modifies a and returns None.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param sx,sy,sz: names of field used as vector components
    :type sx,sy,sz: list of strings
    :rtype: Identical to input

    *Example of use:*

    * `Normalize a vector (array) <Examples/Converter/normalize.py>`_:

    .. literalinclude:: ../build/Examples/Converter/normalize.py

    * `Normalize a vector (pyTree) <Examples/Converter/normalizePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/normalizePT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.magnitude(a, ['sx','sy',sz'])

    Get the magnitude of a vector defined by its 3 vector components for each point.
    The name of created field is composed from the components names.
    For instance 'sx,sy,sz' will create a 'sMagnitude' field.

    Exists also as in place version (_magnitude) that modifies a and returns None.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param sx,sy,sz: names of field used as vector components
    :type sx,sy,sz: list of strings
    :rtype: Identical to input

    *Example of use:*

    * `Compute a vector magnitude (array) <Examples/Converter/magnitude.py>`_:

    .. literalinclude:: ../build/Examples/Converter/magnitude.py

    * `Compute a vector magnitude (pyTree) <Examples/Converter/magnitudePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/magnitudePT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.randomizeVar(a, varName, deltaMin, deltaMax)

    Randomize a fied varName. The modified field is bounded by 
    [f-deltaMin,f+deltaMax] where f is the local field value.

    Exists also as in place version (_randomizeVar) that modifies a and returns None.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param varName: field to randomize
    :type varName: string
    :param deltaMin,deltaMax: range for random
    :type deltaMin,deltaMax: floats
    :rtype: Identical to input

    *Example of use:*

    * `Randomize a field (array) <Examples/Converter/randomizeVar.py>`_:

    .. literalinclude:: ../build/Examples/Converter/randomizeVar.py

    * `Randomize a field (pyTree) <Examples/Converter/randomizeVarPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/randomizeVarPT.py

-----------------------------------------------------------------------------------

.. py:function:: Converter.isFinite(a, var=None)

    Return True if a contains only finite values (no NAN, no INF).

    Exists also as parallel distributed version (C.Mpi.isFinite).

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param var: variable name (optional)
    :type var: string
    :rtype: True or False

    *Example of use:*

    * `Test if fields are finite (array) <Examples/Converter/isFinite.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isFinite.py

    * `Test if fields are finite (pyTree) <Examples/Converter/isFinitePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isFinitePT.py

---------------------------------------------------------------------------


.. Array / PyTree conversion
.. --------------------------------

.. .. py:function:: Converter.convertPyTree2Array(path, t)

    Convert a Python tree node to an array. One have to provide 
    the path of the corresponding node and the Python tree.

    :param path: path of node to convert
    :type path: string
    :param t: pyTree
    :type t: pyTree
    :return: array
    :rtype: Converter array

    *Example of use:*

    * `Node of pyTree to array (array) <Examples/Converter/convertPyTree2Array.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertPyTree2Array.py


Array / PyTree input/output
------------------------------

.. py:function:: Converter.convertFile2Arrays(fileName, format=None, options)

    Read a file and return a list of arrays (array interface).
    For format needing multiple files (for ex: plot3d), multiple 
    files can be specified in file name string as: "file.gbin,file.qbin".
    In file format where variables name are undefined, the following 
    ones are adopted: x, y, z, ro, rou, rov, row, roE, cellN.
    If format is unspecified, the format is guessed from file extension
    or from file header, if possible.
    For a list of available format, see FileFormats_.
    Several options are available to specify the discretization of 
    vector elements (for vector formats such as xfig or svg).
    For a list of available options, see ReadOptions_.

    :param fileName: name of file to read
    :type fileName: string
    :param format: file format (see FileFormats_)
    :type format: string
    :param options: options for vector formats such as svg, xfig (see ReadOptions_)
    :type options: keywords
    :return: list of arrays
    :rtype: list of Converter arrays

    *Example of use:*

    * `Binary tecplot file read <Examples/Converter/convertFile2Arrays.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertFile2Arrays.py


--------------------------------------------------------------------------------
 
.. py:function:: Converter.convertArrays2File(a, fileName, format=None, options)

    Write array or list of arrays to a file (array interface).
    If format is not given, it is guessed from fileName extension.
    For a list of available formats, see FileFormats_.
    For a list of available options, see WriteOptions_.

    :param a: input data
    :type a: [array, list of arrays]
    :param fileName: name of file to read
    :type fileName: string
    :param format: file format (see FileFormats_)
    :type format: string
    :param options: writing options (see WriteOptions_)
    :type options: keywords
    
    *Example of use:*

    * `Binary tecplot file read/write <Examples/Converter/convertArrays2File.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArrays2File.py


-------------------------------------------------------------------------

.. py:function:: Converter.PyTree.convertFile2PyTree(fileName, format=None, options)

    Read a file and return a CGNS pyTree (pyTree interface).
    If format is not given, it is guessed from file header or extension.
    For a list of available format, see FileFormats_.
    For a list of available options, see ReadOptions_.

    :param fileName: name of file to read
    :type fileName: string
    :param format: file format (see FileFormats_) 
    :type format: string
    :param options: reading options (see ReadOptions_)
    :type options: keywords
    :return: a pyTree
    
    *Example of use:*

    * `CGNS file read <Examples/Converter/convertFile2PyTree.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertFile2PyTree.py

-------------------------------------------------------------------------

.. py:function:: Converter.PyTree.convertPyTree2File(t, fileName, format=None, options)

    Write a pyTree to a file (pyTree interface).
    If format is not given, it is guessed from file name extension.
    For a list of available format, see FileFormats_.
    For a list of available options, see WriteOptions_.

    :param t: input data
    :type t: [pyTree, base, zone, list of zones] 
    :param fileName: name of file to read
    :type fileName: string
    :param format: file format (see FileFormats_) 
    :type format: string
    :param options: writing options (see WriteOptions_)
    :type options: keywords

    *Example of use:*

    * `CGNS file write <Examples/Converter/convertPyTree2File.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertPyTree2File.py


--------------------------------------------------------------------------

.. _FileFormats:

    Known formats for read/write functions (array and pyTree interface):

    +------------+-----------+---------------------------------------+
    | Format     | Extension | Description                           |
    +============+===========+=======================================+
    |bin_tp      | .plt      | binary tecplot file                   |
    +------------+-----------+---------------------------------------+
    |fmt_tp      | .dat, .tp | formatted tecplot file                |
    +------------+-----------+---------------------------------------+
    |bin_v3d     | .v3d      | binary v3d file (ONERA)               |
    +------------+-----------+---------------------------------------+
    |fmt_v3d     | .fv3d     | formatted v3d file (ONERA)            |
    +------------+-----------+---------------------------------------+
    |bin_plot3d  | .gbin     | binary plot 3d file (NASA)            |
    +------------+-----------+---------------------------------------+
    |fmt_plot3d  | .gfmt     | formatted plot 3d file (NASA)         |
    +------------+-----------+---------------------------------------+
    |fmt_mesh    | .mesh     | formatted mesh file (INRIA)           |
    +------------+-----------+---------------------------------------+
    |fmt_gmsh    | .msh      | formatted GMSH mesh file (UCL)        |
    +------------+-----------+---------------------------------------+
    |bin_gmsh    | .msh      | binary GMSH mesh file (UCL)           |
    +------------+-----------+---------------------------------------+
    |fmt_su2     | .su2      | formatted SU2 file (STANFORD)         |
    +------------+-----------+---------------------------------------+
    |fmt_foam    | .foam     | formatted FOAM file (OPENFOAM)        |
    +------------+-----------+---------------------------------------+
    |fmt_cedre   | .d        | formatted CEDRE file (ONERA)          |
    +------------+-----------+---------------------------------------+
    |bin_stl     | .bstl     | binary STL file                       |
    +------------+-----------+---------------------------------------+
    |fmt_stl     | .stl .fstl| formatted STL file                    |
    +------------+-----------+---------------------------------------+
    |fmt_pov     | .pov      | formatted povray raytracer file       |
    +------------+-----------+---------------------------------------+
    |fmt_selig   | .selig    | formatted selig file (airfoils)       |
    +------------+-----------+---------------------------------------+
    |fmt_obj     | .obj      | formatted OBJ file (WAVEFRONT)        |
    +------------+-----------+---------------------------------------+
    |bin_gltf    | .gltf     | binary gltf file (KHRONOS, only read) |
    +------------+-----------+---------------------------------------+
    |bin_3ds     | .3ds      | binary 3DS file (3D STUDIO)           |
    +------------+-----------+---------------------------------------+
    |bin_ply     | .ply      | binary PLY file (STANFORD)            |
    +------------+-----------+---------------------------------------+
    |bin_pickle  | .ref      | binary python pickle file             |
    +------------+-----------+---------------------------------------+
    |bin_wav     | .wav      | binary wav 8 bits sound file          |
    +------------+-----------+---------------------------------------+
    |fmt_xfig    | .fig      | formatted XFIG file                   |
    +------------+-----------+---------------------------------------+
    |fmt_svg     | .svg      | formatted SVG file (INKSCAPE)         |
    +------------+-----------+---------------------------------------+
    |bin_png     | .png      | binary PNG file                       |
    +------------+-----------+---------------------------------------+
    |bin_jpg     | .jpg      | binary JPEG file                      |
    +------------+-----------+---------------------------------------+
    |fmt_iges    | .igs      | formatted IGES CAD file               |
    +------------+-----------+---------------------------------------+
    |fmt_step    | .stp      | formatted STEP CAD file               |
    +------------+-----------+---------------------------------------+
    


    Known formats for read/write functions specific to pyTree interface:

    +------------+-----------+---------------------------------------+
    | Format     | Extension | Description                           |
    +============+===========+=======================================+
    |bin_adf     | .adf      | binary CGNS ADF file                  |
    +------------+-----------+---------------------------------------+
    |bin_hdf     | .hdf      | binary CGNS HDF file                  |
    +------------+-----------+---------------------------------------+
    |bin_cgns    | .cgns     | binary CGNS/HDF file                  |
    +------------+-----------+---------------------------------------+
    |bin_tau     | .grid     | binary TAU file                       |
    +------------+-----------+---------------------------------------+
    |bin_fsdm    | .h5       | binary FSDM file                      |
    +------------+-----------+---------------------------------------+


.. _ReadOptions: 

    Options for reading:

    +------------+---------------------------------------------------------------------+-------------------+--------------------------------------+
    | Option     | Description                                                         | Format            | Default value                        |
    +============+=====================================================================+===================+======================================+
    |nptsCurve   | Number of discretization points for curved vector elements          | svg, xfig         | 20                                   |
    +------------+---------------------------------------------------------------------+-------------------+--------------------------------------+
    |nptsLine    | Number of discretization points for lines                           | svg, xfig         | 2                                    |
    +------------+---------------------------------------------------------------------+-------------------+--------------------------------------+
    |density     | Number of discretization points per unit length                     | svg               | -1: not used. If > 0, overides npts  |
    +------------+---------------------------------------------------------------------+-------------------+--------------------------------------+
    |skipTypes   | list of strings (CGNS types) that stop reading when met             | CGNS              | None                                 |
    +------------+---------------------------------------------------------------------+-------------------+--------------------------------------+
    |hausd       | chordal error for CAD discretization                                | igs, stp          | 1. (auto)                            |
    +------------+---------------------------------------------------------------------+-------------------+--------------------------------------+
    |links       | list of list of 4 strings (see after)                               | HDF               | None                                 |
    +------------+---------------------------------------------------------------------+-------------------+--------------------------------------+  

    Links option:

    For hdf format only, when reading, links are always followed but a list of links can be returned.
    If you specify links=[] to convertFile2PyTree, a list of links is returned. Each link is a list
    ['directoryOfPointedFile', 'pointedFile', 'targetNodePath', 'currentNodePath']. The 'directoryOfPointedFile'
    is the directory where the pointed file must be found, 'pointedFile' is the pointed file name, 'targetNodePath' is the path
    of pointed node (in pointed file), 'currentNodePath' is the path of node in current file.

    *Example of use:*

    * `HDF file read with links <Examples/Converter/linksPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/linksPT.py


.. _WriteOptions: 

    Options for writing:

    +------------+--------------------------------------------------------------------------+--------------------+---------------------------------------+-----------------------------------+
    | Option     | Description                                                              | Format             | Possible values                       | Default value                     |
    +============+==========================================================================+====================+=======================================+===================================+
    |isize       | Size of integers                                                         | v3d, p3d, HDF      | 4,8                                   | 8                                 |
    +------------+--------------------------------------------------------------------------+--------------------+---------------------------------------+-----------------------------------+
    |rsize       | Size of reals                                                            | v3d, p3d, HDF      | 4,8                                   | 8                                 |
    +------------+--------------------------------------------------------------------------+--------------------+---------------------------------------+-----------------------------------+
    |endian      | Data endianess                                                           | v3d, p3d           | 'little', 'big'                       | 'big'                             |
    +------------+--------------------------------------------------------------------------+--------------------+---------------------------------------+-----------------------------------+
    |dataFormat  | 'printf' like format for formatted files (%[width][.precision]specifier) | Formatted formats  | '%f', '%.9e', '%16.9e',...            | '%.9e'                            |
    +------------+--------------------------------------------------------------------------+--------------------+---------------------------------------+-----------------------------------+
    |zoneNames   | list of zone names (first struct, the unstruct zones)                    | All                | ['Zone1','Zone2',...]                 | []                                |
    +------------+--------------------------------------------------------------------------+--------------------+---------------------------------------+-----------------------------------+
    |links       | list of list of 4 strings (see after)                                    | HDF                | [['.', 'cart.hdf', '/Base', '/Base']] | None                              |
    +------------+--------------------------------------------------------------------------+--------------------+---------------------------------------+-----------------------------------+

    Links option:
    
    For hdf format only, when writing, link node path can be specified. These nodes are then not written with data but are written as links to a pointed file.
    A link is a list ['directoryOfPointedFile', 'pointedFile', 'targetNodePath', 'currentNodePath']. The 'directoryOfPointedFile'
    is the directory where the pointed file must be found, 'pointedFile' is the pointed file name, 'targetNodePath' is the path
    of pointed node (in pointed file), 'currentNodePath' is the path of node in current file. This function doesn't write 
    the pointed file. You must explicitely write it with another call to convertPyTree2File.

    *Example of use:*

    * `HDF file write with links <Examples/Converter/links2PT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/links2PT.py

---------------------------------------------------------------------------


Preconditionning (hook)
------------------------

Preconditionning is used to create pre-computed opaque search structures 
on zones. These opaque search structures are called hook and are used in
the Geometrical identification functions and other functions of Connector
and Post.

-------------------------------------------------------------------------------

.. py:function:: Converter.createHook(a, functionName)

    Create a hook for use with identification function 'functionName'.
    For "extractMesh" and "adt", input is intended to be a set of zones, otherwise
    a hook is intended to be created on a single zone.

    :param a: input data
    :type a: [array] or [zone]
    :param functionName: function the hook is made for (see functionName_)
    :type functionName: string
    :return: hook
    :rtype: opaque structure

    *Example of use:*

    * `Create hook (array) <Examples/Converter/createHook.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createHook.py

    * `Create hook (pyTree) <Examples/Converter/createHookPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createHookPT.py


.. _functionName:

    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+
    | Function name     | Type of storage                                                                | Usage                                                  |
    +===================+================================================================================+========================================================+
    |'extractMesh'      | Bounding boxes of cells stored in an ADT. Valid for structured and TETRA zones.| Post.extractMesh, Post.extractPoint                    |
    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+
    |'adt'              | Bounding boxes of cells stored in an ADT. Valid for structured and TETRA zones.| Connector.setInterpData, Connector.setIBCData          |
    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+
    |'nodes'            |  Mesh nodes stored in a k-d tree                                               | Converter.identifyNodes, Converter.nearestNodes        |
    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+
    |'faceCenters'      |  Mesh face centers stored in a k-d tree                                        | Converter.identifyFaces, Converter.nearestFaces        |
    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+
    |'elementCenters'   |  Mesh element centers stored in a k-d tree                                     | Converter.identifyElements, Converter.nearestElements  |
    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+


-------------------------------------------------------------------------------

.. py:function:: Converter.createGlobalHook(a, functionName, indir=0)

    Create a global hook (one single search structure) for a set of zones and
    for use with identification function 'functionName' or identifySolutions.
    If indir=1, the function also returns an indirection specifying the
    zone number of each index.

    :param a: input data
    :type a: [arrays] or [zones]
    :param functionName: function the hook is made for (see functionName2_)
    :type functionName: string
    :return: hook
    :rtype: opaque structure

    *Example of use:*

    * `Create global hook (array) <Examples/Converter/createGlobalHook.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createGlobalHook.py

    * `Create global hook (pyTree) <Examples/Converter/createGlobalHookPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createGlobalHookPT.py

.. _functionName2:

    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+
    | Function name     | Type of storage                                                                | Usage                                                  |
    +===================+================================================================================+========================================================+
    |'nodes'            |  Mesh nodes stored in a k-d tree                                               | Converter.identifyNodes, Converter.nearestNodes        |
    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+
    |'faceCenters'      |  Mesh face centers stored in a k-d tree                                        | Converter.identifyFaces, Converter.nearestFaces        |
    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+
    |'elementCenters'   |  Mesh element centers stored in a k-d tree                                     | Converter.identifyElements, Converter.nearestElements  |
    +-------------------+--------------------------------------------------------------------------------+--------------------------------------------------------+


-------------------------------------------------------------------------------

.. py:function:: Converter.freeHook(hook)

    Free a hook created with createHook.

    :param hook: hook
    :type hook: opaque search structure as created by createHook
    
    *Example of use:*

    * `Free hook (array) <Examples/Converter/freeHook.py>`_:

    .. literalinclude:: ../build/Examples/Converter/freeHook.py

    * `Free hook (pyTree) <Examples/Converter/freeHookPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/freeHookPT.py

---------------------------------------------------------------------------


Geometrical identification
----------------------------

.. py:function:: Converter.identifyNodes(hook, a, tol=1.e-11)

    Identify nodes of a with nodes stored in hook. Return the indices of hook 
    corresponding to the nodes of a. If a point is not identified,
    its returned index is -1.

    :param hook: hook
    :type hook: created by createHook
    :param a: input data
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :param tol: matching tolerance
    :type tol: float
    :return: indices of identified points
    :rtype: numpy array or list of numpy arrays

    *Example of use:*

    * `Identify nodes in a hook (array) <Examples/Converter/identifyNodes.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifyNodes.py

    * `Identify nodes in a hook (pyTree) <Examples/Converter/identifyNodesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifyNodesPT.py

    * `Identify nodes in multiple zones (pyTree) <Examples/Converter/identifyNodesMBPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifyNodesMBPT.py

------------------------------------------------------------------------------------------

.. py:function:: Converter.identifyFaces(hook, a, tol=1.e-11)

    Identify face centers of a with points stored in hook. Return the indices of hook 
    corresponding to the faces of a. If a face is not identified,
    its returned index is -1.

    :param hook: hook
    :type hook: created by createHook
    :param a: input data
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :param tol: matching tolerance
    :type tol: float
    :return: indices of identified faces
    :rtype: numpy array or list of numpy arrays

    *Example of use:*

    * `Identify faces in a hook (array) <Examples/Converter/identifyFaces.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifyFaces.py

    * `Identify faces in a hook (pyTree) <Examples/Converter/identifyFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifyFacesPT.py

------------------------------------------------------------------------------------------

.. py:function:: Converter.identifyElements(hook, a, tol=1.e-11)

    Identify element centers of a with points stored in hook. Return the indices of hook 
    corresponding to the elements of a. If a elements is not identified,
    its returned index is -1.

    :param hook: hook
    :type hook: created by createHook
    :param a: input data
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :param tol: matching tolerance
    :type tol: float
    :return: indices of identified elements
    :rtype: numpy array or list of numpy arrays

    *Example of use:*

    * `Identify elements in a hook (array) <Examples/Converter/identifyElements.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifyElements.py

    * `Identify elements in a hook (pyTree) <Examples/Converter/identifyElementsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifyElementsPT.py

------------------------------------------------------------------------------------------

.. py:function:: Converter.identifySolutions(tRcv, tDnr, hookN=None, hookC=None, vars=[], tol=1.e6)

    Set the solution field in tRcv with the nearest point solution of tDnr.
    Hooks must be global hooks on tDnr.

    Exists also as an in-place version (_identifySolutions) which modifies tRcv and returns None.

    :param tRcv: receiver data
    :type tRcv: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :param tDnr: donor data
    :type tDnr: [array,list of arrays] or [pyTree, base, zone, list of zones]
    
    :param hookN: global hook if field on nodes
    :type hookN: created by createGlobalHook
    :param hookC: global hook if field on centers
    :type hookC: created by createGlobalHook
    :param vars: variable names list
    :type vars: list of strings
    :param tol: tolerance for matching
    :type tol: float
    :return: a reference copy of tRcv
    :rtype: identical to input

    *Example of use:*

    * `Identify solutions (array) <Examples/Converter/identifySolutions.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifySolutions.py

    * `Identify solutions (pyTree) <Examples/Converter/identifySolutionsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/identifySolutionsPT.py

------------------------------------------------------------------------------------------


.. py:function:: Converter.nearestNodes(hook, a)

    Find nearest points stored in hook to the nodes of a. 
    Return the indices of hook nearest of a given node of a and
    the corresponing distance.

    :param hook: hook
    :type hook: created by createHook
    :param a: input data
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :return: indices and distance of nearest points
    :rtype: tuple of 2 numpys or list of tuple of 2 numpys

    *Example of use:*

    * `Find nearest nodes in a hook (array) <Examples/Converter/nearestNodes.py>`_:

    .. literalinclude:: ../build/Examples/Converter/nearestNodes.py

    * `Find nearest nodes in a hook (pyTree) <Examples/Converter/nearestNodesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/nearestNodesPT.py

------------------------------------------------------------------------------------------

.. py:function:: Converter.nearestFaces(hook, a)

    Find nearest points stored in hook to the face centers of a. 
    Return the indices of hook nearest of a given face of a and
    the corresponing distance.

    :param hook: hook
    :type hook: created by createHook
    :param a: input data
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :return: indices and distance of nearest points
    :rtype: tuple of 2 numpys or list of tuple of 2 numpys

    *Example of use:*

    * `Find nearest faces in a hook (array) <Examples/Converter/nearestFaces.py>`_:

    .. literalinclude:: ../build/Examples/Converter/nearestFaces.py

    * `Find nearest faces in a hook (pyTree) <Examples/Converter/nearestFacesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/nearestFacesPT.py

------------------------------------------------------------------------------------------

.. py:function:: Converter.nearestElements(hook, a)

    Find nearest points stored in hook to the elements centers of a. 
    Return the indices of hook nearest of a given element of a and
    the corresponing distance.

    :param hook: hook
    :type hook: created by createHook
    :param a: input data
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :return: indices and distance of nearest points
    :rtype: tuple of 2 numpys or list of tuple of 2 numpys

    *Example of use:*

    * `Find nearest elements in a hook (array) <Examples/Converter/nearestElements.py>`_:

    .. literalinclude:: ../build/Examples/Converter/nearestElements.py

    * `Find nearest elements in a hook (pyTree) <Examples/Converter/nearestElementsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/nearestElementsPT.py


------------------------------------------------------------------------------------------

.. py:function:: Converter.createGlobalIndex(a)

    Create a index field corresponding to the vertex number.

    :param a: input data
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :return: input data with a 'globalIndex' field
    :rtype: Identical to input

    *Example of use:*

    * `Create global index (array) <Examples/Converter/createGlobalIndex.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createGlobalIndex.py

    * `Create global index (pyTree) <Examples/Converter/createGlobalIndexPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createGlobalIndexPT.py

------------------------------------------------------------------------------------------

.. py:function:: Converter.recoverGlobalIndex(a, b)

    Push the field of b in a follwing the global index field.

    :param a: input data
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :param b: input data with 'globalIndex' field
    :type b: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :return: modified a with the field of b
    :rtype: Identical to a

    *Example of use:*

    * `Recover global index (array) <Examples/Converter/recoverGlobalIndex.py>`_:

    .. literalinclude:: ../build/Examples/Converter/recoverGlobalIndex.py

    * `Recover global index (pyTree) <Examples/Converter/recoverGlobalIndexPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/recoverGlobalIndexPT.py

---------------------------------------------------------------------------

Client/server to exchange arrays/pyTrees
-------------------------------------------

.. py:function:: Converter.createSockets(nprocs=1, port=15555)

    Create sockets for receiving arrays/pyTrees. If you are sending from a 
    MPI run with nprocs, set nprocs accordingly.

    :param nprocs: the number of process the sender job is running with.
    :type nprocs: int
    :param port: the communication port
    :type port: int
    :return: socket
    :rtype: socket

-------------------------------------------------------------------------------

.. py:function:: Converter.listen(sockets)

    Listen for client sends.

    :param sockets: sockets (created with createSockets)
    :type sockets: sockets
    :return: arrays or pyTrees

    *Example of use:*

    * `Listen for arrays from server (array) <Examples/Converter/listen.py>`_:

    .. literalinclude:: ../build/Examples/Converter/listen.py

    * `Listen for pyTrees from server (pyTree) <Examples/Converter/listenPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/listenPT.py


-------------------------------------------------------------------------------

.. py:function:: Converter.send(a, host='localhost', rank=0, port=15555)
    
    Send data to the server.

    :param a: input data 
    :type a: [array,list of arrays] or [pyTree, base, zone, list of zones]
    :param host: host we are sending to
    :type host: string
    :param rank: rank of sending process
    :type rank: int
    :param port: communication port (must be the same as createSockets)
    :type port: int

    *Example of use:*

    * `Send arrays to server (array) <Examples/Converter/send.py>`_:

    .. literalinclude:: ../build/Examples/Converter/send.py

    * `Send pyTrees to server (pyTree) <Examples/Converter/sendPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/sendPT.py


---------------------------------------------------------------------------


Converter arrays / 3D arrays conversion
-------------------------------------------

In some applications, arrays must be seen as 3D arrays, that is (ni,nj,nk) numpy 
arrays instead of (nfld, ni*nj*nk) arrays. A 3D array is defined 
as [ ['x','y',...],[ ax, ay, ... ] ] where ax is a (ni,nj,nk) numpy array 
corresponding to variable x, and so on...

--------------------------------------------------------------------

.. py:function:: Converter.Array3D.convertArrays2Arrays3D(a)

    Convert arrays to 3D arrays (ni,nj,nk).

    :param a: input data 
    :type a: [list of arrays]
    :return: list of 3D arrays

    *Example of use:*

    * `Create 3D arrays from Converter arrays (array) <Examples/Converter/convertArray2Array3D.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArray2Array3D.py

-------------------------------------------------------------------------------

.. py:function:: Converter.Array3D.convertArrays3D2Arrays(a)

    Convert 3D arrays to Converter arrays.

    :param a: input data 
    :type a: [list of 3D arrays]
    :return: list of arrays

    *Example of use:*

    * `Create Converter arrays from 3D arrays (array) <Examples/Converter/convertArray3D2Array.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertArray3D2Array.py


---------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2


Index
=======

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

