.. elsAProfile documentation master file


Converter.elsAProfile: specific elsA CGNS functions
=====================================================


Preamble
########

This module provides functions to adapt a standard CGNS/python tree
for use with ONERA aerodynamic solver *elsA*. 

To use the module::

    import Converter.elsAProfile as elsAProfile


.. py:module:: Converter.elsAProfile


List of functions
##################

**-- Conversion to elsA CGNS**

.. autosummary::

   Converter.elsAProfile.adaptPeriodicMatch
   Converter.elsAProfile.adaptNearMatch
   Converter.elsAProfile.rmGCOverlap
   Converter.elsAProfile.overlapGC2BC
   Converter.elsAProfile.fillNeighbourList
   Converter.elsAProfile.prefixDnrInSubRegions

**-- Addition of elsA specific nodes**

.. autosummary:: 

    Converter.elsAProfile.addPeriodicDataInSolverParam
    Converter.elsAProfile.addOutput
    Converter.elsAProfile.addOutputForces
    Converter.elsAProfile.addOutputFriction   
    Converter.elsAProfile.addGlobalConvergenceHistory
    Converter.elsAProfile.addReferenceState
    Converter.elsAProfile.addFlowSolution
    Converter.elsAProfile.addFlowSolutionEoR
    Converter.elsAProfile.addTurbulentDistanceIndex
    Converter.elsAProfile.createElsaHybrid


**-- Miscellaneous**

.. autosummary::
    
    Converter.elsAProfile.getCGNSkeys
    Converter.elsAProfile.buildMaskFiles
    Converter.elsAProfile.convert2elsAxdt

.. **-- Macro-functions**

.. .. autosummary::

..     Converter.elsAProfile.convertChimera2elsA
..     Converter.elsAProfile.cleanTree



Contents
###########

Conversion
-----------


.. py:function:: Converter.elsAProfile.adaptPeriodicMatch(t, clean=False)

    Convert periodic grid connectivities to periodic information for elsA. 
    A '.Solver#Property' node is added in the grid connectivity node with
    children nodes of name: 'type', 'jtype', 'jtopo', 'ptype'... 
    A '.Solver#Param' node is added in the zone, providing periodicity information. 
    Children node names are: 'axis_ang_1', 'axis_ang_2',
    'axis_pnt_x', 'axis_pnt_y', 'axis_pnt_z', 'axis_vct_x', 'axis_vct_y', 'axis_vct_z', 'periodic_dir'.
    If 'Periodic_t' nodes are not removed, if RotationAngle exists, it is set in radians if it was
    defined in degrees (regarding the AngleUnits subchild node).
    Exists also as an in-place function (_adaptPeriodicMatch) that modifies t and
    returns None.

    :param t: input data
    :type  t:  pyTree, list of bases, base, list of zones, zone 
    :param clean: if True, removes useless nodes (Periodic_t nodes) for elsA
    :type clean: Boolean
    :return: modified reference copy of t
    :rtype: Identical to input

    *Example of use:*

    * `Adapt periodic match condition for elsA solver (pyTree) <Examples/Converter/adaptPeriodicMatchPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/adaptPeriodicMatchPT.py


---------------------------------------------------------------------------


.. py:function:: Converter.elsAProfile.adaptNearMatch(t)

    Convert the nearmatch grid connectivities compliant with Cassiopee for
    use with elsA.
    In particular, a ".Solver#Property" node is added with all the information required for elsA.
    Exists also as an in-place function (_adaptNearMatch) that modifies t
    and returns None.

    :param t: input data
    :type  t:  pyTree, list of bases, base, list of zones, zone 
    :return: modified reference copy of t
    :rtype: Identical to input

    *Example of use:*

    * `Adapt nearmatch condition for elsA solver (pyTree) <Examples/Converter/adaptNearmatchPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/adaptNearmatchPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.rmGCOverlap(t)

    Remove the overlap grid connectivities described as a GridConnectivity type node 
    of value 'Overset' in the t.
    Exists also as an in-place function (_rmGCOverlap) that modifies t and returns None.

    :param t: input data
    :type  t: pyTree, list of bases, base, list of zones, zone 
    :return: modified reference copy of t 
    :rtype: Identical to input

    *Example of use:*

    * `Remove the Overlap condition described as Grid Connectivity (pyTree) <Examples/Converter/rmGCOverlapPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmGCOverlapPT.py


---------------------------------------------------------------------------


.. py:function:: Converter.elsAProfile.overlapGC2BC(t)
    
    Convert the 'Overset' grid connectivity nodes compliant with Cassiopee to an 'Overlap' BC compliant with elsA.
    The created BC node is attached to a family of BC of name 'Fam_Ovlp'+baseName, where baseName is the name of the parent base. Prefix 'Fam_Ovlp_' of the family BC of 'Overlap' BCs can be redefined
    using Converter.elsAProfile.__FAMOVERLAPBC__ container (e.g. Converter.elsAProfile.__FAMOVERLAPBC__="OVERLAP").
    In case of a doubly defined overlap grid connectivity, the default name of the family of bcs is prefixed by 'Fam_OvlpDD_'.
    In that case, there is one family of classical overlap BCs per receptor base and one family of doubly defined overlap BCs too (if classical and doubly defined overlap bcs exist 
    in that base).
    A '.Solver#Overlap' node is created in each family of BCs and contains the node 'NeighbourList' (where donor zones can be specified). The value is empty here.
    In case of a doubly defined family, a 'doubly_defined' node is added (with value 'active'). 

    This function exists also as an in-place function (_overlapBC2GC) that modifies t 
    and returns None.

    :param t: input tree
    :type t: pyTree
    :return: modified reference copy of t

    *Example of use:*

    * `Convert overlap GCs to overlap BCs with families of BCs attached to each receptor base (pyTree) <Examples/Converter/overlapGC2BCPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/overlapGC2BCPT.py


---------------------------------------------------------------------------


.. py:function:: Converter.elsAProfile.fillNeighbourList(t, sameBase=0)
    
    Fill the NeighbourList node with families of donor zones. Donor zones are obtained by intersection of the receptor base.
    For doubly defined overlap BCs, the NeighbourList node is made of the family of specified donor zones. 
    If sameBase=1, allows for donor zones to be in the same base.

    This function exists also as an in-place function (_fillNeighbourList) that modifies t and
    returns None.

    :param t: input tree
    :type t: pyTree
    :param sameBase: if sameBase=1, donor zones are allowed on the same base as the receptor zone. 
    :type sameBase: integer (0 or 1)
    :return: modified reference copy of t

    *Example of use:*

    * `Fill NeighbourList nodes (pyTree) <Examples/Converter/fillNeighbourListPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/fillNeighbourListPT.py



---------------------------------------------------------------------------



.. py:function:: Converter.elsAProfile.prefixDnrInSubRegions(t)
    
    When Chimera connectivity (computed by Connector.PyTree.setInterpolations) is stored in the pyTree
    as subregion nodes of name 'ID*', corresponding donor zone names are defined by the value of the node 'ID_*'.
    To be read by elsA, the donor zone names are prefixed by their parent base name.

    This function exists also as an in-place function (_prefixDnrInSubRegions) that modifies t and
    returns None.

    :param t: input tree
    :type t: pyTree
    :return: modified reference copy of t

    *Example of use:*

    * `Prefix donor names by their base (pyTree) <Examples/Converter/prefixDnrInSubRegionsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/prefixDnrInSubRegionsPT.py







---------------------------------------------------------------------------

Addition of elsA nodes
----------------------



.. py:function:: Converter.elsAProfile.addPeriodicDataInSolverParam(t, rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.], NAzimutalSectors=0, isChimera=False)

    Add information about periodicity by rotation in '.Solver#Param' node of zones of t.
    Exists also as an in-place function (_addPeriodicDataInSolverParam) that 
    modifies t ans returns None.

    :param t: input data
    :type t: pyTree, list of bases, base, list of zones, zone 
    :param rotationCenter: coordinates of rotation center
    :type rotationCenter: list of 3 floats
    :param rotationAngle: rotation angle along each axis
    :type rotationAngle: list of 3 floats
    :param NAzimutalSectors: number of azimuthal sectors to define 360 degrees 
    :type NAzimutalSectors: integer
    :param isChimera: to say that a zone is involved in a periodic Chimera configuration 
    :type isChimera: Boolean
    :return: modified reference copy of t
    :rtype: same as input

    *Example of use:*

    * 'Add periodic data (rotation) in '.Solver#Param' nodes for elsA <Examples/Converter/addPeriodicDataInSolverParamPT.py>'_:

    .. literalinclude:: ../build/Examples/Converter/addPeriodicDataInSolverParamPT.py


---------------------------------------------------------------------------


.. py:function:: Converter.elsAProfile.addOutput(a, Dict, name='', update=False)

    Add a '.Solver#Output' node suffixed by name in order to perform some extractions defined in a dictionary Dict. 
    Dict is a python dictionary containing all output information. For example:
    >>> Dict={}

    >>> Dict["var"]="convflux_rou convflux_rov convflux_row"
    
    >>> Dict["loc"]=1
    
    >>> Dict["fluxcoef"]=1.
    
    >>> Dict["period"]=10

    Name is optional to suffix the '.Solver#Output' name.
    If the node 'Solver#Output' already exists, it can be cleared and recreated (update=True). 
    In the other case, if the value of the dictionary key is already set in the '.Solver#Output', it is updated.
    Exists also as in place version (_addOutput) that modifies the node a and returns None.

    :param a: input node defining topologically the extraction (e.g. a bc node) 
    :type  a:  standard node (zone, bc,...) 
    :param Dict: dictionary containing all the output information. 
    :type  Dict: python dictionary
    :param name: to suffix the node '.Solver#Output' by name
    :type  name: string
    :param update: if True, removes the existing node of name '.Solver#Output'+name in a
    :type update: Boolean
    :return: modified reference copy of a
    :rtype: identical to input

    *Example of use:*

    * `Add a '.Solver#Output' node for elsA extraction (pyTree) <Examples/Converter/addOutputPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addOutputPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.addOutputForces(a, name="", var=None, loc=4, writingmode=1, period=None, pinf=None, fluxcoef=None, torquecoef=None, xyztorque=None, frame=None, governingEquations="NSTurbulent", xtorque=None, ytorque=None, ztorque=None)

    Add a .Solver#Output:Forces node in order to extract required value. name is a suffix that can be appened to the '.Solver#Output:Forces' name. loc is the location (integer value for elsA) for extraction. writingmode value is the integer corresponding to elsA. period is the extraction frequency

    Exists also as an in-place function (_addOutputForces) that modifies node a and returns None.

    :param a: input node  
    :type  a:  pyTree node
    :param name: suffix node name to add to '.Solver#Output:Forces'
    :type  name: string
    :param var: list of variables names
    :type  var: list of strings
    :param loc: suffix node name to add to '.Solver#Output:Forces'
    :type  loc: string
    :param writingmode: writingmode value for elsA
    :type  writingmode: integer
    :param period: period of extraction
    :type  period: integer
    :param pinf: value of farfield pressure
    :type  pinf: float
    :param fluxcoef: coefficient to correct the forces
    :type  fluxcoef: float
    :param torquecoef: coefficient to correct the torques
    :type  torquecoef: float
    :param xyztorque: coordinates of the torque origin
    :type  xyztorque: list of float
    :param frame: writing frame ('relative','absolute')
    :type  frame: string
    :param governingEquations: governing equations of the model
    :type  governingEquations: string
    :param xtorque: X-coordinate of the torque origin
    :type  xtorque: float
    :param ytorque: Y-coordinate of the torque origin
    :type  ytorque: float
    :param ztorque: Z-coordinate of the torque origin
    :type  ztorque: float
    :return: modified reference copy of a
    :rtype: same as input node type

    *Example of use:*

    * `Add a .Solver#Output:Forces node for elsA extraction (pyTree) <Examples/Converter/addOutputForcesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addOutputForcesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.addOutputFriction(a, name="", var=None, loc=4, writingmode=1, period=None, fluxcoef=None, torquecoef=None, writingframe=None)

    Add a .Solver#Output:Friction node in order to extract required value. name is a suffix that can be appened to the '.Solver#Output:Friction' name. loc is the location (integer value for elsA) for extraction. writingmode value is the integer corresponding to elsA. period is the extraction frequency.
    Exists also as in place version (_addOutputFriction) that modifies t and returns None.

    :param a: input node  
    :type  a: pyTree node 
    :param name: suffix node name to add to '.Solver#Output:Forces'
    :type  name: string
    :param var: list of variables names
    :type  var: list of strings
    :param loc: suffix node name to add to '.Solver#Output:Forces'
    :type  loc: string
    :param writingmode: writingmode value for elsA
    :type  writingmode: integer
    :param period: period of extraction
    :type  period: integer
    :param fluxcoef: coefficient to correct the forces
    :type  fluxcoef: float
    :param torquecoef: coefficient to correct the torques
    :type  torquecoef: float
    :param writingframe: writing frame ('relative','absolute')
    :type  writingframe: string
    :return: modified reference copy of a
    :rtype: identical to a

    *Example of use:*

    * `Add a .Solver#Output:Friction node for elsA solver extraction (pyTree) <Examples/Converter/addOutputFrictionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addOutputFrictionPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.addGlobalConvergenceHistory(t, normValue=0) 

    Add a convergence node in each base of the CGNS/Python tree. The type of norm used in residual computation can be specified (0: L0, 1: L2).
    Exists also as an in-place function (_addGlobalConvergenceHistory) modifying t and returning None.

    :param t: input tree  
    :type  t: pyTree 
    :param normValue: an optional integer, specifying the type of norm as value of the GlobalConvergenceHistory node.
    :type  normValue: integer (0 or 1)
    :return: modified reference copy of t
    :rtype: pyTree

    *Example of use:*

    * `Add a global convergence history node for elsA solver (pyTree) <Examples/Converter/addGlobalConvergenceHistoryPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addGlobalConvergenceHistoryPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.addReferenceState(t, conservative=None, temp=None, turbmod='spalart', name='ReferenceState', comments=None)

    Add a ReferenceState node in each base. Turbulence model, reference variables constants are mandatory.
    Reference temperature can also be defined. 
    Exists also as in place version (_addReferenceState) that modifies t and returns None.
    Depending on the turbulence model, one, two or 7 additional constants must be defined in conservative list, that define turbulence model variables (according to the CGNS standard names):
    
    - For turbmod='spalart':'TurbulentSANuTildeDensity' constant.

    - For turbmod='komega': ['TurbulentEnergyKineticDensity','TurbulentDissipationRateDensity'] constants

    - For turbmod='keps' or 'chien' or 'asm': ['TurbulentEnergyKineticDensity','TurbulentDissipationDensity']

    - For turbmod='smith': ['TurbulentEnergyKineticDensity','TurbulentLengthScaleDensity']

    - For turbmod='kkl' or 'earsm': ['TurbulentEnergyKineticDensity','TurbulentKineticPLSDensity']

    - For turbmod='rsm': ["ReynoldsStressXX","ReynoldsStressXY","ReynoldsStressXZ","ReynoldsStressYY","ReynoldsStressYZ","ReynoldsStressZZ","ReynoldsStressDissipationScale"]

    :param t: input tree  
    :type  t:  pyTree 
    :param conservative: list of constant values of reference variables [conservative + turbulence model variables].
    :type  conservative:  list of floats 
    :param temp: reference temperature (optional).
    :type  temp:  float 
    :param turbmod: turbulence model name.
    :type  turbmod:  string (possible values: 'spalart','komega','keps','chien','asm','smith','kkl','earsm','rsm')
    :param name: name of the ReferenceState node.
    :type  name:  string
    :param comments: optional comments to describe the reference state (added to ReferenceStateDescription node).
    :type  comments:  string
    :return: modified reference copy of t
    :rtype: pyTree

    *Example of use:*

    * `Add a reference state for elsA solver (pyTree) <Examples/Converter/addReferenceStatePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addReferenceStatePT.py


---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.addFlowSolution(t, name='', loc='CellCenter', variables=None, governingEquations=None, writingMode=None, writingFrame='relative', period=None, output=None, addBCExtract=False, protocol="end") 

    Add a FlowSolution node for each zone. name is a suffix that can be appened to the 'FlowSolution' name. loc must be 'Vertex' or 'CellCenter' or 'cellfict'. governingEquation is the value of the 'GoverningEquation' node. output can be optionaly a dictionary specifying the '.Solver#Output' node data. If addBCExtract is true, the boundary windows are also extracted. protocol is an optional string in 'iteration', 'end', 'after', specifying when extraction is performed.
    Note that if governingEquations is set to None, then the governing equations set in the GoverningEquations_t node in the pyTree are chosen.
    If variables=[], it does not create any variable in the FlowSolution node.
    Exists also as in place version (_addFlowSolution) that modifies t and returns None.


    :param t: input tree  
    :type  t:  pyTree 
    :param name: suffix to append to the name 'FlowSolution'
    :type  name: string
    :param loc: location of extraction ('CellCenter', 'Vertex' or 'cellfict')
    :type  loc: string
    :param variables: list of variables (elsA names)
    :type  variables: list of string
    :param governingEquations: kind of governing equations (None or standard governing equations name ('Euler', 'NSTurbulent'...)). Optional: if already defined in the pyTree in the GoverningEquations_t node, its value is chosen.
    :type  governingEquations: string
    :param writingMode: writingmode value for elsA
    :type  writingMode: integer
    :param writingFrame: frame of FlowSolution ('absolute','relative')
    :type  writingFrame: string
    :param period: period of extraction
    :type  period: integer
    :param output: an optional dictionary of node to add to '.Solver#Output' if necessary
    :type  output: dictionary
    :param addBCExtract: True to add extractions from BC windows (pseudo CellFict)
    :type  addBCExtract: Boolean
    :param protocol: Protocol node value ('iteration','end','after') for extractions
    :type  protocol: string
    :return: modified reference copy of t
    :rtype: pyTree: pyTree

    *Example of use:*

    * `Add the FlowSolution nodes for elsA solver (pyTree) <Examples/Converter/addFlowSolutionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addFlowSolutionPT.py

---------------------------------------------------------------------------


.. py:function:: Converter.elsAProfile.addFlowSolutionEoR(t, name='', variables=None, governingEquations=None, writingFrame='relative', addBCExtract=False, protocol="end")

    Add a FlowSolution#EndOfRun node located at cell centers for each zone. Parameter name is a suffix that can be appended to the 'FlowSolution#EndOfRun' name. 
    Parameter governingEquations is the value of the 'GoverningEquations' node. If set to None, then its value is obtained from the 'GoverningEquations_t' node in the pyTree. 
    If variables=[], it does not create any variable in the FlowSolution node.
    If addBCExtract is True, the boundary windows are also extracted. Parameter protocol is an optional string in 'iteration', 'end', 'after', specifying when extraction is performed.
    Exists also as an in-place function (_addFlowSolutionEoR) that modifies t and returns None.

    :param t: input tree
    :type  t:  pyTree
    :param name: suffix to append to 'FlowSolution#EndOfRun'
    :type  name: string
    :param variables: list of variables (elsA names)
    :type  variables: list of strings
    :param governingEquations: kind of governing equations ('Euler','NSLaminar','NSTurbulent')
    :type  governingEquations: string
    :param writingFrame: frame of FlowSolution ('absolute','relative')
    :type  writingFrame: string
    :param addBCExtract: if True, extractions are also achieved from BC windows (pseudo CellFict)
    :type  addBCExtract: Boolean
    :param protocol: Protocol node value ('iteration','end','after') for extractions
    :type  protocol: string
    :return: modified reference copy of t
    :rtype: pyTree: pyTree

    *Example of use:*

    * `Add the FlowSolution#EndOfRun nodes for elsA solver (pyTree) <Examples/Converter/addFlowSolutionEoRPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addFlowSolutionEoRPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.addNeighbours__(t, sameBase=0)

    Fill the NeighbourList nodes if needed with bounding-box domains intersection between bases. If sameBase=1, the intersecting domains are also searched in the current base.

    :param t: input tree  
    :type  t: pyTree
    :param sameBase: choice for keeping the same current base in the bounding-box intersection (1 for yes, 0 else) 
    :type  sameBase:  integer
    :return: modified reference copy of t
    :rtype: pyTree

    *Example of use:*

    * `Fill the NeighbourList nodes with bounding-box domains intersection (pyTree) <Examples/Converter/addNeighboursPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addNeighboursPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.addTurbulentDistanceIndex(t)

    Add a node 'TurbulentDistanceIndex' initialized with -1 (float) in the container of flow solution at centers if TurbulentDistance node exists in the same container. 
    To  set that node in the 'FlowSolution#Init' container, be sure that the 'TurbulentDistance' node is in that container too and set: Internal.__FlowSolutionCenters__='FlowSolution#Init'.
    Exists also as in place version (_addTurbulentDistanceIndex) that modifies t and returns None.

    :param t: input data
    :type  t: pyTree, base, zone, list of zones
    :return: modified reference copy of t
    :rtype: Identical to t

    *Example of use:*

    * `Add the TurbulentDistance index node for elsA solver (pyTree) <Examples/Converter/addTurbulentDistanceIndexPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addTurbulentDistanceIndexPT.py



---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.createElsaHybrid(t, method=0, axe2D=0, methodPE=0) 

    Add nodes required by elsA Hybrid solver for NGON zones.
    For elsA < 3.8.01, use method=0, for newer versions use method=1.
    If the mesh is 2D, use axe2D to precise 2D plane (0: (x,y), 1: (x,z), 2: (y,z)).
    If the mesh has poor quality cells (concave...) use methodPE=1 to build the ParentElement node in a topological manner.
    Exists also as in place version (_createElsaHybrid) that modifies t and returns None.

    :param t: input data
    :type  t: pyTree, base, zone, list of zones
    :param method: 0 (if for elsA < 3.8.01), 1 otherwise
    :type  method: 0 or 1
    :param axe2D: 1 if (y,z), 2 if (x,z), 3 if (x,y)
    :type axe2D: int
    :param methodPE: 0 (for regular mesh), 1 otherwise
    :type  methodPE: 0 or 1
    :return: modified reference copy of t
    :rtype: identical to t

    *Example of use:*

    * `Create specific nodes for elsA Hybrid (pyTree) <Examples/Converter/createElsaHybridPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createElsaHybridPT.py



Miscellaneous
-----------------------------


.. py:function:: Converter.elsAProfile.getCGNSkeys(key, verbose=True)
    
    Return the CGNS name corresponding to an elsA key.

    :param key: elsA key 
    :type key: string
    :return: CGNS name of elsA key
    :rtype: string

    *Example of use:*

    * `Return the CGNS name corresponding to an elsA key <Examples/Converter/getCGNSkeys.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getCGNSkeys.py

---------------------------------------------------------------------------


.. py:function:: Converter.elsAProfile.buildMaskFiles(t, keepOversetHoles=True, fileDir='.', prefixBase=False)

    Write the mask files in bin_v3d format. It can also keep or delete the OversetHoles nodes in the tree. The fileDir variable allows to choose the directory where the v3d hole files will be written.
    OversetHoles nodes can be created first by Connector.PyTree.cellN2OversetHoles.
    Exists also as an in-place function (_buildMaskFiles) that modifies t and returns None.


    :param t: input tree  
    :type  t:  pyTree 
    :param keepOversetHoles: choice for keeping or not the OversetHoles in the tree 
    :type  keepOversetHoles: boolean 
    :param fileDir: path of directory for file writing 
    :type  fileDir: string 
    :param prefixBase: if True, add the base name to the zone name in the created file of name 'hole_*'. Results in 'hole_mybase_myzone.v3d' or 'hole_myzone.v3d' 
    :type prefixBase: Boolean
    :return: modified reference copy of t
    :rtype: pyTree

    *Example of use:*

    * `Build the mask files for elsA solver (pyTree) <Examples/Converter/buildMaskFilesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/buildMaskFilesPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.elsAProfile.convert2elsAxdt(t)

    Macro-function that converts some Cassiopee data stored in the pyTree into elsA-compliant nodes.
    It performs the following functions available as single functions:
    
    - addTurbulentDistanceIndex(t)

    - buildMaskFiles(t)

    - adaptNearMatch(t)

    - adaptPeriodicMatch(t)

    - overlapGC2BC(t)

    - rmGCOverlap(t)

    - fillNeighbourList(t)
    
    - prefixDnrInSubRegions(t) 
   
    Exists also as in place version (_convert2elsAxdt) that modifies t and returns None.

    :param t: input tree  
    :type  t:  pyTree 
    :return: modified reference copy of t
    :rtype: pyTree: pyTree

    *Example of use:*

    * `Convert a tree for elsA solver (pyTree) <Examples/Converter/convert2elsAxdtPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convert2elsAxdtPT.py




---------------------------------------------------------------------------


.. toctree::
   :maxdepth: 2   

Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

