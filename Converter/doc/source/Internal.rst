.. Internal documentation master file


Internal: CGNS/Python tree management
=========================================


Preamble
########

This module provides simple and efficient functions for 
creating/traversing/manipulating CGNS/Python data tree (refered here as *pyTree*).

This module is part of Cassiopee, a free open-source pre- and post-processor 
for CFD simulations.

To use the module::

    import Converter.Internal as Internal

All functions of Cassiopee modules (Converter.PyTree, Transform.PyTree, ...) 
work on 3 types of containers: a container for grid coordinates named __GridCoordinates__, 
a container for node solution named __FlowSolutionNodes__ and a container for 
center solution named __FlowSolutionCenters__. By default::

    Internal.__GridCoordinates__ = 'GridCoordinates'
    Internal.__FlowSolutionNodes__ = 'FlowSolution'
    Internal.__FlowSolutionCenters__ = 'FlowSolution#Centers'


To make the functions operate on another named container for example 
'FlowSolution#Init', start your python script with::

    Internal.__FlowSolutionNodes__ = 'FlowSolution#Init' 

To automatically set the container names used by Cassiopee 
to the ones that exists in your pyTree t, use::

    Internal.autoSetContainers(t)

If multiple containers exists in your pyTree, for instance 'FlowSolution1' and
'FlowSolution2' for 'FlowSolution_t', the the first met is used by Cassiopee.


Through all this documentation, a pyTree node is a python 
list: **['NodeName', value numpy array, [node children list], 'Type_t']** as described by
CGNS/python standard (https://cgns.github.io/CGNS_docs_current/python/sidstopython.pdf).

*'NodeName'* is a string describing the node **name**; *value numpy array* is a 
numpy array corresponding to the **value** stored in node;
*'Type_t'* is a string describing the node **type**; *node children list* is a
list of pyTree nodes that are the children of this node.

.. py:module:: Converter.Internal


List of functions
#################

**-- Node tests**

.. autosummary::

   Converter.Internal.isTopTree
   Converter.Internal.isStdNode
   Converter.Internal.typeOfNode
   Converter.Internal.isType
   Converter.Internal.isName
   Converter.Internal.isValue
   Converter.Internal.isChild

**-- Set/create generic nodes**

.. autosummary::

   Converter.Internal.setName
   Converter.Internal.setType
   Converter.Internal.setValue
   Converter.Internal.createNode
   Converter.Internal.addChild
   Converter.Internal.createChild
   Converter.Internal.createUniqueChild

**-- Access nodes**

.. autosummary::

    Converter.Internal.getName
    Converter.Internal.getType
    Converter.Internal.getChildren
    Converter.Internal.getValue
    Converter.Internal.getVal
    Converter.Internal.getPath
    Converter.Internal.getNodesFromName
    Converter.Internal.getNodeFromName
    Converter.Internal.getByName
    Converter.Internal.getChildFromName
    Converter.Internal.getNodesFromType
    Converter.Internal.getNodeFromType
    Converter.Internal.getByType
    Converter.Internal.getChildFromType
    Converter.Internal.getNodesFromNameAndType
    Converter.Internal.getNodeFromNameAndType
    Converter.Internal.getNodesFromValue

    Converter.Internal.getParentOfNode
    Converter.Internal.getParentFromType
    Converter.Internal.getParentsFromType
    Converter.Internal.getNodePosition

    Converter.Internal.getNodeFromPath
    Converter.Internal.getPathsFromName
    Converter.Internal.getPathsFromType
    Converter.Internal.getPathsFromValue
    Converter.Internal.getPathLeaf
    Converter.Internal.getPathAncestor
    Converter.Internal.getZonePaths

    Converter.Internal.getZones
    Converter.Internal.getZonesPerIteration
    Converter.Internal.getBases
    Converter.Internal.getZoneDim
    Converter.Internal.getZoneType

**-- Check nodes**

.. autosummary::

    Converter.Internal.printTree
    Converter.Internal.getSizeOf
    Converter.Internal.checkPyTree
    Converter.Internal.correctPyTree

**-- Copy nodes**

.. autosummary::

    Converter.Internal.copyRef
    Converter.Internal.copyTree
    Converter.Internal.copyValue
    Converter.Internal.copyNode

**-- Add/remove/move nodes**

.. autosummary::

    Converter.Internal.append
    Converter.Internal.rmNode
    Converter.Internal.rmNodeByPath
    Converter.Internal.rmNodesByName
    Converter.Internal.rmNodesByType
    Converter.Internal.rmNodesByNameAndType
    Converter.Internal.rmNodesByValue
    Converter.Internal.moveNodeFromPaths

**-- Modify nodes**

.. autosummary::
    
    Converter.Internal.merge
    Converter.Internal.renameNode
    Converter.Internal.sortByName
    Converter.Internal.appendBaseName2ZoneName
    Converter.Internal.groupBCByBCType
    
**-- Create specific CGNS nodes**

.. autosummary::

    Converter.Internal.newCGNSTree
    Converter.Internal.newCGNSBase
    Converter.Internal.newZone
    Converter.Internal.newGridCoordinates
    Converter.Internal.newDataArray
    Converter.Internal.newDataClass
    Converter.Internal.newDimensionalUnits
    Converter.Internal.newDimensionalExponents
    Converter.Internal.newDataConversion
    Converter.Internal.newDescriptor
    Converter.Internal.newGridLocation
    Converter.Internal.newIndexArray
    Converter.Internal.newPointList
    Converter.Internal.newPointRange
    Converter.Internal.newRind
    Converter.Internal.newSimulationType
    Converter.Internal.newOrdinal
    Converter.Internal.newDiscreteData
    Converter.Internal.newIntegralData
    Converter.Internal.newElements
    Converter.Internal.newParentElements
    Converter.Internal.newParentElementsPosition
    Converter.Internal.newZoneBC
    Converter.Internal.newBC
    Converter.Internal.newBCDataSet 
    Converter.Internal.newBCData
    Converter.Internal.newBCProperty
    Converter.Internal.newAxiSymmetry
    Converter.Internal.newRotatingCoordinates
    Converter.Internal.newFlowSolution
    Converter.Internal.newZoneGridConnectivity
    Converter.Internal.newGridConnectivity1to1
    Converter.Internal.newGridConnectivity
    Converter.Internal.newGridConnectivityType
    Converter.Internal.newGridConnectivityProperty
    Converter.Internal.newPeriodic
    Converter.Internal.newZoneSubRegion 
    Converter.Internal.newOversetHoles 
    Converter.Internal.newFlowEquationSet
    Converter.Internal.newGoverningEquations
    Converter.Internal.newGasModel
    Converter.Internal.newThermalConductivityModel
    Converter.Internal.newViscosityModel
    Converter.Internal.newTurbulenceClosure
    Converter.Internal.newTurbulenceModel
    Converter.Internal.newThermalRelaxationModel
    Converter.Internal.newChemicalKineticsModel
    Converter.Internal.newEMElectricFieldModel
    Converter.Internal.newEMConductivityModel
    Converter.Internal.newBaseIterativeData
    Converter.Internal.newZoneIterativeData
    Converter.Internal.newRigidGridMotion
    Converter.Internal.newRigidGridMotionType
    Converter.Internal.newReferenceState
    Converter.Internal.newConvergenceHistory
    Converter.Internal.newFamily
    Converter.Internal.newFamilyBC
    Converter.Internal.newGeometryReference
    Converter.Internal.newArbitraryGridMotion
    Converter.Internal.newUserDefinedData
    Converter.Internal.newGravity


Contents
#########

Node tests
-----------

.. py:function:: Converter.Internal.isTopTree(node)

    Return True if node is a top tree node.

    :param node:  Input node
    :type  node:  a pyTree node
    :rtype: Boolean: True or False

    *Example of use:*

    * `Tell if a node is a top tree (pyTree) <Examples/Converter/isTopTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isTopTreePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.isStdNode(node)

    Return 0 if node is a list of standard pyTree nodes, -1 if node is a standard pyTree node,
    -2 otherwise. 

    :param node:  Input node
    :type  node:  a pyTree node
    :rtype: int: 0: node is a list of standard nodes, -1: node is a standard node, -2: otherwise

    *Example of use:*

    * `Tell if a node is a standard one (pyTree) <Examples/Converter/isStdNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isStdNodePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.typeOfNode(node)

    Return 1 if node is a zone node, 2 if node is a list of zones,
    3 if node is a tree, 4 if node is a base, 5 if node is a list
    of bases, -1 otherwise.

    :param node:  Input node
    :type  node:  a pyTree node
    :rtype: int: 1: zone node, 2: list of zones, 3: pyTree, 4: base node, 5: list of bases, -1: otherwise

    *Example of use:*

    * `Tell the type of node (pyTree) <Examples/Converter/typeOfNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/typeOfNodePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.isType(node, ntype)

    Compare given type and node type. Wildcards are accepted for ntype.

    :param node:  Input node
    :type  node:  a pyTree node
    :param ntype:  CGNS type string ('DataArray_t', 'Zone_t', ...)
    :type ntype: string
    :rtype: Boolean: True or False

    *Example of use:*

    * `Tell if a node is of given type (pyTree) <Examples/Converter/isTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isTypePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.isName(node, name)

    Compare given name and node name. Wildcards are accepted for name.

    :param node:  Input node
    :type  node:  a pyTree node
    :param name:  node name to be checked
    :type name: string
    :rtype: Boolean: True or False

    *Example of use:*

    * `Tell if a node has given name (pyTree) <Examples/Converter/isNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isNamePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.isValue(node, value)

    Compare given value and node value. If value is a string, wildcards are accepted.
    Accepts also numpy arrays as value.

    :param node:  Input node
    :type  node:  a pyTree node
    :param value:  node value to be checked
    :type value: int, float, string, numpys
    :rtype: Boolean: True or False

    *Example of use:*

    * `Tell if a node stores given value (pyTree) <Examples/Converter/isValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isValuePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.isChild(start, node)

    Return true if node is a child of start, even at deep levels.
    Exists also as isChild1 and isChild2, limited to 1 or 2 recursivity level.

    :param start:  Input node
    :type  start:  a pyTree node
    :param node:  node to be checked as a child of start
    :type node:  pyTree node
    :rtype: Boolean: True or False

    *Example of use:*

    * `Tell if a node is in the children of another one (pyTree) <Examples/Converter/isChildPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/isChildPT.py

    .. note:: new in version 2.7.
    
Set/create generic nodes
--------------------------

.. py:function:: Converter.Internal.setName(node, name)

    Set the given name to node (node is modified).

    :param node:  Input node
    :type  node:  a pyTree node
    :param name:  node name to be set
    :type name: string
    :rtype: None

    *Example of use:*

    * `Set given name to node (pyTree) <Examples/Converter/setNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/setNamePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.setType(node, ntype)

    Set the given type to node (node is modified).

    :param node:  Input node
    :type  node:  a pyTree node
    :param ntype:  node type to be set
    :type ntype: string
    :rtype: None

    *Example of use:*

    * `Set given type to node (pyTree) <Examples/Converter/setTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/setTypePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.setValue(node, value=None)

    Set the given value in node (node is modified).
    Value can be provided as an int, a float, a string or
    a numpy array but it is always stored as a numpy array in node[1].

    :param node:  Input node
    :type  node:  a pyTree node
    :param value:  node value to be set
    :type value: int, float, string, numpys
    :rtype: None

    *Example of use:*

    * `Set given value in node (pyTree) <Examples/Converter/setValueIPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/setValueIPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.createNode(name, ntype, value=None, children=None, parent=None)

    Create a node with a given name and type and optional value and children.
    If parent is present, created node is also attached to parent node.

    :param name:  created node name
    :type  name:  string
    :param ntype:  created node type ('DataArray_t', 'Zone_t', ...)
    :type ntype:   string
    :param value: created node value
    :type value: int, float, string, numpys
    :param children: a list of children nodes
    :type children: list of pyTree nodes
    :param parent: optional node. If given, created node is attached to parent.
    :type parent: pyTree node
    
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a pyTree node (pyTree) <Examples/Converter/createNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createNodePT.py
 
---------------------------------------------------------------------------

.. py:function:: Converter.Internal.addChild(node, child, pos=-1)

    Insert a child at given index in the children list of node.
    If pos=-1, add it at the end.

    :param node:  modified node
    :type  node:  pyTree node
    :param child:  node added as children of node
    :type child:   pyTree node
    :param pos: position of child in children list
    :type pos: int
    :return: child (identical to input)
    :rtype: pyTree node

    *Example of use:*

    * `Add a child node (pyTree) <Examples/Converter/addChildPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addChildPT.py


---------------------------------------------------------------------------
 
.. py:function:: Converter.Internal.createChild(node, name, ntype, value=None, children=None, pos=-1)

    Create a child node and attach it to a given node.
    Child's node name, type, value and children can be specified.
    Position in children list can also be specified. node is modified and
    newly created child node is returned.
    If a node with the same name already exists in node children list, the
    newly created node is nevertheless added. To avoid creating
    nodes with the same name, consider using createUniqueChild.

    :param node:  modified node
    :type  node:  pyTree node
    :param name:  created node name
    :type  name:  string
    :param ntype:  created node type
    :type ntype:   string
    :param value: created node value
    :type value: int, float, string, numpys
    :param children: optional list of nodes to be attached as children to the created node.
    :type children: list of pyTree nodes
    :param pos: position of child in children list
    :type pos: int
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a child node (pyTree) <Examples/Converter/createChildPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createChildPT.py
 
---------------------------------------------------------------------------
 
.. py:function:: Converter.Internal.createUniqueChild(node, name, ntype, value=None, children=None, pos=-1)

    Same as Internal.createChild except that, if a node with the same name already exists
    in the children list, its value and type are replaced with
    given values (node is modified)
    and newly created or modified child node is returned.

    :param node:  modified node
    :type  node:  pyTree node
    :param name:  created node name
    :type  name:  string
    :param ntype:  created node type
    :type ntype:   string
    :param value: created node value
    :type type: int, float, string, numpys
    :param children: optional list of nodes to be attached as children to the created node.
    :type children: list of pyTree nodes
    :param pos: position of child in children list
    :type pos: int
    :return: created (or modified) node
    :rtype: pyTree node

    *Example of use:*

    * `Create a unique child node (pyTree) <Examples/Converter/createUniqueChildPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createUniqueChildPT.py
 

Acess nodes
------------

.. py:function:: Converter.Internal.getName(node)

    Return the name of node. Completely equivalent to node[0].
   
    :param node:  input node
    :type  node:  pyTree node
    :return: name of node
    :rtype: string

    *Example of use:*

    * `Return name of pyTree node (pyTree) <Examples/Converter/getNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNamePT.py

---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getType(node)

    Return the type of node. Completely equivalent to node[3].
   
    :param node:  input node
    :type  node:  pyTree node
    :return: type of node
    :rtype: string

    *Example of use:*

    * `Return type of pyTree node (pyTree) <Examples/Converter/getTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getTypePT.py

---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getChildren(node)

    Return the type of node. Completely equivalent to node[2].
   
    :param node:  input node
    :type  node:  pyTree node
    :return: children nodes
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return children of pyTree node (pyTree) <Examples/Converter/getChildrenPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getChildrenPT.py

---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getValue(node)

    Return value of node. Depending of stored value, return string, int, float or
    numpy array. It differs from node[1], which is always a numpy array.

   
    :param node:  input node
    :type  node:  pyTree node
    :return: node value
    :rtype: string, int, float, numpy array

    *Example of use:*

    * `Return value of pyTree node (pyTree) <Examples/Converter/getValueIPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getValueIPT.py

---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getVal(node)

    Return value of node always as a numpy. Completely equivalent to node[1].
   
    :param node:  input node
    :type  node:  pyTree node
    :return: node value
    :rtype: numpy array

    *Example of use:*

    * `Return numpy value of pyTree node (pyTree) <Examples/Converter/getValPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getValPT.py

    .. note:: new in version 3.2.

---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getPath(t, node, pyCGNSLike=False)

    Return path of node in t. Path is a string describing the node names from
    t to node. For instance: Base/Zone/GridConnectivity.
    If node is not found, return None.
   
    :param t:  starting node
    :type  t:  pyTree node
    :param node:  input node
    :type  node:  pyTree node
    :param pyCGNSLike:  if True, paths don't start with CGNSTree
    :type  pyCGNSLike:  boolean
    :return: path of node in t
    :rtype: string

    *Example of use:*

    * `Return the path of node (pyTree) <Examples/Converter/getPathPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getPathPT.py


---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodesFromName(t, name)

    Return a list of nodes matching given name (wildcards allowed).
    To accelerate search, this function can be limited to 1, 2 or 3 levels from 
    the starting node using getNodesFromName1, getNodesFromName2, getNodesFromName3. 
    Then wildcards are then no longer allowed.

    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :param name:  nodes name we are looking for 
    :type  name:  string
    :return: list of nodes that matches name (shared with t)
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return list of nodes specified by a name (pyTree) <Examples/Converter/getNodesFromNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodesFromNamePT.py
 
-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodeFromName(t, name)

    Return the first node found of given name.
    Starting node must be a standard pyTree node. 
    Wildcards are NOT accepted. 
    To accelerate search, this function can be limited to 1, 2 or 3 levels from 
    the starting node using getNodeFromName1, getNodeFromName2, getNodeFromName3. 
    If not found, it returns None. This is a fast routine.

    :param t:  starting node
    :type  t:  pyTree node
    :param name:  node name we are looking for
    :type  name:  string
    :return: node that matches name (shared with t)
    :rtype: pyTree node

    *Example of use:*

    * `Return the first found node specified by a name (pyTree) <Examples/Converter/getNodeFromNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodeFromNamePT.py
    
 
-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getByName(t, name, recursive=-1)

    Return a standard node containing the list of matching name nodes as children. 
    Wildcards are accepted.
   
    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :param name:  node name we are looking for
    :type  name:  string
    :param recursive: if set to 1,2,3, search is limited to 1,2,3 levels. If set to -1, the search is not limited.
    :type recursive: int
    :return: standard node with children list set as list of nodes that match name
    :rtype: pyTree node

    *Example of use:* 

    * `Return list of nodes specified by name as a standard node (pyTree) <Examples/Converter/getByNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getByNamePT.py

---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getChildFromName(node)

    Return the first child of node matching given name (one level search). 
    If not found, return None.
   
    :param node:  input node
    :type  node:  pyTree node
    :return: found node or None
    :rtype: pyTree node

    *Example of use:*

    * `Return first child of name (pyTree) <Examples/Converter/getChildFormNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getChildFromNamePT.py

    .. note:: new in version 3.2.

-----------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodesFromType(t, ntype)

    Return a list of nodes matching given type.
    To accelerate search, this function can be limited to 1, 2 or 3 levels from 
    the starting node using getNodesFromType1, getNodesFromType2, getNodesFromType3. 

    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :param ntype:  node type we are looking for
    :type  ntype:  string
    :return: list of nodes that matches given type (shared with t)
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return list of nodes specified by a type (pyTree) <Examples/Converter/getNodesFromTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodesFromTypePT.py
 
-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodeFromType(t, ntype)

    Return the first node found of given type.
    Starting node must be a standard pyTree node. 
    To accelerate search, this function can be limited to 1, 2 or 3 levels from 
    the starting node using getNodeFromType1, getNodeFromType2, getNodeFromType3.
    If not found, it returns None. This is a fast routine.

    :param t:  starting node
    :type  t:  pyTree node
    :param ntype:  node type we are looking for
    :type  ntype:  string
    :return: node that matches type (shared with t)
    :rtype: pyTree node

    *Example of use:*

    * `Return node specified by a type (pyTree) <Examples/Converter/getNodeFromTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodeFromTypePT.py
    

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getByType(t, ntype, recursive=-1)

    Return a standard node containing the list of matching type nodes as children.
   
    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :param ntype:  node type we are looking for
    :type  ntype:  string
    :param recursive: if set to 1,2,3, search is limited to 1,2,3 levels. If set to -1, the search is not limited.
    :type recursive: int
    :return: standard node with children list set as list of nodes that match type
    :rtype: pyTree node

    *Example of use:*

    * `Return list of nodes specified by type as a standard node (pyTree) <Examples/Converter/getByTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getByTypePT.py

---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getChildFromType(node)

    Return the first child of node matching given type (one level search). 
    If not found, return None.
   
    :param node:  input node
    :type  node:  pyTree node
    :return: found node or None
    :rtype: pyTree node

    *Example of use:*

    * `Return first child of type (pyTree) <Examples/Converter/getChildFormTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getChildFromTypePT.py

    .. note:: new in version 3.2.

-----------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodesFromNameAndType(t, name, ntype)

    Return a list of nodes matching given name and type.
    Wildcards are accepted for name and type.
    
    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :param name: node name we are looking for
    :type name: string
    :param ntype: node type we are looking for
    :type  ntype: string
    :return: list of nodes that matches given name and type (shared with t)
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return list of nodes specified by a name and type (pyTree) <Examples/Converter/getNodesFromNameAndTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodesFromNameAndTypePT.py
 
    .. note:: new in version 2.7.

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodeFromNameAndType(t, name, ntype)

    Return the first node found of given name and type.
    Starting node must be a standard pyTree node. 
    If not found, it returns None.

    :param t:  starting node
    :type  t:  pyTree node
    :param name: node name we are looking for
    :type name: string
    :param ntype: node type we are looking for
    :type  ntype: string
    :return: node that matches name and type (shared with t)
    :rtype: pyTree node

    *Example of use:*

    * `Return node specified by a name and type (pyTree) <Examples/Converter/getNodeFromNameAndTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodeFromNameAndTypePT.py
    
    .. note:: new in version 2.7.

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodesFromValue(t, value)

    Return a list of nodes matching given value. if value is a string, wildcards
    are accepted.
   
    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :param value:  node value we are looking for
    :type  value:  string, int, float, numpy
    :return: list of nodes that match given value (shared with t)
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return list of nodes specified by value (pyTree) <Examples/Converter/getNodesFromValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodesFromValuePT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getParentOfNode(t, node)

    Return the parent node of given node in t. t must be a higher node in the
    tree. It returns p (the parent node) and c (the position number of node
    in the parent's children list).
    If t is a node of a pyTree, then p[2][c] = node.
    If t is a list of standard nodes, then p == t and p[c] = node.
    If node is not found, then p is None.
    This in an expansive routine. Prefer top-down tree traversal, if
    possible.
   
    :param t:  higher node
    :type  t:  pyTree node
    :param node:  node the parent of which is looked for
    :type  node:  pyTree node

    :return: parent node and position in children list
    :rtype: tuple (pyTree node, c)

    *Example of use:*

    * `Return parent node (pyTree) <Examples/Converter/getParentOfNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getParentOfNodePT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getParentFromType(t, node, parentType)

    Return the parent node of given node in t. The parent can be any level upper.
    t must be a higher node in the tree. 
    It returns the first parent node matching that type.
   
    :param t:  higher node
    :type  t:  pyTree node
    :param node:  node the parent of which is looked for
    :type  node:  pyTree node
    :param parentType: type of parent
    :type parentType: string

    :return: parent node 
    :rtype: pyTree node

    *Example of use:*

    * `Return parent node of given type (pyTree) <Examples/Converter/getParentFromTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getParentFromTypePT.py

    .. note:: new in version 2.7.

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getParentsFromType(t, node, parentType)

    Return a list of parent nodes of given node in t of given type. The parent can be any level upper.
    t must be a higher node in the tree. 
   
    :param t:  higher node
    :type  t:  pyTree node
    :param node:  node the parent of which is looked for
    :type  node:  pyTree node
    :param parentType: type of parent
    :type parentType: string

    :return: parent node list 
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return parent nodes of given type (pyTree) <Examples/Converter/getParentsFromTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getParentsFromTypePT.py

    .. note:: new in version 2.7.

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodePosition(node, parent)

    Return the position of node in parent children list.
    If node is not found, return -1.

    :param node:  node which position need to be found 
    :type  node:  pyTree node
    :param parent:  the parent of node
    :type  parent:  pyTree node
    
    :return: node position in parent children 
    :rtype: int

    *Example of use:*

    * `Return position of a node (pyTree) <Examples/Converter/getNodePositionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodePositionPT.py

    .. note:: new in version 2.9.

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getNodeFromPath(t, path)

    Return a node from its path string. Note that node path is
    relative to input node. If not found, it returns None.
   
    :param t:  starting node
    :type  t:  pyTree node
    :param path:  path ('Base/Zone')
    :type  path:  string

    :return: found node (shared with t)
    :rtype: pyTree node

    *Example of use:*

    * `Return node from path (pyTree) <Examples/Converter/getNodeFromPathPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getNodeFromPathPT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getPathsFromName(node, name, pyCGNSLike=False)

    Return a list of paths corresponding to a given name. The 
    paths are built from node. Wildcards are accepted for name.
   
    :param node:  starting node
    :type  node:  pyTree node or list of pyTree nodes
    :param name:  node name we are looking for
    :type  type:  string
    :param pyCGNSLike:  if True, paths don't start with CGNSTree
    :type  pyCGNSLike:  boolean
    :return: list of paths of nodes matching given name
    :rtype: list of strings

    *Example of use:*

    * `Return list of paths specified by name (pyTree) <Examples/Converter/getPathsFromNamesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getPathsFromNamesPT.py

    .. note:: new in version 2.5.

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getPathsFromType(node, ntype, pyCGNSLike=False)

    Return a list of paths corresponding to a given type. The 
    paths are built from node.
   
    :param node:  starting node
    :type  node:  pyTree node or list of pyTree nodes
    :param ntype:  node type we are looking for
    :type  ntype:  string
    :param pyCGNSLike:  if True, paths don't start with CGNSTree
    :type  pyCGNSLike:  boolean
    :return: list of paths of nodes matching given type
    :rtype: list of strings

    *Example of use:*

    * `Return list of paths specified by type (pyTree) <Examples/Converter/getPathsFromTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getPathsFromTypePT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getPathsFromValue(node, value, pyCGNSLike=False)

    Return a list of paths corresponding to a given value. The 
    paths are built from node. If value is a string, wildcards are accepted.
   
    :param node:  starting node
    :type  node:  pyTree node or list of pyTree nodes
    :param value:  node value we are looking for
    :type  value:  string, int, float, numpy
    :param pyCGNSLike:  if True, paths don't start with CGNSTree
    :type  pyCGNSLike:  boolean
    :return: list of paths of nodes matching given type
    :rtype: list of strings

    *Example of use:*

    * `Return list of paths specified by value (pyTree) <Examples/Converter/getPathsFromValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getPathsFromValuePT.py

    .. note:: New in version 2.5

------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getPathLeaf(path)

    Return the last term of path.
   
    :param path:  input path
    :type  path:  string
    :return: last term of path
    :rtype: string

    *Example of use:*

    * `Return leaf of path (pyTree) <Examples/Converter/getPathLeafPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getPathLeafPT.py

    .. note:: New in version 2.5

------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getPathAncestor(path, level=-1)

    Return the path ancestor of path.
   
    :param path:  input path
    :type  path:  string
    :param level: number of level to go up
    :type level: int
    :return: path of ancestor
    :rtype: string

    *Example of use:*

    * `Return ancestor path (pyTree) <Examples/Converter/getPathAncestorPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getPathAncestorPT.py

    .. note:: New in version 2.5

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getZonePaths(t, pyCGNSLike=False)

    Return the list of paths of Zone_t nodes. Paths are relative
    to node t.

    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :param pyCGNSLike:  if True, paths don't start with CGNSTree
    :type  pyCGNSLike:  boolean
    :return: list of Zone_t paths
    :rtype: list of strings

    *Example of use:*

    * `Return zone paths of tree (pyTree) <Examples/Converter/getZonePathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getZonePathsPT.py

    .. note:: New in version 2.5

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getZones(t)

    Return the list of Zone_t nodes.
   
    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :return: list of Zone_t nodes (shared with t)
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return zone nodes of tree (pyTree) <Examples/Converter/getZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getZonesPT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getZonesPerIteration(t, iteration=None, time=None)

    Return the list of Zone_t nodes matching a given iteration.
    A BaseIterativeData node must exist in input pyTree.
    If iteration is provided, return a list of zones matching iteration.
    If time is provided, return a list of zones mathcing time
    If none is set, return a list of list of zones for each iterations.

    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :param iteration: iteration number
    :type  iteration: int
    :param time: desired time
    :type  time: float
    :return: list of Zone_t nodes (shared with t) or list of list of zones
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return zone nodes of tree of given iteration (pyTree) <Examples/Converter/getZonesPerIterationPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getZonesPerIterationPT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getBases(t)

    Return the list of CGNSBase_t nodes.
   
    :param t:  starting node
    :type  t:  pyTree node or list of pyTree nodes
    :return: list of CGNSBase_t nodes (shared with t)
    :rtype: list of pyTree nodes

    *Example of use:*

    * `Return base nodes of tree (pyTree) <Examples/Converter/getBasesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getBasesPT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getZoneDim(zone) 

    Return the dimension of a zone node.
    Return
    ['Structured', ni, nj, nk, celldim] for structured grids, return
    ['Unstructured', np, ne, eltsName, celldim] for unstructured grids.
    np is the number of points, ne the number of elements, celldim is
    0, 1, 2, 3 depending on element dimension, eltsName is one of NODE,
    BAR, TRI, QUAD, TETRA, PYRA, PENTA, NGON, MULTIPLE.
   
    :param zone:  a 'Zone_t' node
    :type  zone:  pyTree node of type 'Zone_t'
    :return: ['Structured', ni, nj, nk, celldim] or ['Unstructured', np, ne, eltsName, celldim]
    :rtype: list of data

    *Example of use:*

    * `Return zone dimension information (pyTree) <Examples/Converter/getZoneDimPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getZoneDimPT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getZoneType(zone) 

    Return the nature of zone.
    Return 1 if zone is structured, 2 if unstructured, 0 if failed.

    :param zone:  a 'Zone_t' node
    :type  zone:  pyTree node of type 'Zone_t'
    :return: 1 or 2
    :rtype: int

    *Example of use:*

    * `Return zone nature (pyTree) <Examples/Converter/getZoneTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getZoneTypePT.py


Check nodes
------------


.. py:function:: Converter.Internal.printTree(node, file=None, stdOut=None, editor=None, color=False) 

    Pretty print a pyTree or a pyTree node to screen or in file.
   
    :param node:  input node
    :type  node:  pyTree node of list of pyTree nodes
    :param file:  file name (optional)
    :type file: string
    :param stdOut: a file object (as created by open)(optional)
    :type stdOut: file object
    :param editor: an editor name (optional)
    :type editor: string
    :param color: trigger colored output with ANSI codes (optional)
    :type color: True or False
    
    *Example of use:*

    * `Pretty print pyTree nodes (pyTree) <Examples/Converter/printTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/printTreePT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.getSizeOf(node) 

    Return the size of input node and attached nodes in octets.

    :param node:  input node
    :type  node:  pyTree node of list of pyTree nodes
    :return: size of node in octets
    :rtype: int

    *Example of use:*

    * `Get size of node (pyTree) <Examples/Converter/getSizeOfPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getSizeOfPT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.checkPyTree(t, level=-20) 

    Check pyTree t following level (0: valid version node,
    1: node conformity, 2: unique base name,
    3: unique zone name, 4: unique BC name, 5: valid BC range, 6: valid
    opposite BC range for match and nearmatch, 7: referenced familyZone and
    familyBCs must be defined in bases, 8: valid CGNS types, 9: valid connectivity,
    10: valid CGNS flowfield name and dimension).
    If level=-20, all previous checks are performed.
    
    Return a list of pairs of invalid nodes and error message.
   
    :param t:  input pyTree
    :type  t:  pyTree node or list of pyTree nodes
    :param level:  level of check
    :type  level:  int
    :return: [(node, 'errror message')]
    :rtype: list of tuple (node, 'error')

    *Example of use:*

    * `Check pyTree (pyTree) <Examples/Converter/checkPyTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/checkPyTreePT.py

-----------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.correctPyTree(t, level=-20) 

    Correct a pyTree t following level (0: valid version node,
    1: node conformity, 2: unique base name,
    3: unique zone name, 4: unique BC name, 5: valid BC range, 6: valid
    opposite BC range for match and nearmatch, 7: referenced familyZone and
    familyBCs must be defined in bases, 8: valid CGNS types, 9: valid connectivity,
    10: valid CGNS flowfield name and dimension). 

    Generally invalid nodes are suppressed.
    If level=-20, all previous checks are performed.

    Exists also as in place version (_correctPyTree) that modifies t and returns None.
   
    :param t: input pyTree
    :type  t: pyTree node or list of pyTree nodes
    :param level:  level of check
    :type  level:  int 
    :return: modified reference copy of t
    :rtype: same as input

    *Example of use:*

    * `Correct pyTree (pyTree) <Examples/Converter/correctPyTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/correctPyTreePT.py

Copy nodes
-----------

.. py:function:: Converter.Internal.copyRef(node) 

    Copy recursively sharing node values (in particular data numpys are shared).

    :param node:  input data
    :type  node:  pyTree node or list of pyTree nodes
    :return: copy of input data
    :rtype: same as input

    *Example of use:*

    * `Copy pyTree with references (pyTree) <Examples/Converter/copyRefPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/copyRefPT.py

--------------------------------------------------------------------------

.. py:function:: Converter.Internal.copyTree(node) 

    Fully copy a tree. Node values (in particular data numpys) are copied.

    :param node:  input data
    :type  node:  pyTree node or list of pyTree nodes
    :return: copy of input data
    :rtype: same as input

    *Example of use:*

    * `Copy pyTree (pyTree) <Examples/Converter/copyTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/copyTreePT.py

--------------------------------------------------------------------------

.. py:function:: Converter.Internal.copyValue(node, byName=None, byType=None) 

    Copy recursively the value of nodes specified by byName or byType string. 
    Wildcards are accepted.

    :param node:  input data
    :type  node:  pyTree node or list of pyTree nodes
    :param byName: name of nodes to be copied or None
    :type byName: string
    :param byType: type of nodes to be copied or None
    :type byType: string
    :return: copy of input data
    :rtype: same as input

    *Example of use:*

    * `Copy value (pyTree) <Examples/Converter/copyValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/copyValuePT.py

--------------------------------------------------------------------------

.. py:function:: Converter.Internal.copyNode(node) 

    Copy only this node (no recursion). Node value (in particular data numpys) is copied.

    :param node:  input node
    :type  node:  pyTree node or list of pyTree nodes
    :return: copy of input node
    :rtype: same as input

    *Example of use:*

    * `Copy node (pyTree) <Examples/Converter/copyNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/copyNodePT.py



Add/remove node
---------------

.. py:function:: Converter.Internal.append(t, node, path) 

    Append a node by specifying its path. Note that the path is relative to t.

    Exists also as in place version (_append) that modifies t and returns None.
   

    :param t: starting node
    :type  t: pyTree node
    :param node:  node to be added
    :type  node:  pyTree node
    :param path: the path where node should be added ('Base/')
    :return: modified reference copy of t
    :rtype: same as t

    *Example of use:*

    * `Append a node by its path (pyTree) <Examples/Converter/appendPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/appendPT.py

--------------------------------------------------------------------------

.. py:function:: Converter.Internal.rmNode(t, node)

    Remove given node in t. t is modified.

    :param t: starting node
    :type  t: pyTree node
    :param node:  node to be removed
    :type  node:  pyTree node
    :rtype: None

    *Example of use:*

    * `Remove a node (pyTree) <Examples/Converter/rmNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmNodePT.py

--------------------------------------------------------------------------

.. py:function:: Converter.Internal.rmNodeByPath(t, path) 

    Remove node that corresponds to path in t.

    Exists also as in place version (_rmNodeByPath) that modifies t and returns None.

    :param t: input node
    :type  t: pyTree node
    :param path:  path of node to be removed (ex: 'Base/Zone0')
    :type  path:  string
    :return: reference copy of t with node removed
    :rtype: pyTree node

    *Example of use:*

    * `Remove node by its path (pyTree) <Examples/Converter/rmNodeByPathPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmNodeByPathPT.py

----------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.rmNodesByName(t, name) 

    Remove all nodes in t that match given name.

    Exists also as in place version (_rmNodesByName) that modifies t and returns None.
    Exists also as in place with search limited to 1 or 2 levels as _rmNodesByName1
    and _rmNodesByName2.

    :param t: input node
    :type  t: pyTree node or list of pyTree nodes
    :param name:  name of nodes to be removed (wildcards are accepted)
    :type  name:  string
    :return: reference copy of t with nodes removed
    :rtype: same as t

    *Example of use:*

    * `Remove nodes that match given name (pyTree) <Examples/Converter/rmNodesByNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmNodesByNamePT.py

----------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.rmNodesByType(t, ntype) 

    Remove all nodes in t that match given type.

    Exists also as in place version (_rmNodesByType) that modifies t and returns None.
    Exists also as in place with search limited to 1 or 2 levels as _rmNodesByType1
    and _rmNodesByType2.

    :param t: input node
    :type  t: pyTree node or list of pyTree nodes
    :param ntype:  type of nodes to be removed
    :type  ntype:  string
    :return: reference copy of t with nodes removed
    :rtype: same as t

    *Example of use:*

    * `Remove nodes that match given type (pyTree) <Examples/Converter/rmNodesByTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmNodesByTypePT.py

----------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.rmNodesByNameAndType(t, name, ntype) 

    Remove all nodes that match given type and given name at the same time.

    Exists also as in place version (_rmNodesByNameAndType) that modifies t and returns None.

    :param t: input node
    :type  t: pyTree node or list of pyTree nodes
    :param name:  name of nodes to be removed (wildcards accepted)
    :type  name:  string
    :param ntype:  type of nodes to be removed
    :type  ntype:  string
    :return: reference copy of t with nodes removed
    :rtype: same as t

    *Example of use:*

    * `Remove nodes that match given name and type (pyTree) <Examples/Converter/rmNodesByNameAndTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmNodesByNameAndTypePT.py

----------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.rmNodesByValue(t, value) 

    Remove all nodes in t that match given value.

    Exists also as in place version (_rmNodesByValue) that modifies t and returns None.

    :param t: input node
    :type  t: pyTree node or list of pyTree nodes
    :param value:  value of nodes to be removed
    :type  value:  string, number...
    :return: reference copy of t with nodes removed
    :rtype: same as t

    *Example of use:*

    * `Remove nodes that match given value (pyTree) <Examples/Converter/rmNodesByValuePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmNodesByValuePT.py


----------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.moveNodeFromPaths(t, path1, path2) 

    Move node located at path1 in t to path2.

    Exists also as in place version (_moveNodeFromPaths) that modifies t and returns None.

    :param t: input node
    :type  t: pyTree node
    :param path1: initial path of node to move
    :type  path1: string
    :param path2: destination path
    :type  path2: string
    :return: reference copy of t with node moved
    :rtype: same as t

    .. note:: new in version 3.2.

    *Example of use:*

    * `Move node from path1 to path2 (pyTree) <Examples/Converter/moveNodeFromPathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/moveNodeFromPathsPT.py


Modify nodes
-------------

.. py:function:: Converter.Internal.merge(A)

    Merge a list of trees defined in [t1, t2,...]. Return a merged tree.
    If a node appears more than once in different trees of the list,
    the first found node is kept.

    :param A: list of trees to be merged
    :type  A: list of pyTrees
    :return: merged tree (reference copy)
    :rtype: pyTree

    *Example of use:*

    * `Merge trees (pyTree) <Examples/Converter/mergePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/mergePT.py

----------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.renameNode(t, name, newName) 

    Rename node named 'name' with 'newName' in t. Occurances of 
    name elsewhere in the tree t is replaced with 'newName'.
    Widlcards are accepted for name. 

    Exists also as in place version (_renameNode) that modifies t and returns None.

    :param t: input node 
    :type  t: pyTree node or list of pyTree nodes
    :return: reference copy of t
    :rtype: same as t

    *Example of use:*

    * `Rename node (pyTree) <Examples/Converter/renameNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/renameNodePT.py

----------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.sortByName(t, recursive=True) 

    Sort nodes by their name (alphabetical order). If recursive=True,
    sort also chidren nodes.

    Exists also as in place version (_sortByName) that modifies t and returns None.

    :param t: starting pyTree node
    :type  t: pyTree node or list of pyTree nodes
    :param recursive: if True, perform sort also on all children nodes
    :type  recursive: boolean
    :return: reference copy of t
    :rtype: same as t

    *Example of use:*

    * `Sort nodes by name (pyTree) <Examples/Converter/sortByNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/sortByNamePT.py

------------------------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.appendBaseName2ZoneName(t, updateRef=True, separator='_', trailing='') 

    Append base name to zone name (resulting name is baseName_zoneName for each
    zone). If updateRef=True, reference to zone names in BC,... are replaced.
    Separator between base and zone name can be set with separator, a trailing
    string can also be added.

    Exists also as in place version (_appendBaseName2ZoneName) that modifies t and returns None.

    :param t: input date
    :type  t: top pyTree node or list of Base nodes
    :param updateRef: True if reference to zone name has to be changed
    :type updateRef: Boolean
    :param separator: separator between base and zone name
    :type separator: string
    :param trailing: trailing to be optionally added to merge name
    :type trailing: string
    :return: reference copy of input
    :rtype: same as input

    *Example of use:*

    * `Append base name to zone name (pyTree) <Examples/Converter/appendBaseName2ZoneNamePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/appendBaseName2ZoneNamePT.py

----------------------------------------------------------------------------------------

.. py:function:: Converter.Internal.groupBCByBCType(t, btype='BCWall', name='FamWall') 

    For each base, gather all BCs of same BCType in a family named FamilyName.
    The family node is created and the BC is tagged with BC family.

    Exists also as in place version (_groupBCByBCType) that modifies t and returns None.

    :param t: input data
    :type  t: pyTree or list of Base nodes
    :param btype: Type of BC to be grouped in a family
    :type btype: string
    :param name: name of family to be created
    :type name: string
    :return: reference copy of input
    :rtype: same as input

    *Example of use:*

    * `Group BC from their types (pyTree) <Examples/Converter/groupBCByBCTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/groupBCByBCTypePT.py

    .. note:: New in version 2.4


Create specific CGNS nodes
----------------------------

.. py:function:: Converter.Internal.newCGNSTree()

    Create a tree node with the CGNS version node.
    
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new tree node (pyTree) <Examples/Converter/newCGNSTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newCGNSTreePT.py

---------------------------------------------------------------------------------

.. py:function:: Converter.Internal.newCGNSBase(name='Base', cellDim=3, physDim=3, parent=None) 

    Create a base node. cellDim is the dimension of zone cells of this
    base, physDim is the
    dimension of physical space. For example, triangle zones in 3D space will
    correspond to cellDim=2 and physDim=3. If parent is not None, attach
    it to parent node.

    :param name: name of base
    :type name: string
    :param cellDim: dimension of zone cells in Base (1, 2, 3)
    :type cellDim: int
    :param physDim: physical dimension of space (1, 2, 3)
    :type physDim: int
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new base node (pyTree) <Examples/Converter/newCGNSBasePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newCGNSBasePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newZone(name='Zone', zsize=None, ztype='Structured', family=None, parent=None) 

    Create a zone node. zsize is the dimension of zone, ztype is the type of
    zone ('Structured' or 'Unstructured'). If family is not None, attach a
    zone family node. If parent is not None, attach it to parent node.

    :param name: name of zone
    :type name: string
    :param zsize: number of points, number of elements, 0
    :type zsize: list of integers
    :param ztype: 'Structured' or 'Unstructured'
    :type ztype: string
    :param family: optional family name
    :type family: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new zone node (pyTree) <Examples/Converter/newZonePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newZonePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGridCoordinates(name=Internal.__GridCoordinates__, parent=None) 

    Create a GridCoordinates node.
    If parent is not None, attach it to parent node.

    :param name: name of GridCoordinates container
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GridCoordinates node (pyTree) <Examples/Converter/newGridCoordinatesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGridCoordinatesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newDataArray(name='Data', value=None, parent=None) 

    Create a DataArray node. value can be a string, an int, a float,
    a numpy of ints, a numpy of floats.
    If parent is not None, attach it to parent node.

    :param name: name of node
    :type name: string
    :param value: value to put in node
    :type value: string, int, float, numpy
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new DataArray node (pyTree) <Examples/Converter/newDataArrayPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newDataArrayPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newDataClass(value='UserDefined', parent=None) 

    Create a DataClass node.
    If parent is not None, attach it to parent node.

    :param value: value to put in node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new DataClass node (pyTree) <Examples/Converter/newDataClassPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newDataClassPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newDimensionalUnits(massUnit='Kilogram', lengthUnit='Meter', timeUnit='Second', temperatureUnit='Kelvin', angleUnit='Radian', parent=None) 

    Create a DimensionalUnits node. Arguments describe the units of the problem. 
    If parent is not None, attach it to parent node.

    :param massUnit: mass unit you want
    :type massUnit: string in 'Null', 'UserDefined', 'Kilogram', 'Gram', 'Slug', 'PoundMass',
    :param lengthUnit: length unit you want
    :type lengthUnit: string in 'Null', 'UserDefined', 'Meter', 'Centimeter', 'Millimeter', 'Foot', 'Inch'
    :param timeUnit: time unit you want
    :type timeUnit: string in 'Null', 'UserDefined', 'Second'
    :param temperatureUnit: temprature unit you want
    :type temperatureUnit: string in 'Null', 'UserDefined', 'Kelvin', 'Celsius', 'Rankine'
    :param angleUnit: angle unit you want
    :type angleUnit: string in 'Null', 'UserDefined', 'Radian', 'Degree'
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new DimensionalUnits node (pyTree) <Examples/Converter/newDimensionalUnitsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newDimensionalUnitsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newDimensionalExponents(massExponent='Kilogram', lengthExponent='Meter', timeExponent='Second', temperatureExponent='Kelvin', angleExponent='Radian', parent=None) 

    Create a DimensionalExponents node. Arguments describe the unit exponent of the problem. 
    If parent is not None, attach it to parent node.

    :param massExponent: exponent for mass
    :type massExponent: float
    :param lengthExponent: exponent for length
    :type lengthExponent: float
    :param timeExponent: exponent for time
    :type timeExponent: float
    :param temperatureExponent: exponent for temperature
    :type temperatureExponent: float
    :param angleExponent: exponent for angle
    :type angleExponent: float
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new DimensionalExponents node (pyTree) <Examples/Converter/newDimensionalExponentsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newDimensionalExponentsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newDataConversion(conversionScale=1., conversionOffset=0., parent=None) 

    Create a DataConversion node. Arguments describe the conversion factors. 
    Data(raw) = Data(nondimensional)*ConversionScale + ConversionOffset
    If parent is not None, attach it to parent node.

    :param conversionScale: scale for conversion
    :type conversionScale: float
    :param conversionOffset: offset for conversion
    :type conversionScale: float
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new DataConversion node (pyTree) <Examples/Converter/newDataConversionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newDataConversionPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newDescriptor(name='Descriptor', value='', parent=None) 

    Create a Descriptor node. 
    If parent is not None, attach it to parent node.

    :param name: name of node
    :type name: string
    :param value: value to put in node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new Descriptor node (pyTree) <Examples/Converter/newDescriptorPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newDescriptorPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGridLocation(value='CellCenter', parent=None)

    Create a GridLocation node. 
    If parent is not None, attach it to parent node.

    :param value: value to put in node ('CellCenter', 'Vertex'...)
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GridLocation node (pyTree) <Examples/Converter/newGridLocationPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGridLocationPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newIndexArray(name='Index', value=None, parent=None) 

    Create a indexArray node. 
    If parent is not None, attach it to parent node.

    :param name: name of node
    :type name: string
    :param value: integer values of indices
    :type value: list of integers
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new IndexArray node (pyTree) <Examples/Converter/newIndexArrayPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newIndexArrayPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newPointList(name='PointList', value=None, parent=None) 

    Create a PointList node. 
    If parent is not None, attach it to parent node.

    :param name: name of node
    :type name: string
    :param value: list of point indices
    :type value: list of integers
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new PointList node (pyTree) <Examples/Converter/newPointListPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newPointListPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newPointRange(name='PointRange', value=None, parent=None) 

    Create a PointRange node. 
    If parent is not None, attach it to parent node.

    :param name: name of node
    :type name: string
    :param value: list of point indices ([imin,imax,jmin,jmax,kmin,kmax] for a structured point range)
    :type value: list of integers
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new PointRange node (pyTree) <Examples/Converter/newPointRangePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newPointRangePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newRind(value=None, parent=None) 

    Create a Rind node. Rind contains the number of ghost cells for a structured block.
    If parent is not None, attach it to parent node.

    :param value: list of integers (6 for a structured block)
    :type value: list of integers
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new Rind node (pyTree) <Examples/Converter/newRindPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newRindPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newSimulationType(value='TimeAccurate', parent=None) 

    Create a SimulationType node. 
    If parent is not None, attach it to parent node.

    :param value: value of the node ('TimeAccurate','NonTimeAccurate')
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new SimulationType node (pyTree) <Examples/Converter/newSimulationTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newSimulationTypePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newOrdinal(value=0, parent=None) 

    Create an Ordinal node. 
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: integer
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new Ordinal node (pyTree) <Examples/Converter/newOrdinalPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newOrdinalPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newDiscreteData(name='DiscreteData', parent=None) 

    Create a DiscreteData node. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new DiscreteData node (pyTree) <Examples/Converter/newDiscreteDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newDiscreteDataPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newIntegralData(name='IntegralData', parent=None) 

    Create an IntegralData node. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new IntegralData node (pyTree) <Examples/Converter/newIntegralDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newIntegralDataPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newElements(name='Elements', etype='UserDefined', econnectivity=None, erange=None , eboundary=0, parent=None) 

    Create a Elements node. etype is the element type ('BAR', 'TRI', 'QUAD', ...), econnectivity is the connectivity numpy array and eboundary ... 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param etype: type of the elements ('BAR', 'TRI', 'QUAD', 'TETRA', 'HEXA', 'PENTA', 'NODE', 'PYRA', 'NGON', 'NFACE')
    :type etype: string
    :param econnectivity: connectivity of elements
    :type econnectivity: numpy array
    :param erange: element range
    :type erange: list of two integers
    :param eboundary: number of elements at boundary (0 by default)
    :type eboundary: integer
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new Elements node (pyTree) <Examples/Converter/newElementsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newElementsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newParentElements(value=None, parent=None) 

    Create a ParentElements node. value is a numpy array. 
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: numpy array
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ParentElements node (pyTree) <Examples/Converter/newParentElementsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newParentElementsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newParentElementsPosition(value=None, parent=None) 

    Create a ParentElementsPosition node. value is a numpy array. 
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: numpy array
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ParentElementsPosition node (pyTree) <Examples/Converter/newParentElementsPositionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newParentElementsPositionPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newZoneBC(parent=None)  

    Create a ZoneBC node. 
    If parent is not None, attach it to parent node.

    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ZoneBC node (pyTree) <Examples/Converter/newZoneBCPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newZoneBCPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newBC(name='BC', pointRange=None, pointList=None, btype='Null', family=None, parent=None)  

    Create a BC node. It can be defined by a pointRange for structured zones or a pointList of faces for unstructured zones. btype specifies the BC type. A BC familyName can also be defined. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param pointRange: list of point indices ([imin,imax,jmin,jmax,kmin,kmax] for a structured point range)
    :type pointRange: list of integers
    :param pointList: list of point indices (for unstructured grid)
    :type pointList: list of integers
    :param btype: type of BC (BCWall, BCFarfield, FamilySpecified...)
    :type btype: string
    :param family: name of the family for a family specified BC
    :type family: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new BC node (pyTree) <Examples/Converter/newBCPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newBCPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newBCDataSet(name='BCDataSet', value='Null', gridLocation=None, parent=None)  

    Create a BCDataSet node. value must be a BCType. GridLocation ('FaceCenter', 'Vertex') can be specified. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param value: value of node (UserDefined, BCWall...)
    :type value: string
    :param gridLocation: location of the grid points (Vertex, FaceCenter, CellCenter...)
    :type gridLocation: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new BCDataSet node (pyTree) <Examples/Converter/newBCDataSetPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newBCDataSetPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newBCData(name='BCData', parent=None)   

    Create a BCData node. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new BCData node (pyTree) <Examples/Converter/newBCDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newBCDataPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newBCProperty(wallFunction='Null', area='Null', parent=None)   

    Create a BCProperty node. 
    If parent is not None, attach it to parent node.

    :param wallFunction: type of wall function (Null, UserDefined, Generic)
    :type wallFunction: string
    :param area: type of area (Null, UserDefined, BleedArea, CaptureArea)
    :type area: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new BCProperty node (pyTree) <Examples/Converter/newBCPropertyPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newBCPropertyPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newAxiSymmetry(referencePoint=[0.,0.,0.], axisVector=[0.,0.,0.], parent=None)   

    Create a AxiSymmetry node. 
    If parent is not None, attach it to parent node.

    :param referencePoint: coordinates of the axis point
    :type referencePoint: list of 3 floats
    :param axisVector: axis vector
    :type axisVector: list of 3 floats
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new AxiSymmetry node (pyTree) <Examples/Converter/newAxiSymmetryPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newAxiSymmetryPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newRotatingCoordinates(rotationCenter=[0.,0.,0.], rotationRateVector=[0.,0.,0.], parent=None)   

    Create a RotatingCoordinates node. 
    If parent is not None, attach it to parent node.

    :param rotationCenter: coordinates of the rotation center
    :type rotationCenter: list of 3 floats
    :param rotationRateVector: vector of rotation rate
    :type rotationRateVector: list of 3 floats
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new RotatingCoordinates node (pyTree) <Examples/Converter/newRotatingCoordinatesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newRotatingCoordinatesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newFlowSolution(name=__FlowSolutionNodes__, gridLocation='Vertex', parent=None)   

    Create a newFlowSolution node. 
    If parent is not None, attach it to parent node.

    :param name: name or the container of the flow solution node
    :type name: string or container
    :param gridLocation: location of the node grid points
    :type gridLocation: string (Vertex, CellCenter...)
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new FlowSolution node (pyTree) <Examples/Converter/newFlowSolutionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newFlowSolutionPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newZoneGridConnectivity(name='ZoneGridConnectivity', parent=None)   

    Create a newZoneGridConnectivity node. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ZoneGridConnectivity node (pyTree) <Examples/Converter/newZoneGridConnectivityPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newZoneGridConnectivityPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGridConnectivity1to1(name='Match', donorName=None, pointRange=None, pointList=None, pointRangeDonor=None, pointListDonor=None, transform=None, parent=None)   

    Create a newGridConnectivity1to1 node. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param donorName: name of donor zone
    :type donorName: string
    :param pointRange: list of point indices of the local zone ([imin,imax,jmin,jmax,kmin,kmax] for a structured point range)
    :type pointRange: list of integers
    :param pointList: list of point indices of the local zone (for unstructured grid)
    :type pointList: list of integers
    :param pointRangeDonor: list of point indices of the donor zone ([imin,imax,jmin,jmax,kmin,kmax] for a structured point range)
    :type pointRangeDonor: list of integers
    :param pointListDonor: list of point indices of the donor zone (for unstructured grid)
    :type pointListDonor: list of integers
    :param transform: transformation of the orientation between local and donor zones
    :type transform: list of integers ([1,2,3], [-2,-1,-3]...)
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GridConnectivity1to1 node (pyTree) <Examples/Converter/newGridConnectivity1to1PT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGridConnectivity1to1PT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGridConnectivity(name='Overlap', donorName=None, ctype='Overset', parent=None)   

    Create a newGridConnectivity node. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param donorName: name of donor zone
    :type donorName: string
    :param ctype: connectivity type ('Overset')
    :type ctype: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GridConnectivity node (pyTree) <Examples/Converter/newGridConnectivityPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGridConnectivityPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGridConnectivityType(ctype='Overset',parent=None)   

    Create a newGridConnectivityType node. 
    If parent is not None, attach it to parent node.

    :param ctype: connectivity type ('Overset')
    :type ctype: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GridConnectivityType node (pyTree) <Examples/Converter/newGridConnectivityTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGridConnectivityTypePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGridConnectivityProperty(parent=None)   

    Create a newGridConnectivityProperty node. 
    If parent is not None, attach it to parent node.

    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GridConnectivityProperty node (pyTree) <Examples/Converter/newGridConnectivityPropertyPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGridConnectivityPropertyPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newPeriodic(rotationCenter=[0.,0.,0.], rotationAngle=[0.,0.,0.], translation=[0.,0.,0.], parent=None)   

    Create a Periodic node. 
    If parent is not None, attach it to parent node.

    :param rotationCenter: coordinates of the rotation center
    :type rotationCenter: list of 3 floats
    :param rotationAngle: angles of rotation
    :type rotationAngle: list of 3 floats
    :param translation: translation vector
    :type translation: list of 3 floats
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new Periodic node (pyTree) <Examples/Converter/newPeriodicPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newPeriodicPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newZoneSubRegion(name='SubRegion', pointRange=None, pointList=None, bcName=None, gcName=None, gridLocation=None, parent=None)

    Create a ZoneSubRegion node. 
    If parent is not None, attach it to parent node.

    :param name: name of node
    :type name: string
    :param pointRange: list of point indices ([imin,imax,jmin,jmax,kmin,kmax] for a structured point range)
    :type pointRange: list of integers
    :param pointList: list of point indices (for unstructured grid)
    :type pointList: list of integers
    :param bcName:  name of the BC node to which is connected the ZoneSubRegion
    :type bcName: string
    :param gcName: name of the GridConnectivity node to which is connected the ZoneSubRegion
    :type gcName: string
    :param gridLocation: location of zoneSubRegion data  (Vertex, FaceCenter, CellCenter…)
    :type gridLocation: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ZoneSubRegion node (pyTree) <Examples/Converter/newZoneSubRegionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newZoneSubRegionPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newOversetHoles(name='OversetHoles', pointRange=None, pointList=None, parent=None)   

    Create a OversetHoles node. 
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param pointRange: list of point indices ([imin,imax,jmin,jmax,kmin,kmax] for a structured point range)
    :type pointRange: list of integers
    :param pointList: list of point indices (for unstructured grid)
    :type pointList: list of integers
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new OversetHoles node (pyTree) <Examples/Converter/newOversetHolesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newOversetHolesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newFlowEquationSet(parent=None)

    Create a FlowEquationSet node. 
    If parent is not None, attach it to parent node.

    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new FlowEquationSet node (pyTree) <Examples/Converter/newFlowEquationSetPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newFlowEquationSetPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGoverningEquations(value='Euler', parent=None)

    Create a GoverningEquations node. value is the equation type in 'FullPotential', 'Euler', 'NSLaminar','NSTurbulent', 'NSLaminarIncompressible', 'NSTurbulentIncompressible'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GoverningEquations node (pyTree) <Examples/Converter/newGoverningEquationsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGoverningEquationsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGasModel(value='Ideal', parent=None)

    Create a GasModel node. value is the model type in 'Ideal', 'VanderWaals', 'CaloricallyPerfect', 'ThermallyPerfect', 'ConstantDensity', 'RedlichKwong'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GasModel node (pyTree) <Examples/Converter/newGasModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGasModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newThermalConductivityModel(value='Null', parent=None)

    Create a ThermalConductivityModel node. value is the model type in 'ConstantPrandtl', 'PowerLaw', 'SutherlandLaw'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ThermalConductivityModel node (pyTree) <Examples/Converter/newThermalConductivityModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newThermalConductivityModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newViscosityModel(value='Null', parent=None)

    Create a ViscosityModel node. value is the model type in 'Constant', 'PowerLaw', 'SutherlandLaw'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ViscosityModel node (pyTree) <Examples/Converter/newViscosityModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newViscosityModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newTurbulenceClosure(value='Null', parent=None)

    Create a TurbulenceClosure node. value is the closure type in 'EddyViscosity', 'ReynoldStress', 'ReynoldsStressAlgebraic'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new TurbulenceClosure node (pyTree) <Examples/Converter/newTurbulenceClosurePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newTurbulenceClosurePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newTurbulenceModel(value='Null', parent=None) 

    Create a TurbulenceModel node. value is the model type in 'Algebraic_BaldwinLomax', 'Algebraic_CebeciSmith', 'HalfEquation_JohnsonKing', 'OneEquation_BaldwinBarth', 'OneEquation_SpalartAllmaras', 'TwoEquation_JonesLaunder', 'TwoEquation_MenterSST', 'TwoEquation_Wilcox'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new TurbulenceModel node (pyTree) <Examples/Converter/newTurbulenceModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newTurbulenceModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newThermalRelaxationModel(value='Null', parent=None)

    Create a ThermalRelaxationModel node. value is the model type in 'Frozen', 'ThermalEquilib', 'ThermalNonequilb'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ThermalRelaxationModel node (pyTree) <Examples/Converter/newThermalRelaxationModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newThermalRelaxationModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newChemicalKineticsModel(value='Null', parent=None)

    Create a ChemicalKineticsModel node. value is the model type in 'Frozen', 'ChemicalEquilibCurveFit', 'ChemicalEquilibMinimization', 'ChemicalNonequilib'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ChemicalKineticsModel node (pyTree) <Examples/Converter/newChemicalKineticsModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newChemicalKineticsModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newEMElectricFieldModel(value='Null', parent=None)

    Create a EMElectricFieldModel node. value is the model type in 'Constant', 'Frozen', 'Interpolated', 'Voltage'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new EMElectricFieldModel node (pyTree) <Examples/Converter/newEMElectricFieldModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newEMElectricFieldModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newEMMagneticFieldModel(value='Null', parent=None)

    Create a EMMagneticFieldModel node. value is the model type in 'Constant', 'Frozen', 'Interpolated'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new EMMagneticFieldModel node (pyTree) <Examples/Converter/newEMMagneticFieldModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newEMMagneticFieldModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newEMConductivityModel(value='Null', parent=None)

    Create a EMConductivityModel node. value is the model type in 'Constant', 'Frozen', 'Equilibrium_LinRessler', 'Chemistry_LinRessler'.
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new EMConductivityModel node (pyTree) <Examples/Converter/newEMConductivityModelPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newEMConductivityModelPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newBaseIterativeData(name='BaseIterativeData', parent=None)

    Create a BaseIterativeData node.
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new BaseIterativeData node (pyTree) <Examples/Converter/newBaseIterativeDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newBaseIterativeDataPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newZoneIterativeData(name='ZoneIterativeData', parent=None)

    Create a ZoneIterativeData node.
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ZoneIterativeData node (pyTree) <Examples/Converter/newZoneIterativeDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newZoneIterativeDataPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newRigidGridMotion(name='Motion', origin=[0.,0.,0.], mtype='Null', parent=None)

    Create a RigidGridMotion node. mtype is the motion type ('ConstantRate', 'VariableRate').
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param origin: coordinates of the origin point of the motion
    :type origin: list of 3 floats
    :param mtype: motion type
    :type mtype: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new RigidGridMotion node (pyTree) <Examples/Converter/newRigidGridMotionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newRigidGridMotionPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newRigidGridMotionType(value='ConstantRate', parent=None)

    Create a RigidGridMotionType node. value is the motion type ('ConstantRate', 'VariableRate').
    If parent is not None, attach it to parent node.

    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new RigidGridMotionType node (pyTree) <Examples/Converter/newRigidGridMotionTypePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newRigidGridMotionTypePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newReferenceState(name='ReferenceState', parent=None)

    Create a ReferenceState node.
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ReferenceState node (pyTree) <Examples/Converter/newReferenceStatePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newReferenceStatePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newConvergenceHistory(name='GlobalConvergenceHistory', value=0, parent=None)

    Create a ConvergenceHistory node. value is an iteration number.
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param value: value of the node (iteration number)
    :type value: integer
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ConvergenceHistory node (pyTree) <Examples/Converter/newConvergenceHistoryPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newConvergenceHistoryPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newFamily(name='Family', parent=None)

    Create a zone Family node. 
    If parent is not None, attach it to parent node.

    :param name: name of the family
    :type name: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new zone Family node (pyTree) <Examples/Converter/newFamilyPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newFamilyPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newFamilyBC(value='UserDefined', parent=None)

    Create a FamilyBC node. value is a BC type string.
    If parent is not None, attach it to parent node.

    :param value: BC type ('BCWall','UserDefined'...)
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new FamilyBC node (pyTree) <Examples/Converter/newFamilyBCPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newFamilyBCPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGeometryReference(value='Null', file='MyCAD.iges', parent=None)

    Create a GeometryReference node. value is a type of CAD ('NASA-IGES', 'SDRC', 'Unigraphics', 'ProEngineer', 'ICEM-CFD').
    If parent is not None, attach it to parent node.

    :param value: CAD type
    :type value: string
    :param file: name of CAD file
    :type file: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new GeometryReference node (pyTree) <Examples/Converter/newGeometryReferencePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGeometryReferencePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newArbitraryGridMotion(name='Motion', value='Null', parent=None)

    Create a ArbitraryGridMotion node. value is the type of motion ('NonDeformingGrid', 'DeformingGrid').
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param value: type of motion
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new ArbitraryGridMotion node (pyTree) <Examples/Converter/newArbitraryGridMotionPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newArbitraryGridMotionPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newUserDefinedData(name='UserDefined', value=None, parent=None)

    Create a UserDefinedData node to store user specific data.
    If parent is not None, attach it to parent node.

    :param name: name of the node
    :type name: string
    :param value: value of the node
    :type value: string
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new UserDefinedData node (pyTree) <Examples/Converter/newUserDefinedDataPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newUserDefinedDataPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Internal.newGravity(value=[0.,0.,9.81], parent=None)

    Create a Gravity node. value is the gravity vector.
    If parent is not None, attach it to parent node.

    :param value: gravity vector
    :type value: list of 3 floats
    :param parent: optional parent node
    :type parent: pyTree node
    :return: created node
    :rtype: pyTree node

    *Example of use:*

    * `Create a new Gravity node (pyTree) <Examples/Converter/newGravityPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/newGravityPT.py


.. toctree::
   :maxdepth: 2   

Index
######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

