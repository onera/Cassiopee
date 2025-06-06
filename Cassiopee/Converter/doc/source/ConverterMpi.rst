.. Converter.Mpi documentation master file

:tocdepth: 2


Converter.Mpi: distributed pyTree services
=============================================


Preamble
######## 

This module provides services to deal with distributed pyTrees.

A distributed pyTree is a tree where
zones are distributed over different processes.

The number of a process is called the rank. 

Three new concepts are introduced in addition to
standard pyTrees: the **skeleton tree**, the **loaded skeleton tree** and
the **partial tree**.

A **skeleton tree (S)** is a full pyTree where numpy arrays  
of DataArray_t type nodes are replaced by None.

A **loaded skeleton tree (LS)** is a skeleton tree for which zones 
attributed to the current rank are fully loaded.

A **partial tree (P)** is a pyTree with only zones attributed
to the current rank fully loaded and no skeleton zones. 
It can be viewed as a *loaded skeleton tree* with skeleton zones suppressed.

Generally, Cassiopee functions will operate seamlessly on
a partial tree, if the function doesn't requires data exchanges.
If the function requires exchanges, then a specific version exists
in the Module.Mpi module. For instance, for Converter, only center2Node
requires exchanges and is redefined in Converter.Mpi.

To use the module::

    import Converter.Mpi as Cmpi 

To run a python script in parallel with two processes::

    mpirun -np 2 python script.py 


.. py:module:: Converter.Mpi


List of functions
##################

**-- Input/output**

.. autosummary::
   :nosignatures:

   Converter.Mpi.convertFile2SkeletonTree
   Converter.Mpi.convertFile2PyTree
   Converter.Mpi.readZones
   Converter.Mpi.writeZones
   Converter.Mpi.convertPyTree2File

**-- Conversion**

.. autosummary::
   :nosignatures:

    Converter.Mpi.convert2PartialTree
    Converter.Mpi.convert2SkeletonTree
    Converter.Mpi.createBBoxTree

**-- Communication Graphs**

.. autosummary::
   :nosignatures:

    Converter.Mpi.getProc
    Converter.Mpi.setProc
    Converter.Mpi.getProcDict
    Converter.Mpi.computeGraph

**-- Exchanges**

.. autosummary::
   :nosignatures:

    Converter.Mpi.setCommunicator
    Converter.Mpi.addXZones
    Converter.Mpi.rmXZones
    Converter.Mpi.allgatherTree

**-- Actions**

.. autosummary::
   :nosignatures:

    Converter.Mpi.trace
    Converter.Mpi.center2Node


Contents
#########

Input/output
-------------

.. py:function:: Converter.Mpi.convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5, maxDepth=-1, links=None)

    Read a skeleton tree (**S**) from file (adf or hdf file format only). The loaded in
    memory skeleton tree is identical on all processors.

    If float data array size of DataArray_t type nodes is lower than maxFloatSize then the 
    array is loaded. Otherwise it is set to None. 
    If maxDepth is specified, load is limited to maxDepth levels. 

    :param fileName: file name to read from
    :type fileName: string
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param maxFloatSize: the maxSize of float array to load
    :type maxFloatSize: int
    :param maxDepth: max depth of load
    :type maxDepth: int
    :param links: if not None, return a list of links in file
    :type links: list of list of 4 strings
    :return: Skeleton tree
    :rtype: pyTree node

    *Example of use:*

    * `Read skeleton tree (pyTree) <Examples/Converter/convertFile2SkeletonTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertFile2SkeletonTreePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.convertFile2PyTree(fileName, format=None, proc=None)

    If proc=None, read a full tree. The fully loaded in
    memory tree is identical on all processors.
    
    If proc is given, read a partial tree of zones corresponding to proc.

    :param fileName: file name to read from
    :type fileName: string
    :param format: any converter format (optional)
    :type format: string
    :param proc: None or rank number 
    :type proc: None or int
    :return: fully loaded or partial tree
    :rtype: pyTree node

    *Example of use:*

    * `Read full tree (pyTree) <Examples/Converter/convertFile2PyTreeMPIPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertFile2PyTreeMPIPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.readZones(t, fileName, format=None, rank=None, zoneNames=None)

    Fill the data of skeleton zones of t following zone rank or zone name (adf or hdf).

    If rank is not None, zone must have been attributed to ranks either with
    Distributor2.distribute or setProc. If the zone rank corresponds to 
    process rank, the zone is filled with data read from file.

    If zoneNames is not None, zone with corresponding name are filled with data read from file.

    Exists also as in place version (_readZones) that modifies t
    and returns None.

    :param t: input data
    :type t: [pyTree]
    :param fileName: file name to read from
    :type fileName: string
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param rank: the processor of zones to read
    :type rank: int
    :param zoneNames: paths of zones to read (if rank is not set)
    :type zoneNames: list of strings
    :return: modified reference copy of t
    :rtype: Identical to t

    *Example of use:*

    * `Read some zones in a skeleton tree (pyTree) <Examples/Converter/readZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/readZonesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.writeZones(t, fileName, format=None, rank=None, zoneNames=None)

    Write some zones in an existing file (adf or hdf) according to zone rank or zone name.
    If by rank, zone must have been attributed to processors either with
    Distributor2.distribute or setProc. 

    :param t: input data
    :type t: [pyTree]
    :param fileName: file name to write to
    :type fileName: string
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param rank: the processor of written zones
    :type rank: int
    :param zoneNames: paths of written zones (if rank is not set)
    :type zoneNames: list of strings
    
    *Example of use:*

    * `Write some zones in an existing file (pyTree) <Examples/Converter/writeZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/writeZonesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.convertPyTree2File(t, fileName, format=None, links=[], ignoreProcNodes=False)

   Write a skeleton tree (**S**), a loaded skeleton tree (**LS**) or a 
   partial tree (**P**) to a file (adf or hdf).

   :param t: input data
   :type t: [pyTree]
   :param fileName: file name to write to
   :type fileName: string
   :param format: bin_cgns, bin_adf, bin_hdf (optional)
   :type format: string
   :param links: optional list of links to be written
   :type links: list of list of 4 strings
   :param ignoreProcNodes: if true, only write zones with procNode set to rank, else write all proc zones
   :type ignoreProcNodes: boolean
    
   *Example of use:*

   * `Write a tree to a file (pyTree) <Examples/Converter/convertPyTree2FileMPI.py>`_:

   .. literalinclude:: ../build/Examples/Converter/convertPyTree2FileMPI.py


---------------------------------------------------------------------------

Conversions
----------------

.. py:function:: Converter.Mpi.convert2SkeletonTree(t, maxSize=6)

    Convert a tree (**LS** or **P**) to a skeleton tree (S). In a skeleton tree,
    numpys in DataArray_t nodes are replaced by None.

    Exists also as in place version (_convert2SkeletonTree) that modifies t
    and returns None.

    :param t: input data
    :type t: [pyTree, base, zone, list of zones]
    :param maxSize: numpy bigger than maxSize are remplace by None
    :type maxSize: int
    :rtype: Identical to t

    *Example of use:*

    * `Convert a tree to a skeleton tree (pyTree) <Examples/Converter/convert2SkeletonTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convert2SkeletonTreePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.convert2PartialTree(t, rank=-1)

    Convert a loaded skeleton tree (**LS**) to a partial tree (**P**). 
    If rank=-1, all skeleton zones are suppressed.
    If rank>=0, zones with proc != rank are suppressed.

    Exists also as in place version (_convert2PartialTree) that modifies t
    and returns None.

    :param t: input data
    :type t: [pyTree, base, zone, list of zones]
    :param rank: if rank=-1: suppress all skeleton zones, if rank>=0: suppress zones with proc=rank 
    :rtype: Identical to t

    *Example of use:*

    * `Convert a tree to a partial tree (pyTree) <Examples/Converter/convert2PartialTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convert2PartialTreePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.createBBoxTree(t, method='AABB')

    From a partial tree (**P**) or a loaded skeleton tree (**LS**), create a full tree containing 
    the bbox of zones. A bbox is a structured grid made of 8 points englobing
    zone. The returned tree is identical on all processors. 
    Argument method can be 'AABB' (axis aligned bbox) or 'OBB' (oriented bbox). 

    :param t: input data
    :type t: [pyTree, base, zone, list of zones]
    :param method: 'AABB': axis aligned bbox, 'OBB': oriented bbox
    :rtype: Identical to t

    *Example of use:*

    * `Create a BBox tree (pyTree) <Examples/Converter/createBBoxTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/createBBoxTreePT.py

---------------------------------------------------------------------------

Graphs
------------


.. py:function:: Converter.Mpi.getProc(z)

    Get the rank of zone. 
    It only returns the value stored in .Solver#Param/proc node.
    You can use setProc or Distributor2 to create the .Solver#Param/proc node.

    :param z: input zone
    :type z: [zone]
    :rtype: int

    *Example of use:*

    * `Get the proc of a zone (pyTree) <Examples/Converter/getProcPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getProcPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.setProc(t, rank)

    Set rank in t. t can be a skeleton (**S**), partial (**P**) or full tree.
    It only creates a .Solver#Param/proc node for the zones of t.

    Exists also as in place version (_setProc) that modifies t and returns None.

    :param t: input data
    :type t: [pyTree, base, zone, list of zones]
    :param rank: the rank value to set
    :type rank: int
    :return: modified reference copy of t
    :rtype: Identical to t

    *Example of use:*

    * `Set the rank in zones (pyTree) <Examples/Converter/setProcPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/setProcPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.getProcDict(t)

    Return the rank information stored in .Solver#Param/proc as 
    a dictionary proc['zoneName'].

    :param t: input data
    :type t: [pyTree, base, zone, list of zones]
    :rtype: Dictionary

    *Example of use:*

    * `Get the dictionary of zone procs (pyTree) <Examples/Converter/getProcDictPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getProcDictPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.computeGraph(t, type='bbox', t2=None, procDict=None, rank=0, intersectionsDict=None)

    Compute a communication graph. The graph is a dictionary such that graph[proc1][proc2] contains the names of zones of proc1 that are "connected" to at least one zone on proc2.

    - If type='bbox', a zone is connected to another if their bbox intersects. A must be a bbox tree.
    - If type='bbox2', a zone is connected to another if their bbox intersects and are not in the same base. A must be a bbox tree.
    - If type='bbox3', a zone is connected to another of another tree if their bbox intersects. A and t2 must be a bbox tree.
    - If type='match' (S/LS/P), a zone is connected to another if they have a match between them. A can be a skeleton, loaded skeleton or a partial tree.
    - If type='ID' (S/LS/P), a zone is connected to another if they have interpolation data between them. A can be a skeleton, a loaded skeleton or a partial tree.
    - If type='IBCD' (S/LS/P), a zone is connected to another if they have IBC data between them. A can be a skeleton, a loaded skeleton or a partial tree.
    - If type='ALLD' (S/LS/P), a zone is connected to another if they have Interpolation or IBC data between them. A can be a skeleton, a loaded skeleton or a partial tree.
    - If type='proc', a zone is attributed to another proc than the one it is loaded on. A can be a skeleton, a loaded skeleton tree or a partial tree.
    - If type='POST', t defines the donor tree, where the interpolation data is stored and t2 the receptor tree, as they do not define the same zones. Requires procDict and procDict2

    :param t: input data
    :type t: [pyTree, base, zone, list of zones]
    :param type: type of graph
    :type type: string in 'bbox', 'bbox2', 'bbox3', 'match', 'ID', 'IBCD', 'ALLD', 'proc'
    :param t2: optional second tree for type='bbox3'
    :type t2: pyTree
    :param procDict: if provided, used for zone affectation
    :type procDict: dictionary
    :param intersectionDict: dictionary of intersections
    :type intersectionDict: python dictionary
    :rtype: python dictionary of communications

    *Example of use:*

    * `Return the communication graph (pyTree) <Examples/Converter/computeGraphPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/computeGraphPT.py

---------------------------------------------------------------------------

Exchanges
------------

.. py:function:: Converter.Mpi.setCommunicator(com)

    Set the MPI communicator for Cassiopee exchanges. By default,
    it is set to MPI.COMM_WORLD.

    :param com: communicator to set
    :type com: MPI communicator


---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.addXZones(t, graph)

    For a partial tree, add zones loaded on a different process that are 
    connected to local zones through the graph.

    Exists also as in place version (_addXZones) that modifies t and returns None.

    :param t: input tree
    :type t: [pyTree]
    :param graph: communication graph as defined by computeGraph
    :type graph: dictionary
    :return: modified reference copy of t
    :rtype: tree with connected zones added

    *Example of use:*

    * `Add connected zones to a partial tree (pyTree) <Examples/Converter/addXZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/addXZonesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.rmXZones(t)

    For a partial tree, remove zones created by addXZones.

    Exists also as in place version (_rmXZones) that modifies t and returns None.

    :param t: input tree
    :type t: [pyTree]
    :rtype: tree with connected zones suppressed

    *Example of use:*

    * `Remove connected zones from a partial tree (pyTree) <Examples/Converter/rmXZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/rmXZonesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.allgatherTree(t)

    Gather a distributed tree on all processors.
    All processors then see the same tree.

    :param t: input tree
    :type t: [pyTree]
    :rtype: merged and gathered tree (identical on all processors)

    *Example of use:*

    * `Gather a distributed tree (pyTree) <Examples/Converter/allgatherTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/allgatherTreePT.py


---------------------------------------------------------------------------


Actions
-------------

.. py:function:: Converter.Mpi.trace(text, cpu=True, mem=True, fileName="stdout", method=0)

    Enable to monitor CPU usage and memory usage for each node/process.
    If fileName="stdout", information is written to standard output with the processor number. 
    If fileName is another name, information is written in "filenameXXX" files, one for each process. If filename does not have a format for the process number and/or an file extension, it is appended "%03d" and ".out", respectively.
    If cpu=True, time elapsed since the previous call to "trace" by this node is written.
    If mem=True, the current usage of memory of each process is written.
    Method used to find memory: 0: Rss in smaps, 1: tracemalloc.

    :param text: text to write
    :type text: string
    :param cpu: True to write cpu usage information
    :type cpu: boolean
    :param mem: True to write memory usage information
    :type mem: boolean
    :param fileName: name of the file to write information to or "stdout"
    :type fileName: string
    :param method: method used to check memory
    :type method: 0: Rss, 1: tracemalloc

    *Example of use:*

    * `Write a trace (pyTree) <Examples/Converter/tracePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/tracePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Mpi.center2Node(t, var=None, cellNType=0, graph=None)

    Perform a center to node conversion for a distributed tree.
    If var is set, var can be a field name or a container name 
    (Internal.__FlowSolutionNodes__, Internal.__FlowSolutionCenters__, ...).
    Then, center2Node is performed in the given field.
    Otherwise, the zone with its coordinates is moved to node.

    :param t: input tree
    :type t: [pyTree]
    :rtype: merged and gathered tree (identical on all processors)

    *Example of use:*

    * `Perform center2Node for a distributed tree (pyTree) <Examples/Converter/center2NodeMpiPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/center2NodeMpiPT.py


---------------------------------------------------------------------------


.. toctree::
   :maxdepth: 2   

Index
###################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

