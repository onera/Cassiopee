.. Filter documentation master file

Filter: partial reading/writing of files
=========================================


Preamble
########

This module provides services for partially reading or writing cgns files (HDF or ADF).

To use the module::

    import Converter.Filter as Filter

.. py:module:: Converter.Filter

List of functions
##################

**-- Low level layer**

.. autosummary::

   Converter.Filter.convertFile2SkeletonTree
   Converter.Filter.readNodesFromPaths
   Converter.Filter.readNodesFromFilter
   Converter.Filter.readPyTreeFromPaths
   Converter.Filter.writeNodesFromPaths
   Converter.Filter.writePyTreeFromPaths
   Converter.Filter.writePyTreeFromFilter
   Converter.Filter.deletePaths

**-- High level layer**

.. autosummary::

    Converter.Filter.Handle
    Converter.Filter.Handle.loadSkeleton
    Converter.Filter.Handle.getVariables
    Converter.Filter.Handle.loadZones
    Converter.Filter.Handle.loadZonesWoVars
    Converter.Filter.Handle.loadVariables
    Converter.Filter.Handle.writeZones
    Converter.Filter.Handle.writeZonesWoVars
    Converter.Filter.Handle.writeVariables
    Converter.Filter.Handle.loadFromProc
    Converter.Filter.Handle.loadAndDistribute
    Converter.Filter.Handle.loadAndSplit
    

Contents
#########

Low level layer
----------------

.. py:function:: Converter.Filter.convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5, maxDepth=-1, links=None)

    Read a skeleton tree. 
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

    For documentation on links, see Converter read options. 

    *Example of use:*

    * `Read skeleton tree (pyTree) <Examples/Converter/convertFile2SkeletonTreeFPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertFile2SkeletonTreeFPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.readNodesFromPaths(fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1, dataShape=None, skipTypes=None, com=None)

    Read nodes specified by their paths.
    If maxFloatSize=-1, all data are loaded, otherwise data are loaded
    only if the number of elements is lower that maxFloatSize.
    If maxDepth=-1, the read is fully recursive. Otherwise, load is limited
    to maxDepth levels.
    If skipTypes is specified, children of nodes of given type are not loaded (HDF only).

    :param fileName: file name to read from
    :type fileName: string
    :param paths: paths to read
    :type paths: list of strings
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param maxFloatSize: the maxSize of float array to load
    :type maxFloatSize: int
    :param maxDepth: max depth of load
    :type maxDepth: int
    :param dataShape: dictionary of returned data shapes if not None
    :type dataShape: None or dictionary of shapes of data
    :param skipTypes: list of CGNS types to skip
    :type skipTypes: None or list of strings
    :param com: optional MPI communicator. If set, triggers parallel IO
    :type com: MPI communicator
    :return: read nodes
    :rtype: pyTree node list

    *Example of use:*

    * `Read nodes from file (pyTree) <Examples/Converter/readNodesFromPathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/readNodesFromPathsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.readNodesFromFilter(fileName, filter, format='bin_hdf', com=None)

    Partially read nodes from a file specified by a filter.
    Filter is a dictionary for each path to be read.
    It specifies the slice of array you want to load from file and write to memory.
    Each slice information is made of [ offset, stride, count, block ].
    Offset specifies the starting index of slice.
    Stride specifies the stride of slice (for example 1 means all elements, 2 means one over two, ...).
    Count is the number of element to read.
    Block is always 1.

    For example, for structured grids: [[imin,jmin,kmin], [1,1,1], [imax-imin+1,jmax-jmin+1,kmax-kmin+1], [1,1,1]].

    For example, for unstructured grids: [[istart], [1], [iend-imax+1], [1]].
    
    Only for HDFfile format.

    :param fileName: file name to read from
    :type fileName: string
    :param filter: paths and indices to be read 
    :type filter: dictionary of lists
    :param format: bin_hdf
    :type format: string
    :param com: optional MPI communicator. If set, tirggers parralel IO
    :type com: MPI communicator
    :return: dictionary of read node data
    :rtype: dictionary of numpys

    *Example of use:*

    * `Partially read nodes from file (pyTree) <Examples/Converter/readNodesFromFilterPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/readNodesFromFilterPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.readPyTreeFromPaths(t, fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1, dataShape=None, skipTypes=None, com=None)

    Read nodes of t specified by their paths.
    Exists also as in place function (_readPyTreeFromPaths) that modifies t
    and returns None.

    :param t: input tree
    :type t: pyTree
    :param fileName: file name to read from
    :type fileName: string
    :param paths: paths to read
    :type paths: list of strings
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param maxFloatSize: the maxSize of float array to load
    :type maxFloatSize: int
    :param maxDepth: max depth of load
    :type maxDepth: int
    :param dataShape: dictionary of returned data shapes if not None
    :type dataShape: None or dictionary of shapes of data
    :param skipTypes: list of CGNS types to skip
    :type skipTypes: None or list of strings
    :param com: optional MPI communicator. If set, triggers parallel IO
    :type com: MPI communicator
    :rtype: modified tree

    *Example of use:*

    * `Read nodes from file and modify tree (pyTree) <Examples/Converter/readPyTreeFromPathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/readPyTreeFromPathsPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.writeNodesFromPaths(fileName, paths, nodes, format=None, maxDepth=-1, mode=0)

    Write given nodes to specified paths in file.
    If mode=0 (append), nodes are appened to path location. Nevertheless, if a node with identical name already 
    exists in the path node children, it will be replaced by the appended ones.
    If mode=1 (replace), nodes are replaced to path location. 
    If maxDepth>0, replace mode kill children of replaced node.
    If maxDepth=0, replace mode replaces value and type of node (not the name).

    :param fileName: file name to write to
    :type fileName: string
    :param paths: paths to write to
    :type paths: list of strings
    :param nodes: nodes to write
    :type nodes: list of pyTree nodes
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param maxDepth: max depth to write
    :type maxDepth: int
    :param mode: writing mode (0: append, 1: replace)
    :type mode: int

    *Example of use:*

    * `Write nodes from paths (pyTree) <Examples/Converter/writeNodesFromPathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/writeNodesFromPathsPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.writePyTreeFromPaths(t, fileName, paths, format=None, maxDepth=-1)

    Write given paths of tree to the same specified paths in file.

    :param t: input pyTree
    :type t: pyTree node
    :param fileName: file name to write to
    :type fileName: string
    :param paths: paths to write to
    :type paths: list of strings
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param maxDepth: max depth to write
    :type maxDepth: int

    *Example of use:*

    * `Write pyTree from paths (pyTree) <Examples/Converter/writePyTreeFromPathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/writePyTreeFromPathsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.writePyTreeFromFilter(t, fileName, filter, format='bin_hdf', com=None)

    Write partial data from tree t to file specified by a filter.

    :param t: input pyTree
    :type t: pyTree node
    :param fileName: file name to write to
    :type fileName: string
    :param filter: paths and indices to be read 
    :type filter: dictionary of lists
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param com: optional MPI communicator. If set, triggers parallel IO
    :type com: MPI communicator
    
    *Example of use:*

    * `Write pyTree from filter (pyTree) <Examples/Converter/writePyTreeFromFilterPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/writePyTreeFromFilterPT.py

---------------------------------------------------------------------------


.. py:function:: Converter.Filter.deletePaths(fileName, paths, format=None)

    Delete paths in file.

    :param fileName: file name to read from
    :type fileName: string
    :param paths: paths to read (relative to a)
    :type paths: list of strings
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string

    *Example of use:*

    * `Delete some nodes in file from paths (pyTree) <Examples/Converter/deletePathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/deletePathsPT.py


High level layer
----------------

.. py:function:: Converter.Filter.Handle(fileName)

    Create a handle on a file to enable partial reading.
    The file must be a CGNS/ADF or CGNS/HDF file.

    :param fileName: file name to read from
    :type fileName: string
    :rtype: handle class 
    
    *Example of use:*

    * `Create a handle on a file (pyTree) <Examples/Converter/handlePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/handlePT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.loadSkeleton(maxDepth=3, readProcNode=False)

    Load a skeleton tree from file (a tree of depth maxDepth where no data are loaded).

    :param maxDepth: the depth you want to load (-1 means all tree)
    :type maxDepth: int
    :param readProcNode: if true, force reading of proc node. Usefull for maxDepth <= 4.
    :type readProcNode: boolean
    :rtype: a skeleton pyTree

    *Example of use:*

    * `Load a skeleton tree from file (pyTree) <Examples/Converter/loadSkeletonPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/loadSkeletonPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.getVariables(cont=None)

    Get the names of variables contained in file. 
    This function must be called after loadSkeleton.

    :param cont: container name. Can be a CGNS name 'FlowSolution', ... or 'centers' or 'nodes'
    :return: list of variable names contained in file
    :rtype: list of strings

    *Example of use:*

    * `Read variable list from file (pyTree) <Examples/Converter/getVariablesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getVariablesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.loadZones(a, znp=None)

    Fully load specified zones (coordinates, fields, grid connectivity, boundary conditions) in tree.
    This function must be called after loadSkeleton.

    :param a: modified pyTree 
    :type a: pyTree
    :param znp: paths of zones to load (must be a list of 'BaseName/ZoneName')
    :type znp: list of strings
    
    *Example of use:*

    * `Fully load zones (pyTree) <Examples/Converter/loadZonesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/loadZonesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.loadZonesWoVars(a, znp=None, bbox=None)

    Load specified zones (coordinates, grid connectivity, boundary conditions) in tree.
    If bbox=[xmin,ymin,zmin,xmax,ymax,zmax] is specified, load only zones
    intersecting this bbox.
    This function must be called after loadSkeleton.

    :param a: modified pyTree 
    :type a: pyTree
    :param znp: path of zones to load from (starting from top)
    :type znp: list of strings
    :param bbox: optional bbox
    :type bbox: list of 6 floats

    *Example of use:*

    * `Load zones without variable (pyTree) <Examples/Converter/loadZonesWoVarsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/loadZonesWoVarsPT.py

---------------------------------------------------------------------------
        
.. py:function:: Converter.Filter.Handle.loadVariables(a, var, znp=None)

    Load specified variables in tree.
    This function must be called after loadSkeleton.

    :param a: modified pyTree 
    :type a: pyTree
    :param var: variables to load
    :type var: string or list of strings
    :param znp: path of zones to load from (starting from top)
    :type znp: list of strings

    *Example of use:*

    * `Load given variables (pyTree) <Examples/Converter/loadVariablesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/loadVariablesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.writeZones(a, fileName=None, znp=None)

    Fully write specified zones (coordinates, fields, grid connectivity, boundary conditions) in file.
    
    :param a: input pyTree 
    :type a: pyTree or list of zones
    :param fileName: file name if different of handle name
    :type fileName: string
    :param znp: path of zones to write to in file (starting from root)
    :type znp: list of strings
    
    *Example of use:*

    * `Fully write zones (pyTree) <Examples/Converter/writeZonesFPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/writeZonesFPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.writeZonesWoVars(a, fileName=None, znp=None)

    Write specified zones without fields (coordinates, grid connectivity, boundary conditions) in file.
    
    :param a: input pyTree 
    :type a: pyTree or list of zones
    :param fileName: file name if different of handle name
    :type fileName: string
    :param znp: path of zones to write to in file (starting from root)
    :type znp: list of strings
    
    *Example of use:*

    * `Write zones without fields (pyTree) <Examples/Converter/writeZonesWoVarsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/writeZonesWoVarsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.writeVariables(a, var, fileName=None, znp=None)

    Write specified variables in file.
    
    :param a: input pyTree 
    :type a: pyTree or list of zones
    :param var: variables to write
    :type var: string or list of strings
    :param fileName: file name if different of handle name
    :type fileName: string
    :param znp: path of zones to write to in file (starting from root)
    :type znp: list of strings
    
    *Example of use:*

    * `Write variables (pyTree) <Examples/Converter/writeVariablesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/writeVariablesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.loadFromProc(loadVariables=True)

    Load on each processor the zones with the corresponding proc node.
    The zones in file must have a .Solver#Param/proc node.

    :param loadVariables: If true, load all variables in file. Otherwise load only coordinates
    :type loadVariables: Boolean
    :rtype: partial tree on each processor

    *Example of use:*

    * `Load tree from proc node (pyTree) <Examples/Converter/loadFromProcPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/loadFromProcPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.loadAndDistribute(strategy=None, algorithm='graph', loadVariables=True)

    Load and distribute zones of file on the different processors.
    
    :param strategy: strategy for distribution. Can be None (only the number of points of block is considered), 'match' (use matching boundaries to optimize distribution)
    :type  strategy: string
    :param algorithm: algorithm for distribution. Can be 'graph', 'fast', 'gradient'. See Distributor2 documentation.
    :type algorithm: string
    :param loadVariables: If true, load all variables in file. Otherwise load only coordinates
    :type loadVariables: Boolean
    :rtype: partial tree on each processor

    *Example of use:*

    * `Load and distribute tree (pyTree) <Examples/Converter/loadAndDistributePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/loadAndDistributePT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.loadAndSplit()

    Load and split zones of file on the different processors (only for structured zones).
    
    :rtype: partial tree on each processor

    *Example of use:*

    * `Load and split tree (pyTree) <Examples/Converter/loadAndSplitPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/loadAndSplitPT.py



.. toctree::
   :maxdepth: 2   

Index
######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
