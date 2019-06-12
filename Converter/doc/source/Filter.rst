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
   Converter.Filter.deletePaths

**-- High level layer**

.. autosummary::

    Converter.Filter.Handle
    Converter.Filter.Handle.loadSkeleton
    Converter.Filter.Handle.getVariables
    Converter.Filter.Handle.loadZonesWoVars
    Converter.Filter.Handle.loadVariables


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

    *Example of use:*

    * `Read skeleton tree (pyTree) <Examples/Converter/convertFile2SkeletonTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertFile2SkeletonTreePT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.readNodesFromPaths(fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1, skipTypes=None)

    Read nodes specified by their paths.
    If maxFloatSize=-1, all data are loaded, otherwise data are loaded
    only if the number of elements is lower that maxFloatSize.
    If maxDepth=-1, the read is fully recursive. Otherwise, load is limited
    to maxDepth levels.
    If skipTypes is specified, load is stopped when given node type is met (HDF only).

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
    :param skipTypes: list of CGNS types to skip
    :type skipTypes: None or list of strings
    :return: read nodes
    :rtype: pyTree node list

    *Example of use:*

    * `Read nodes from file (pyTree) <Examples/Converter/readNodesFromPathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/readNodesFromPathsPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.readNodesFromFilter(fileName, filter, format='bin_hdf', com=None)

    Partially read nodes specified by a filter.
    Filter is a dictionary for each path to be read.
    For structured grids: [[imin,jmin,kmin], [1,1,1], [imax,jmax,kmax], [1,1,1]]
    For unstructured grids: [[istart], [1], [iend], [1]].
    Only for HDFfile format.

    :param fileName: file name to read from
    :type fileName: string
    :param filter: paths and indices to be read 
    :type filter: dictionary of lists
    :param format: bin_hdf
    :type format: string
    :param com: communicator if run with mpi
    :type com: int
    :return: dictionary of read nodes
    :rtype: dictionary of numpys

    *Example of use:*

    * `Partially read nodes from file (pyTree) <Examples/Converter/readNodesFromFilterPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/readNodesFromFilterPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.readPyTreeFromPaths(t, fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1)

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
    :rtype: modified tree

    *Example of use:*

    * `Read nodes from file and modify tree (pyTree) <Examples/Converter/readPyTreeFromPathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/readPyTreeFromPathsPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.writeNodesFromPaths(fileName, paths, nodes, format=None, maxDepth=-1, mode=0)

    Write given nodes to specified paths in file.
    If mode=0 (append), nodes are appened to path location.
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


.. py:function:: Converter.Filter._readPyTreeFromPaths(a, fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1)

    Read nodes specified by their paths and stored it in a (in place).

    :param a: input data
    :type a: [pyTree, base, Zone, list of Zones]
    :param fileName: file name to read from
    :type fileName: string
    :param paths: paths to read (relative to a)
    :type paths: list of strings
    :param format: bin_cgns, bin_adf, bin_hdf (optional)
    :type format: string
    :param maxFloatSize: the maxSize of float array to load
    :type maxFloatSize: int
    :param maxDepth: max depth of load
    :type maxDepth: int

    *Example of use:*

    * `Read nodes from file and modify pyTree (pyTree) <Examples/Converter/readPyTreeFromPathsPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/readPyTreeFromPathsPT.py

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

.. py:function:: Converter.Filter.Handle.loadSkeleton()

    Load a skeleton tree from file (a tree where no data are loaded).

    :rtype: a skeleton pyTree     

    *Example of use:*

    * `Load a skeleton tree from file (pyTree) <Examples/Converter/loadSkeletonPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/loadSkeletonPT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.getVariables()

    Get the variables contained in file. This function minimal reads from file
    and store variable names in handle. This function must be called after loadSkeleton.

    :return: list of variables contained in file
    :rtype: list of strings

    *Example of use:*

    * `Read variable list from file (pyTree) <Examples/Converter/getVariablesPT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/getVariablesPT.py

---------------------------------------------------------------------------

.. py:function:: Converter.Filter.Handle.loadZonesWoVars(a, znp=None, bbox=None)

    Load specified zones (coordinates, grid connectivity, boundary conditions).
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

    Load specified variables.
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


.. toctree::
   :maxdepth: 2   

Index
######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
