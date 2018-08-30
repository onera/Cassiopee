.. Filter documentation master file

Filter: partial reading/writing of files
=========================================


Preamble
########

This module provides services for partially reading or writing cgns files (HDF or ADF).

To use the module::

    import Converter.Filter as Filter

List of functions
##################

**-- Low layer**

.. autosummary::

   Converter.Filter.convertFile2SkeletonTree
   Converter.Filter.readNodesFromPaths
   Converter.Filter.readNodesFromFilter
   Converter.Filter.writeNodesFromPaths
   Converter.Filter.deletePaths


**-- High layer**



Contents
#########

Per node reading/writing
------------------------

.. py:function:: Converter.Filter.convertFile2SkeletonTree(fileName, format=None, maxFloatSize=5, maxDepth=-1)

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
    :return: Skeleton tree
    :rtype: pyTree node

    *Example of use:*

    * `Read skeleton tree (pyTree) <Examples/Converter/convertFile2SkeletonTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Converter/convertFile2SkeletonTreePT.py


---------------------------------------------------------------------------

.. py:function:: Converter.Filter.readNodesFromPaths(fileName, paths, format=None, maxFloatSize=-1, maxDepth=-1)

    Read nodes specified by their paths.

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

.. py:function:: Converter.Filter.writeNodesFromPaths(fileName, paths, nodes, format=None, maxDepth=-1, mode=0)

    Write given nodes to specified paths in file.
    If mode=0 (append), nodes are appened to path location.
    If mode=1 (replace), nodes are replaced to path location. 
    If maxDepth>0, replace mode kill children of replaced node.
    If maxDepth=0, replace mode replace value and type of node (not the name).

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

    * `Write nodes from file (pyTree) <Examples/Converter/writeNodesFromPathsPT.py>`_:

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

.. toctree::
   :maxdepth: 2   

Index
######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
