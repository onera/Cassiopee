.. Compressor documentation master file


Compressor: Field compression module
=====================================


Preamble
########

Compressor enables fields compression for arrays/pyTrees.

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

To use the module with the Compressor array interface::

   import Compressor

To use the module with the CGNS/Python interface::

    import Compressor.PyTree as Compressor


.. py:module:: Compressor

List of functions
#################

**-- Index field compression**

.. autosummary::

   Compressor.deltaIndex

**-- Object serializer/compression**

.. autosummary::

   Compressor.pack
   Compressor.unpack

**-- CGNS Zones/tree compression**

.. autosummary::

   Compressor.PyTree.compressCartesian
   Compressor.PyTree.uncompressCartesian
   Compressor.PyTree.compressCellN
   Compressor.PyTree.compressCoords
   Compressor.PyTree.compressFields
   Compressor.PyTree.compressElements
   Compressor.PyTree.uncompressAll

Contents
########

Index field compression
------------------------

.. py:function:: Compressor.deltaIndex(a, ref)

    Compress a list of indices using delta algorithm. 
    The return Delta contains
    the number of added indices in a when compared to ref, 
    the list of added indices, the
    number of suppressed indices, the list of suppressed indices.

    :param a: input indices
    :type a: numpy of ints
    :param ref: compared indices
    :type ref: numpy
    :return: list of added indices, the number of supressed indices, list of suppress indices
    :rtype: (numpy, int, numpy)

    * `Compression by delta (numpy) <Examples/Compressor/deltaIndex.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/deltaIndex.py


---------------------------------------

Object serialize/compression
-----------------------------

.. py:function:: Compressor.pack(a)

    Serialize/compress a python object a. For now, this is only a general interface
    to pickle module.

    :param a: any python object
    :type a: python object
    :return: serialized stream

    * `Object serialization (numpy) <Examples/Compressor/pack.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/pack.py
    

---------------------------------------

.. py:function:: Compressor.unpack(a)

    Deserialize/decompress a serialized stream b. For now, this is only a general interface
    to pickle module.

    :param a: a serialized stream as produced by pack
    :type a: serialized stream
    :return: python object

    * `Object deserialization (numpy) <Examples/Compressor/unpack.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/unpack.py


---------------------------------------


.. py:function:: Compressor.PyTree.compressCartesian(a)

    Compress zones if they are regular Cartesian grids. Create a
    CaresianData node containing the 6 floats
    corresponding to first point and steps in 3 directions.

    Exists also as an in-place version (_compressCartesian) which modifies a and returns None.

    :param a: input data
    :type a: [zone, list of zones, base, pyTree]
    :return: identical to input

    * `Cartesian compression (pyTree) <Examples/Compressor/compressCartesianPT.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/compressCartesianPT.py

---------------------------------------

.. py:function:: Compressor.PyTree.uncompressCartesian(a)

    Uncompress zones that has been compressed with compressCartesian.
    Exists also as an in-place version (_uncompressCartesian) which modifies a and returns None.

    :param a: input data
    :type a: [zone, list of zones, base, pyTree]
    :return: identical to input

    * `Cartesian compression (pyTree) <Examples/Compressor/uncompressCartesianPT.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/uncompressCartesianPT.py


---------------------------------------


.. py:function:: Compressor.PyTree.compressCellN(a)

    Compress cellN fields (valued 0,1,2).

    Exists also as an in-place version (_compressCellN) which modifies a and returns None.

    :param a: input data
    :type a: [zone, list of zones, base, pyTree]
    :return: identical to input

    * `CellN compression (pyTree) <Examples/Compressor/compressCellNPT.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/compressCellNPT.py


---------------------------------------


.. py:function:: Compressor.PyTree.compressCoords(a, tol=1.e-8, ctype=0)

    Compress zone coordinates with sz library with a relative tolerance tol.

    Exists also as an in-place version (_compressCoords) which modifies a and returns None.

    :param a: input data
    :type a: [zone, list of zones, base, pyTree]
    :param tol: control relative error on output
    :type tol: float
    :param ctype: compression algorithm
    :type ctype: 0 (sz), 1 (zfp)
    :return: identical to input

    * `Coordinates compression (pyTree) <Examples/Compressor/compressCoordsPT.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/compressCoordsPT.py


---------------------------------------

.. py:function:: Compressor.PyTree.compressFields(a, tol=1.e-8)

    Compress zone fields with sz library with a relative tolerance tol.

    Exists also as an in-place version (_compressFields) which modifies a and returns None.

    :param a: input data
    :type a: [zone, list of zones, base, pyTree]
    :param tol: control relative error on output
    :type tol: float
    :param ctype: compression algorithm
    :type ctype: 0 (sz), 1 (zfp)
    :return: identical to input

    * `Field compression (pyTree) <Examples/Compressor/compressFieldsPT.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/compressFieldsPT.py

---------------------------------------


.. py:function:: Compressor.PyTree.compressElements(a)

    Compress zone elements (connectivity).

    Exists also as an in-place version (_compressElements) which modifies a and returns None.

    :param a: input data
    :type a: [zone, list of zones, base, pyTree]
    :return: identical to input

    * `Element compression (pyTree) <Examples/Compressor/compressElementsPT.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/compressElementsPT.py


---------------------------------------

.. py:function:: Compressor.PyTree.uncompressAll(a, tol=1.e-8)

    Uncompress zones compressed with the previous compressors.

    Exists also as an in-place version (_uncompressAll) which modifies a and returns None.

    :param a: input data
    :type a: [zone, list of zones, base, pyTree]
    :return: identical to input

    * `Zone decompression (pyTree) <Examples/Compressor/uncompressAllPT.py>`_:

    .. literalinclude:: ../build/Examples/Compressor/uncompressAllPT.py

---------------------------------------

.. toctree::
   :maxdepth: 2   


Index
#######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

