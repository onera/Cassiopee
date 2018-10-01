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

.. toctree::
   :maxdepth: 2   


Index
#######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

