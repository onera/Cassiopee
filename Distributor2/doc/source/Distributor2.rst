.. Distributor2 documentation master file

Distributor2: distribution of grids on processors
=================================================

Preamble
########

This module provides functions to distribute blocks on a given 
number of processors. At the end of the process, each block
will have a number corresponding to the processor it must be affected
to for a balanced computation, depending on given criterias. This
module doesn't perform splitting (see the Transform module for
that).

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use with the array interface, you have to import Distributor2 module::

   import Distributor2 as D2

For use with the pyTree interface::

    import Distributor2.PyTree as D2


.. py:module:: Distributor2


List of functions
##################

**-- Automatic load balance**

.. autosummary::

   Distributor2.distribute
   Distributor2.PyTree.distribute

**-- Various operations**

.. autosummary::

   Distributor2.PyTree.addProcNode
   Distributor2.PyTree.getProc
   Distributor2.PyTree.getProcDict
   Distributor2.PyTree.getProcList
   Distributor2.PyTree.copyDistribution
   Distributor2.PyTree.printProcStats
   Distributor2.Mpi.redispatch

   


Contents
#########

.. py:function:: Distributor2.distribute(A, NProc, prescribed=[], perfo=[], weight=[], com=[], algorithm='graph', mode='nodes', nghost=0)

    Distribute automatically the blocks amongst NProc processors.

    - prescribed is a list of blocks that are forced to be on a given processor.

    For instance, prescribed[2] = 0 means that block 2 MUST be affected to processor 0.

    - perfo is a tuple or a tuple list for each processor. 

    Each tuple describes the relative weight of solver CPU time regarding the communication speed and latence (solverWeight, latenceWeight, comSpeedWeight).

    - weight is a list of weight for each block indicating the relative cost for solving each block.
    
    - com is a ixj matrix describing the volume of points exchanged between bloc i and bloc j.

    - algorithm can be chosen in: 'gradient', 'genetic', 'fast', 'graph'

    - mode='node', 'cells': optimize distribution of block considering node (cells) numbers.

    - nghost: take into account ghost cells (only for structured grids)

    :param a: Input data
    :type  a: [array, list of arrays]
    :param N: number of processors
    :type N: int
    :param prescribed: list of prescribed blocks
    :type prescribed: list of ints
    :param perfo: list of performance for each processor
    :type perfo: list of tuples
    :param weight: list of weight for each block
    :type weight: list of ints
    :param algorithm: ['gradient', 'genetic', 'fast', 'graph']
    :type algorithm: string
    :param nghost: number of ghost cells present in the mesh
    :type nghost: int

    The function output is a stats dictionary.
    stat['distrib'] is a vector describing the attributed processor for each 
    block, stats['meanPtsPerProc'] is the mean number of points per proc,
    stats['varMin'] is the minimum variation of number of points,
    stats['varMax'] is the maximum variation of number of points,
    stats['varRMS'] is the mean variation of number of points,
    stats['nptsCom'] is the number of points exchanged between processors for
    communication, stats['comRatio'] is the ratio between the number of points exchanged between 
    processors in this configuration divided by the total number of matching/overlap boundary points,
    stats['adaptation'] is the value of the optimized function.

    
    
    *Example of use:*

    * `Distribute arrays (array) <Examples/Distributor2/distribute.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/distribute.py

    

===================================================================

.. py:function:: Distributor2.PyTree.distribute(A, NProc, prescribed={}, perfo=[], weight=[], useCom='all', algorithm='graph', mode='nodes', nghost=0)

    Distribute automatically the blocks amongst NProc processors.

    With the pyTree interface, the user-defined node .Solver#Param/proc 
    is updated with the attributed processor number.

    If useCom=0, only the grid number of points is taken into account.
    If useCom='all', matching and overlap communications are taken into account. 
    If useCom='match', only match connectivity are taken into account.
    if useCom='overlap', only overlap connectivity are taken into account.
    if useCom='bbox', overlap between zone bbox is taken into account.
    if useCom='ID', ID (interpolation or match) and IBCD (IBM points) are taken into account.

    When using distributed trees, prescribed must be a dictionary containing 
    the zones names as key, and the prescribed proc as value.  
    weight is also a dictionary where the keys are the zone names and the weight as the value.
    It is not mandatory to assign a weight to all the zones of the pyTree. Default value is assumed 1,
    only different weight values can be assigned to zones. 
    t can be either a skeleton or a loaded skeleton pyTree for useCom=0 or useCom='match', 
    but must be a loaded skeleton tree only for the other settings.
    
    :param a: Input data
    :type  a: [pyTree, base, zone, list of zones]
    :param N: number of processors
    :type N: int
    :param prescribed: dictionary of prescribed block (optional)
    :type prescribed: dictionary
    :param perfo: list of perfo for each processor (optional)
    :type perfo: list of tuples
    :param weight: dictionary of weights for block (optional)
    :type weight: dictionary
    :param useCom: tell what to use to measure communication volumes
    :type useCom: ['0, 'all', 'match', 'overlap', 'bbox', 'ID']
    :param algorithm: ['gradient', 'genetic', 'fast', 'graph']
    :type algorithm: string
    :param nghost: number of ghost cells present in the mesh
    :type nghost: int

    *Example of use:*

    * `Distribute zones (pyTree) <Examples/Distributor2/distributePT.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/distributePT.py

    
Various operations
---------------------

.. py:function:: Distributor2.PyTree.addProcNode(a, NProc)

    Add a "proc" node to all zones of a with given value.
    Exists also as in place version (_addProcNode) that modifies
    a and returns None.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param NProc: proc to be set
    :type NProc: int
    :return: reference copy of a
    :rtype: identical to input

    *Example of use:*

    * `Add a proc node (pyTree) <Examples/Distributor2/addProcNodePT.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/addProcNodePT.py
    
----------------------------------------------------------------

.. py:function:: Distributor2.PyTree.getProc(a)

    Return the proc value of a zone or a list of zones.
    
    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :return: the affected proc of zone
    :rtype: int or list of ints (for multiple zones)
    
    *Example of use:*

    * `Get proc of a zone (pyTree) <Examples/Distributor2/getProcPT.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/getProcPT.py

----------------------------------------------------------------

.. py:function:: Distributor2.PyTree.getProcDict(a, prefixByBase=False)

    Return a dictionary where procDict['zoneName'] is the no
    of proc affected to zone 'zoneName'.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param prefixByBase: if true, add base prefix to zone name
    :type prefixByBase: boolean
    :return: the dictionary of zone/proc.
    :rtype: dictionary
    
    *Example of use:*

    * `Return the dictionary of zones/proc (pyTree) <Examples/Distributor2/getProcDictPT.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/getProcDictPT.py

----------------------------------------------------------------

.. py:function:: Distributor2.PyTree.getProcList(a, NProc=None)

    Return procList where procList[proc] is a list of zone names
    attributed to the proc processor.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :return: the affected proc of zone
    :rtype: int or list of ints
    
    *Example of use:*

    * `Return the list of zones affected to a proc (pyTree) <Examples/Distributor2/getProcListPT.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/getProcListPT.py    

----------------------------------------------------------------

    
.. py:function:: Distributor2.PyTree.copyDistribution(a, b)

    Copy the distribution from b to a matching zones by their name.
    Exists also as in place version (_copyDistribution) that modifies
    a and returns None.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param b: original data
    :type b: [pyTree, base, zone, list of zones]
    :return: modifie reference copy of a
    :rtype: same as input data

    *Example of use:*

    * `Copy proc distribution from one tree to another (pyTree) <Examples/Distributor2/copyDistributionPT.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/copyDistributionPT.py    

    

----------------------------------------------------------------

.. py:function:: Distributor2.redispatch(a)

    Redispatch a tree where a new distribution is defined in the node 'proc'.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :return: modifie reference copy of a
    :rtype: same as input data

    *Example of use:*

    * `Redispatch a tree (pyTree) <Examples/Distributor2/redispatchPT.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/redispatchPT.py    
    

----------------------------------------------------------------

.. py:function:: Distributor2.PyTree.printProcStats(a, stats=None, NProc=None)

    Print statistics for each processor: number of points and list of zones names.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param stats: dictionary obtained from Distributor2.distribute 
    :type stats: Python dictionary
    :param NProc: number of processors 
    :type NProc: integer
    :return: None

    *Example of use:*

    * `Print procs statistics after distribution of a tree (pyTree) <Examples/Distributor2/printProcStatsPT.py>`_:

    .. literalinclude:: ../build/Examples/Distributor2/printProcStatsPT.py    
    
    .. note:: new in version 2.7.



.. toctree::
   :maxdepth: 2   


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

