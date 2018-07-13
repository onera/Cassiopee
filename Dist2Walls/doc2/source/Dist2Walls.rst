.. Dist2Walls documentation master file


Dist2Walls: wall distance computation
=========================================

Preamble
########

Dist2Walls gathers efficient algorithms for computing the distance fields 
for arrays (as defined in Converter documentation) or 
for CGNS/python tree (pyTrees).

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use with the array interface, you have to import Dist2Walls module::

   import Dist2Walls

For use with the pyTree interface::

    import Dist2Walls.PyTree


.. py:module:: Dist2Walls

List of functions
##################

**-- Wall distance computation**

.. autosummary::

    Dist2Walls.distance2Walls


Contents
#########

Wall distance computation
--------------------------

.. py:function:: Dist2Walls.distance2Walls(a, bodies, type='ortho', loc='centers', signed=0, dim=3)

    Computes the distance field from a set of bodies.
    compute the distance field located at nodes or centers of zone a (or zones in A), provided a list 
    of surfaces defining the bodies to which the distance is computed.

    Two algorithms are available:

    - type='ortho' means a distance computed by an orthogonal projection to the surface faces defined by bodies.
    - type='mininterf' returns the minimum distance of the point to the vertices of bodies.

    If loc='nodes', returns a distance computed at nodes of a (A), else if loc='centers, distance is computed at cell centers
    of a (A).

    Parameter 'signed'=1 enables to compute a signed distance (negative inside bodies). 
    When using signed distances, each body in bodies list must be a closed and watertight surface.
    In array version, cellnbodies provides the 'cellN' field for any vertex in bodies. Default value is 1.
    The algorithm 'ortho' does not take into account a body face if cellN=0 for all the vertices of that face.
    The algorithm 'mininterf' does not compute the distance to a vertex of cellN=0.
    
    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param bodies: body definition
    :type bodies: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param type: type of wall distance computation in ['ortho', 'mininterf']
    :type type: string
    :param loc: location of distance field in ['nodes', 'centers'] 
    :type loc: string
    :param signed: if 0 absolut distance, if 1 signed distance (negative inside)
    :type signed: int

    In the pyTree version, 'cellN' variable must be stored in bodies directly.
    If loc='nodes', the distance field is stored as a 'TurbulentDistance' field located at nodes, and 
    if loc='centers', it is stored in nodes located at centers.
    
    *Example of use:*

    * `Compute distance to walls (array) <Examples/Dist2Walls/distance2Walls.py>`_:

    .. literalinclude:: ../build/Examples/Dist2Walls/distance2Walls.py

    * `Compute distance to walls (pyTree) <Examples/Dist2Walls/distance2WallsPT.py>`_:

    .. literalinclude:: ../build/Examples/Dist2Walls/distance2WallsPT.py
    


.. toctree::
   :maxdepth: 2   


Indices and tables
###################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

