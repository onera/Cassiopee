.. Modeler documentation master file

:tocdepth: 2

Modeler: some basic models
==========================

Preamble
########

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use with the array interface, you have to import Modeler module::

   import Modeler.Models as Models

For use with the pyTree interface::

    import Modeler.PyTree as Models


.. py:module:: Modeler

List of functions
##################

**-- Boxes**

.. autosummary::

    Modeler.Models.skySphere
    


Boxes
-------


.. py:function:: Modeler.Models.skySphere(Xc, R)

    
    Create a sky sphere of center Xc and radius R. The sphere is uv mapped.
    This sphere can be used to map 360 images obtained with equirectangular projection.

    :param Xc: center of sphere
    :type Xc: tuple of 3 floats
    :param R: radius of sphere
    :type R: float
    
    *Example of use:*

    * `Create a sky sphere (array) <Examples/Modeler/skySphere.py>`_:

    .. literalinclude:: ../build/Examples/Modeler/skySphere.py

    * `Create a sky sphere (pyTree) <Examples/Modeler/skySpherePT.py>`_:

    .. literalinclude:: ../build/Examples/Modeler/skySpherePT.py


.. toctree::
   :maxdepth: 2


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

