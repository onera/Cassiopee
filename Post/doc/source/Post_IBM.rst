.. Post IBM documentation master file


IBM: Immersed boundary method specific post-processing 
=======================================================

Specific post-processing for immersed boundaries.

These functions work with a solution tree "t", a geometry tree "tb", and/or a connectivity tree "tc".


.. py:module:: Post.IBM

List of functions
#################


**-- Post-processing of IBs**

.. autosummary::

    #Post.IBM.computeSkinVariables
    Post.IBM.extractConvectiveTerms
    Post.IBM.extractIBMInfo
    Post.IBM.extractPressureHO
    Post.IBM.extractPressureHO2
    Post.IBM.loads
    #Post.IBM.prepareSkinReconstruction
    #Post.IBM.unsteadyLoads


Contents
########

.. py:function:: Post.IBM.extractConvectiveTerms(tc)

    Computes the convective terms required for the thin boundary layers equations (TBLE) and stores them in the tc.
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]

    *Example of use:*

    * `Compute the convective terms (pyTree) <Examples/Post/extractConvectiveTermsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractConvectiveTermsPT.py


---------------------------------------

.. py:function:: Post.IBM.extractIBMInfo(tc,filename='IBMInfo.cgns')

    Extracts the geometrical information required for the IBM (i.e. wall points, target points, and image points).
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]

    *Example of use:*

    * `Extract the IBM geometrical information (pyTree) <Examples/Post/extractIBMInfoPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractIBMInfoPT.py


---------------------------------------

.. py:function:: Post.IBM.extractPressureHO(tc)

    1st order extrapolation of the pressure at the immersed boundary (IB).
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]

    *Example of use:*

    * `1st order extrapolation of the pressure at the IB (pyTree) <Examples/Post/extractPressureHOPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressureHOPT.py


---------------------------------------

.. py:function:: Post.IBM.extractPressureHO2(tc)

    2nd order extrapolation of the pressure at the immersed boundary (IB).
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]

    *Example of use:*

    * `2nd order extrapolation of the pressure at the IB (pyTree) <Examples/Post/extractPressureHO2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressureHO2PT.py


---------------------------------------

.. py:function:: Post.IBM.loads(t_case, tc_in=None, tc2_in=None, wall_out=None, alpha=0., beta=0., gradP=False, order=1, Sref=None, famZones=[])

    Computes the viscous and pressure forces on the IB. If tc_in=None, t_case must also contain the projection of the flow field solution onto the IB.

    :param t_case: geometry tree
    :type  t_case: [zone, list of zones, base, tree]
    :param tc_in: connectivity tree 
    :type  tc_in: [zone, list of zones, base, tree, or None]
    :param tc2_in: connectivity tree of second image point (if present)
    :type  tc2_in: [zone, list of zones, base, tree, or None]
    :param wall_out: file name for the output of the forces at the wall and at the cell centers
    :type wall_out: string or None
    :param alpha: Angle with respect to (0,Z) axe (in degrees)
    :type alpha: float
    :param beta: Angle with respect to (0,Y) axe (in degrees)
    :type beta: float
    :param gradP: calculate the pressure gradient?
    :type gradP: boolean
    :param order: pressure extrapolation order
    :type order: integer
    :param Sref: reference surface area
    :type Sref: float or None
    :param famZones: name of familys for which IBM data is extracted
    :type famZones: list of strings or None

    *Example of use:*

    * `Computes the viscous and pressure forces on an IB (pyTree) <Examples/Post/loadsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/loadsPT.py


---------------------------------------

.. py:function:: Post.IBM.unsteadyloads(tb, Sref=None, alpha=0., beta=0.)

    Computes the viscous and pressure forces on the IB during the computation of the solution. 

    :param tb: geometry tree with solution projected onto it
    :type  tb: [zone, list of zones, base, tree]
    :param Sref: reference surface area
    :type Sref: float or None
    :param alpha: Angle with respect to (0,Z) axe (in degrees)
    :type alpha: float
    :param beta: Angle with respect to (0,Y) axe (in degrees)
    :type beta: float

    *Example of use:*

    * `Computes the viscous and pressure forces on an IB during the computation of the solution (pyTree) <Examples/Post/unsteadyloadsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/unsteadyloadsPT.py
