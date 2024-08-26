.. Post IBM documentation master file

:tocdepth: 2


Post.IBM: immersed boundary method specific post-processing 
===========================================================

Specific post-processing for immersed boundaries (IB).

These functions work with a solution tree "t", a geometry tree "tb", and/or a connectivity tree "tc".

All of these functions can be executed in both sequential and parallel contexts.


.. py:module:: Post.IBM

List of functions
#################


**-- Extract additional IBM information in the connectivity tree**

.. autosummary::
    :nosignatures:

    Post.IBM.extractPressureGradients
    Post.IBM.extractPressureHighOrder
    Post.IBM.extractYplusAtImagePoints

**-- IBM post-processing**

.. autosummary::
    :nosignatures:

    Post.IBM.prepareSkinReconstruction
    Post.IBM.computeSkinVariables

    .. Post.IBM.computeExtraVariables

    Post.IBM.computeAerodynamicLoads
    Post.IBM.computeAerodynamicCoefficients


Contents
########
.. py:function:: Post.IBM.extractPressureGradients(t, tc, secondOrder=False)

    Computes and extracts the pressure gradients from the solution tree (t) to the image points. The pressure gradients are computed within  the fluid domain using a second order technique. The pressure gradients are finally obtained at the image points with the same interpolation technique as for the immersed boundary treatment used during the simulation. The updated information is stored in the connectivity tree (tc) within the IBCD* sub-regions. The solution tree (t) remains unchanged.
    
    Also available as in-place (_extractPressureGradients).
    
    :param t: solution tree
    :type  t: [zone, list of zones, base, tree]
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :param secondOrder: if True, also computes second order pressure gradients
    :type  secondOrder: boolean
    :return: same as input with updated pressure gradient solution in each IBCD zone
    
    *Example of use:*

    * `extracts the pressure gradients at the image points (pyTree) <Examples/Post/extractPressureGradientsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressureGradientsPT.py

---------------------------------------

.. py:function:: Post.IBM.extractPressureHighOrder(tc, order=1)

    Extrapolates the wall pressure (1st or 2nd order) at the immersed boundaries and stores the solution in the connectivity tree. This function requires the pressure gradient information stored in the IBCD* sub-regions in the x, y, and z directions (see Post.IBM.extractPressureGradients).
    
    Also available as in-place (_extractPressureHighOrder).

    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :param order: int
    :type  order: 1 or 2
    :return: same as input with updated pressure solution in each IBCD zone

    *Example of use:*

    * `extrapolates the wall pressure (1st or 2nd order) at the immersed boundaries (pyTree) <Examples/Post/extractPressureHighOrderPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressureHighOrderPT.py

---------------------------------------

.. py:function:: Post.IBM.extractYplusAtImagePoints(tc)

    Extracts the yplus values at the image points and stores the solution in the connectivity tree. This function uses the yplus information located in the IBCD* subregions corresponding to the yplus values calculated at the target points during the simulation.

    Also available as in-place (_extractYplusAtImagePoints).

    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :return: same as input with updated yplus solution (yplusIP) in each IBCD zone

    *Example of use:*

    * `extracts the yplus values at the image points (pyTree) <Examples/Post/extractYplusAtImagePointsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractYplusAtImagePointsPT.py

---------------------------------------

.. py:function:: Post.IBM.prepareSkinReconstruction(tb, tc, dimPb=3, ibctypes=[])

    Prepares the flow solution extraction at immersed boundaries. This function extracts the IBM wall points originally stored in the connectivity tree to create a cloud of points that can be projected onto the surface. This function then pre-calculates and stores the intepolation data for the MLS projection (Moving Least Square, 3rd order), which can be particularly useful when more than one on-the-fly force extraction is requested during the simulation.

    When run in parallel, this function also automatically splits the case tree into NP parts (where NP is the number of MPI processes) and dispatches them between procs.

    If ibctypes is not empty, only the information associated with some immersed boundary conditions is calculated. 

    :param tb: surface mesh (TRI-type) with density, pressure, utau, and velocity variables
    :type tb: [zone, list of zones, base, tree]
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param ibctypes: list of IBC conditions
    :type dimPb: list of integers
    :return: Communication graph for IBM post-processing (for parallel use) and surface tree with interpolation data stored in each zone

    *Example of use:*

    * `prepares the flow solution extraction at immersed boundaries (pyTree) <Examples/Post/prepareSkinReconstructionPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/prepareSkinReconstructionPT.py

---------------------------------------

.. py:function:: Post.IBM.computeSkinVariables(ts, tc, graphIBCDPost, dimPb=3, ibctypes=[])

    Computes the surface flow solution at the wall using the IBM information updated and stored in the connectivity tree at each iteration of the flow simulation.
    This function operates in conjunction with Post.IBM.prepareSkinReconstruction, which must be called beforehand.

    Density, pressure and velocity vector (VelocityX, VelocityY, VelocityZ) information is updated in each zone of the surface tree.
    Depending on the immersed boundary condition, additional variables may also be updated, such as utau and yplus when using wall models.

    Also available as in-place (_computeSkinVariables).

    :param ts: surface mesh (TRI-type) with interpolation data stored in each zone
    :type ts: [zone, list of zones, base, tree]
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :param graphIBCDPost: communication graph for IBM post-processing
    :type graphIBCDPost: python dictionary
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param ibctypes: list of IBC conditions
    :type dimPb: list of integers
    :return: surface tree with updated flow information

    *Example of use:*

    * `computes the surface flow solution at the wall (pyTree) <Examples/Post/computeSkinVariablesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeSkinVariablesPT.py

---------------------------------------

.. py:function:: Post.IBM.computeAerodynamicLoads(ts, ts2=None, dimPb=3, famZones=[], Pref=None, center=(0.,0.,0.), verbose=0)

    Computes the aerodynamic loads acting on the immersed boundaries (forces and moments). This function computes the pressure and friction contributions separately. If an additional surface tree solution (ts2) is specified, the pressure solution is extracted from ts2 to ts, since ts2 is usually used to extract the pressure solution further away from the wall. A list of family names can also be specified to integrate only some parts of the geometry.

    This function returns a list of four lists containing the integration information in the body frame (aeroLoads = [forcePressure, forceFriction, momentPressure, momentFriction]).

    If Pref is omitted, its value is extracted from the reference state stored in ts.

    Warning: the pressure and friction coefficients, as well as the aerodynamic loads are still dimensionalized.

    :param ts: surface mesh (TRI-type) with updated flow solution
    :type ts: [zone, list of zones, base, tree]
    :param ts2: optional second surface mesh (TRI-type) with updated pressure solution
    :type ts2: [zone, list of zones, base, tree]
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param famZones: list of family names
    :type famZones: list of strings
    :param Pref: reference pressure for the pressure coefficient
    :type Pref: float
    :param center: reference center for the integration of the aerodynamic moments
    :type center: tuple of three floats
    :param verbose: if verbose > 0, print integration information
    :type verbose: 0 or 1
    :return: surface tree with updated flow information (Cp, Cf, etc.) and aeroLoads

    *Example of use:*

    * `computes the aerodynamic loads at the wall (pyTree) <Examples/Post/computeAerodynamicLoadsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeAerodynamicLoadsPT.py

---------------------------------------

.. py:function:: Post.IBM.computeAerodynamicCoefficients(ts, aeroLoads, dimPb=3, Sref=None, Lref=None, Qref=None, alpha=0., beta=0., verbose=0)

    Normalizes aerodynamic coefficients and places integration information in the wind frame. This function uses the a priori information obtained with Post.IBM.computeAerodynamicLoads.

    If Qref is omitted, its value is extracted from the reference state stored in ts.

    If Sref is omitted, its value is calculated as the area of the geometry surface.

    If Lref is omitted, it is set to Lref = 1 by default.

    :param ts: surface mesh (TRI-type) with extra information (dimensionalized pressure and friction coefficients)
    :type ts: [zone, list of zones, base, tree]
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param Sref: reference surface
    :type Sref: float
    :param Lref: reference length
    :type Lref: float
    :param Qref: reference dynamic pressure
    :type Qref: float
    :param alpha: angle of attack (x-z plane)
    :type alpha: float
    :param beta: angle of sideslip (x-y plane)
    :type beta: float
    :param verbose: if verbose > 0, print integration information
    :type verbose: 0 or 1
    :return: Surface tree with normalized Cp & Cf and updated aeroLoads with adimensionalized integration information placed in the wind frame.

    *Example of use:*

    * `normalizes aerodynamic coefficients and places integration information in the wind frame (pyTree) <Examples/Post/computeAerodynamicCoefficientsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeAerodynamicCoefficientsPT.py

---------------------------------------

.. .. py:function:: Post.IBM.computeExtraVariables(tb, PInf, QInf, variables=['Cp','Cf','frictionX','frictionY','frictionZ','frictionMagnitude','ShearStress'])

..     Computes additional variables required for the IBM post-processing. Uses density, pressure, utau, and velocity variables located at the vertices of the surface stored in the geometry tree (tb).

..     Possible extra variables are 'Cp', 'Cf', 'frictionX', 'frictionY', 'frictionZ', 'frictionMagnitude', 'ShearStress', 'gradnP' and 'gradtP'.

..     :param tb: surface mesh (TRI-type) with density, pressure, utau, and velocity variables.
..     :type tb: [zone, list of zones, base, tree]
..     :param PInf: reference pressure to compute Cp
..     :type PInf: real
..     :param QInf: reference dynamic pressure
..     :type QInf: real
..     :param variables: list of variables to be computed
..     :type variables: list of strings
..     :return: surface tree with additional variables located at the cell centers

..     *Example of use:*

..     * `computes additional variables required for the IBM post-processing (pyTree) <Examples/Post/computeExtraVariablesIBMPT.py>`_:

..     .. literalinclude:: ../build/Examples/Post/computeExtraVariablesIBMPT.py

.. ---------------------------------------

