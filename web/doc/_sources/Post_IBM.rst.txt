.. Post IBM documentation master file

:tocdepth: 2


Post.IBM: immersed boundary method specific post-processing 
===========================================================

Specific post-processing for immersed boundaries (IB).

These functions work with a solution tree "t", a geometry tree "tb", and/or a connectivity tree "tc".


.. py:module:: Post.IBM

List of functions
#################


**-- Post-processing for IB**

.. autosummary::

    Post.IBM.extractIBMWallFields
    Post.IBM.extractShearStress
    Post.IBM.extractLocalPressureGradients
    Post.IBM.extractYplusIP
    Post.IBM.extractPressureHO
    Post.IBM.extractPressureHO2
    Post.IBM.extractConvectiveTerms
    Post.IBM.computeExtraVariables
    Post.IBM.loads
    
    .. Post.IBM.computeSkinVariables
    .. Post.IBM.prepareSkinReconstruction
    .. Post.IBM.unsteadyLoads


Contents
########

.. py:function:: Post.IBM.extractIBMWallFields(tc, tb=None, coordRef='wall', famZones=[], IBCNames="IBCD_*", extractIBMInfo=False)

    Projects the solution computed and stored at IBM points onto the vertices of the surface.

    If tb is None, returns the cloud of IBM points.
    Else, the solution is projected onto the bodies, using a third-order accurate Moving Least Squares interpolation.
    
    Returns density, pressure, utau, yplus, velocity components. (Optional: yplus at image points, pressure gradients, curvature coefficient, temperature)

    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :param tb: surface mesh (TRI-type)
    :type tb: [zone, list of zones, base, tree]
    :param coordRef: reference coordinates for the cloud of IBM points (default is IBM wall points)
    :type coordRef: 'wall','target' or 'image'
    :param famZones: list of IBC families to be projected
    :type famZones: list of family names
    :param extractIBMInfo: if True, extracts all IBM point coordinates (wall, target and image points)
    :type extractIBMInfo: boolean
    :return: surface tree with the flow solution (density, pressure, friction velocity, yplus)

    *Example of use:*

    * `Projects the solution at IBM wall points onto the vertices of the surface (pyTree) <Examples/Post/extractIBMWallFieldsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractIBMWallFieldsPT.py

---------------------------------------

.. py:function:: Post.IBM.extractShearStress(tb)

    Computes the shear stress on the immersed boundary surface using utau values. Exists also as in-place (_extractShearStress).

    :param tb: surface mesh (TRI-type) with density, velocity, utau variable
    :type tb: [zone, list of zones, base, tree]
    :return: surface tree with the shear stress variables located at the cell centers ("ShearStressXX", "ShearStressYY", "ShearStressZZ", "ShearStressXY", "ShearStressXZ", "ShearStressYZ")

    *Example of use:*

    * `Computes the shear stress on the immersed boundary surface using utau values. (pyTree) <Examples/Post/computeShearStressPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeShearStressPT.py

---------------------------------------

.. py:function:: Post.IBM.extractLocalPressureGradients(tb)

    Computes the pressure gradients in the wall normal/tangent direction. Exists also as in-place (_extractLocalPressureGradients).

    :param tb: surface mesh (TRI-type) with pressure gradients in x, y and z directions
    :type tb: [zone, list of zones, base, tree]
    :return: surface tree with the normal and tangential pressure gradient variables located at the cell centers ("gradtP" and "gradnP")

    *Example of use:*

---------------------------------------

.. py:function:: Post.IBM.extractYplusIP(tc)

    Computes yplus values at image points and stores them in the tc. Exists also as in-place (_extractYplusIP).

    These new yplus values require yplus information located at target points as well as all IBM point coordinates (wall, target and image points)

    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :return: same as input with yplusIP field in each IBCD zone.

    *Example of use:*

---------------------------------------

.. py:function:: Post.IBM.extractPressureHO(tc, extractDensity=False)

    Extrapolates the wall pressure (1st order) at the immersed boundaries and stores the solution in the tc. Exists also as in-place (_extractPressureHO).

    Requires pressure gradient information in the x, y and z directions.
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :param extractDensity: if True, modifies the density solution using perfect gas law and the updated pressure solution
    :type  extractDensity: boolean
    :return: same as input with updated pressure solution in each IBCD zone.
    
    *Example of use:*

    * `1st order extrapolation of the pressure at the IB (pyTree) <Examples/Post/extractPressureHOPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressureHOPT.py


---------------------------------------

.. py:function:: Post.IBM.extractPressureHO2(tc, extractDensity=False)

    Extrapolates the wall pressure (2nd order) at the immersed boundaries and stores the solution in the tc. Exists also as in-place (_extractPressureHO2).

    Requires first and second order pressure gradient information in the x, y and z directions.

    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :param extractDensity: if True, modifies the density solution using perfect gas law and the updated pressure solution
    :type  extractDensity: boolean
    :return: same as input with updated pressure solution in each IBCD zone.
	     
    *Example of use:*

    * `2nd order extrapolation of the pressure at the IB (pyTree) <Examples/Post/extractPressureHO2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressureHO2PT.py

---------------------------------------

.. py:function:: Post.IBM.extractConvectiveTerms(tc)

    Computes the convective terms required for the thin boundary layers equations (TBLE) and stores them in the tc.

    Requires velocity gradient information in the x, y and z directions.
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :return:  same as input with convective terms in each IBCD zone (conv1: u*(du/dx) and conv2: v*(du/dy)).

    *Example of use:*

    * `Computes the convective terms (pyTree) <Examples/Post/extractConvectiveTermsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractConvectiveTermsPT.py

---------------------------------------

.. py:function:: Post.IBM.computeExtraVariables(tb, PInf, QInf, variables=['Cp','Cf','frictionX','frictionY','frictionZ','frictionMagnitude','ShearStress'])

    Computes additional variables required for the IBM post-processing. Uses density, pressure, utau, and velocity variables located at the vertices of tb.

    Possible extra variables are 'Cp', 'Cf', 'frictionX', 'frictionY', 'frictionZ', 'frictionMagnitude', 'ShearStress', 'gradnP' and 'gradtP'.

    :param tb: surface mesh (TRI-type) with density, pressure, utau, and velocity variables.
    :type tb: [zone, list of zones, base, tree]
    :param PInf: reference pressure to compute Cp
    :type PInf: real
    :param QInf: reference dynamic pressure
    :type QInf: real
    :param variables: list of variables to be computed
    :type variables: list of strings
    :return: surface tree with additional variables located at the cell centers

    *Example of use:*

    * `Computes variables using variables density, pressure, utau, and velocity at vertices of tb (pyTree) <Examples/Post/computeExtraVariablesIBMPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeExtraVariablesIBMPT.py

.. ---------------------------------------

.. py:function:: Post.IBM.loads(tb_in, tc_in=None, tc2_in=None, wall_out=None, alpha=0., beta=0., Sref=None, order=1, gradP=False, famZones=[])

    Computes the viscous and pressure forces on the immersed boundaries (IB). If tc_in=None, tb_in must also contain the projection of the flow field solution onto the surface.

    if tc and tc2 are not None, uses the pressure information at second image points.

    :param tb_in: geometry tree
    :type  tb_in: [zone, list of zones, base, tree]
    :param tc_in: connectivity tree 
    :type  tc_in: [zone, list of zones, base, tree, or None]
    :param tc2_in: connectivity tree containing IBM information for the second image point (optional)
    :type  tc2_in: [zone, list of zones, base, tree, or None]
    :param wall_out: file name for the output
    :type wall_out: string or None
    :param alpha: angle of attack (x-y plane) (in degrees)
    :type alpha: float
    :param beta: angle of attack (x-z plane) (in degrees)
    :type beta: float
    :param gradP: if True, extracts the wall pressure using pressure gradient information (see extractPressureHO() or extractPressureHO2())
    :type gradP: boolean
    :param order: pressure extrapolation order (when gradP is active)
    :type order: 1 or 2
    :param Sref: reference surface area for calculating the aerodynamic coefficients (CD/CL). if Sref is None, Sref is computed as the surface area
    :type Sref: float or None
    :param famZones: List of IBC families on which loads are computed.
    :type famZones: list of strings or None
    :return: surface tree with the flow solution. Lists of CD/CL per base.
       
    *Example of use:*

    * `Computes the viscous and pressure forces on an IB (pyTree) <Examples/Post/loadsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/loadsPT.py

---------------------------------------

.. .. py:function:: Post.IBM.unsteadyloads(tb, Sref=None, alpha=0., beta=0.)

    Computes the viscous and pressure forces on the IB during the computation of the solution. 

    :param tb: geometry tree with solution projected onto it
    :type  tb: [zone, list of zones, base, tree]
    :param Sref: reference surface area
    :type Sref: float or None
    :param alpha: Angle with respect to (0,Z) axe (in degrees)
    :type alpha: float
    :param beta: Angle with respect to (0,Y) axe (in degrees)
    :type beta: float
    :return: tree with the solution at the IB and the viscous and pressure loads
	     
    *Example of use:*

    * `Computes the viscous and pressure forces on an IB during the computation of the solution (pyTree) <Examples/Post/unsteadyloadsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/unsteadyloadsPT.py


.. ---------------------------------------

.. .. py:function:: Post.IBM.extractMassFlowThroughSurface(tb, t, famZones=[])

    Returns massflow through a surface defined by tb and returns tb. If famZones is a list of families, then only the
    zones of tb where the
    Currently: only sequential mode!

    :param tb: geometry tree
    :type  tb: [zone, list of zones, base, tree]
    :param t: solution tree with (Density,VelocityX, VelocityY, VelocityZ) stored at cell centers.
    :type t: pyTree 
    :param famZones: list of names of families of zones of tb where the massflow must be computed.
    :type famZones: list of strings
	     
    *Example of use:*

    * `Computes the massflow through an inlet surface (pyTree) <Examples/Post/extractMassFlowThroughSurfacePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractMassFlowThroughSurfacePT.py
