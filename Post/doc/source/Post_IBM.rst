.. Post IBM documentation master file


Post.IBM: immersed boundary method specific post-processing 
===========================================================

Specific post-processing for immersed boundaries (IB).

These functions work with a solution tree "t", a geometry tree "tb", and/or a connectivity tree "tc".


.. py:module:: Post.IBM

List of functions
#################


**-- Post-processing of IBs**

.. autosummary::

    .. Post.IBM.computeSkinVariables
    Post.IBM.extractIBMWallFields
    Post.IBM.extractIBMInfo
    Post.IBM.extractShearStress
    Post.IBM.computeExtraVariables
    .. Post.IBM.loads
    Post.IBM.extractPressureHO
    Post.IBM.extractPressureHO2
    Post.IBM.extractConvectiveTerms
    .. Post.IBM.prepareSkinReconstruction
    .. Post.IBM.unsteadyLoads


Contents
########

.. py:function:: Post.IBM.extractIBMInfo(tc,filename='IBMInfo.cgns')

    Extracts the geometrical information required for the IBM (i.e. wall points, target points, and image points).
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :return: tree with geometrical information required for the IBM

    *Example of use:*

    * `Extract the IBM geometrical information (pyTree) <Examples/Post/extractIBMInfoPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractIBMInfoPT.py

---------------------------------------

.. py:function:: Post.IBM.extractIBMWallFields(tc, tb=None, coordRef='wall', famZones=[], front=1)

    Project the solution at IBM wall points onto the vertices of the surface.
    If tb is None, returns the cloud of points, 
    else interpolate by Moving Least Squares interpolation on tb.
    Returns density, pressure, utau, yplus, velocity components.

    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :param tb: surface mesh (TRI-type)
    :type tb: [zone, list of zones, base, tree]
    :param coordRef: coordinates of IBM points to be projected (default is IBM wall points) :'wall','cible','image'
    :type coordRef: string
    :param famZones: list of families of IBCs to be projected
    :type famZones: list of names of families
    :return: surface tree with solution (density, pressure, friction velocity, yplus)

    *Example of use:*

    * `Project the solution at IBM wall points onto the vertices of the surface (pyTree) <Examples/Post/extractIBMWallFieldsPT.py>`_:

    .. .. literalinclude:: ../build/Examples/Post/extractIBMWallFieldsPT.py

---------------------------------------

.. py:function:: Post.IBM.extractShearStress(tb)

    Computes the shear stress using utau values at vertices of the surface mesh

    :param tb: surface mesh (TRI-type) with density, velocity, utau variable
    :type tb: [zone, list of zones, base, tree]
    :return: surface tree with shear stress variables "ShearStressXX", "ShearStressYY","ShearStressZZ","ShearStressXY", "ShearStressXZ", "ShearStressYZ"

    *Example of use:*

    * `Computes the shear stress using utau values at vertices of the surface mesh (pyTree) <Examples/Post/computeShearStressPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeShearStressPT.py

---------------------------------------

.. py:function:: Post.IBM.computeExtraVariables(tb, PInf, QInf, variables=['Cp','Cf','frictionX','frictionY','frictionZ', 'frictionMagnitude','ShearStress'])

    Computes variables using variables density, pressure, utau, and velocity at vertices of tb.
    Solution is located at cell centers.

    :param tb: surface mesh (TRI-type) with density, velocity, utau variable
    :type tb: [zone, list of zones, base, tree]
    :param PInf: reference pressure to compute Cp
    :type PInf: real
    :param QInf: reference dynamic pressure
    :type QInf: real
    :param variables: list of variables to be computed.
    :type variables: list of strings
    :return: surface tree with additional variables.

    *Example of use:*

    * `Computes variables using variables density, pressure, utau, and velocity at vertices of tb (pyTree) <Examples/Post/computeExtraVariablesIBMPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeExtraVariablesIBMPT.py

.. ---------------------------------------

.. .. py:function:: Post.IBM.loads(t_case, tc_in=None, tc2_in=None, wall_out=None, alpha=0., beta=0., gradP=False, order=1, Sref=None, famZones=[])

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
    :param famZones: name of families of immersed boundaries on whih loads are computed
    :type famZones: list of strings or None
    :return: tree with the solution at the IB and the viscous and pressure loads
       
    *Example of use:*

    * `Computes the viscous and pressure forces on an IB (pyTree) <Examples/Post/loadsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/loadsPT.py

---------------------------------------

.. py:function:: Post.IBM.extractConvectiveTerms(tc)

    Computes the convective terms required for the thin boundary layers equations (TBLE) and stores them in the tc.
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :return: same as input

    *Example of use:*

    * `Compute the convective terms (pyTree) <Examples/Post/extractConvectiveTermsPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractConvectiveTermsPT.py

---------------------------------------

.. py:function:: Post.IBM.extractPressureHO(tc)

    1st order extrapolation of the pressure at the IB.
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :return: same as input
    
    *Example of use:*

    * `1st order extrapolation of the pressure at the IB (pyTree) <Examples/Post/extractPressureHOPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressureHOPT.py


---------------------------------------

.. py:function:: Post.IBM.extractPressureHO2(tc)

    2nd order extrapolation of the pressure at the IB.
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, base, tree]
    :return: same as input
	     
    *Example of use:*

    * `2nd order extrapolation of the pressure at the IB (pyTree) <Examples/Post/extractPressureHO2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressureHO2PT.py


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
