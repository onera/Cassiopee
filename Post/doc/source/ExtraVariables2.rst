.. Post.ExtraVariables2 documentation master file


ExtraVariables2: derived fields from primitive variables 
=========================================================

This module compute derived fields from primitive variables.


.. py:module:: Post.ExtraVariables2


List of functions
##################


**-- Volume fields**

.. autosummary::

    Post.ExtraVariables2.extractTree
    Post.ExtraVariables2.computeVorticity2
    Post.ExtraVariables2.computeVorticityMagnitude2
    Post.ExtraVariables2.computeQCriterion2
    Post.ExtraVariables2.computeLambda2
    Post.ExtraVariables2.extractPressure
    Post.ExtraVariables2.extractVelocityMagnitude
    Post.ExtraVariables2.extractMach
    Post.ExtraVariables2.extractViscosityMolecular
    Post.ExtraVariables2.extractViscosityMolecular
    Post.ExtraVariables2.extractMutSurMu

**-- Surface fields**

.. autosummary::

    Post.ExtraVariables2.extractShearStress
    Post.ExtraVariables2.extractTaun
    Post.ExtraVariables2.extractPn
    Post.ExtraVariables2.extractForce
    Post.ExtraVariables2.extractFrictionVector
    Post.ExtraVariables2.extractFrictionMagnitude
    Post.ExtraVariables2.extractUTau

**-- 1D profiles**
    Post.ExtraVariables2.extractProfile
    Post.ExtraVariables2.extractyplus


Contents
#########

Volume fields
--------------------


.. py:function:: Post.ExtraVariables2.extractTree(t, vars=['centers:Density', 'centers:VelocityX', 'centers:VelocityY', 'centers:VelocityZ', 'centers:Temperature', 'centers:TurbulentSANuTilde'])

    Keep only some variables from tree. This is just a reference tree (no extra memory is used).

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param vars: list of vars to keep in returned tree
    :type vars: list of strings
    :return: tree with selected variables
    :rtype: identical to input

    *Example of use:*

    * `Extract tree (pyTree) <Examples/Post/extractTreePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractTreePT.py

--------------------

.. py:function:: Post.ExtraVariables2.computeVorticity2(t, ghostCells=False)

    Compute vorticity on t from Velocity field in centers. 
    If t contains ghost cells, set argument to True.
    Exists also as in place function (_computeVoriticity2) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param ghostCells: must be true if t contains ghost cells
    :type ghostCells: boolean
    :return: tree with "VorticityX,"VorticityY","VorticityZ" in centers
    :rtype: identical to input

    *Example of use:*

    * `Compute vorticity (pyTree) <Examples/Post/computeVorticity2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeVorticity2PT.py

--------------------

.. py:function:: Post.ExtraVariables2.computeVorticityMagnitude2(t, ghostCells=False)

    Compute vorticity magnitude on t from Velocity field in centers. 
    If t contains ghost cells, set argument to True.
    Exists also as in place function (_computeVoriticityMagnitude2) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param ghostCells: must be true if t contains ghost cells
    :type ghostCells: boolean
    :return: tree with "VorticityMagnitude" in centers
    :rtype: identical to input

    *Example of use:*

    * `Compute vorticity magnitude (pyTree) <Examples/Post/computeVorticityMagnitude2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeVorticityMagnitude2PT.py

--------------------

.. py:function:: Post.ExtraVariables2.computeQCriterion2(t, ghostCells=False)

    Compute Q criterion on t from Velocity field in centers. 
    If t contains ghost cells, set argument to True.
    Exists also as in place function (_computeQCriterion2) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param ghostCells: must be true if t contains ghost cells
    :type ghostCells: boolean
    :return: tree with "QCriterion" in centers
    :rtype: identical to input

    *Example of use:*

    * `Compute Q criterion (pyTree) <Examples/Post/computeQCriterion2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeQCriterion2PT.py

--------------------

.. py:function:: Post.ExtraVariables2.computeLambda2(t, ghostCells=False)

    Compute lambda2 on t from Velocity field in centers. 
    If t contains ghost cells, set argument to True.
    Exists also as in place function (_computeLambda2) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param ghostCells: must be true if t contains ghost cells
    :type ghostCells: boolean
    :return: tree with "lambda2" in centers
    :rtype: identical to input

    *Example of use:*

    * `Compute lambda2 (pyTree) <Examples/Post/computeLambda2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeLambda2PT.py


--------------------

.. py:function:: Post.ExtraVariables2.extractPressure(t)

    Compute Pressure on t from Temperature and Density field in centers with P = ro r T. 
    The tree t must have a ReferenceState node.
    Cv and Gamma are taken from ReferenceState and r = Cv * (Gamma-1).
    Exists also as in place function (_extractPressure) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :return: tree with "Pressure" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract pressure (pyTree) <Examples/Post/extractPressurePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPressurePT.py

-------------------------------------------

.. py:function:: Post.ExtraVariables2.extractVelocityMagnitude(t)

    Compute velocity magnitude on t from Velocity field in centers. 
    Exists also as in place function (_extractVelocityMagnitude) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :return: tree with "VelocityMagnitude" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract velocity magnitude (pyTree) <Examples/Post/extractVelocityMagnitudePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractVelocityMagnitudePT.py


--------------------

.. py:function:: Post.ExtraVariables2.extractMach(t)

    Compute Mach on t from Velocity, Temperature and Density field in centers with M = u/sqrt(gamma p/ro) and p = ro r T. 
    The tree t must have a ReferenceState node.
    Cv and Gamma are taken from ReferenceState and r = Cv * (Gamma-1).
    Exists also as in place function (_extractMach) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :return: tree with "Mach" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract mach (pyTree) <Examples/Post/extractMachPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractMachPT.py


--------------------

.. py:function:: Post.ExtraVariables2.extractViscosityMolecular(t)

    Compute ViscosityMolecular on t from Temperature field in centers with Sutherland law. 
    The tree t must have a ReferenceState node.
    Cs, Mus, Ts are taken from ReferenceState.
    Exists also as in place function (_extractViscosityMolecular) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :return: tree with "ViscosityMolecular" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract viscosity molecular (pyTree) <Examples/Post/extractViscosityMolecularPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractViscosityMolecularPT.py

--------------------

.. py:function:: Post.ExtraVariables2.extractViscosityEddy(t)

    Compute ViscosityEddy on t from TurbulentSANuTilde, ViscosityMolecular and Density field in centers with 
    kappa = ro * nutilde / mu
    and mut = ro * nutilde * kappa^3 / (kappa^3 + 7.1^3). 
    Exists also as in place function (_extractViscosityEddy) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :return: tree with "ViscosityEddy" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract viscosity eddy (pyTree) <Examples/Post/extractViscosityEddyPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractViscosityEddyPT.py

--------------------

.. py:function:: Post.ExtraVariables2.extractMutSurMu(t)

    Compute ViscosityEddy divided by ViscosityMolecular on t 
    from ViscosityEddy and ViscosityMolecular in centers. 
    Exists also as in place function (_extractMutSurMu) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :return: tree with "MutSurMu" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract Mut over Mu (pyTree) <Examples/Post/extractMutSurMuPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractMutSurMuPT.py




Surface fields
--------------------

.. py:function:: Post.ExtraVariables2.extractShearStress(teff)

    Compute ShearStress on teff 
    from ViscosityMolecular and gradxVelocityX,... in centers. 
    Exists also as in place function (_extractShearStress) that modifies t and returns None.

    :param teff: input tree
    :type  teff: [zone, list of zones, base, tree]
    :return: tree with "ShearStressXX,XY,XZ,YY,YZ,ZZ" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract shearStress (pyTree) <Examples/Post/extractShearStressPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractShearStressPT.py


---------------------------

.. py:function:: Post.ExtraVariables2.extractTaun(teff)

    Compute tau.n on teff from ShearStress in centers. 
    Exists also as in place function (_extractTaun) that modifies t and returns None.

    :param teff: input tree
    :type  teff: [zone, list of zones, base, tree]
    :return: tree with "taunx,y,z" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract tau.n (pyTree) <Examples/Post/extractTaunPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractTaunPT.py


---------------------------

.. py:function:: Post.ExtraVariables2.extractPn(teff)

    Compute P.n on teff from Pressure in centers. 
    Exists also as in place function (_extractPn) that modifies t and returns None.

    :param teff: input tree
    :type  teff: [zone, list of zones, base, tree]
    :return: tree with "Pnx,y,z" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract P.n (pyTree) <Examples/Post/extractPnPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractPnPT.py


---------------------------

.. py:function:: Post.ExtraVariables2.extractForce(teff, withPInf=None)

    Compute the force field on teff from Pressure and ShearStress in centers. 
    If withPinf is None: F = -p.n + tau.n
    Else: F = -(p-pinf).n + tau.n
    Exists also as in place function (_extractForce) that modifies t and returns None.

    :param teff: input tree
    :type  teff: [zone, list of zones, base, tree]
    :param withPinf: None or infinite field pressure
    :type withPinf: None or float
    :return: tree with "Fx,y,z" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract Force (pyTree) <Examples/Post/extractForcePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractForcePT.py

---------------------------

.. py:function:: Post.ExtraVariables2.extractFrictionVector(teff)

    Compute the friciton vector on teff from ShearStress in centers
    with taut = tau.n - (n. tau.n) n.
    Exists also as in place function (_extractFrictionVector) that modifies t and returns None.

    :param teff: input tree
    :type  teff: [zone, list of zones, base, tree]
    :return: tree with "FrictionX,FrictionY,FrictionZ" in centers
    :rtype: identical to input

    *Example of use:*

    * `Extract Force (pyTree) <Examples/Post/extractFrictionVectorPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractFrictionVectorPT.py



