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

    Compute vorticity on t from Velocity in centers. 
    If t contains ghost cells, set argument to True.
    Exists also as in place function (_computeVoriticity2) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param ghostCells: must be true if t contains ghost cells
    :type ghostCells: boolean
    :return: tree with vorticity
    :rtype: identical to input

    *Example of use:*

    * `Compute vorticity (pyTree) <Examples/Post/computeVorticity2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeVorticity2PT.py

--------------------

.. py:function:: Post.ExtraVariables2.computeVorticityMagnitude2(t, ghostCells=False)

    Compute vorticity magnitude on t from Velocity in centers. 
    If t contains ghost cells, set argument to True.
    Exists also as in place function (_computeVoriticityMagnitude2) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param ghostCells: must be true if t contains ghost cells
    :type ghostCells: boolean
    :return: tree with vorticity magnitude
    :rtype: identical to input

    *Example of use:*

    * `Compute vorticity magnitude (pyTree) <Examples/Post/computeVorticityMagnitude2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeVorticityMagnitude2PT.py

--------------------

.. py:function:: Post.ExtraVariables2.computeQCriterion2(t, ghostCells=False)

    Compute Q criterion on t from Velocity in centers. 
    If t contains ghost cells, set argument to True.
    Exists also as in place function (_computeQCriterion2) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param ghostCells: must be true if t contains ghost cells
    :type ghostCells: boolean
    :return: tree with Q criterion
    :rtype: identical to input

    *Example of use:*

    * `Compute Q criterion (pyTree) <Examples/Post/computeQCriterion2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeQCriterion2PT.py

--------------------

.. py:function:: Post.ExtraVariables2.computeLambda2(t, ghostCells=False)

    Compute lambda2 on t from Velocity in centers. 
    If t contains ghost cells, set argument to True.
    Exists also as in place function (_computeLambda2) that modifies t and returns None.

    :param t: input tree
    :type  t: [zone, list of zones, base, tree]
    :param ghostCells: must be true if t contains ghost cells
    :type ghostCells: boolean
    :return: tree with lambda2 criterion
    :rtype: identical to input

    *Example of use:*

    * `Compute lambda2 (pyTree) <Examples/Post/computeLambda2PT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeLambda2PT.py
