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

