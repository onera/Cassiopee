.. RigidMotion documentation master file

RigidMotion: compute/define rigid motions
=========================================

Preamble
########

RigidMotion enables to define or compute rigid motions for arrays (as defined in Converter documentation) or for CGSN/Python trees (pyTrees).

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use with the array interface, you have to import RigidMotion module::

   import RigidMotion

For use with the pyTree interface::

    import RigidMotion.PyTree as RigidMotion


.. py:module:: RigidMotion


List of functions
##################

**-- Prescribed motions**

.. autosummary::

   RigidMotion.setPrescribedMotion1
   RigidMotion.setPrescribedMotion2
   RigidMotion.setPrescribedMotion3

**-- General functions**

.. autosummary::

   RigidMotion.evalPosition

   


Contents
#########

.. py:function:: RigidMotion.setPrescribedMotion1(a, motionName, tx="0", ty="0", tz="0", cx="0", cy="0", cz="0", ex="0", ey="0", ez="0", angle="0")

    Set a precribed motion defined by a translation of the origin (tx,ty,tz), the center of a rotation (cx,cy,cz), the second point of the rotation axis (ex,ey,ez) and the rotation angle in degrees. They can depend on time {t}.

    :param a: Input data
    :type  a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tx: translation in x motion string
    :type tx: string
    :param ty: translation in y motion string
    :type ty: string
    :param tz: translation in z motion string
    :type tz: string
    :param cx: rotation center x coordinate motion string
    :type cx: string
    :param cy: rotation center y coordinate motion string
    :type cy: string
    :param cz: rotation center z coordinate motion string
    :type cz: string
    :param ex: rotation axis x coordinate motion string
    :type ex: string
    :param ey: rotation axis y coordinate motion string
    :type ey: string
    :param ez: rotation axis z coordinate motion string
    :type ez: string
    :param angle: rotation angle motion string
    :type angle: string
    
    *Example of use:*

    * `Set a prescribed motion of type 1 (pyTree) <Examples/RigidMotion/setPrescribedMotion1PT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/setPrescribedMotion1PT.py

------------------------------------------------------------------------------------------------

.. py:function:: RigidMotion.setPrescribedMotion2(a, motionName, transl_speed, psi0, pis0_b, alp_pnt, alp_vct, alp0, rot_pnt, rot_vct, rot_omg, del_pnt, del_vct, del0, delc, dels, bet_pnt, bet_vct, bet0, betc, bets, tet_pnt, tet_vct, tet0, tetc, tets, span_vct, pre_lag_pnt, pre_lag_vct, pre_lag_ang, pre_con_pnt, pre_con_vct, pre_con_ang)

    Set a precribed motion defined by a elsA rotor motion. Arguments are identical to elsA rotor motion. 

    :param a: Input data
    :type  a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tx: translation in x motion string
    :type tx: string
    :param ty: translation in y motion string
    :type ty: string
    :param tz: translation in z motion string
    :type tz: string
    :param cx: rotation center x coordinate motion string
    :type cx: string
    :param cy: rotation center y coordinate motion string
    :type cy: string
    :param cz: rotation center z coordinate motion string
    :type cz: string
    :param ex: rotation axis x coordinate motion string
    :type ex: string
    :param ey: rotation axis y coordinate motion string
    :type ey: string
    :param ez: rotation axis z coordinate motion string
    :type ez: string
    :param angle: rotation angle motion string
    :type angle: string
    
    *Example of use:*

    * `Set a prescribed motion of type 2 (pyTree) <Examples/RigidMotion/setPrescribedMotion2PT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/setPrescribedMotion2PT.py

------------------------------------------------------------------------------------------------

------------------------------------------------------------------------------------------------

.. py:function:: RigidMotion.setPrescribedMotion3(a, motionName, transl_speed, axis_pnt, axis_vct, omega)

    Set a precribed motion defined by a constant speed rotation and translation. 
    omega is in rad/time unit.

    :param a: Input data
    :type  a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param tx: translation in x motion string
    :type tx: string
    :param ty: translation in y motion string
    :type ty: string
    :param tz: translation in z motion string
    :type tz: string
    :param cx: rotation center x coordinate motion string
    :type cx: string
    :param cy: rotation center y coordinate motion string
    :type cy: string
    :param cz: rotation center z coordinate motion string
    :type cz: string
    :param ex: rotation axis x coordinate motion string
    :type ex: string
    :param ey: rotation axis y coordinate motion string
    :type ey: string
    :param ez: rotation axis z coordinate motion string
    :type ez: string
    :param angle: rotation angle motion string
    :type angle: string
    
    *Example of use:*

    * `Set a prescribed motion of type 3 (pyTree) <Examples/RigidMotion/setPrescribedMotion3PT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/setPrescribedMotion3PT.py

------------------------------------------------------------------------------------------------




General functions
---------------------

.. py:function:: RigidMotion.PyTree.evalPosition(a, time)

    Evaluate the position at time t according to a motion.
    If the motion is defined in a with setPrescribedMotion. 

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param time: evaluation time
    :type time: float
    :return: reference copy of a
    :rtype: identical to input

    *Example of use:*

    * `Evaluate position (pyTree) <Examples/RigidMotion/evalPositionPT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/evalPositionPT.py
    

    Evaluate position at given time, when motion is described by a function. F(t) is a 
    function describing motion. F(t) = (centerAbs(t), centerRel(t), rot(t)), 
    where centerAbs(t) are the coordinates of the rotation center in the absolute frame, centerRel(t) are the coordinates of the rotation center in the relative 
    (that is array's) frame and rot(t), the rotation matrix.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param time: evaluation time
    :type time: float
    :param F: motion function
    :type F: python function
    :return: reference copy of a
    :rtype: identical to input

    *Example of use:*

    * `Evaluate position with function (pyTree) <Examples/RigidMotion/evalPosition2PT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/evalPosition2PT.py
    

.. toctree::
   :maxdepth: 2   


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

