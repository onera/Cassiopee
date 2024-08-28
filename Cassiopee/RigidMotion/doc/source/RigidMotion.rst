.. RigidMotion documentation master file

:tocdepth: 2

RigidMotion: compute/define rigid motions
=========================================

Preamble
########

RigidMotion enables to define or compute rigid motions for arrays (as defined in Converter documentation) or for CGNS/Python trees (pyTrees).

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
   :nosignatures:

   RigidMotion.PyTree.setPrescribedMotion1
   RigidMotion.PyTree.setPrescribedMotion2
   RigidMotion.PyTree.setPrescribedMotion3

**-- General functions**

.. autosummary::
   :nosignatures:

   RigidMotion.PyTree.evalPosition
   RigidMotion.PyTree.evalGridSpeed
   RigidMotion.PyTree.copyGrid2GridInit
   RigidMotion.PyTree.copyGridInit2Grid
   

Contents
#########

.. py:function:: RigidMotion.PyTree.setPrescribedMotion1(a, motionName, tx="0", ty="0", tz="0", cx="0", cy="0", cz="0", ex="0", ey="0", ez="0", angle="0")

    Set a prescribed motion defined by a translation of the origin (tx,ty,tz), the center of a rotation (cx,cy,cz), the second point of the rotation axis (ex,ey,ez) and the rotation angle in degrees. They can depend on time {t}.
    
    Exists also as an in-place version (_setPrescribedMotion1) which modifies a and returns None.
    
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

.. py:function:: RigidMotion.PyTree.setPrescribedMotion2(a, motionName, transl_speed, psi0, pis0_b, alp_pnt, alp_vct, alp0, rot_pnt, rot_vct, rot_omg, del_pnt, del_vct, del0, delc, dels, bet_pnt, bet_vct, bet0, betc, bets, tet_pnt, tet_vct, tet0, tetc, tets, span_vct, pre_lag_pnt, pre_lag_vct, pre_lag_ang, pre_con_pnt, pre_con_vct, pre_con_ang)

    Set a prescribed motion defined by a rigid rotor motion. Arguments are identical to elsA rotor motion. 

    Exists also as an in-place version (_setPrescribedMotion2) which modifies a and returns None.

    :param a: Input data
    :type  a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param transl_speed: translation speed
    :type transl_speed: a 3-tuple of floats
    :param psi0: initial pitch angle (in degrees)
    :type psi0: float
    :param psi0_b: angle for blade position wrt leading blade (in degrees)
    :type psi0_b: float
    :param alp_pnt: origin of rotor shaft
    :type alp_pnt: a 3-tuple of floats
    :param alp_vct: axis of rotor shaft
    :type alp_vct: a 3-tuple of floats
    :param alp0: rotor shaft angle (in degrees)
    :type alp0: float
    :param rot_pnt: rotation center 
    :type rot_pnt: 3-tuple of floats
    :param rot_vct: rotation axis
    :type rot_vct: 3-tuple of floats
    :param rot_omg: rotor angular velocity (in radians per sec)
    :type rot_omg:  float

    :param del_pnt: origin of lead-lag 
    :type del_pnt: 3-tuple of floats
    :param del_vct: lead-lag axis
    :type del_vct:  3-tuple of floats
    :param del0: lead-lag angle (in degrees)
    :type del0: float
    :param delc: cosine part of harmonics for lead-lag
    :type delc: tuple of floats
    :param dels: sine part of harmonics for lead-lag
    :type dels: tuple of floats

    :param bet_pnt: origin of flapping motion
    :type bet_pnt: 3-tuple of floats
    :param bet_vct: flapping axis
    :type bet_vct:  3-tuple of floats
    :param bet0: flapping angle (in degrees)
    :type bet0: float
    :param betc: cosine part of harmonics for conicity
    :type betc: tuple of floats
    :param bets: sine part of harmonics for conicity
    :type bets: tuple of floats

    :param tet_pnt: origin of pitching motion
    :type tet_pnt: 3-tuple of floats
    :param tet_vct: pitching axis
    :type tet_vct:  3-tuple of floats
    :param tet0: collective pitch angle (in degrees)
    :type tet0: float
    :param tetc: cyclic pitch cosine part
    :type tetc: tuple of floats
    :param tets: cyclic pitch sine part
    :type tets: tuple of floats
 
    :param span_vct: reference blade spanwise axis
    :type span_vct: 3-tuple of floats
    
    :param pre_lag_pnt: origin of pre-lag 
    :type pre_lag_pnt: 3-tuple of floats
    :param pre_lag_vct: pre-lag axis
    :type pre_lag_vct: 3-tuple of floats
    :param pre_lag_ang: pre-lag angle (in degrees)
    :type pre_lag_ang: float
    :param pre_con_pnt: origin of pre-conicity
    :type pre_con_pnt: 3-tuple of floats
    :param pre_con_vct: pre-conicity axis
    :type pre_con_vct: 3-tuple of floats
    :param pre_con_ang: pre-conicity angle (in degrees)
    :type pre_con_ang: float

    *Example of use:*

    * `Set a prescribed motion of type 2 (pyTree) <Examples/RigidMotion/setPrescribedMotion2PT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/setPrescribedMotion2PT.py

------------------------------------------------------------------------------------------------

.. py:function:: RigidMotion.PyTree.setPrescribedMotion3(a, motionName, transl_speed, axis_pnt, axis_vct, omega)

    Set a precribed motion defined by a constant speed rotation and constant translation vector. 
    omega is in rad/time unit.
    Since rotation is applied before translation, the center of rotation (axis_pnt) is
    moving with translation speed also.

    Exists also as an in-place version (_setPrescribedMotion3) which modifies a and returns None.

    :param a: Input data
    :type  a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param transl_speed: translation vector
    :type transl_speed: tuple of 3 floats
    :param axis_pnt: rotation axis (constant in translated frame)
    :type axis_pnt: tuple of 3 floats
    :param axis_vect: vector axis (constant in traslated frame)
    :type axis_vect: tuple of 3 floats
    :param omega: constant rotation speed
    :type omega: float 
    
    *Example of use:*

    * `Set a prescribed motion of type 3 (pyTree) <Examples/RigidMotion/setPrescribedMotion3PT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/setPrescribedMotion3PT.py

------------------------------------------------------------------------------------------------




General functions
---------------------

.. py:function:: RigidMotion.PyTree.evalPosition(a, time)

    Evaluate the position at time t according to a motion.
    The motion must be defined in a with setPrescribedMotion.
    If GridCoordinates#Init is present, it is used to compute position. 
    Otherwise, Grid coordinates in a must be the coordinates at time=0.

    Exists also as an in-place version (_evalPosition) which modifies a and returns None.

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
    

------------------------------------------------------------------------------------------------

.. py:function:: RigidMotion.PyTree.evalGridSpeed(a, time)

    Evaluate grid speed at given time.
    The position must already have been evaluated at this time.

    Exists also as an in-place version (_evalGridSpeed) which modifies a and returns None.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param time: evaluation time
    :type time: float
    :return: reference copy of a
    :rtype: identical to input

    *Example of use:*

    * `Evaluate speed (pyTree) <Examples/RigidMotion/evalGridSpeedPT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/evalGridSpeedPT.py
    

------------------------------------------------------------------------------------------------

.. py:function:: RigidMotion.PyTree.copyGrid2GridInit(a, mode=0)

    Copy GridCoordinates to GridCoordinates#Init.
    If mode=0, only if grid has a TimeMotion node.
    If mode=1, always copy.

    Exists also as an in-place version (_copyGrid2GridInit) which modifies a and returns None.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param mode: behaviour
    :type mode: 0 or 1
    :return: reference copy of a
    :rtype: identical to input

    *Example of use:*

    * `Copy GridCoordinates to GridCoordinates#Init (pyTree) <Examples/RigidMotion/copyGrid2GridInitPT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/copyGrid2GridInitPT.py

------------------------------------------------------------------------------------------------

.. py:function:: RigidMotion.PyTree.copyGridInit2Grid(a)

    Copy GridCoordinates#Init to GridCoordinates if it exists.
    
    Exists also as an in-place version (_copyGridInit2Grid) which modifies a and returns None.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :return: reference copy of a
    :rtype: identical to input

    *Example of use:*

    * `Copy GridCoordinates#Init to GridCoordinates (pyTree) <Examples/RigidMotion/copyGridInit2GridPT.py>`_:

    .. literalinclude:: ../build/Examples/RigidMotion/copyGridInit2GridPT.py


---------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

