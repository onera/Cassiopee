.. Post documentation master file


Rotor: rotor specific post-processing 
======================================

Specific post-processing for rotors or propellers.

Those function work with a solution tree "t" or a blade surface tree "teff".


.. py:module:: Post.Rotor

List of functions
##################


**-- Force extractions**

.. autosummary::

    Post.Rotor.extractSlices
    Post.Rotor.computeZb
    Post.Rotor.computeThrustAndTorque


Contents
#########

Force extractions
------------------


.. py:function:: Post.Rotor.extractSlices(teff, bladeName, psi, radius, ROINF, PINF, ASOUND, MTIP, AR, CHORD, MU, relativeShaft=0., localFrame=True, delta=0.05, accumulatorSlices=None, accumulatorCnM2=None, accumulatorCmM2=None)

    Extract slices on blades. Export Cp, Cf on those slices. 
    Compute CnM2 and CmM2 on those slices.

    :param teff: surface stress tree
    :type  teff: [zone, list of zones, base, tree]
    :param bladeName: name of the blade base to work with
    :type bladeName: string
    :param psi: angular angle of blade in teff (in degree)
    :type psi: float
    :param radius: list of radius to extract
    :type radius: list of floats
    :param ROINF: infinite flow density
    :type ROINF: float
    :param PINF: infinite flow pressure
    :type PINF: float
    :param ASOUND: infinite flow sound speed
    :type ASOUND: float
    :param MTIP: blade mach tip
    :type MTIP: float
    :param AR: blade length in m
    :type AR: float
    :param CHORD: blade chord
    :type CHORD: float
    :param MU: advance ratio
    :type MU: float
    :param relativeShaft: relative shaft angle if the mesh is not in the wind frame
    :type relativeShaft: float
    :param localFrame: if True, return CnM2 and CmM2 in relative (blade section) frame
    :type localFrame: boolean
    :param accumulatorSlices: if not None, accumulate slices
    :type accumulatorSlices: dictionary
    :param accumulatorCnM2: if not None, accumulate CnM2
    :type accumulatorCnM2: dictionary
    :param accumulatorCmM2: if not None, accumulate CnM2
    :type accumulatorCmM2: dictionary
    :return: list of slices (one for each radius)
    :rtype: list of zones

    *Example of use:*

    * `Extract slices (pyTree) <Examples/Post/extractSlicesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractSlicesPT.py

---------------------------------------

.. py:function:: Post.Rotor.computeZb(teff, psi, ROINF, ASOUND, MTIP, AR, SIGMA, relativeShaft=0., accumulatorZb=None)

    Compute Zb in the wind frame.
    
    :param teff: surface stress tree
    :type  teff: [zone, list of zones, base, tree]
    :param psi: angular angle of blade in teff (in degree)
    :type psi: float
    :param ROINF: infinite flow density
    :type ROINF: float
    :param ASOUND: infinite flow sound speed
    :type ASOUND: float
    :param MTIP: blade mach tip
    :type MTIP: float
    :param AR: blade length in m
    :type AR: float
    :param SIGMA: rotor solidity (= Nb*c / pi*AR)
    :type SIGMA: float
    :param relativeShaft: relative shaft angle if the mesh is not in the wind frame
    :type relativeShaft: float
    :param accumulatorZb: if not None, accumulate Zb
    :type accumulatorZb: dictionary    
    :return: [Xb,Yb,Zb]
    :rtype: list of 3 floats

    *Example of use:*

    * `Compute Zb (pyTree) <Examples/Post/computeZbPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeZbPT.py


---------------------------------------

.. py:function:: Post.Rotor.computeThrustAndTorque(teff, psi, PINF, center=(0,0,0), relativeShaft=0., accumulatorThrust=None)

    Compute Thrust in the rotor frame (that is orthogonal to rotor).

    :param teff: surface stress tree
    :type  teff: [zone, list of zones, base, tree]
    :param psi: angular angle of blade in teff (in degree)
    :type psi: float
    :param PINF: infinite flow pressure
    :type PINF: float
    :param center: center for momentum computations
    :type center: list of 3 floats
    :param relativeShaft: relative shaft angle if the mesh is not in the rotor frame
    :type relativeShaft: float
    :param accumulatorThrust: if not None, accumulate thrust and torque
    :type accumulatorThrust: dictionary    
    :return: thrust=[tx,ty,tz] and torque=[mx,my,mz]
    :rtype: 2 lists of 3 floats

    *Example of use:*

    * `Compute thrust and torque (pyTree) <Examples/Post/computeThrustAndTorquePT.py>`_:

    .. literalinclude:: ../build/Examples/Post/computeThrustAndTorquePT.py
