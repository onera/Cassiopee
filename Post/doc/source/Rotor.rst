.. Post documentation master file

:tocdepth: 2


Rotor: rotor specific post-processing 
======================================

Specific post-processing for rotors or propellers.

Those functions work with a solution tree "t" or a blade surface tree "teff".


.. py:module:: Post.Rotor

List of functions
##################


**-- Force extractions**

.. autosummary::

    Post.Rotor.extractSlices
    Post.Rotor.computeZb
    Post.Rotor.computeThrustAndTorque

**-- Accumulator export**

.. autosummary::

    Post.Rotor.exportAccumulatorPerPsi
    Post.Rotor.exportAccumulatorPerRadius
    Post.Rotor.exportAccumulatorMap


Contents
#########

Force extractions
------------------


.. py:function:: Post.Rotor.extractSlices(teff, bladeName, psi, radius, RoInf, PInf, ASOUND, Mtip, AR, CHORD, MU, relativeShaft=0., localFrame=True, delta=0.05, accumulatorSlices=None, accumulatorCnM2=None, accumulatorCmM2=None)

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
    :param RoInf: infinite flow density
    :type RoInf: float
    :param PInf: infinite flow pressure
    :type PInf: float
    :param ASOUND: infinite flow sound speed
    :type ASOUND: float
    :param Mtip: blade mach tip
    :type Mtip: float
    :param AR: blade length
    :type AR: float
    :param CHORD: blade mean chord
    :type CHORD: float
    :param MU: advance ratio
    :type MU: float
    :param relativeShaft: relative shaft angle if the mesh is not in the wind frame
    :type relativeShaft: float
    :param localFrame: if True, return CnM2 and CmM2 in relative (blade section) frame
    :type localFrame: boolean
    :param delta: mean mesh step on blade in the span wise direction
    :type delta: float
    :param accumulatorSlices: if not None, accumulate slices
    :type accumulatorSlices: dictionary
    :param accumulatorCnM2: if not None, accumulate CnM2
    :type accumulatorCnM2: dictionary
    :param accumulatorCmM2: if not None, accumulate CmM2
    :type accumulatorCmM2: dictionary
    :return: list of slices, list of CnM2, list of CmM2 (one for each radius)
    :rtype: list of zones, list of 3 floats, list of 3 floats

    *Example of use:*

    * `Extract slices (pyTree) <Examples/Post/extractSlicesPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/extractSlicesPT.py

---------------------------------------

.. py:function:: Post.Rotor.computeZb(teff, psi, RoInf, ASOUND, Mtip, AR, SIGMA, relativeShaft=0., accumulatorZb=None)

    Compute Zb in the wind frame.
    
    :param teff: surface stress tree
    :type  teff: [zone, list of zones, base, tree]
    :param psi: angular angle of blade in teff (in degree)
    :type psi: float
    :param RoInf: infinite flow density
    :type RoInf: float
    :param ASOUND: infinite flow sound speed
    :type ASOUND: float
    :param Mtip: blade mach tip
    :type Mtip: float
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

.. py:function:: Post.Rotor.computeThrustAndTorque(teff, psi, PInf, center=(0,0,0), relativeShaft=0., accumulatorThrust=None)

    Compute Thrust in the rotor frame (that is orthogonal to rotor).

    :param teff: surface stress tree
    :type  teff: [zone, list of zones, base, tree]
    :param psi: angular angle of blade in teff (in degree)
    :type psi: float
    :param PInf: infinite flow pressure
    :type PInf: float
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

Accumulator export
-------------------

.. py:function:: Post.Rotor.exportAccumulatorPerPsi(accumulator, psi=0., vars=['F1','F2'])

    Export a given psi of an accumulator (psi,rad) in a 1D zone.
    For distributed computations, the exported zone is identical on all processors.

    :param accumulator: (psi,rad) accumulator
    :type  accumulator: dictionary
    :param psi: angular angle to be extracted (in degree)
    :type psi: float
    :param vars: the name of variables stored in accumulator
    :type vars: list of strings
    :return: a single Zone with vars corresponding to psi
    :rtype: Zone

    *Example of use:*

    * `Export accumulator for given psi (pyTree) <Examples/Post/exportAccumulatorPerPsiPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/exportAccumulatorPerPsiPT.py

---------------------------------------

.. py:function:: Post.Rotor.exportAccumulatorPerRadius(accumulator, rad=0., vars=['F1','F2'])

    Export a given radius of an accumulator (psi,rad) in a 1D zone.
    For distributed computations, the exported zone is identical on all processors.

    :param accumulator: (psi,rad) accumulator
    :type  accumulator: dictionary
    :param rad: radius to be extracted
    :type rad: float
    :param vars: the name of variables stored in accumulator
    :type vars: list of strings
    :return: a single Zone with vars corresponding to rad
    :rtype: Zone

    *Example of use:*

    * `Export accumulator for given rad (pyTree) <Examples/Post/exportAccumulatorPerRadiusPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/exportAccumulatorPerRadiusPT.py


---------------------------------------

.. py:function:: Post.Rotor.exportAccumulatorMap(accumulator, vars=['Fx','Fy','Fz'])

    Export accumulator (psi,rad) to a 2D zone.
    For distributed computations, the exported zone is identical on all processors.

    :param accumulator: (psi,rad) accumulator
    :type  accumulator: dictionary
    :param vars: the name of variables stored in accumulator
    :type vars: list of strings
    :return: a single Zone with fields
    :rtype: Zone

    *Example of use:*

    * `Export accumulator to a map (pyTree) <Examples/Post/exportAccumulatorMapPT.py>`_:

    .. literalinclude:: ../build/Examples/Post/exportAccumulatorMapPT.py
