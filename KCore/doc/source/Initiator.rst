.. Initiator documentation master file

Initiator: solution Initialization Module
=========================================

Preamble
########

Initiator module works on arrays (as defined in Converter) or on CGNS/python trees (pyTrees) 
containing grid information (coordinates must be defined).

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use with the array interface, you have to import Initiator module::

   import Initiator as I

For use with the pyTree interface::

    import Initiator.PyTree as I


.. py:module:: Initiator

List of functions
##################

**-- CFD field initialisations**

.. autosummary::

   initConst
   initLamb
   initVisbal
   initScully
   initYee
   overlayField

**-- Adimensioning**

.. autosummary::

   Adim.adim1
   Adim.adim2
   Adim.adim3
   Adim.dim1
   Adim.dim2
   Adim.dim3


Contents
#########

CFD field initialisations
--------------------------

The input data defines a grid or a set of grids on which the solution has to be initialized.

The created field variables are the variables defined in the input data. 
If the five conservative variables are not present, then the default output variables are the coordinates and the conservative variables.

---------------------------------------

.. py:function:: Initiator.initConst(a, adim='adim1', MInf=None, alphaZ=0., alphaY=0., ReInf=1.e8, loc='nodes')

    Initialization by a constant field, given Mach number, incident flow angles and Reynolds number.
    
    Exists also as in place version (_initConst) that modifies a and returns None.

    :param a:  input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param adim:   Name of adim ('adim1', 'adim2', 'dim1', ...) - see Adimensioning section
    :type  adim:   string
    :param MInf:   freestream Mach number
    :type  MInf:   float
    :param alphaZ: Angle with respect to (0,Z) axe (in degrees)
    :type  alphaZ: float
    :param alphaY: Angle with respect to (0,Y) axe (in degrees)
    :type  alphaY: float
    :param ReInf:  Reynolds Number
    :type  ReInf:  float
    :param loc: created field localisation ('nodes' or 'centers') - only for pyTree interface
    :type loc: string
    :rtype: Identical to a

    *Example of use:*

    * `Constant field initialization (array) <Examples/Initiator/initConst.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/initConst.py

    * `Constant field initialization (pyTree) <Examples/Initiator/initConstPT.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/initConstPT.py

---------------------------------------

.. py:function:: Initiator.initLamb(a, (x0,y0), Gamma=2., MInf=0.5, loc='nodes')

    Initialization of conservative variables by a 2D Lamb vortex at position (x0,y0), intensity Gamma
    and infinite Mach number MInf.
    
    Exists also as in place version (_initLamb) that modifies a and returns None.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param Gamma:  Intensity of vortex
    :type  Gamma:   float
    :param MInf:   Mach number
    :type  MInf:   float
    :param loc: created field localisation ('nodes' or 'centers')  - only for pyTree interface
    :rtype: Identical to a

    *Example of use:*

    * `Field initialization with a lamb vortex (array) <Examples/Initiator/lamb.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/lamb.py

    * `Field initialization with a lamb vortex (pyTree) <Examples/Initiator/lambPT.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/lambPT.py

---------------------------------------

.. py:function:: Initiator.initVisbal(a, (x0,y0), Gamma=2., MInf=0.5, loc='nodes')

    Initialization of conservative variables by a 2D Visbal vortex at position (x0,y0), intensity Gamma
    and infinite Mach number MInf.
    
    Exists also as in place version (_initVisbal) that modifies a and returns None.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param Gamma:  Intensity of vortex
    :type  Gamma:   float
    :param MInf:   Mach number
    :type  MInf:   float
    :param loc: field localisation ('nodes' or 'centers')  - only for pyTree interface
    :rtype: Identical to a

    *Example of use:*

    * `Field initialization with a visbal vortex (array) <Examples/Initiator/visbal.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/visbal.py

    * `Field initialization with a visbal vortex (pyTree) <Examples/Initiator/visbalPT.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/visbalPT.py

---------------------------------------

.. py:function:: Initiator.initScully(a, (x0,y0), Gamma=2., coreRadius=1., MInf=0.5, loc='nodes')

    Initialization of conservative variables by a 2D Scully vortex at position (x0,y0), intensity Gamma
    and infinite Mach number MInf.
    
    Exists also as in place version (_initScully) that modifies a and returns None.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param Gamma:  Intensity of vortex
    :type  Gamma:   float
    :param coreRadius: radius of vortex core
    :type coreRadius: float
    :param MInf:   Mach number
    :type  MInf:   float
    :param loc: field localisation ('nodes' or 'centers')  - only for pyTree interface
    :rtype: Identical to a

    *Example of use:*

    * `Field initialization with a scully vortex (array) <Examples/Initiator/scully.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/scully.py

    * `Field initialization with a scully vortex (pyTree) <Examples/Initiator/scullyPT.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/scullyPT.py

---------------------------------------

.. py:function:: Initiator.initYee(a, (x0,y0), Gamma=2., MInf=0.5, loc='nodes')

    Initialization of conservative variables by a 2D Yee vortex at position (x0,y0), intensity Gamma
    and infinite Mach number MInf.
    
    Exists also as in place version (_initYee) that modifies a and returns None.

    :param a:  Input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param Gamma:  Intensity of vortex
    :type  Gamma:   float
    :param MInf:   Mach number
    :type  MInf:   float
    :param loc: field localisation ('nodes' or 'centers')  - only for pyTree interface
    :rtype: Identical to a
 
    *Example of use:*

    * `Field initialization with a Yee vortex (array) <Examples/Initiator/yee.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/yee.py

    * `Field initialization with a Yee vortex (pyTree) <Examples/Initiator/yeePT.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/yeePT.py

---------------------------------------

.. py:function:: Initiator.overlayField(a, b, MInf=0.5, loc='nodes')

    Overlay the field of two solutions defined on the same grid (defined in a and b). 
    Only density, velocity and energy stagnation density are overlaid. 
    Both fields must use the same adimensioning, which corresponds to the 
    same infinite Mach number MInf. 
    
    Exists also as in place version (_overlayField) that modifies a and returns None.

    :param a:  First input data
    :type  a:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param b:  Second input data
    :type  b:  [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param MInf:   Mach number
    :type  MInf:   float
    :param loc: field localisation ('nodes' or 'centers')  - only for pyTree interface
    :rtype: Identical to a

    *Example of use:*

    * `Field initialization by overlaying two solutions (array) <Examples/Initiator/overlay.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/overlay.py

    * `Field initialization by overlaying two solutions (pyTree) <Examples/Initiator/overlayPT.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/overlayPT.py


Adimensioning
--------------


.. py:function:: Initiator.Adim.adim1(MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8, MutSMuInf=0.2, TurbLevelInf=1.e-4)

    Return a reference state corresponding to an adimensioning by density, sound speed and temperature.
    Returned state is adimensioned. When using this state, the mesh must also be adimensioned.

    :param MInf:   Mach number
    :type  MInf:   float
    :param alphaZ: Angle with respect to (0,Z) axe (in degrees)
    :type  alphaZ: float
    :param alphaY: Angle with respect to (0,Y) axe (in degrees)
    :type  alphaY: float
    :param ReInf:  Reynolds Number
    :type  ReInf:  float
    :param MutSMuInf: ratio of mut/mu (turbulence viscosity/molecular viscosity) at infinity 
    :type MutSMuInf: float
    :param TurbLevelInf: turbulence level at infinity
    :return: [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
              ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
              Mus, Cs, Ts, Pr]
    :rtype: list

    *Example of use:*

    * `Get adimensioned state <Examples/Initiator/adim1.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/adim1.py

    .. note:: New in version 2.5

---------------------------------------

.. py:function:: Initiator.Adim.adim2(MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8, MutSMuInf=0.2, TurbLevelInf=1.e-4)

    Return a reference state corresponding to an adimensioning by density, fluid velocity and temperature.
    Returned state is adimensioned. When using this state, the mesh must also be adimensioned.

    :param MInf:   Mach number
    :type  MInf:   float
    :param alphaZ: Angle with respect to (0,Z) axe (in degrees)
    :type  alphaZ: float
    :param alphaY: Angle with respect to (0,Y) axe (in degrees)
    :type  alphaY: float
    :param ReInf:  Reynolds Number
    :type  ReInf:  float
    :param MutSMuInf: ratio of mut/mu (turbulence viscosity/molecular viscosity) at infinity 
    :type MutSMuInf: float
    :param TurbLevelInf: turbulence level at infinity
    :rtype: [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

    *Example of use:*

    * `Get adimensioned state 2 <Examples/Initiator/adim2.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/adim2.py

    .. note:: New in version 2.5

---------------------------------------

.. py:function:: Initiator.Adim.adim3(MInf=0.5, alphaZ=0., alphaY=0., ReInf=1.e8, LInf=1., MutSMuInf=0.2, TurbLevelInf=1.e-4)

    Return a reference state corresponding to an adimensioning by density, sound speed and temperature.
    Returned state is adimensioned. When using this state, the mesh doesn't need to be adimensioned and
    has a characteritic length of LInf.

    :param MInf:   Mach number
    :type  MInf:   float
    :param alphaZ: Angle with respect to (0,Z) axe (in degrees)
    :type  alphaZ: float
    :param alphaY: Angle with respect to (0,Y) axe (in degrees)
    :type  alphaY: float
    :param ReInf:  Reynolds Number
    :type  ReInf:  float
    :param LInf:  Characteristic length of mesh (on which ReInf is based)
    :type  LInf:  float
    :param MutSMuInf: ratio of mut/mu (turbulence viscosity/molecular viscosity) at infinity 
    :type MutSMuInf: float
    :param TurbLevelInf: turbulence level at infinity
    :rtype: [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

    *Example of use:*

    * `Get adimensioned state 3 <Examples/Initiator/adim3.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/adim3.py

    .. note:: New in version 2.5

---------------------------------------

.. py:function:: Initiator.Adim.dim1(UInf=2.7777, TInf=298.15, PInf=101325., LInf=1., alphaZ=0., alphaY=0., MutSMuInf=0.2, TurbLevelInf=1.e-4)

    Return a dimensioned reference state corresponding to dry air, considered as a perfect gaz.
    Returned a dimensioned state (USI).

    :param UInf:   Fluid velocity in m/s
    :type  UInf:   float
    :param TInf:   Temperature in K (0 degree=273.15K)
    :type  TInf:   float
    :param PInf:   Pressure in Pa
    :type  PInf:   float
    :param LInf:    Reference length (m). Usefull to compute Reynolds.
    :type  LInf:   float
    :param MutSMuInf: ratio of mut/mu (turbulence viscosity/molecular viscosity) at infinity 
    :type MutSMuInf: float
    :param TurbLevelInf: turbulence level at infinity
    :rtype: [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

    *Example of use:*

    * `Get dimensioned state <Examples/Initiator/dim1.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/dim1.py

    .. note:: New in version 2.5

---------------------------------------

.. py:function:: Initiator.Adim.dim2(UInf=2.7777, TInf=298.15, RoInf=1.225, LInf=1., alphaZ=0., alphaY=0., MutSMuInf=0.2, TurbLevelInf=1.e-4)

    Return a dimensioned reference state corresponding to dry air, considered as a perfect gaz.
    Only the input is different from dim1.
    Returned a dimensioned state (USI).

    :param UInf:   Fluid velocity in m/s
    :type  UInf:   float
    :param TInf:   Temperature in K (0 degree=273.15K)
    :type  TInf:   float
    :param RoInf:  Input density (kg/m3)
    :type  RoInf:   float
    :param LInf:    Reference length (m). Usefull to compute Reynolds.
    :type  LInf:   float
    :param MutSMuInf: ratio of mut/mu (turbulence viscosity/molecular viscosity) at infinity 
    :type MutSMuInf: float
    :param TurbLevelInf: turbulence level at infinity
    :rtype: [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

    *Example of use:*

    * `Get dimensioned state 2 <Examples/Initiator/dim2.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/dim2.py

    .. note:: New in version 2.5

---------------------------------------

.. py:function:: Initiator.Adim.dim3(UInf=2.7777, PInf=101325., RoInf=1.225, LInf=1., alphaZ=0., alphaY=0., MutSMuInf=0.2, TurbLevelInf=1.e-4)

    Return a dimensioned reference state corresponding to dry air, considered as a perfect gaz.
    Only the input is different from dim1.
    Returned a dimensioned state (USI).

    :param UInf:   Fluid velocity in m/s
    :type  UInf:   float
    :param PInf:   Pressure in Pa
    :type  PInf:   float
    :param RoInf:  Input density (kg/m3)
    :type  RoInf:   float
    :param LInf:    Reference length (m). Usefull to compute Reynolds.
    :type  LInf:   float
    :param MutSMuInf: ratio of mut/mu (turbulence viscosity/molecular viscosity) at infinity 
    :type MutSMuInf: float
    :param TurbLevelInf: turbulence level at infinity
    :rtype: [RoInf, RouInf, RovInf, RowInf, RoEInf, PInf, TInf, cvInf, MInf,
            ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
            Mus, Cs, Ts, Pr]

    *Example of use:*

    * `Get dimensioned state 3 <Examples/Initiator/dim3.py>`_:

    .. literalinclude:: ../build/Examples/Initiator/dim3.py

    .. note:: New in version 2.5 

.. toctree::
   :maxdepth: 2   


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

