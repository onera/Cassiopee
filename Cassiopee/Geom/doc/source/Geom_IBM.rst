.. Geom IBM documentation master file

:tocdepth: 2


Geom.IBM: specific geometry modifications for IBMs 
==================================================

Specific geometry modification functions for immersed boundaries.

These functions require a geometry tree "tb" or a connectivity tree "tc".

Notes on IBCTypes
#################
Table outlining the various IBCs currently supported. Please note that the "Name" are case sensitive (e.g. Slip is not supported)

+-------------------------------------------------------+---------------+-------------------+
| IBC Type                                              | Name          | Integer Identifier|
+=======================================================+===============+===================+
| Wall slip                                             | slip          | 0                 |
+-------------------------------------------------------+---------------+-------------------+
| Wall no slip                                          | noslip        | 1                 |
+-------------------------------------------------------+---------------+-------------------+
| Wall model: Logarithmic                               | Log           | 2                 |
+-------------------------------------------------------+---------------+-------------------+
| Wall model: Musker                                    | Musker        | 3                 |
+-------------------------------------------------------+---------------+-------------------+
| Outflow pressure                                      | outpress      | 4                 |
+-------------------------------------------------------+---------------+-------------------+
| Injection                                             | inj           | 5                 |
+-------------------------------------------------------+---------------+-------------------+
| Wall model: Thin Boundary Layer Equations (TBLE)      | TBLE          | 6                 |
+-------------------------------------------------------+---------------+-------------------+
| Wall no slip with rotation (limited capabilities)     | MuskerMob     | 7                 |
+-------------------------------------------------------+---------------+-------------------+
| Wall model: Pohlhausen                                | Pohlhausen    | 8                 |
+-------------------------------------------------------+---------------+-------------------+
| Wall model: Thwaites                                  | Thwaites      | 9                 |
+-------------------------------------------------------+---------------+-------------------+
| Wall model: Mafzal                                    | Mafzal        | 10                |
+-------------------------------------------------------+---------------+-------------------+
| Wall model: Full TBLE                                 | TBLE_FULL     | 11                |
+-------------------------------------------------------+---------------+-------------------+
| Wall no slip with curvature radius                    | slip_cr       | 100               |
+-------------------------------------------------------+---------------+-------------------+

Note: Wall no slip with rotation with an integer identifier of 7 assumes the x-axis (i.e. x=0) is the axis of rotation. A more complete IBM with rotation is in the process of being developed and will be shortly be available making the current IBC decrepit.

.. py:module:: Geom.IBM

List of functions
#################


**-- Setting Snear & Dfar**

.. autosummary::
   :nosignatures:

    Geom.IBM.setSnear
    Geom.IBM.setDfar
    Geom.IBM.snearFactor
    .. Geom.IBM.getDfarOpt
    
**-- Setting IBC Type**

.. autosummary::
   :nosignatures:

    Geom.IBM.setIBCType
    Geom.IBM.changeIBCType
    Geom.IBM.initOutflow
    Geom.IBM.initInj
    Geom.IBM.setFluidInside


Contents
########

Note that all the functions have an in-place version, modifying directly the data without copy.
The function names must be prefixed by an '_' (e.g. _setSnear for the in-place version of setSnear)

Setting Snear & Dfar
--------------------


.. py:function:: Geom.IBM.setSnear(tb, snear)

    Set the snear for a geometry defined by tb. Exists also as in-place (_setSnear). Snear is the local Cartesian spacing close to  cells intersected by the immersed boundary.

    :param tb: geometry tree
    :type  tb: [zone, list of zones, tree]
    :param snear: snear value
    :type snear: float
    :return: same as input

    *Example of use:*
    
    * `Set the values of the snears (pyTree) <Examples/Geom/setSnearPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/setSnearPT.py

---------------------------------------


.. py:function:: Geom.IBM.setDfar(tb, dfar)

    Set the dfar for a geometry defined by tb. Exists also as in-place (_setDfar). Dfar is the distance from the center of the bounding box of the immersed boundary to the edge of the domain.

    :param tb: geometry tree
    :type  tb: [zone, list of zones, tree]
    :param dfar: dfar value
    :type dfar: float    
    :return: same as input

    *Example of use:*
    
    * `Set the value of the dfar (pyTree) <Examples/Geom/setDfarPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/setDfarPT.py

---------------------------------------



    .. .. py:function:: Geom.IBM.getDfarOpt(tb, vmin, snear, factor=10, nlevel=-1)

       Computes the optimal dfar to get the exact snear.

       :param tb: geometry tree
       :type  tb: [zone, list of zones, tree]
       :param vmin: number of points per elementary structured octree zone
       :type vmin: integer
       :param snear: targeted snear (smallest value)
       :type snear: float
       :param factor: factor*L is the extent of the octree, where L is the maximum length of the bounding box of tb.
       :type factor: float
       :param nlevel: number of levels prescribed for the octree (-1 means that the number of refinement levels are computed according to the parameter factor).
       :type nlevel: integer
       :return: float

       *Example of use:*

       * `Get the value of the dfar to respect exactly snear (pyTree) <Examples/Geom/getDfarOptPT.py>`_:

       .. literalinclude:: ../build/Examples/Geom/getDfarOptPT.py

---------------------------------------


.. py:function:: Geom.IBM.snearFactor(tb, sfactor)

    Multiply the snears in the geometry defined by tb by a factor. Exists also as in-place (_snearFactor).

    :param tb: geometry tree
    :type  tb: [zone, list of zones, tree]
    :param sfactor: multiplying factor
    :type sfactor: float
    :return: same as input

    *Example of use:*
    
    * `Modifying the value of snears (pyTree) <Examples/Geom/snearFactorPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/snearFactorPT.py
			


Setting IBC Type
----------------

.. py:function:: Geom.IBM.setIBCType(tb, ibctype)

    Set the type of IBC for the geometry defined by tb. Exists also as in-place (_setIBCType). See the table in "Notes on IBCTypes" for the IBCs currently supported.

    :param tb: geometry tree
    :type  tb: [zone, list of zones, tree]
    :param ibctype: name of the type of IBC
    :type ibctype: string
    :return: same as input

    *Example of use:*
    
    * `Set the type of IBC (pyTree) <Examples/Geom/setIBCTypePT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/setIBCTypePT.py
			
---------------------------------------

.. py:function:: Geom.IBM.changeIBCType(tc,oldBCType,newBCType)

    Change the IBC type in a connectivity tree. Exists also as in-place (_changeIBCType). Please refer to the table in "Notes on IBCTypes" for details on the integer identifies for the various IBC types.

    :param tc: connectivity tree
    :type  tc: [zone, list of zones, tree]
    :param oldBCType: type of ibc
    :type oldBCType: integer
    :param newBCType: type of ibc
    :type newBCType: integer
    :return: same as input

    *Example of use:*
    
    * `Change the type of IBC (pyTree) <Examples/Geom/changeIBCTypePT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/changeIBCTypePT.py

---------------------------------------

.. py:function:: Geom.IBM.setFluidInside(tb)

    Define the fluid inside a surface defined by tb. In that case, the IBM mesh will be defined inside tb. Exists also as in-place (_setFluidInside).

    :param tb: geometry tree 
    :type  tb: [zone, list of zones, tree]
    :return: same as input

    *Example of use:*
    
    * `Define fluid inside surface (pyTree) <Examples/Geom/setFluidInsidePT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/setFluidInsidePT.py

---------------------------------------

.. py:function:: Geom.IBM.symmetrizePb(tb)

    Add a symmetry plane and symmetrize the body base; input surface must be y>=0 if the symmetry plane is at y=0.
    Creates the geometry for the symmetry plane, the snear and slip IBC type.
    :param tb: geometry tree 
    :type  tb: [zone, list of zones, tree]
    :param bodyNameSym: name of the base where
    :type bodyNameSym: string
    :param snear_sym: mesh resolution at the symmetry plane
    :type snear_sym: float
    :param dir_sym: direction of the symmetry plane
    :type dir_sym: integer (1:x,2:y,3:z)                    
    :return: same as input

    *Example of use:*
    
    * `Symmetrize the surface (pyTree) <Examples/Geom/symetrizePbPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/symetrizePbPT.py
                                
---------------------------------------

.. py:function:: Geom.IBM.initOutflow(tc, familyName, Pstatic, InterpolPlane=None, PressureVar=0,isDensityConstant=False)

    Set the value of the static pressure Pstatic for the outflow pressure IBC with family name familyName. A plane InterpolPlane may also be provided with only static pressure variable or various variables with static pressure as the PressureVar (e.g. 2nd) variable).
    Exists also as in-place (_initOutflow).
    
    :param tc: connectivity tree
    :type  tc: [zone, list of zones, tree]
    :param familyName: familyName
    :type familyName: string
    :param Pstatic: static pressure
    :type Pstatic: float
    :param PInterpolPlane: interpolation plane
    :type PInterpolPlane: [zone, list of zones, tree]
    :param PressureVar: variable number of static pressure
    :type PressureVar: integer
    :param isDensityConstant: is density at the outflow constant and equal to the reference density
    :type isDensityConstant: boolean
    :return: same as input

    *Example of use:*
    
    * `Set outflow IBC (pyTree) <Examples/Geom/initOutflowPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/initOutflowPT.py


---------------------------------------

.. py:function:: Geom.IBM.initInj(tc, familyName, Ptot, Htot, injDir=[1.,0.,0.], InterpolPlane=None, PressureVar=0, EnthalpyVar=0)

    Set the total pressure Ptot, total enthalpy Htot, and direction of the flow injDir (w.r.t the absolute frame of reference) for the injection IBC with family name familyName. A plane InterpolPlane may also be provided with at least the total pressure and total enthalpy variables with the former and latter as the PressureVar (e.g. 2nd) and EnthalpyVar (e.g. 4th) variables, respectively.
    Exists also as in-place (_initInj). 

    :param tc: connectivity tree
    :type  tc: [zone, list of zones, tree]
    :param familyName: familyName
    :type familyName: string
    :param Ptot: total pressure
    :type Ptot: float
    :param Htot: total enthalpy
    :type Htot: float
    :param injDir: direction of the injection w.r.t to the reference coordinate axis
    :type injDir: float list
    :param PInterpolPlane: interpolation plane
    :type PInterpolPlane: [zone, list of zones, tree]
    :param PressureVar: variable number of total pressure			  
    :type PressureVar: integer
    :param EnthalpyVar: variable number of total enthalpy		       
    :type EnthalpyVar: integer
    :return: same as input

    *Example of use:*
    
    * `Set injection IBC (pyTree) <Examples/Geom/initInjPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Geom/initInjPT.py   

---------------------------------------


.. toctree::
   :maxdepth: 2   


Index
#######

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
		
