.. Connector IBM documentation master file

:tocdepth: 2


Connector.IBM: immersed boundary method grid connectivity module
================================================================

.. py:module:: Connector.IBM


List of functions
#################


**-- IBM Connectivity**

.. autosummary::
   :nosignatures:

   <!--Connector.IBM._computeFrictionVelocity-->
   <!--Connector.IBM._blankClosestTargetCells-->
   <!--Connector.IBM._removeBlankedGrids-->
   <!--Connector.IBM.blankByIBCBodies-->
   <!--Connector.IBM._addBCOverlaps-->
   <!--Connector.IBM._addExternalBCs-->
   <!--Connector.IBM.getIBMFront-->
   <!--Connector.IBM._pushBackImageFront2-->
   <!--Connector.IBM._smoothImageFront-->
   <!--Connector.IBM._smoothImageFrontBackward-->
   <!--Connector.IBM.gatherFront-->
   <!--Connector.IBM.doInterp-->
   <!--Connector.IBM._extractIBMInfo_param-->
   <!--Connector.IBM.extractIBMInfo-->
   <!--Connector.IBM.extractIBMInfo2-->

   <!--Connector.IBM.getAllIBMPoints-->
   <!--Connector.IBM.prepareIBMData-->
   <!--Connector.IBM.createWallAdapt-->
   <!--Connector.IBM.createIBMWZones-->
   <!--Connector.IBM._computeKcurvParameter-->
   <!--Connector.IBM._signDistance-->

   Connector.IBM.dist2wallIBM
   Connector.IBM.blankingIBM
   Connector.IBM.buildFrontIBM
   Connector.IBM.setInterpDataIBM
   Connector.IBM.initializeIBM

Contents
###########

Main functions
--------------

.. py:function:: Connector.IBM.dist2wallIBM(t, tb, dimPb=3, frontType=1, Reynolds=1.e6, yplus=100, Lref=1., correctionMultiCorpsF42=False, heightMaxF42=-1.)

    Compute the wall distance for IBM pre-processing. Exists also as in-place (_dist2wallIBM). 
    
    The Reynolds number, yplus and Lref are used to calculate a custom modeling height when using frontType 42.

    :param t: computational tree
    :type t: tree
    :param tb: geometry tree
    :type tb: tree
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param frontType: type of IBM front
    :type frontType: 0, 1, 2 or 42
    :param Reynolds: Reynolds number (frontType 42)
    :type Reynolds: float
    :param yplus: estimated yplus at the first computed cells (frontType 42)
    :type yplus: float
    :param Lref: characteristic length of the geometry (frontType 42)
    :type Lref: float
    :param correctionMultiCorpsF42: if True, computes the wall distance w.r.t each body that is not a symmetry plane (frontType 42)
    :type correctionMultiCorpsF42: boolean
    :param heightMaxF42: if heightMaxF42 > 0: uses a maximum modeling height to speed up individual wall distance calculations when correctionMultiCorpsF42 is active (frontType 42)
    :type heightMaxF42: float

    *Example of use:*
    
    * `Compute the wall distance for IBM (pyTree) <Examples/Connector/dist2wallIBMPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Connector/dist2wallIBMPT.py

---------------------------------------

.. py:function:: Connector.IBM.blankingIBM(t, tb, dimPb=3, frontType=1, IBCType=1, depth=2, Reynolds=1.e6, yplus=100, Lref=1., twoFronts=False, correctionMultiCorpsF42=False, blankingF42=False, wallAdaptF42=None, heightMaxF42=-1.)

    Blank the computational tree by IBC bodies for IBM pre-processing. Exists also as in-place (_blankingIBM). 
    
    The Reynolds number, yplus and Lref are used to calculate a custom modeling height when using frontType 42. 
    
    The wallAdaptF42 file must be obtained with Connector.IBM.createWallAdapt().

    :param t: computational tree
    :type t: tree
    :param tb: geometry tree
    :type tb: tree
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param frontType: type of IBM front
    :type frontType: 0, 1, 2 or 42
    :param IBCType: type of IBM, -1: IB target points are located inside the solid, 1: IB target points are located in the fluid
    :type IBCType: -1 or 1
    :param depth: depth of overlapping regions
    :type depth: int
    :param Reynolds: Reynolds number (frontType 42)
    :type Reynolds: float
    :param yplus: estimated yplus at the first computed cells (frontType 42)
    :type yplus: float
    :param Lref: characteristic length of the geometry (frontType 42)
    :type Lref: float
    :param twoFronts: if True, performs the IBM pre-processing for an additional image point positioned farther away
    :type twoFronts: boolean
    :param correctionMultiCorpsF42: if True, ensures that there are calculated points between the immersed bodies by using individual wall distances (frontType 42)
    :type correctionMultiCorpsF42: boolean
    :param blankingF42: if True, reduces as much as possible the number of IB target points inside the boundary layer (frontType 42)
    :type blankingF42: boolean
    :param wallAdaptF42: use a previous computation to adapt the positioning of IB target points around the geometry according to a target yplus (frontType 42)
    :type wallAdaptF42: cloud of IB target points with yplus information
    :param heightMaxF42: if heightMaxF42 > 0: maximum modeling height for the location of IB target points around the geometry (frontType 42)
    :type heightMaxF42: float

    *Example of use:*
    
    * `Blanking for IBM (pyTree) <Examples/Connector/blankingIBMPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Connector/blankingIBMPT.py

---------------------------------------
 
.. py:function:: Connector.IBM.buildFrontIBM(t, tc, dimPb=3, frontType=1, cartesian=False, twoFronts=False, check=False)

    Build the IBM front for IBM pre-processing.

    :param t: computational tree
    :type t: tree
    :param tc: connectivity tree
    :type tc: tree
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param frontType: type of IBM front
    :type frontType: 0, 1, 2 or 42
    :param cartesian: if True, activates optimized algorithms for Cartesian meshes
    :type cartesian: boolean
    :param twoFronts: if True, performs the IBM pre-processing for an additional image point positioned farther away
    :type twoFronts: boolean
    :param check: if True, saves front.cgns (and front2.cgns if twoFronts is active)
    :type check: boolean

    *Example of use:*
    
    * `Build the IBM front (pyTree) <Examples/Connector/buildFrontIBMPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Connector/buildFrontIBMPT.py

---------------------------------------
 
.. py:function:: Connector.IBM.setInterpDataIBM(t, tc, tb, front, front2=None, dimPb=3, frontType=1, IBCType=1, depth=2, Reynolds=1.e6, yplus=100, Lref=1., cartesian=False, twoFronts=False)

    Compute the transfer coefficients and data for IBM pre-processing. The information are stored in the connectivity tree (IBCD* zones). Exists also as in-place (_setInterpDataIBM). 

    The Reynolds number, yplus and Lref are used to calculate a custom modeling height when using frontType 42. 

    front and front2 must be obtained with Connector.IBM.buildFrontIBM().

    :param t: computational tree
    :type t: tree
    :param tc: connectivity tree
    :type tc: tree
    :param front: front of image points
    :type front: tree
    :param front2: front of second image points (optional)
    :type front2: tree
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param frontType: type of IBM front
    :type frontType: 0, 1, 2 or 42
    :param IBCType: type of IBM, -1: IB target points are located inside the solid, 1: IB target points are located in the fluid
    :type IBCType: -1 or 1
    :param depth: depth of overlapping regions
    :type depth: int
    :param Reynolds: Reynolds number (frontType 42)
    :type Reynolds: float
    :param yplus: estimated yplus at the first computed cells (frontType 42)
    :type yplus: float
    :param Lref: characteristic length of the geometry (frontType 42)
    :type Lref: float
    :param cartesian: if True, activates optimized algorithms for Cartesian meshes
    :type cartesian: boolean
    :param twoFronts: if True, performs the IBM pre-processing for an additional image point positioned farther away
    :type twoFronts: boolean

    *Example of use:*
    
    * `Compute IBM coefficients (pyTree) <Examples/Connector/setInterpDataIBMPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Connector/setInterpDataIBMPT.py

---------------------------------------

.. py:function:: Connector.IBM.initializeIBM(t, tc, tb, tinit=None, dimPb=3, twoFronts=False)

    Initialize the computational and connectivity trees for IBM pre-processing.

    tinit might be used to initialize the flow solution in t.
    
    :param t: computational tree
    :type t: tree
    :param tc: connectivity tree
    :type tc: tree
    :param tb: geometry tree
    :type tb: tree
    :param tinit: computational tree from previous computation
    :type tinit: tree
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param twoFronts: if True, creates a new connectivity tree that contains second image points information
    :type twoFronts: boolean