.. Connector IBM documentation master file

:tocdepth: 2


Connector.IBM: immersed boundary method grid connectivity module
================================================================

.. py:module:: Connector.IBM


List of functions
#################


**-- IBM Connectivity**

.. autosummary::

   Connector.IBM._computeFrictionVelocity
   Connector.IBM._blankClosestTargetCells
   Connector.IBM._removeBlankedGrids
   Connector.IBM.blankByIBCBodies
   Connector.IBM._addBCOverlaps
   Connector.IBM._addExternalBCs
   Connector.IBM._modifPhysicalBCs__
   Connector.IBM.getIBMFront
   Connector.IBM._pushBackImageFront2
   Connector.IBM._smoothImageFront
   Connector.IBM._smoothImageFrontBackward
   Connector.IBM.gatherFront
   Connector.IBM.doInterp
   Connector.IBM.doInterp2
   Connector.IBM.doInterp3
   Connector.IBM._extractIBMInfo_param
   Connector.IBM.extractIBMInfo
   Connector.IBM.extractIBMInfo2

   Connector.IBM.getAllIBMPoints
   Connector.IBM.prepareIBMData
   Connector.IBM.prepareIBMData2
   Connector.IBM.createWallAdapt
   Connector.IBM.createIBMWZones
   Connector.IBM._computeKcurvParameter
   Connector.IBM._signDistance

   Connector.IBM.dist2wallIBM
   Connector.IBM.blankingIBM
   Connector.IBM.buildFrontIBM
   Connector.IBM.setInterpDataIBM
   

Contents
###########

Main functions
--------------

.. py:function:: Connector.IBM._dist2wallIBM(t, tb, dimPb=3, correctionMultiCorpsF42=False, frontType=1, yplus=100, Reynolds=1.e6, Lref=1., heightMaxF42=-1.)

    Compute wall distance for IBM pre-processing.

    :param t: computation tree
    :type t: tree
    :param tb: immersed body geometry
    :type tb: tree
    :param dimPb: problem dimension
    :type dimPb: 2 or 3
    :param correctionMultiCorpsF42: ??
    :type correctionMultiCorpsF42: boolean
    :param frontType: type of IBM front
    :type frontType: 1, 2, 42
    :param yplus: desired yplus (front 42)
    :type yplus: float
    :param Reynolds: Reynolds on body
    :type Reynolds: float
    :param Lref: Reference length of body
    :type Lref: float
    :param heightMaxF42: max height (front42)
    :type heightMaxF42: float

    *Example of use:*
    
    * `Compute wall distance for IBM (pyTree) <Examples/Connector/dist2wallIBMPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Connector/dist2wallIBMPT.py

---------------------------------------

.. py:function:: Connector.IBM._blankingIBM(t, tb, dimPb=3, frontType=1, IBCType=1, DEPTH=2, yplus=100, Reynolds=1.e6, Lref=1., heightMaxF42=-1., correctionMultiCorpsF42=False, wallAdaptF42=None, blankingF42=False, twoFronts=False)

    Blank t by bodies tb for IBM pre-processing.

    *Example of use:*
    
    * `Blanking for IBM (pyTree) <Examples/Connector/blankingIBMPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Connector/blankingIBMPT.py

---------------------------------------
 
.. py:function:: Connector.IBM.buildFrontIBM(t, tc, dimPb=3, frontType=1, interpDataType=0, cartesian=False, twoFronts=False, check=False)

    Build the IBM front.

    *Example of use:*
    
    * `Build IBM front (pyTree) <Examples/Connector/buildFrontIBMPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Connector/buildFrontIBMPT.py

---------------------------------------
 
.. py:function:: Connector.IBM.setInterpDataIBM(t, tc, tb, front, front2=None, dimPb=3, frontType=1, DEPTH=2, IBCType=1, interpDataType=0, Reynolds=1.e6, yplus=100, Lref=1., twoFronts=False, cartesian=False)

    Compute transfer coefficients and data for IBM and store them in tc.

    *Example of use:*
    
    * `Compute IBM coefficients (pyTree) <Examples/Connector/setInterpDataIBMPT.py>`_:
    
    .. literalinclude:: ../build/Examples/Connector/setInterpDataIBMPT.py
