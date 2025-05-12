.. CPlot documentation master file

:tocdepth: 2

CPlot: a light plotter for arrays/pyTree
=========================================

Preamble
########

CPlot is a simple plotter for arrays (as defined in Converter 
documentation) or for CGNS/python trees (pyTrees as defined in Converter/Internal 
documentation).

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use with the array interface, you have to import CPlot module::

   import CPlot

For use with the pyTree interface::

   import CPlot.PyTree as CPlot


.. py:module:: CPlot

List of functions
##################

**-- Actions**

.. autosummary::
   :nosignatures:

   CPlot.display
   CPlot.render
   CPlot.delete
   CPlot.add
   CPlot.replace
   CPlot.finalizeExport
   CPlot.PyTree.display360

**-- Set / Get functions**

.. autosummary::
   :nosignatures:

   CPlot.getState
   CPlot.getSelectedZone
   CPlot.getSelectedZones
   CPlot.getSelectedStatus
   CPlot.getActiveZones
   CPlot.getActiveStatus
   CPlot.getActivePoint
   CPlot.getActivePointIndex
   CPlot.getMouseState
   CPlot.getKeyboard
   CPlot.resetKeyboard
   CPlot.changeVariable
   CPlot.changeStyle
   CPlot.changeBlanking
   CPlot.setState
   CPlot.setMode
   CPlot.setSelectedZones
   CPlot.unselectAllZones

**-- Camera setting and motion**

.. autosummary::
   :nosignatures:

    CPlot.lookFor
    CPlot.moveCamera
    CPlot.travelLeft

**-- Set rendering informations in pyTree**

.. autosummary::
   :nosignatures:

    CPlot.PyTree.addRender2Zone
    CPlot.PyTree.addRender2PyTree
    CPlot.PyTree.loadView


Contents
#########

Keys in CPlot window
---------------------------

    Keys must be pressed when CPlot window is active.

    + **f**: fit view to data.
    + **Ctrl+f**: switch between full screen and windowed mode.
    + **Left/right Arrows** or **mouse drag**: rotate model.
    + **Up/down Arrows** or **mouse drag**: rotate model.
    + **Ctrl + Up/down Arrows** or **mouse wheel**: zoom in and out
    + **Shift + Arrows** or **right mouse drag**: translate model.
    + **Ctrl + Arrows** or **Ctrl + right mouse drag**: tilt model.

    + **Shift + left mouse click**: select zone.
    + **Shift + Ctrl + left mouse click**: multiple select.
    + **Ctrl + left mouse click**: Accurate select (click on nearest mesh node)
    + **Shift + right mouse click**: deactivate (hide) zone.
    + **Shift + double left mouse click**: center view on clicked point.

    + **1** or **Shift+1**: toggle between mesh/solid modes.
    + **2** or **Shift+2**: display fields (switch variable - next and previous).
    + **3**: set render mode.
    + **Space bar**: toggle select all zones.
    + **\< \>**: toggle ambiguous selections.
    + **m** or **M**: toggle between 2D and 3D mode.
    + **z** or **Z**: select zones one by one.
    + **a** or **A**: activate(show)/deactivate(hide) a selected zone.
    + **l**: look for selected zone.
    + **i** or **I** or **Ctrl+i** or **Ctrl+I**: change displayed i plane (structured zones).
    + **j** or **J** or **Ctrl+j** or **Ctrl+J**: change displayed j plane.
    + **k** or **K** or **Ctrl+k** or **Ctrl+K**: change displayed k plane.
    + **q**: quit.
    
---------------------------------------------------------------------------

Actions
--------------------------


.. py:function:: CPlot.display(a, ...)

    .. A2.O0.D1
    
    Display entity.
    Display function has a lot of optional options that can be specified as arguments.
    In offscreen mode, you can render with openGL if you have a GPU or with osmesa
    if you only have a CPU.

    :param a: input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param dim: dimension of data. 1: 1D, 2: 2D, 3: 3D (default: 3)
    :type dim: int
    :param mode: display mode. 0 or 'Mesh': mesh, 1 or 'Solid': solid, 2 or 'Render': render, 3 or 'Scalar': scalar field, 4 or 'Vector': vector field (default: 0)
    :type mode: int or string
    :param scalarField: scalar field number or scalar field name (ex:'Density')
    :type scalarField: int or string
    :param vectorField1,2,3: vector field number or vector field name
    :type vectorField1,2,3: int or string
    :param displayInfo: 0 means no info display (default: 1)
    :type displayInfo: int
    :param displayIsoLegend: 0 means no iso legend display (default: 0)
    :type displayIsoLegend: int
    :param meshStyle: 0: white solid and red wireframe, 1: colored wireframe, 2: colored solid and wireframe, 3: cyan solid and black wireframe, 4: colored solid and black wireframe (default: 2)
    :type meshStyle: int
    :param solidStyle: 0: blue, 1: colored by zone, 3: white, 4: colored by zone outlined (default: 0)
    :type solidStyle: int
    :param scalarStyle: 0: banded, 1: banded+mesh, 2: lines, 3: lines+mesh (default: 0)
    :type scalarStyle: int
    :param vectorStyle: 0: RGB, 1: arrows, 2: lines (default: 0)
    :type vectorStyle: int
    :param vectorDensity: the density of vectors (default: 0.)
    :type vectorDensity: float
    :param vectorScale: scale of vector in % (default: 100.)
    :type vectorScale: float in 0-100
    :param vectorNormalize: if 1, displayed vectors are normalized (default: 0)
    :type vectorNormalize: 0 or 1
    :param vectorShowSurface: if 1, display surface in vector mode (vectorStyle=1) (default: 1)
    :type vectorShowSurface: 0 or 1
    :param vectorShape: type of arrows for vectors (vectorStyle=1) (default: 0)
    :type vectorShape: 0 (3D arrows), 1 (Flat arrows), 2 (Tetra arrows)
    :param vectorProjection: if 1, vectors are projected on surface (default: 0)
    :type vectorProjection: 0 or 1
    :param colormap: 0-1: Blue2Red, 2-3: BiColorRGB, 4-5: BiColorHSV, 6-7: TriColorRGB, 8-9: TriColorHSV, 10-11: MultiColorRGB, 12-13: MultiColorHSV, 14-15: Diverging, 16-17: Viridis, 18-19: Inferno, 20-21: Magma, 22-23: Plasma, 24-25: Jet, 26-27: Greys, 28-29: Nice Blue, 30-31: Greens (default: 0)
    :type colormap: int (upper number activates light)
    :param colormapC1: Hexa string for starting color of bi/tri colors colormaps (ex: #FFFFFF)
    :type colormapC1: string
    :param colormapC2: Hexa string for ending color of bi/tri colors colormaps (ex: #FFFFFF)
    :type colormapC2: string
    :param colormapC3: Hexa string for mid color of tri colors colormaps (ex: #FFFFFF)
    :type colormapC3: string
    :param niso: number of isos (default: 25)
    :type niso: int
    :param isoEdges: width of iso edges for scalar display (default: -1)
    :type isoEdges: float
    :param isoScales: list of min and max of a variable ([varName, niso, min, max] or [varName, niso, min, max, cutmin, cutmax])(default: [])
    :type isoScales: list of [string, int, float, float] or [string, int, float, float, float, float]. Additional colormap number can be added.
    :param win: (sizeWinX, sizeWinY) window size (default: 700,700)
    :type win: tuple of 2 ints
    :param posCam: (x,y,z) camera position
    :type posCam: tuple of 3 floats
    :param posEye: (x,y,z) eye position
    :type posEye: tuple of 3 floats
    :param dirCam: (x,y,z) camera direction (default: 0,0,1)
    :type dirCam: tuple of 3 floats
    :param viewAngle: camera angle in degrees (default: 50.)
    :type viewAngle: float
    :param bgColor: background color. 0-13 (default: 0)
    :type bgColor: int
    :param backgroundFile: name of a image png file (default: 0)
    :type backgroundFile: string
    :param shadow: display shadows. 0-1 (default: 0)
    :type shadow: int
    :param dof: depth of field smoothing. 0-1 (default: 0)
    :type dof: int
    :param dofPower: power of depth of field effect (default: 3.)
    :type dofPower: float 
    :param lighOffset: offset to default light position (default: (0,0))
    :type lightOffset: tuple of two floats
    :param gamma: gamma correction (default: 1.)
    :type gamma: float 
    :param toneMapping: none (0), ACE (1), Filmic (2), Uchimura (3)
    :type toneMapping: int     
    :param stereo: 1 or 2 means red/cyan anaglyph (default: 0)
    :type stereo: int
    :param stereoDist: distance between eyes for stereo
    :type stereoDist: float
    :param export: file name for export (.png)
    :type export: string
    :param exportResolution: resolution for export ("1920x1080")
    :type exportResolution: string
    :param zoneNames: optional list of zone names (same size as arrays, struct zones, then unstruct zones)
    :type zoneNames: list of strings
    :param renderTags: optional list of render tags (same size as arrays, struct zones, then unstruct zones)
    :type renderTags: list of strings
    :param frameBuffer: the number of the frame buffer we are rendering to for offscreen rendering
    :type frameBuffer: int
    :param offscreen: 0: direct rendering, 1: offscreen rendering (osmesa), 2: offscreen rendering (openGL), 3: partial composite offscreen rendering (openGL), 4: final composite offscreen rendering (openGL), 5: partial composite offscreen rendering (osmesa), 6: final composite offscreen rendering (osmesa), 7: parallel offscreen rendering (osmesa) (default: 0)
    :type offscreen: int

    *Example of use:*

    * `Display (array) <Examples/CPlot/display.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/display.py

    * `Display (pyTree) <Examples/CPlot/displayPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/displayPT.py
    
    * `Display offscreen with openGL (pyTree) <Examples/CPlot/displayOffScreen2PT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/displayOffScreen2PT.py
    
    * `Display offscreen with osmesa (pyTree) <Examples/CPlot/displayOffScreenPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/displayOffScreenPT.py

------------------------------------------

.. py:function:: CPlot.PyTree.display360(a, stereo=0, offscreen=1, ...)
    
    Display a 360-degree view of the pyTree. This function is useful for creating
    panoramic views of the data.
    Accepts the same argument as CPlot.display.
    If stereo=0, no stereo.
    If stereo=1, ODS stereo (only for offscreen=1 or 7).

    :param a: input pyTree
    :type a: pyTree
    :param stereo: stereo type.
    :type stereo: 0, 1, 2
    :param offscreen: renderer: 1 (osmesa), 2 (openGL)
    :type offscreen: 1, 2

    *Example of use:*

    * `Display 360 view (pyTree) <Examples/CPlot/display360PT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/display360PT.py

    * `Display 360 view stereo (pyTree) <Examples/CPlot/display360sPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/display360sPT.py

------------------------------------------

.. py:function:: CPlot.render()

    .. A2.O0.D1
    
    Force rendering. Must be used after functions that don't render (ex: add, delete, ...).

    *Example of use:*

    * `Render (array) <Examples/CPlot/render.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/render.py

    * `Render (pyTree) <Examples/CPlot/renderPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/renderPT.py

------------------------------------------

    
.. py:function:: CPlot.delete(list)

    .. A2.O0.D1
    
    Delete zones from plotter. This function does not render. Argument is either
    a list of zone numbers (struct zones then unstructured zones order) or a list 
    of zone names if zoneNames argument has been provided before to display function.

    :param list: list of zone number or zone names
    :type list: list of int or strings

    *Example of use:*

    * `Delete zones (array) <Examples/CPlot/delete.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/delete.py

    * `Delete zones (pyTree) <Examples/CPlot/deletePT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/deletePT.py    
     

-------------------------------------------

.. py:function:: CPlot.add(A, ..., a)

    .. A2.O0.D1
    
    Add/insert one array/zone in plotter. This function does not render. 
    
    For array interface:
    ::
    
        CPlot.add(A, no, a)
    
    no is the position of insertion of a in A.
    Replace also in A.    

    For the pyTree interface:
    ::
    
        CPlot.add(A, nob, noz, a)
    
    Insert or append a to the base nob, at position
    noz in the zone list. If noz=-1, append to end of list.

    :param A: input data
    :type A: arrays, pyTree or list of zones
    :param no: position of zone to add in A
    :type no: int
    :param nob: position of base of zone to add in A
    :type nob: int
    :param noz: position of zone to add in A
    :type noz: int
    :param a: data to add
    :type a: array or zone
    
    *Example of use:*

    * `Add a zone (array) <Examples/CPlot/add.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/add.py

    * `Add zone (pyTree) <Examples/CPlot/addPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/addPT.py    


---------------------------------------------

.. py:function:: CPlot.replace(A, ..., a)

    .. A2.O0.D1

    For array interface::
    
        CPlot.replace(A, no, a)
        
    Performs A[no]=a, keeping plotter coherent.
    
    For pyTree interface::
    
        CPlot.replace(A, nob, noz, a)
        
    Performs t[2][nob][2][noz]=a, keeping plotter coherent. 
    This function does not render. 

    :param A: input data
    :type A: arrays, pyTree or zones
    :param no: position of zone to add in A
    :type no: int
    :param nob: position of base of zone to add in A
    :type nob: int
    :param noz: position of zone to add in A
    :type noz: int
    :param a: data to add
    :type a: array or zone
    
    *Example of use:*

    * `Replace a zone (array) <Examples/CPlot/replace.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/replace.py

    * `Replace a zone (pyTree) <Examples/CPlot/replacePT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/replacePT.py


-----------------------------------------------
    

.. py:function:: CPlot.finalizeExport(action=0)

    .. A2.O0.D0
    
    Finalize an export. Wait for the end of file writing.
    
    Argument must be identical to the offscreen argument of display function.
    
    When doing offscreen rendering with openGL (offscreen=2,3,4), finalize export is needed for
    offscreen=2 (wait until file is written) and offscreen=4 (wait until file is written and clear buffers for next image).
    When doing offscreen rendering with osmesa (offscreen=1,5,6,7), finalize export is needed
    for offscreen=6 (clear buffer for next image).

    if action=-1, close the mpeg file.

    :param action: identical to offscreen argument of display
    :type action: int

    *Example of use:*

    * `Finalize an export (array) <Examples/CPlot/displayOffScreen2.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/displayOffScreen2.py


---------------------------------------------------------------------------

Set / Get functions
--------------------------

.. py:function:: CPlot.getState(stateName)

    .. A2.O0.D0
    
    Return the specified state value as stored in plotter.
    Available stateName are the same as the display
    function arguments.

    :param stateName: name of state to be retrieved
    :type stateName: string

    *Example of use:*

    * `Get state from plotter (array) <Examples/CPlot/getState.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getState.py

    * `Get state from plotter (pyTree) <Examples/CPlot/getStatePT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getStatePT.py

---------------------------------------

.. py:function:: CPlot.getSelectedZone()

    .. A2.O0.D0
    
    Return the currently selected zone. If none is selected, return -1. If
    multiple zones are selected, return the last selected zone.

    :return: no of selected zone
    :rtype: int    


    *Example of use:*

    * `Get selected zone (array) <Examples/CPlot/getSelectedZone.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getSelectedZone.py    

---------------------------------------

.. py:function:: CPlot.getSelectedZones()

    .. A2.O0.D0
    
    Return the list of selected zones. If none is selected, return [].

    *Example of use:*

    * `Get selected zones (array) <Examples/CPlot/getSelectedZones.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getSelectedZones.py    

---------------------------------------


.. py:function:: CPlot.getSelectedStatus(nz)

    .. A2.O0.D0
    
    Return the selected status (1: selected, 0: not selected) of zone nz.

    :param nz: zone number
    :type nz: int
    :return: status of zone (1: selected, 0: not selected)
    :rtype: int

    *Example of use:*

    * `Get selected status of zone (array) <Examples/CPlot/getSelectedStatus.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getSelectedStatus.py    


---------------------------------------


.. py:function:: CPlot.getActiveZones()

    .. A2.O0.D0
    
    Return the list of active (visible) zones.

    :return: list of zone numbers
    :rtype: list of ints

    *Example of use:*

    * `Get active zones (array) <Examples/CPlot/getActiveZones.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getActiveZones.py    

-----------------------------------------------

.. py:function:: CPlot.getActiveStatus(nz)

    .. A2.O0.D0
    
    Return the active status (1: active, 0: inactive) of zone nz.
    Active zones are visible, unactive zones are hidden.
    
    :param nz: number of zone
    :type nz: int
    :return: active status of zone
    :rtype: int

    *Example of use:*

    * `Get active status of zone (array) <Examples/CPlot/getActiveStatus.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getActiveStatus.py    


-----------------------------------------------

.. py:function:: CPlot.getActivePoint()

    .. A2.O0.D0
    
    Return the last clicked point position (coordinates in 3D world) as 
    a list of three coordinates.

    :return: active point position
    :rtype: tuple of 3 floats

    *Example of use:*

    * `Get active point (array) <Examples/CPlot/getActivePoint.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getActivePoint.py    

-----------------------------------------------

.. py:function:: CPlot.getActivePointIndex()
    
    .. A2.O0.D0
    
    Return the active point index. 
    For structured grids, return [ind, indc, i,j,k], where ind is the global index of the nearest node
    to active point, indc is the global index of the nearest center to
    active point and i,j,k are the indices of nearest node. 
    For unstructured grids, return [ind, indc, no, 0, 0],
    where ind is the global index of nearest node, indc 
    is the nearest center to the clicked point and no is the number of ind
    in the center2node connectivity of indc.
    If there is no active point, return [].

    :return: active point index
    :rtype: list of 4 ints

    *Example of use:*

    * `Get active point index (array) <Examples/CPlot/getActivePointIndex.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getActivePointIndex.py

-----------------------------------------------


.. py:function:: CPlot.getMouseState()
    
    .. A2.O0.D0
    
    Return the current button state
    of mouse (0: left pressed, 1: middle pressed, 2: right pressed, 
    5: not pressed) and the current mouse position (if pressed). Use 
    it when dragging.
    
    :return: mouse state and mouse position 
    :rtype: tuple of 2 ints

    *Example of use:*

    * `Get mouse state (array) <Examples/CPlot/getMouseState.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getMouseState.py


-----------------------------------------------

.. py:function:: CPlot.getKeyboard()

    .. A2.O0.D0
    
    Return the pressed keys as a string.

    :return: keys pressed in CPlot window
    :rtype: string

    *Example of use:*

    * `Get keyboard keys (array) <Examples/CPlot/getKeyboard.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getKeyboard.py


-----------------------------------------------
    
.. py:function:: CPlot.resetKeyboard()

    .. A2.O0.D0 
    
    Reset the pressed key string stored in plotter.
    
-----------------------------------------------

.. py:function:: CPlot.changeVariable()

    .. A2.O0.D0
    
    Change displayed variable.
    
    *Example of use:*

    * `Change variable (array) <Examples/CPlot/changeVariable.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/changeVariable.py


-----------------------------------------------

.. py:function:: CPlot.changeStyle()

    .. A2.O0.D0
    
    Change CPlot display style.

    *Example of use:*

    * `Change Style (array) <Examples/CPlot/changeStyle.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/changeStyle.py
    
-----------------------------------------------

.. py:function:: CPlot.changeBlanking()

    Change the blanking procedure.
    
    *Example of use:*

    * `Change the blanking procedure (array) <Examples/CPlot/changeBlanking.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/changeBlanking.py

-----------------------------------------------

.. py:function:: CPlot.setState(dim, mode, ...)

    .. A2.O0.D0
    
    Set a CPlot state. The same keywords as display can be used.

    Additional keywords are:
    
    + **ghostifyDeactivatedZones**: 1 means deactivated zones will appear blended.
    + **edgifyActivatedZones**: 1 means activated zones will appear as edges.
    + **edgifyDeactivatedZones**: 1 means deactivated zones will appear as edges.
    + **message**: "A string" or "Clear"
    + **viewAngle**: the camera angle (default: 50 degrees).
    + **cursor**: mouse cursor type (0: normal, 1: cross, 2: wait).
    + **sobelThreshold**: sobel threshold for zone outlines (default: -0.5).
    + **sharpenPower**: power of sharpening image post-processing (default: -0.5).
    + **selectionStyle**: style for selection (default: 0).
    + **activateShortCut**: if False, deactivate shortCut keys (def: True).
    + **billBoards**: list of billboard image files ['file.png',1,1] (default: None).
    + **billBoardSize**: size of billboard. If -1, use distance to fit billboards (default: -1).
    + **materials**: list of material image files used in textured rendering ['mat.png'] (default: None).

    *Example of use:*
    
    * `Set a state in CPlot (array) <Examples/CPlot/setState.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/setState.py    


-----------------------------------------------


.. py:function:: CPlot.setMode(mode)

    .. A2.O0.D0
    
    Set CPlot display mode (0 or 'Mesh': means mesh, 1 or 'Solid': means solid, 2 or 'Render': means render, 
    3 or 'Scalar' means scalar field, 4 or 'Vector' means vector fields)

    :param mode: mode to set
    :type mode: int or string

    *Example of use:*

    * `Set a mode in CPlot (array) <Examples/CPlot/setmode.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/setMode.py    
    

-----------------------------------------------

.. py:function:: CPlot.setSelectedZones(list)

    Set the selected zone status (1: selected, 0: not selected) by zone global
    number.

    :param list: list of zone number and status
    :type  list: list of tuples of 2 ints

    *Example of use:*

    * `Set the selected status for a list of zones (array) <Examples/CPlot/setSelectedZones.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/setSelectedZones.py    


-----------------------------------------------

.. py:function:: CPlot.unselectAllZones()

    .. A2.O0.D0
    
    Unselect all zones.

    *Example of use:*

    * `Unselect all zones (array) <Examples/CPlot/unselectAllZones.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/unselectAllZones.py    


-----------------------------------------------

.. py:function:: CPlot.setActiveZones(list)

    Set the active (visible) status for given zones.

    :param list: list of zone number and status
    :type  list: list of tuples of 2 ints

    *Example of use:*

    * `Set the active status for a list of zones (array) <Examples/CPlot/setActiveZones.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/setActiveZones.py


-----------------------------------------------

.. py:function:: CPlot.setZoneNames(list)

    Set the specified zone names.

    :param list: list of zone number and zone names
    :type  list: list of tuples of int and string

    *Example of use:*

    * `Set the zone names for a list of zones (array) <Examples/CPlot/setZoneNames.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/setZoneNames.py

Camera setting and motion
-----------------------------------------------

.. py:function:: CPlot.lookFor()

    Look for selected zone. It positions the camera for a clear view
    on the currently selected zone.

    *Example of use:*

    * `Look for selected zone (array) <Examples/CPlot/lookFor.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/lookFor.py

-----------------------------------------------

.. py:function:: CPlot.moveCamera(posCams, posEyes=None, dirCams=None, moveEye=False, N=100, speed=1., pos=-1)

    Move camera along camera points. The list of points specifies the path of the camera.
    The camera will move along this path, making N steps. pos will position the 
    camera at step number pos along the path. 
    If posEyes is specified, posEye (that
    is the position the camera is looking to) will follow this path.
    If posEyes is not specified and moveEye is true, the posEye will follow the path
    Otherwise, the posEye will stay at initial posEye.

    :param posCams: coordinates of camera points
    :type posCams: list of tuple of 3 floats or 1D STRUCT Zone
    :param posEyes: coordinates of eyes points
    :type posEyes: list of tuple of 3 floats or 1D STRUCT Zone
    :param dirCams: camera directions
    :type dirCams: list of tuple of 3 floats or 1D STRUCT Zone
    :param moveEye: if True, the eye follow camera points
    :type moveEye: Boolean
    :param speed: speed of camera motion
    :type speed: float
    :param N: number of camera positions
    :type N: int
    :param pos: position in posCams (in 0,N)
    :type pos: int
    :return: current posCam, posEye, dirCam
    :rtype: 3 lists of 3 floats

    *Example of use:*

    * `Move camera along check points (array) <Examples/CPlot/moveCamera.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/moveCamera.py

    * `Move camera along check points (pyTree) <Examples/CPlot/moveCameraPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/moveCameraPT.py

    * `Move camera along check points in offscreen mode (pyTree) <Examples/CPlot/moveCameraOffScreenPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/moveCameraOffScreenPT.py

-----------------------------------------------


.. py:function:: CPlot.travelLeft(xr, N=100)

    .. A2.O0.D0
    
    Travel camera left/Right/Up/Down/In/Out. 
    Xr is the range (in 0.,1.). 
    N is the number of check points.

    :return: final posCam, posEye, dirCam
    :rtype: 3 lists of 3 floats

    *Example of use:*

    * `Travel camera (array) <Examples/CPlot/travel.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/travel.py

    * `Travel camera (pyTree) <Examples/CPlot/travel.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/travel.py

---------------------------------------------------------------------------

Set rendering informations in pyTree
--------------------------------------

.. py:function:: CPlot.PyTree.addRender2Zone(a, material=None, color=None, blending=None, meshOverlay=None, shaderParameters=None)

    Add rendering info to a zone. Info are added in a .RenderInfo user
    defined node. Use Render mode in display for rendering effects.
    Exists also as in place version (_addRender2Zone) that modifies a
    and returns None.
    Shader parameters are described in shaderSettings_.

    :param a: input zone
    :type a: zone node
    :param material: material to set (in 'Solid', 'Flat', 'Glass', 'Chrome', 'Metal', 'Wood', 'Marble', 'Granite', 'Brick', 'XRay', 'Cloud', 'Gooch', 'Sphere', 'Texmat')
    :type material: string
    :param color: color to set (in 'White', 'Grey', ... or '#FFFF')
    :type color: string
    :param blending: opacity factor (in [0.,1.])
    :type blending: float
    :param meshOverlay: if 1 then overlay the mesh
    :type meshOverlay: 0 or 1
    :param shaderParameters: two float that parametrize shaders
    :type shaderParameters: list of two floats in [0.,2.]

    *Example of use:*

    * `Add render information to zone (pyTree) <Examples/CPlot/addRender2ZonePT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/addRender2ZonePT.py
    
-----------------------------------------------

.. py:function:: CPlot.PyTree.addRender2PyTree(a, slot=0, posCam=None, posEye=None, dirCam=None, mode=None, scalarField=None, niso=None, isoScales=None, isoEdges=None, isoLight=None, isoLegend=None, colormap=None, materials=None, bumpMaps=None, billBoards=None)

    Add rendering info to a tree. Info are added in a .RenderInfo user
    defined node. To load the settings to the view, call explicitely CPlot.loadView.
    Exists also as in place version (_addRender2PyTree) that modifies a and
    returns None.

    :param a: input tree
    :type a: pyTree
    :param slot: the slot number
    :type slot: int
    :param posCam: camera position
    :type posCam: list of 3 floats
    :param posEye: eye position
    :type posEye: list of 3 floats
    :param dirCam: camera direction
    :type dirCam: list of 3 floats
    :param mode: displayed mode ('Mesh', 'Solid', 'Scalar', 'Vector', 'Render') 
    :type mode: string
    :param scalarField: scalar field to display in Scalar mode
    :type scalarField: string
    :param niso: number of isos to display in Scalar mode
    :type niso: int
    :param isoScales: list of min and max of a variable ([varName, niso, min, max] or [varName, niso, min, max, cutmin, cutmax])(default: [])
    :type isoScales: list of [string, int, float, float] or [string, int, float, float, float, float]
    :param isoEdges: size of edges in Scalar mode
    :type isoEdges: float
    :param isoLight: set to 1 if light is used in Scalar mode
    :type isoLight: 0 or 1
    :param isoLegend: set to 1 if legend is displayed in Scalar mode
    :type isoLegend: 0 or 1
    :param colormap: name of the colormap for Scalar mode ('Blue2Red', ...)
    :type colormap: string
    :type materials: list of image file names used for texture mapping
    :param materials: list of strings
    :type bumpMaps: list of image file names used for bump mapping
    :param bumpMaps: list of strings
    :type billBoards: list of image file names used for billboarding
    :param billboards: list of strings

    *Example of use:*

    * `Add render information to pyTree (pyTree) <Examples/CPlot/addRender2PyTreePT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/addRender2PyTreePT.py
    
-----------------------------------------------

.. py:function:: CPlot.PyTree.loadView(a, slot=0)

    Load a view defined in a slot to the plotter.
    A view must already have been stored in pyTree a using CPlot.addRender2PyTree.
    
    :param a: input tree
    :type a: pyTree
    :param slot: number of slot to load
    :type slot: int

    *Example of use:*

    * `Load (pyTree) <Examples/CPlot/loadViewPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/loadViewPT.py


---------------------------------------------------------------------------

Shader settings
---------------------------
.. _shaderSettings:

    Shaders can be adjusted with the two parameters of the addRender2Zone function.

    Here are the meaning of each parameters for each shader:

    - Solid: [1] specularFactor; [2] diffuseFactor
    - Flat: [1] None; [2] None
    - Glass: [1] mix; [2] mix2
    - Chrome: [1] mix color/envmap; [2] base color intensity
    - Metal: [1] reflection intensity; [2] bump size
    - Wood: [1] ray size, [2] bump height
    - Marble: [1] ray size, [2] None
    - Granite: [1] bump size; [2] bump height
    - Brick: [1] brick size; [2] color width
    - XRay: [1] fall off; [2] color modifier
    - Cloud: [1] size; [2] None
    - Gooch : [1] width; [2] shininess
    - Sphere: [1] size of sphere; [2] type of billboard
    - TexMat: [1] specularFactor; [2] texture number



---------------------------------------------------------------------------

.. toctree::
   :maxdepth: 2


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

