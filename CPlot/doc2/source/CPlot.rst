.. CPlot documentation master file

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

   CPlot.display
   CPlot.render
   CPlot.delete
   CPlot.add
   CPlot.replace
   CPlot.pressKey
   CPlot.finalizeExport


**-- Set / Get functions**

.. autosummary::

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

    CPlot.lookFor
    CPlot.moveCamera
    CPlot.travelLeft

**-- Set rendering informations in pyTree**

.. autosummary::

    CPlot.PyTree.addRender2Zone
    CPlot.PyTree.addRender2PyTree
    CPlot.PyTree.loadView


Contents
#########

Keys in CPlot window
---------------------------

    Must be pressed when CPlot window is active.

    + **f**: fit view to data.
    + **Ctrl+f**: switch between full screen and windowed mode.
    + **Arrows** or **left mouse**: move around.
    + **Shift + Arrows** or **right mouse**: strafe.
    + **Ctrl + Arrows** or **Ctrl + left mouse**: move your head.
    + **Ctrl + right mouse**: tilt your head.
    + **Shift + left mouse**: select zone.
    + **Shift + Ctrl + left mouse**: multiple select.
    + **Ctrl + left mouse**: Accurate select (click on node)
    + **Shift + right mouse**: deactivate zone.
    + **Shift + double left mouse**: center view.
    + **o** or **left mouse**: move up.
    + **p** or **left mouse**: move down.
    + **1** or **Shift+1**: display fields (switch variable next and previous).
    + **Space bar**: display mesh.
    + **Shift+Space bar**: display solid.
    + **m** or **M**: toggle between 2D and 3D mode.
    + **z** or **Z**: select zones one by one.
    + **a** or **A**: activate(show)/deactivate(hide) a selected zone.
    + **l**: look for selected zone.
    + **i** or **I** or **Ctrl+i** or **Ctrl+I**: change displayed i plane (structured zones).
    + **j** or **J** or **Ctrl+j** or **Ctrl+J**: change displayed j plane.
    + **k** or **K** or **Ctrl+k** or **Ctrl+K**: change displayed k plane.
    + **q**: quit.
    
Actions
--------------------------


.. py:function:: CPlot.display(a, ...)

    Display entity.
    Display function has a lot of optional options that can be specified as arguments.

    :param a: Input data
    :type a: [array, list of arrays] or [pyTree, base, zone, list of zones]
    :param dim: dimension of data. 1: 1D, 2: 2D, 3: 3D (default: 3)
    :type dim: int
    :param mode: display mode. 0 or 'Mesh': mesh, 1 or 'Solid': solid, 2 or 'Render': render, 3 or 'Scalar': scalar field, 4 or 'Vector': vector field (default: 0)
    :type mode: int or string
    :param scalarField: scalar field number or scalar field name
    :type scalarField: int or string
    :param vectorField1,2,3: vector field number or vector Field name
    :type vectorField1,2,3: int or string
    :param displayBB: 0: bounding box display (default: 1)
    :type displayBB: int
    :param displayInfo: 0 means no info display (default: 1)
    :type displayInfo: int
    :param displayIsoLegend: 0 means no iso legend display (default: 0)
    :type displayIsoLegend: int
    :param meshStyle: 0: white solid and red wireframe, 1: colored wireframe, 2: colored solid and wireframe, 3: cyan solid and black wireframe (default: 2)
    :type meshStyle: int
    :param solidStyle: 0: blue, 1: colored by zone, 3: cyan (default: 1)
    :type solidStyle: int
    :param scalarStyle: 0: banded, 1: banded+mesh, 2: lines, 3: lines+mesh (default: 0)
    :type scalarStyle: int
    :param vectorStyle: 0: RGB, 1: triangles, 2: lines, 3: arrows, 4: uniform lines
    :type vectorStyle: int
    :param vectorDensity: the density of vectors
    :type vectorDensity: float
    :param vectorNormalize: if 1, displayed vectors are normalized
    :type vectorNormalize: 0 or 1
    :param colormap: 0: Blue2Red, 2: Green2Red, 4: Black2White, 6: White2Black, 8: Diverging (default: 0)
    :type colormap: int
    :param niso: number of isos (default: 25)
    :type niso: int
    :param isoEdges: width of iso edges for scalar display (default: -1)
    :type isoEdges: float
    :param isoScales: list of min and max of a variable [novar, niso, min, max] (default: [])
    :type isoScales: list of set of 4 values
    :param win: (sizeWinX, sizeWinY) window size (default: 700,700)
    :type win: tuple of 2 ints
    :param posCam: (x,y,z) camera position
    :type posCam: tuple of 3 floats
    :param posEye: (x,y,z) eye position
    :type posEye: tuple of 3 floats
    :param dirCam: (x,y,z) camera direction (default: 0,0,1)
    :type dirCam: tuple of 3 floats
    :param viewAngle: camera angle in degrees (default: 50)
    :type viewAngle: float
    :param bgColor: background color. 0-10 (default: 0)
    :type bgColor: int
    :param shadow: display shadows. 0-1 (default: 0)
    :type shadow: int
    :param dof: depth of field smoothing. 0-1 (default: 0)
    :type dof: int
    :param stereo: 1 or 2 means red/cyan anaglyph (default: 0)
    :type stereo: int
    :param stereoDist: distance between eyes for stereo
    :type stereoDist: float
    :param export: file name for export
    :type export: string
    :param exportResolution: resolution for export ("1200x900")
    :type exportResolution: string
    :param zoneNames: optional list of zone names (same size as arrays, struct zones, then unstruct zones)
    :type zoneNames: list of strings
    :param renderTags: optional list of render tags (same size as arrays, struct zones, then unstruct zones)
    :type renderTags: list of strings
    :param offscreen: 1 means offscreen rendering (mesa), 2 means offscreen rendering (openGL) (default: 0)
    :type offscreen: int

    *Example of use:*

    * `Display (array) <Examples/CPlot/display.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/display.py

    * `Display (pyTree) <Examples/CPlot/displayPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/displayPT.py

------------------------------------------

.. py:function:: CPlot.render()

    Force rendering.

    *Example of use:*

    * `Render (array) <Examples/CPlot/render.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/render.py

    * `Render (pyTree) <Examples/CPlot/renderPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/renderPT.py

------------------------------------------

    
.. py:function:: CPlot.delete(list)

    Delete zones from plotter. This function does not render. Argument is either
    a list of zone numbers (struct zones then unstructured zones order) or a list 
    of zone names if zoneNames arg has been provided before to display function:

    :param list: list of zone number or zone names
    :type list: list of int or strings

    *Example of use:*

    * `Delete zones (array) <Examples/CPlot/delete.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/delete.py

    * `Delete zones (pyTree) <Examples/CPlot/deletePT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/deletePT.py    
     

-------------------------------------------

.. py:function:: CPlot.add(A, no, a)

    Add/insert one array/zone in plotter. This function does not render. 
    For array interface, no is the position of insertion of a in A.
    Replace also in A

    For the pyTree interface, insert or append a to the base nob, at position
    noz in the zone list. If noz=-1, append to end of list.

    :param A: Initial data
    :type A: arrays, pyTree or zones
    :param no: position of zone to add in A
    :type no: int
    :param a: data to add
    :type a: array or zone
    
    *Example of use:*

    * `Add a zone (array) <Examples/CPlot/add.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/add.py

    * `Add zone (pyTree) <Examples/CPlot/addPT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/addPT.py    


---------------------------------------------

.. py:function:: CPlot.replace(A, no, a)

    Performs A[no]=a, keeping plotter coherent.
    Performs t[2][nob][2][noz]=a, keeping plotter coherent. 
    This function does not render. 

    :param A: Initial data
    :type A: arrays, pyTree or zones
    :param no: position of zone to add in A
    :type no: int
    :param a: data to add
    :type a: array or zone
    
    *Example of use:*

    * `Replace a zone (array) <Examples/CPlot/replace.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/replace.py

    * `Replace a zone (pyTree) <Examples/CPlot/replacePT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/replacePT.py


-----------------------------------------------

.. py:function:: CPlot.pressKey()

    Wait for a key to be pressed.

-----------------------------------------------
    

.. py:function:: CPlot.finalizeExport(action=0)

    Finalize an export. Wait for the end of file writing (action=0)
    or in mpeg export, close the mpeg file (action=1).
    
    :param action: if 0, means wait until file is written, if 1, means close mpeg file
    :type action: int

    *Example of use:*

    * `Finalize an export (array) <Examples/CPlot/displayOffScreen2.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/displayOffScreen2.py


Set / Get functions
--------------------------

.. py:function:: CPlot.getState(stateName)

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

    Return the currently selected zone. If none is selected, return -1. If
    multiple zones are selected, return the last selected zone.

    :return: no of selected zone
    :rtype: int    


    *Example of use:*

    * `Get selected zone (array) <Examples/CPlot/getSelectedZone.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getSelectedZone.py    

---------------------------------------

.. py:function:: CPlot.getSelectedZones()

    Return the list of selected zones. If none is selected, return [].

    *Example of use:*

    * `Get selected zones (array) <Examples/CPlot/getSelectedZones.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getSelectedZones.py    

---------------------------------------


.. py:function:: CPlot.getSelectedStatus(nz)

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

    Return the list of active (displayed) zones.

    :return: list of zone numbers
    :rtype: list of ints

    *Example of use:*

    * `Get active zones (array) <Examples/CPlot/getActiveZones.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getActiveZones.py    

-----------------------------------------------

.. py:function:: CPlot.getActiveStatus(nz)

    Return the active status (1: active, 0: inactive) of zone nz.
    
    :param nz: number of zone
    :type nz: int
    :return: active status of zone
    :rtype: int

    *Example of use:*

    * `Get active status of zone (array) <Examples/CPlot/getActiveStatus.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getActiveStatus.py    


-----------------------------------------------

.. py:function:: CPlot.getActivePoint()

    Return the last clicked point position (coordinates in 3D world) as 
    a list of three coordinates.

    :return: active point position
    :rtype: tuple of 3 floats

    *Example of use:*

    * `Get active point (array) <Examples/CPlot/getActivePoint.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getActivePoint.py    

-----------------------------------------------

.. py:function:: CPlot.getActivePointIndex()
    
    Return the active point index. For structured grids, return 
    [ind, indc, i,j,k], where ind is the global index of the nearest node
    to active point, indc is the global index of the nearest center to
    active point and i,j,k are the indices of nearest node. For
    unstructured grids, return [ind, indc, no, 0, 0],
    where ind is the global index of nearest node, indc 
    is the nearest center to the clicked point and no is the number of ind
    in the center2node connectivity of indc.

    :return: active point index
    :rtype: list of 4 ints

    *Example of use:*

    * `Get active point index (array) <Examples/CPlot/getActivePointIndex.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getActivePointIndex.py

-----------------------------------------------


.. py:function:: CPlot.getMouseState()
    
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

    Return the pressed keys as a string.

    :return: keys pressed in CPlot window
    :rtype: string

    *Example of use:*

    * `Get keyboard keys (array) <Examples/CPlot/getKeyboard.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/getKeyboard.py


-----------------------------------------------
    
.. py:function:: CPlot.resetKeyboard()

    Reset the pressed key string stored in plotter.
    
-----------------------------------------------

.. py:function:: CPlot.changeVariable()

    Change displayed variable.
    
    *Example of use:*

    * `Change variable (array) <Examples/CPlot/changeVariable.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/changeVariable.py


-----------------------------------------------

.. py:function:: CPlot.changeStyle()

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

    Set a CPlot state. The same keywords as display can be used.

    Additional keywords are:
    
    + **ghostifyDeactivatedZones**: 1 means deactivated zones will appear blended.
    + **edgifyActivatedZones**: 1 means activated zones will appear as edges.
    + **edgifyDeactivatedZones**: 1 means deactivated zones will appear as edges.
    + **message**: "A string" or "Clear"
    + **viewAngle**: the camera angle (default: 50 degrees).
    + **cursor**: mouse cursor type (0: normal, 1: cross, 2: wait).
    + **lightOffset**: offset to default light position (default: (0,0)).
    + **dofPower**: power of depth of field effect (default: 6.).
    + **selectionStyle**: style for selection (default: 0).

    *Example of use:*
    
    * `Set a state in CPlot (array) <Examples/CPlot/setState.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/setState.py    


-----------------------------------------------


.. py:function:: CPlot.setMode(mode)

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

    Unselect all zones.

    *Example of use:*

    * `Unselect all zones (array) <Examples/CPlot/unselectAllZones.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/unselectAllZones.py    


-----------------------------------------------

.. py:function:: CPlot.setActiveZones(list)


    Set the active (displayed) zones.

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

    Look for selected zone. It positions the camera for a clear wiew
    on the currently selected zone.

    *Example of use:*

    * `Look for selected zone (array) <Examples/CPlot/lookFor.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/lookFor.py

-----------------------------------------------

.. py:function:: CPlot.moveCamera(list, moveEye=False, N=100, speed=50.)

    Move camera along check points.

    :param listOfPts: coordinates of check points
    :type listOfPts: list of tuple of 3 floats
    :param moveEye: if True, the eye follow check points
    :type moveEye: Boolean
    :param speed: speed of camera motion
    :type speed: float

    *Example of use:*

    * `Move camera along check points (array) <Examples/CPlot/moveCamera.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/moveCamera.py

-----------------------------------------------


.. py:function:: CPlot.travelLeft(xr, N=100)

    Travel camera left/Right/Up/Down/In/Out. 
    Xr is the range (in 0.,1.). 
    N is the number of check points.

    *Example of use:*

    * `Travel camera (array) <Examples/CPlot/travel.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/travel.py



Set rendering informations in pyTree
--------------------------------------

.. py:function:: CPlot.PyTree.addRender2Zone(a, material, color, blending, meshOverlay, shaderParameters)

    Add rendering info to a zone. Info are added in a .RenderInfo user
    defined node. Use Render mode.
    Exists also as in place version (_addRender2Zone) that modifies a
    and returns None.

    :param a: input zone
    :type a: zone node
    :param material: material to set (in 'Solid', 'Flat', 'Glass', 'Chrome', 'Metal', 'Wood', 'Marble', 'Granite', 'Brick', 'XRay', 'Cloud', 'Gooch', 'Sphere')
    :type material: string
    :param color: color to set (in 'White', 'Grey', ... or '#FFFF')
    :type color: string
    :param blending: opacity factor (in [0.,1.])
    :type blending: float
    :param meshOverlay: if 1 then overlay the mesh
    :type meshOverlay: 0 or 1
    :param shaderParameters: two float that parametrize shaders
    :type shaderParameters: list of two floats

    *Example of use:*

    * `Add render information to zone (pyTree) <Examples/CPlot/addRender2ZonePT.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/addRender2ZonePT.py
    
-----------------------------------------------

.. py:function:: CPlot.PyTree.addRender2PyTree(a, slot, posCam, posEye, dirCam, mode, scalarField, niso, isoScales, isoEdges, isoLight, colormap)

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


.. toctree::
   :maxdepth: 2


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

