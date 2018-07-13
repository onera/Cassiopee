.. tkCassiopee documentation master file


Preamble
########
    
tkCassiopee is the graphical user interface to Cassiopee modules.

The application is made from different applets. Applets are stored in tabs (Tree, 
State, Edge, Surf, ...).

To call the application from shell::
    
    cassiopee <file>
    


General guidelines
###################
    
Each applet can be triggered from top menu **Apps** or
by right clicking on a tab. It can be
discarded by left clicking on its name at top of frame or by pressing
CTRL+c.

The default settings of an applet can be modified in the preference file
by right clicking on the applet name and chosing 'save'.
Settings can be reset by chosing 'reset' in the menu.
An applet can be pinned, meaning that it will automatically
be opened at next restart.

The graphical window is called the CPlot window. You can rotate
your model using left mouse button. You can move your model
using left mouse button. You can tilt the model by pressing
CTRL+right mouse button. You can move the way your looking without
modifying your position by pressing CTRL+left mouse button.
View can be quickly centered by left double clicking.

A zone can be selected by using 
SHIFT + left mouse button, can be deactivated using 
SHIFT + right mouse button. Multiple selection is available
by pressing SHIFT + CTRL + left mouse button or by draging
mouse with SHIFT pressed.

Some CPlot functions are available through the **CPlot**
top menu. Other functions are described in the CPlot module
documentation.

Generally, functions of the applet work on selected zones.

One level-undo is available in the **CPlot** top menu or by pressing 
CTRL+z or by pressing the back arrow in the tool bar.

You can chose the location of displayed data (node/center) for
structured zones by clicking on 'change displayed location'
in the **CPlot** top menu.

Selecting 'Quit' in **File** top menu exits cassiopee. You can
also exits by pressing 'q' when in CPlot window.

Cassiopee:ToolBar
------------------

    The tool bar is situated just under the menu bar. It contains:

    - The save icon: click for quick save,
    - The undo icon: click for undo (only one level), 
    - The delete block icon: click will delete selected blocks,
    - The copy block icon: click will copy selected blocks,
    - The fit view icon: click to fit view to selection,
    - The select all blocks icon: select all zones, 
    - The view deactivated zones icon: enables to view deactivated zones as ghost grids,
    - The main tree icon: force view on main tree.    


Cassiopee:Tree
---------------
    
This tab gathers applets related to pyTree management.

    - **tkTree**: enables the visualization of pyTrees as text trees.
    In tkTree window, left mouse button
    selects zone, right mouse button deactivates/activates zone,
    double clicking on a node name enables to change its name,
    'Suppr' key deletes selected node, CTRL+e extend the window size,
    CTRL+r shrink the window.

    - **tkTreeOps**: enables to move selection to another base, reorder tree nodes and edit node values.

    - **tkCheckPyTree**: enables to check your pyTree with different level of checking (node conformity, zone names, boundary condition consistency...).

    - **tkFilter**: enables to select/activate zones that match the filter rule. 
    For filtering by name, you can use standard python regexp. 
    For instance, Zone.[5-12] will filter zones between Zone.5 and Zone.12.
    You can also filter zones by size (number of points), by multigrid level,
    by affected processor, by Chimera priority, or by formula. With this,
    you can specify a rule like '{Density}>0.1': zones that match the rule
    for at least one grid point will be selected.

    - **tkFamily**: enables to create zone families or BC families.

Cassiopee:State
----------------
    
This menu gathers applets related to state modification.

    - **tkState**: modify the global state of your pyTree (dimension, equations,...). 
    Pressing 'Set state' save state in your pyTree.

    - **tkPrefs**: enables to change the mesh display style, the solid display style, ... 
    You can add applets that open automatically
    each time cassiopee starts in the 'Auto open field'. For instance:
    'tkTree; tkBlock' will automatically open those two applets at next start.
    Save your preferences.

    - **tkPerfo**: settings to simplify the display and set the number of threads.

    - **tkContainers**: change containers. 
    By default, all functions work on 'GridCoordinates',
    'FlowSolution' and 'FlowSolution#Centers' containers. If you want
    functions to operate on other containers, change their name here.

    - **tkRuler**: measure your distance on your model using this applet. 
    First click on model set first point, second click indicates the distance.
    Click again on 'Measure mode' to end.

    - **tkFind**: find a given index in mesh.

    
Cassiopee:Edge
-----------------
    
    - **tkCanvas**: enables to create a canvas for mesh positioning or drawing. The position and size of the canvas can be modified.

    - **tkPoint**: enables to draw points.

    - **tkDraw**: enables to draw basic shapes.

    - **tkExtractEdges**: enables to extract edges generally from surfaces.

    - **tkMapEdge**: enables points redistribution (remap) on edges.

Cassiopee:Surf
---------------
    
This menu gathers applets related to surface creation/modification.

    - **tkBasicSurfs**: create basic surfaces (sphere, tetra, ...)

    - **tkText**: create a text.

    - **tkFixer2**: fix gaps in surfaces. 
    This applet enables manual closure of holes in surfaces. Select the contour of your hole, then click on 'Fix gap in contour'. You can bump the generated surface using the slider. This applet enables also automatic closure of all holes in a model. You must set surfaces defining a component into one zone, then click on 'Fix gap in patches'.

    - **tkBoolean**: perform boolean operation between surfaces.

    - **tkSculpt**: very basic sculpting tool.

    - **tkPaint**: change field values by painting.

    - **tkMapSurf**: perform surface remeshing by projecting an octree on surface.

    - **tkFilterSurfs**: enables to filter or inflate a surface.

    - **tkSurfaceWalk**: create meshes by walking on surfaces.

    - **tkProjection**: project a surface on another surface.
 
    
Cassiopee:Mesh
---------------
    
This menu gathers applets related to mesh creation/modification.

    - **tkCells**: enables mesh cell modification (suppress/refine cells).
    Click on a mode, then click on your mesh. Click again on the mode
    button when done.

    - **tkStretch**: for structured grids, stretch, refine or uniformize the grid in given direction. 
    For stretching, you must specify a grid step. 
    This step is enforced where you have last clicked.

    - **tkExtrusion**: create meshes by extrusion.

    - **tkTetraMesher**: create tetra meshes.

    - **tkTFI**: create meshes by transfinite interpolation.

    - **tkSmooth**: provides mesh smoothing.

    - **tkOctree**: enables octree mesh generation.

    - **tkCollarMesh**: create collar mesh (for junction bteween two solid bodies).

    - **tkBlader**: a blade dedicated mesher.

    - **tkMeshQual**: enables to compute various grid quality map and check for negative volume cells.

    - **tkMeshInfo**: enables to get various informations from your model (number of points, min/max of values...).

Cassiopee:Block
----------------
    
This menu gathers applets related to block creation/modification.

    - **tkBlock**: enables basic block operations (remove, copy...).
    Exterior faces returns the exterior faces of a zone as an unstructured
    zone. Close merge points in a mesh that are closer than epsilon, the
    resulting mesh connectivity is cleaned.

    - **tkTransform**: enables basic transformation of blocks (rotation, translation, ...). 
    When clicking on 'Translate by clicking', you
    must then click on a point of the zone to translate, then on 
    the destination point.

    - **tkNGon**: preforms NGon (polyedral) operations.

    - **tkSplit**: enables splitting or join operations on a block. 
    'Splitsize' splits each zone in order to get the required number of points.

    'SplitMP' eliminates multiple point junction in a structured mesh.
    
    'SplitConnexity' identifies connex parts in an unstructured block.

    - **tkReorder**: enables to reorder a zone. 
    Unstructured zones are reordered in order to have normals with the same orientation on each zone. 
    Structured zones are reordered by exchanging i- and j- numerotation.

Cassiopee:BC
--------------

This menu gathers applets related to boundary conditions 
creation/modification.

    - **tkBC**: enables to set interactively the boundary conditions.
    'View Mesh/BC' enables to visualize the boundary conditions
    of a certain type. 'View undefined BC' shows boundary conditions
    that are lacking. By clicking on BC, you can then set a
    boundary to a certain type using the 'setBCWith' button.
    'ConnectMatch' automatically computes the matching boundary condition
    in your model.

    - **tkChimera**: perform hole cutting, overlap optimization between overset grids.

    - **tkIBC**: create data for immersed boundary conditions.

    - **tkExtractBC**: extract a certain type of BC to a zone.

Cassiopee:Motion
------------------

    - **tkRigidMotion**: enables definition of rigid motions.

    - **tkTime**: manage time for motion visualization.

    
Cassiopee:Solver
-----------------
    
This menu gathers applets related to solvers.

    - **tkInit**: initialize solutions or wall distance.

    - **tkDistributor**: distributes blocks over processors. 
    Enter the number of processors and the weight of communication cost relative to
    solver cost per iteration. Click on 'Distribute tree'. Distribution stats
    are available when clicking on stats. 'Set Proc Field' creates a field
    in zones containing the attributed processor number for each zone.

    - **tkDist2Walls**: computes wall distance.

    - **tkCassiopeeSolver**: export CGNS files suitable for Cartesian solver.

    - **tkElsaSolver**: export CGNS files suitable for elsAxdt.

Cassiopee:Post
---------------
    
This menu gathers applets related to post-processing.

    - **tkVariables**: add/rm variables in your solution pyTree.
    - **tkExtractMesh**: interpolate solution from one mesh to another.
    - **tkStream**: extract streamlines.
    - **tkIsoLine**: extract an iso-line or a set of isolines to a zone.
    - **tkIsoSurf**: extract an iso-surface to a zone.
    - **tkInteg**: perform field integration on surfaces.

Cassiopee:Visu
----------------
    
This menu gathers applets related to pyTree visualization.

    - **tkView**: perform view settings.
    - **tkPlot**: perform 1D plot of data.
    - **tkSlice**: extract/view slices in mesh.
    - **tkCelln**: enables to display the location of interpolated, 
    blanked points described by a 'cellN' or 'cellNF' field.
    - **tkBackground**: add a background.

Cassiopee:Render
-----------------

    - **tkRenderSet**: enables to chose the color and material of each zone.
    - **tkStereo**: enable the stereo anaglyph mode.
    - **tkEffects**: enable special effects for view such as shadow, DOF.
    - **tkDemo**: rotate or move around your model automatically. Chose speed in slider.
    - **tkPovRay**: if povray is installed on your computer, you can use this applet to raytrace your scene using povray. 
    Chose the name of your output (used for image and pov output), chose your background
    and the size of output image, then click on 'Render scene'.
    - **tkLuxRender**: exports file for LuxRender.

For gurus
############
    
tkCassiopee is a fully modular GUI. You can add your own 
applet easily. First copy the file CPlot/apps/tkPersonalSample.py
to tkMyApplet.py. Then add your own buttons and functions in this file.
Finally, in cassiopee, go to top menu **Tools** and
**Add a personal app**. Enter the file name. Your applet
will appear in this menu.

