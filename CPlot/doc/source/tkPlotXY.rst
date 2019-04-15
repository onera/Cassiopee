.. tkPlotXY documentation master file

tkPlotXY : plot of curves
==========================

Preamble
########

tkPlotXY is a 2D plotting library based on Matplotlib. The aim of tkPlotXY is to provide
to users an easier scriptable interface and a useful graphical interface in the mean time.
This documentation focuses only on the scriptable interface. To know more about its graphical interface,
some tutos will soon be available.

tkPlotXY uses preferentially 1D-data from pyTrees but in the scriptable interface, some other ways to
define datas are available and will be exposed in this document.

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use in a python script, you have to import tkPlotXY module::

   import tkPlotXY

.. py:module:: tkPlotXY

List of classes
##################

tkPlotXY is based on classes. Some of them are internal classes used for display. They are not documented here. It has to be remarked that some classes have a 'TK' suffix at the end
of their name. These classes are equivalent to the one without suffix, but they have been developped to work inside the tkInter context. It means that for python scripting
only classes without the suffix 'TK' should be used.

**-- Classes**

.. autoclass:: GraphEditor
.. autoclass:: Desktop
.. autoclass:: Graph
.. autoclass:: Axis
.. autoclass:: DirAxis
.. autoclass:: Legend
.. autoclass:: AxisGrid
.. autoclass:: Grid
.. autoclass:: LevelGrid
.. autoclass:: Curve
.. autoclass:: SubPlotParams
.. autoclass:: TightLayout
.. autoclass:: Movie

GraphEditor
###########
.. py:class:: tkPlotXY.GraphEditor

An object of class GraphEditor allows you to create a Desktop. Accessing to the Desktop will give you the possibility to plot all your graphs. Moreover, the Desktop contains all the data that can be used to generate plots.
For python scripting interface, the first step is to create an object of this class GraphEditor and to access its Desktop. This is performed by the function openGraphEditor that directly returns the Desktop. Then the data can be added to this Desktop object. Finally it is used to generate all the graphs.

The first step is to create a graphEditor and to get its Desktop using :

.. py:function:: tkPlotXY.openGraphEditor

.. autosummary:: tkPlotXY.openGraphEditor

.. code-block:: python

   import tkPlotXY as tkP
   # Create a graphEditor
   graphDesktop = tkP.openGraphEditor(None)

Desktop
###########
.. py:class:: tkPlotXY.Desktop

The Desktop deals with the data management and the graph plotting.

Data management
----------------

The data can be loaded by a Desktop object from a pyTree or from a dictionnary. Only the 1D-array data from a pyTree will be loaded while the data loaded from a dictionnary has to be compliant with the following structure :
{Base/Zone (string) : {Variable name (string): data (array)}}

Several methods are available to set, update or even remove data :



.. automethod:: tkPlotXY.Desktop.addZone
.. automethod:: tkPlotXY.Desktop.setData
.. automethod:: tkPlotXY.Desktop.replaceZone
.. automethod:: tkPlotXY.Desktop.deleteZoneFromData


.. code-block:: python

    import numpy as np
    import tkPlotXY as tkP
    # Create a graphEditor
    graphDesktop = tkP.openGraphEditor(None)
    # Generate data
    t = np.arange(0., 5., 0.002)
    dataFromDict = {'Zone1':
                            {
                                'Iteration':t,
                                'Residual':np.sin(t),
                                'Cf':np.sin(t/2),
                                'Debit':t*t
                            }
                    }
    # Set data
    graphDesktop.setData(dataFromDict)
    # Display data
    for zone in graphDesktop.data.keys():
        for var in graphDesktop.data[zone].keys():
            print zone, ' : ', var, ' : ',graphDesktop.data[zone][var]


.. code-block:: python

    import tkPlotXY as tkP
    # Create a graphEditor
    graphDesktop = tkP.openGraphEditor(None)
    # Generate data with a Lamb vortex (pyTree) -
    import Generator.PyTree as G
    import Initiator.PyTree as I
    import Post.PyTree as P
    import Converter.PyTree as C

    NI = 100; NJ = 100
    HI = 50./(NI-1); HJ = 50./(NJ-1)
    tree = G.cart((0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
    tree = I.initLamb(tree, position=(7.,7.), Gamma=2., MInf=0.8, loc='centers')
    tree = P.isoSurfMC(tree, 'CoordinateZ', 0.5)
    tree = P.isoSurfMC(tree, 'CoordinateY', 0.7)

    # Save generated data as cgns
    C.convertPyTree2File(tree,'vortex_slice.hdf')
    # Load data as pyTree
    tree=C.convertFile2PyTree('./vortex_slice.hdf')

    # Set data
    graphDesktop.setData(tree)
    # Display data
    for zone in graphDesktop.data.keys():
        for var in graphDesktop.data[zone].keys():
            print zone, ' : ', var, ' : ',graphDesktop.data[zone][var]

Graph creation
--------------

Once data have been loaded into the Desktop, you can create as many graphs as you need with the Desktop object. After that, all the drawing will be driven by the graph object itself. This is performed by the method :

.. automethod:: tkPlotXY.Desktop.createGraph

.. code-block:: python

    # Create First Graph
    graph_0 = graphDesktop.createGraph('MyFirstGraph','1:1')
    # Create Second Graph
    graph_1 = graphDesktop.createGraph('MySecondGraph','2:1')

Graph
######

.. autoclass:: Graph

Creating a Graph object will automatically generate an Axis, a Grid and a Legend objects for each plots on the graph.
Only curves have to be created and then attached to a given graph.

Then each object can be configured. To do so, it is mandatory to access these objects thanks to the graph.

All of these actions are described in the concerned item section (Curve, Axis, Grid or Legend).

Configure the Graph object
---------------------------

Some times, using a matricial Graph (for instance '2:2') will provide you an unacceptable drawing. Indeed, the axis label of a plot may be overlapped by the plot below. For all these reasons, you may be interested in advanced configuration for your Graph object such as positionning, padding, margin ...

Two possibilities are available: TightLayout or SubPlotParams.

In order to improve your drawing using SubPlotParams, please use the method:

.. automethod:: tkPlotXY.Graph.updateSubPlotParams

where *params* is a dictionnary such that:

.. code-block:: python

   params = {'left':...,'right':...,'top':...,'bottom':...,'hspace':...,'wspace':...,'isActive':True}

For example, let us create a Graph object *graph_3* with 4 plots inside ('2:2') and let us try to improve the poisitionnment of this graph with SubPlotParams:

.. code-block:: python

   # Graph creation
   graph_3 = graphDesktop.createGraph('MyThirdGraph','2:2')
   # Improving drawing of the Graph thanks to SubPlotParams
   graph_3.updateSubPlotParams({'isActive':True,'right':0.97,'top':0.97,'wspace':0.3})

The other way to improve this kind of drawing is to use TightLayout:

.. automethod:: tkPlotXY.Graph.updateTightLayout

where:

.. code-block:: python

   params = {'isActive':True, 'pad':..., 'hpad':..., 'wpad':...}

For example, let us create a Graph object *graph_3* with 4 plots inside ('2:2') and let us try to improve the poisitionnment of this graph with TightLayout.

.. code-block:: python

   # Graph creation
   graph_3 = graphDesktop.createGraph('MyThirdGraph','2:2')
   # Improving drawing of the Graph thanks to TightLayout
   graph_3.updateTightLayout({'isActive':True,'pad':1.1,'hpad':0.1,'wpad':0.1})

Curve
######

.. autoclass:: Curve

Creating a curve
------------------

To create a curve, one has just to create an object of class Curve. All the settings, listed in the section 'Editing a curve', can be already configured during the creation of the object. For instance:

.. code-block:: python

   curve_0 = Curve(zone=['Base/cart'],varx='CoordinateX',vary='Density_FlowSolution' , line_color='#7f00ff' , marker_face_color='#7f00ff' ,marker_edge_color='#7f00ff' )

Editing a curve
-----------------

To edit a curve, you can use the method :

.. automethod:: tkPlotXY.Curve.setValue

with *variable* according to the following tab:


.. csv-table:: **Available variables to set a curve**
   :header: "Variable", "Allowed values", "Description"
   :widths: 30, 30, 30

   "zone                 ","`List of zones`","`List of zones that should be plotted`"
   "varx                 ","`Variable name (string)`", "`X-coordinate variable name (LaTeX available with $$)`"
   "vary                 ","`Variable name (string)`", "`Y-coordinate variable name (LaTeX available with $$)`"
   "line_color           ","`Html color code (string)`", "`Color used to plot the line`"
   "line_style           ","`'solid', 'dashed', 'dashdot', 'dotted', 'None'`", "`Style of line to use`"
   "line_width           ","`(float)`", "`Width of the line`"
   "marker_style         ","`'none', 'plus', 'star', 'pixel', 'point', 'star3_down', 'star3_up', 'star3_left', 'star3_right', 'triangle_left', 'triangle_right', 'diamond', 'hexagon2', 'triangle_up', 'hline', 'thin_diamond', 'hexagon1', 'circle', 'pentagon', 'square', 'triangle_down', 'x'`", "`Type of marker to use`"
   "marker_size          ","`(float)`", "`Size of the marker`"
   "marker_edge_color    ","`Html color code (string)`","`Color of the edge of the marker`"
   "marker_edge_width    ","`(float)`","`Width of the edge of the marker`"
   "marker_face_color    ","`Html color code (string)`","`Color of the face marker`"
   "marker_sampling_start","`(int)`","`Index on data to start plotting markers`"
   "marker_sampling_end  ","`(int)`","`Index on data to stop plotting markers`"
   "marker_sampling_step ","`(int)`","`Step between index to plot markers`"
   "legend_label         ","`(string)`","`Name of the curve to display in the legend`"
   "legend_display       ","`(bool)`","`Display the current curve in the legend`"
   "visible              ","`(bool)`","`Hide or show the curve`"
   "axis                 ","`(int)`","`Axis in which the curve has to be plotted`"

For instance to edit a curve to a dashed curve:

.. code-block:: python

   curve_0.setValue('line_style','dashed')

A curve can be edited all the time. The graph has just to be updated after the modification of the curve property.

Adding a curve to a given plot on a given graph
-------------------------------------------------

To attach a curve to a given plot inside a given graph, use the method:

.. automethod:: tkPlotXY.Graph.addCurve

where *iCurSubGraph* identifies the plot inside the graph thanks to its matricial position. For instance, to add a curve on the second line and first column of the graph *graph_1*:

.. code-block:: python

   graph_1.addCurve('2:1',curve_0)

Axis
######

.. autoclass:: tkPlotXY.Axis
.. autoclass:: tkPlotXY.DirAxis

While a Graph object  is created, an axis system (X and Y DirAxis) is generated for each plot on the graph. This system of X and Y DirAxis is an Axis object.


Access the Axis system
--------------------------

To get the Axis system of  given plot on a given graph, use the method on your graph:

.. automethod:: tkPlotXY.Graph.getAxis

where *iCurSubGraph* identifies the plot inside the graph thanks to its matricial position. Moreover, in case of a multiple axis plot (2 or more Y DirAxis on the same plot for instance, see *Multiple axis system section*), you can specify the number identifying your axis system using *ind*. Note that the original axis system has the index 0 and then the index is increased for each new axis system.

.. code-block:: python

   axis_2 = graph_1.getAxis('2:1',ind=0) # Equivalent to axis_2 = graph_1.getAxis('2:1')

Multiple axis system
----------------------

On the same plot, you can use multiple axis system. You can decide to twin your X or Y DirAxis or even to create a new independant axis system. to generate your new axis system (twin or independant), use the method:

.. automethod:: tkPlotXY.Graph.addAxis

where *iCurSubGraph* identifies the plot inside the graph thanks to its matricial position, *shared* can take the value : 'x','X','y','Y' or None. If None is used, then an independant axis system will be created. If an other value is used, then *axis* allows you to specify the index of the axis system you want to clone. Remember, this index starts at 0 for each plot and is then locally increased for each new axis system in a plot. Or you can directly give the object Axis you want to clone by using the parameter 'axis'. This function returns the newly created axis object.

.. code-block:: python

   # Get axis_2
   ind_axis_2 = axis_2.getInd() # returns the index of the current axis
   axis_2 = graph_1.getAxis('2:1', ind_axis_2) # Equivalent to axis_2 = graph_1.getAxis('2:1') because here ind_axis_2 = 0 !

   # Twin X DirAxis of axis_2
   axis_3 = graph_1.addAxis('2:1',shared='x',axis=axis_2) # equivalent to "axis_3 = graph_1.addAxis('2:1',shared='x',ind=ind_axis_2)"
   ind_axis_3 = axis_3.getInd()

Changing the Axis of a curve
----------------------------

Once a curve has been added to a given plot on a given graph and that this plot is composed of several axis, then it is possible to change the axis where the curve will be drawn into the given plot. To do so, you just need to edit the attribute *axis* of your Curve object. You can either use the axis object itself or its index.

.. code-block:: python

   curve_3.setValue('axis',axis_3)
   # equivalent to
   curve_3.setValue('ind_axis',ind_axis_3)

Editing the Axis system
------------------------

You can edit it by accessing to the Axis object and using the method:

.. automethod:: tkPlotXY.Axis.setValue

where *axis* has to be 'x' or 'y'.

Or you can directly access the X (resp. Y) DirAxis object by using the attribute *x* (resp. *y*) of the class Axis that will return the X DirAxis (resp. Y DirAxis) and then you can use the method:

.. automethod:: tkPlotXY.DirAxis.setValue

with *variable* according to the following tab:



.. csv-table:: **Available variables to set a DirAxis**
   :header: "Variable", "Allowed values", "Description"
   :widths: 30, 30, 30

   "axis_logscale        ","`(bool)`","`Use logscale for selected DirAxis`"
   "axis_autoscale       ","`(bool)`","`Use auto-scaling for selected DirAxis`"
   "axis_min             ","`(float)`","`Minimum range to plot for the selected DirAxis`"
   "axis_max             ","`(float)`","`Maximum range to plot for the selected DirAxis`"
   "axis_label           ","`Label name (string)`","`Label for the selected DirAxis (LaTeX formula are available with $$)`"
   "axis_inverted        ","`(bool)`","`Invert the orientation for the selected DirAxis`"
   "axis_visible         ","`(bool)`","`Show or hide the selected DirAxis`"
   "axis_position        ","`For X-DirAxis : 'top','bottom','both' and for Y-DirAxis : 'left','right','both'`","`Position where to plot the axis line, its ticks and the  label for selected DirAxis`"
   "axis_offset          ","`(float)`","`Introduce an offset for the axis line, its ticks and the label for selected DirAxis`"
   "axis_label_fontsize  ","`(float)`","`Set the size of the label font for the selected DirAxis`"

For instance, setting logscale on the Y-DirAxis of the *axis_3* previously defined as a twin X-DirAxis of *axis_2* (See *Multiple axis system* section)

.. code-block:: python

   axis_3.y.setValue('axis_logscale',True) # equivalent to axis_3.setValue('y','axis_inverted',True)

Grid
#######

.. autoclass:: tkPlotXY.Grid
.. autoclass:: LevelGrid
.. autoclass:: AxisGrid

Once an axis system is created, a Grid object is attached to this new axis system. It means that for each Axis object there exists a unique associated Grid object. This Grid object is composed of two LevelGrid objects : *major* and *minor* which corresponds to the main grid and the second grid. Each LevelGrid object contains two AxisGrid : *X* and *Y*. To put it in a nutshell, a Grid object describes 4 AxisGrid objects : *major X*, *major Y*, *minor X* and *minor Y*.

Since Grid objects are automatically generated during the creation of Axis object, there is no need to create Grid object. It is just needed to be able to access it.

Access a Grid object
---------------------

To access the associated Grid object to a given Axis on a given plot inside a given Graph, use the method:

.. automethod:: tkPlotXY.Graph.getGrid

where *iCurSubGraph* identifies the plot inside the graph thanks to its matricial position. Moreover, in case of a multiple axis plot (2 or more Y DirAxis on the same plot for instance, see *Multiple axis system section*), you can specify the number identifying your axis system using *ind* or directly specify the axis object using *axis*. Note that the original axis system has the index 0 and then the index is increased for each new axis system.

For instance, to get the Grid associated to the Axis *axis_3* previously defined:

.. code-block:: python

   grid_3 = graph_1.getGrid('2:1',ind=1)
   # is equivalent to
   grid_3 = graph_1.getGrid('2:1',axis=axis_3)

If needed, you can access directly the LevelGrid object using the attributes *minor* and *major* of class Grid and then you can get the *AxisGrid* using the attributes *x* and *y* of class LevelGrid.

.. code-block:: python

   grid_3 = graph_1.getGrid('2:1',axis=axis_3)
   grid_3_majorX = grid_3.major.x
   grid_3_majorY = grid_3.major.y
   grid_3_minorX = grid_3.minor.x
   grid_3_minorY = grid_3.minor.y

Editing a Grid object
---------------------

You can edit a Grid object using the method:

.. automethod:: tkPlotXY.Grid.setValue

where *level* ('major' or 'minor' expected) and *direction* ('x' or 'y' expected) identify the AxisGrid to edit

or the LevelGrid object using the method:

.. automethod:: tkPlotXY.LevelGrid.setValue

where *direction* ('x' or 'y' expected) identifies the AxisGrid to edit

or directly the AxisGrid object using the method:

.. automethod:: tkPlotXY.AxisGrid.setValue

Authorized *variable* and *value* are described in the following tab:

.. csv-table:: **Available variables to set an AxisGrid**
   :header: "Variable", "Allowed values", "Description"
   :widths: 30, 30, 30

   "display          ","`(bool)`",`Show or hide the AxisGrid`"
   "grid_color       ","`Html color code (string)`",`Modify the color of the AxisGrid`"
   "grid_style       ","`'solid', 'dashed', 'dashdot', 'dotted', 'None'`",`Modify the line style for the AxisGrid`"
   "grid_width       ","`(float)`",`Modify the line width for the AxisGrid`"
   "grid_tick_number ","`(int)`",`Change the number of ticks on the AxisGrid`"
   "grid_tick_size   ","`(float)`",`Change the size of the ticks`"

For instance, to add the major grid as dashed lines for X and Y axis on the Grid grid_3 previously defined:

.. code-block:: python

   grid_3.setValue('major','x','display',True)
   grid_3.setValue('major','x','grid_style','dashed')

   grid_3.setValue('major','y','display',True)
   grid_3.setValue('major','y','grid_style','dashed')

Legend
#######

.. autoclass:: tkPlotXY.Legend

While a Graph object  is created, a Legend object is associated to each plot of the Graph. There is no need to create manually a Legend object. You just need to access it in order to edit it.

Accessing a Legend object
----------------------------

To access a Legend object of a given plot inside a given Graph, use the method:

.. automethod:: tkPlotXY.Graph.getLegend

where *iCurSubGraph* identifies the plot inside the graph thanks to its matricial position.

For example, to get the Legend of the plot on the second line of graph_1:

.. code-block:: python

   legend_2 = graph_1.getLegend('2:1')

Editing a Legend object
------------------------

To edit a Legend object, you can use the method:

.. automethod:: tkPlotXY.Legend.setValue

Authorized *variable* and *value* are described in the following tab:

.. csv-table:: **Available variables to set a Legend**
   :header: "Variable", "Allowed values", "Description"
   :widths: 30, 30, 30

   "legend_display                 ","`(bool)`","`Show or hide the legend box`"
   "legend_title                   ","`(string)`","`Title of the legend`"
   "legend_border_width            ","`(float)`","`Width of the legend box border`"
   "legend_border_color            ","`Html color code (string)`","`Color of the border of the legend box`"
   "legend_background_color        ","`Html color code (string)`","`Background color of the legend box`"
   "legend_background_color_active ","`(bool)`","`Use or not transparency as background for the box of the legend`"
   "legend_position                ","`'best', 'upper left', 'upper center', 'upper right', 'center left', 'center', 'center right', 'lower left', 'lower center', 'lower    right'`","`Position of the legend box`"
   "legend_ncol                    ","`(int)`","`Number of columns to display the legend`"
   "legend_label_weight            ","`'normal','bold'`","`Use bold font for the curves name`"
   "legend_label_style             ","`'normal','italic'`","`Use italic font for the curves name`"
   "legend_label_size              ","`(float)`","`Font size for the curves name`"
   "legend_label_color             ","`Html color code (string)`","`Font color used for the name of the curves`"
   "legend_title_weight            ","`'normal','bold'`","`Use bold font for the legend title`"
   "legend_title_style             ","`'normal','italic'`","`Use italic font for the legend title`"
   "legend_title_size              ","`(float)`","`Font size for the legend title`"
   "legend_title_color             ","`Html color code (string)`","`Font color used for the legend title`"

For example, the following code set the title of Legend legend_2, that has been previously defined, with a bold font:

.. code-block:: python

   legend_2.setValue('legend_title_weight','bold')

Update, View and Save your graph
#################################

Now that ve have listed all elements that can be used to configure your plots, let us address the main objective of your python script: visualize your plot. First of all, you will need to update them to take into account all the modifications that have been added. Then you can either display on screen or save the figure.

Update figures
------------------

To update a given plot on a given Graph, use the method:

.. automethod:: tkPlotXY.Graph.updateGraph

where *iCurSubGraph* identifies the plot inside the graph thanks to its matricial position.

For instance, in order to update the plot on second line of *graph_2*:

.. code-block:: python

   graph_2.updateGraph('2:1')

If you want to update all the plots of a given Graph, you can use the method:

.. automethod:: tkPlotXY.Graph.drawFigure

For example, to update the two plots on our *graph_2*:

.. code-block:: python

   graph_2.drawFigure()

Please, do not be confused, the method **drawFigure** does not display the figure !

Display on screen
------------------

In order to display on screen a Graph, use the method:

.. automethod:: tkPlotXY.Graph.showFigure

Do not forget to add a call to **time.sleep()** after this command in order to let the display active more longer than just a pop-up !

.. code-block:: python

   import time
   ####
   #...
   ####
   graph_2.showFigure()
   time.sleep(5.)  # stop for 5 sec. here

Save figures
--------------

Use the method **save** of the class Graph to save your drawing in a file:

.. automethod:: tkPlotXY.Graph.save

where *path* is the output file path and format is the format you wish to use to save your figure (available formats are : 'emf', 'eps', 'pdf', 'png', 'ps', 'raw', 'rgba', 'svg', 'svgz').

.. code-block:: python

   graph_2.save('/home/User/Example/myNiceGraph.png')


Extra usages
#############

Load & Save configurations script
---------------------------------

Configurations scripts are python scripts automatically generated by the GUI. They can be loaded to ease the proccess of tunning your plot. For instance, if you have a plot that you often draw, instead of re-creating each time your drawing, you can simply load your configuration and use it as a base. Moreover, this configurations scripts are created by the graphical user interface and you will only need to adapt a few elements (removing the 'TK' suffix of all the classes basically) to use them directly as script without the GUI.


Complete example
#################

.. code-block:: python

   import tkPlotXY as tkP
   # Create a graphEditor
   graphDesktop = tkP.openGraphEditor(None)
   # Generate data with a Lamb vortex (pyTree) -
   import Generator.PyTree as G
   import Initiator.PyTree as I
   import Post.PyTree as P
   import Converter.PyTree as C
   import time,os
   cwd = os.getcwd()

   DEBUG_CHECKDATA = True

   # Creating a vortex
   NI = 100; NJ = 100
   HI = 50./(NI-1); HJ = 50./(NJ-1)
   tree = G.cart((0.,0.,0.), (HI,HJ,1.), (NI,NJ,2))
   tree = I.initLamb(tree, position=(7.,7.), Gamma=2., MInf=0.8, loc='centers')
   tree = P.isoSurfMC(tree, 'CoordinateZ', 0.5)
   tree = P.isoSurfMC(tree, 'CoordinateY', 0.7)

   # Save generated data as cgns
   C.convertPyTree2File(tree,'vortex_slice.hdf')
   # Load data as pyTree
   tree=C.convertFile2PyTree('./vortex_slice.hdf')

   ########################### Set data
   graphDesktop.setData(tree)
   if DEBUG_CHECKDATA:
       for z in graphDesktop.data.keys():
           print '*-'*15
           print 'Zone : ',z
           for k in graphDesktop.data[z].keys():
               print '---> Var : ',k

   ########################### Graph creation
   # Create First Graph
   graph_0 = graphDesktop.createGraph('MyFirstGraph','1:1')
   # Create Second Graph
   graph_1 = graphDesktop.createGraph('MySecondGraph','2:1')

   ########################### Curve creation
   # Create the first curve
   curve_0 = tkP.Curve(zone=['Base/cart'],varx='CoordinateX',vary='Density@FlowSolution' , line_color='#7f00ff' , marker_face_color='#7f00ff' ,marker_edge_color='#7f00ff' )

   # Create the second curve
   curve_1 = tkP.Curve(zone=['Base/cart'],varx='CoordinateX',vary='MomentumZ@FlowSolution' , line_color='#0404B4' , marker_face_color='#0404B4' ,marker_edge_color='#0404B4' )

   # Create the third curve
   curve_2 = tkP.Curve(zone=['Base/cart'],varx='CoordinateX',vary='MomentumX@FlowSolution' , line_color='#FF00FF' , marker_face_color='#FF00FF' ,marker_edge_color='#FF00FF' )

   # Create the fourth curve
   curve_3 = tkP.Curve(zone=['Base/cart'],varx='CoordinateX',vary='MomentumY@FlowSolution' , line_color='#FFBF00' , marker_face_color='#FFBF00' ,marker_edge_color='#FFBF00' )

   ########################### Attaching curves to graph
   # First curve to graph_0
   graph_0.addCurve('1:1',curve_0)
   # Second curve to graph_1 first line
   graph_1.addCurve('1:1',curve_1)
   # Third curve to graph_1 second line
   graph_1.addCurve('2:1',curve_2)
   # Fourth curve to graph_1 second line
   graph_1.addCurve('2:1',curve_3)

   ########################### Editing curves
   ## Name for the legend
   curve_0.setValue('legend_label','Density')
   curve_1.setValue('legend_label','MomentumZ')
   curve_2.setValue('legend_label','MomentumX')
   curve_3.setValue('legend_label','MomentumY')

   ## curve_3 : dashed
   curve_3.setValue('line_style','dashed')

   ## curve_2 : add markers
   curve_2.setValue('marker_style','circle')
   curve_2.setValue('marker_sampling_step',20) # 1 marker over 20

   ########################### Axis properties
   ## 1/- Get axis
   axis_0 = graph_0.getAxis('1:1') # ind = 0
   ind_axis_0 = axis_0.getInd()
   # #
   axis_1 = graph_1.getAxis('1:1') # ind = 0
   ind_axis_1 = axis_1.getInd()
   #
   axis_2 = graph_1.getAxis('2:1') # ind = 0
   ind_axis_2 = axis_2.getInd()
   ## 2/- Twining X axis for plot 2:1 on graph_1
   axis_3 = graph_1.addAxis('2:1',shared='x',axis=axis_2) # equivalent to "axis_3 = graph_1.addAxis('2:1',shared='x',ind=ind_axis_2)""
   ind_axis_3 = axis_3.getInd()

   # Set the position of axis label
   axis_2.setValue('y','axis_position','left')
   axis_3.setValue('y','axis_position','right')
   # Change the label text
   axis_1.setValue('y','axis_label','$\\rho W$')
   axis_2.setValue('y','axis_label','$\\rho U$')
   axis_3.setValue('y','axis_label','$\\rho V$')

   # ########################### Changing axis for curve_3
   curve_3.setValue('axis',axis_3) # equivalent to "curve_3.setValue('ind_axis',ind_axis_3)"

   ########################### Editing Grid properties
   # Get the grid objects
   grid_0 = graph_0.getGrid('1:1',ind=ind_axis_0)   # equivalent to "grid_0 = graph_0.getGrid('1:1',axis=axis_0)"
   grid_1 = graph_1.getGrid('1:1',ind=ind_axis_1)   # equivalent to "grid_1 = graph_1.getGrid('1:1',axis=axis_1)"
   grid_2 = graph_1.getGrid('2:1',axis=axis_2)   # equivalent to "grid_2 = graph_1.getGrid('2:1',ind=ind_axis_2)"
   grid_3 = graph_1.getGrid('2:1',axis=axis_3)   # equivalent to "grid_3 = graph_1.getGrid('2:1',ind=ind_axis_3)"
   # Display a solid grid for the major grids on X & Y axis_2
   grid_2.major.x.setValue('grid_style','solid')
   grid_2.major.y.setValue('grid_style','solid')

   ########################### Editing Legend properties
   # Get the legend objects
   legend_0 = graph_0.getLegend('1:1')
   legend_1 = graph_1.getLegend('1:1')
   legend_2 = graph_1.getLegend('2:1')
   # Hide legend_1
   legend_1.setValue('legend_display',False)
   # Reduce legend label size
   legend_2.setValue('legend_label_size',8)
   # Increase legend title size and set its font as bold
   legend_2.setValue('legend_title_size',10)
   legend_2.setValue('legend_title_weight','bold')
   # Position legend in the lower right corner
   legend_2.setValue('legend_position','lower right')

   ########################### Use SubPlotParams
   params = {'left':0.12,'right':0.87,'top':0.90,'bottom':0.12,'isActive':True,'hspace':0.3}
   graph_1.updateSubPlotParams(params)

   ########################### Update, view & save
   # Update
   graph_1.drawFigure()
   # Display
   graph_1.showFigure()
   # Wait
   time.sleep(5.)
   # Save
   graph_1.save(os.path.join(cwd,'MyNiceFigure.png'))

.. toctree::
   :maxdepth: 2


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
