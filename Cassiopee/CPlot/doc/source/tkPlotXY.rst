.. tkPlotXY documentation master file

:tocdepth: 2

tkPlotXY: plot of curves
=========================

Preamble
########

tkPlotXY is a 2D plotting library based on Matplotlib. The aim of tkPlotXY is to provide
to users an easier scriptable interface and a useful graphical interface for plotting data in 
arrays or pyTrees.

tkPlotXY uses preferentially 1D-data from pyTrees but in the scriptable interface, some other ways to
define datas are available and will be exposed in this document.

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

For use in a python script, you have to import tkPlotXY module::

   import tkPlotXY

.. py:module:: tkPlotXY


One line plot function
######################

tkPlotXY can be used with a single function.

.. py:function:: tkPlotXY.plot(a, varx, vary, rangex=None, rangey=None, export=None, ...)

    Plot 1D zones of a.

    :param a: input data
    :type a: [pyTree, base, zone, list of zones]
    :param varx: name of variable to plot in X
    :type varx: string
    :param vary: name of variable to plot in Y
    :type vary: string
    :param rangex: if not None, range for x variable. If None, automatic setting.
    :type rangex: None or list of two floats
    :param rangey: if not None, range for y variable. If None, automatic setting.
    :type rangey: None or list of two floats
    :param xlabel: if not None, name to display on x axis.
    :type xlabel: None or string
    :param ylabel: if not None, name to display on y axis.
    :type ylabel: None or string
    :param xformat: if not None, format for x axis values.
    :type xformat: None or string like '9.3f'
    :param yformat: if not None, format for y axis values.
    :type yformat: None or string like '9.3f'
    :param legends: if not None, the name of curves to be displayed in legend.
    :type legends: None or list of strings
    :param export: if None, interactive plot, otherwise name of export file.
    :type export: None or name of export file
    :param lineWidth: width of plot lines
    :type lineWidth: float (default: 1.5) or list of floats for each zone
    :param lineColor: color of plot lines. Html string (#FFFFFF) or color name ('black).
    :type lineColor: string or list of strings for each zone
    :param markerStyle: style of marker, 'none', '+', 'o', ...
    :type markerStyle: string or list of strings for each zone
    :param markerWidth: width of markers
    :type markerWidth: float (default: 6.5) or list of floats for each zone
    :param markerFaceColor: face color of markers. Html string (#FFFFFF) or color name ('black).
    :type markerFaceColor: string or list of string for each zone
    :param markerEdgeColor: edge color of markers. Html string (#FFFFFF) or color name ('black).
    :type markerEdgeColor: string

    *Example of use:*

    * `Plot 1d zones (pyTree) <Examples/CPlot/tkPlotXYSingleLine.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/tkPlotXYSingleLine.py

---------------------------------------------------------------------------


Usage with classes
###################

In case of the previous function doesnt provide enough flexibility, you can use 
tkPlotXY with classes that provide access to a full customization. Here goes 
two simple examples, an interactive one and a batch one.

    *Example of use:*

    * `Interactive tkPlotXY classes usage (pyTree) <Examples/CPlot/tkPlotXYInteractive.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/tkPlotXYInteractive.py

    * `Batch tkPlotXY classes usage (pyTree) <Examples/CPlot/tkPlotXYBatch.py>`_:

    .. literalinclude:: ../build/Examples/CPlot/tkPlotXYBatch.py


---------------------------------------------------------------------------


.. toctree::
   :maxdepth: 2


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
