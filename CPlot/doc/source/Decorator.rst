.. CPlot.Decorator documentation master file

CPlot.Decorator: decoration of CPlot images using matplotlib
============================================================

Preamble
########

CPlot can generate images of meshes or flow fields. Those image can be
further enhanced with any additional matplotlib items using this module. 

This module is part of Cassiopee, a free open-source pre- and post-processor for CFD simulations.

To import CPlot.Decorator module::

   import CPlot.Decorator as Decorator


.. py:module:: CPlot.Decorator


List of functions
##################

**-- Actions**

.. autosummary::

   CPlot.display
   CPlot.render
   CPlot.delete
   CPlot.add
   CPlot.replace
   CPlot.finalizeExport




Contents
#########

    
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
    :param meshStyle: 0: white solid and red wireframe, 1: colored wireframe, 2: colored solid and wireframe, 3: cyan solid and black wireframe (default: 2)
    :type meshStyle: int
    :param solidStyle: 0: blue, 1: colored by zone, 3: white, 4: colored by zone outlined (default: 1)
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
    :param vectorProjection: 1 of vectors are projected on surface (default: 0)
    :type vectorProjection: 0 or 1
    :param colormap: 0-1: Blue2Red, 2-3: Green2Red, 4-5: BiColorRGB, 6-7: BiColorHSV, 8-9: Diverging, 10-11: TriColorRGB, 12-13: TriColorHSV, 14-15: Grey2White, 16-17: Viridis, 18-19: Inferno, 20-21: magma, 22-23: plasma, 24-25: nice blue (default: 0)
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
    :type isoScales: list of [string, int, float, float] or [string, int, float, float, float, float]
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
    :param exportResolution: resolution for export ("1920x1080")
    :type exportResolution: string
    :param zoneNames: optional list of zone names (same size as arrays, struct zones, then unstruct zones)
    :type zoneNames: list of strings
    :param renderTags: optional list of render tags (same size as arrays, struct zones, then unstruct zones)
    :type renderTags: list of strings
    :param frameBuffer: the number of the frame buffer we are rendering to for offscreen rendering
    :type frameBuffer: int
    :param offscreen: 0:  rendering, 1: offscreen rendering (osmesa), 2: offscreen rendering (openGL), 3: partial composite offscreen rendering (openGL), 4: final composite offscreen rendering (openGL), 5: partial composite offscreen rendering (osmesa), 6: final composite offscreen rendering (osmesa), 7: parallel offscreen rendering (osmesa) (default: 0)
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


.. toctree::
   :maxdepth: 2


Index
########

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

