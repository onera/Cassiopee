"""Simple plotter for CFD.
"""
__version__ = '4.2'
__author__ = "Christophe Benoit, Stephanie Peron, Pascal Raud, Matthieu Soismier, Bertrand Michel"

# cplot and cplotOSMesa must not be imported at the same time

import KCore.kcore as KCore
import KCore.Vector as Vector
from . import ColorMaps

import os.path
import time, numpy
from sys import version_info

import importlib
cplotm = None
# Get the CPlot module (cplot or cplotOSMesa)
# IN: offscreen: offscreen rendering mode (0=default, 1/5/6/7=OSMesa modes, 2/3/4=OpenGL)
# OUT: None
def getModule(offscreen=0):
    global cplotm
    if cplotm is None:
        if offscreen == 1 or offscreen == 5 or offscreen == 6 or offscreen == 7:
            cplotm = importlib.import_module("CPlot.cplotOSMesa")
        else:
            cplotm = importlib.import_module("CPlot.cplot")
    return None

__timeStep__ = 0.01
__slot__ = None

#==============================================================================
# -- configuration --
#==============================================================================

# Configure the type of renderer used internally
# IN: useRender: renderer type (direct rendering, display lists, or VBO)
# IN: offscreen: offscreen rendering mode (default=0)
# OUT: None
def configure(useRender, offscreen=0):
    """Configure CPlot for direct rendering (cplot.useDirect), display Lists (cplot.useDL)
        or VBO (cplot.useVBO)"""
    getModule(offscreen)
    cplotm.configure(useRender)

# Check if direct rendering is available on host
# IN: None
# OUT: True if direct rendering is available, False otherwise
def hasDirectRendering():
    """Detect if direct rendering is available on host (may not work on windows)"""
    try:
        import subprocess
        out = subprocess.check_output('glxinfo | grep "direct rendering"', shell=True, stderr=subprocess.STDOUT)
        if "Yes" in out: return True
        else: return False
    except: return False

#==============================================================================
# -- display --
#==============================================================================

# Main display function for CFD arrays
# IN: arrays: list of arrays to display
# IN: dim: dimension (3=3D, 2=2D, 1=1D, -1=auto)
# IN: mode: display mode (-1=auto)
# IN: scalarField: scalar field index for coloring
# IN: vectorField1/2/3: vector field components
# IN: displayBB: display bounding box (0/1)
# IN: displayInfo: display info overlay (0/1)
# IN: displayIsoLegend: display isosurface legend (0/1)
# IN: meshStyle: mesh display style
# IN: solidStyle: solid display style
# IN: scalarStyle: scalar display style
# IN: vectorStyle: vector display style
# IN: vectorScale: vector scaling factor
# IN: vectorDensity: vector density
# IN: vectorNormalize: normalize vectors (0/1)
# IN: vectorShowSurface: show surface with vectors (0/1)
# IN: vectorShape: vector shape type
# IN: vectorProjection: vector projection mode
# IN: colormap: colormap index
# IN: colormapC1/C2/C3: colormap color stops
# IN: colormapC: colormap color array
# IN: niso: number of isosurfaces
# IN: isoEdges: display isosurface edges (0/1)
# IN: isoScales: isosurface scale factors
# IN: win: window position (x,y)
# IN: posCam: camera position (x,y,z)
# IN: posEye: camera eye position (x,y,z)
# IN: dirCam: camera direction (x,y,z)
# IN: viewAngle: camera view angle
# IN: bgColor: background color
# IN: backgroundFile: background image file
# IN: shadow: enable shadows (0/1)
# IN: lightOffset: light offset (x,y)
# IN: dof: depth of field (0/1)
# IN: dofPower: depth of field power
# IN: gamma: gamma correction
# IN: toneMapping: tone mapping mode
# IN: stereo: stereo rendering mode
# IN: stereoDist: stereo separation distance
# IN: panorama: panorama mode
# IN: export: export mode
# IN: exportResolution: export resolution
# IN: exportAA: anti-aliasing level for export
# IN: zoneNames: zone names list
# IN: renderTags: render tags list
# IN: frameBuffer: frame buffer index
# IN: offscreen: offscreen mode
# IN: posCamList: list of camera positions for animation
# IN: posEyeList: list of eye positions for animation
# IN: dirCamList: list of directions for animation
# OUT: None
def display(arrays,
            dim=-1, mode=-1,
            scalarField=-1,
            vectorField1=-1, vectorField2=-1, vectorField3=-1,
            displayBB=-1, displayInfo=-1, displayIsoLegend=-1,
            meshStyle=-1, solidStyle=-1, scalarStyle=-1,
            vectorStyle=-1, vectorScale=-1., vectorDensity=-1.,
            vectorNormalize=-1, vectorShowSurface=-1, vectorShape=-1,
            vectorProjection=-1,
            colormap=-1, colormapC1="", colormapC2="", colormapC3="",
            colormapC=None,
            niso=-1, isoEdges=-1, isoScales=[],
            win=(-1,-1),
            posCam=(-999,-999,-999),
            posEye=(-999,-999,-999),
            dirCam=(-999,-999,-999),
            viewAngle=-1.,
            bgColor=-1, backgroundFile="None",
            shadow=-1, lightOffset=(-999,-999),
            dof=-1, dofPower=-1, gamma=-1, toneMapping=-1,
            stereo=-1, stereoDist=-1., panorama=0,
            export="None", exportResolution="None", exportAA=-1,
            zoneNames=[], renderTags=[], frameBuffer=-1,
            offscreen=0,
            posCamList=None, posEyeList=None, dirCamList=None):
    """Display arrays.
    Usage: display(arrays)"""
    if arrays != [] and not isinstance(arrays[0], list): arrays = [arrays]
    global __slot__
    getModule(offscreen)
    if offscreen == 1 or offscreen == 5 or offscreen == 6 or offscreen == 7:
        if __slot__ is None:
            shaderPath = os.path.dirname(__file__)+'/OSMESA/'
            cplotm.setShaderPath(shaderPath)
            cplotm.displayNew(arrays, dim, mode, scalarField, vectorField1,
                              vectorField2, vectorField3, displayBB, displayInfo,
                              displayIsoLegend, meshStyle, solidStyle,
                              scalarStyle, vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                              vectorShowSurface, vectorShape, vectorProjection,
                              colormap, colormapC1, colormapC2, colormapC3, colormapC,
                              niso, isoEdges, isoScales, win,
                              posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                              shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                              stereo, stereoDist, panorama,
                              export, exportResolution, exportAA,
                              zoneNames, renderTags, frameBuffer, offscreen,
                              posCamList, posEyeList, dirCamList)
            __slot__ = 1
        else:
            cplotm.displayAgain(arrays, dim, mode, scalarField, vectorField1,
                                vectorField2, vectorField3, displayBB, displayInfo,
                                displayIsoLegend, meshStyle, solidStyle,
                                scalarStyle, vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                                vectorShowSurface, vectorShape, vectorProjection,
                                colormap, colormapC1, colormapC2, colormapC3, colormapC,
                                niso, isoEdges, isoScales, win,
                                posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                                shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                                stereo, stereoDist, panorama,
                                export, exportResolution, exportAA,
                                zoneNames, renderTags, frameBuffer, offscreen,
                                posCamList, posEyeList, dirCamList)
        return
    if __slot__ is None:
        shaderPath = os.path.dirname(__file__)+'/'
        cplotm.setShaderPath(shaderPath)
        displayNew__(arrays, dim, mode, scalarField, vectorField1,
                     vectorField2, vectorField3, displayBB, displayInfo,
                     displayIsoLegend, meshStyle, solidStyle,
                     scalarStyle, vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                     vectorShowSurface, vectorShape, vectorProjection,
                     colormap, colormapC1, colormapC2, colormapC3, colormapC,
                     niso, isoEdges, isoScales, win,
                     posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                     shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                     stereo, stereoDist, panorama,
                     export, exportResolution, exportAA,
                     zoneNames, renderTags, frameBuffer, offscreen,
                     posCamList, posEyeList, dirCamList)
        __slot__ = 1
    else:
        displayAgain__(arrays, dim, mode, scalarField, vectorField1,
                       vectorField2, vectorField3, displayBB, displayInfo,
                       displayIsoLegend, meshStyle, solidStyle,
                       scalarStyle, vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                       vectorShowSurface, vectorShape, vectorProjection,
                       colormap, colormapC1, colormapC2, colormapC3, colormapC,
                       niso, isoEdges, isoScales, win,
                       posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                       shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                       stereo, stereoDist, panorama,
                       export, exportResolution, exportAA,
                       zoneNames, renderTags, frameBuffer, offscreen,
                       posCamList, posEyeList, dirCamList)

# Force rendering of the current scene
# IN: None
# OUT: None
def render():
    """Force render.
    Usage: render()"""
    cplotm.render()

# Delete zones from the plotter
# IN: zlist: list of zone indices to delete
# OUT: None
def delete(zlist):
    """Delete zones from plotter.
    Usage: delete([i1,i2,...])"""
    if __slot__ is None: return
    cplotm.delete(zlist)

#==============================================================================
# Add one array/zone to plotter
# IN: arrays: list of arrays
# IN: no: index position to insert (-1=append)
# IN: array: array to add
# IN: zoneName: name of the zone
# IN: renderTag: render tag associated with zone
# OUT: None
def add(arrays, no, array, zoneName=None, renderTag=None):
    """Add one zone to plotter.
    Usage: add(arrays, no, array, zoneName, renderTag)"""
    nzs = 0; nzu = 0
    for i in range(no):
        if len(arrays[i]) == 5: nzs += 1
        else: nzu += 1
    arrays.insert(no, array)
    if __slot__ is None:
        if no == -1: no = len(arrays)
        arrays.insert(no, array)
        if len(arrays) == 1:
            display(arrays, zoneNames=[zoneName], renderTags=[renderTag])
        else: display(arrays)
    else:
        cplotm.add(array, (nzs, nzu), zoneName, renderTag)

#==============================================================================
# Replace an array/zone in plotter
# IN: arrays: list of arrays
# IN: no: index of array to replace
# IN: array: new array to replace with
# IN: zoneName: name of the zone
# IN: renderTag: render tag associated with zone
# OUT: None
def replace(arrays, no, array, zoneName=None, renderTag=None):
    """Replace arrays[no] by array.
    Usage: replace(arrays, no, array, zoneName, renderTag)"""
    zone = arrays[no]
    if len(zone) == 5: oldType = 1
    else: oldType = 2
    nzs = 0; nzu = 0
    for i in range(no):
        if len(arrays[i]) == 5: nzs += 1
        else: nzu += 1
    arrays[no] = array
    if __slot__ is None:
        if len(arrays) == 1:
            display(arrays, zoneNames=[zoneName], renderTags=[renderTag])
        else: display(arrays)
    else:
        cplotm.replace(array, (nzs, nzu, oldType), zoneName, renderTag)

#==============================================================================
# Display 1D curves/plots
# IN: arrays: list of arrays or numpy arrays [x,y]
# IN: slot: slot index for multi-grid display
# IN: gridPos: grid position (row, col)
# IN: gridSize: grid size (rows, cols) (-1,-1=auto)
# IN: bgBlend: background blend factor
# IN: var1: first variable name (x-axis)
# IN: var2: second variable name (y-axis)
# IN: r1: range for var1 (min, max)
# IN: r2: range for var2 (min, max)
# IN: offscreen: offscreen mode
# OUT: None
def display1D(arrays, slot=0, gridPos=(0,0), gridSize=(-1,-1),
              bgBlend=1., var1='x', var2='y', r1=None, r2=None, offscreen=0):
    """Display 1D plots.
    Usage: display1D(arrays, slot, ....)"""
    import Converter
    import numpy
    if len(arrays) > 0:
        if isinstance(arrays[0], numpy.ndarray): # numpy [x,y]
            if len(arrays) < 2:
                raise ValueError('display1D: requires at least two numpys [x,y]')
            x = arrays[0]; y = arrays[1]
            n = x.size
            array = Converter.array(var1+','+var2,n,1,1)
            array[1][0,:] = x[:]
            array[1][1,:] = y[:]
            arrays = [array]
        elif not isinstance(arrays[0], list): arrays = [arrays]
    if __slot__ is None: display([])
    try:
        arrays = Converter.convertArray2Hexa(arrays)
        minv1 = Converter.getMinValue(arrays, var1)
        maxv1 = Converter.getMaxValue(arrays, var1)
        minv2 = Converter.getMinValue(arrays, var2)
        maxv2 = Converter.getMaxValue(arrays, var2)
        if r1 is None: r1 = (minv1, maxv1)
        if r2 is None: r2 = (minv2, maxv2)
    except: pass
    if r1 is None: r1 = (0.,1.)
    if r2 is None: r2 = (0.,1.)
    getModule(offscreen)
    cplotm.display1D(arrays, slot, gridPos, gridSize,
                     bgBlend, var1, var2, r1, r2)
    time.sleep(__timeStep__)

#==============================================================================
# Wait for user to press a key in plotter window
# IN: None
# OUT: None
def pressKey():
    """Wait for user to press a key.
    Usage: pressKey()"""
    cplotm.pressKey()

#==============================================================================
# -- get/set functions --
#==============================================================================

# Get plotter state
# IN: mode: state mode to retrieve (e.g., "posCam", "colormap", etc.)
# OUT: state value
def getState(mode):
    """Return a state in plotter.
    Usage: n = getState(mode)"""
    return cplotm.getState(mode)

# Get selected zone
# IN: None
# OUT: currently selected zone index
def getSelectedZone():
    """Return the selected zone in plotter.
    Usage: n = getSelectedZone()"""
    return cplotm.getSelectedZone()

# Get selected zones
# IN: None
# OUT: list of currently selected zone indices
def getSelectedZones():
    """Return the selected zones in plotter.
    Usage: list = getSelectedZones()"""
    return cplotm.getSelectedZones()

# Get selected status
# IN: zone: zone index
# OUT: selection status (0=not selected, 1=selected)
def getSelectedStatus(zone):
    """Return the selected status of a zone in plotter.
    Usage: status = getSelectedStatus(zone)"""
    return cplotm.getSelectedStatus(zone)

# Get active zones
# IN: None
# OUT: list of active (displayed) zone indices
def getActiveZones():
    """Return the active (displayed) zones in plotter.
    Usage: list = getActiveZones()"""
    return cplotm.getActiveZones()

# Get active status
# IN: zone: zone index
# OUT: active status (0=not active, 1=active)
def getActiveStatus(zone):
    """Return the active status of a zone in plotter.
    Usage: status = getActiveStatus(zone)"""
    return cplotm.getActiveStatus(zone)

# Get active point
# IN: None
# OUT: active (clicked) point coordinates (x,y,z)
def getActivePoint():
    """Return the active (clicked) point in plotter.
    Usage: n = getActivePoint()"""
    return cplotm.getActivePoint()

# Get active point index
# IN: None
# OUT: active point index in mesh
def getActivePointIndex():
    """Return the active (clicked) point index.
    Usage: n = getActivePointIndex()"""
    return cplotm.getActivePointIndex()

# Get active point field values
# IN: None
# OUT: field values at active point
def getActivePointF():
    """Return the active (clicked) point field values.
    Usage: f = getActivePointF()"""
    return cplotm.getActivePointF()

# Get mouse state
# IN: None
# OUT: mouse state (position and button state)
def getMouseState():
    """Return mouse state (mouse position and button state)."""
    return cplotm.getMouseState()

# Get keyboard state
# IN: None
# OUT: pressed keys string
def getKeyboard():
    """Return the pressed keys.
    Usage: n = getKeyboard()"""
    return cplotm.getKeyboard()

# Reset keyboard state
# IN: None
# OUT: None
def resetKeyboard():
    """Reset the keyboard string.
    Usage: resetKeyboard()"""
    return cplotm.resetKeyboard()

# Convert r,g,b tuple in [0,1] to hexa color string
# IN: T: tuple of (r, g, b) values in range [0,1]
# OUT: hexa color string (e.g., "#FF0000")
def convertRGB2String(T):
    """Convert r,g,b tuple to hexa color string."""
    r,g,b = T
    a1 = int(r*255.)
    a2 = int(g*255.)
    a3 = int(b*255.)
    return f"#{a1:02X}{a2:02X}{a3:02X}"

# Filter and convert colormap values to internal format
# IN: values: list of [colormap, colormapC1, colormapC2, colormapC3, colormapC]
# OUT: list of filtered colormap values in internal format
def filterColormap(values):
    [colormap, colormapC1, colormapC2, colormapC3, colormapC] = values
    shift = colormap % 2
    if colormap == 16 or colormap == 17: # Viridis
        colormapC = ColorMaps.cmapDict['Viridis']
        colormap = 10+shift
    elif colormap == 18 or colormap == 19: # Inferno
        colormapC = ColorMaps.cmapDict['Inferno']
        colormap = 10+shift
    elif colormap == 20 or colormap == 21: # Magma
        colormapC = ColorMaps.cmapDict['Magma']
        colormap = 10+shift
    elif colormap == 22 or colormap == 23: # Plasma
        colormapC = ColorMaps.cmapDict['Plasma']
        colormap = 10+shift
    elif colormap == 24 or colormap == 25: # Jet
        colormapC = ColorMaps.cmapDict['Jet']
        colormap = 10+shift
    elif colormap == 26 or colormap == 27: # Greys
        colormapC = ColorMaps.cmapDict['Greys']
        colormap = 10+shift
    elif colormap == 28 or colormap == 29: # NiceBlue
        colormapC1='#000000'; colormapC2='#FFFFFF'; colormapC3='#0061A5'
        colormap = 6+shift
    elif colormap == 30 or colormap == 31: # Greens
        colormapC = ColorMaps.cmapDict['Greens']
        colormap = 10+shift

    # if colormapC are given in rgb, convert to hexa string
    if not isinstance(colormapC1, str):
        colormapC1 = convertRGB2String(colormapC1)
    if not isinstance(colormapC2, str):
        colormapC2 = convertRGB2String(colormapC2)
    if not isinstance(colormapC1, str):
        colormapC3 = convertRGB2String(colormapC3)

    return [colormap, colormapC1, colormapC2, colormapC3, colormapC]

# Get filtered colormap index from current state
# Detects the actual colormap used when getState returns 10 or 11
# IN: None
# OUT: colormap index (16-31 for named colormaps)
def getFilteredColormap():
    colormap = getState("colormap")
    if colormap != 10 and colormap != 11: return colormap
    if colormap == 10: shift = 0
    else: shift = 1
    colormapC = getState("colormapC")
    CC = numpy.array([colormapC[i] for i in range(0,5)])
    C1 = numpy.array([ColorMaps.Jet[i] for i in range(0,5)])
    if numpy.allclose(CC, C1): return 24+shift
    C1 = numpy.array([ColorMaps.Magma[i] for i in range(0,5)])
    if numpy.allclose(CC, C1): return 20+shift
    C1 = numpy.array([ColorMaps.Inferno[i] for i in range(0,5)])
    if numpy.allclose(CC, C1): return 18+shift
    C1 = numpy.array([ColorMaps.Viridis[i] for i in range(0,5)])
    if numpy.allclose(CC, C1): return 16+shift
    C1 = numpy.array([ColorMaps.Plasma[i] for i in range(0,5)])
    if numpy.allclose(CC, C1): return 22+shift
    C1 = numpy.array([ColorMaps.Greys[i] for i in range(0,5)])
    if numpy.allclose(CC, C1): return 26+shift
    C1 = numpy.array([ColorMaps.Greens[i] for i in range(0,5)])
    if numpy.allclose(CC, C1): return 30+shift
    return colormap

# Set CPlot state parameters
# IN: dim: dimension (3=3D, 2=2D, 1=1D)
# IN: mode: display mode
# IN: scalarField: scalar field index
# IN: vectorField1/2/3: vector field components
# IN: displayBB: display bounding box (0/1)
# IN: displayInfo: display info overlay (0/1)
# IN: displayIsoLegend: display isosurface legend (0/1)
# IN: meshStyle: mesh display style
# IN: solidStyle: solid display style
# IN: scalarStyle: scalar display style
# IN: vectorStyle: vector display style
# IN: vectorScale: vector scaling factor
# IN: vectorDensity: vector density
# IN: vectorNormalize: normalize vectors (0/1)
# IN: vectorShowSurface: show surface with vectors (0/1)
# IN: vectorShape: vector shape type
# IN: vectorProjection: vector projection mode
# IN: colormap: colormap index
# IN: colormapC1/C2/C3: colormap color stops
# IN: colormapC: colormap color array
# IN: niso: number of isosurfaces
# IN: isoEdges: display isosurface edges (0/1)
# IN: isoScales: isosurface scale factors
# IN: win: window position (x,y)
# IN: posCam: camera position (x,y,z)
# IN: posEye: camera eye position (x,y,z)
# IN: dirCam: camera direction (x,y,z)
# IN: viewAngle: camera view angle
# IN: lightOffset: light offset (x,y)
# IN: bgColor: background color
# IN: backgroundFile: background image file
# IN: shadow: enable shadows (0/1)
# IN: dof: depth of field (0/1)
# IN: dofPower: depth of field power
# IN: gamma: gamma correction
# IN: toneMapping: tone mapping mode
# IN: sobelThreshold: sobel threshold for edge detection
# IN: sharpenPower: sharpening power
# IN: ssaoPower: screen space ambient occlusion power
# IN: ghostifyDeactivatedZones: ghostify deactivated zones (0/1)
# IN: edgifyActivatedZones: edgify activated zones (0/1)
# IN: edgifyDeactivatedZones: edgify deactivated zones (0/1)
# IN: simplifyOnDrag: simplify rendering on drag (0/1)
# IN: export: export mode
# IN: exportResolution: export resolution
# IN: exportAA: anti-aliasing level for export
# IN: continuousExport: continuous export mode (0/1)
# IN: envmap: environment map file
# IN: message: message to display
# IN: stereo: stereo rendering mode
# IN: stereoDist: stereo separation distance
# IN: cursor: cursor type
# IN: gridSize: grid size (rows, cols)
# IN: timer: timer value
# IN: selectionStyle: selection style
# IN: activateShortCuts: activate keyboard shortcuts (0/1)
# IN: billBoards: billboard settings
# IN: billBoardSize: billboard size
# IN: materials: material settings
# IN: bumpMaps: bump map settings
# IN: frameBuffer: frame buffer index
# IN: offscreen: offscreen mode
# OUT: None
def setState(dim=-1,
             mode=-1,
             scalarField=-1,
             vectorField1=-1, vectorField2=-1, vectorField3=-1,
             displayBB=-1,
             displayInfo=-1,
             displayIsoLegend=-1,
             meshStyle=-1,
             solidStyle=-1,
             scalarStyle=-1,
             vectorStyle=-1,
             vectorScale=-1.,
             vectorDensity=-1.,
             vectorNormalize=-1,
             vectorShowSurface=-1,
             vectorShape=-1,
             vectorProjection=-1,
             colormap=-1,
             colormapC1="", colormapC2="", colormapC3="", colormapC=None,
             niso=-1,
             isoEdges=-1,
             isoScales=[],
             win=(-1,-1),
             posCam=(-999,-999,-999),
             posEye=(-999,-999,-999),
             dirCam=(-999,-999,-999),
             viewAngle=-1,
             lightOffset=(-999,-999),
             bgColor=-1,
             backgroundFile="None",
             shadow=-1,
             dof=-1, dofPower=-1,
             gamma=-1,
             toneMapping=-1,
             sobelThreshold=-1,
             sharpenPower=-1,
             ssaoPower=-1,
             ghostifyDeactivatedZones=-1,
             edgifyActivatedZones=-1,
             edgifyDeactivatedZones=-1,
             simplifyOnDrag=-1,
             export="None",
             exportResolution="None",
             exportAA=-1,
             continuousExport=-1,
             envmap="None", message="None",
             stereo=-1, stereoDist=-1.,
             cursor=-1,
             gridSize=(-1, -1),
             timer=-1,
             selectionStyle=-1,
             activateShortCuts=-1,
             billBoards=None,
             billBoardSize=-1,
             materials=None, bumpMaps=None,
             frameBuffer=-1,
             offscreen=0):
    """Set CPlot state.
    Usage: setState(posCam=(12,0,0))"""
    getModule(offscreen)
    if colormap != -1:
        [colormap, colormapC1, colormapC2, colormapC3, colormapC] = filterColormap( [colormap, colormapC1, colormapC2, colormapC3, colormapC] )
    if offscreen == 1 or offscreen == 5 or offscreen == 6 or offscreen == 7: # must set in the right module
        from . import cplotOSMesa as cplot
    else:
        from . import cplot
    cplotm.setState(dim, mode, scalarField, vectorField1, vectorField2,
                    vectorField3, displayBB, displayInfo, displayIsoLegend,
                    meshStyle, solidStyle, scalarStyle,
                    vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                    vectorShowSurface, vectorShape, vectorProjection,
                    colormap, colormapC1, colormapC2, colormapC3, colormapC,
                    niso, isoEdges, isoScales, win,
                    posCam, posEye, dirCam, viewAngle, lightOffset,
                    bgColor, backgroundFile, shadow,
                    dof, dofPower, gamma, toneMapping,
                    sobelThreshold, sharpenPower, ssaoPower,
                    ghostifyDeactivatedZones, edgifyActivatedZones,
                    edgifyDeactivatedZones, simplifyOnDrag,
                    export, exportResolution, exportAA, continuousExport,
                    envmap, message, stereo, stereoDist,
                    cursor, gridSize, timer, selectionStyle,
                    activateShortCuts, billBoards, billBoardSize,
                    materials, bumpMaps, frameBuffer, offscreen)

# Set CPlot display mode
# IN: mode: display mode (0=solid, 1=wireframe, etc.)
# OUT: None
def setMode(mode):
    """Set CPlot display mode.
    Usage: setMode(0)"""
    cplotm.setMode(mode)

# Change displayed variable interactively
# IN: None
# OUT: None
def changeVariable():
    """Change displayed variable.
    Usage: changeVariable()"""
    cplotm.changeVariable()

# Change CPlot display style interactively
# IN: None
# OUT: None
def changeStyle():
    """Change CPlot display style.
    Usage: changeStyle()"""
    cplotm.changeStyle()

# Change CPlot info display style interactively
# IN: None
# OUT: None
def changeInfoDisplay():
    """Change CPlot info display style.
    Usage: changeInfoDisplay()"""
    cplotm.changeInfoDisplay()

# Change the blanking procedure interactively
# IN: None
# OUT: None
def changeBlanking():
    """Change the blanking procedure.
    Usage: changeBlanking()"""
    cplotm.changeBlanking()

# Set CPlot display dimension
# IN: dim: dimension (3=3D, 2=2D, 1=1D)
# OUT: None
def setDim(dim):
    """Set CPlot display dim 3, 2 or 1.
    Usage: setDim(2)"""
    cplotm.setDim(dim)

# Set the active (clicked) point in plotter
# IN: x: x coordinate
# IN: y: y coordinate
# IN: z: z coordinate
# OUT: None
def setActivePoint(x,y,z):
    """Set the active (clicked) point in plotter.
    Usage: setActivePoint(x,y,z)"""
    return cplotm.setActivePoint(x,y,z)

# Set selected zones
# IN: zlist: list of zone indices to select [(slot, zone), ...]
# OUT: None
def setSelectedZones(zlist):
    """Set selected zones.
    Usage: setSelectedZones([(0,1),(1,1),...])"""
    cplotm.setSelectedZones(zlist)

# Unselect all zones
# IN: None
# OUT: None
def unselectAllZones():
    """Unselect all zones.
    Usage: unselectAllZones()"""
    cplotm.unselectAllZones()

# Set active (displayed) zones
# IN: zlist: list of zone indices to activate [(slot, zone), ...]
# OUT: None
def setActiveZones(zlist):
    """Set active (displayed) zones.
    Usage: setActiveZones([(0,1)])"""
    cplotm.setActiveZones(zlist)

# Set zone names
# IN: zlist: list of (zone_index, name) tuples
# OUT: None
def setZoneNames(zlist):
    """Set zone names.
    Usage: setZoneNames([(0,'myZone')])"""
    cplotm.setZoneNames(zlist)

#==============================================================================
# Look for (focus on) selected zones
# IN: None
# OUT: None
def lookFor():
    """Look for selected zones.
    Usage: lookFor()"""
    cplotm.lookFor()

# Fit the view to objects in scene
# IN: None
# OUT: None
def fitView():
    """Fit the view to objects.
    Usage: fitView()"""
    cplotm.fitView()

# Finalize export operation
# IN: action: export action type (0=default, 1/5/6/7=special modes)
# OUT: None
def finalizeExport(action=0):
    """Finalize export."""
    if action == 1 or action == 5 or action == 6 or action == 7:
        cplotm.finalizeExport(action)
        return
    while cplotm.isDisplayRunning() == 0: pass
    cplotm.finalizeExport(action)

# Hide the CPlot window
# IN: None
# OUT: None
def hide():
    """Hide window."""
    cplotm.hide()

# Show the CPlot window if hidden
# IN: None
# OUT: None
def show():
    """Show window if it has been hidden with flush."""
    cplotm.show()

# Set the CPlot window title
# IN: file: file name to display
# IN: path: file path to display
# OUT: None
def setWindowTitle(file, path):
    """Set the CPlot window title."""
    getModule()
    cplotm.setWindowTitle(file, path)

#==============================================================================
# camera
#==============================================================================
# Move camera following check points (animation path)
# IN: posCams: list of camera positions for path
# IN: posEyes: list of eye positions for path (optional)
# IN: dirCams: list of camera directions for path (optional)
# IN: moveEye: move eye position along path (0/1)
# IN: N: number of interpolation points
# IN: speed: animation speed factor
# IN: pos: start position index (-1=start from beginning)
# OUT: tuple of (posCam, posEye, dirCam) final camera state
def moveCamera(posCams, posEyes=None, dirCams=None, moveEye=False, N=100, speed=1., pos=-1):
    """Move posCam and posEye following check points."""
    # Set d, array of posCams and N nbre of points
    import Geom
    if len(posCams) == 5 and isinstance(posCams[0], str): # input struct array
        N = posCams[2]
        d = posCams
        pinc = KCore.isNamePresent(posCams, 'incEye')
        if pinc >= 0: pinc = posCams[1][pinc]
        else: pinc = None
    else: # list
        N = max(N, 3)
        Np = len(posCams)
        pOut = []
        P0 = posCams[0]
        pOut.append(P0)
        for i in range(1,Np):
            P1 = posCams[i]
            sub = Vector.sub(P1,P0)
            if Vector.norm(sub)>1.e-10: pOut.append(P1)
            P0 = P1
        if len(pOut) == 1:
            d = Geom.polyline(pOut*N)
        elif len(pOut) == 2:
            p = Geom.polyline(pOut)
            d = Geom.spline(p, 2, N)
        else:
            p = Geom.polyline(pOut)
            d = Geom.spline(p, 3, N)
        pinc = None

    # Set e, array of posEye of N pts
    if posEyes is not None:
        if len(posEyes) == 5 and isinstance(posEyes[0], str): # input struct array
            Neye = posEyes[2]
            if Neye != N:
                import Generator
                dis = Geom.getDistribution(d)
                posEyes = Generator.map(posEyes, dis, 1)
            e = posEyes
        else: # list
            Np = len(posEyes)
            pOut = []
            P0 = posEyes[0]
            pOut.append(P0)
            for i in range(1,Np):
                P1 = posEyes[i]
                sub = Vector.sub(P1,P0)
                if Vector.norm(sub)>1.e-10: pOut.append(P1)
                P0 = P1
            if len(pOut) == 1:
                e = Geom.polyline(pOut*N)
            elif len(pOut) == 2:
                p = Geom.polyline(pOut)
                e = Geom.spline(p, 2, N)
            else:
                p = Geom.polyline(pOut)
                e = Geom.spline(p, 3, N)
    else: e = None

    # Set dc, array of dirCams of N pts
    if dirCams is not None:
        if len(dirCams) == 5 and isinstance(dirCams[0], str): # input struct array
            Ndc = dirCams[2]
            if Ndc != N:
                import Generator
                dis = Geom.getDistribution(d)
                dirCams = Generator.map(dirCams, dis, 1)
            dc = dirCams
        else: # list
            Np = len(dirCams)
            pOut = []
            P0 = dirCams[0]
            pOut.append(P0)
            for i in range(1,Np):
                P1 = dirCams[i]
                sub = Vector.sub(P1,P0)
                if Vector.norm(sub)>1.e-10: pOut.append(P1)
                P0 = P1
            if len(pOut) == 1:
                dc = Geom.polyline(pOut*N)
            elif len(pOut) == 2:
                p = Geom.polyline(pOut)
                dc = Geom.spline(p, 2, N)
            else:
                p = Geom.polyline(pOut)
                dc = Geom.spline(p, 3, N)
    else: dc = None

    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')

    if pos == -1:
        i = 0
        while i < N-1:
            time.sleep(__timeStep__*speed*0.06)
            if i > N-11: inc = N-i-1
            else: inc = 10
            posCam = (d[1][0,i],d[1][1,i],d[1][2,i])
            if e is not None:
                posEye = (e[1][0,i],e[1][1,i],e[1][2,i])
                if dc is not None:
                    dirCam = (dc[1][0,i],dc[1][1,i],dc[1][2,i])
                    setState(posCam=posCam, posEye=posEye, dirCam=dirCam)
                else:
                    setState(posCam=posCam, posEye=posEye)
            elif moveEye:
                posEye = (d[1][0,i+inc],d[1][1,i+inc],d[1][2,i+inc])
                setState(posCam=posCam, posEye=posEye)
            else: setState(posCam=posCam)
            i += 1
    else:
        i = pos; i = min(pos, N-1)
        if pinc is not None: inc = int(pinc[i])
        else: inc = 10
        inc = min(inc, N-i-1)
        posCam = (d[1][0,i],d[1][1,i],d[1][2,i])
        if e is not None:
            posEye = (e[1][0,i],e[1][1,i],e[1][2,i])
            if dc is not None:
                dirCam = (dc[1][0,i],dc[1][1,i],dc[1][2,i])
                setState(posCam=posCam, posEye=posEye, dirCam=dirCam)
            else:
                setState(posCam=posCam, posEye=posEye)
        elif moveEye:
            posEye = (d[1][0,i+inc],d[1][1,i+inc],d[1][2,i+inc])
            setState(posCam=posCam, posEye=posEye)
        else: setState(posCam=posCam)

    return posCam, posEye, dirCam

# Travel camera right (orbit around target)
# IN: xr: travel distance as fraction of camera distance
# IN: N: number of interpolation points
# OUT: tuple of (posCam, posEye, dirCam) final camera state
def travelRight(xr=0.1, N=100):
    """Travel camera right."""
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = 2.*3.14*R*xr*0.5
    d2 = Vector.cross(d1, dirCam)
    d2 = Vector.normalize(d2)
    d2 = Vector.mul(L,d2)
    P2 = Vector.sub(posCam, d2)
    d3 = Vector.sub(posEye, P2)
    d4 = Vector.cross(d3, dirCam)
    d4 = Vector.normalize(d4)
    d4 = Vector.mul(L,d4)
    P3 = Vector.sub(P2, d4)
    checkPoints = [posCam,tuple(P2),tuple(P3)]
    return moveCamera(checkPoints, N=N)

# Travel camera left (orbit around target)
# IN: xr: travel distance as fraction of camera distance
# IN: N: number of interpolation points
# OUT: tuple of (posCam, posEye, dirCam) final camera state
def travelLeft(xr=0.1, N=100):
    """Travel camera left."""
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = 2*3.14*R*xr*0.5
    d2 = Vector.cross(d1, dirCam)
    d2 = Vector.normalize(d2)
    d2 = Vector.mul(L,d2)
    P2 = Vector.add(posCam, d2)
    d3 = Vector.sub(posEye, P2)
    d4 = Vector.cross(d3, dirCam)
    d4 = Vector.normalize(d4)
    d4 = Vector.mul(L,d4)
    P3 = Vector.add(P2, d4)
    checkPoints = [posCam,tuple(P2),tuple(P3)]
    return moveCamera(checkPoints, N=N)

# Travel camera up (move along view direction)
# IN: xr: travel distance as fraction of camera distance
# IN: N: number of interpolation points
# OUT: tuple of (posCam, posEye, dirCam) final camera state
def travelUp(xr=0.1, N=100):
    """Travel camera up."""
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = 2*3.14*R*xr*1.
    d2 = Vector.normalize(dirCam)
    d2 = Vector.mul(L,d2)
    P2 = Vector.add(posCam, d2)
    checkPoints = [posCam,tuple(P2)]
    return moveCamera(checkPoints, N=N)

# Travel camera down (move opposite to view direction)
# IN: xr: travel distance as fraction of camera distance
# IN: N: number of interpolation points
# OUT: tuple of (posCam, posEye, dirCam) final camera state
def travelDown(xr=0.1, N=100):
    """Travel camera down."""
    posCam = getState('posCam')
    posEye = getState('posEye')
    dirCam = getState('dirCam')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = 2*3.14*R*xr*1.
    d2 = Vector.normalize(dirCam)
    d2 = Vector.mul(L,d2)
    P2 = Vector.sub(posCam, d2)
    checkPoints = [posCam,tuple(P2)]
    return moveCamera(checkPoints, N=N)

# Zoom camera in (move closer to target)
# IN: xr: zoom factor as fraction of camera distance
# IN: N: number of interpolation points
# OUT: tuple of (posCam, posEye, dirCam) final camera state
def travelIn(xr=0.1, N=100):
    """Zoom camera in."""
    posCam = getState('posCam')
    posEye = getState('posEye')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = R*xr
    d2 = Vector.mul(L,d1)
    P2 = Vector.add(posCam, d2)
    checkPoints = [posCam,tuple(P2)]
    return moveCamera(checkPoints, N=N)

# Zoom camera out (move away from target)
# IN: xr: zoom factor as fraction of camera distance
# IN: N: number of interpolation points
# OUT: tuple of (posCam, posEye, dirCam) final camera state
def travelOut(xr=0.1, N=100):
    """Zoom camera out."""
    posCam = getState('posCam')
    posEye = getState('posEye')
    d1 = Vector.sub(posEye,posCam)
    R = Vector.norm(d1)
    L = R*xr
    d2 = Vector.mul(L,d1)
    P2 = Vector.sub(posCam, d2)
    checkPoints = [posCam,tuple(P2)]
    return moveCamera(checkPoints, N=N)

#==============================================================================
# image as array
#==============================================================================

# Blur an image array
# IN: a: image array or list of arrays to blur
# IN: blurSigma: blur sigma value (default=0.8)
# OUT: None (modifies array in place)
def blur(a, blurSigma=0.8):
    """Blur an image array."""
    getModule()
    if isinstance(a[0], list):
        for i in a: cplotm.blur(i, blurSigma)
    else: cplotm.blur(a, blurSigma)
    return None

#==============================================================================
# -- Internal functions --
#==============================================================================

# Set internal file name for CPlot
# IN: name: file name to set
# OUT: None
def setFileName__(name):
    cplotm.setFileName(name)

# Internal function: Create new display window (threaded)
def displayNew__(arrays, dim, mode, scalarField, vectorField1, vectorField2,
                 vectorField3, displayBB, displayInfo, displayIsoLegend,
                 meshStyle, solidStyle, scalarStyle, vectorStyle,
                 vectorScale, vectorDensity, vectorNormalize, vectorShowSurface,
                 vectorShape, vectorProjection,
                 colormap, colormapC1, colormapC2, colormapC3, colormapC,
                 niso, isoEdges, isoScales, win,
                 posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                 shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                 stereo, stereoDist, panorama,
                 export, exportResolution, exportAA,
                 zoneNames, renderTags, frameBuffer, offscreen,
                 posCamList, posEyeList, dirCamList):
    global __slot__
    import threading
    if colormap != -1:
        [colormap, colormapC1, colormapC2, colormapC3, colormapC] = filterColormap( [colormap, colormapC1, colormapC2, colormapC3, colormapC] )
    if offscreen > 0: daemon = True
    else: daemon = False
    if version_info[0] == 2: # python 2 - no daemon
        a = threading.Thread(None, cplotm.displayNew, None,
                             (arrays, dim, mode, scalarField, vectorField1,
                              vectorField2, vectorField3, displayBB, displayInfo,
                              displayIsoLegend,
                              meshStyle, solidStyle, scalarStyle,
                              vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                              vectorShowSurface, vectorShape, vectorProjection,
                              colormap, colormapC1, colormapC2, colormapC3, colormapC,
                              niso, isoEdges, isoScales,
                              win, posCam, posEye, dirCam, viewAngle,
                              bgColor, backgroundFile,
                              shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                              stereo, stereoDist, panorama,
                              export, exportResolution, exportAA,
                              zoneNames, renderTags, frameBuffer, offscreen,
                              posCamList, posEyeList, dirCamList), {})
    else: # python3 - daemon
        a = threading.Thread(None, cplotm.displayNew, None,
                             (arrays, dim, mode, scalarField, vectorField1,
                              vectorField2, vectorField3, displayBB, displayInfo,
                              displayIsoLegend,
                              meshStyle, solidStyle, scalarStyle,
                              vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                              vectorShowSurface, vectorShape, vectorProjection,
                              colormap, colormapC1, colormapC2, colormapC3, colormapC,
                              niso, isoEdges, isoScales,
                              win, posCam, posEye, dirCam, viewAngle,
                              bgColor, backgroundFile,
                              shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                              stereo, stereoDist, panorama,
                              export, exportResolution, exportAA,
                              zoneNames, renderTags, frameBuffer, offscreen,
                              posCamList, posEyeList, dirCamList), {}, daemon=daemon)
    a.start()
    __slot__ = a

# Internal function: Update existing display (threaded)
def displayAgain__(arrays, dim, mode, scalarField, vectorField1, vectorField2,
                   vectorField3, displayBB, displayInfo, displayIsoLegend,
                   meshStyle, solidStyle, scalarStyle, vectorStyle,
                   vectorScale, vectorDensity, vectorNormalize, vectorShowSurface,
                   vectorShape, vectorProjection,
                   colormap, colormapC1, colormapC2, colormapC3, colormapC,
                   niso, isoEdges, isoScales,
                   win, posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                   shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                   stereo, stereoDist, panorama,
                   export, exportResolution, exportAA,
                   zoneNames, renderTags, frameBuffer, offscreen,
                   posCamList, posEyeList, dirCamList):
    if colormap != -1:
        [colormap, colormapC1, colormapC2, colormapC3, colormapC] = filterColormap( [colormap, colormapC1, colormapC2, colormapC3, colormapC] )
    cplotm.displayAgain(arrays, dim, mode, scalarField, vectorField1,
                        vectorField2, vectorField3, displayBB, displayInfo,
                        displayIsoLegend,
                        meshStyle, solidStyle, scalarStyle, vectorStyle,
                        vectorScale, vectorDensity, vectorNormalize, vectorShowSurface,
                        vectorShape, vectorProjection,
                        colormap, colormapC1, colormapC2, colormapC3, colormapC,
                        niso, isoEdges, isoScales,
                        win, posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                        shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                        stereo, stereoDist, panorama,
                        export, exportResolution, exportAA,
                        zoneNames, renderTags, frameBuffer, offscreen,
                        posCamList, posEyeList, dirCamList)
    time.sleep(__timeStep__)
