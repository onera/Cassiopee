"""Simple plotter for CFD.
"""
__version__ = '4.0'
__author__ = "Christophe Benoit, Stephanie Peron, Pascal Raud, Matthieu Soismier, Bertrand Michel"

# cplot and cplotOSMesa must not be imported at the same time

import KCore.kcore as KCore
import KCore.Vector as Vector
from . import ColorMaps

import os.path
import time
from sys import version_info

__timeStep__ = 0.01
__slot__ = None

#==============================================================================
# -- configuration --
def configure(useRender):
    """Configure CPlot for direct rendering (cplot.useDirect), display Lists (cplot.useDL)
        or VBO (cplot.useVBO)"""
    from . import cplot
    cplot.configure(useRender)

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
            export="None", exportResolution="None",
            zoneNames=[], renderTags=[], frameBuffer=-1,
            offscreen=0,
            posCamList=None, posEyeList=None, dirCamList=None):
    """Display arrays.
    Usage: display(arrays)"""
    if arrays != [] and not isinstance(arrays[0], list): arrays = [arrays]
    if offscreen == 1 or offscreen == 5 or offscreen == 6 or offscreen == 7:
        from . import cplotOSMesa
        global __slot__
        if __slot__ is None:
            shaderPath = os.path.dirname(__file__)+'/OSMESA/'
            cplotOSMesa.setShaderPath(shaderPath)
            cplotOSMesa.displayNew(arrays, dim, mode, scalarField, vectorField1,
                                   vectorField2, vectorField3, displayBB, displayInfo,
                                   displayIsoLegend, meshStyle, solidStyle,
                                   scalarStyle, vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                                   vectorShowSurface, vectorShape, vectorProjection,
                                   colormap, colormapC1, colormapC2, colormapC3, colormapC,
                                   niso, isoEdges, isoScales, win,
                                   posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                                   shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                                   stereo, stereoDist, panorama,
                                   export, exportResolution, zoneNames, renderTags,
                                   frameBuffer, offscreen,
                                   posCamList, posEyeList, dirCamList)
            __slot__ = 1
        else:
            cplotOSMesa.displayAgain(arrays, dim, mode, scalarField, vectorField1,
                                     vectorField2, vectorField3, displayBB, displayInfo,
                                     displayIsoLegend, meshStyle, solidStyle,
                                     scalarStyle, vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                                     vectorShowSurface, vectorShape, vectorProjection,
                                     colormap, colormapC1, colormapC2, colormapC3, colormapC,
                                     niso, isoEdges, isoScales, win,
                                     posCam, posEye, dirCam, viewAngle, bgColor, backgroundFile,
                                     shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                                     stereo, stereoDist, panorama,
                                     export, exportResolution, zoneNames, renderTags,
                                     frameBuffer, offscreen,
                                     posCamList, posEyeList, dirCamList)
        return
    from . import cplot
    if __slot__ is None:
        shaderPath = os.path.dirname(__file__)+'/'
        cplot.setShaderPath(shaderPath)
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
                     export, exportResolution, zoneNames, renderTags,
                     frameBuffer, offscreen,
                     posCamList, posEyeList, dirCamList)
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
                       export, exportResolution,
                       zoneNames, renderTags, frameBuffer, offscreen,
                       posCamList, posEyeList, dirCamList)

#==============================================================================
def render():
    """Force render.
    Usage: render()"""
    from . import cplot
    cplot.render()

#==============================================================================
def delete(zlist):
    """Delete zones from plotter.
    Usage: delete([i1,i2,...])"""
    if __slot__ is None: return
    from . import cplot
    cplot.delete(zlist)

#==============================================================================
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
        from . import cplot
        cplot.add(array, (nzs, nzu), zoneName, renderTag)

#==============================================================================
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
        from . import cplot
        cplot.replace(array, (nzs, nzu, oldType), zoneName, renderTag)

#==============================================================================
def display1D(arrays, slot=0, gridPos=(0,0), gridSize=(-1,-1),
              bgBlend=1., var1='x', var2='y', r1=None, r2=None):
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
    from . import cplot
    cplot.display1D(arrays, slot, gridPos, gridSize,
                    bgBlend, var1, var2, r1, r2)
    time.sleep(__timeStep__)

#==============================================================================
def pressKey():
    """Wait for user to press a key.
    Usage: pressKey()"""
    from . import cplot
    cplot.pressKey()

#==============================================================================
# -- get/set functions --
def getState(mode):
    """Return a state in plotter.
    Usage: n = getState(mode)"""
    from . import cplot
    return cplot.getState(mode)

def getSelectedZone():
    """Return the selected zone in plotter.
    Usage: n = getSelectedZone()"""
    from . import cplot
    return cplot.getSelectedZone()

def getSelectedZones():
    """Return the selected zones in plotter.
    Usage: list = getSelectedZones()"""
    from . import cplot
    return cplot.getSelectedZones()

def getSelectedStatus(zone):
    """Return the selected status of a zone in plotter.
    Usage: status = getSelectedStatus(zone)"""
    from . import cplot
    return cplot.getSelectedStatus(zone)

def getActiveZones():
    """Return the active (displayed) zones in plotter.
    Usage: list = getActiveZones()"""
    from . import cplot
    return cplot.getActiveZones()

def getActiveStatus(zone):
    """Return the active status of a zone in plotter.
    Usage: status = getActiveStatus(zone)"""
    from . import cplot
    return cplot.getActiveStatus(zone)

def getActivePoint():
    """Return the active (clicked) point in plotter.
    Usage: n = getActivePoint()"""
    from . import cplot
    return cplot.getActivePoint()

def getActivePointIndex():
    """Return the active (clicked) point index.
    Usage: n = getActivePointIndex()"""
    from . import cplot
    return cplot.getActivePointIndex()

def getActivePointF():
    """Return the active (clicked) point field values.
    Usage: f = getActivePointF()"""
    from . import cplot
    return cplot.getActivePointF()

def getMouseState():
    """Return mouse state (mouse position and button state)."""
    from . import cplot
    return cplot.getMouseState()

def getKeyboard():
    """Return the pressed keys.
    Usage: n = getKeyboard()"""
    from . import cplot
    return cplot.getKeyboard()

def resetKeyboard():
    """Reset the keyboard string.
    Usage: resetKeyboard()"""
    from . import cplot
    return cplot.resetKeyboard()

# Ajoute des colormaps indirectes
def filterColormap(values):
    [colormap, colormapC1, colormapC2, colormapC3, colormapC] = values
    if colormap == 16 or colormap == 17: # Viridis
        colormapC = ColorMaps.Viridis
        if colormap == 16: colormap = 10
        elif colormap == 17: colormap = 11
    elif colormap == 18 or colormap == 19: # Inferno
        colormapC = ColorMaps.Inferno
        if colormap == 18: colormap = 10
        elif colormap == 19: colormap = 11
    elif colormap == 20 or colormap == 21: # Magma
        colormapC = ColorMaps.Magma
        if colormap == 20: colormap = 10
        elif colormap == 21: colormap = 11
    elif colormap == 22 or colormap == 23: # Plasma
        colormapC = ColorMaps.Plasma
        if colormap == 22: colormap = 10
        elif colormap == 23: colormap = 11
    elif colormap == 24 or colormap == 25: # Jet
        colormapC = ColorMaps.Jet2
        if colormap == 24: colormap = 10
        elif colormap == 25: colormap = 11
    elif colormap == 26 or colormap == 27: # Greys
        colormapC = ColorMaps.Greys
        if colormap == 26: colormap = 10
        elif colormap == 27: colormap = 11
    elif colormap == 28 or colormap == 29: # NiceBlue
        colormapC1='#000000'; colormapC2='#FFFFFF'; colormapC3='#0061A5'
        if colormap == 28: colormap = 6
        elif colormap == 29: colormap = 7
    elif colormap == 30 or colormap == 31: # Greys
        colormapC = ColorMaps.Greens
        if colormap == 30: colormap = 10
        elif colormap == 31: colormap = 11

    return [colormap, colormapC1, colormapC2, colormapC3, colormapC]

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
    if colormap != -1:
        [colormap, colormapC1, colormapC2, colormapC3, colormapC] = filterColormap( [colormap, colormapC1, colormapC2, colormapC3, colormapC] )
    if offscreen == 1 or offscreen == 5 or offscreen == 6 or offscreen == 7: # must set in the right module
        from . import cplotOSMesa as cplot
    else:
        from . import cplot
    cplot.setState(dim, mode, scalarField, vectorField1, vectorField2,
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
                   export, exportResolution, continuousExport,
                   envmap, message,
                   stereo, stereoDist,
                   cursor, gridSize, timer, selectionStyle,
                   activateShortCuts, billBoards, billBoardSize,
                   materials, bumpMaps, frameBuffer, offscreen)

def setMode(mode):
    """Set CPlot display mode.
    Usage: setMode(0)"""
    from . import cplot
    cplot.setMode(mode)

def changeVariable():
    """Change displayed variable.
    Usage: changeVariable()"""
    from . import cplot
    cplot.changeVariable()

def changeStyle():
    """Change CPlot display style.
    Usage: changeStyle()"""
    from . import cplot
    cplot.changeStyle()

def changeInfoDisplay():
    """Change CPlot info display style.
    Usage: changeInfoDisplay()"""
    from . import cplot
    cplot.changeInfoDisplay()

def changeBlanking():
    """Change the blanking procedure.
    Usage: changeBlanking()"""
    from . import cplot
    cplot.changeBlanking()

def setDim(dim):
    """Set CPlot display dim 3, 2 or 1.
    Usage: setDim(2)"""
    from . import cplot
    cplot.setDim(dim)

def setActivePoint(x,y,z):
    """Set the active (clicked) point in plotter.
    Usage: setActivePoint(x,y,z)"""
    from . import cplot
    return cplot.setActivePoint(x,y,z)

def setSelectedZones(zlist):
    """Set selected zones.
    Usage: setSelectedZones([(0,1),(1,1),...])"""
    from . import cplot
    cplot.setSelectedZones(zlist)

def unselectAllZones():
    """Unselect all zones.
    Usage: unselectAllZones()"""
    from . import cplot
    cplot.unselectAllZones()

def setActiveZones(zlist):
    """Set active (displayed) zones.
    Usage: setActiveZones([(0,1)])"""
    from . import cplot
    cplot.setActiveZones(zlist)

def setZoneNames(zlist):
    """Set zone names.
    Usage: setZoneNames([(0,'myZone')])"""
    from . import cplot
    cplot.setZoneNames(zlist)

#==============================================================================
def lookFor():
    """Look for selected zones.
    Usage: lookFor()"""
    from . import cplot
    cplot.lookFor()

def fitView():
    """Fit the view to objects.
    Usage: fitView()"""
    from . import cplot
    cplot.fitView()

def finalizeExport(action=0):
    """Finalize export."""
    if action == 1 or action == 5 or action == 6 or action == 7:
        from . import cplotOSMesa
        cplotOSMesa.finalizeExport(action)
        return
    from . import cplot
    while cplot.isDisplayRunning() == 0: pass
    cplot.finalizeExport(action)

def hide():
    """Hide window."""
    from . import cplot
    cplot.hide()

def show():
    """Show window if it has been hidden with flush."""
    from . import cplot
    cplot.show()

def setWindowTitle(file, path):
    """Set the CPlot window title."""
    from . import cplot
    cplot.setWindowTitle(file, path)

#==============================================================================
# camera
#==============================================================================
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
# image
#==============================================================================
def blur(a, blurSigma=0.8):
    """Blur an image array."""
    from . import cplot
    if isinstance(a[0], list):
        for i in a: cplot.blur(i, blurSigma)
    else: cplot.blur(a, blurSigma)
    return None

#==============================================================================
# -- Internal functions --
def setFileName__(name):
    from . import cplot
    cplot.setFileName(name)

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
                 export, exportResolution, zoneNames, renderTags, frameBuffer, offscreen,
                 posCamList, posEyeList, dirCamList):
    global __slot__
    import threading
    if colormap != -1:
        [colormap, colormapC1, colormapC2, colormapC3, colormapC] = filterColormap( [colormap, colormapC1, colormapC2, colormapC3, colormapC] )
    if offscreen > 0: daemon = True
    else: daemon = False
    from . import cplot
    if version_info[0] == 2: # python 2 - no daemon
        a = threading.Thread(None, cplot.displayNew, None,
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
                              export, exportResolution,
                              zoneNames, renderTags, frameBuffer, offscreen,
                              posCamList, posEyeList, dirCamList), {})
    else: # python3 - daemon
        a = threading.Thread(None, cplot.displayNew, None,
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
                              export, exportResolution,
                              zoneNames, renderTags, frameBuffer, offscreen,
                              posCamList, posEyeList, dirCamList), {}, daemon=daemon)
    a.start()
    __slot__ = a

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
                   export, exportResolution, zoneNames, renderTags, frameBuffer, offscreen,
                   posCamList, posEyeList, dirCamList):
    if colormap != -1:
        [colormap, colormapC1, colormapC2, colormapC3, colormapC] = filterColormap( [colormap, colormapC1, colormapC2, colormapC3, colormapC] )
    from . import cplot
    cplot.displayAgain(arrays, dim, mode, scalarField, vectorField1,
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
                       export, exportResolution, zoneNames, renderTags,
                       frameBuffer, offscreen,
                       posCamList, posEyeList, dirCamList)
    time.sleep(__timeStep__)
