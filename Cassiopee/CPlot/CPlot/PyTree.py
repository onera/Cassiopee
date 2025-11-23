"""Simple plotter for CFD.
"""
import numpy
import os.path
from . import CPlot

# Separateur intra-nom
SEP1 = '/'

__version__ = CPlot.__version__
# Delay between display
__timeStep__ = CPlot.__timeStep__
# CPlot display fields localized in __LOCATION__ ('nodes' or 'centers')
__LOCATION__ = 'nodes'

# FILE name for default decorator output
decorator = '.decorator.png'

try:
    import Converter.PyTree as C
    import Converter.Internal as Internal
except:
    raise ImportError("CPlot.PyTree: requires Converter.PyTree module.")

#==============================================================================
def display(t,
            dim=-1,
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
            vectorScale=-1., vectorDensity=-1., vectorNormalize=-1,
            vectorShowSurface=-1, vectorShape=-1, vectorProjection=-1,
            colormap=-1,
            colormapC1="", colormapC2="", colormapC3="", colormapC=None,
            niso=25,
            isoEdges=-0.5,
            isoScales=[],
            win=(-1,-1),
            posCam=(-999,-999,-999),
            posEye=(-999,-999,-999),
            dirCam=(-999,-999,-999),
            viewAngle=-1.,
            bgColor=-1,
            backgroundFile="None",
            shadow=-1, lightOffset=(-999,-999),
            dof=-1, dofPower=-1, gamma=-1, toneMapping=-1,
            stereo=-1, stereoDist=-1., panorama=0,
            export="None", exportResolution="None", exportAA=-1,
            location="unchanged",
            frameBuffer=-1,
            offscreen=0,
            posCamList=None, posEyeList=None, dirCamList=None):
    """Display pyTrees.
    Usage: display(t)"""
    global __LOCATION__
    if location != 'unchanged': __LOCATION__ = location
    if __LOCATION__ == 'nodes':
        t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    else: t = C.node2Center(t)
    zoneNames = C.getZoneNames(t)
    renderTags = getRenderTags(t)
    arrays = C.getAllFields(t, 'nodes', api=3)
    CPlot.display(arrays, dim, mode, scalarField, vectorField1,
                  vectorField2, vectorField3, displayBB, displayInfo,
                  displayIsoLegend, meshStyle,
                  solidStyle, scalarStyle, vectorStyle,
                  vectorScale, vectorDensity, vectorNormalize,
                  vectorShowSurface, vectorShape, vectorProjection,
                  colormap, colormapC1, colormapC2, colormapC3, colormapC,
                  niso, isoEdges, isoScales, win,
                  posCam, posEye, dirCam, viewAngle,
                  bgColor, backgroundFile,
                  shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                  stereo, stereoDist, panorama,
                  export, exportResolution, exportAA,
                  zoneNames, renderTags, frameBuffer, offscreen,
                  posCamList, posEyeList, dirCamList)

# CB: temporaire. Raw data direct display.
def displayRaw(t,
               dim=-1,
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
               vectorScale=-1., vectorDensity=-1., vectorNormalize=-1,
               vectorShowSurface=-1, vectorShape=-1, vectorProjection=-1,
               colormap=-1,
               colormapC1="", colormapC2="", colormapC3="", colormapC=None,
               niso=25,
               isoEdges=-0.5,
               isoScales=[],
               win=(-1,-1),
               posCam=(-999,-999,-999),
               posEye=(-999,-999,-999),
               dirCam=(-999,-999,-999),
               viewAngle=-1.,
               bgColor=-1,
               backgroundFile="None",
               shadow=-1, lightOffset=(-999,-999),
               dof=-1, dofPower=-1, gamma=-1, toneMapping=-1,
               stereo=-1, stereoDist=-1., panorama=0,
               export="None", exportResolution="None", exportAA=-1,
               location="unchanged",
               frameBuffer=-1,
               offscreen=0):
    zoneNames = C.getZoneNames(t)
    renderTags = getRenderTags(t)
    arrays = C.getAllFields(t, 'nodes', api=3)
    CPlot.display(arrays, dim, mode, scalarField, vectorField1,
                  vectorField2, vectorField3, displayBB, displayInfo,
                  displayIsoLegend, meshStyle,
                  solidStyle, scalarStyle, vectorStyle,
                  vectorScale, vectorDensity, vectorNormalize,
                  vectorShowSurface, vectorShape, vectorProjection,
                  colormap, colormapC1, colormapC2, colormapC3, colormapC,
                  niso, isoEdges, isoScales, win,
                  posCam, posEye, dirCam, viewAngle,
                  bgColor, backgroundFile,
                  shadow, lightOffset, dof, dofPower, gamma, toneMapping,
                  stereo, stereoDist, panorama,
                  export, exportResolution, exportAA,
                  zoneNames, renderTags, frameBuffer, offscreen)

#==============================================================================
def render():
    """Force render.
    Usage: render()"""
    CPlot.render()

#==============================================================================
# This function doesnt remove zone from t
def delete(zlist):
    """Delete zones from plotter."""
    CPlot.delete(zlist)

#==============================================================================
def add(t, nob, noz, zone):
    """Add/insert a zone to plotter.
    Usage: add(t, nob, noz, zone)"""
    zoneName = t[2][nob][0]+SEP1+zone[0]
    renderTag = getRenderTags__(zone, [])[0]

    if noz == -1: noz = len(t[2][nob][2]) # insere a la fin
    t[2][nob][2].insert(noz, zone)
    if CPlot.__slot__ is None: display(t); return
    (nzs, nzu) = getNzs(t, zone)

    if __LOCATION__ == 'nodes': # ajouter le toptree
        zone = C.center2Node(zone, Internal.__FlowSolutionCenters__)
    else: zone = C.node2Center(zone)
    array = C.getAllFields(zone, 'nodes', api=3)[0]

    CPlot.cplotm.add(array, (nzs, nzu), zoneName, renderTag)

#==============================================================================
def replace(t, nob, noz, zone):
    """Replace a zone in plotter.
    Usage: replace(t, nob, noz, zone)"""
    zoneName = t[2][nob][0]+SEP1+zone[0]
    renderTag = getRenderTags__(zone, [])[0]
    oldType = Internal.getZoneType(t[2][nob][2][noz])
    t[2][nob][2][noz] = zone
    if CPlot.__slot__ is None: display(t); return
    (nzs, nzu) = getNzs(t, zone)

    if __LOCATION__ == 'nodes': # ajouter le toptree
        zone = C.center2Node(zone, Internal.__FlowSolutionCenters__)
    else: zone = C.node2Center(zone)
    array = C.getAllFields(zone, 'nodes', api=3)[0]

    CPlot.cplotm.replace(array, (nzs, nzu, oldType), zoneName, renderTag)

#==============================================================================
def display1D(t, slot=0, gridPos=(0,0), gridSize=(-1,-1),
              bgBlend=1., var1='', var2='',
              r1=None, r2=None):
    """Display 1D plots.
    Usage: display1D(t, slot, ....)"""
    import Converter
    if len(t) > 0 and isinstance(t[0], numpy.ndarray): # as numpys
        if len(t) < 2:
            raise ValueError('display1D: requires at least two numpys [x,y]')
        x = t[0]; y = t[1]
        n = x.size
        array = Converter.array(var1+','+var2,n,1,1)
        array[1][0,:] = x[:]
        array[1][1,:] = y[:]
        arrays = [array]
    else:
        arrays = C.getAllFields(t, 'nodes', api=1) # as pyZone
    CPlot.display1D(arrays, slot, gridPos, gridSize, bgBlend,
                    var1, var2, r1, r2)

#==============================================================================
def pressKey():
    """Wait for user to press a key.
    Usage: pressKey()"""
    CPlot.pressKey()

#==============================================================================
# -- get functions --
#==============================================================================
def getState(mode):
    """Return the state in plotter.
    Usage: n = getState(mode)"""
    return CPlot.getState(mode)

def getSelectedZone():
    """Return the selected zone in plotter.
    Usage: n=getSelectedZone()"""
    return CPlot.getSelectedZone()

def getSelectedZones():
    """Return the selected zones in plotter.
    Usage: list=getSelectedZones()"""
    return CPlot.getSelectedZones()

def getSelectedStatus(zone):
    """Return the selected status of a zone in plotter.
    Usage: status=getSelectedStatus(zone)"""
    return CPlot.getSelectedStatus(zone)

def getActiveZones():
    """Return the active (displayed) zones in plotter.
    Usage: list=getActiveZones()"""
    return CPlot.getActiveZones()

def getActiveStatus(zone):
    """Return the active status of a zone in plotter.
    Usage: status=getActiveStatus(zone)"""
    return CPlot.getActiveStatus(zone)

def getActivePoint():
    """Return the active point coordinates in plotter.
    Usage: n=getActivePoint()"""
    return CPlot.getActivePoint()

def getActivePointIndex():
    """Return the active point index.
    Usage: n = getActivePointIndex()"""
    return CPlot.getActivePointIndex()

def getActivePointF():
    """Return the active point field values.
    Usage: n = getActivePointF()"""
    return CPlot.getActivePointF()

def getMouseState():
    """Return mouse state (mouse position and button state)."""
    return CPlot.getMouseState()

def getKeyboard():
    """Return the pressed keys.
    Usage: n = getKeyboard()"""
    return CPlot.getKeyboard()

def resetKeyboard():
    """Reset the keyboard string.
    Usage: resetKeyboard()"""
    return CPlot.resetKeyboard()

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
             colormapC1="",
             colormapC2="",
             colormapC3="",
             colormapC=None,
             niso=-1,
             isoEdges=-1,
             isoScales=[],
             win=(-1,-1),
             posCam=(-999,-999,-999),
             posEye=(-999,-999,-999),
             dirCam=(-999,-999,-999),
             viewAngle=-1.,
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
             gridSize=(-1,-1),
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
    CPlot.setState(dim, mode, scalarField, vectorField1, vectorField2,
                   vectorField3, displayBB, displayInfo, displayIsoLegend,
                   meshStyle, solidStyle, scalarStyle,
                   vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                   vectorShowSurface, vectorShape, vectorProjection,
                   colormap, colormapC1, colormapC2, colormapC3, colormapC,
                   niso, isoEdges, isoScales, win,
                   posCam, posEye, dirCam, viewAngle, lightOffset,
                   bgColor, backgroundFile,
                   shadow, dof, dofPower, gamma, toneMapping,
                   sobelThreshold, sharpenPower, ssaoPower,
                   ghostifyDeactivatedZones, edgifyActivatedZones,
                   edgifyDeactivatedZones, simplifyOnDrag,
                   export, exportResolution, exportAA,
                   continuousExport, envmap, message,
                   stereo, stereoDist,
                   cursor, gridSize, timer, selectionStyle,
                   activateShortCuts, billBoards, billBoardSize,
                   materials, bumpMaps, frameBuffer, offscreen)

def setMode(mode):
    """Set CPlot display mode.
    Usage: setMode(0)"""
    CPlot.setMode(mode)

def changeVariable():
    """Change displayed variable.
    Usage: changeVariable()"""
    CPlot.changeVariable()

def changeStyle():
    """Change CPlot display style.
    Usage: changeStyle()"""
    CPlot.changeStyle()

def changeInfoDisplay():
    """Change CPlot info display style.
    Usage: changeInfoDisplay()"""
    CPlot.changeInfoDisplay()

def changeBlanking():
    """Change CPlot blanking procedure.
    Usage: changeBlanking()"""
    CPlot.changeBlanking()

def setDim(dim):
    """Set CPlot display dim 3, 2 or 1.
    Usage: setDim(2)"""
    CPlot.setDim(dim)

def setActivePoint(x,y,z):
    """Set the active (clicked) point in plotter.
    Usage: setActivePoint(x,y,z)"""
    return CPlot.setActivePoint(x,y,z)

def setSelectedZones(zlist):
    """Set selected zones.
    Usage: setSelectedZones([(0,1),(1,1),...])"""
    CPlot.setSelectedZones(zlist)

def unselectAllZones():
    """Unselect all zones.
    Usage: unselectAllZones()"""
    CPlot.unselectAllZones()

def setActiveZones(zlist):
    """Set active (displayed) zones.
    Usage: setActiveZones([(0,1)])"""
    CPlot.setActiveZones(zlist)

def setZoneNames(zlist):
    """Set zone names.
    Usage: setZoneNames([(0,'myZone')])"""
    CPlot.setZoneNames(zlist)

def lookFor():
    """Look for selected zones.
    Usage: lookFor()"""
    CPlot.lookFor()

def fitView():
    """Fit view to all zones.
    Usage: firView()"""
    CPlot.fitView()

def finalizeExport(action=0):
    """Finalize export (barrier, end movie, stop continuous export."""
    CPlot.finalizeExport(action)

def hide():
    """Hide window."""
    CPlot.hide()

def show():
    """Show window if it has been hidden with flush."""
    CPlot.show()

def moveCamera(posCams, posEyes=None, dirCams=None, moveEye=False, N=100, speed=1., pos=-1):
    """Move camera.
    Usage: moveCamera(checkPoints, moveEye, N, speed, pos)."""
    if isinstance(posCams[0], str): # zone
        posCams = C.getAllFields(posCams, 'nodes', api=1)[0]
    if posEyes is not None and isinstance(posEyes[0], str): # zone
        posEyes = C.getAllFields(posEyes, 'nodes', api=1)[0]
    if dirCams is not None and isinstance(dirCams[0], str): # zone
        dirCams = C.getAllFields(dirCams, 'nodes', api=1)[0]
    ret = CPlot.moveCamera(posCams, posEyes, dirCams, moveEye, N, speed, pos)
    return ret

def travelRight(xr=0.1, N=100):
    """Travel camera right."""
    ret = CPlot.travelRight(xr, N)
    return ret

def travelLeft(xr=0.1, N=100):
    """Travel camera left."""
    ret = CPlot.travelLeft(xr, N)
    return ret

def travelUp(xr=0.1, N=100):
    """Travel camera up."""
    ret = CPlot.travelUp(xr, N)
    return ret

def travelDown(xr=0.1, N=100):
    """Travel camera down."""
    ret = CPlot.travelDown(xr, N)
    return ret

def travelIn(xr=0.1, N=100):
    """Zoom camera in."""
    ret = CPlot.travelIn(xr, N)
    return ret

def travelOut(xr=0.1, N=100):
    """Zoom camera out."""
    ret = CPlot.travelOut(xr, N)
    return ret

def setWindowTitle(file, path):
    """Set CPlot window title."""
    CPlot.setWindowTitle(file, path)

#==============================================================================
# -- Numbering functions --
#==============================================================================
# Retourne la numerotation (base, zone) dans l'arbre
# a partir de la numerotation zone de CPlot
#==============================================================================
def updateCPlotNumbering(t):
    """Update the CPlot numbering."""
    Nz = []; Nb = []
    nodes = Internal.getZones(t)
    nstrf = 0
    for z in nodes:
        Nz.append(1); Nb.append(1)
        ztype = Internal.getZoneType(z)
        if ztype == 1: nstrf += 1

    nstr = 0; nunstr = 0
    nb = -1
    for b in t[2]:
        if b[3] == 'CGNSBase_t':
            nz = 0
            for z in b[2]:
                if z[3] == 'Zone_t':
                    ztype = Internal.getZoneType(z)
                    if ztype == 1:
                        Nz[nstr] = nz
                        Nb[nstr] = nb; nstr += 1
                    else:
                        Nz[nunstr+nstrf] = nz
                        Nb[nunstr+nstrf] = nb; nunstr += 1
                nz += 1
        nb += 1
    return (Nb, Nz)

#==============================================================================
# Retourne le numero global nz de la basename/zonename dans CPlot
#==============================================================================
def updateCPlotGlobalNumbering(t):
    dnz = {}
    # count total struct zones
    nstr = 0
    for b in t[2]:
        if b[3] == 'CGNSBase_t':
            for z in b[2]:
                if z[3] == 'Zone_t':
                    type = Internal.getZoneType(z)
                    if type == 1: nstr += 1

    nzs = 0 # global cplot numbering
    nzu = nstr
    for b in t[2]:
        if b[3] == 'CGNSBase_t':
            dnz[b[0]] = {}
            for z in b[2]:
                if z[3] == 'Zone_t':
                    type = Internal.getZoneType(z)
                    if type == 1: dnz[b[0]][z[0]] = nzs; nzs += 1
                    else: dnz[b[0]][z[0]] = nzu; nzu += 1
    return dnz

#==============================================================================
# Retourne la numerotation CPlot de la zone baseName/zoneName de t
# IN: t: un noeud pyTree (et pas autre chose)
# IN: baseName: le nom de la base
# IN: zoneName: le nom de la zone
#==============================================================================
def getCPlotNumber(t, baseName, zoneName):
    """Return the CPlot number of zone defined by baseName/zoneName."""
    base = Internal.getChildFromName(t, baseName) # must be unique
    z = Internal.getChildFromName(base, zoneName) # must be unique
    type = Internal.getZoneType(z)
    bases = Internal.getBases(t)

    if type == 1: # zoneName is structured
        nstr = 0
        for b in bases:
            nodes = Internal.getNodesFromType1(b, 'Zone_t')
            if b[0] != baseName:
                for z in nodes:
                    type = Internal.getZoneType(z)
                    if type == 1: nstr += 1
            else:
                for z in nodes:
                    type = Internal.getZoneType(z)
                    if type == 1: nstr += 1
                    if z[0] == zoneName: return nstr-1
        return nstr-1
    else: # zoneName is unstructured
        nstr = 0
        for b in bases:
            nodes = Internal.getNodesFromType1(b, 'Zone_t')
            for z in nodes:
                type = Internal.getZoneType(z)
                if type == 1: nstr += 1
        unstr = nstr
        for b in bases:
            nodes = Internal.getNodesFromType1(b, 'Zone_t')
            if b[0] != baseName:
                for z in nodes:
                    type = Internal.getZoneType(z)
                    if type == 2: unstr += 1
            else:
                for z in nodes:
                    type = Internal.getZoneType(z)
                    if type == 2: unstr += 1
                    if z[0] == zoneName: return unstr-1
        return unstr-1

#==============================================================================
# Retourne le nombre de zones structurees (nzs) et le nombre de zones
# non structurees (nzu) qui sont avant zone dans l'arbre t (sans compter zone)
#==============================================================================
def getNzs(t, zone):
    """Get the number of structured zones and untructured zones in t before zone."""
    nzs = 0; nzu = 0
    nodes = Internal.getZones(t)
    for z in nodes:
        if id(z) == id(zone): return (nzs, nzu)
        type = Internal.getZoneType(z)
        if type == 1: nzs += 1
        else: nzu += 1
    return (nzs, nzu)

#==============================================================================
# Delete la selection (CPlot) dans t
# IN: t: pyTree
# IN: Nb, Nz: CPlot numbering
# IN: nzs: selection CPlot
# OUT: t: avec les zones correspondant a la selection supprimees
#==============================================================================
def deleteSelection(t, Nb, Nz, nzs):
    """Delete a selection of zones."""
    nbases = len(t[2])
    bases = []
    for i in range(nbases): bases.append([])
    for nz in nzs:
        nob = Nb[nz]+1
        noz = Nz[nz]
        bases[nob].append(noz)

    for i in range(nbases):
        l = bases[i]
        for a in range(len(l)):
            for b in range(a+1, len(l)):
                if l[a] < l[b]:
                    temp = l[a]; l[a] = l[b]; l[b] = temp

    for i in range(1, nbases):
        l = bases[i]
        for a in range(len(l)):
            del t[2][i][2][l[a]]
    return t

#==============================================================================
# Determine si la selection correspond a une base complete
# Retourne le no de la base si c'est le cas
# Retourne -1 sinon
#==============================================================================
def isSelAFullBase(t, Nb, nzs):
    """Return true if selection corresponds to a full base."""
    if nzs == []: return -1
    fullBase = 0
    for nz in nzs:
        nob = Nb[nz]+1
        if fullBase == 0: fullBase = nob
        if fullBase != nob: fullBase = -1
    if fullBase > 0:
        zones = Internal.getNodesFromType1(t[2][fullBase], 'Zone_t')
        if len(nzs) != len(zones): fullBase = -1
    return fullBase

#==============================================================================
# getRenderTags
# Retourne les tags de render dans l'arbre, trie ainsi:
# En premier, les zones structurees
# En second, les zones non structurees
# Un tag est de la forme:
# color:material:blending:meshOverlay:meshColor:meshWidth:shaderParameters,
# si renderInfo n'est pas present, on retourne None:None:None:None:None:None:None.
#==============================================================================
def getRenderTags(t):
    """Return the render tags of zones."""
    bases = Internal.getBases(t)
    renderTags = []
    if bases == []:
        zones = Internal.getNodesFromType1(t, 'Zone_t')
        for z in zones:
            type = Internal.getZoneType(z)
            if type == 1: renderTags = getRenderTags__(z, renderTags)
        for z in zones:
            type = Internal.getZoneType(z)
            if type == 2: renderTags = getRenderTags__(z, renderTags)
    else:
        for b in bases:
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                type = Internal.getZoneType(z)
                if type == 1: renderTags = getRenderTags__(z, renderTags)
        for b in bases:
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            for z in zones:
                type = Internal.getZoneType(z)
                if type == 2: renderTags = getRenderTags__(z, renderTags)
    return renderTags

#==============================================================================
# IN: z: zone node
#==============================================================================
def getRenderTags__(z, renderTags):
    ri = Internal.getNodeFromName1(z, '.RenderInfo')
    if ri is None: renderTags.append('None:None:None:None:None:None:None:None')
    else:
        rc = Internal.getNodeFromName1(ri, 'Color') # Zone color
        if rc is None: color = 'None'
        else: color = Internal.getValue(rc)
        rm = Internal.getNodeFromName1(ri, 'Material') # Zone material
        if rm is None: material = 'None'
        else: material = Internal.getValue(rm)
        rm = Internal.getNodeFromName1(ri, 'Blending') # Zone blending
        if rm is None: blending = 'None'
        else: blending = str(Internal.getValue(rm))
        rm = Internal.getNodeFromName1(ri, 'MeshOverlay') # Mesh overlay on/off
        if rm is None: meshOverlay = 'None'
        else: meshOverlay = str(Internal.getValue(rm))
        rm = Internal.getNodeFromName1(ri, 'MeshColor') # Mesh overlay color
        if rm is None: meshColor = 'None'
        else: meshColor = Internal.getValue(rm)
        rm = Internal.getNodeFromName1(ri, 'MeshWidth') # Mesh overlay width
        if rm is None: meshWidth = 'None'
        else: meshWidth = str(Internal.getValue(rm))
        rm = Internal.getNodeFromName1(ri, 'ShaderParameters')
        if rm is None: shaderParameters = 'None:None'
        else:
            v = rm[1]
            if isinstance(v, numpy.ndarray):
                shaderParameters = ''
                lgt = v.shape[0]
                for i in range(lgt):
                    shaderParameters += str(v[i])
                    if i < lgt-1: shaderParameters += ':'
            else: shaderParameters = 'None:None'
        renderTags.append(color+':'+material+':'+blending+':'+meshOverlay+':'+meshColor+':'+meshWidth+':'+shaderParameters)
    return renderTags

# -- addRender2Zone
def addRender2Zone(t, material=None, color=None, blending=None,
                   meshOverlay=None, meshColor=None, meshWidth=None,
                   shaderParameters=None):
    """Add a renderInfo node to a zone node.
    Usage: addRender2Zone(zone, renderInfo)"""
    tp = Internal.copyRef(t)
    _addRender2Zone(tp, material, color, blending,
                    meshOverlay, meshColor, meshWidth, shaderParameters)
    return tp

def _addRender2Zone(a, material=None, color=None, blending=None,
                    meshOverlay=None, meshColor=None, meshWidth=None,
                    shaderParameters=None):
    """Add a renderInfo node to a zone node.
    Usage: addRender2Zone(zone, renderInfo)"""
    zones = Internal.getZones(a)
    for z in zones:

        ri = Internal.createUniqueChild(z, '.RenderInfo', 'UserDefinedData_t')

        if material is not None:
            Internal._createUniqueChild(ri, 'Material', 'DataArray_t', value=material)

        if color is not None:
            Internal._createUniqueChild(ri, 'Color', 'DataArray_t', value=color)

        if blending is not None:
            Internal._createUniqueChild(ri, 'Blending', 'DataArray_t', value=blending)

        if meshOverlay is not None:
            Internal._createUniqueChild(ri, 'MeshOverlay', 'DataArray_t', value=meshOverlay)

        if meshColor is not None:
            Internal._createUniqueChild(ri, 'MeshColor', 'DataArray_t', value=meshColor)

        if meshWidth is not None:
            Internal._createUniqueChild(ri, 'MeshWidth', 'DataArray_t', value=meshWidth)

        if shaderParameters is not None:
            Internal._createUniqueChild(ri, 'ShaderParameters', 'DataArray_t', value=shaderParameters)
    return None

# -- addRender2PyTree
def addRender2PyTree(t, slot=0, posCam=None, posEye=None, dirCam=None,
                     mode=None, scalarField=None, niso=None, isoScales=None,
                     isoEdges=None, isoLight=None, isoLegend=None,
                     colormap=None, colormapC1=None, colormapC2=None, colormapC3=None, colormapC=None,
                     materials=None, bumpMaps=None, billBoards=None,
                     shadow=None, lightOffsetX=None, lightOffsetY=None,
                     dof=None, tone=None, gamma=None, dofPower=None,
                     sharpenPower=None, ssaoPower=None):
    """Add a renderInfo node to a tree.
    Usage: addRender2PyTree(t, slot, renderInfo)"""
    a = Internal.copyRef(t)
    _addRender2PyTree(a, slot, posCam, posEye, dirCam,
                      mode, scalarField, niso, isoScales,
                      isoEdges, isoLight, isoLegend,
                      colormap, colormapC1, colormapC2, colormapC3, colormapC,
                      materials, bumpMaps, billBoards,
                      shadow, lightOffsetX, lightOffsetY,
                      dof, tone, gamma, dofPower,
                      sharpenPower, ssaoPower)
    return a

def _addRender2PyTree(a, slot=0, posCam=None, posEye=None, dirCam=None,
                      mode=None, scalarField=None, niso=None, isoScales=None,
                      isoEdges=None, isoLight=None, isoLegend=None,
                      colormap=None, colormapC1=None, colormapC2=None, colormapC3=None, colormapC=None,
                      materials=None, bumpMaps=None, billBoards=None,
                      shadow=None, lightOffsetX=None, lightOffsetY=None,
                      dof=None, tone=None, gamma=None, dofPower=None,
                      sharpenPower=None, ssaoPower=None):
    """Add a renderInfo node to a tree.
    Usage: addRender2PyTree(t, renderInfo)"""
    if a[3] != 'CGNSTree_t': return None
    exist = 0
    for i in a[2]:
        if i[0] == '.RenderInfo': break
        exist += 1
    if exist < len(a[2]): ri = a[2][exist]
    else: ri = Internal.createChild(a, '.RenderInfo', 'UserDefinedData_t')

    # find slot
    sl = Internal.getNodeFromName1(ri, 'Slot%d'%slot)
    if sl is None:
        sl = Internal.createChild(ri, 'Slot%d'%slot, 'UserDefinedData_t')

    if posCam is not None:
        rt = Internal.createUniqueChild(sl, 'posCam', 'DataArray_t', value=posCam)

    if posEye is not None:
        rt = Internal.createUniqueChild(sl, 'posEye', 'DataArray_t', value=posEye)

    if dirCam is not None:
        rt = Internal.createUniqueChild(sl, 'dirCam', 'DataArray_t', value=dirCam)

    if mode is not None:
        rt = Internal.createUniqueChild(sl, 'mode', 'DataArray_t', value=mode)

    if scalarField is not None:
        rt = Internal.createUniqueChild(sl, 'scalarField', 'DataArray_t', value=scalarField)

    if niso is not None:
        rt = Internal.createUniqueChild(sl, 'niso', 'DataArray_t', value=niso)

    if isoScales is not None: # must be a list or a list of list
        if len(isoScales) > 0 and isinstance(isoScales[0], list): # list of list
            for iso in isoScales:
                n = len(iso)-1
                v = numpy.empty(n, numpy.float64)
                for i in range(n): v[i] = float(iso[i+1])
                Internal.createUniqueChild(sl, 'isoScales[%s]'%str(iso[0]), 'DataArray_t', value=v)
        else:
            n = len(isoScales)-1
            v = numpy.empty(n, numpy.float64)
            for i in range(n): v[i] = float(isoScales[i+1])
            Internal.createUniqueChild(sl, 'isoScales[%s]'%str(isoScales[0]), 'DataArray_t', value=v)

    if isoEdges is not None:
        rt = Internal.createUniqueChild(sl, 'isoEdges', 'DataArray_t', value=isoEdges)

    if isoLight is not None:
        rt = Internal.createUniqueChild(sl, 'isoLight', 'DataArray_t', value=isoLight)

    if isoLegend is not None:
        rt = Internal.createUniqueChild(sl, 'isoLegend', 'DataArray_t', value=isoLegend)

    if colormap is not None:
        rt = Internal.createUniqueChild(sl, 'colormap', 'DataArray_t', value=colormap)
    if colormapC1 is not None:
        rt = Internal.createUniqueChild(sl, 'colormapC1', 'DataArray_t', value=colormapC1)
    if colormapC2 is not None:
        rt = Internal.createUniqueChild(sl, 'colormapC2', 'DataArray_t', value=colormapC2)
    if colormapC3 is not None:
        rt = Internal.createUniqueChild(sl, 'colormapC3', 'DataArray_t', value=colormapC3)
    if colormapC is not None:
        rt = Internal.createUniqueChild(sl, 'colormapC', 'DataArray_t', value=colormapC)

    # Under .RenderInfo
    if materials is not None:
        rt = Internal.createUniqueChild(ri, 'materials', 'UserDefinedData_t')
        prevValues = [Internal.getValue(i) for i in rt[2]]
        cnt = len(prevValues)
        li = prevValues
        for i in materials:
            if i not in prevValues: li.append(i)
        for c, f in enumerate(li):
            Internal._createUniqueChild(rt, 'file%d'%c, 'DataArray_t', value=f)

    if bumpMaps is not None:
        rt = Internal.createUniqueChild(ri, 'bumpMaps', 'UserDefinedData_t')
        prevValues = [Internal.getValue(i) for i in rt[2]]
        cnt = len(prevValues)
        li = prevValues
        for i in bumpMaps:
            if i not in prevValues: li.append(i)
        for c, f in enumerate(li):
            Internal._createUniqueChild(rt, 'file%d'%c, 'DataArray_t', value=f)

    if billBoards is not None:
        rt = Internal.createUniqueChild(ri, 'billBoards', 'UserDefinedData_t')
        prevValues = [Internal.getValue(i) for i in rt[2]]
        cnt = len(prevValues)
        li = prevValues
        for i in billBoards:
            if i not in prevValues: li.append(i)
        for c, f in enumerate(li):
            Internal._createUniqueChild(rt, 'file%d'%(c+cnt), 'DataArray_t', value=f)

    if shadow is not None:
        rt = Internal.createUniqueChild(sl, 'shadow', 'UserDefinedData_t', value=int(shadow))

    if lightOffsetX is not None:
        rt = Internal.createUniqueChild(sl, 'lightOffsetX', 'UserDefinedData_t', value=lightOffsetX)

    if lightOffsetY is not None:
        rt = Internal.createUniqueChild(sl, 'lightOffsetY', 'UserDefinedData_t', value=lightOffsetY)

    if dof is not None: # dof activate general post-processing of frame buffer
        rt = Internal.createUniqueChild(sl, 'dof', 'UserDefinedData_t', value=dof)

    if tone is not None:
        rt = Internal.createUniqueChild(sl, 'tone', 'UserDefinedData_t', value=tone)

    if gamma is not None:
        rt = Internal.createUniqueChild(sl, 'gamma', 'UserDefinedData_t', value=gamma)

    if dofPower is not None:
        rt = Internal.createUniqueChild(sl, 'dofPower', 'UserDefinedData_t', value=dofPower)

    if sharpenPower is not None:
        rt = Internal.createUniqueChild(sl, 'sharpenPower', 'UserDefinedData_t', value=sharpenPower)

    if ssaoPower is not None:
        rt = Internal.createUniqueChild(sl, 'ssaoPower', 'UserDefinedData_t', value=ssaoPower)

    return None

#==============================================================================
# IN: nom de la colormap
# IN: light: 0 ou 1
# OUT: style (entier)
#==============================================================================
def colormap2Style(colormapName, light=0):
    style = 0
    if colormapName == 'Blue2Red': style = 0
    elif colormapName == 'BiColorRGB': style = 2
    elif colormapName == 'BiColorHSV': style = 4
    elif colormapName == 'TriColorRGB': style = 6
    elif colormapName == 'TriColorHSV': style = 8
    elif colormapName == 'Diverging': style = 14
    elif colormapName == 'Viridis': style = 16
    elif colormapName == 'Inferno': style = 18
    elif colormapName == 'Magma': style = 20
    elif colormapName == 'Plasma': style = 22
    elif colormapName == 'Jet': style = 24
    elif colormapName == 'Greys': style = 26
    elif colormapName == 'NiceBlue': style = 28
    elif colormapName == 'Greens': style = 30
    if light == 1: style += 1
    return style

def style2Colormap(style):
    light = style % 2
    if style == 0 or style == 1: colormap = 'Blue2Red'
    elif style == 2 or style == 3: colormap = 'BiColorRGB'
    elif style == 4 or style == 5: colormap = 'BiColorHSV'
    elif style == 6 or style == 7: colormap = 'TriColorRGB'
    elif style == 8 or style == 9: colormap = 'TriColorHSV'
    elif style == 14 or style == 15: colormap = 'Diverging'
    elif style == 16 or style == 17: colormap = 'Viridis'
    elif style == 18 or style == 19: colormap = 'Inferno'
    elif style == 20 or style == 21: colormap = 'Magma'
    elif style == 22 or style == 23: colormap = 'Plasma'
    elif style == 24 or style == 25: colormap = 'Jet'
    elif style == 26 or style == 27: colormap = 'Greys'
    elif style == 28 or style == 29: colormap = 'NiceBlue'
    elif style == 30 or style == 31: colormap = 'Greens'
    return colormap, light

def getFilteredColormap():
    return CPlot.getFilteredColormap()

#==============================================================================
# loadView from slot
#==============================================================================
def loadView(t, slot=0):
    """Load a view stored in slot."""
    renderInfo = Internal.getNodeFromName1(t, '.RenderInfo')
    if renderInfo is None: return
    slot = Internal.getNodeFromName1(renderInfo, 'Slot%d'%slot)
    if slot is None: return
    pos = Internal.getNodeFromName1(slot, 'posCam')
    if pos is None:
        n = pos[1]; CPlot.setState(posCam=(n[0], n[1], n[2]))
    pos = Internal.getNodeFromName1(slot, 'posEye')
    if pos is not None:
        n = pos[1]; CPlot.setState(posEye=(n[0], n[1], n[2]))
    pos = Internal.getNodeFromName1(slot, 'dirCam')
    if pos is not None:
        n = pos[1]; CPlot.setState(dirCam=(n[0], n[1], n[2]))
    pos = Internal.getNodeFromName1(slot, 'mode')
    if pos is not None:
        mode = Internal.getValue(pos)
        if mode == 'Mesh': imode = 0
        elif mode == 'Solid': imode = 1
        elif mode == 'Render': imode = 2
        elif mode == 'Scalar': imode = 3
        elif mode == 'Vector': imode = 4
        else: imode = 0
        CPlot.setState(mode=imode)
    pos = Internal.getNodeFromName1(slot, 'scalarField')
    if pos is not None:
        field = Internal.getValue(pos)
        avars = C.getVarNames(t)[0]
        ifield = 0
        for i in avars:
            if i == field: break
            if i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ':
                ifield += 1
        CPlot.setState(scalarField=ifield)
    pos = Internal.getNodeFromName1(slot, 'niso')
    if pos is not None:
        niso = Internal.getValue(pos); CPlot.setState(niso=int(niso))
    pos = Internal.getNodeFromName1(slot, 'isoEdges')
    if pos is not None:
        edge = Internal.getValue(pos); CPlot.setState(isoEdges=edge)
    pos = Internal.getNodeFromName1(slot, 'isoLight')
    if pos is not None: light = Internal.getValue(pos)
    else: light = 1
    pos = Internal.getNodeFromName1(slot, 'isoLegend')
    if pos is not None: legend = Internal.getValue(pos)
    else: legend = 0
    pos = Internal.getNodeFromName1(slot, 'colormap')
    if pos is not None: colormap = Internal.getValue(pos)
    else: colormap = 'Blue2Red'
    pos = Internal.getNodeFromName1(slot, 'colormapC1')
    if pos is not None: colormapC1 = Internal.getValue(pos)
    else: colormapC1 = '#000000'
    pos = Internal.getNodeFromName1(slot, 'colormapC2')
    if pos is not None: colormapC2 = Internal.getValue(pos)
    else: colormapC2 = '#FFFFFF'
    pos = Internal.getNodeFromName1(slot, 'colormapC3')
    if pos is not None: colormapC3 = Internal.getValue(pos)
    else: colormapC3 = '#777777'
    pos = Internal.getNodeFromName1(slot, 'colormapC')
    if pos is not None: colormapC = Internal.getValue(pos)
    else: colormapC = []
    style = colormap2Style(colormap, light)

    if legend == 1: CPlot.setState(displayIsoLegend=legend)

    CPlot.setState(colormap=style)
    pos = Internal.getNodesFromName(slot, 'isoScales*')
    scales = []
    for p in pos:
        name = p[0]
        name = name.replace('isoScales[', '')
        name = name[0:-2]
        #try: name = int(name)
        #except: pass
        out = [name]+p[1].tolist()
        scales.append(out)
    if scales != []: CPlot.setState(isoScales=scales)

    # RenderInfo
    pos = Internal.getNodeFromName1(renderInfo, 'materials')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i))
        CPlot.setState(materials=out)
    pos = Internal.getNodeFromName1(renderInfo, 'bumpMaps')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i))
        CPlot.setState(bumpMaps=out)
    pos = Internal.getNodeFromName1(renderInfo, 'billBoards')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i))
        CPlot.setState(billBoards=out)

#==============================================================================
# loadGlobalFiles (material, bumpmaps, billboards)
#==============================================================================
def loadImageFiles(t, offscreen=0):
    """Load image files (texture, billboards, bumpmaps) in CPlot."""
    if t == []: return
    renderInfo = Internal.getNodeFromName1(t, '.RenderInfo')
    if renderInfo is None: return None
    pos = Internal.getNodeFromName1(renderInfo, 'materials')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i))
        CPlot.setState(materials=out, offscreen=offscreen)
    pos = Internal.getNodeFromName1(renderInfo, 'bumpMaps')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i))
        CPlot.setState(bumpMaps=out, offscreen=offscreen)
    pos = Internal.getNodeFromName1(renderInfo, 'billBoards')
    if pos is not None:
        out = []
        for i in pos[2]: out.append(Internal.getValue(i))
        CPlot.setState(billBoards=out, offscreen=offscreen)
    return None

#==============================================================================
# subfunction of display 360. Display 6 views rotating over posCam.
#==============================================================================
def display360__(t, posCam, posEye, dirCam, offscreen, exportRez, kwargs):
    import KCore.Vector as Vector

    # resolution for the square view images
    locRez = exportRez.split('x')[1]
    locRez = int(locRez)//2
    locRez = max(locRez, 100) # minimum 100 pixels
    locRez = min(locRez, 8192) # maximum 8192 pixels, generally the max texture size
    locRez = "%dx%d"%(locRez, locRez)

    # Compute all view vectors
    v1 = Vector.sub(posEye, posCam) # view vector
    vz = Vector.normalize(dirCam)
    v2 = Vector.cross(vz, v1) # second view vector
    n = Vector.norm(v1)
    v3 = Vector.mul(n, vz) # third view vector

    lkwargs = kwargs.copy()
    fov = 90.

    exportRoot = kwargs.get('export', 'export.png')
    exportRoot = os.path.dirname(exportRoot)
    if exportRoot == '': exportRoot = '.'

    # right
    posCam0 = posCam
    posEye0 = Vector.sub(posCam, v2)
    dirCam0 = dirCam
    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = fov
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_right.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # left
    posCam0 = posCam
    posEye0 = Vector.add(posCam, v2)
    dirCam0 = dirCam
    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = fov
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_left.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # front
    posCam0 = posCam
    posEye0 = posEye
    dirCam0 = dirCam
    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = fov
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_front.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # back
    posCam0 = posCam
    posEye0 = Vector.sub(posCam, v1)
    dirCam0 = dirCam
    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = fov
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_back.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # top
    posEye0 = Vector.add(posCam, v3)
    dirCam0 = Vector.mul(-1, v1)
    lkwargs['posCam'] = posCam
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = fov
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_top.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # bottom
    posEye0 = Vector.sub(posCam, v3)
    dirCam0 = Vector.mul(+1, v1)
    lkwargs['posCam'] = posCam
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = fov
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_bottom.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    return None

# subfunction of display 360. Display the 6 views with rotating posCam (wrong stereo)
# requires that nothing is in the corners
def display360WS__(t, posCam, posEye, dirCam, offscreen, exportRez, stereoShift, kwargs):
    import KCore.Vector as Vector

    lkwargs = kwargs.copy()

    # resolution for the square view images
    locRez = exportRez.split('x')[1]
    locRez = int(locRez)//2
    locRez = max(locRez, 100) # minimum 100 pixels
    locRez = min(locRez, 8192) # maximum 8192 pixels, generally the max texture size
    locRez = "%dx%d"%(locRez, locRez)

    exportRoot = kwargs.get('export', 'export.png')
    exportRoot = os.path.dirname(exportRoot)
    if exportRoot == '': exportRoot = '.'

    # Compute all front view vectors
    v1 = Vector.sub(posEye, posCam) # view vector
    vz = Vector.normalize(dirCam) # third view vector
    v2 = Vector.cross(v1, vz) # second view vector
    v2 = Vector.normalize(v2)

    import Geom.PyTree as D
    import Transform.PyTree as T

    # front image
    theta = 0.
    point = D.point(v1)
    point = T.rotate(point, (0,0,0), vz, theta)
    v1p = C.getValue(point, 'GridCoordinates', 0)
    point = D.point(v2)
    point = T.rotate(point, (0,0,0), vz, theta)
    v2p = C.getValue(point, 'GridCoordinates', 0)
    dv = Vector.mul(stereoShift, v2p)

    posCam0 = Vector.add(posCam, dv)
    posEye0 = Vector.add(v1p, posCam)
    dirCam0 = dirCam

    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = 90.
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_front.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # top image
    point = D.point(v1p)
    point = T.rotate(point, (0,0,0), v2p, +90)
    v1z = C.getValue(point, 'GridCoordinates', 0)
    point = D.point(dirCam)
    point = T.rotate(point, (0,0,0), v2p, +90)
    v2z = C.getValue(point, 'GridCoordinates', 0)

    posEye0 = Vector.add(v1z, posCam)
    dirCam0 = v2z

    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = 90.
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_top.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # bottom image
    point = D.point(v1p)
    point = T.rotate(point, (0,0,0), v2p, -90)
    v1z = C.getValue(point, 'GridCoordinates', 0)
    point = D.point(dirCam)
    point = T.rotate(point, (0,0,0), v2p, -90)
    v2z = C.getValue(point, 'GridCoordinates', 0)

    posEye0 = Vector.add(v1z, posCam)
    dirCam0 = v2z

    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = 90.
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_bottom.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # right image
    theta = -90.
    point = D.point(v1)
    point = T.rotate(point, (0,0,0), vz, theta)
    v1p = C.getValue(point, 'GridCoordinates', 0)
    point = D.point(v2)
    point = T.rotate(point, (0,0,0), vz, theta)
    v2p = C.getValue(point, 'GridCoordinates', 0)
    dv = Vector.mul(stereoShift, v2p)

    posCam0 = Vector.add(posCam, dv)
    posEye0 = Vector.add(v1p, posCam)
    dirCam0 = dirCam

    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = 90.
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_right.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # left image
    theta = +90.
    point = D.point(v1)
    point = T.rotate(point, (0,0,0), vz, theta)
    v1p = C.getValue(point, 'GridCoordinates', 0)
    point = D.point(v2)
    point = T.rotate(point, (0,0,0), vz, theta)
    v2p = C.getValue(point, 'GridCoordinates', 0)
    dv = Vector.mul(stereoShift, v2p)

    posCam0 = Vector.add(posCam, dv)
    posEye0 = Vector.add(v1p, posCam)
    dirCam0 = dirCam

    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = 90.
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_left.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    # back image
    theta = +180.
    point = D.point(v1)
    point = T.rotate(point, (0,0,0), vz, theta)
    v1p = C.getValue(point, 'GridCoordinates', 0)
    point = D.point(v2)
    point = T.rotate(point, (0,0,0), vz, theta)
    v2p = C.getValue(point, 'GridCoordinates', 0)
    dv = Vector.mul(stereoShift, v2p)

    posCam0 = Vector.add(posCam, dv)
    posEye0 = Vector.add(v1p, posCam)
    dirCam0 = dirCam

    lkwargs['posCam'] = posCam0
    lkwargs['posEye'] = posEye0
    lkwargs['dirCam'] = dirCam0
    lkwargs['viewAngle'] = 90.
    lkwargs['exportResolution'] = locRez
    lkwargs['export'] = exportRoot+'/run/cube_back.png'
    display(t, **lkwargs)
    finalizeExport(offscreen)

    return None

# subfunction of display 360. Display the n views with rotating posCam
def display360ODS__(t, posCam, posEye, dirCam, offscreen, exportRez, stereoShift, kwargs):
    import Converter.Mpi as Cmpi
    import KCore.Vector as Vector

    lkwargs = kwargs.copy()

    exportRoot = kwargs.get('export', 'export.png')
    exportRoot = os.path.dirname(exportRoot)
    if exportRoot == '': exportRoot = '.'

    # number of images, 1 per pixel
    nangles = exportRez.split('x')[0]
    nangles = int(nangles)
    # fov of each image
    fov = 90.

    # locrez of each image
    locRez = exportRez.split('x')[1]
    locRez1 = 2
    locRez2 = int(locRez)
    locRez = "%dx%d"%(locRez1, locRez2)

    # Compute all front view vectors
    v1 = Vector.sub(posEye, posCam) # view vector
    vz = Vector.normalize(dirCam) # third view vector
    v2 = Vector.cross(v1, vz) # second view vector
    v2 = Vector.normalize(v2)

    import Geom.PyTree as D
    import Transform.PyTree as T

    # start from -pi to pi and rotate left
    for i in range(nangles):

        # simple parallel hack
        #if i%Cmpi.size != Cmpi.rank: continue

        theta = 180. - i*360./nangles

        point = D.point(v1)
        point = T.rotate(point, (0,0,0), vz, theta)
        v1p = C.getValue(point, 'GridCoordinates', 0)
        point = D.point(v2)
        point = T.rotate(point, (0,0,0), vz, theta)
        v2p = C.getValue(point, 'GridCoordinates', 0)
        dv = Vector.mul(stereoShift, v2p)

        # front image
        posCam0 = Vector.add(posCam, dv)
        posEye0 = Vector.add(v1p, posCam)
        dirCam0 = dirCam
        print('front %d / %d'%(i,nangles))

        lkwargs['posCam'] = posCam0
        lkwargs['posEye'] = posEye0
        lkwargs['dirCam'] = dirCam0
        lkwargs['viewAngle'] = fov
        lkwargs['exportResolution'] = locRez
        lkwargs['export'] = exportRoot+'/run/front_%05d.png'%i
        display(t, **lkwargs)
        finalizeExport(offscreen)

        # top image
        point = D.point(v1p)
        point = T.rotate(point, (0,0,0), v2p, +90)
        v1z = C.getValue(point, 'GridCoordinates', 0)
        point = D.point(dirCam)
        point = T.rotate(point, (0,0,0), v2p, +90)
        v2z = C.getValue(point, 'GridCoordinates', 0)

        posEye0 = Vector.add(v1z, posCam)
        dirCam0 = v2z
        print('top %d / %d'%(i,nangles))

        lkwargs['posCam'] = posCam0
        lkwargs['posEye'] = posEye0
        lkwargs['dirCam'] = dirCam0
        lkwargs['viewAngle'] = fov
        lkwargs['exportResolution'] = locRez
        lkwargs['export'] = exportRoot+'/run/top_%05d.png'%i
        display(t, **lkwargs)
        finalizeExport(offscreen)

        # bottom image
        point = D.point(v1p)
        point = T.rotate(point, (0,0,0), v2p, -90)
        v1z = C.getValue(point, 'GridCoordinates', 0)
        point = D.point(dirCam)
        point = T.rotate(point, (0,0,0), v2p, -90)
        v2z = C.getValue(point, 'GridCoordinates', 0)

        posEye0 = Vector.add(v1z, posCam)
        dirCam0 = v2z
        print('bot %d / %d'%(i,nangles))

        lkwargs['posCam'] = posCam0
        lkwargs['posEye'] = posEye0
        lkwargs['dirCam'] = dirCam0
        lkwargs['viewAngle'] = fov
        lkwargs['exportResolution'] = locRez
        lkwargs['export'] = exportRoot+'/run/bottom_%05d.png'%i
        display(t, **lkwargs)
        finalizeExport(offscreen)

    Cmpi.barrier() # wait for completion
    return None

# subfunction of display 360. Display the n views with rotating posCam
# only for osmesa
def display360ODS2__(t, posCam, posEye, dirCam, offscreen, exportRez, stereoShift, kwargs):

    import KCore.Vector as Vector
    lkwargs = kwargs.copy()

    posCamList = []; posEyeList = []; dirCamList = []

    # number of images, 1 per pixel
    nangles = exportRez.split('x')[0]
    nangles = int(nangles)
    # fov of each image
    fov = 90.

    # locrez of each image
    locRez = exportRez.split('x')[1]
    locRez1 = 2
    locRez2 = int(locRez)
    locRez = "%dx%d"%(locRez1, locRez2)

    lkwargs['viewAngle'] = fov
    lkwargs['exportResolution'] = locRez

    # Compute all front view vectors
    v1 = Vector.sub(posEye, posCam) # view vector
    vz = Vector.normalize(dirCam) # third view vector
    v2 = Vector.cross(v1, vz) # second view vector
    v2 = Vector.normalize(v2)

    import Geom.PyTree as D
    import Transform.PyTree as T

    # start from -pi to pi and rotate left
    for i in range(nangles):

        theta = 180. - i*360./nangles

        point = D.point(v1)
        point = T.rotate(point, (0,0,0), vz, theta)
        v1p = C.getValue(point, 'GridCoordinates', 0)
        point = D.point(v2)
        point = T.rotate(point, (0,0,0), vz, theta)
        v2p = C.getValue(point, 'GridCoordinates', 0)
        dv = Vector.mul(stereoShift, v2p)

        # front image
        posCam0 = Vector.add(posCam, dv)
        posEye0 = Vector.add(v1p, posCam)
        dirCam0 = dirCam
        #print('front %d / %d'%(i,nangles))

        posCamList += posCam0
        posEyeList += posEye0
        dirCamList += dirCam0

        # top image
        point = D.point(v1p)
        point = T.rotate(point, (0,0,0), v2p, +90)
        v1z = C.getValue(point, 'GridCoordinates', 0)
        point = D.point(dirCam)
        point = T.rotate(point, (0,0,0), v2p, +90)
        v2z = C.getValue(point, 'GridCoordinates', 0)

        posEye0 = Vector.add(v1z, posCam)
        dirCam0 = v2z
        #print('top %d / %d'%(i,nangles))

        posCamList += posCam0
        posEyeList += posEye0
        dirCamList += dirCam0

        # bottom image
        point = D.point(v1p)
        point = T.rotate(point, (0,0,0), v2p, -90)
        v1z = C.getValue(point, 'GridCoordinates', 0)
        point = D.point(dirCam)
        point = T.rotate(point, (0,0,0), v2p, -90)
        v2z = C.getValue(point, 'GridCoordinates', 0)

        posEye0 = Vector.add(v1z, posCam)
        dirCam0 = v2z
        #print('bot %d / %d'%(i,nangles))

        posCamList += posCam0
        posEyeList += posEye0
        dirCamList += dirCam0

        lkwargs['posCamList'] = posCamList
        lkwargs['posEyeList'] = posEyeList
        lkwargs['dirCamList'] = dirCamList

    display(t, **lkwargs)

    return None

#==============================================================================
# display360 (offscreen=1, 2 or 7)
# type360=0 (360 degres), =1 (180 degres)
#==============================================================================
def display360(t, type360=0, **kwargs):
    """Display for 360 images."""
    import KCore.Vector as Vector
    import Converter.Mpi as Cmpi
    posCam = kwargs.get("posCam", (0,0,0))
    posEye = kwargs.get("posEye", (1,0,0))
    dirCam = kwargs.get("dirCam", (0,0,1))

    # get export resolution (final image), offscreen mode, stereo and stereo dist
    export = kwargs.get("export", "image360.png")
    exportRez = kwargs.get("exportResolution", "4096x2048")
    offscreen = kwargs.get("offscreen", 1)
    stereo = kwargs.get("stereo", 0)
    stereoDist = kwargs.get("stereoDist", 0.07) # stereoDist is in real world distance
    kwargs['stereo'] = 0 # force no anaglyph

    # orthogonalisation de v1 et dirCam si ils ne sont pas orthos
    v1 = Vector.sub(posEye, posCam) # view vector
    vz = Vector.normalize(dirCam)
    s = Vector.dot(v1, vz)
    v1 = Vector.sub(v1, Vector.mul(s, vz))
    posEye = Vector.add(posCam, v1)

    # display
    if stereo == 0:
        # display 6 views
        display360__(t, posCam, posEye, dirCam, offscreen, exportRez, kwargs)
        # Create the 360 image from cube images
        if Cmpi.rank == 0:
            panorama(export, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

    elif stereo == 1: # stereo ODS internal (only offscreen=1 or 7)

        if offscreen != 1 and offscreen != 7: raise ValueError("display360: stereo=1 only for osmesa."); return None

        export1 = export.rsplit('.', 1)
        if len(export1) == 2: export1 = export1[0]+'_1.'+export1[1]
        else: export1 = export+'_1'
        export2 = export.rsplit('.', 1)
        if len(export2) == 2: export2 = export2[0]+'_2.'+export2[1]
        else: export2 = export+'_2'

        # right eye
        #stereoDist = 0. # forced to 0 for debug
        kwargs['export'] = export1
        display360ODS2__(t, posCam, posEye, dirCam, offscreen, exportRez, stereoDist/2., kwargs)

        # left eye
        kwargs['export'] = export2
        display360ODS2__(t, posCam, posEye, dirCam, offscreen, exportRez, -stereoDist/2., kwargs)

        # stitch
        if Cmpi.rank == 0:
            panoramaStereo(export, export1, export2, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

    elif stereo == 2: # stereo = 2 (wrong stereo)
        export1 = export.rsplit('.', 1)
        if len(export1) == 2: export1 = export1[0]+'_1.'+export1[1]
        else: export1 = export+'_1'
        export2 = export.rsplit('.', 1)
        if len(export2) == 2: export2 = export2[0]+'_2.'+export2[1]
        else: export2 = export+'_2'

        # right eye
        display360WS__(t, posCam, posEye, dirCam, offscreen, exportRez, stereoDist/2., kwargs)
        if Cmpi.rank == 0:
            panorama(export1, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

        # left eye
        display360WS__(t, posCam, posEye, dirCam, offscreen, exportRez, -stereoDist/2., kwargs)
        if Cmpi.rank == 0:
            panorama(export2, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

        # stitch
        if Cmpi.rank == 0:
            panoramaStereo(export, export1, export2, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

    elif stereo == 3: # stereo = 3 (wrong stereo 2)
        export1 = export.rsplit('.', 1)
        if len(export1) == 2: export1 = export1[0]+'_1.'+export1[1]
        else: export1 = export+'_1'
        export2 = export.rsplit('.', 1)
        if len(export2) == 2: export2 = export2[0]+'_2.'+export2[1]
        else: export2 = export+'_2'

        v1 = Vector.sub(posEye, posCam) # view vector
        vz = Vector.normalize(dirCam) # third view vector
        v2 = Vector.cross(v1, vz) # second view vector
        v2 = Vector.normalize(v2)
        dv = Vector.mul(stereoDist/2, v2)

        # right eye
        posCam0 = Vector.add(posCam, dv)
        display360__(t, posCam0, posEye, dirCam, offscreen, exportRez, kwargs)
        if Cmpi.rank == 0:
            panorama(export1, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

        # left eye
        posCam0 = Vector.sub(posCam, dv)
        display360__(t, posCam0, posEye, dirCam, offscreen, exportRez, kwargs)
        if Cmpi.rank == 0:
            panorama(export2, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

        # stitch
        if Cmpi.rank == 0:
            panoramaStereo(export, export1, export2, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

    elif stereo == 4: # stereo ODS external

        export1 = export.rsplit('.', 1)
        if len(export1) == 2: export1 = export1[0]+'_1.'+export1[1]
        else: export1 = export+'_1'
        export2 = export.rsplit('.', 1)
        if len(export2) == 2: export2 = export2[0]+'_2.'+export2[1]
        else: export2 = export+'_2'

        # right eye
        #stereoDist = 0. # forced to 0 for debug
        kwargs['export'] = export1
        display360ODS__(t, posCam, posEye, dirCam, offscreen, exportRez, stereoDist/2., kwargs)
        if Cmpi.rank == 0:
            panoramaODS(export1, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

        # left eye
        kwargs['export'] = export2
        display360ODS__(t, posCam, posEye, dirCam, offscreen, exportRez, -stereoDist/2., kwargs)
        if Cmpi.rank == 0:
            panoramaODS(export2, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

        # stitch
        if Cmpi.rank == 0:
            panoramaStereo(export, export1, export2, exportRez, type360=type360)
        Cmpi.barrier() # wait for completion

    return None

# assemble 6 cube images en une image panoramique
# type360=0 -> 360, type360=1 -> 180
def panorama(export, exportRez, type360=0):
    res = exportRez.split('x')
    if type360 == 0: resx = int(res[0]); resy = int(res[1])
    else: resx = int(res[1]); resy = int(res[1])
    import Generator.PyTree as G
    import CPlot.cplot
    exportRoot = os.path.dirname(export)
    if exportRoot == '': exportRoot = '.'

    a1 = C.convertFile2PyTree(exportRoot+'/run/cube_left.png')
    a1 = C.getFields('nodes', a1, api=3)[0]
    a2 = C.convertFile2PyTree(exportRoot+'/run/cube_right.png')
    a2 = C.getFields('nodes', a2, api=3)[0]
    a3 = C.convertFile2PyTree(exportRoot+'/run/cube_bottom.png')
    a3 = C.getFields('nodes', a3, api=3)[0]
    a4 = C.convertFile2PyTree(exportRoot+'/run/cube_top.png')
    a4 = C.getFields('nodes', a4, api=3)[0]
    a5 = C.convertFile2PyTree(exportRoot+'/run/cube_back.png')
    a5 = C.getFields('nodes', a5, api=3)[0]
    a6 = C.convertFile2PyTree(exportRoot+'/run/cube_front.png')
    a6 = C.getFields('nodes', a6, api=3)[0]
    a7 = G.cart((0,0,0), (1,1,1), (resx, resy,1))
    C._addVars(a7, ['r','g','b','a'])
    a7f = C.getFields('nodes', a7, api=3)[0]
    CPlot.cplot.panorama(a1, a2, a3, a4, a5, a6, a7f, type360)
    C.convertPyTree2File(a7, export)
    return a7

# assemble 2 images panoramiques en une image panoramique stereo
# export1: right, export2: left
def panoramaStereo(export, export1, export2, exportRez, type360=0):
    import Generator.PyTree as G
    # assemble 2 panoramic images in a single stereo image
    a1 = C.convertFile2PyTree(export1)
    a2 = C.convertFile2PyTree(export2)
    a1 = Internal.getZones(a1)[0]
    a2 = Internal.getZones(a2)[0]
    a1[0] = "right"; a2[0] = "left"
    locRez = exportRez.split('x')
    if type360 == 0: # 360
        ni = int(locRez[0]); nj = int(locRez[1])
        a = G.cart((0,0,0), (1,1,1), (ni,2*nj,1))
        C._addVars(a, ['r','g','b','a'])
        for v in ['r','g','b','a']:
            pr = Internal.getNodeFromName2(a, v)[1]
            pr1 = Internal.getNodeFromName2(a1, v)
            pr2 = Internal.getNodeFromName2(a2, v)
            if pr1 is not None and pr2 is not None:
                pr1 = pr1[1]; pr2 = pr2[1]
                pr[0:ni,0:nj] = pr1[0:ni, 0:nj]
                pr[0:ni,nj:2*nj] = pr2[0:ni, 0:nj]
            else:
                if v == 'a': pr[0:ni, 0:2*nj] = 255
                else: pr[0:ni, 0:2*nj] = 0
    else: # 180
        ni = int(locRez[1]); nj = int(locRez[1])
        a = G.cart((0,0,0), (1,1,1), (2*ni,nj,1))
        C._addVars(a, ['r','g','b','a'])
        for v in ['r','g','b','a']:
            pr = Internal.getNodeFromName2(a, v)[1]
            pr1 = Internal.getNodeFromName2(a1, v)
            pr2 = Internal.getNodeFromName2(a2, v)
            if pr1 is not None and pr2 is not None:
                pr1 = pr1[1]; pr2 = pr2[1]
                pr[0:ni,0:nj] = pr1[0:ni,0:nj]
                pr[ni:2*ni,0:nj] = pr2[0:ni,0:nj]
            else:
                if v == 'a': pr[0:2*ni,0:nj] = 255
                else: pr[0:2*ni,0:nj] = 0

    C.convertPyTree2File(a, export) # finale

# assemble n cube images en une image panoramique
# type360=0 -> 360, type360=1 -> 180
def panoramaODS(export, exportRez, type360=0):

    res = exportRez.split('x')
    if type360 == 0: resx = int(res[0]); resy = int(res[1])
    else: resx = int(res[1]); resy = int(res[1])
    import Generator.PyTree as G
    import CPlot.cplot

    exportRoot = os.path.dirname(export)
    if exportRoot == '': exportRoot = '.'

    nangles = exportRez.split('x')[0]
    nangles = int(nangles)

    front = []; top = []; bottom = []
    for i in range(nangles):
        a1 = C.convertFile2PyTree(exportRoot+'/run/front_%05d.png'%i)
        a1 = C.getFields('nodes', a1, api=3)[0]
        front.append(a1)
        a1 = C.convertFile2PyTree(exportRoot+'/run/top_%05d.png'%i)
        a1 = C.getFields('nodes', a1, api=3)[0]
        top.append(a1)
        a1 = C.convertFile2PyTree(exportRoot+'/run/bottom_%05d.png'%i)
        a1 = C.getFields('nodes', a1, api=3)[0]
        bottom.append(a1)

    a7 = G.cart((0,0,0), (1,1,1), (resx, resy,1))
    C._addVars(a7, ['r','g','b','a'])
    a7f = C.getFields('nodes', a7, api=3)[0]
    CPlot.cplot.panoramaODS(front, top, bottom, a7f, type360)
    C.convertPyTree2File(a7, export)

# blur an image
def _blur(t, blurSigma=0.8):
    """Blur an image."""
    arrays = C.getAllFields(t, 'nodes', api=3)
    CPlot.blur(arrays, blurSigma)
    return None