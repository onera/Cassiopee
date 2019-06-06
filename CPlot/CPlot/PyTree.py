"""Plotter functions for pyTrees.
"""
import numpy
import CPlot
try: import cplot
except: ImportError("CPlot: is partially installed (no display).")

try: range = xrange
except: pass

# Separateur intra-nom
SEP1 = '/'

__version__ = CPlot.__version__
# Delay between display
__timeStep__ = CPlot.__timeStep__
# CPlot display fields localized in __LOCATION__ ('nodes' or 'centers')
__LOCATION__ = 'nodes'

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
            vectorScale=-1.,
            vectorDensity=-1.,
            vectorNormalize=-1,
            vectorShowSurface=-1,
            vectorShape=-1,
            vectorProjection=-1,
            colormap=-1,
            niso=25,
            isoEdges=-0.5,
            isoScales=[],
            win=(-1,-1),
            posCam=(-999,-999,-999),
            posEye=(-999,-999,-999),
            dirCam=(-999,-999,-999),
            viewAngle=-1.,
            bgColor=-1,
            shadow=-1,
            dof=-1,
            stereo=-1,
            stereoDist=-1.,
            export="None",
            exportResolution="None",
            location='unchanged',
            offscreen=0):
    """Display pyTrees.
    Usage: display(t)"""
    global __LOCATION__
    if location != 'unchanged': __LOCATION__ = location
    if __LOCATION__ == 'nodes':
        t = C.center2Node(t, Internal.__FlowSolutionCenters__)
    else: t = C.node2Center(t)
    zoneNames = C.getZoneNames(t)
    renderTags = getRenderTags(t)
    arrays = C.getAllFields(t, 'nodes', api=2)
    CPlot.display(arrays, dim, mode, scalarField, vectorField1,
                  vectorField2, vectorField3, displayBB, displayInfo,
                  displayIsoLegend, meshStyle,
                  solidStyle, scalarStyle, vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                  vectorShowSurface, vectorShape, vectorProjection,
                  colormap, niso, isoEdges, isoScales, win,
                  posCam, posEye, dirCam, viewAngle,
                  bgColor, shadow, dof, stereo, stereoDist,
                  export, exportResolution, 
                  zoneNames, renderTags, offscreen)

#==============================================================================
def render():
    """Force render.
    Usage: render()"""
    CPlot.render()

#==============================================================================
def delete(list):
    """Delete zones from plotter.
    Usage: delete([i1,i2,...])"""
    CPlot.delete(list)
    
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
    array = C.getAllFields(zone, 'nodes', api=2)[0]
    
    cplot.add(array, (nzs, nzu), zoneName, renderTag)

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
    array = C.getAllFields(zone, 'nodes', api=2)[0]
   
    cplot.replace(array, (nzs, nzu, oldType), zoneName, renderTag)
    
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
        arrays = C.getAllFields(t, 'nodes') # as pyZone
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
             shadow=-1,
             dof=-1, dofPower=-1,
             gamma=-1,
             sobelThreshold=-1,
             ghostifyDeactivatedZones=-1,
             edgifyActivatedZones=-1,
             edgifyDeactivatedZones=-1,
             export="None",
             exportResolution="None",
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
             materials=None, bumpMaps=None):
    """Set CPlot state.
    Usage: setState(posCam=(12,0,0))"""
    CPlot.setState(dim, mode, scalarField, vectorField1, vectorField2,
                   vectorField3, displayBB, displayInfo, displayIsoLegend,
                   meshStyle, solidStyle, scalarStyle,
                   vectorStyle, vectorScale, vectorDensity, vectorNormalize, 
                   vectorShowSurface, vectorShape, vectorProjection, colormap,
                   niso, isoEdges, isoScales, win,
                   posCam, posEye, dirCam, viewAngle, lightOffset,
                   bgColor, shadow, dof, dofPower, gamma, sobelThreshold,
                   ghostifyDeactivatedZones, edgifyActivatedZones,
                   edgifyDeactivatedZones,
                   export, exportResolution, continuousExport,
                   envmap, message,
                   stereo, stereoDist, cursor, gridSize, timer, selectionStyle,
                   activateShortCuts, billBoards, billBoardSize, 
                   materials, bumpMaps)
    
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

def setSelectedZones(list):
    """Set selected zones.
    Usage: setSelectedZones([(0,1),(1,1),...])"""
    CPlot.setSelectedZones(list)

def unselectAllZones():
    """Unselect all zones.
    Usage: unselectAllZones()"""
    CPlot.unselectAllZones()

def setActiveZones(list):
    """Set active (displayed) zones.
    Usage: setActiveZones([(0,1)])"""
    CPlot.setActiveZones(list)

def setZoneNames(list):
    """Set zone names.
    Usage: setZoneNames([(0,'myZone')])"""
    CPlot.setZoneNames(list)

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
        posCams = C.getAllFields(posCams, 'nodes')
    if posEyes is not None and isinstance(posEyes[0], str): # zone
        posEyes = C.getAllFields(posEyes, 'nodes')
    if dirCams is not None and isinstance(dirCams[0], str): # zone
        dirCams = C.getAllFields(dirCams, 'nodes')
    CPlot.moveCamera(posCams, posEyes, dirCams, moveEye, N, speed, pos)

def travelRight(xr=0.1, N=100):
    CPlot.travelRight(xr, N)

def travelLeft(xr=0.1, N=100):
    CPlot.travelLeft(xr, N)

def travelUp(xr=0.1, N=100):
    CPlot.travelUp(xr, N)

def travelDown(xr=0.1, N=100):
    CPlot.travelDown(xr, N)

def travelIn(xr=0.1, N=100):
    CPlot.travelIn(xr, N)

def travelOut(xr=0.1, N=100):
    CPlot.travelOut(xr, N)

#==============================================================================
# -- Numbering functions --
#==============================================================================
# Retourne la numerotation (base, zone) dans l'arbre
# a partir de la numerotation zone de CPlot
#==============================================================================
def updateCPlotNumbering(t):
    Nz = []; Nb = []
    nodes = Internal.getZones(t)
    nstr = 0
    for z in nodes:
        Nz.append(1); Nb.append(1)
        type = Internal.getZoneType(z)
        if type == 1: nstr += 1

    bases = t[2][1:]
    str = 0; unstr = 0
    nb = 0
    for b in bases:
        nodes = b[2]
        nz = 0
        for z in nodes:
            if z[3] == 'Zone_t':
                type = Internal.getZoneType(z)
                if type == 1:
                    Nz[str] = nz
                    Nb[str] = nb; str += 1
                else:
                    Nz[unstr+nstr] = nz
                    Nb[unstr+nstr] = nb; unstr += 1
            nz += 1
        nb += 1
    return (Nb, Nz)

#==============================================================================
# Retourne la numerotation CPlot de la zone baseName/zoneName de t
# IN: t: un noeud pyTree (et pas autre chose)
# IN: baseName: le nom de la base
# IN: zoneName: le nom de la zone
#==============================================================================
def getCPlotNumber(t, baseName, zoneName):
    base = Internal.getNodeFromName1(t, baseName) # must be unique
    z = Internal.getNodeFromName1(base, zoneName) # must be unique
    type = Internal.getZoneType(z)
    bases = Internal.getBases(t)
   
    if type == 1:
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
    else:
        nstr = 0
        for b in bases:
            nodes = Internal.getNodesFromType1(b, 'Zone_t')
            for z in nodes:
                type = Internal.getZoneType(z)
                if type == 1: nstr += 1
        unstr = nstr
        for b in bases:
            nodes = Internal.getNodesFromType1(b, 'Zone_t')
            if (b[0] != baseName):
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
# Un tag est de la forme color:material:blending:meshOverlay:shaderParameters,
# si renderInfo n'est pas present, on retourne None:None:None:None:None.
#==============================================================================
def getRenderTags(t):
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
    if ri is None: renderTags.append('None:None:None:None:None:None')
    else:
        rc = Internal.getNodeFromName1(ri, 'Color')
        if rc is None: color = 'None'
        else:
            v = Internal.getValue(rc)
            color = v
        rm = Internal.getNodeFromName1(ri, 'Material')
        if rm is None: material = 'None'
        else:
            v = Internal.getValue(rm)
            material = v
        rm = Internal.getNodeFromName1(ri, 'Blending')
        if rm is None: blending = 'None'
        else:
            v = Internal.getValue(rm)
            blending = str(v)
        rm = Internal.getNodeFromName1(ri, 'MeshOverlay')
        if rm is None: meshOverlay = 'None'
        else:
            v = Internal.getValue(rm)
            meshOverlay = str(v)
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
        renderTags.append(color+':'+material+':'+blending+':'+meshOverlay+':'+shaderParameters)
    return renderTags

# -- addRender2Zone
def addRender2Zone(t, material=None, color=None, blending=None,
                   meshOverlay=None, shaderParameters=None):
  """Add a renderInfo node to a zone node.
  Usage: addRender2Zone(zone, renderInfo)"""
  tp = Internal.copyRef(t)
  _addRender2Zone(tp, material, color, blending,
                  meshOverlay, shaderParameters)
  return tp

def _addRender2Zone(a, material=None, color=None, blending=None,
                    meshOverlay=None, shaderParameters=None):
  """Add a renderInfo node to a zone node.
  Usage: addRender2Zone(zone, renderInfo)"""
  if material is None and color is None and blending is None and meshOverlay is None and shaderParameters is None: return None
  
  zones = Internal.getZones(a)
  if zones == []: raise TypeError("addRender2Zone: zone must be a zone node.")
  for z in zones:

    ri = Internal.createUniqueChild(z, '.RenderInfo', 'UserDefinedData_t')

    if material is not None:
      rt = Internal.createUniqueChild(ri, 'Material', 'DataArray_t', value=material)

    if color is not None:
      rt = Internal.createUniqueChild(ri, 'Color', 'DataArray_t', value=color)

    if blending is not None:
      rt = Internal.createUniqueChild(ri, 'Blending', 'DataArray_t', value=blending)

    if meshOverlay is not None:
      rt = Internal.createUniqueChild(ri, 'MeshOverlay', 'DataArray_t', value=meshOverlay)

    if shaderParameters is not None:
      rt = Internal.createUniqueChild(ri, 'ShaderParameters', 'DataArray_t', value=shaderParameters)
  return None

# -- addRender2PyTree
def addRender2PyTree(t, slot=0, posCam=None, posEye=None, dirCam=None,
                     mode=None, scalarField=None, niso=None, isoScales=None,
                     isoEdges=None, isoLight=None,
                     colormap=None, materials=None, bumpMaps=None, billBoards=None):
  """Add a renderInfo node to a tree.
  Usage: addRender2PyTree(t, slot, renderInfo)"""
  a = Internal.copyRef(t)
  _addRender2PyTree(a, slot, posCam, posEye, dirCam,
                    mode, scalarField, niso, isoScales,
                    isoEdges, isoLight, colormap, 
                    materials, bumpMaps, billBoards)
  return a

def _addRender2PyTree(a, slot=0, posCam=None, posEye=None, dirCam=None,
                     mode=None, scalarField=None, niso=None, isoScales=None,
                     isoEdges=None, isoLight=None,
                     colormap=None, materials=None, bumpMaps=None, billBoards=None):
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

  if isoScales is not None: # must be a list
    rt = Internal.getNodeFromName1(sl, 'isoScales')
    n = len(isoScales)
    
    if rt is None:
      v = numpy.empty((n), numpy.float64)
      for i in range(n): v[i] = float(isoScales[i])
      sl[2].append(['isoScales', v, [], 'DataArray_t'])
    elif rt[1] is not None:
      old = rt[1]
      dict = {}
      l = old.shape[0]
      for i in range(0,l,4): dict[old[i]] = [old[i+1], old[i+2], old[i+3]]
      for i in range(0,n,4): dict[isoScales[i]] = [isoScales[i+1], isoScales[i+2], isoScales[i+3]]
      k = dict.keys(); l = len(k)
      v = numpy.empty((4*l), numpy.float64)
      c = 0
      for i in k:
        v[c] = i; v[c+1] = dict[i][0]; v[c+2] = dict[i][1]; v[c+3] = dict[i][2]
        c += 4
      rt[1] = v

  if isoEdges is not None:
    rt = Internal.createUniqueChild(sl, 'isoEdges', 'DataArray_t', value=isoEdges)

  if isoLight is not None:
    rt = Internal.createUniqueChild(sl, 'isoLight', 'DataArray_t', value=isoLight)

  if colormap is not None:
    rt = Internal.createUniqueChild(sl, 'colormap', 'DataArray_t', value=colormap)

  # Under .RenderInfo
  if materials is not None:
    rt = Internal.createUniqueChild(ri, 'materials', 'UserDefinedData_t')
    prevValues = [Internal.getValue(i) for i in rt[2]]
    cnt = len(prevValues)
    li = prevValues
    for i in materials:
      if i not in prevValues: li.append(i)
    for c, f in enumerate(li):
        rt = Internal.createUniqueChild(rt, 'file%d'%c, 'DataArray_t', value=f)

  if bumpMaps is not None:
    rt = Internal.createUniqueChild(ri, 'bumpMaps', 'UserDefinedData_t')
    prevValues = [Internal.getValue(i) for i in rt[2]]
    cnt = len(prevValues)
    li = prevValues
    for i in bumpMaps:
      if i not in prevValues: li.append(i)
    for c, f in enumerate(li):
        rt = Internal.createUniqueChild(rt, 'file%d'%c, 'DataArray_t', value=f)

  if billBoards is not None:
    rt = Internal.createUniqueChild(ri, 'billBoards', 'UserDefinedData_t')
    prevValues = [Internal.getValue(i) for i in rt[2]]
    cnt = len(prevValues)
    li = prevValues
    for i in billBoards:
      if i not in prevValues: li.append(i)
    for c, f in enumerate(li):
        re = Internal.createUniqueChild(rt, 'file%d'%(c+cnt), 'DataArray_t', value=f)
  return None

#==============================================================================
# loadView form slot
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
        vars = C.getVarNames(t)[0]
        ifield = 0
        for i in vars:
            if i == field: break
            if (i != 'CoordinateX' and i != 'CoordinateY' and i != 'CoordinateZ'):
                ifield += 1
        CPlot.setState(scalarField=ifield)
    pos = Internal.getNodeFromName1(slot, 'niso')
    if pos is not None:
        niso = Internal.getValue(pos); CPlot.setState(niso=int(niso))
    pos = Internal.getNodeFromName1(slot, 'isoEdges')
    if pos is not None:
        edge = Internal.getValue(pos); CPlot.setState(isoEdges=edge)
    pos = Internal.getNodeFromName1(slot, 'isoLight')
    if pos is not None:
        light = Internal.getValue(pos)
    else: light = 1
    pos = Internal.getNodeFromName1(slot, 'colormap')
    if pos is not None:
        colormap = Internal.getValue(pos)
    else: colormap = 'Blue2Red'
    style = 0
    if colormap == 'Blue2Red': style = 0
    elif colormap == 'Green2Red': style = 2
    elif colormap == 'Black2White': style = 4
    elif colormap == 'White2Black': style = 6
    elif colormap == 'Diverging': style = 8
    if light == 1: style += 1
    CPlot.setState(colormap=style)
    pos = Internal.getNodeFromName1(slot, 'isoScales')
    if pos is not None:
        CPlot.setState(isoScales=pos[1].tolist())
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
def loadImageFiles(t):
    """Load image files (texture, billboards, bumpmaps) in CPlot."""
    if t == []: return
    renderInfo = Internal.getNodeFromName1(t, '.RenderInfo')
    if renderInfo is None: return None
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
    return None
