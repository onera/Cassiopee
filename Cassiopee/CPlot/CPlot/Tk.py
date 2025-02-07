"""Tkinter interface for CPlot"""
try: import tkinter as TK
except ImportError: import Tkinter as TK
from . import Ttk as TTK
import Converter.PyTree as C
import Converter
import Converter.Internal as Internal
from . import CPlot as CP
import Transform
import Post
from . import PyTree as CPlot
from . import Panels
from . import iconics
import os, os.path
from sys import version_info
import re, fnmatch

# Set this to suppress CPlot firewalls (fixNGon, breakConnect) and enables direct v4
FIREWALL = True

#==============================================================================
# Variables globales partagees entre toutes les apps tk
#==============================================================================
# Fichier de donnees lu (en cours), nom du fichier et/ou Filter handle et/ou CAD hook
FILE = ''
HANDLE = None
CADHOOK = None
# Fichier pour les exports
EXPORTFILE = ''
# Update du menu Apps (semble poser pb sur certains linux)
UPDATEAPPSMENU = False
# pyTree contenant les donnees affichee
t = []
# pyTree contenant l'arbre precedent (pour le undo)
tp = []
# pyTree transcient (pour les donnees temporaires)
dt = []
# Numerotation globale entre CPlot et le pyTree t
Nb = []; Nz = []
# Widget des sorties textes (commun a toutes les apps)
TXT = None
# TkTree module
TKTREE = None
# TkPlotXY module
TKPLOTXY = None
# Other modules dictionary
TKMODULES = {} ; TKMODULEFRAMES = {}
# Pref dictionary (contains only strings)
PREFS = {}
# BUSY=True quand l'app est occupee
__BUSY__ = False
# localisation des donnees affichees par CPlot
__LOCATION__ = 'nodes'
# Undo=True si le undo est actif
__UNDO__ = True
# - CPlot filters -
# Points affiches
# 0: exclusive - tous les elements en non structure
# 1: all pts - mais affiche seulement les faces externes
# 2 et plus: one over n + faces externes affichees
__ONEOVERN__ = 1
# Champs transmis a CPlot
__FIELD__ = '__all__'
# Status de l'arbre (MAIN=1 ou <=0 - voir ci-dessous)
__MAINTREE__ = 1
# Les differents status pour l'arbre
MAIN=1; DEFINEDBC=-1; TIME=-2; SLICE=-2; CELLN=-3; MESHQUAL=-4;
UNDEFINEDBC=-5; IJK=-6; MONITOR=-7;UNDEFINEDIBC=-8;

# Sauvegarde des zones actives de main (avant de basculer sur un arbre
# temporaire)
__MAINACTIVEZONES__ = []

# Garde la trace de certains widgets
WIDGETS = {}

#==============================================================================
# Fonts
#==============================================================================
def defineFonts(FONTTYPE, FONTSIZE):
    GENERALFONT = (FONTTYPE, FONTSIZE)
    FRAMEFONT = (FONTTYPE, FONTSIZE, 'bold')
    LABELFONT = (FONTTYPE, FONTSIZE)
    MENUFONT = (FONTTYPE, FONTSIZE+1)
    BUTTONFONT = (FONTTYPE, FONTSIZE)
    TEXTFONT = (FONTTYPE, FONTSIZE+1)
    MSGDFONT = (FONTTYPE, FONTSIZE+2)
    FIXEDFONT = ('Courier', FONTSIZE, 'bold')
    return (GENERALFONT, FRAMEFONT, LABELFONT, MENUFONT, BUTTONFONT,
            TEXTFONT, MSGDFONT, FIXEDFONT)

FONTTYPE = 'Helvetica'
FONTSIZE = 10
(GENERALFONT, FRAMEFONT, LABELFONT, MENUFONT, BUTTONFONT, TEXTFONT, MSGDFONT, FIXEDFONT) = defineFonts(FONTTYPE, FONTSIZE)

#==============================================================================
# Frame style
#==============================================================================
FRAMESTYLE = TK.GROOVE

#==============================================================================
# Les types de fichiers lisibles par Converter
fileTypes=[('converter', '*.plt'),
           ('bin tecplot', '*.plt'),
           ('converter','*.tp'),
           ('fmt tecplot','*.tp'),
           ('converter','*.dat'),
           ('fmt tecplot','*.dat'),
           ('converter','*.v3d'),
           ('bin v3d','*.v3d'),
           ('converter','*.fv3d'),
           ('fmt v3d','*.fv3d'),
           ('converter','*.su2'),
           ('fmt su2','*.su2'),
           ('converter','*.cgns'),
           ('bin adf cgns','*.cgns'),
           ('converter','*.adf'),
           ('bin adf cgns','*.adf'),
           ('converter','*.hdf'),
           ('bin hdf cgns','*.hdf'),
           ('converter','*.hdf5'),
           ('bin hdf cgns','*.hdf5'),
           ('converter','*.mesh'),
           ('fmt mesh','*.mesh'),
           ('converter','*.msh'),
           ('fmt gmsh','*.msh'),
           ('converter','*.stl'),
           ('fmt stl','*.stl'),
           ('converter','*.bstl'),
           ('bin stl','*.bstl'),
           ('converter','*.iges'),
           ('fmt iges','*.iges'),
           ('converter','*.igs'),
           ('fmt iges','*.igs'),
           ('converter','*.step'),
           ('fmt step','*.step'),
           ('converter','*.stp'),
           ('fmt step','*.stp'),
           ('converter', '*.obj'),
           ('fmt obj', '*.obj'),
           ('converter', '*.3ds'),
           ('bin 3ds', '*.3ds'),
           ('converter', '*.ply'),
           ('bin ply', '*.ply'),
           ('converter', '*.pov'),
           ('fmt pov', '*.pov'),
           ('converter', '*.wav'),
           ('bin wav', '*.wav'),
           ('converter', '*.fig'),
           ('fmt fig', '*.fig'),
           ('converter', '*.svg'),
           ('fmt svg', '*.svg'),
           ('converter', '*.gts'),
           ('fmt gts', '*.gts'),
           ('converter', '*.png'),
           ('bin png', '*.png'),
           ('converter', '*.jpg'),
           ('bin jpg', '*.jpg'),
           ('converter', '*.jpeg'),
           ('bin jpg', '*.jpeg'),
           ('converter', '*.PLT'),
           ('bin tecplot', '*.PLT'),
           ('converter','*.TP'),
           ('fmt tecplot','*.TP'),
           ('converter','*.DAT'),
           ('fmt tecplot','*.DAT'),
           ('converter','*.V3D'),
           ('bin v3d','*.V3D'),
           ('converter','*.FV3D'),
           ('fmt v3d','*.FV3D'),
           ('converter','*.SU2'),
           ('fmt su2','*.SU2'),
           ('converter','*.CGNS'),
           ('bin adf cgns','*.CGNS'),
           ('converter','*.ADF'),
           ('bin adf cgns','*.ADF'),
           ('converter','*.HDF'),
           ('bin hdf cgns','*.HDF'),
           ('converter','*.HDF5'),
           ('bin hdf cgns','*.HDF5'),
           ('converter','*.MESH'),
           ('fmt mesh','*.MESH'),
           ('converter','*.STL'),
           ('bin stl','*.STL'),
           ('converter','*.BSTL'),
           ('bin stl','*.BSTL'),
           ('converter','*.IGES'),
           ('fmt iges','*.IGES'),
           ('converter','*.IGS'),
           ('fmt iges','*.IGS'),
           ('converter', '*.OBJ'),
           ('fmt obj', '*.OBJ'),
           ('converter', '*.3DS'),
           ('bin 3ds', '*.3DS'),
           ('converter', '*.PLY'),
           ('bin ply', '*.PLY'),
           ('converter', '*.POV'),
           ('fmt pov', '*.POV'),
           ('converter', '*.WAV'),
           ('bin wav', '*.WAV'),
           ('converter', '*.FIG'),
           ('fmt fig', '*.FIG'),
           ('converter', '*.SVG'),
           ('fmt svg', '*.SVG'),
           ('converter', '*.GTS'),
           ('fmt gts', '*.GTS'),
           ('converter', '*.PNG'),
           ('bin png', '*.PNG'),
           ('converter', '*.JPG'),
           ('bin jpg', '*.JPG'),
           ('converter', '*.JPEG'),
           ('bin jpg', '*.JPEG'),
           ('All files', '*')
           ]

#==============================================================================
# Retourne le module si il est loade, le load et cree l'App sinon
#==============================================================================
def getModule(app):
    global TKMODULES
    if app not in TKMODULES: return None
    if TKMODULES[app] is None:
        #try:
        module = __import__(app)
        TKMODULES[app] = module
        frame = TKMODULEFRAMES[app]
        module.createApp(frame)
        #except: print('Warning: module %s can not be loaded.'%app)
    return TKMODULES[app]

def openApp(app):
    module = getModule(app)
    if module is not None: module.showApp()

# 0: None, 2: watch, 3: cross
def setCursor(cursor, B=None, C=None, D=None):
    win = WIDGETS['masterWin']
    if cursor == 2: win.config(cursor="watch")
    elif cursor == 1: win.config(cursor="cross")
    elif cursor == 0: win.config(cursor="")
    win.update()
    if B is not None:
        if cursor == 2: B.config(cursor="watch")
        elif cursor == 1: B.config(cursor="cross")
        elif cursor == 0: B.config(cursor="")
        B.update()
    if C is not None:
        if cursor == 2: C.config(cursor="watch")
        elif cursor == 1: C.config(cursor="cross")
        elif cursor == 0: C.config(cursor="")
        C.update()
    if D is not None:
        if cursor == 2: D.config(cursor="watch")
        elif cursor == 1: D.config(cursor="cross")
        elif cursor == 0: D.config(cursor="")
        D.update()

#==============================================================================
# Extrait des arrays de a pour les envoyer au plotter
# - transforme la solution en centres en noeuds
# - prend 1 pt sur N si besoin
# - prend les elts exterieurs pour les Tetra et les Hexa
#==============================================================================
def buildCPlotArrays(a, topTree=[]):
    if CPlot.__LOCATION__ == 'nodes':
        if __FIELD__ == '__all__':
            a = C.center2Node(a, Internal.__FlowSolutionCenters__)
        else:
            v = __FIELD__.split(':')
            if len(v) == 2 and v[0] == 'centers':
                a = C.center2Node(a, __FIELD__)
    else: a = C.node2Center(a)

    if FIREWALL: api = 1
    else: api = 3

    if __FIELD__ == '__all__':
        arrays = C.getAllFields(a, 'nodes', api=api)
    elif __FIELD__ == '__none__':
        arrays = C.getFields(Internal.__GridCoordinates__, a, api=api)
    else:
        arrays = C.getFields(Internal.__GridCoordinates__, a, api=api)
        v = __FIELD__.split(':')
        if len(v) == 2: v = v[1]
        else: v = __FIELD__
        arrays2 = C.getField(v, a, api=api)
        for i, b in enumerate(arrays):
            if b != []: Converter._addVars([arrays[i], b])

    if __ONEOVERN__ > 1:
        for i, b in enumerate(arrays):
            if len(b) == 5:
                arrays[i] = Transform.oneovern(b, (__ONEOVERN__,__ONEOVERN__,__ONEOVERN__))

    # Transmet les maillages contenant les borders elts pour HEXA, TETRA,
    # PYRA, PENTA, NGON
    if __ONEOVERN__ > 0:
        for i, b in enumerate(arrays):
            if FIREWALL:
                if b[3] == 'TETRA' or b[3] == 'HEXA' or b[3] == 'PYRA' or b[3] == 'PENTA':
                    arrays[i] = Post.exteriorElts(b)
                if b[3] == 'NGON' and b[2][0,2] > 2:
                    arrays[i] = Post.exteriorElts(b)
    return arrays

#==============================================================================
# Surcharge de la fonction display de CPlot pour prendre en compte des
# filtres supplementaires (ONEOVERN, ...)
#==============================================================================
def display(t, dim=-1,
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
            viewAngle=-1.,
            bgColor=-1,
            backgroundFile="None",
            shadow=-1, lightOffset=(-999,-999),
            dof=-1, dofPower=-1, gamma=-1, toneMapping=-1,
            stereo=-1, stereoDist=-1., panorama=0,
            export='None', exportResolution='None',
            location='unchanged',
            mainTree=1):
    """Display pyTrees.
    Usage: display(t)"""
    global __MAINTREE__, __MAINACTIVEZONES__
    if mainTree <= 0 and __MAINTREE__ == 1:
        __MAINACTIVEZONES__ = CP.getActiveZones()
    if location != 'unchanged': CPlot.__LOCATION__ = location
    zoneNames = C.getZoneNames(t)
    renderTags = CPlot.getRenderTags(t)
    arrays = buildCPlotArrays(t)
    CP.display(arrays, dim, mode, scalarField, vectorField1, vectorField2,
               vectorField3, displayBB, displayInfo, displayIsoLegend,
               meshStyle, solidStyle, scalarStyle, vectorStyle, vectorScale, vectorDensity, vectorNormalize,
               vectorShowSurface, vectorShape, vectorProjection,
               colormap, colormapC1, colormapC2, colormapC3, colormapC,
               niso, isoEdges, isoScales,
               win, posCam, posEye, dirCam, viewAngle,
               bgColor, backgroundFile,
               shadow, lightOffset, dof, dofPower, gamma, toneMapping,
               stereo, stereoDist, panorama,
               export, exportResolution,
               zoneNames, renderTags)
    if mainTree == 1 and __MAINTREE__ <= 0:
        __MAINTREE__ = 1
        active = [(i,0) for i in range(len(arrays))]
        for i in __MAINACTIVEZONES__: active[i] = (i,1)
        CP.setActiveZones(active)
    __MAINTREE__ = mainTree

#==============================================================================
# Surcharge de la fonction replace de CPlot
# Les filtres sont ajoutes
#==============================================================================
def replace(t, nob, noz, zone):
    zoneName = t[2][nob][0]+Internal.SEP1+zone[0]
    renderTag = CPlot.getRenderTags__(zone, [])[0]

    oldType = Internal.getZoneType(t[2][nob][2][noz])
    t[2][nob][2][noz] = zone
    if CP.__slot__ is None: display(t); return
    (nzs, nzu) = CPlot.getNzs(t, zone)
    array = buildCPlotArrays(zone, topTree=t)[0]
    from . import cplot
    cplot.replace(array, (nzs, nzu, oldType), zoneName, renderTag)

#==============================================================================
# Surcharge de la fonction add de CPlot
# Les filtres sont ajoutes
#==============================================================================
def add(t, nob, noz, zone):
    zoneName = t[2][nob][0]+Internal.SEP1+zone[0]
    renderTag = CPlot.getRenderTags__(zone, [])[0]

    if noz == -1: noz = len(t[2][nob][2]) # insere a la fin
    t[2][nob][2].insert(noz, zone)
    if CP.__slot__ is None: display(t); return
    (nzs, nzu) = CPlot.getNzs(t, zone)
    array = buildCPlotArrays(zone, topTree=t)[0]
    from . import cplot
    cplot.add(array, (nzs, nzu), zoneName, renderTag)

#==============================================================================
def displayMainTree():
    if t == []: return
    display(t)

#==============================================================================
def viewDeactivatedZones():
    if t == []: return
    state = CPlot.getState('ghostifyDeactivatedZones')
    if state == 1: CPlot.setState(ghostifyDeactivatedZones=0)
    else: CPlot.setState(ghostifyDeactivatedZones=1)

#==============================================================================
# Montre la selection dans le tkTree
#==============================================================================
def showSelectionInTkTree(event=None):
    if t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        TXT.insert('START', 'Selection is empty.\n')
        TXT.insert('START', 'Error: ', 'Error'); return
    if __MAINTREE__ == 1:
        nz = nzs[0]
        nob = Nb[nz]+1
        noz = Nz[nz]
        baseName = t[2][nob][0]
        zoneName = t[2][nob][2][noz][0]
        TKTREE.WIDGETS['tree'].focusOnGivenZone(baseName, zoneName)
    elif __MAINTREE__ == DEFINEDBC: # id par baseName/ZoneName/bcName
        (Nob, Noz) = CPlot.updateCPlotNumbering(dt)
        nz = nzs[0]
        nob = Nob[nz]+1
        noz = Noz[nz]
        baseName = dt[2][nob][0]
        zoneName = dt[2][nob][2][noz][0]
        zoneName = zoneName.split(Internal.SEP1,2)
        if len(zoneName) == 3:
            baseName = zoneName[0]
            bcName = zoneName[2]
            zoneName = zoneName[1]
            TKTREE.WIDGETS['tree'].focusOnGivenZone(baseName, zoneName, bcName)
    elif __MAINTREE__ == UNDEFINEDBC: # id par baseName/zoneName/c
        (Nob, Noz) = CPlot.updateCPlotNumbering(dt)
        nz = nzs[0]
        nob = Nob[nz]+1
        noz = Noz[nz]
        baseName = dt[2][nob][0]
        zoneName = dt[2][nob][2][noz][0]
        # zoneName a ete etendu par appendBaseName2ZoneName
        zoneName = zoneName.split(Internal.SEP1,2)
        if len(zoneName) == 3:
            baseName = zoneName[0]
            zoneName = zoneName[1]
            TKTREE.WIDGETS['tree'].focusOnGivenZone(baseName, zoneName)
    elif __MAINTREE__ == SLICE or __MAINTREE__ == CELLN: # id by BaseName/ZoneName
        (Nob, Noz) = CPlot.updateCPlotNumbering(dt)
        nz = nzs[0]
        nob = Nob[nz]+1
        noz = Noz[nz]
        baseName = dt[2][nob][0]
        zoneName = dt[2][nob][2][noz][0]
        zoneName = zoneName.split(Internal.SEP1,1)
        if len(zoneName) == 2:
            baseName = zoneName[0]
            zoneName = zoneName[1]
            TKTREE.WIDGETS['tree'].focusOnGivenZone(baseName, zoneName)
    elif __MAINTREE__ == UNDEFINEDIBC: # id par baseName/zoneName/c
        (Nob, Noz) = CPlot.updateCPlotNumbering(dt)
        nz = nzs[0]
        nob = Nob[nz]+1
        noz = Noz[nz]
        baseName = dt[2][nob][0]
        zoneName = dt[2][nob][2][noz][0]
        # zoneName a ete etendu par appendBaseName2ZoneName
        zoneName = zoneName.split(Internal.SEP1,2)
        if len(zoneName) == 3:
            baseName = zoneName[0]
            zoneName = zoneName[1]
            TKTREE.WIDGETS['tree'].focusOnGivenZone(baseName, zoneName)

#==============================================================================
# Upgrade a pyTree to work flawlessly with tkCassiopee
# tkCassiopee requires: variables for all zones must be the same
# coordinates must exists, version node is the first node
#==============================================================================
def upgradeTree(t):
    Internal.autoSetContainers(t)
    Internal._correctPyTree(t, level=0) # version node
    #t = Internal.correctPyTree(t, level=9) # connectivity
    if FIREWALL: Internal._fixNGon(t) # suppressed in v4
    try:
        if C.isNamePresent(t, 'CoordinateX') <= 0: C._addVars(t, 'CoordinateX')
        if C.isNamePresent(t, 'CoordinateY') <= 0: C._addVars(t, 'CoordinateY')
        if C.isNamePresent(t, 'CoordinateZ') <= 0: C._addVars(t, 'CoordinateZ')
    except: pass
    return t

#==============================================================================
# tkFile dialog (multiple files) retourne:
# Sur windows, une string avec les noms encadres par {}
# Sous unix, retourne une liste de strings
# Sous unix, il faut parfois laisser tomber le premier de la liste
# qui est le fichier par defaut (initFile)
# Cette fonction rend toujours une liste de strings
#==============================================================================
def fixFileString__(files, initFile=None):
    try:
        import platform
        system = platform.uname()[0]
    except: system = 'unix'

    if version_info[0] == 2 and isinstance(files, unicode): # windows old bug (single unicode)
        import sys
        encoding = sys.getfilesystemencoding()
        # try to find { and }
        out = []
        while len(files) > 0:
            c = 0
            pos1 = files.find(u'{', c)
            if pos1 == -1: break
            c = pos1+1
            pos2 = files.find(u'}', c)
            if pos2 == -1: break
            c = pos2+1
            if pos2 > pos1: out.append(files[pos1+1:pos2])
            files = files[:pos1] + files[pos2+1:]

        # split les autres
        sp = files.split(u' ')
        for s in sp:
            s = s.encode(encoding)
            if s != ' ' and s != '': out.append(s)
    elif system == 'Windows': # doesnt return initfile
        if initFile != '' and initFile is not None:
            if len(files) == 0: out = [initFile]
            else: out = files
        else: out = files
    else: # unix
        # unix retourne aussi initFile en premier
        if initFile != '' and initFile is not None:
            if len(files) == 0: out = [initFile]
            else: out = files[1:]
        else: out = files
    if version_info[0] == 2:
        out = [o.encode('utf-8') for o in out] # Force utf-8
    return out

#==============================================================================
# Pour tkFileDialog (simple), fix l'encoding de la chaine si
# necessaire
#==============================================================================
def fixFileString2__(file):
    #if isinstance(file, unicode):
    #    import sys
    #    encoding = sys.getfilesystemencoding()
    #    s = file.encode(encoding)
    #    return s
    #else: return file
    if version_info[0] == 2: return file.encode('utf-8')
    else: return file

#==============================================================================
# Load a file par un dialog
# OUT: FILE: nom du fichier choisi
# OUT: t: l'arbre pyTree loade
# OUT: Nb; Nz: la correspondance de zones entre CPlot et le pyTree
#==============================================================================
def loadFile(event=None):
    global FILE; global t; global Nb; global Nz; global __FIELD__
    try: import tkinter.filedialog as tkFileDialog
    except: import tkFileDialog
    files = tkFileDialog.askopenfilenames(
        filetypes=fileTypes, initialfile=FILE, multiple=1)
    if files == '' or files is None or files == (): # user cancel
        return
    files = fixFileString__(files, FILE)
    try:
        FILE = files[0]
        t = []
        for file in files:
            t2 = C.convertFile2PyTree(file, density=1.)
            if t == []: t = t2
            else: t = C.mergeTrees(t, t2)
        t = upgradeTree(t)
        (Nb, Nz) = CPlot.updateCPlotNumbering(t); TKTREE.updateApp()
        if 'tkContainers' in TKMODULES: TKMODULES['tkContainers'].updateApp()
        if TKPLOTXY is not None: TKPLOTXY.updateApp()
        Panels.updateRenderPanel()
        fileName = os.path.split(FILE)[1]
        filePath = os.path.split(FILE)[0]
        TXT.insert('START', 'File %s read.\n'%fileName)
        changeWindowTitle(fileName, filePath)
        errors = Internal.checkPyTree(t, level=5)
        if errors == []: display(t)
        else:
            t = Internal.correctPyTree(t, level=5)
            display(t); Panels.displayErrors(errors, header='Checking pyTree')
    except:
        TXT.insert('START', 'Cannot read file '+os.path.split(files[0])[1]+'.\n')
        TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# Add a file par un dialog
# OUT: t: l'arbre pyTree loade
# OUT: Nb; Nz: la correspondance de zones entre CPlot et le pyTree
#==============================================================================
def addFile():
    global t; global Nb; global Nz
    try: import tkinter.filedialog as tkFileDialog
    except: import tkFileDialog
    files = tkFileDialog.askopenfilenames(
        filetypes=fileTypes, initialfile=FILE, multiple=1)
    if files == '' or files is None or files == (): # user cancel
        return
    files = fixFileString__(files, FILE)

    try:
        saveTree()
        for f in files:
            t2 = C.convertFile2PyTree(f, density=1.)
            # Fusion des bases de t et t2
            if t == []: t = t2
            else:
                #t = C.mergeTrees(t, t2)
                if Internal.getBases(t2) != []: # merge trees
                    t = C.mergeTrees(t, t2)
                else: # ajoute le noeud sous la selection de tkTree
                    node = TKTREE.getCurrentSelectedNode()
                    node[2] += t2[2][1:]
        t = upgradeTree(t)
        (Nb, Nz) = CPlot.updateCPlotNumbering(t); TKTREE.updateApp()
        if 'tkContainers' in TKMODULES: TKMODULES['tkContainers'].updateApp()
        if TKPLOTXY is not None: TKPLOTXY.updateApp()
        Panels.updateRenderPanel()
        TXT.insert('START', 'File '+os.path.split(files[0])[1]+' added.\n')
        errors = Internal.checkPyTree(t, level=5)
        if errors == []: display(t)
        else:
            t = Internal.correctPyTree(t, level=5)
            display(t); Panels.displayErrors(errors, header='Checking pyTree')
    except:
        TXT.insert('START', 'Can not add file '+os.path.split(files[0])[1]+'.\n')
        TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# Save a file par un dialog
# OUT: FILE: le fichier choisi
#==============================================================================
def saveFile():
    global FILE
    try: import tkinter.filedialog as tkFileDialog
    except ImportError: import tkFileDialog
    ret = tkFileDialog.asksaveasfilename(filetypes=fileTypes, initialfile=FILE)
    if ret == '' or ret is None or ret == (): # user cancel
        return
    try:
        FILE = fixFileString2__(ret)
        C.convertPyTree2File(t, FILE)
        fileName = os.path.split(FILE)[1]
        filePath = os.path.split(FILE)[0]
        changeWindowTitle(fileName, filePath)
        TXT.insert('START', 'File '+fileName+' saved.\n')
    except:
        TXT.insert('START', 'Can not save file '+os.path.split(ret)[1]+'.\n')

#==============================================================================
# Quick save file
#==============================================================================
def quickSaveFile(event=None):
    global FILE
    try: import tkinter.filedialog as tkFileDialog
    except: import tkFileDialog
    if FILE == '':
        ret = tkFileDialog.asksaveasfilename(filetypes=fileTypes)
        if ret == '' or ret is None or ret == (): # user cancel
            return
        FILE = fixFileString2__(ret)
    try:
        C.convertPyTree2File(t, FILE)
        TXT.insert('START', 'File '+os.path.split(FILE)[1]+' saved.\n')
    except:
        TXT.insert('START', 'Can not save file '+os.path.split(FILE)[1]+'.\n')

#==============================================================================
# Quick reload file
#==============================================================================
def quickReloadFile(event=None):
    global t; global Nb; global Nz; global __FIELD__
    try:
        t = C.convertFile2PyTree(FILE, density=1.)
        t = upgradeTree(t)
        (Nb, Nz) = CPlot.updateCPlotNumbering(t); TKTREE.updateApp()
        if 'tkContainers' in TKMODULES: TKMODULES['tkContainers'].updateApp()
        if TKPLOTXY is not None: TKPLOTXY.updateApp()
        Panels.updateRenderPanel()
        fileName = os.path.split(FILE)[1]
        filePath = os.path.split(FILE)[0]
        TXT.insert('START', 'File %s reloaded.\n'%fileName)
        changeWindowTitle(fileName, filePath)
        errors = Internal.checkPyTree(t, level=5)
        if errors == []: display(t)
        else:
            t = Internal.correctPyTree(t, level=5)
            display(t); Panels.displayErrors(errors, header='Checking pyTree')
    except:
        TXT.insert('START', 'Can not reload file '+os.path.split(FILE)[1]+'.\n')

#==============================================================================
# Quick reload skeleton from file
#==============================================================================
def quickReloadSkeleton(event=None):
    global t; global HANDLE
    try:
        import Converter.Filter as Filter
        h = Filter.Handle(FILE)
        t = h.loadSkeleton(maxDepth=-1)
        HANDLE = h
        TKTREE.updateApp()
        fileName = os.path.split(FILE)[1]
        filePath = os.path.split(FILE)[0]
        TXT.insert('START', 'File %s reloaded.\n'%fileName)
        changeWindowTitle(fileName, filePath)
    except:
        TXT.insert('START', 'Can not reload file '+os.path.split(FILE)[1]+'.\n')

#==============================================================================
# Load skeleton from file
#==============================================================================
def loadFileSkeleton(event=None):
    global FILE; global t
    try: import tkinter.filedialog as tkFileDialog
    except: import tkFileDialog
    files = tkFileDialog.askopenfilenames(
        filetypes=fileTypes, initialfile=FILE, multiple=1)
    if files == '' or files is None or files == (): # user cancel
        return
    files = fixFileString__(files, FILE)
    try:
        FILE = files[0]
        t = []
        for file in files:
            import Converter.Filter as Filter
            h = Filter.Handle(file)
            t2 = h.loadSkeleton(maxDepth=-1)
            if t == []: t = t2
            else: t = C.mergeTrees(t, t2)
        HANDLE = Filter.Handle(FILE)
        TKTREE.updateApp()
        fileName = os.path.split(FILE)[1]
        filePath = os.path.split(FILE)[0]
        TXT.insert('START', 'File %s read.\n'%fileName)
        changeWindowTitle(fileName, filePath)
    except:
        TXT.insert('START', 'Cannot read file '+os.path.split(files[0])[1]+'.\n')
        TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# Save selected zones to a file. Open a dialog.
#==============================================================================
def saveSelZones2File():
    if t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        TXT.insert('START', 'Selection is empty.\n')
        TXT.insert('START', 'Error: ', 'Error'); return
    try: import tkinter.filedialog as tkFileDialog
    except: import tkFileDialog
    ret = tkFileDialog.asksaveasfilename(filetypes=fileTypes)
    if ret == '' or ret is None or ret == (): # user cancel
        return

    base1 = []; base2 = []; base3 = []
    for nz in nzs:
        nob = Nb[nz]+1
        noz = Nz[nz]
        z = t[2][nob][2][noz]
        dim = Internal.getZoneDim(z)
        if dim[4] <= 1: base1.append(z)
        elif dim[4] == 2: base2.append(z)
        else: base3.append(z)
    t2 = C.newPyTree()
    if len(base1) > 0:
        t2 = C.addBase2PyTree(t2, 'Base1', cellDim=1)
        b = Internal.getNodeFromName1(t2, 'Base1')
        b[2] += base1
    if len(base2) > 0:
        t2 = C.addBase2PyTree(t2, 'Base2', cellDim=2)
        b = Internal.getNodeFromName1(t2, 'Base2')
        b[2] += base2
    if len(base3) > 0:
        t2 = C.addBase2PyTree(t2, 'Base3', cellDim=3)
        b = Internal.getNodeFromName1(t2, 'Base3')
        b[2] += base3
    try:
        C.convertPyTree2File(t2, ret)
        fileName = os.path.split(ret)[1]
        TXT.insert('START', 'Selected zones saved to '+fileName+'.\n')
    except:
        TXT.insert('START', 'Can not save file '+os.path.split(ret)[1]+'.\n')

#==============================================================================
# Save selected node to a file. Open a dialog.
#==============================================================================
def saveNode2File():
    if t == []: return

    # Get selected node from tkTree
    if TKTREE is None: return

    node = TKTREE.getCurrentSelectedNode()
    if node is None:
        TXT.insert('START', 'No selected node.\n')
        return

    try: import tkinter.filedialog as tkFileDialog
    except: import tkFileDialog
    ret = tkFileDialog.asksaveasfilename(filetypes=fileTypes)
    if ret == '' or ret is None or ret == (): # user cancel
        return

    try:
        tp = C.newPyTree([])
        tp[2].append(node)
        C.convertPyTree2File(tp, ret)
        fileName = os.path.split(ret)[1]
        TXT.insert('START', 'Selected node saved to '+fileName+'.\n')
    except:
        TXT.insert('START', 'Can not save file '+os.path.split(ret)[1]+'.\n')

#==============================================================================
# Load file that stores *Cassiopee* preferences
# load $HOME/.cassiopee or $HOME/.cassiopee/config
# Remplit le dictionnaire PREFS
#==============================================================================
def loadPrefFile():
    global PREFS
    homePath = os.path.expanduser('~')
    if homePath is None: homePath = os.getenv('HOME')
    if homePath is None: homePath = os.getenv('USERPROFILE')
    if homePath is None: homePath = ''
    kdir = homePath+'/.cassiopee'
    exist = os.path.exists(kdir)
    if not exist: savePrefFile(); return []
    isdir = os.path.isdir(kdir)
    if not isdir: file = open(homePath+'/.cassiopee', 'r')
    else:
        try: file = open(homePath+'/.cassiopee/config', 'r')
        except: savePrefFile(); return []
    plist = file.read()
    file.close()
    plist = plist.split('\n')
    for i in plist:
        if i != '':
            a = i.split(':',1)
            if len(a) == 2:
                prefName = a[0]
                val = a[1]
                PREFS[prefName] = val
    return

#==============================================================================
# Essai d'importer ttk
#==============================================================================
def importTtk():
    try: import tkinter.ttk as ttk
    except ImportError:
        try: import ttk
        except ImportError: ttk = None
    return ttk

#==============================================================================
# Set preferences
# IN: PREFS dict
# L'etat du logiciel depend de :
# L'etat de CPlot (modifiable par des CPlot.get et CPlot.set)
# Les variables globales de CPlot.Tk
#==============================================================================
def setPrefs():
    global __UNDO__, __ONEOVERN__, FONTTYPE, FONTSIZE
    global GENERALFONT, FRAMEFONT, LABELFONT, MENUFONT, BUTTONFONT
    global TEXTFONT, MSGDFONT, FIXEDFONT, FRAMESTYLE
    for i in PREFS:
        val = PREFS[i]
        if i == 'tkViewMode':
            if val == 'Mesh': CPlot.setState(mode=0)
            elif val == 'Solid': CPlot.setState(mode=1)
            elif val == 'Render': CPlot.setState(mode=2)
            elif val == 'Scalar': CPlot.setState(mode=3)
            elif val == 'Vector': CPlot.setState(mode=4)
        elif i == 'tkViewMeshStyle':
            if val == 'Red wires+solid': style = 0
            elif val == 'Multicolor wireframes': style = 1
            elif val == 'Multicolor wires+solid': style = 2
            elif val == 'Black wires+solid': style = 3
            elif val == 'Black wires+solid2': style = 4
            else: style = 0
            CPlot.setState(meshStyle=style)
        elif i == 'tkViewSolidStyle':
            if val == 'Monocolor/1-side': style = 0
            elif val == 'Multicolor/2-sides': style = 1
            elif val == 'White/2-sides': style = 3
            elif val == 'Multicolor/outlined': style = 4
            else: style = 0
            CPlot.setState(solidStyle=style)
        elif i == 'tkEffectsAngle':
            val = float(val)
            CPlot.setState(viewAngle=val)
        elif i == 'tkEffectsShadow':
            if val == '1': CPlot.setState(shadow=1)
            else: CPlot.setState(shadow=0)
        elif i == 'tkEffectsDOF':
            if val == '1': CPlot.setState(dof=1)
            else: CPlot.setState(dof=0)
        elif i == 'tkViewNiso':
            val = int(val)
            CPlot.setState(niso=val)
        elif i == 'tkViewColormap':
            if val == 'Blue2Red': style = 0
            elif val == 'Green2Red': style = 2
            elif val == 'BiColorRGB': style = 4
            elif val == 'BiColorHSV': style = 6
            elif val == 'Diverging': style = 8
            elif val == 'TriColorRGB': style = 10
            elif val == 'TriColorHSV': style = 12
            else: style = 0
            if 'tkViewIsoLight' in PREFS:
                if PREFS['tkViewIsoLight'] == 'IsoLight on': style += 1
            else: style += 1
            CPlot.setState(colormap=style)
        elif i == 'tkViewLegend':
            val = int(val)
            CPlot.setState(displayIsoLegend=val)
        elif i == 'tkViewDim':
            if val == '2D': CPlot.setDim(2)
            elif val == '3D': CPlot.setDim(3)
        elif i == 'tkViewEdgeA':
            if val == '1': CPlot.setState(edgifyActivatedZones=1)
            else: CPlot.setState(edgifyActivatedZones=0)
        elif i == 'tkViewEdgeD':
            if val == '1': CPlot.setState(edgifyDeactivatedZones=1)
            else: CPlot.setState(edgifyDeactivatedZones=0)
        elif i == 'tkViewScalarStyle':
            if val == 'Bands': style = 0
            elif val == 'Bands+mesh': style = 1
            CPlot.setState(scalarStyle=style)
        elif i == 'selectionStyle':
            if val == 'Blue': style = 0
            elif val == 'Alpha': style = 1
            CPlot.setState(selectionStyle=style)
        elif i == 'simplifyOnDrag':
            if val == 'None': style = 0
            elif val == 'BBox': style = 1
            CPlot.setState(simplifyOnDrag=style)

        elif i == 'undo':
            val = int(val)
            if val == 0: __UNDO__ = False
            else: __UNDO__ = True
        elif i == 'displayInfo':
            val = int(val)
            CPlot.setState(displayInfo=val)
            CPlot.setState(displayBB=val)
        elif i == 'bgColor':
            val = int(val)
            CPlot.setState(bgColor=val)
        elif i == 'backgroundFile':
            CPlot.setState(backgroundFile=val)
        elif i == 'envmap':
            CPlot.setState(envmap=val)
        elif i == 'fontType':
            FONTTYPE = val
            (GENERALFONT, FRAMEFONT, LABELFONT, MENUFONT, BUTTONFONT,
             TEXTFONT, MSGDFONT, FIXEDFONT) = defineFonts(FONTTYPE, FONTSIZE)
        elif i == 'fontSize':
            val = int(val)
            FONTSIZE = val
            (GENERALFONT, FRAMEFONT, LABELFONT, MENUFONT, BUTTONFONT,
             TEXTFONT, MSGDFONT, FIXEDFONT) = defineFonts(FONTTYPE, FONTSIZE)
        elif i == 'frameStyle':
            if val == 'GROOVE': FRAMESTYLE = TK.GROOVE
            elif val == 'RIDGE': FRAMESTYLE = TK.RIDGE
            elif val == 'SUNKEN': FRAMESTYLE = TK.SUNKEN
            elif val == 'RAISED': FRAMESTYLE = TK.RAISED
        elif i == 'viewAngle':
            val = float(val)
            CPlot.setState(viewAngle=val)
        elif i == 'GridCoordinatesContainer':
            Internal.__GridCoordinates__ = val
        elif i == 'FlowSolutionNodesContainer':
            Internal.__FlowSolutionNodes__ = val
        elif i == 'FlowSolutionCentersContainer':
            Internal.__FlowSolutionCenters__ = val
        elif i == 'tkPerfoPoints':
            if val == 'All points': __ONEOVERN__ = 1
            elif val == 'One over 2': __ONEOVERN__ = 2
            elif val == 'One over 3': __ONEOVERN__ = 3
            elif val == 'One over 4': __ONEOVERN__ = 4
            elif val == 'Exclusive': __ONEOVERN__ = 0

#==============================================================================
# Save pref file that stores personal modules paths and prefs
# Save $HOME/.cassiopee or $HOME/.cassiopee/config
# IN: pathList: la liste des chemins des modules additionnels
#==============================================================================
def savePrefFile():
    homePath = os.path.expanduser('~')
    if homePath is None: homePath = os.getenv('HOME')
    if homePath is None: homePath = os.getenv('USERPROFILE')
    if homePath is None: homePath = ''
    kdir = homePath+'/.cassiopee'
    exist = os.path.exists(kdir)
    if not exist:
        os.makedirs(kdir)
        file = open(homePath+'/.cassiopee/config', 'w')
    else:
        isdir = os.path.isdir(kdir) # ok
        if not isdir:
            file = open(homePath+'/.cassiopee', 'w')
            #os.rename(kdir, homePath+'/.cassiopee_save')
            #os.makedirs(kdir)
            #os.rename(homePath+'/.cassiopee_save', homePath+'/.cassiopee/config')
        else:
            file = open(homePath+'/.cassiopee/config', 'w')
    for i in PREFS: file.write(i+':'+PREFS[i]+'\n')
    file.close()

#==============================================================================
# Undo (1 niveau)
#==============================================================================
def undo(event=None):
    if not __UNDO__: return
    global t, Nb, Nz
    if t == []: return
    if __MAINTREE__ <= 0: CPlot.display(t); return
    if tp == []: return
    t = tp
    (Nb, Nz) = CPlot.updateCPlotNumbering(t); TKTREE.updateApp()
    TXT.insert('START', 'Undo performed.\n')
    display(t)

#==============================================================================
# Sauvegarde l'arbre courant dans l'arbre precedent
# Doit etre appele avant de modifier t
#==============================================================================
def saveTree():
    if not __UNDO__: return
    global dt, tp
    dt = [] # pour economiser de la memoire
    if t == []: tp = C.newPyTree(['Base'])
    else: tp = Internal.copyRef(t)

#==============================================================================
def Quit():
    #WIDGETS['masterWin'].destroy()
    #WIDGETS['masterWin'].quit()
    del iconics.PHOTO
    os._exit(0)

#==============================================================================
def setCPlotMode0(): # mode mesh
    CPlot.setMode(0)

#==============================================================================
def setCPlotMode1(): # mode solid
    CPlot.setMode(1)

#==============================================================================
def setCPlotMode2(): # mode render
    CPlot.setMode(2)
#==============================================================================
def setCPlotMode3(): # mode scalar
    CPlot.changeVariable()

#==============================================================================
def changeCPlotStyle():
    CPlot.changeStyle()

#==============================================================================
def setLocCenters():
    global __LOCATION__
    __LOCATION__ = 'centers'
    menu = WIDGETS['cplotMenu']
    menu.entryconfig(0, label='Display Nodes')
    menu.entryconfig(1, label='Display Centers*')
    TXT.insert('START', 'Centers displayed.\n')
    display(t, location=__LOCATION__)

#==============================================================================
def setLocNodes():
    global __LOCATION__
    __LOCATION__ = 'nodes'
    menu = WIDGETS['cplotMenu']
    menu.entryconfig(0, label='Display Nodes*')
    menu.entryconfig(1, label='Display Centers')
    TXT.insert('START', 'Nodes displayed.\n')
    display(t, location=__LOCATION__)

#==============================================================================
def lookFor():
    if t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []: CPlot.fitView(); return
    CPlot.lookFor()
    showSelectionInTkTree()

#==============================================================================
def changeCPlotBlanking():
    CPlot.changeBlanking()
    b = CPlot.getState("blanking")
    if b == 0: TXT.insert('START', 'Blanking turned off.\n')
    else: TXT.insert('START', 'Blanking turned on.\n')

#==============================================================================
def cplotExport():
    global EXPORTFILE
    try: import tkinter.filedialog as tkFileDialog
    except: import tkFileDialog
    ret = tkFileDialog.asksaveasfilename(
        title='Export as...',
        initialfile=EXPORTFILE,
        filetypes=[('Portable Network Graphics', '*.png'),
                   ('Portable pixmap', '*.ppm'),
                   ('Bitmap Postscript', '*.ps'),
                   ('Mpeg movie', '*.mpeg')])
    if ret == '' or ret is None or ret == (): # user cancel
        return
    try: exportResolution = PREFS['exportResolution']
    except: exportResolution = 'None'
    try:
        fileString = fixFileString2__(ret)
        CPlot.setState(exportResolution=exportResolution)
        CPlot.setState(export=fileString)
        myFile = os.path.split(fileString)[1]
        ext = os.path.splitext(myFile)[1]
        if ext == '.mpeg' or ext == '.MPEG':
            CPlot.setState(continuousExport=1)
        EXPORTFILE = fileString
        TXT.insert('START', 'File '+os.path.split(EXPORTFILE)[1]+
                   ' exported.\n')
    except:
        TXT.insert('START', 'Can not export to file '+
                   os.path.split(fileString)[1]+'.\n')
    return

#==============================================================================
def finalizeExport():
    CPlot.finalizeExport()

#==============================================================================
def rmBlock():
    global t, Nb, Nz
    if t == []: return
    if __MAINTREE__ <= 0:
        TXT.insert('START', 'Fail on a temporary tree.\n')
        TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []: return

    deletedZoneNames = []
    for nz in nzs:
        nob = Nb[nz]+1
        noz = Nz[nz]
        deletedZoneNames.append(t[2][nob][0]+Internal.SEP1+t[2][nob][2][noz][0])

    saveTree()
    t = CPlot.deleteSelection(t, Nb, Nz, nzs)

    TXT.insert('START', 'Selected zones deleted.\n')
    (Nb, Nz) = CPlot.updateCPlotNumbering(t); TKTREE.updateApp()
    CPlot.delete(deletedZoneNames)
    CPlot.render()

#==============================================================================
def copyBlock():
    global t, Nb, Nz
    if t == []: return
    if __MAINTREE__ <= 0:
        TXT.insert('START', 'Fail on a temporary tree.\n')
        TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []: return

    saveTree()
    t = C.addBase2PyTree(t, 'COPY', 3)
    base = Internal.getNodeFromName1(t, 'COPY')
    for nz in nzs:
        nob = Nb[nz]+1
        noz = Nz[nz]
        z = t[2][nob][2][noz]
        zp = Internal.copyRef(z)
        zp[0] = C.getZoneName(z[0]+'.dup')
        base[2].append(zp)

    TXT.insert('START', 'Selected zones duplicated.\n')
    (Nb, Nz) = CPlot.updateCPlotNumbering(t); TKTREE.updateApp()
    display(t)

#==============================================================================
def unselectAll():
    if t == []: return
    if __MAINTREE__ <= 0: CPlot.display(t)
    #bases = Internal.getBases(t)
    selected = []
    s = 0

    nodes = Internal.getZones(t)
    for no in range(len(nodes)): selected.append( (no, s) )

    TXT.insert('START', 'Tree unselected.\n')
    CPlot.setSelectedZones(selected)

#==============================================================================
def toggleSelectAll():
    if t == []: return
    if __MAINTREE__ <= 0: CPlot.display(t)
    bases = Internal.getBases(t)
    selected = []
    s = -1

    for b in bases:
        baseName = b[0]
        nodes = Internal.getNodesFromType1(b, 'Zone_t')
        if nodes != []:
            noz = CPlot.getCPlotNumber(t, baseName, nodes[0][0])
            if s == -1:
                sp = CPlot.getSelectedStatus(noz)
                if sp == 0: s = 1
                else: s = 0
                break

    nodes = Internal.getZones(t)
    for no in range(len(nodes)): selected.append( (no, s) )

    if s == 0: TXT.insert('START', 'Tree unselected.\n')
    elif s == 1: TXT.insert('START', 'Tree selected.\n')
    CPlot.setSelectedZones(selected)

#==============================================================================
# Inverse les zones activee et desactivees
def revertActivated():
    if t == []: return
    nz = len(Internal.getZones(t))
    nzs = CPlot.getActiveZones()
    active = [(i,1) for i in range(nz)]
    for n in nzs: active[n] = (n,0)
    CPlot.setActiveZones(active)

#==============================================================================
class infoBulle(TK.Toplevel):
    # IN: text: text to be displayed
    # IN: temps: delai for display
    # IN: btype=0: standard, btype=1: applet menu info
    def __init__(self, parent=None, text='', temps=1000, btype=0,
                 textVariable=None):
        TK.Toplevel.__init__(self, parent, bd=1, bg='black')
        self.tps = temps
        self.parent = parent
        self.withdraw()
        self.overrideredirect(1)
        self.transient()
        self.btype = btype
        if btype == 1: # menu
            if textVariable is not None:
                l = TK.Label(self, textvariable=textVariable, bg="white",
                             justify='right', takefocus=0)
            else:
                l = TK.Label(self, text=text, bg="white", justify='right',
                             takefocus=0)
        else: # std label
            if textVariable is not None:
                l = TK.Label(self, textvariable=textVariable, bg="yellow",
                             justify='left', takefocus=0)
            else:
                l = TK.Label(self, text=text, bg="yellow", justify='left',
                             takefocus=0)
        l.update_idletasks()
        l.pack()
        l.update_idletasks()
        self.label = l
        self.tipwidth = l.winfo_width()
        self.tipheight = l.winfo_height()
        if btype == 1: # menu
            self.parent.bind('<Enter>', self.delai2)
            #self.parent.bind('<Button-1>', self.efface)
            self.parent.bind('<Leave>', self.efface)
        else:  # std
            self.parent.bind('<Enter>', self.delai)
            #self.parent.bind('<Button-1>', self.efface)
            self.parent.bind('<Leave>', self.efface)
    def delai(self, event): # pour le cas std
        self.action = self.parent.after(self.tps, self.affiche)
    def delai2(self, event): # pour les menus
        self.parent.focus_set() # frame est le parent
        self.action = self.parent.after(self.tps, self.affiche)
    def affiche(self):
        self.update_idletasks()
        if self.btype == 1: # menu
            posX = self.parent.winfo_rootx()+self.parent.winfo_width()
            posY = self.parent.winfo_rooty()
        else:
            posX = self.parent.winfo_rootx()+self.parent.winfo_width()
            posY = self.parent.winfo_rooty()+self.parent.winfo_height()

        if posX > self.winfo_screenwidth():
            posX = self.winfo_screenwidth() - self.tipwidth
        elif posX + self.tipwidth > self.winfo_screenwidth():
            posX -= self.winfo_width()
        if posY + self.tipheight > self.winfo_screenheight():
            posY -= self.parent.winfo_height() + self.winfo_height()

        self.geometry('+%d+%d'%(posX,posY))
        self.deiconify()
    def efface(self, event):
        self.withdraw()
        self.parent.after_cancel(self.action)

def stt__(s, mode=0):
    s = s.replace(' ', '')
    if mode == 0:
        s = s.replace('(', '')
        s = s.replace(')', '')
        s = s.replace(',', ';')
    s = s.replace('|', ';')
    c = 0; out = ''; o = 0
    for i in s:
        if i == ';':
            o += 1
            if o == 1: out += i
        else:
            o = 0
            out += i
        c += 1
    s = out.split(';')
    #s = [i for i in s if i]
    return s

#==============================================================================
# Extrait les variables de chaines recuperes des widgets
# ces chaines sont de la forme :
# type0,1,2: 'var1; var2; var3' ou 'var1'
# type3: 'ind;(ind);(i,j,k)'
# On retourne une liste de variables.
# type0,1,2: [var1,var2,var3]
# type3: [ind, ind, (i,j,k)]
# Si la liste est vide, c'est que la conversion a echoue.
# IN: varString: chaine de type 'var1 ; vars2'
# IN: type: type de conversion tentee
# OUT: une liste contenant les valeurs correspondants [v1, v2, ...] ou une
# liste d'indices
# type=0 -> string; type=1 -> float; type=2 -> int; type=3 -> indices
#==============================================================================
def varsFromWidget(varString, type=0):
    if type == 0: # string
        return stt__(varString)
    elif type == 1: # float
        s = stt__(varString)
        ret = []
        for i in s:
            try: val = float(i)
            except: val = -1.
            ret.append(val)
        return ret
    elif type == 2: # int
        s = stt__(varString)
        ret = []
        for i in s:
            try: val = int(i)
            except:
                try: val = int(float(i))
                except: val = -1
            ret.append(val)
        return ret
    elif type == 3: # indices
        s = stt__(varString, mode=1)
        ret = []

        #from ast import literal_eval
        #for i in s:
        #    try: val = literal_eval(i); ret.append(val)
        #    except: ret.append(0)

        for i in s:
            try:
                val = eval(i)
                if isinstance(val, list): ret += val
                else: ret.append(val)
            except: ret.append(0)

        return ret
    return []

#==============================================================================
# Retourne un liste de zones
# En mode MAINTREE:
# Si une ou plusieurs zones sont selectionnees: retourne ces zones
# Sinon: retourne les zones actives
# En mode tree temporaire: retourne []
#==============================================================================
def getValidZones():
    if __MAINTREE__ == 1:
        nzs = CPlot.getSelectedZones()
        if nzs != []:
            out = []
            zones = Internal.getZones(t)
            for nz in nzs:
                nob = Nb[nz]+1
                noz = Nz[nz]
                z = t[2][nob][2][noz]
                out.append(z)
            return out
        nzs = CPlot.getActiveZones()
        if nzs != []:
            out = []
            zones = Internal.getZones(t)
            for nz in nzs:
                nob = Nb[nz]+1
                noz = Nz[nz]
                z = t[2][nob][2][noz]
                out.append(z)
            return out
        zones = Internal.getZones(t)
        return zones
    else: return []

#==============================================================================
# Tool bar
# Create a tool bar on top of win
# mode = 0: for Cassiopee
# mode = 1: for kcgnsview
#==============================================================================
def toolBar(win, mode=0):
    # Change l'icone de la fenetre
    win.tk.call('wm', 'iconphoto', win._w, iconics.PHOTO[13])

    frame = TTK.Frame(win)
    frame.rowconfigure(0, weight=0)

    width, height = 30,30

    if mode == 0: # Cassiopee
        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[0], borderwidth=0,
                        command=quickSaveFile)
        BB = infoBulle(parent=B, text='Save.')
        B.grid(row=0, column=0)

        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[1], borderwidth=0,
                        command=undo)
        BB = infoBulle(parent=B, text='Undo.')
        B.grid(row=0, column=1)

        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[11], borderwidth=0,
                        command=quickReloadFile)
        BB = infoBulle(parent=B, text='Reload current file.')
        B.grid(row=0, column=2)

        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[3], borderwidth=0,
                        command=cplotExport)
        BB = infoBulle(parent=B, text='Export image.')
        B.grid(row=0, column=3)

        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[4], borderwidth=0,
                        command=lookFor)
        BB = infoBulle(parent=B, text='Fit view to selection\nor fit to full size.')
        B.grid(row=0, column=4)

        # B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
        #             image=iconics.PHOTO[5], borderwidth=0,
        #             command=toggleSelectAll)
        # BB = infoBulle(parent=B, text='Select/unselect all zones.')
        # B.grid(row=0, column=5)
        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[14], borderwidth=0,
                        command=copyBlock)
        BB = infoBulle(parent=B, text='Duplicate selected zones.')
        B.grid(row=0, column=5)

        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[2], borderwidth=0,
                        command=rmBlock)
        BB = infoBulle(parent=B, text='Remove selected zones.')
        B.grid(row=0, column=6)

        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[6], borderwidth=0,
                        command=viewDeactivatedZones)
        BB = infoBulle(parent=B, text='View deactivated zones.')
        B.grid(row=0, column=7)

        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[12], borderwidth=0,
                        command=revertActivated)
        BB = infoBulle(parent=B, text='Switch active <-> inactive zones.')
        B.grid(row=0, column=8)

        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[7], borderwidth=0,
                        command=displayMainTree)
        BB = infoBulle(parent=B, text='Display main tree.')
        B.grid(row=0, column=9)
    else: # kcgnsview
        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[11], borderwidth=0,
                        command=quickReloadSkeleton)
        BB = infoBulle(parent=B, text='Reload current file.')
        B.grid(row=0, column=0)
        B = TTK.Button2(frame, compound=TK.TOP, width=width, height=height,
                        image=iconics.PHOTO[1], borderwidth=0,
                        command=undo)
        BB = infoBulle(parent=B, text='Undo.')
        B.grid(row=0, column=1)

    frame.grid(sticky=TK.NSEW, columnspan=2)

#==============================================================================
# Minimum application: win, menu, txt
# IN: title: application title
# IN: show: create and show if true
# IN: mode: 0 (Cassiopee), 1 (kcgnsview)
#==============================================================================
def minimal(title, show=True, mode=0):

    global TXT, TKTREE
    win = TK.Tk()
    WIDGETS['masterWin'] = win
    if not show: win.withdraw()
    win.title(title)

    # ttk style
    TTK.installLocalThemes(win)
    TTK.setTheme(PREFS.get('guitheme', 'None'))

    win.protocol("WM_DELETE_WINDOW", Quit)
    win.bind('<Control-Key-s>', quickSaveFile)
    win.bind('<Control-Key-o>', loadFile)
    win.bind('<Control-Key-p>', Panels.openLoadPanel)
    win.bind('<Control-Key-z>', undo)
    win.option_add('*Font', GENERALFONT)
    win.option_add('*Label.font', LABELFONT)
    win.option_add('*Menu.font', MENUFONT)
    win.option_add('*Button.font', BUTTONFONT)
    win.option_add('*Dialog.msg.font', MSGDFONT)
    win.columnconfigure(0, weight=1)
    win.columnconfigure(1, weight=1)
    win.rowconfigure(0, weight=0)
    win.rowconfigure(1, weight=0)
    win.rowconfigure(2, weight=0)
    win.rowconfigure(3, weight=1)
    win.resizable(0,0) # fixed size window
    #win.minsize(325, 325)

    menu = TK.Menu(win, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
    win.config(menu=menu)

    # menu file
    if mode == 0: # Cassiopee
        file = TK.Menu(menu, tearoff=0, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
        menu.add_cascade(label='File', menu=file)
        file.add_command(label='Open', accelerator='Ctrl+o', command=loadFile)
        file.add_command(label='Add', command=addFile)
        file.add_command(label='Open load panel', accelerator='Ctrl+p', command=Panels.openLoadPanel)
        file.add_separator()
        file.add_command(label='Save', accelerator='Ctrl+s', command=quickSaveFile)
        file.add_command(label='Save as...', command=saveFile)
        file.add_command(label='Save selected zones', command=saveSelZones2File)
        file.add_command(label='Save selected node', command=saveNode2File)
        file.add_separator()
        file.add_command(label='Quit', command=Quit)
    else: # kcgnsview
        file = TK.Menu(menu, tearoff=0, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
        menu.add_cascade(label='File', menu=file)
        file.add_command(label='Open', accelerator='Ctrl+o', command=loadFileSkeleton)
        file.add_separator()
        file.add_command(label='Quit', command=Quit)

    # menu CPlot
    if mode == 0:
        cplot = TK.Menu(menu, tearoff=0, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
        WIDGETS['cplotMenu'] = cplot
        menu.add_cascade(label='CPlot', menu=cplot)
        cplot.add_command(label='Display Nodes*', command=setLocNodes)
        cplot.add_command(label='Display Centers', command=setLocCenters)
        cplot.add_command(label='Toggle blanking (on/off)', command=changeCPlotBlanking)
        cplot.add_separator()
        cplot.add_command(label='Export image/movie', command=cplotExport)
        cplot.add_command(label='Finalize movie', command=finalizeExport)
        cplot.add_command(label='E-mail image/report bug', command=mail2Friends)
        cplot.add_command(label='Save image in doc/blog', command=save2Doc)
        cplot.add_separator()
        cplot.add_command(label='Undo change', accelerator='Ctrl+z', command=undo)
    else:
        cplot = TK.Menu(menu, tearoff=0, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
        WIDGETS['cplotMenu'] = cplot
        menu.add_cascade(label='CPlot', menu=cplot)

    # Menu specific tools
    if mode == 0:
        tools = TK.Menu(menu, tearoff=0, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
        menu.add_cascade(label='Tools', menu=tools)
        tools.add_command(label='Save prefs', command=savePrefFile)
        tools.add_separator()
        #tools.add_command(label='Activate key', command=Panels.activation)
    else:
        tools = TK.Menu(menu, tearoff=0, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
        menu.add_cascade(label='Tools', menu=tools)

    # Menu Help
    help = TK.Menu(menu, tearoff=0, bg=TTK.BACKGROUNDCOLOR, fg=TTK.FOREGROUNDCOLOR)
    menu.add_cascade(label='Help', menu=help)
    help.add_command(label='Online User doc', command=getOnlineDoc)
    help.add_command(label='Online Help forum', command=getOnlineForum)
    help.add_command(label='Online Tutorials', command=getOnlineTutorials)
    help.add_separator()
    help.add_command(label='About', command=Panels.about)

    win.config(menu=menu)

    toolBar(win, mode=mode)

    # Text only
    #TXT = TK.Text(win, width=30, height=1, background='White', font=TEXTFONT)
    #TXT.tag_config("Error", foreground="red")
    #TXT.tag_config("Warning", foreground="green")
    #TXT.mark_set('START', TK.INSERT)
    #TXT.mark_gravity('START', TK.LEFT)
    #TXT.grid(sticky=TK.EW, columnspan=2)

    # Text + search entry pour Cortano
    F = TTK.Frame(win, width=30, height=1, takefocus=1)
    F.columnconfigure(0, weight=1)
    TXT = TTK.Text(F, width=30, height=1, background='White', font=TEXTFONT)
    TXT.tag_config("Error", foreground="red")
    TXT.tag_config("Warning", foreground="green")
    TXT.mark_set('START', TK.INSERT)
    TXT.mark_gravity('START', TK.LEFT)
    TXT.grid(sticky=TK.EW)
    from . import tkSearchBar
    E = tkSearchBar.createSearchBar2(F)
    E.grid(row=1, sticky=TK.NSEW)
    F.grid(sticky=TK.NSEW, columnspan=2)

    try: TKTREE = __import__('tkTree'); TKTREE.createApp(win)
    except: TKTREE = None

    return (win, menu, file, tools)

#==============================================================================
# Minimum application: tktree, notebook, frames, menu, txt
# IN: title: application title
# IN: show: create and show if true
# IN: mode: 0 (Cassiopee), 1 (kcgnsview)
#==============================================================================
def minimal2(title, show=True, mode=0):
    global TKTREE
    (win, menu, file, tools) = minimal(title, show, mode)

    # Frame container
    F = TTK.Frame(win)
    F.grid(columnspan=2, sticky=TK.NSEW)
    F.columnconfigure(0, weight=0)
    F.columnconfigure(1, weight=1)
    F.rowconfigure(0, weight=0)
    F.rowconfigure(1, weight=0)
    F.rowconfigure(2, weight=1)
    F.rowconfigure(3, weight=0)

    # Cree le TkTree (colonne1)
    TKTREE = __import__('tkTree')
    TKTREE.createApp(F); TKTREE.showApp()

    # Cree le notebook avec boutons
    nb = noteBook(F, TK.LEFT, menu)
    frames = []; menus = []; buttons = []

    # Frame pour tree
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    tree = TK.Menu(menu, tearoff=0); menus.append(tree)
    bt = nb.add_screen(frame, 'Tree', tree)
    frames.append(frame); buttons.append(bt)
    WIDGETS['TreeNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['TreeNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour state
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    state = TK.Menu(menu, tearoff=0); menus.append(state)
    bt = nb.add_screen(frame, 'State', state)
    frames.append(frame); buttons.append(bt)
    WIDGETS['StateNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['StateNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour edge
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    edge = TK.Menu(menu, tearoff=0); menus.append(edge)
    bt = nb.add_screen(frame, 'Edge', edge)
    frames.append(frame); buttons.append(bt)
    WIDGETS['EdgeNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['EdgeNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour surf
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    surf = TK.Menu(menu, tearoff=0); menus.append(surf)
    bt = nb.add_screen(frame, 'Surf', surf)
    frames.append(frame); buttons.append(bt)
    WIDGETS['SurfNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['SurfNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour mesh
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    mesh = TK.Menu(menu, tearoff=0); menus.append(mesh)
    bt = nb.add_screen(frame, 'Mesh', mesh)
    frames.append(frame); buttons.append(bt)
    WIDGETS['MeshNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['MeshNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour block
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    block = TK.Menu(menu, tearoff=0); menus.append(block)
    bt = nb.add_screen(frame, 'Block', block)
    frames.append(frame); buttons.append(bt)
    WIDGETS['BlockNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['BlockNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour BC
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    bc = TK.Menu(menu, tearoff=0); menus.append(bc)
    bt = nb.add_screen(frame, 'BC', bc)
    frames.append(frame); buttons.append(bt)
    WIDGETS['BCNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['BCNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour motion
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    bc = TK.Menu(menu, tearoff=0); menus.append(bc)
    bt = nb.add_screen(frame, 'Motion', bc)
    frames.append(frame); buttons.append(bt)
    WIDGETS['MotionNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['MotionNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour Solver
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    solver = TK.Menu(menu, tearoff=0)
    bt = nb.add_screen(frame, 'Solver', solver); menus.append(solver)
    frames.append(frame); buttons.append(bt)
    WIDGETS['SolverNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['SolverNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour Post
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    post = TK.Menu(menu, tearoff=0); menus.append(post)
    bt = nb.add_screen(frame, 'Post', post)
    frames.append(frame); buttons.append(bt)
    WIDGETS['PostNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['PostNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour Visu
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    visu = TK.Menu(menu, tearoff=0); menus.append(visu)
    bt = nb.add_screen(frame, 'Visu', visu)
    frames.append(frame); buttons.append(bt)
    WIDGETS['VisuNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['VisuNoteBook'].grid(sticky=TK.NSEW)

    # Frame pour Render
    frame = TTK.Frame(nb())
    frame.columnconfigure(0, weight=1)
    render = TK.Menu(menu, tearoff=0); menus.append(render)
    bt = nb.add_screen(frame, 'Render', render)
    frames.append(frame); buttons.append(bt)
    WIDGETS['RenderNoteBook'] = TTK.NoteBook(frame)
    WIDGETS['RenderNoteBook'].grid(sticky=TK.NSEW)

    # Sauvegarde
    WIDGETS['noteBook'] = nb
    WIDGETS['noteBookFrames'] = frames # chaque frame est le conteneur des applets
    WIDGETS['noteBookMenus'] = menus # petits menus
    WIDGETS['noteBookButtons'] = buttons # boutons du noteBook

    return (win, frames, menu, menus, file, tools)

#==============================================================================
# IN: app: nom d'applet (tkChose ou --- ou Sub/tkChose)
# si app=---, un separateur est affiche
# si app=Sub/tkChose, un sous-menu Sub est ajoute
# IN: menu: menu ou on ajoute l'app
# IN: frame: la frame de l'app
# IN: submenus: dict to keep trace of allready created submenus
# OUT: TKMODULES: le dictionnaire des modules importes
#==============================================================================
def addMenuItem(app, menu, frame, submenus, auto):
    global MODULES, MODULEFRAMES
    app = app.split('/')
    if len(app) == 2: submenu = app[0]; app = app[1]
    else: submenu = None; app = app[0]

    if submenu is None:
        if app == '---': menu.add_separator()
        else:
            TKMODULES[app] = None
            TKMODULEFRAMES[app] = frame
            name = app; name = '  '+name
            menu.add_command(label=name, command=lambda x=app:openApp(x))
            if auto[app] == 1: openApp(app)

    else: # submenu
        if submenu in submenus:
            myMenu = submenus[submenu]
        else:
            myMenu = TK.Menu(menu, tearoff=0)
            submenus[submenu] = myMenu
            menu.add_cascade(label=submenu, menu=myMenu)
        if app == '---': myMenu.add_separator()
        else:
            TKMODULES[app] = None
            TKMODULEFRAMES[app] = frame
            name = app; name = '  '+name
            myMenu.add_command(label=name, command=lambda x=app:openApp(x))
            if auto[app] == 1: openApp(app)

#==============================================================================
class noteBook:
    # initialization. receives the master widget
    # reference and the notebook orientation
    def __init__(self, master, side=TK.LEFT, menu=None):
        self.active_fr = None
        self.count = 0
        self.choice = TK.IntVar(value=0)
        self.menu = menu

        # allows the TOP and BOTTOM radiobuttons' positioning.
        if side in (TK.TOP, TK.BOTTOM): self.side = TK.LEFT
        else: self.side = TK.TOP

        # creates notebook's frames structure
        self.rb_fr = TTK.Frame(master, borderwidth=0, relief=TK.RIDGE)
        self.rb_fr.grid(sticky='ewn', column=0, row=0, rowspan=2)
        self.screen_fr = TTK.Frame(master, borderwidth=1, relief=TK.GROOVE)
        self.screen_fr.columnconfigure(0, weight=1)
        self.screen_fr.grid(column=1, row=1, sticky=TK.NSEW)

    # Return a master frame reference for the external frames (screens)
    def __call__(self):
        return self.screen_fr

    # add a new frame (screen) to the (bottom/left of the) notebook
    def add_screen(self, fr, title, menu_fr=None):
        if menu_fr is None:
            b = TTK.Radiobutton(self.rb_fr, text=title,
                                #image=image, compound=TK.TOP, pady=0, border=TK.ROUND,
                                offrelief=TK.GROOVE,
                                indicatoron=False, selectcolor='#ffffff',
                                variable=self.choice, value=self.count)
        else:
            b = TTK.Radiobutton(self.rb_fr, text=title,
                                #image=image, compound=TK.TOP, pady=0, border=TK.ROUND,
                                offrelief=TK.GROOVE,
                                indicatoron=False, selectcolor='#ffffff',
                                variable=self.choice, value=self.count,
                                command=lambda: self.display(fr, menu_fr))
            b.bind('<ButtonRelease-3>', lambda event: self.displayMenu(event, fr, menu_fr, b))
            menu_fr.config(bg=TTK.BACKGROUNDCOLOR)
            menu_fr.config(fg=TTK.FOREGROUNDCOLOR)
            b.grid(sticky=TK.EW)

        # ensures the first frame will be the first selected/enabled
        if not self.active_fr:
            fr.grid(sticky=TK.NSEW)
            self.active_fr = fr
            if self.menu is not None and menu_fr is not None:
                if UPDATEAPPSMENU:
                    self.menu.insert_cascade(index=3, label='Apps', menu=menu_fr)

        self.count += 1

        # returns a reference to the newly created
        # radiobutton (allowing its configuration/destruction)
        return b

    # hides the former active frame and shows
    # another one, keeping its reference
    def display(self, fr, menu_fr):
        self.active_fr.grid_forget()
        fr.grid(sticky=TK.NSEW)
        self.active_fr = fr
        if self.menu is not None and menu_fr is not None:
            if UPDATEAPPSMENU:
                self.menu.delete(3)
                self.menu.insert_cascade(index=3, label='Apps', menu=menu_fr)
                #self.menu.entryconfigure(3, menu=menu_fr)

    # Display the gleize little menu
    def displayMenu(self, event, fr, menu_fr, b):
        menu_fr.tk_popup(event.x_root+50, event.y_root, 0)
        TTK.selectRadioButton(b)
        self.display(fr, menu_fr)

#==============================================================================
def getOnlineDoc():
    try:
        import webbrowser
        TXT.insert('START', 'Opening online documentation.\n')
        webbrowser.open('https://cassiopee.onera.fr/Userguide.html')
    except:
        TXT.insert('START', 'Can not open online documentation.\n')
        TXT.insert('START', 'Error: ', 'Error')
    return

#==============================================================================
def getOnlineForum():
    try:
        import webbrowser
        TXT.insert('START', 'Opening online documentation.\n')
        webbrowser.open('https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/cassiopee-community')
    except:
        TXT.insert('START', 'Can not open online forum.\n')
        TXT.insert('START', 'Error: ', 'Error')
    return

#==============================================================================
def getOnlineTutorials():
    try:
        import webbrowser
        TXT.insert('START', 'Opening online tutorials.\n')
        webbrowser.open('https://cassiopee.onera.fr/Tutorials/Tutorials.html')
    except:
        TXT.insert('START', 'Can not open online tutorials.\n')
        TXT.insert('START', 'Error: ', 'Error')
    return

#==============================================================================
# Verifie si l'applet est dans les auto-open (PREFS['auto'])
# Retourne False: non, retourne True: oui
#==============================================================================
def isAppAutoOpen(name):
    if 'auto' not in PREFS: return False
    auto = PREFS['auto']
    if auto.find(name) != -1: return True
    else: return False

#==============================================================================
# Interne: appele pour changer le status d'une applet (auto-open ou non)
#==============================================================================
def toggleAutoOpen(name, m):
    opened = isAppAutoOpen(name)
    if opened:
        # Le retire de la liste
        auto = PREFS['auto']
        auto = auto.replace(name+';', '')
        auto = auto.replace(name, '')
        PREFS['auto'] = auto
        m.entryconfigure('Unpin ', image=iconics.PHOTO[9])
        m.entryconfigure('Unpin ', label='Pin ')
    else:
        if 'auto' not in PREFS: PREFS['auto'] = name+';'
        else:
            auto = PREFS['auto']
            auto += name+';'
            PREFS['auto'] = auto
        m.entryconfigure('Pin ', image=iconics.PHOTO[10])
        m.entryconfigure('Pin ', label='Unpin ')
    savePrefFile()

#==============================================================================
# Ajoute a un TkMenu no m un pin menu
#==============================================================================
def addPinMenu(m, name):
    opened = isAppAutoOpen(name)
    if opened: txt='Unpin '; icon = iconics.PHOTO[10]
    else: txt='Pin '; icon = iconics.PHOTO[9]
    m.add_command(
        label=txt, image=icon, compound=TK.RIGHT,
        command=lambda name=name, m=m: toggleAutoOpen(name, m))

#==============================================================================
# Fonction generale d'interface
# Applique la fonction zone a zone
#==============================================================================
def GIF(Function, functionName='myFunction', *args):
    global Nb, Nz
    if t == []: return
    # Ne marche pas sur un arbre temporaire
    if __MAINTREE__ <= 0:
        TXT.insert('START', 'Fail on a temporary tree.\n')
        TXT.insert('START', 'Error: ', 'Error'); return

    # Fonction appliquee a CTK.t
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        TXT.insert('START', 'Selection is empty.\n')
        TXT.insert('START', 'Error: ', 'Error'); return

    # Sauve l'arbre
    saveTree()
    # Applique la fonction zone a zone
    fail = False; errors = []
    for nz in nzs:
        nob = Nb[nz]+1; noz = Nz[nz]
        try:
            a = Function(t[2][nob][2][noz], *args)
            replace(t, nob, noz, a)
        except Exception as e: fail = True; errors += [0, str(e)]
    if not fail: TXT.insert('START', functionName+' succeeds.\n')
    else:
        Panels.displayErrors(errors, header='Error: '+functionName)
        TXT.insert('START', functionName+' fails for at least one zone.\n')
        TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(t)
    (Nb, Nz) = CPlot.updateCPlotNumbering(t)
    TKTREE.updateApp()
    CPlot.render()

    # Log
    import Log
    if Log.LOGGING:
        sel = Log.getSelectedZones()
        if sel is not None:
            Log.LOG += sel
        Log.displayLog()

#===============================================================================
# Mail current image to friends
#===============================================================================
def mail2Friends():
    # Open panel for 'mail'+'message'
    Panels.openMailWindow()
    #Panels.mailData['mailWindow'].master.wait_window(Panels.mailData['mailWindow'])

#===============================================================================
# Save image to document (odt)
#===============================================================================
def save2Doc():
    # Open panel for doc
    Panels.openDocWindow()
    #Panels.docData['docWindow'].master.wait_window(Panels.docData['docWindow'])

#==============================================================================
# Change title : change le titre dans la fenetre CPlot (if top window) + Tk
#==============================================================================
def changeWindowTitle(fileName, filePath="."):
    if fileName == '': return
    CPlot.setWindowTitle(fileName, filePath)
    win = WIDGETS['masterWin']
    win.title('Cassiopee'+C.__version__+' - '+fileName)

# Fonction pour matcher un filter dans une chaine
def matchString(filter, myString):
    if ('*' in filter)|('?' in filter)|('[' in filter):
        ret = fnmatch.fnmatch(myString, filter)
    else:
        ret = False
        if re.search(filter, myString, re.IGNORECASE) is not None: ret = True
    return ret

#==============================================================================
# Meta load function of multiple files
# Si partial: load un CTK.t squelette + HANDLE
# Si full: CTK.t full
# mode = 'full', 'partial', 'auto'
# Cette fonction ne fait pas upgrade et updateApps
#==============================================================================
def tkLoadFile(files, mode='full'):
    global FILE; global HANDLE; global t
    if mode == 'auto':
        try:
            size = 0
            for f in files: size += os.path.getsize(f) # en octets
        except:
            print('Error: convertFile2PyTree: fail to read file %s.'%files[0])
            return
        if size > 1000000000: print('INFO: size: %f Gb'%(size/1000000000))
        elif size > 1000000: print('INFO: size: %f Mb'%(size/1000000))
        else: print('INFO: size: %f kb'%(size/1000))
        maxSize = PREFS.get('maxFileSizeForLoad', 6.) # en Gb
        maxSize = maxSize * 100000000
        if size > maxSize: mode = 'partial'
        else: mode = 'full'

    if mode == 'partial':
        fileName = files[0]
        try: format = Converter.checkFileType(fileName)
        except:
            print('Error: convertFile2PyTree: fail to read file %s.'%fileName)
            return
        if format != 'bin_adf' and format != 'bin_hdf': mode = 'full'

    if mode == 'partial': # partial load
        import Converter.Filter as Filter
        HANDLE = Filter.Handle(files[0])
        t = HANDLE.loadSkeleton()
        HANDLE._loadTreeExtras(t)
        Filter._convert2PartialTree(t)
        HANDLE.getVariables()

    if mode == 'full': # full load of multiple files
        t = []
        for file in files:
            try:
                t2 = C.convertFile2PyTree(file, density=1.)
                if t == []: t = t2
                else: t = C.mergeTrees(t, t2)
            except:
                print('Error: convertFile2PyTree: fail to read file %s.'%file)

    # common part
    FILE = files[0]
