# - fix CAD app -
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import CPlot.iconics as iconics
import math

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def readCAD(event=None):
    import OCC.PyTree as OCC
    fileName = VARS[0].get()
    fileFmt = VARS[1].get()
    if CTK.CADHOOK is not None:
        hook = CTK.CADHOOK
        OCC.freeHook(hook)
        # free the previous hook
        edges = Internal.getNodeFromName1(CTK.t, 'EDGES')
        edges[2] = []
        faces = Internal.getNodeFromName1(CTK.t, 'FACES')
        faces[2] = []
    hook = OCC.readCAD(fileName, fileFmt)
    # Previous hmax, hausd?
    [hmin, hmax, hausd] = OCC.getCADcontainer(CTK.t)
    if hmax is None or (hmax < 0 and hausd < 0):
        (hmax,hmin,hausd) = OCC.occ.analyseEdges(hook)
    OCC._setCADcontainer(CTK.t, fileName, fileFmt, hmin, hmax, hausd)
    CTK.CADHOOK = hook
    # remesh and redisplay
    CTK.setCursor(2, WIDGETS['frame'])
    OCC._meshAllEdges(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd)
    OCC._meshAllFacesTri(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd)
    CTK.setCursor(0, WIDGETS['frame'])
    CTK.display(CTK.t)
    CTK.TXT.insert('START', 'CAD loaded from %s.\n'%fileName)

#==============================================================================
def writeCAD(event=None):
    import OCC.PyTree as OCC
    if CTK.CADHOOK is None: return
    hook = CTK.CADHOOK
    fileName = VARS[0].get()
    fileFmt = VARS[1].get()
    OCC.writeCAD(hook, fileName, fileFmt)
    CTK.TXT.insert('START', 'CAD written in %s.\n'%fileName)

#==============================================================================
def sewCAD(event=None):
    import OCC.PyTree as OCC
    if CTK.CADHOOK is None: return
    tol = CTK.varsFromWidget(VARS[2].get(), 1)[0]
    hook = CTK.CADHOOK
    [hmin, hmax, hausd] = OCC.getCADcontainer(CTK.t)
    faces = []
    nzs = CPlot.getSelectedZones()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        b = CTK.t[2][nob]
        if b[0] != 'FACES': continue
        z = CTK.t[2][nob][2][noz]
        CAD = Internal.getNodeFromName1(z, 'CAD')
        if CAD is not None:
            no = Internal.getNodeFromName1(CAD, 'no')
            no = Internal.getValue(no)
            faces.append(no)

    CTK.setCursor(2, WIDGETS['frame'])
    CTK.setCursor(2, WIDGETS['sewingButton'])

    OCC._sewing(hook, faces, tol)
    # remesh CAD and redisplay
    edges = Internal.getNodeFromName1(CTK.t, 'EDGES')
    edges[2] = []
    faces = Internal.getNodeFromName1(CTK.t, 'FACES')
    faces[2] = []
    OCC._meshAllEdges(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd) # loose manual remeshing
    OCC._meshAllFacesTri(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd)

    CTK.setCursor(0, WIDGETS['frame'])
    CTK.setCursor(0, WIDGETS['sewingButton'])

    NL = OCC.getNbLonelyEdges(CTK.t)
    VARS[4].set('Lonely edges: %d'%NL)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    CTK.TXT.insert('START', 'CAD sewed with %g.\n'%tol)

#==============================================================================
def filletCAD(event=None):
    import OCC.PyTree as OCC
    if CTK.CADHOOK is None: return
    radius = CTK.varsFromWidget(VARS[3].get(), 1)[0]
    hook = CTK.CADHOOK
    [hmin, hmax, hausd] = OCC.getCADcontainer(CTK.t)
    # Get selected edges
    nzs = CPlot.getSelectedZones()
    edges = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        b = CTK.t[2][nob]
        if b[0] != 'EDGES': continue
        z = CTK.t[2][nob][2][noz]
        CAD = Internal.getNodeFromName1(z, 'CAD')
        if CAD is not None:
            no = Internal.getNodeFromName1(CAD, 'no')
            no = Internal.getValue(no)
            edges.append(no)
    if edges == []:
        CTK.TXT.insert('START', 'No valid edges in selection.\n')
        return

    CTK.setCursor(2, WIDGETS['frame'])
    CTK.setCursor(2, WIDGETS['filletButton'])

    try:
        OCC._addFillet(hook, edges, radius)
    except:
        CTK.setCursor(0, WIDGETS['frame'])
        CTK.setCursor(0, WIDGETS['filletButton'])
        CTK.TXT.insert('START', 'Add fillet fails.\n')
        return

    # remesh CAD and redisplay
    edges = Internal.getNodeFromName1(CTK.t, 'EDGES')
    edges[2] = []
    faces = Internal.getNodeFromName1(CTK.t, 'FACES')
    faces[2] = []
    OCC._meshAllEdges(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd) # loose manual remeshing...
    OCC._meshAllFacesTri(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd)
    CTK.setCursor(0, WIDGETS['frame'])
    CTK.setCursor(0, WIDGETS['filletButton'])

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    CTK.TXT.insert('START', 'Fillet added.\n')

#==============================================================================
def removeFaces(event=None):
    import OCC.PyTree as OCC
    if CTK.CADHOOK is None: return
    hook = CTK.CADHOOK
    #[hmin, hmax, hausd] = OCC.getCADcontainer(CTK.t)
    # Get selected faces
    nzs = CPlot.getSelectedZones()
    faces = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        b = CTK.t[2][nob]
        if b[0] != 'FACES': continue
        z = CTK.t[2][nob][2][noz]
        CAD = Internal.getNodeFromName1(z, 'CAD')
        if CAD is not None:
            no = Internal.getNodeFromName1(CAD, 'no')
            no = Internal.getValue(no)
            faces.append(no)
    if faces == []:
        CTK.TXT.insert('START', 'No valid faces in selection.\n')
        return

    CTK.setCursor(2, WIDGETS['frame'])
    CTK.setCursor(2, WIDGETS['removeFacesButton'])

    # old style (full remesh)
    #edgeMap = []; faceMap = []
    #OCC._removeFaces(hook, faces, edgeMap, faceMap)
    #edges = Internal.getNodeFromName1(CTK.t, 'EDGES')
    #edges[2] = []
    #faces = Internal.getNodeFromName1(CTK.t, 'FACES')
    #faces[2] = []
    #OCC._meshAllEdges(hook, CTK.t, hmax=hmax, hausd=hausd)
    #OCC._meshAllFacesTri(hook, CTK.t, hmax=hmax, hausd=hausd)

    # new style (no remesh)
    nbEdges = OCC.getNbEdges(hook)
    nbFaces = OCC.getNbFaces(hook)
    new2OldEdgeMap = []; new2OldFaceMap = []
    OCC._removeFaces(hook, faces, new2OldEdgeMap, new2OldFaceMap)
    OCC._updateTree(CTK.t, nbEdges, nbFaces, new2OldEdgeMap, new2OldFaceMap)

    CTK.setCursor(0, WIDGETS['frame'])
    CTK.setCursor(0, WIDGETS['removeFacesButton'])
    NL = OCC.getNbLonelyEdges(CTK.t)
    VARS[4].set('Lonely edges: %d'%NL)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

    CTK.TXT.insert('START', 'Faces removed from CAD.\n')

#==============================================================================
# filter lonelyEdges based on color
def filterLonelyEdges(event=None):
    EDGES = Internal.getNodeFromName1(CTK.t, 'EDGES')
    if EDGES is None: return
    FACES = Internal.getNodeFromName1(CTK.t, 'FACES')
    activeMode = 1
    activeList = []
    if FACES is not None:
        zones = Internal.getZones(FACES)
        if len(zones) > 0:
            noz = CPlot.getCPlotNumber(CTK.t, 'FACES', zones[0][0])
            activeMode = CPlot.getActiveStatus(noz)

        if activeMode == 0: # previously deactivated
            for z in zones:
                noz = CPlot.getCPlotNumber(CTK.t, 'FACES', z[0])
                #CPlot.setActiveZones([(noz,1)])
                activeList.append((noz,1))
        else:
            for z in zones:
                noz = CPlot.getCPlotNumber(CTK.t, 'FACES', z[0])
                #CPlot.setActiveZones([(noz,0)])
                activeList.append((noz,0))

    if activeMode == 1:
        for z in Internal.getZones(EDGES):
            noz = CPlot.getCPlotNumber(CTK.t, 'EDGES', z[0])
            p = Internal.getNodeFromPath(z, '.RenderInfo/Color')
            active = 0
            if p is not None:
                if Internal.getValue(p) == 'Red': active = 1
            #CPlot.setActiveZones([(noz,active)])
            activeList.append((noz,active))

    else:
        for z in Internal.getZones(EDGES):
            noz = CPlot.getCPlotNumber(CTK.t, 'EDGES', z[0])
            #CPlot.setActiveZones([(noz,1)])
            activeList.append((noz,1))

    CPlot.setActiveZones(activeList)

    CTK.TXT.insert('START', 'Lonely edges displayed.\n')

#==============================================================================
def fillHole(event=None):
    import OCC.PyTree as OCC
    if CTK.CADHOOK is None: return
    hook = CTK.CADHOOK
    [hmin, hmax, hausd] = OCC.getCADcontainer(CTK.t)
    # Get selected edges
    nzs = CPlot.getSelectedZones()
    edges = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        b = CTK.t[2][nob]
        if b[0] != 'EDGES': continue
        z = CTK.t[2][nob][2][noz]
        CAD = Internal.getNodeFromName1(z, 'CAD')
        if CAD is not None: edges.append(z)
    if edges == []:
        CTK.TXT.insert('START', 'No valid edges in selection.\n')
        return

    CTK.setCursor(2, WIDGETS['frame'])
    CTK.setCursor(2, WIDGETS['fillHoleButton'])

    edges = OCC.orderEdgeList(edges)
    #print('edgeList', edges, flush=True)
    try:
        OCC._fillHole(hook, edges)
    except:
        CTK.setCursor(0, WIDGETS['frame'])
        CTK.setCursor(0, WIDGETS['fillHoleButton'])
        CTK.TXT.insert('START', 'Fill hole fails.\n')
        return

    # remesh CAD and redisplay
    edges = Internal.getNodeFromName1(CTK.t, 'EDGES')
    edges[2] = []
    faces = Internal.getNodeFromName1(CTK.t, 'FACES')
    faces[2] = []
    OCC._meshAllEdges(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd) # loose manual remeshing...
    OCC._meshAllFacesTri(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd)
    CTK.setCursor(0, WIDGETS['frame'])
    CTK.setCursor(0, WIDGETS['fillHoleButton'])

    NL = OCC.getNbLonelyEdges(CTK.t)
    VARS[4].set('Lonely edges: %d'%NL)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    CTK.TXT.insert('START', 'Fill hole in CAD.\n')

#==============================================================================
def setTrimFace1():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    selected = ''
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        b = CTK.t[2][nob]
        if b[0] != 'FACES': continue
        z = CTK.t[2][nob][2][noz]
        CAD = Internal.getNodeFromName1(z, 'CAD')
        if CAD is not None:
            no = Internal.getNodeFromName1(CAD, 'no')
            no = Internal.getValue(no)
            selected += str(no)+' '
    #print(selected, flush=True)
    VARS[5].set(selected)

#==============================================================================
def trimFaces(event=None):
    import OCC.PyTree as OCC
    if CTK.CADHOOK is None: return
    hook = CTK.CADHOOK
    [hmin, hmax, hausd] = OCC.getCADcontainer(CTK.t)
    # Get selected faces
    nzs = CPlot.getSelectedZones()
    faces1 = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        b = CTK.t[2][nob]
        if b[0] != 'FACES': continue
        z = CTK.t[2][nob][2][noz]
        CAD = Internal.getNodeFromName1(z, 'CAD')
        if CAD is not None:
            no = Internal.getNodeFromName1(CAD, 'no')
            no = Internal.getValue(no)
            faces1.append(no)
    if faces1 == []:
        CTK.TXT.insert('START', 'No valid faces in selection.\n')
        return

    faces2 = []
    st = VARS[5].get()
    st = st.split(' ')
    for s in st:
        try:
            val = int(s)
            faces2.append(val)
        except: pass
    if faces2 == []:
        CTK.TXT.insert('START', 'No valid faces in selection.\n')
        return

    CTK.setCursor(2, WIDGETS['frame'])
    CTK.setCursor(2, WIDGETS['trimFacesButton'])

    OCC.occ.trimFaces(hook, faces1, faces2)

    # remesh CAD and redisplay
    edges = Internal.getNodeFromName1(CTK.t, 'EDGES')
    edges[2] = []
    faces = Internal.getNodeFromName1(CTK.t, 'FACES')
    faces[2] = []
    OCC._meshAllEdges(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd) # loose manual remeshing...
    OCC._meshAllFacesTri(hook, CTK.t, hmin=hmin, hmax=hmax, hausd=hausd)

    CTK.setCursor(0, WIDGETS['frame'])
    CTK.setCursor(0, WIDGETS['trimFacesButton'])
    NL = OCC.getNbLonelyEdges(CTK.t)
    VARS[4].set('Lonely edges: %d'%NL)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

    CTK.TXT.insert('START', 'Faces trimmed in CAD.\n')

#==============================================================================
def checkWatertight(event=None):
    import OCC.PyTree as OCC
    if CTK.t == []: return
    b = Internal.getNodeFromName1(CTK.t, 'FACES')
    if b is None:
        CTK.TXT.insert('START', 'No FACES base. Check not performed.\n')
        return
    import Transform.PyTree as T
    import Post.PyTree as P
    import Generator.PyTree as G

    CTK.setCursor(2, WIDGETS['frame'])
    CTK.setCursor(2, WIDGETS['checkWatertight'])

    CTK.saveTree()
    p = Internal.getNodeFromName1(CTK.t, 'LEAKS')
    if p is not None:
        zones = []
        for z in Internal.getZones(p): zones.append('LEAKS'+Internal.SEP1+z[0])
        CPlot.delete(zones)
        gnob = C.getNobOfBase(p, CTK.t)
        del CTK.t[2][gnob]

    tol = 1.e-6
    CAD = Internal.getNodeFromName1(CTK.t, 'CAD')
    if CAD is not None:
        hmax = Internal.getNodeFromName1(CAD, 'hmax')
        hmax = Internal.getValue(hmax)
        hmin = Internal.getNodeFromName1(CAD, 'hmin')
        hmin = Internal.getValue(hmin)
        tol = (hmin+hmax)*0.2

    CTK.t = C.addBase2PyTree(CTK.t, 'LEAKS', 1)
    p = Internal.getNodeFromName1(CTK.t, 'LEAKS')
    gnob = C.getNobOfBase(p, CTK.t)

    ef = OCC.getComponents(CTK.t)

    VARS[6].set('Components: %d'%(len(ef)))

    isWatertight = False
    #try:
    isWatertight = True
    for f in ef:
        try:
            ext = P.exteriorFaces(f)
            ext = T.splitConnexity(ext)
            for i in ext: CTK.add(CTK.t, gnob, -1, i)
            if len(ext) != 0: isWatertight = False; break
        except: isWatertight = True

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

    CTK.setCursor(0, WIDGETS['frame'])
    CTK.setCursor(0, WIDGETS['checkWatertight'])

    if  isWatertight: TTK.setButtonGreen(WIDGETS['checkWatertight'])
    else: TTK.setButtonRed(WIDGETS['checkWatertight'])

    if isWatertight: CTK.TXT.insert('START', 'CAD is watertight with tol %g.\n'%tol)
    else: CTK.TXT.insert('START', 'CAD is not watertight with tol %g.\n'%tol)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCADFix  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Fix CAD.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkCADFix')
    WIDGETS['frameMenu'] = FrameMenu

    #- VARS -
    tol = 1.e-6; NL = 0
    fileName = ''; fileFmt = 'fmt_step'
    if CTK.CADHOOK is not None:
        import OCC.PyTree as OCC
        fileName, fileFmt = OCC.getFileAndFormat(CTK.CADHOOK)
        CAD = Internal.getNodeFromName1(CTK.t, 'CAD')
        if CAD is not None:
            NL = OCC.getNbLonelyEdges(CTK.t)
            tol = Internal.getNodeFromName1(CAD, 'hmax')
            tol = Internal.getValue(tol)
            power = math.floor(math.log10(abs(tol)))
            tol = round(tol, -power)

    # -0- CAD file name -
    V = TK.StringVar(win); V.set(fileName); VARS.append(V)
    # -1- CAD file format -
    V = TK.StringVar(win); V.set(fileFmt); VARS.append(V)
    # -2- Sewing tolerance -
    V = TK.StringVar(win); V.set(str(tol)); VARS.append(V)
    # -3- Fillet radius -
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)
    # -4- Bilan des edges lonely
    V = TK.StringVar(win); V.set('Lonely edges: %d'%NL); VARS.append(V)
    # -5- List of face to trim (compound 1)
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -6- Number of components
    V = TK.StringVar(win); V.set('Components: 0'); VARS.append(V)

    # CAD file name
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=15)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='CAD file name.')

    # CAD file format
    B = TTK.OptionMenu(Frame, VARS[1], 'fmt_step', 'fmt_iges')
    B.grid(row=0, column=1, sticky=TK.EW)

    # Read/write CAD file
    B = TTK.Button(Frame, text="Read", command=readCAD)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Read CAD in tree.')

    B = TTK.Button(Frame, text="Write", command=writeCAD)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Write CAD to file.')

    # Lonely edges
    B = TTK.Label(Frame, textvariable=VARS[4])
    B.bind("<Button-1>", filterLonelyEdges)
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of lonely edges.')

    # Number of components
    B = TTK.Label(Frame, textvariable=VARS[6])
    B.grid(row=2, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of connex watertight components.')

    # Check watertight
    B = TTK.Button(Frame, text="Check watertight", command=checkWatertight)
    B.grid(row=3, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Check if CAD is watertight.')
    WIDGETS['checkWatertight'] = B

    # Sewing
    B = TTK.Button(Frame, text="Sew", command=sewCAD)
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Sew CAD to fix multiple edges.')
    WIDGETS['sewingButton'] = B

    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=10)
    B.grid(row=4, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Sewing tolerance.')
    B.bind('<Return>', sewCAD)

    # Fillet
    B = TTK.Button(Frame, text="Fillet", command=filletCAD)
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Make a fillet from edges selection.')
    WIDGETS['filletButton'] = B

    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=10)
    B.grid(row=5, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fillet radius.')
    B.bind('<Return>', filletCAD)

    # Remove faces
    B = TTK.Button(Frame, text="Remove faces", command=removeFaces)
    B.grid(row=6, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Remove selected faces from CAD.')
    WIDGETS['removeFacesButton'] = B

    # Fill hole
    B = TTK.Button(Frame, text="Fill hole", command=fillHole)
    B.grid(row=6, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fill hole from CAD edges.')
    WIDGETS['fillHoleButton'] = B

    # Trim
    B = TTK.Button(Frame, text="Trim faces", command=trimFaces)
    B.grid(row=7, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Trim faces.')
    WIDGETS['trimFacesButton'] = B
    B = TTK.Button(Frame, command=setTrimFace1,
                   image=iconics.PHOTO[8], padx=0, pady=0)
    BB = CTK.infoBulle(parent=B, text='Set the faces for first compound.')
    B.grid(row=7, column=1, columnspan=1, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    try: CTK.WIDGETS['TreeNoteBook'].add(WIDGETS['frame'], text='tkCADFix')
    except: pass
    CTK.WIDGETS['TreeNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    CTK.WIDGETS['TreeNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkCADFix '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
