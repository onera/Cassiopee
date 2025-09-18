# - tkFixer2 -
"""Fix holes in meshes."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T
import Intersector.PyTree as XOR
import Post.PyTree as P
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Fix one gap
#==============================================================================
def fixGap():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # get bump factor
    factor = 2*WIDGETS['bump'].get() / 100.

    # Contours
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    contours = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        contours.append(CTK.t[2][nob][2][noz])

    contours = C.convertArray2Tetra(contours)
    contours = T.join(contours)
    contours = G.close(contours)
    contours = T.splitManifold(contours)

    fail = False; out = []
    for contour in contours:
        try:
            p = G.fittingPlaster(contour, bumpFactor=factor)
            b = G.gapfixer(contour, p)
            out.append(b)
        except Exception as e:
            fail = True
            Panels.displayErrors([0,str(e)], header='Error: gapfixer on %s.'%contour[0])

    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'FIXED', 2)
    base = Internal.getNodesFromName1(CTK.t, 'FIXED')[0]
    nob = C.getNobOfBase(base, CTK.t)
    for b in out:
        CTK.add(CTK.t, nob, -1, b)
        #CTK.add(CTK.t, nob, -1, p)
    if not fail:
        CTK.TXT.insert('START', 'Gap fixed.\n')
    else:
        CTK.TXT.insert('START', 'Gap fixing fails.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Fix multiple gaps, eliminate overlap in patches
#==============================================================================
def fixGaps():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # get mode
    modeString = VARS[0].get()
    mode = 2; coplanar = 0 # unknown
    if modeString == 'Centers': mode = 1; coplanar = 0
    elif modeString == 'Nodes': mode = 0; coplanar = 0
    elif modeString == 'Slice': mode = 0; coplanar = 1

    # Patches
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    patches = []
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        patches.append(CTK.t[2][nob][2][noz])
    nob0 = CTK.Nb[nzs[0]]+1
    patches = C.convertArray2Tetra(patches)
    patches = G.close(patches)
    b = G.gapsmanager(patches, mode=mode, coplanar=coplanar)
    CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
    CTK.t[2][nob0][2] += b

    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Gaps fixed.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def conformUnstr():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    tol = CTK.varsFromWidget(VARS[1].get(), type=1)
    if len(tol) != 1:
        CTK.TXT.insert('START', 'Tolerance is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = tol[0]

    split = VARS[2].get()
    if split == '1': split = True
    else: split = False

    CTK.setCursor(2, WIDGETS['conformUnstr'])
    C._fillMissingVariables(CTK.t)
    patches = []
    CTK.saveTree()

    if len(nzs) == 2: # two by two
        nz = nzs[0]
        nob1 = CTK.Nb[nz]+1
        noz1 = CTK.Nz[nz]
        z1 = CTK.t[2][nob1][2][noz1]
        nz = nzs[1]
        nob2 = CTK.Nb[nz]+1
        noz2 = CTK.Nz[nz]
        z2 = CTK.t[2][nob2][2][noz2]
        if split:
            #zo1 = XOR.conformUnstr(z1, z2, tol=tol, itermax=1)
            #zo2 = XOR.conformUnstr(z2, z1, tol=tol, itermax=1)
            z = T.join([z1, z2])
            try:
                z = XOR.conformUnstr(z, tol=tol, itermax=1)
            except Exception as e:
                CTK.setCursor(0, WIDGETS['conformUnstr'])
                Panels.displayErrors([0,str(e)], header='Error: conformUnstr')
                CTK.TXT.insert('START', 'ConformUnstr failed\n'); return

            zones = T.splitManifold(z); lz = len(zones)
            if lz > 0: CTK.replace(CTK.t, nob1, noz1, zones[0])
            if lz > 1: CTK.replace(CTK.t, nob2, noz2, zones[1])
            for i in zones[2:]: CTK.add(CTK.t, nob1, -1, i)
        else:
            try:
                zo1 = XOR.conformUnstr(z1, z2, tol=tol, itermax=1)
                zo2 = XOR.conformUnstr(z2, z1, tol=tol, itermax=1)
            except Exception as e:
                CTK.setCursor(0, WIDGETS['conformUnstr'])
                Panels.displayErrors([0,str(e)], header='Error: conformUnstr')
                CTK.TXT.insert('START', 'ConformUnstr failed\n'); return

            CTK.replace(CTK.t, nob1, noz1, zo1)
            CTK.replace(CTK.t, nob2, noz2, zo2)
    else:
        if split:
            zones = []
            for nz in nzs:
                nob = CTK.Nb[nz]+1
                noz = CTK.Nz[nz]
                zones.append(CTK.t[2][nob][2][noz])
            z = T.join(zones)
            z = XOR.conformUnstr(z, tol=tol, itermax=1)
            zones = T.splitManifold(z); lz = len(zones)
            if lz <= len(nzs):
                for i in range(lz):
                    nz = nzs[i]
                    nob = CTK.Nb[nz]+1
                    noz = CTK.Nz[nz]
                    CTK.replace(CTK.t, nob, noz, zones[i])
            nob = CTK.Nb[0]+1
            for i in zones[lz:]: CTK.add(CTK.t, nob, -1, i)
        else:
            for nz in nzs:
                nob = CTK.Nb[nz]+1
                noz = CTK.Nz[nz]
                z = CTK.t[2][nob][2][noz]
                z = XOR.conformUnstr(z, tol=tol, itermax=1)
                CTK.replace(CTK.t, nob, noz, z)

    CTK.TXT.insert('START', 'Surface conformized.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['conformUnstr'])

#==============================================================================
def forceMatch():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    tol = CTK.varsFromWidget(VARS[3].get(), type=1)
    if len(tol) != 1:
        CTK.TXT.insert('START', 'Tolerance is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = tol[0]

    #P1 = CTK.varsFromWidget(VARS[4].get(), type=1)
    #P2 = CTK.varsFromWidget(VARS[5].get(), type=1)

    #if len(nzs) != 2:
    #    CTK.TXT.insert('START', 'Need two patches.\n')
    #    CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.setCursor(2, WIDGETS['forceMatch'])
    CTK.saveTree()

    nz = nzs[0]
    nob1 = CTK.Nb[nz]+1
    noz1 = CTK.Nz[nz]
    z1 = CTK.t[2][nob1][2][noz1]
    if len(nzs) > 1:
        nz = nzs[1]
        nob2 = CTK.Nb[nz]+1
        noz2 = CTK.Nz[nz]
        z2 = CTK.t[2][nob2][2][noz2]
    else: z2 = None

    #z1 = []
    #name = VARS[6].get()
    #names = name.split(';')
    #for v in names:
    #    v = v.lstrip(); v = v.rstrip()
    #    sname = v.split('/', 1)
    #    base = Internal.getNodeFromName1(CTK.t, sname[0])
    #    if base is not None:
    #        nodes = Internal.getNodesFromType1(base, 'Zone_t')
    #        for z in nodes:
    #            if z[0] == sname[1]: z1.append(z)
    #if len(z1) > 0: z1 = z1[0]
    #else: z1 = None

    c1 = []
    name = VARS[7].get()
    names = name.split(';')
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        base = Internal.getNodeFromName1(CTK.t, sname[0])
        if base is not None:
            nodes = Internal.getNodesFromType1(base, 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: c1.append(z)
    if c1 != []: c1 = T.join(c1)

    #z2 = []
    #name = VARS[8].get()
    #names = name.split(';')
    #for v in names:
    #    v = v.lstrip(); v = v.rstrip()
    #    sname = v.split('/', 1)
    #    base = Internal.getNodeFromName1(CTK.t, sname[0])
    #    if base is not None:
    #        nodes = Internal.getNodesFromType1(base, 'Zone_t')
    #        for z in nodes:
    #            if z[0] == sname[1]: z2.append(z)
    #if len(z2) > 0: z2 = z2[0]
    #else: z2 = None

    c2 = []
    name = VARS[9].get()
    names = name.split(';')
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        base = Internal.getNodeFromName1(CTK.t, sname[0])
        if base is not None:
            nodes = Internal.getNodesFromType1(base, 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: c2.append(z)
    if c2 != []: c2 = T.join(c2)

    if c1 == [] and c2 == []:
        G._close(z1, tol)
        if z2 is not None: G._close(z2, tol)
    elif c1 != [] and c2 != []:
        G._forceMatch(z1, z2, C1=c1, C2=c2)
    else:
        G._forceMatch(z1, z2, tol=tol)

    # Try merge
    if z2 is not None: z = T.join(z1, z2)
    else: z = z1
    z = G.close(z)

    z = T.reorder(z, (1,))
    CTK.replace(CTK.t, nob1, noz1, z)
    if z2 is not None:
        CPlot.delete([CTK.t[2][nob2][0]+Internal.SEP1+z2[0]])
        del CTK.t[2][nob2][2][noz2]

    if VARS[10].get() == '1': displayFix()

    CTK.TXT.insert('START', 'Surfaces match.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['forceMatch'])

#==============================================================================
def setPointCoordinates(V):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    point = CPlot.getActivePoint()
    if point != []:
        V.set(str(point[0])+';'+str(point[1])+';'+str(point[2]))

# V: VARS, W: Widget
def setSel(V, W):
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
        z = CTK.t[2][nob][2][noz]
        selected += CTK.t[2][nob][0]+'/'+z[0]+';'
    selected = selected[0:-1]
    V.set(selected); W.xview(len(selected)-1)

# Cree une base avec les courbes a fixer
# a partir de toutes les zones TRI de t
def getFix(t):
    zones = Internal.getZones(t)
    zonesF = []
    for z in zones:
        dim = Internal.getZoneDim(z)
        if dim[3] == 'TRI': zonesF.append(z)
    try: ext = P.exteriorFaces(zonesF)
    except: ext = []
    if ext != []:
        ext = T.splitManifold(ext)
        #ext = T.splitSharpEdges(ext)
    return ext

def displayFix(event=None):
    if CTK.t == []: return
    # Find base
    b = Internal.getNodeFromName1(CTK.t, 'TOFIX')
    if b is None: b = Internal.newCGNSBase('TOFIX', parent=CTK.t)
    # delete previous zones
    size = len(b[2])
    for c in range(size):
        z = b[2][size-c-1]
        if z[3] == 'Zone_t':
            CPlot.delete([b[0]+Internal.SEP1+z[0]])
        del b[2][size-c-1]
    # add new zones
    ext = getFix(CTK.t)
    nob = C.getNobOfBase(b, CTK.t)
    for z in ext: CTK.add(CTK.t, nob, -1, z)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    return None

def removeFix(event=None):
    # Find base
    b = Internal.getNodeFromName1(CTK.t, 'TOFIX')
    if b is None: b = Internal.newCGNSBase('TOFIX', parent=CTK.t)
    # delete previous zones
    size = len(b[2])
    for c in range(size):
        z = b[2][size-c-1]
        if z[3] == 'Zone_t':
            CPlot.delete([b[0]+Internal.SEP1+z[0]])
        del b[2][size-c-1]
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    return None

def toggleFix(event=None):
    if VARS[10].get() == '0': removeFix()
    else: displayFix()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFixer2  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Fix holes in surfaces.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=0)
    Frame.columnconfigure(2, weight=1)
    Frame.columnconfigure(3, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkFixer2')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- mode pour gapsmanager -
    V = TK.StringVar(win); V.set('Nodes'); VARS.append(V)
    # -1- tol pour conformUnstr -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -2- checkbox if split in conformization
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -3- tolerance for forcematch
    V = TK.StringVar(win); V.set('1.e-6'); VARS.append(V)
    # -4- point 1 coordinates -
    V = TK.StringVar(win); V.set('0;0;0'); VARS.append(V)
    # -5- point 2 coordinates -
    V = TK.StringVar(win); V.set('0;0;0'); VARS.append(V)
    # -6- Bloc 1 -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -7- Courbe 1 -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -8- Bloc 2 -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -9- Courbe 2 -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -10- fix mode : display curves to be fixed
    V = TK.StringVar(win); V.set('0'); VARS.append(V)

    # - Slider -
    B = TTK.Scale(Frame, from_=-50, to=50, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, value=0)
    WIDGETS['bump'] = B
    B.grid(row=0, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Bump.')

    # - Fix gap in contour -
    B = TTK.Button(Frame, text="Fix gap in contour", command=fixGap)
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fix gap in a contour.')

    # - Fix gaps in a set of patches (with overlap) -
    B = TTK.Button(Frame, text="Fix gaps in patches", command=fixGaps)
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fix gaps in patches (even with overlap).')
    B = TTK.OptionMenu(Frame, VARS[0], 'Nodes', 'Centers', 'Unknown', 'Slice')
    B.grid(row=1, column=2, columnspan=2, sticky=TK.EW)

    # - conformUnstr -
    B = TTK.Button(Frame, text="conformUnstr", command=conformUnstr)
    WIDGETS['conformUnstr'] = B
    B.grid(row=2, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Conformize a TRI surface.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Tolerance used in conformUnstr operations.\n0. means automatic setting.')
    B = TTK.Checkbutton(Frame, text='', variable=VARS[2])
    B.grid(row=2, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split while conformizing.')

    # - forceMatch -
    #B = TTK.Entry(Frame, textvariable=VARS[6], background='White', width=15)
    #B.grid(row=3, column=0, columnspan=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Source1.')
    #B = TTK.Button(Frame, image=iconics.PHOTO[8],
    #               command=lambda x=VARS[6] : setSel(x), padx=0)
    #B.grid(row=3, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[7], background='White', width=15)
    B.grid(row=3, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Courbe1.')
    WIDGETS['courbe1'] = B
    B = TTK.Button(Frame, image=iconics.PHOTO[8],
                   command=lambda x=VARS[7], y=WIDGETS['courbe1'] : setSel(x, y), padx=0)
    B.grid(row=3, column=1, sticky=TK.EW)

    #B = TTK.Entry(Frame, textvariable=VARS[8], background='White', width=15)
    #B.grid(row=4, column=0, columnspan=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Source1.')
    #B = TTK.Button(Frame, image=iconics.PHOTO[8],
    #               command=lambda x=VARS[8] : setSel(x), padx=0)
    #B.grid(row=4, column=2, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[9], background='White', width=15)
    B.grid(row=3, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Courbe2.')
    WIDGETS['courbe2'] = B
    B = TTK.Button(Frame, image=iconics.PHOTO[8],
                   command=lambda x=VARS[9],y=WIDGETS['courbe2'] : setSel(x,y), padx=0)
    B.grid(row=3, column=3, sticky=TK.EW)

    B = TTK.Button(Frame, text="forceMatch", command=forceMatch)
    WIDGETS['forceMatch'] = B
    B.grid(row=4, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Force two patch to match.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=5)
    B.grid(row=4, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Tolerance used in forceMatch operations.')
    B = TTK.Checkbutton(Frame, text='', variable=VARS[10], command=toggleFix)
    BB = CTK.infoBulle(parent=B, text='Show curves to fix.')
    B.grid(row=4, column=3, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkFixer2')
    except: pass
    CTK.WIDGETS['SurfNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['SurfNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkFixer '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
