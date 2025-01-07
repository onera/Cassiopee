# - tkExtrusion -
"""Generate mesh by extrusion."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Generator.PyTree as G
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Geom.PyTree as D
import Converter.Internal as Internal
import Transform.PyTree as T
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#=============================================================================
def revolve():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    # - axis -
    name = VARS[4].get()
    names = name.split(';')
    curve = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: curve.append(z)
    if len(curve) == 0:
        CTK.TXT.insert('START', 'Axis is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    teta = CTK.varsFromWidget(VARS[5].get(), type=1)
    if len(teta) != 1:
        CTK.TXT.insert('START', 'Revolve angle is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    teta = teta[0]

    Nteta = CTK.varsFromWidget(VARS[6].get(), type=2)
    if len(Nteta) != 1:
        CTK.TXT.insert('START', 'Number of points for revolve is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    Nteta = Nteta[0]

    # - Extrait 2 pts de l'axe -
    axis = curve[0]
    [xo,yo,zo] = C.getValue(axis, Internal.__GridCoordinates__, 0)

    dim = Internal.getZoneDim(axis)
    if dim[0] == 'Structured': np = dim[1]*dim[2]*dim[3]
    else: np = dim[1]
    [x1,y1,z1] = C.getValue(axis, Internal.__GridCoordinates__, np-1)

    ntx = x1-xo; nty = y1-yo; ntz = z1-zo

    # - Extrait rmod si il y en a
    if len(curve) > 1: rmod = curve[1]
    else: rmod = None

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            z = D.axisym(z, (xo,yo,zo), (ntx,nty,ntz), teta, Nteta, rmod)
            CTK.replace(CTK.t, nob, noz, z)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Mesh revolved with axis.\n')
    else:
        Panels.displayErrors(errors, header='Error: revolve')
        CTK.TXT.insert('START', 'Revolution fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def setAxis():
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
    VARS[4].set(selected)

#==============================================================================
def setCurve():
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
    VARS[3].set(selected)

#=============================================================================
def lineDrive():
    extrudeWCurve(mode=0)

def orthoDrive():
    extrudeWCurve(mode=1)

#=============================================================================
def extrudeWCurve(mode=0):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    # - curve -
    name = VARS[3].get()
    names = name.split(';')
    curve = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: curve.append(z)
    if len(curve) == 0:
        CTK.TXT.insert('START', 'Curve is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            if mode == 0: z = D.lineDrive(z, curve)
            else: z = D.orthoDrive(z, curve)
            CTK.replace(CTK.t, nob, noz, z)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Mesh extruded with curve(s).\n')
    else:
        Panels.displayErrors(errors, header='Error: extrusion')
        CTK.TXT.insert('START', 'Extrusion fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# AddNormalLayers
#==============================================================================
def addLayers():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    h = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(h) != 1:
        CTK.TXT.insert('START', 'Layer height is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    h = h[0]
    N = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(N) != 1:
        CTK.TXT.insert('START', 'Number of layers is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    N = N[0]
    smooth = CTK.varsFromWidget(VARS[2].get(), type=2)
    if len(smooth) != 1:
        CTK.TXT.insert('START', 'Smoothing iterations number is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    smooth = smooth[0]

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.setCursor(2, WIDGETS['addLayers'])
    CTK.saveTree()
    zlist = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        zlist.append(z)

    d = G.cart((0,0,0), (h,1,1), (N+1,1,1))
    fail = False; errors = []
    try:
        zlist = G.addNormalLayers(zlist, d, niter=smooth)
    except Exception as e:
        fail = True; errors += [0,str(e)]

    for z in zlist: z[0] = C.getZoneName(z[0]) # unique name

    CTK.t = C.addBase2PyTree(CTK.t, 'MESHES')
    base = Internal.getNodeFromName1(CTK.t, 'MESHES')

    if not fail:
        nob = C.getNobOfBase(base, CTK.t)
        for i in zlist: CTK.add(CTK.t, nob, -1, i)
        CTK.TXT.insert('START', 'Normal layers added.\n')
    else:
        Panels.displayErrors(errors, header='Error: addNormalLayers')
        CTK.TXT.insert('START', 'Add normal layers fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['addLayers'])

#==============================================================================
def addkplanes():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    N = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(N) != 1:
        CTK.TXT.insert('START', 'Number of layers is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    N = N[0]
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    fail = False
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        #try:
        z = T.addkplane(CTK.t[2][nob][2][noz], N=N)
        CTK.replace(CTK.t, nob, noz, z)
        #except Exception as e: fail = True
    if not fail: CTK.TXT.insert('START', 'K planes added.\n')
    else:
        CTK.TXT.insert('START', 'add K planes failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkExtrusion  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Extrude a mesh.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkExtrusion')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Hauteur de chaque maille -
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)
    if 'tkExtrusionHeight' in CTK.PREFS:
        V.set(CTK.PREFS['tkExtrusionHeight'])
    # -1- Nombre de layers a ajouter
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    if 'tkExtrusionNLayers' in CTK.PREFS:
        V.set(CTK.PREFS['tkExtrusionNLayers'])
    # -2- Nombre d'iterations de lissage
    V = TK.StringVar(win); V.set('50'); VARS.append(V)
    if 'tkExtrusionSmooth' in CTK.PREFS:
        V.set(CTK.PREFS['tkExtrusionSmooth'])
    # -3- Driving curve
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -4- Axis for revolve
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -5- Angle for revolve
    V = TK.StringVar(win); V.set('360.'); VARS.append(V)
    if 'tkExtrusionRevAngle' in CTK.PREFS:
        V.set(CTK.PREFS['tkExtrusionRevAngle'])
    # -6- Npts for revolve
    V = TK.StringVar(win); V.set('30'); VARS.append(V)
    if 'tkExtrusionNpts' in CTK.PREFS:
        V.set(CTK.PREFS['tkExtrusionNpts'])

    # - AddNormalLayers -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Height of each layer.')

    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of layers.')

    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=0, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of smoothing iterations.')

    B = TTK.Button(Frame, text="Add layers", command=addLayers)
    WIDGETS['addLayers'] = B
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Add normal layers to surface.')

    B = TTK.Button(Frame, text="Add k planes", command=addkplanes)
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Add k planes.')

    # - Extrusion suivant une courbe -
    B = TTK.Button(Frame, text="Curve(s)", command=setCurve,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=2, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set driving curve(s) for extrusion.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White')
    B.grid(row=2, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Driving curve(s) for extrusion.')

    B = TTK.Button(Frame, text="LineDrive", command=lineDrive)
    B.grid(row=3, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extrude following linearly curve(s).')

    B = TTK.Button(Frame, text="OrthoDrive", command=orthoDrive)
    B.grid(row=3, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extrude following orthogonally curve(s).')

    # - Revolve -
    B = TTK.Button(Frame, text="Axis", command=setAxis,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=4, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set axis for revolve.')
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White')
    B.grid(row=4, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Axis for revolve or\nAxis + rmod curve for revolve.')

    B = TTK.Entry(Frame, textvariable=VARS[5], background='White', width=5)
    B.grid(row=5, column=0,  columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Angle for revolve.')
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White', width=5)
    B.grid(row=5, column=1,  columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of points in revolve.')

    B = TTK.Button(Frame, text="Revolve", command=revolve)
    B.grid(row=5, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Revolve following axis.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MeshNoteBook'].add(WIDGETS['frame'], text='tkExtrusion')
    except: pass
    CTK.WIDGETS['MeshNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['MeshNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkExtrusionHeight'] = VARS[0].get()
    CTK.PREFS['tkExtrusionNLayers'] = VARS[1].get()
    CTK.PREFS['tkExtrusionSmooth'] = VARS[2].get()
    CTK.PREFS['tkExtrusionRevAngle'] = VARS[5].get()
    CTK.PREFS['tkExtrusionNpts'] = VARS[6].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('1.e-1')
    VARS[1].set('10')
    VARS[2].set('50')
    VARS[5].set('360.')
    VARS[6].set('30')
    CTK.PREFS['tkExtrusionHeight'] = VARS[0].get()
    CTK.PREFS['tkExtrusionNLayers'] = VARS[1].get()
    CTK.PREFS['tkExtrusionSmooth'] = VARS[2].get()
    CTK.PREFS['tkExtrusionRevAngle'] = VARS[5].get()
    CTK.PREFS['tkExtrusionNpts'] = VARS[6].get()
    CTK.savePrefFile()

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
    (win, menu, file, tools) = CTK.minimal('tkExtrusion '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
