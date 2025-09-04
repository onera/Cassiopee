# - tkSurfaceWalk -
"""Generate meshes by walking on surface."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def walkIn():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    # Constraints
    name = VARS[0].get()
    names = name.split(';')
    constraints = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): constraints.append(z)

    # surfaces
    name = VARS[1].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: surfaces.append(z)

    # - Hauteur de chaque maille -
    dhloc = CTK.varsFromWidget(VARS[2].get(), type=1); dhloc = dhloc[0]
    N = CTK.varsFromWidget(VARS[3].get(), type=2); N = N[0]
    dh = G.cart((0.,0.,0.),(dhloc,1.,1.),(N+1,1,1))

    # - nb d'iterations de lissage -
    nit = CTK.varsFromWidget(VARS[4].get(), type=2); nit = nit[0]
    # - contour -
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    contours = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        contour = C.convertBAR2Struct(z)
        contours.append(contour)
    # surfaceWalk
    zlist = []
    fail = False; errors = []
    for c in contours:
        try:
            z = G.surfaceWalk(surfaces, c, dh, constraints=constraints,
                              niter=nit, check=1)
            zlist.append(z)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    # Ajout dans la base SURFACES
    CTK.t = C.addBase2PyTree(CTK.t, 'SURFACES')
    bases = Internal.getNodesFromName1(CTK.t, 'SURFACES')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in zlist: CTK.add(CTK.t, nob, -1, i)
    #C._fillMissingVariables(CTK.t)
    if not fail: CTK.TXT.insert('START', 'Surface walk done.\n')
    else:
        Panels.displayErrors(errors, header='Error: surfaceWalk')
        CTK.TXT.insert('START', 'Surface walk fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def walkOut():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    # Constraints
    name = VARS[0].get()
    names = name.split(';')
    constraints = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if (bases != []):
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): constraints.append(z)

    # surfaces
    name = VARS[1].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if (bases != []):
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): surfaces.append(z)

    # - Hauteur de chaque maille -
    dhloc = CTK.varsFromWidget(VARS[2].get(), type=1); dhloc = dhloc[0]
    N = CTK.varsFromWidget(VARS[3].get(), type=2); N = N[0]
    dh = G.cart((0.,0.,0.),(dhloc,1.,1.),(N+1,1,1))

    # - nb d'iterations de lissage -
    nit = CTK.varsFromWidget(VARS[4].get(), type=2); nit = nit[0]
    # - contour -
    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    contours = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        contour = C.convertBAR2Struct(z)
        contour = T.reorder(contour,(-1,2,3))
        contours.append(contour)
    # surfaceWalk
    zlist = []
    fail = False; errors = []
    for c in contours:
        try:
            z = G.surfaceWalk(surfaces, c, dh, constraints=constraints,
                              niter=nit, check=1)
            zlist.append(z)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    # Ajout dans la base SURFACES
    CTK.t = C.addBase2PyTree(CTK.t, 'SURFACES')
    bases = Internal.getNodesFromName1(CTK.t, 'SURFACES')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in zlist: CTK.add(CTK.t, nob, -1, i)
    #C._fillMissingVariables(CTK.t)
    if not fail: CTK.TXT.insert('START', 'Surface walk done.\n')
    else:
        Panels.displayErrors(errors, header='Error: surfaceWalk')
        CTK.TXT.insert('START', 'Surface walk fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def setConstraintContour():
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
    VARS[0].set(selected)

#==============================================================================
def setSurface():
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
    VARS[1].set(selected)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkSurfaceWalk  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Generate meshes by orthogonal\nwalk on surfaces.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkSurfaceWalk')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Constraint contour -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- Projection surface -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -2- Hauteur de chaque maille -
    V = TK.StringVar(win); V.set('1.e-1'); VARS.append(V)
    if 'tkSurfaceWalkHeight' in CTK.PREFS:
        V.set(CTK.PREFS['tkSurfaceWalkHeight'])
    # -3- Nombre de layers a ajouter
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    if 'tkSurfaceWalkNLayers' in CTK.PREFS:
        V.set(CTK.PREFS['tkSurfaceWalkNLayers'])
    # -4- Nombre d'iterations de lissage
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkSurfaceWalkSmooth' in CTK.PREFS:
        V.set(CTK.PREFS['tkSurfaceWalkSmooth'])

    # - Surface -
    B = TTK.Button(Frame, text="Proj. surf", command=setSurface,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set projection surfaces.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White')
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Projection surfaces.')

    # - Contours -
    B = TTK.Button(Frame, text="Constraints", command=setConstraintContour,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=1, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set constraint contours.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Constraint contours.')

    # - Walk
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Height of each layer.')

    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=5)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of layers.')

    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=5)
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of smoothing iterations.')

    B = TTK.Button(Frame, text="Walk in", command=walkIn)
    B.grid(row=3, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Walk in.')

    B = TTK.Button(Frame, text="Walk out", command=walkOut)
    B.grid(row=3, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Walk out.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkSurfaceWalk')
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
def updateApp():
    return

#==============================================================================
def saveApp():
    CTK.PREFS['tkSurfaceWalkHeight'] = VARS[2].get()
    CTK.PREFS['tkSurfaceWalkNLayers'] = VARS[3].get()
    CTK.PREFS['tkSurfaceWalkSmooth'] = VARS[4].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[2].set('1.e-1')
    VARS[3].set('1')
    VARS[4].set('0')
    CTK.PREFS['tkSurfaceWalkHeight'] = VARS[0].get()
    CTK.PREFS['tkSurfaceWalkNLayers'] = VARS[1].get()
    CTK.PREFS['tkSurfaceWalkSmooth'] = VARS[1].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if (__name__ == "__main__"):
    import sys
    if (len(sys.argv) == 2):
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkSurfaceWalk '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
