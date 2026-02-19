# - tkCollarMesh -
"""Volume collar mesh generation."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Generator.PyTree as G
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setConstraintContour1():
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
    VARS[8].set(selected)

#==============================================================================
def setConstraintContour2():
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
    VARS[9].set(selected)

#==============================================================================
def setSurface1():
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
def setSurface2():
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
def unionCollarMesh():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    # - surfaces1 -
    name = VARS[0].get()
    names = name.split(';')
    surface1 = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): surface1.append(z)

    # - surfaces2 -
    name = VARS[1].get()
    names = name.split(';')
    surface2 = []
    for v in names:
        v = v.lstrip(); v= v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: surface2.append(z)

    # - constraints1 -
    constraints1 = []
    name = VARS[8].get()
    names = name.split(';')
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: constraints1.append(z)

    # - constraints2 -
    constraints2 = []
    name = VARS[9].get()
    names = name.split(';')
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: constraints2.append(z)

    # - Hauteur de chaque maille en j -
    dhj = CTK.varsFromWidget(VARS[2].get(), type=1); dhj = dhj[0]
    Nj = CTK.varsFromWidget(VARS[3].get(), type=2); Nj = Nj[0]
    distribj = G.cart((0.,0.,0.),(dhj,1.,1.),(Nj+1,1,1))

    # - Hauteur de chaque maille en k -
    dhk = CTK.varsFromWidget(VARS[5].get(), type=1); dhk = dhk[0]
    Nk = CTK.varsFromWidget(VARS[6].get(), type=2); Nk = Nk[0]
    distribk = G.cart((0.,0.,0.),(dhk,1.,1.),(Nk+1,1,1))

    # - Nb d'iterations de lissage -
    niterj = CTK.varsFromWidget(VARS[4].get(), type=2); niterj = niterj[0]
    niterk = CTK.varsFromWidget(VARS[7].get(), type=2); niterk = niterk[0]

    # - Contour -
    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    contours = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        contour = C.convertBAR2Struct(z)
        contours.append(contour)
    zlist = []
    if contours != []:
        for c in contours:
            zones = G.collarMesh(surface1, surface2, distribj, distribk,
                                 niterj=niterj, niterk=niterk, alphaRef=180.,
                                 type='union', contour=c,
                                 constraints1=constraints1,
                                 constraints2=constraints2)
            zlist += zones
    else:
        zones = G.collarMesh(surface1, surface2, distribj, distribk,
                             niterj=niterj, niterk=niterk,
                             type='union', constraints1=constraints1,
                             constraints2=constraints2)
        zlist += zones

    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'COLLAR')
    bases = Internal.getNodesFromName1(CTK.t, 'COLLAR')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in zlist: CTK.add(CTK.t, nob, -1, i)
    CTK.t = C.fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Collar grids created.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def differenceCollarMesh():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    # surfaces1
    name = VARS[0].get()
    names = name.split(';')
    surface1 = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: surface1.append(z)

    # surfaces2
    name = VARS[1].get()
    names = name.split(';')
    surface2 = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: surface2.append(z)
    # - constraints1 -
    constraints1 = []
    name = VARS[8].get()
    names = name.split(';')
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: constraints1.append(z)

    # - constraints2 -
    constraints2 = []
    name = VARS[9].get()
    names = name.split(';')
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: constraints2.append(z)

    # - Hauteur de chaque maille en j -
    dhj = CTK.varsFromWidget(VARS[2].get(), type=1); dhj = dhj[0]
    Nj = CTK.varsFromWidget(VARS[3].get(), type=2); Nj = Nj[0]
    distribj = G.cart((0.,0.,0.),(dhj,1.,1.),(Nj+1,1,1))

    # - Hauteur de chaque maille en k -
    dhk = CTK.varsFromWidget(VARS[5].get(), type=1); dhk = dhk[0]
    Nk = CTK.varsFromWidget(VARS[6].get(), type=2); Nk = Nk[0]
    distribk = G.cart((0.,0.,0.),(dhk,1.,1.),(Nk+1,1,1))

    # - Nb d'iterations de lissage -
    niterj = CTK.varsFromWidget(VARS[4].get(), type=2); niterj = niterj[0]
    niterk = CTK.varsFromWidget(VARS[7].get(), type=2); niterk = niterk[0]

    # - Contour -
    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    contours = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        contour = C.convertBAR2Struct(z)
        contours.append(contour)
    zlist = []
    if contours != []:
        for c in contours:
            zones = G.collarMesh(surface1, surface2, distribj, distribk,
                                 niterj=niterj, niterk=niterk, alphaRef=180.,
                                 type='difference', contour=c,
                                 constraints1=constraints1,
                                 constraints2=constraints2)
            zlist += zones
    else:
        zones = G.collarMesh(surface1, surface2, distribj, distribk,
                             niterj=niterj, niterk=niterk,
                             type='difference', constraints1=constraints1,
                             constraints2=constraints2)
        zlist += zones
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'COLLAR')
    bases = Internal.getNodesFromName1(CTK.t, 'COLLAR')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in zones: CTK.add(CTK.t, nob, -1, i)
    CTK.t = C.fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Collar grids created.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCollarMesh  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Generate collar meshes.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=2)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=1)
    Frame.columnconfigure(3, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkCollarMesh')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Surface1-
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- Surface2-
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -2- hauteur de la maille dans la direction j
    V = TK.StringVar(win); V.set('1.e-1'); VARS.append(V)
    if 'tkCollarMeshHj' in CTK.PREFS: V.set(CTK.PREFS['tkCollarMeshHj'])
    # -3- nombre de layers a ajouter en j
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    if 'tkCollarMeshNj' in CTK.PREFS: V.set(CTK.PREFS['tkCollarMeshNj'])
    # -4- nombre d'iteration de lissage dans la direction j
    V = TK.StringVar(win); V.set('50'); VARS.append(V)
    if 'tkCollarMeshSj' in CTK.PREFS: V.set(CTK.PREFS['tkCollarMeshSj'])
    # -5- hauteur de la maille dans la direction k
    V = TK.StringVar(win); V.set('1.e-1'); VARS.append(V)
    if 'tkCollarMeshHk' in CTK.PREFS: V.set(CTK.PREFS['tkCollarMeshHk'])
    # -6- nombre de layers a ajouter en k
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    if 'tkCollarMeshNk' in CTK.PREFS: V.set(CTK.PREFS['tkCollarMeshNk'])
    # -7- nombre d'iteration de lissage dans la direction k
    V = TK.StringVar(win); V.set('50'); VARS.append(V)
    if 'tkCollarMeshSk' in CTK.PREFS: V.set(CTK.PREFS['tkCollarMeshSk'])
    # -8- Constraints for surface1
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -9- Constraints for surface2
    V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Surface1 -
    B = TTK.Button(Frame, text="Surf1", command=setSurface1,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=0, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set collar surfaces 1.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)

    # - Surface2 -
    B = TTK.Button(Frame, text="Surf2", command=setSurface2,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=1, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set collar surfaces 2.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White')
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)

    # - Constraints1 -
    B = TTK.Button(Frame, text="Constr.1", command=setConstraintContour1,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=2, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set constraint contours for Surf1.')
    B = TTK.Entry(Frame, textvariable=VARS[8], background='White')
    B.grid(row=2, column=0,  columnspan=2, sticky=TK.EW)

    # - Constraints2 -
    B = TTK.Button(Frame, text="Constr.2", command=setConstraintContour2,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=3, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set constraint contours for Surf2.')
    B = TTK.Entry(Frame, textvariable=VARS[9], background='White')
    B.grid(row=3, column=0,  columnspan=2, sticky=TK.EW)

    # - Walk in j -
    B = TTK.Label(Frame, text="j-dir", width=2)
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=4, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Height of each j-layer.')

    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=5)
    B.grid(row=4, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of j-layers.')

    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=5)
    B.grid(row=4, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of j-smoothing iterations.')

    # - Walk in k -
    B = TTK.Label(Frame, text="k-dir", width=2)
    B.grid(row=5, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White', width=5)
    B.grid(row=5, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Height of each k-layer.')

    B = TTK.Entry(Frame, textvariable=VARS[6], background='White', width=5)
    B.grid(row=5, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of k-layers.')

    B = TTK.Entry(Frame, textvariable=VARS[7], background='White', width=5)
    B.grid(row=5, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of k-smoothing iterations.')

    # - Union -
    B = TTK.Button(Frame, text="Union", command=unionCollarMesh)
    B.grid(row=6, column=0, columnspan=2,sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Collar mesh for union assembly.')

    # - Difference -
    B = TTK.Button(Frame, text="Difference", command=differenceCollarMesh)
    B.grid(row=6, column=2, columnspan=2,sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Collar mesh for difference assembly.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MeshNoteBook'].add(WIDGETS['frame'], text='tkCollarMesh')
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
    CTK.PREFS['tkCollarMeshHj'] = VARS[2].get()
    CTK.PREFS['tkCollarMeshNj'] = VARS[3].get()
    CTK.PREFS['tkCollarMeshSj'] = VARS[4].get()
    CTK.PREFS['tkCollarMeshHk'] = VARS[5].get()
    CTK.PREFS['tkCollarMeshNk'] = VARS[6].get()
    CTK.PREFS['tkCollarMeshSk'] = VARS[7].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[2].set('1.e-1')
    VARS[3].set('10')
    VARS[4].set('50')
    VARS[5].set('1.e-1')
    VARS[6].set('10')
    VARS[7].set('50')
    CTK.PREFS['tkCollarMeshHj'] = VARS[2].get()
    CTK.PREFS['tkCollarMeshNj'] = VARS[3].get()
    CTK.PREFS['tkCollarMeshSj'] = VARS[4].get()
    CTK.PREFS['tkCollarMeshHk'] = VARS[5].get()
    CTK.PREFS['tkCollarMeshNk'] = VARS[6].get()
    CTK.PREFS['tkCollarMeshSk'] = VARS[7].get()
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
    (win, menu, file, tools) = CTK.minimal('tkCollarMesh '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
