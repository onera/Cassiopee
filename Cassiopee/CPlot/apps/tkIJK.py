# - tkIJK -
"""Show only sub-indices of structured grids."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Transform.PyTree as T
import Post.PyTree as P
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Generator.PyTree as G
import Converter.Internal as Internal

# VIEWDATA
XYZVIEWDATA = None
IJKVIEWDATA = None

# local widgets list
WIDGETS = {}; VARS = []

def setPlusI(event=None):
    set__(VARS[0], VARS[1], VARS[6], +1)

def setMoinsI(event=None):
    set__(VARS[0], VARS[1], VARS[6], -1)

def setFullI(event=None):
    setFull__(VARS[0], VARS[1])

def setPlusJ(event=None):
    set__(VARS[2], VARS[3], VARS[7], +1)

def setMoinsJ(event=None):
    set__(VARS[2], VARS[3], VARS[7], -1)

def setFullJ(event=None):
    setFull__(VARS[2], VARS[3])

def setPlusK(event=None):
    set__(VARS[4], VARS[5], VARS[8], +1)

def setMoinsK(event=None):
    set__(VARS[4], VARS[5], VARS[8], -1)

def setFullK(event=None):
    setFull__(VARS[4], VARS[5])

def set__(V0, V1, Voff, off):
    v = CTK.varsFromWidget(Voff.get(), type=2)
    if len(v) > 0: off = v[0]*off
    else: return
    v = CTK.varsFromWidget(V0.get(), type=2)
    if len(v) > 0: imin = v[0]
    else: return
    v = CTK.varsFromWidget(V1.get(), type=2)
    if len(v) > 0: imax = v[0]
    else: return
    m = min(imin, imax)+off
    m = max(m,1)
    V0.set(str(m))
    V1.set(str(m))
    show()

def setFull__(V0, V1):
    V0.set(str(1))
    V1.set(str(-1))
    show()

#==============================================================================
# show
#==============================================================================
def show(event=None):
    if CTK.t == []: return
    # Get indices
    v = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(v) > 0: imin = v[0]
    v = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(v) > 0: imax = v[0]
    v = CTK.varsFromWidget(VARS[2].get(), type=2)
    if len(v) > 0: jmin = v[0]
    v = CTK.varsFromWidget(VARS[3].get(), type=2)
    if len(v) > 0: jmax = v[0]
    v = CTK.varsFromWidget(VARS[4].get(), type=2)
    if len(v) > 0: kmin = v[0]
    v = CTK.varsFromWidget(VARS[5].get(), type=2)
    if len(v) > 0: kmax = v[0]

    if CTK.__MAINTREE__ == CTK.MAIN:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()
    active = []
    tp = Internal.appendBaseName2ZoneName(CTK.t, updateRef=False,
                                          separator=Internal.SEP1)
    for z in CTK.__MAINACTIVEZONES__: active.append(tp[2][CTK.Nb[z]+1][2][CTK.Nz[z]])

    temp = C.newPyTree(['Base']); temp[2][1][2] += active

    if CTK.__MAINTREE__ == CTK.SLICE:
        exts = Internal.getNodeFromName1(CTK.dt, 'Edges')
        if exts is None: exts = []
        else: exts = exts[2]

    # Subzone les zones actives
    CTK.dt = C.newPyTree(['Base','Edges'])
    for z in Internal.getZones(temp):
        try:
            zp = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax)); zp[0] = z[0]
            CTK.dt[2][1][2].append(zp)
        except: pass
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp(CTK.dt)

    # Ajoute les edges pour les grilles structurees
    if CTK.__MAINTREE__ == CTK.MAIN:
        exts = []
        for z in active:
            ztype = Internal.getZoneType(z)
            if ztype == 1:
                zp = P.exteriorFacesStructured(z)
                exts += zp
            else:
                #zp = P.exteriorFaces(z)
                #zp = P.sharpEdges(zp)
                zp = []
                exts += zp
    CTK.dt[2][2][2] += exts

    # deactivate les zones exts
    lenZ = len(CTK.dt[2][1][2]); lenExts = len(CTK.dt[2][2][2])
    activeExt = [(i,1) for i in range(lenZ+lenExts)]
    for i in range(lenZ): activeExt[i] = (i,1)
    for i in range(lenExts): activeExt[i+lenZ] = (i+lenZ,0)

    CTK.display(CTK.dt, mainTree=CTK.SLICE)
    CPlot.setActiveZones(activeExt)
    CPlot.setState(edgifyDeactivatedZones=1)

#==============================================================================
# extract : ajoute l'entite dans CTK.t/EXTRACT
#==============================================================================
def extract(event=None):
    if CTK.t == []: return
    # Get indices
    v = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(v) > 0: imin = v[0]
    v = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(v) > 0: imax = v[0]
    v = CTK.varsFromWidget(VARS[2].get(), type=2)
    if len(v) > 0: jmin = v[0]
    v = CTK.varsFromWidget(VARS[3].get(), type=2)
    if len(v) > 0: jmax = v[0]
    v = CTK.varsFromWidget(VARS[4].get(), type=2)
    if len(v) > 0: kmin = v[0]
    v = CTK.varsFromWidget(VARS[5].get(), type=2)
    if len(v) > 0: kmax = v[0]

    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()
    active = []
    tp = Internal.appendBaseName2ZoneName(CTK.t, updateRef=False,
                                          separator=Internal.SEP1)
    for z in CTK.__MAINACTIVEZONES__: active.append(tp[2][CTK.Nb[z]+1][2][CTK.Nz[z]])

    temp = C.newPyTree(['Base']); temp[2][1][2] += active

    C._addBase2PyTree(CTK.t, 'EXTRACT')
    base = Internal.getNodeFromName1(CTK.t, 'EXTRACT')

    # Subzone les zones actives
    for z in Internal.getZones(temp):
        try:
            zp = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax))
            zp[0] = C.getZoneName(z[0])
            base[2].append(zp)
        except: pass

    CTK.display(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.TXT.insert('START', 'Zone extracted in base EXTRACT.\n')

# Met i,j,k dans Coordinates des zones
def setIJK2Coordinates(z):
    dim = Internal.getZoneDim(z)
    zc = G.cart((0,0,0), (1,1,1), (dim[1], dim[2], dim[3]))
    FC = Internal.getNodeFromName1(z, Internal.__FlowSolutionCenters__)
    FN = Internal.getNodeFromName1(z, Internal.__FlowSolutionNodes__)
    if FN is not None: zc[2].append(FN)
    if FC is not None: zc[2].append(FC)
    return zc

#==============================================================================
# viewIJK : affiche la selection en coordonnees IJK
#==============================================================================
def viewIJK(event=None):
    global XYZVIEWDATA
    if CTK.t == []: return
    if CTK.__MAINTREE__ != 1: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    sel = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        if Internal.getZoneType(z) == 1:
            zp = setIJK2Coordinates(z)
            sel.append(zp)
    CTK.dt = C.newPyTree(['Base'])
    CTK.dt[2][1][2] += sel
    XYZVIEWDATA = CPlot.getState('posCam') + CPlot.getState('posEye') + CPlot.getState('dirCam')
    if IJKVIEWDATA is not None:
        posCam = IJKVIEWDATA[0:3]
        posEye = IJKVIEWDATA[3:6]
        dirCam = IJKVIEWDATA[6:9]
        CTK.display(CTK.dt, posCam=posCam, posEye=posEye, dirCam=dirCam, mainTree=CTK.IJK)
    else:
        CTK.display(CTK.dt, mainTree=CTK.IJK)
        CPlot.fitView()
    CTK.TXT.insert('START', 'Viewing IJK of selection.\n')

def backFromIJK(event=None):
    global IJKVIEWDATA
    if CTK.t == []: return
    if CTK.__MAINTREE__ == CTK.IJK:
        IJKVIEWDATA = CPlot.getState('posCam') + CPlot.getState('posEye') + CPlot.getState('dirCam')

    if CTK.__MAINTREE__ != 1 and XYZVIEWDATA is not None:
        posCam = XYZVIEWDATA[0:3]
        posEye = XYZVIEWDATA[3:6]
        dirCam = XYZVIEWDATA[6:9]
        CTK.display(CTK.t, posCam=posCam, posEye=posEye, dirCam=dirCam, mainTree=CTK.MAIN)
        CTK.TXT.insert('START', 'Viewing full mesh.\n')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkIJK  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Display mesh informations.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkIJK')
    WIDGETS['frameMenu'] = FrameMenu

    #- VARS -
    # -0- Imin -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -1- Imax -
    V = TK.StringVar(win); V.set('-1'); VARS.append(V)
    # -2- Jmin -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -3- Jmax -
    V = TK.StringVar(win); V.set('-1'); VARS.append(V)
    # -4- Kmin -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -5- Kmax -
    V = TK.StringVar(win); V.set('-1'); VARS.append(V)
    # -6- stepI -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -7- stepJ -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -8- stepK -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)

    # - Subzone indices -
    F = TTK.LabelFrame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE,
                       text='I', font=CTK.FRAMEFONT, takefocus=1)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=1)
    F.columnconfigure(2, weight=0)
    F.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text='+', width=5, command=setPlusI)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(F, text='-', width=5, command=setMoinsI)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text='X', width=1, command=setFullI)
    BB = CTK.infoBulle(parent=B, text='Reset to full zone.')
    B.grid(row=0, column=2, sticky=TK.EW)

    B = TTK.Entry(F, textvariable=VARS[0], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='imin.')
    B.bind('<Return>', show)
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[1], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='imax.')
    B.bind('<Return>', show)
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[6], background='White', width=1)
    BB = CTK.infoBulle(parent=B, text='i step.')
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)

    F = TTK.LabelFrame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE,
                       text='J', font=CTK.FRAMEFONT, takefocus=1)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=1)
    F.columnconfigure(2, weight=0)
    F.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(F, text='+', width=5, command=setPlusJ)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(F, text='-', width=5, command=setMoinsJ)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text='X', width=1, command=setFullJ)
    BB = CTK.infoBulle(parent=B, text='Reset to full zone.')
    B.grid(row=0, column=2, sticky=TK.EW)

    B = TTK.Entry(F, textvariable=VARS[2], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='jmin.')
    B.bind('<Return>', show)
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[3], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='jmax.')
    B.bind('<Return>', show)
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[7], background='White', width=1)
    BB = CTK.infoBulle(parent=B, text='j step.')
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)

    F = TTK.LabelFrame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE,
                       text='K', font=CTK.FRAMEFONT, takefocus=1)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=1)
    F.columnconfigure(2, weight=0)
    F.grid(row=0, column=2, sticky=TK.EW)
    B = TTK.Button(F, text='+', width=5, command=setPlusK)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(F, text='-', width=5, command=setMoinsK)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text='X', width=1, command=setFullK)
    BB = CTK.infoBulle(parent=B, text='Reset to full zone.')
    B.grid(row=0, column=2, sticky=TK.EW)

    B = TTK.Entry(F, textvariable=VARS[4], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='kmin.')
    B.bind('<Return>', show)
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[5], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='kmax.')
    B.bind('<Return>', show)
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[8], background='White', width=1)
    BB = CTK.infoBulle(parent=B, text='k step.')
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)

    B = TTK.Button(Frame, text="Extract", command=extract)
    B.grid(row=2, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract subzones to main tree.')

    B = TTK.Button(Frame, text="View IJK", command=viewIJK)
    B.grid(row=3, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='View selection as IJK grid.')

    B = TTK.Button(Frame, text="Back", command=backFromIJK)
    B.grid(row=3, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Back to main view.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['VisuNoteBook'].add(WIDGETS['frame'], text='tkIJK')
    except: pass
    CTK.WIDGETS['VisuNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['VisuNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkIJJK '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
