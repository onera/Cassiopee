# - mapping/remeshing -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Generator.PyTree as G
import Geom.PyTree as D
import Transform.PyTree as T
import Post.PyTree as P
import KCore.Vector as Vector

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def mapCurvature():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    dir = VARS[1].get()
    if dir == 'i-indices': dir = 1
    elif dir == 'j-indices': dir = 2
    elif dir == 'k-indices': dir = 3
    else: dir = 0

    power = VARS[4].get()
    try: power = float(power)
    except: power = 1.

    power2 = VARS[5].get()
    try: power2 = float(power2)
    except: power2 = 0.5

    CTK.saveTree()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dims = Internal.getZoneDim(z)
        if dims[3] == 'BAR':
            z = C.convertBAR2Struct(z)
            dims = Internal.getZoneDim(z)
        if dir == 1:
            zp = G.mapCurvature(z, int(power*dims[1]), power2, 1)
        elif dir == 2:
            zp = G.mapCurvature(z, int(power*dims[2]), power2, 2)
        elif dir == 3:
            zp = G.mapCurvature(z, int(power*dims[3]), power2, 3)
        else:
            zp = G.mapCurvature(z, int(power*dims[1]), power2, 1)
            if dims[2] > 1:
                zp = G.mapCurvature(zp, int(power*dims[2]), power2, 2)
            if dims[3] > 1:
                zp = G.mapCurvature(zp, int(power*dims[3]), power2, 3)
        CTK.replace(CTK.t, nob, noz, zp)
    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'MapCurvature done.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def uniformizeMesh(z, N, dir):
    N = int(N)
    distrib = G.cart( (0,0,0), (1./(N-1),1,1), (N,1,1) )
    zp = G.map(z, distrib, dir)
    return zp

#==============================================================================
def uniformize():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    dir = VARS[1].get()
    if dir == 'i-indices': dir = 1
    elif dir == 'j-indices': dir = 2
    elif dir == 'k-indices': dir = 3
    else: dir = 0

    power = VARS[3].get()
    try: power = float(power)
    except: power = 1.

    CTK.saveTree()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dims = Internal.getZoneDim(z)
        if dims[3] == 'BAR':
            z = C.convertBAR2Struct(z)
            dims = Internal.getZoneDim(z)
        if dir == 1:
            zp = uniformizeMesh(z, power*dims[1], 1)
        elif dir == 2:
            zp = uniformizeMesh(z, power*dims[2], 2)
        elif dir == 3:
            zp = uniformizeMesh(z, power*dims[3], 3)
        else:
            zp = uniformizeMesh(z, power*dims[1], 1)
            if dims[2] > 1:
                zp = uniformizeMesh(zp, power*dims[2], 2)
            if dims[3] > 1:
                zp = uniformizeMesh(zp, power*dims[3], 3)
        CTK.replace(CTK.t, nob, noz, zp)
    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Uniformize done.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def refine():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    dir = VARS[1].get()
    if dir == 'i-indices': dir = 1
    elif dir == 'j-indices': dir = 2
    elif dir == 'k-indices': dir = 3
    else: dir = 0

    power = VARS[2].get()
    try: power = float(power)
    except: power = 2.

    CTK.saveTree()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dim = Internal.getZoneDim(z)
        if dim[3] == 'TRI': # raffinement TRI
            iter = int(round(power / 2.,0))
            for i in xrange(iter):
                 #z = C.initVars(z, 'centers:__indic__', 1)
                 #z = P.refine(z, 'centers:__indic__')
                 P._refine(z, w=1./16.) # butterfly
                 #z = C.rmVars(z, 'centers:__indic__')
            CTK.replace(CTK.t, nob, noz, z)

        elif (dim[3] == 'BAR' or dim[0] == 'Structured'):
            if dim[3] == 'BAR': z = C.convertBAR2Struct(z)
            z = G.refine(z, power, dir)
            CTK.replace(CTK.t, nob, noz, z)

    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Refine done.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def enforceMesh(z, dir, N, width, ind, h):
    i1 = ind[2]; j1 = ind[3]; k1 = ind[4]
    if (dir == 1):
        zt = T.subzone(z, (1,j1,k1), (N,j1,k1))
    elif (dir == 2):
        zt = T.subzone(z, (i1,1,k1), (i1,N,k1))
        zt = T.reorder(zt, (2,1,3))
    else:
        zt = T.subzone(z, (i1,j1,1), (i1,j1,N))
        zt = T.reorder(zt, (3,2,1))

    l = D.getLength(zt)
    zt = D.getCurvilinearAbscissa(zt)
    distrib = C.cpVars(zt, 's', zt, 'CoordinateX')
    distrib = C.initVars(distrib, 'CoordinateY', 0.)
    distrib = C.initVars(distrib, 'CoordinateZ', 0.)
    distrib = C.rmVars(distrib, 's')

    Nr = int(width*N);
    #print h, l, Nr, i1, j1, k1

    if dir == 1:
        val = C.getValue(zt, 's', i1-1)
        indl = (i1,j1,k1); inp1 = (i1+1,j1,k1); indm1 = (i1-1,j1,k1)
    elif dir == 2:
        val = C.getValue(zt, 's', j1-1)
        indl = (i1,j1,k1); inp1 = (i1,j1+1,k1); indm1 = (i1,j1-1,k1)
    else:
        val = C.getValue(zt, 's', k1-1)
        indl = (i1,j1,k1); inp1 = (i1,j1,k1+1); indm1 = (i1,j1,k1-1)

    Xc = CPlot.getActivePoint()
    valf = val
    Pind = C.getValue(z, 'GridCoordinates', indl)
    if (ind < N-1): # cherche avec indp1
        Pindp1 = C.getValue(z, 'GridCoordinates', indp1)
        v1 = Vector.sub(Pindp1, Pind)
        v2 = Vector.sub(Xc, Pind)
        if Vector.dot(v1,v2) >= 0:
            val2 = C.getValue(zt, 's', i1)
            alpha = Vector.norm(v2)/Vector.norm(v1)
            valf = val+alpha*(val2-val)
    if (ind > 0 and val == valf): # cherche avec indm1
        Pindm1 = C.getValue(z, 'GridCoordinates', indm1)
        v1 = Vector.sub(Pindm1, Pind)
        v2 = Vector.sub(Xc, Pind)
        if Vector.dot(v1,v2) >= 0:
            val2 = C.getValue(zt, 's', i1-2)
            alpha = Vector.norm(v2)/Vector.norm(v1)
            valf = val+alpha*(val2-val)

    if (h < 0): distrib = G.enforcePoint(distrib, valf)
    elif (dir == 1):
        if (i1 == 1):
            distrib = G.enforcePlusX(distrib, h/l, Nr, Nr)
        elif (i1 == N):
            distrib = G.enforceMoinsX(distrib, h/l, Nr, Nr)
        else:
            distrib = G.enforceX(distrib, valf, h/l, Nr, Nr)
    elif (dir == 2):
        if (j1 == 1):
            distrib = G.enforcePlusX(distrib, h/l, Nr, Nr)
        elif (j1 == N):
            distrib = G.enforceMoinsX(distrib, h/l, Nr, Nr)
        else:
            distrib = G.enforceX(distrib, valf, h/l, Nr, Nr)
    elif (dir == 3):
        if (k1 == 1):
            distrib = G.enforcePlusX(distrib, h/l, Nr, Nr)
        elif (k1 == N):
            distrib = G.enforceMoinsX(distrib, h/l, Nr, Nr)
        else:
            distrib = G.enforceX(distrib, valf, h/l, Nr, Nr)
    try:
        z1 = G.map(z, distrib, dir); return z1
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: enforce')
        return None

#==============================================================================
def enforce():
    if (CTK.t == []): return
    if (CTK.__MAINTREE__ <= 0):
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    width = WIDGETS['slider'].get() / 100.

    h = VARS[0].get()
    try: h = float(h)
    except: h = 1.e-6
    xdir = VARS[1].get()

    ind = CPlot.getActivePointIndex()

    fail = False
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dims = Internal.getZoneDim(z)
        if (dims[0] == 'Structured'):
            if (xdir == 'i-indices'):
                dir = 1; N = dims[1]
                zp = enforceMesh(z, dir, N, width, ind, h)
            elif (xdir == 'j-indices'):
                dir = 2; N = dims[2]
                zp = enforceMesh(z, dir, N, width, ind, h)
            elif (xdir == 'k-indices'):
                dir = 3; N = dims[3]
                zp = enforceMesh(z, dir, N, width, ind, h)
            else: # les 3 dirs
                zp = enforceMesh(z, 1, dims[1], width, ind, h)
                if (zp is not None and dims[2] > 1):
                    zp = enforceMesh(zp, 2, dims[2], width, ind, h)
                if (zp is not None and dims[3] > 1):
                    zp = enforceMesh(zp, 3, dims[3], width, ind, h)
        if zp is not None: CTK.replace(CTK.t, nob, noz, zp)
        else: fail = True

    if fail:
        CTK.TXT.insert('START', 'Stretch failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    else: CTK.TXT.insert('START', 'Stretch done.\n')
    #C._fillMissingVariables(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def setMeshWidth(event=None):
    val = WIDGETS['slider'].get()
    VARS[6].set('Width of mesh concerned with remeshing (%.2f %% of points).'%(val / 100.))

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkStretch', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Stretch meshes locally.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=4)
    Frame.columnconfigure(2, weight=4)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkStretch')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- enforced height -
    V = TK.StringVar(win); V.set('1.e-6'); VARS.append(V)
    if 'tkStretchHeight' in CTK.PREFS:
        V.set(CTK.PREFS['tkStretchHeight'])
    # -1- direction pour remap
    V = TK.StringVar(win); V.set('i-j-k indices'); VARS.append(V)
    # -2- refinement power for refine
    V = TK.StringVar(win); V.set('2.'); VARS.append(V)
    # -3- refinement factor for uniformize
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -4- refinement factor for mapCurvature
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -5- power for mapCurvature
    V = TK.StringVar(win); V.set('0.5'); VARS.append(V)
    if 'tkStretchCurvPower' in CTK.PREFS:
        V.set(CTK.PREFS['tkStretchCurvPower'])
    # -6- Width of mesh info bulle
    V = TK.StringVar(win); V.set('Width of mesh concerned with remeshing (% of points).'); VARS.append(V)

    # - Slider -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  command=setMeshWidth, borderwidth=1, value=50)
    WIDGETS['slider'] = B
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[6])
    
    # - Direction -
    B = TTK.OptionMenu(Frame, VARS[1], 'i-indices', 'j-indices', 'k-indices',
                       'i-j-k indices')
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Index direction for stretch or refine.')

    # - enforce -
    B = TTK.Button(Frame, text="Enforce", command=enforce)
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Enforce given height in mesh following dir.')

    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=3)
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Enforced height in mesh.')

    # - refine -
    B = TTK.Button(Frame, text="Refine", command=refine)
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refine a given mesh keeping original distribution.')
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=3)
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refinement factor.')

    # - uniformize -
    B = TTK.Button(Frame, text="Uniformize", command=uniformize)
    B.grid(row=3, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Uniformize the mesh distribution.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=3)
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refinement factor.')

    # - mapCurvature -
    B = TTK.Button(Frame, text="mapCurvature", command=mapCurvature)
    B.grid(row=4, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Remesh following curvature.')
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=3)
    B.grid(row=4, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refinement factor.')
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White', width=3)
    B.grid(row=4, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Proportionality to curvature.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW)

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkStretchHeight'] = VARS[0].get()
    CTK.PREFS['tkStretchCurvPower'] = VARS[5].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('1.e-1')
    VARS[5].set('0.5')
    CTK.PREFS['tkStretchHeight'] = VARS[0].get()
    CTK.PREFS['tkStretchCurvPower'] = VARS[5].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)
    
#==============================================================================
if (__name__ == "__main__"):
    import sys
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkStretch '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
