# - tkCartWrap -
"""Surface remesh with cartesian wrapper."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setConstraint():
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

#==============================================================================
def mapSurf(z, density, smoothIter, eltType, constraints, strength):
    import Generator as G
    import Converter
    import Connector as X
    import Post as P
    import Transform as T

    a = C.getFields(Internal.__GridCoordinates__, z)[0]

    # Grille cartesienne (reguliere)
    BB = G.bbox([a])
    xmin = BB[0]; ymin = BB[1]; zmin = BB[2]
    xmax = BB[3]; ymax = BB[4]; zmax = BB[5]
    ni = density*(xmax-xmin); nj = density*(ymax-ymin);
    nk = density*(zmax-zmin)
    if ni < 2: ni = 2
    if nj < 2: nj = 2
    if nk < 2: nk = 2
    hi = (xmax-xmin)/(ni-1); hj = (ymax-ymin)/(nj-1); hk = (zmax-zmin)/(nk-1)
    h = min(hi, hj); h = min(h, hk)
    ni = int((xmax-xmin)/h)+7; nj = int((ymax-ymin)/h)+7
    nk = int((zmax-zmin)/h)+7
    b = G.cart( (xmin-3*h, ymin-3*h, zmin-3*h), (h, h, h), (ni,nj,nk) )
    celln = Converter.array('cellN', ni-1, nj-1, nk-1)
    celln = Converter.initVars(celln, 'cellN', 1.)

    # Masquage
    xrayDim1 = 2000; xrayDim2 = 2000
    cellno = X.blankCells([b], [celln], [a], blankingType=1, delta=h*0.05, 
                          dim=3, XRaydim1=xrayDim1, XRaydim2=xrayDim2)

    # Selection du front
    b1 = P.selectCells2(b, cellno[0], strict=1)
    #Converter.convertArrays2File([b1,a], 'select.plt')
    bext = P.exteriorFaces(b1)
    #Converter.convertArrays2File([bext], 'ext.plt')
    bexts = T.splitConnexity(bext)
    f = bexts[1]; f = Converter.initVars(f, 'cellN', 1)
    f = T.reorder(f, (-1,))
    #Converter.convertArrays2File([f], 'front.plt')

    # Lissage du front
    f = T.smooth(f, eps=0.5, niter=10) # lissage du front quad

    if eltType == 'TRI' or eltType == 'QUAD':
        # Projection du front
        if smoothIter == -1:
            f = T.projectOrtho(f, [a])
        elif smoothIter == 0:
            f = T.projectOrtho(f, [a])
            if constraints != []: 
                f = Converter.initVars(f, 'indic', 0)
                f = G.snapSharpEdges(f, constraints, step=0.3*h, angle=30.)
        else:
            f = T.projectOrtho(f, [a])
            #if constraints != []:
            #    f = G.snapSharpEdges(f, constraints, step=0.3*h, angle=30.)
            for smooth in range(smoothIter):
                f = T.smooth(f, eps=0.5, niter=2, 
                             fixedConstraints=constraints,
                             #projConstraints=constraints,
                             delta=strength)
                f = T.projectOrtho(f, [a])
            if constraints != []:
                f = Converter.initVars(f, 'indic', 0)
                f = G.snapSharpEdges(f, constraints, step=0.3*h, angle=30.)
        proj = f
        if eltType == 'TRI': proj = Converter.converter.convertQuad2Tri(proj)
        return C.convertArrays2ZoneNode('remaped', [proj])

    if eltType == 'STRUCT':
        N = 2
        distrib = G.cart((0,0,0), (1./(N-1),1,1), (N,1,1))
        if smoothIter == 0: f = T.projectOrtho(f, [a])
        for smooth in range(smoothIter):
            f = T.projectOrtho(f, [a])
            f = T.smooth(f, eps=0.5, niter=2, fixedConstraints=constraints,
                         delta=strength)
        r = G.fillWithStruct(f, Vmin=10)

        for smooth in range(smoothIter):
            print('Smoothing iteration %d...'%smooth)
            r = T.smooth(r, eps=0.5, niter=2, fixedConstraints=constraints,
                         delta=strength)
            print('Projection iteration %d...'%smooth)
            r = T.projectOrtho(r, [a])

        ret = []
        for z in r:
            ret.append(C.convertArrays2ZoneNode('remaped', [z]))
        return ret

#==============================================================================
def remap(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    sel = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        sel.append(z)
    if len(sel) != 1:
        CTK.TXT.insert('START', 'Remap can only map one surface.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Elt type
    eltType = VARS[2].get()

    # density
    density = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(density) != 1:
        CTK.TXT.insert('START', 'Density is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    density = density[0]

    # smooth
    smooth = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(smooth) != 1:
        CTK.TXT.insert('START', 'Smooth iter is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    smooth = smooth[0]

    # Constraint strength
    strength = CTK.varsFromWidget(VARS[4].get(), type=1)
    if len(strength) != 1:
        CTK.TXT.insert('START', 'Constraint strength is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    strength = strength[0]

    # Constraint
    name = VARS[3].get()
    names = name.split(';')
    constraints = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]:
                    coord = C.getFields(Internal.__GridCoordinates__, z)[0]
                    constraints.append(coord)
    CTK.saveTree()

    nob = CTK.Nb[nzs[0]]+1
    noz = CTK.Nz[nzs[0]]
    z = CTK.t[2][nob][2][noz]
    zp = mapSurf(z, density, smooth, eltType, constraints, strength)
    C._rmVars(zp, 'nodes:cellN')
    if eltType == 'TRI' or eltType == 'QUAD':
        CTK.replace(CTK.t, nob, noz, zp)
    else:
        CTK.replace(CTK.t, nob, noz, zp[0])
        for i in zp[1:]: CTK.add(CTK.t, nob, -1, i)

    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Surface remapped.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCartWrap  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Remap surfaces with cartesian wrapper.\nCtrl+w to close applet.', temps=0, btype=1)
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
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkCartWrap')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Point density -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkCartWrapDensity' in CTK.PREFS: 
        V.set(CTK.PREFS['tkCartWrapDensity'])
    # -1- Smoother power -
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    if 'tkCartWrapSmooth' in CTK.PREFS: 
        V.set(CTK.PREFS['tkCartWrapSmooth'])
    # -2- Elt type -
    V = TK.StringVar(win); V.set('QUAD'); VARS.append(V)
    if 'tkCartWrapElts' in CTK.PREFS: 
        V.set(CTK.PREFS['tkCartWrapElts'])
    # -3- Constraint
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -4- Constraint strength
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkCartWrapConsStrength' in CTK.PREFS: 
        V.set(CTK.PREFS['tkCartWrapConsStrength'])

    # - Point density -
    B = TTK.Label(Frame, text="Point density")
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=8)
    B.grid(row=0, column=1, sticky=TK.EW)
    B.bind('<Return>', remap)
    BB = CTK.infoBulle(parent=B, text='Point density (npts/length unit).')

    # - Smoother power -
    B = TTK.Label(Frame, text="Smooth iter")
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=8)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of smoother iterations.')

    # - Constraint strength -
    B = TTK.Label(Frame, text="Constraint strength")
    B.grid(row=2, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=8)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Constraint strength.')

    # - Smoother constraint -
    B = TTK.Button(Frame, text="Constraint", command=setConstraint,
                   image=iconics.PHOTO[8], compound=TK.RIGHT, padx=0, pady=0)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the constraint curve for smoother.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=8)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Constraint curve for smoother.')

    # - Remap -
    B = TTK.Button(Frame, text="Remap", command=remap)
    BB = CTK.infoBulle(parent=B, text='Remap a surface.')
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[2], 'QUAD', 'TRI', 'STRUCT')
    BB = CTK.infoBulle(parent=B, text='Elements in remap surface.')
    B.grid(row=4, column=1, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkCartWrap')
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
def saveApp():
    CTK.PREFS['tkCartWrapDensity'] = VARS[0].get()
    CTK.PREFS['tkCartWrapSmooth'] = VARS[1].get()
    CTK.PREFS['tkCartWrapElts'] = VARS[2].get()
    CTK.PREFS['tkCartWrapConsStrength'] = VARS[4].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('1.')
    VARS[1].set('10')
    VARS[2].set('QUAD')
    VARS[4].set('1.')
    CTK.PREFS['tkCartWrapDensity'] = VARS[0].get()
    CTK.PREFS['tkCartWrapSmooth'] = VARS[1].get()
    CTK.PREFS['tkCartWrapElts'] = VARS[2].get()
    CTK.PREFS['tkCartWrapConsStrength'] = VARS[4].get()
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
    (win, menu, file, tools) = CTK.minimal('tkCartWrap '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
