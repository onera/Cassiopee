# - Octree mesher app -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Generator.PyTree as G
import Connector.PyTree as X
import Converter.Internal as Internal
import numpy
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

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
    VARS[6].set(selected)

#==============================================================================
def expandLayer():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    level = CTK.varsFromWidget(VARS[4].get(), type=2)
    if level == []:
        CTK.TXT.insert('START', 'Invalid level.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    else: level = level[0]
    CTK.saveTree()
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            z = G.expandLayer(z, level=level)
            CTK.replace(CTK.t, nob, noz, z)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    #C._fillMissingVariables(CTK.t)
    if not fail:
        CTK.TXT.insert('START', 'Level %d expanded.\n'%level)
    else:
        Panels.displayErrors(errors, header='Error: expandLayers')
        CTK.TXT.insert('START', 'Expand layer fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def octree2Struct():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    vmin = CTK.varsFromWidget(VARS[7].get(), type=2)
    if vmin == []:
        CTK.TXT.insert('START', 'Invalid number of points per structured zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    else: vmin = vmin[0]

    type = VARS[8].get()
    ext = 0; AMR = 0; optimized = 1; merged = 1
    if (type == 'Match'): ext = 0; AMR = 0
    elif (type == 'Overlap1'): ext = 1; AMR = 0
    elif (type == 'Overlap2'): ext = 2; AMR = 0
    elif (type == 'AMR0'): ext = 0; AMR = 1
    elif (type == 'AMR1'): ext = 1; AMR = 1
    elif (type == 'AMR2'): ext = 2; AMR = 2

    CTK.saveTree()
    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    fail = False
    zlist = []; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            zlist = G.octree2Struct(z, vmin=vmin, ext=ext, optimized=optimized, 
                                    merged=merged, AMR=AMR)    
        except Exception as e: 
            fail = True; errors += [0,str(e)]

    CTK.t = C.addBase2PyTree(CTK.t, 'CARTESIAN')
    bases = Internal.getNodesFromName1(CTK.t, 'CARTESIAN')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in zlist: CTK.add(CTK.t, nob, -1, i)

    #C._fillMissingVariables(CTK.t)
    if fail == False:
        CTK.TXT.insert('START', 'Structured octree generated.\n')
    else:
        Panels.displayErrors(errors, header='Error: octree2Struct')
        CTK.TXT.insert('START', 'octree2Struct fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def bodyFit():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # EquationDimension
    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); dim = 3

    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Blanking
    tp = C.newPyTree(['Base'])
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        tp[2][1][2].append(z)

    # surfaces des corps
    name = VARS[6].get()
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
    if (len(surfaces) == 0):
        CTK.TXT.insert('START', 'Invalid body surface.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Blank
    BM = numpy.zeros((1, 1), dtype=Internal.E_NpyInt); BM[0,0] = 1
    tp = X.blankCells(tp, [surfaces], blankingMatrix=BM,
                      blankingType='node_in', dim=dim)

    CTK.saveTree()

    # Back to tree
    c = 0
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = tp[2][1][2][c]
        CTK.t[2][nob][2][noz] = z
        c += 1

    # Snap
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            z = G.snapFront(z, surfaces, optimized=2)
            CTK.t[2][nob][2][noz] = z
        except Exception as e:
            fail = True; errors += [0,str(e)]

    #C._fillMissingVariables(CTK.t)
    if fail == False:
        CTK.TXT.insert('START', 'Snapped to body.\n')
    else:
        Panels.displayErrors(errors, header='Error: snapFront')
        CTK.TXT.insert('START', 'Snap fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def adaptInsideOctree():
    if CTK.t == []: return
    if (CTK.__MAINTREE__ <= 0):
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # EquationDimension
    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); dim = 3

    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    snearsarg = VARS[0].get()
    try: volmin = float(snearsarg); volmin = volmin**dim
    except:
        CTK.TXT.insert('START', 'snear is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

        snearsarg = snearsarg.split(';')
        for s in snearsarg: snears.append(float(s))

    # Octree meshes
    tp = C.newPyTree(['Base'])
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        tp[2][1][2].append(z)

     # surfaces des corps
    name = VARS[6].get()
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
    if (len(surfaces) == 0):
        CTK.TXT.insert('START', 'Invalid body surface.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Blank   
    CTK.saveTree()

    BM = numpy.zeros((1, 1), dtype=Internal.E_NpyInt); BM[0,0] = 1
    end = 0
    count = 0
    while end == 0:
        tp = X.blankCells(tp, [surfaces], blankingMatrix=BM,
                          blankingType='center_in', dim=dim)
        tp = G.getVolumeMap(tp)
        C._initVars(tp,'{centers:indicator}=({centers:cellN}<1.)*({centers:vol}>%20.16g)'%volmin)
        end = 1
        if C.getMaxValue(tp,'centers:indicator') == 1.: end = 0
        for noz in range(len(tp[2][1][2])):
            tp[2][1][2][noz] = G.adaptOctree(tp[2][1][2][noz],'centers:indicator', balancing=1)
        count += 1

    # Back to tree
    c = 0
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = tp[2][1][2][c]
        CTK.t[2][nob][2][noz] = z
        c += 1
    fail = False
    #C._fillMissingVariables(CTK.t)
    if fail == False:
        CTK.TXT.insert('START', 'Adapt octree computed.\n')
    else: 
        CTK.TXT.insert('START', 'Adapt octree failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()    

#==============================================================================
def hexaOctree():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    snearsarg = VARS[0].get()
    snears = []
    try: snears.append(float(snearsarg))
    except:
        snearsarg = snearsarg.split(';')
        for s in snearsarg: snears.append(float(s))

    dfar = VARS[1].get()
    try: dfar = float(dfar)
    except:
        CTK.TXT.insert('START', 'dfar is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    balancing = VARS[2].get()
    if balancing == 'Balanced': balancing = 1
    else: balancing = 0

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    surfs = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dim = Internal.getZoneDim(z)
        if (dim[4] == 2 or dim[4] == 1): surfs.append(z) # surfaces ou courbes

    if len(snears) == 1 :
        snears = [snears[0] for x in range(len(surfs))]
    elif len(snears) != len(surfs):
        CTK.TXT.insert('START', 'snear is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'OCTREE')
    b = Internal.getNodesFromName1(CTK.t, 'OCTREE')[0]
    fail = False
    try:
        o = G.octree(surfs, snears, dfar=dfar, balancing=balancing)
        nob = C.getNobOfBase(b[0], CTK.t)
        CTK.add(CTK.t, nob, -1, o)
    except Exception as e:
        fail = True; print('Error: octree: %s.'%str(e))
    #C._fillMissingVariables(CTK.t)
    if fail == False:
        CTK.TXT.insert('START', 'Hexa octree computed.\n')
    else: 
        CTK.TXT.insert('START', 'Hexa octree failed.\n')
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
                           text='tkOctre  [ + ]  e', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Create octrees.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkOctree')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Snears -
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)
    if 'tkOctreeSnear' in CTK.PREFS: V.set(CTK.PREFS['tkOctreeSnear'])
    # -1- Dfar -
    V = TK.StringVar(win); V.set('10.'); VARS.append(V)
    if 'tkOctreeDfar' in CTK.PREFS: V.set(CTK.PREFS['tkOctreeDfar'])
    # -2- Balanced/unbalanced -
    V = TK.StringVar(win); V.set('Balanced'); VARS.append(V)
    if 'tkOctreeBalance' in CTK.PREFS: V.set(CTK.PREFS['tkOctreeBalance'])
    # -3- Vmins
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    if 'tkOctreeVmin' in CTK.PREFS: V.set(CTK.PREFS['tkOctreeVmin'])
    # -4- Level to expand
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkOctreeExpand' in CTK.PREFS: V.set(CTK.PREFS['tkOctreeExpand'])
    # -5- Type of body fitting
    V = TK.StringVar(win); V.set('Snap'); VARS.append(V)
    # -6- Body surfaces for body fitting
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -7- Vmin (number of points for each structured block -
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    # -8- Type of extension between structured blocks -
    V = TK.StringVar(win); V.set('Overlap2'); VARS.append(V)

    # - Dfar -
    rown=0
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=3)
    B.grid(row=rown, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Distance to far boundaries.')

    # - Snears -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=3)
    B.grid(row=rown, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Local size near bodies (global or list).')

    # - Balanced/unbalanced -
    B = TTK.OptionMenu(Frame, VARS[2], 'Balanced', 'Unbalanced')
    B.grid(row=rown, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Make a balanced octree or not.')

    # - Create hexa octree -
    rown += 1
    B = TTK.Button(Frame, text="Create hexa octree", command=hexaOctree)
    B.grid(row=rown, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create an hexa octree around surfaces.')

    # - Adapt hexa octree -
    rown+=1
    B = TTK.Button(Frame, text="Surfaces", command=setSurface,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=rown, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Define refinement surfaces.')
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White')
    B.grid(row=rown, column=0, columnspan=2, sticky=TK.EW)
    rown+=1
    B = TTK.Button(Frame, text="Adapt octree inside surfaces", command=adaptInsideOctree)
    B.grid(row=rown, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Adapt octree inside surfaces.')

    # - Expand layer -
    rown+=1
    B = TTK.Button(Frame, text="Expand layer", command=expandLayer)
    B.grid(row=rown, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Expand one level.')
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=5)
    B.grid(row=rown, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Level to expand.')

    # - octree2Struct -
    rown+=1
    B = TTK.Button(Frame, text="Octree2Struct", command=octree2Struct)
    B.grid(row=rown, column=0, columnspan=1, sticky=TK.EW)
    B = TK.Entry(Frame, textvariable=VARS[7], background='White', width=5)
    B.grid(row=rown, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of points in each structured blocs.')
    B = TTK.OptionMenu(Frame, VARS[8], 'Match', 'Overlap1', 'Overlap2', 'AMR0', 'AMR1', 'AMR2')
    B.grid(row=rown, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of structured blocks.')

    # - Body surfaces -
    rown+=1
    B = TTK.Button(Frame, text="Bodies", command=setSurface,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=rown, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Define bodies.')
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White')
    B.grid(row=rown, column=0, columnspan=2, sticky=TK.EW)

    # - Adapt to body
    rown+=1
    B = TTK.Button(Frame, text="Body fit", command=bodyFit)
    B.grid(row=rown, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fit the octree to bodies.')
    B = TTK.OptionMenu(Frame, VARS[5], 'Snap')
    B.grid(row=rown, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of body fitting.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MeshNoteBook'].add(WIDGETS['frame'], text='tkOctree')
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
    CTK.PREFS['tkOctreeSnear'] = VARS[0].get()
    CTK.PREFS['tkOctreeDfar'] = VARS[1].get()
    CTK.PREFS['tkOctreeBalance'] = VARS[2].get()
    CTK.PREFS['tkOctreeVmin'] = VARS[3].get()
    CTK.PREFS['tkOctreeExpand'] = VARS[4].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('0.1')
    VARS[1].set('10.')
    VARS[2].set('Balanced')
    VARS[3].set('10')
    VARS[4].set('0')
    CTK.PREFS['tkOctreeSnear'] = VARS[0].get()
    CTK.PREFS['tkOctreeDfar'] = VARS[1].get()
    CTK.PREFS['tkOctreeBalance'] = VARS[2].get()
    CTK.PREFS['tkOctreeVmin'] = VARS[3].get()
    CTK.PREFS['tkOctreeExpand'] = VARS[4].get()
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
    (win, menu, file, tools) = CTK.minimal('tkOctree '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
