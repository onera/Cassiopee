# - Chimera app -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.PyTree as C
import Connector.PyTree as X
import Converter.Internal as Internal
import numpy
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setPriority(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    priority = VARS[9].get(); priority = int(priority)

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        b = CTK.t[2][nob]
        C._addChimera2Base(b, 'Priority', priority)
    CTK.TXT.insert('START', 'Prioriy set to %d.\n'%priority)
    
#==============================================================================
def setSurface():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
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
# Init le cellN a 1 pour la selection ou tout l'arbre
#==============================================================================
def initCellN():
    if CTK.t == []: return
    type = VARS[8].get()
    CTK.saveTree()
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        if type == 'node_in':
            C._initVars(CTK.t, 'cellN', 1.)
        else:
            C._initVars(CTK.t, 'centers:cellN', 1.)
    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            if type == 'node_in':
                C._initVars(CTK.t[2][nob][2][noz], 'cellN', 1)
            else:
                C._initVars(CTK.t[2][nob][2][noz], 'centers:cellN', 1.)
        #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'cellN variable init to 1.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    
#==============================================================================
# Applique les BC overlaps pour tout l'arbre
#==============================================================================
def applyOverlap():
    if CTK.t == []: return
    CTK.saveTree()
    depth = VARS[7].get(); depth = int(depth)
    CTK.t = X.applyBCOverlaps(CTK.t, depth=depth)
    CTK.t = X.setDoublyDefinedBC(CTK.t, depth=depth)
    CTK.TXT.insert('START', 'cellN variable modified near overlap BCs.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    
#==============================================================================
# Applique setHoleInterpolatedPoints pour tout l'arbre
#==============================================================================
def setHoleInterpolatedPoints():
    if CTK.t == []: return
    CTK.saveTree()
    depth = VARS[7].get(); depth = int(depth)
    CTK.t = X.setHoleInterpolatedPoints(CTK.t, depth=depth)
    CTK.TXT.insert('START', 'cellN variable modified near holes.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# blanking
# Blank la selection avec la surface fournie 
#==============================================================================
def blank():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # type
    type = VARS[8].get()

    # EquationDimension
    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); dim = 3

    # Blanking surfaces
    name = VARS[6].get()
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

    # Reglages XRay
    delta = float(VARS[1].get())
    tol = float(VARS[2].get())

    # Creation de l'arbre temporaire
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    t = C.newPyTree(['Base'])
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        t[2][1][2].append(z)
    # Create blanking Matrix
    BM = numpy.zeros((1, 1), numpy.int32); BM[0,0] = 1

    # BlankCells
    CTK.saveTree()
    t = X.blankCells(t, [surfaces], blankingMatrix=BM,
                     blankingType=type,
                     delta=delta, dim=dim, tol=tol)

    # Back to tree
    c = 0
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = t[2][1][2][c]
        CTK.t[2][nob][2][noz] = z
        c += 1

    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Blanking done.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    
#==============================================================================
# Optimize Overlap
#==============================================================================
def optimize():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    dw = 1
    if VARS[0].get() == 'inactive': dw = 0

    depth = VARS[7].get(); depth = int(depth)
    
    # Creation de l'arbre temporaire
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
    # copie du premier niveau sans les enfants
    t = ['tree', None, [], 'CGNSTree_t']
    for i in CTK.t[2]: t[2].append([i[0], i[1], [], i[3]])

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        t[2][nob][2].append(z)

    # Recuperation des priorites dans l'arbre
    bases = Internal.getBases(CTK.t)
    prios = []
    for b in bases:
        cont = Internal.getNodesFromName1(b, '.Solver#Chimera')
        if cont != []:
            prio = Internal.getNodesFromName3(cont[0], 'Priority')
            if prio != []: prios.append(b[0]); prios.append(int(prio[0][1][0,0]))

    # Optimisation du recouvrement
    CTK.saveTree()
    t = X.optimizeOverlap(t, double_wall=dw, priorities=prios)
    t = X.maximizeBlankedCells(t, depth=depth)

    c = [0 for x in xrange(len(bases))]
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        CTK.t[2][nob][2][noz] = t[2][nob][2][c[nob-1]]
        c[nob-1] += 1
        
    CTK.TXT.insert('START', 'Overlapping optimized.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# createOversetHoles
#==============================================================================
def createOversetHoles():
    if CTK.t == []: return
    CTK.t = X.cellN2OversetHoles(CTK.t)
    CTK.TXT.insert('START', 'OversetHoles nodes created.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkChimera', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Chimera connectivity.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    Frame.columnconfigure(3, weight=1)
    Frame.columnconfigure(4, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkChimera')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- double wall: active ? -
    V = TK.StringVar(win); V.set('inactive'); VARS.append(V)
    # -1- delta XRay -
    V = TK.StringVar(win); V.set('1.e-10'); VARS.append(V)    
    # -2- tolerance XRay -
    V = TK.StringVar(win); V.set('1.e-8'); VARS.append(V)
    # -3- Base name 1
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -4- Relation
    V = TK.StringVar(win); V.set('+'); VARS.append(V)
    # -5- Base name 2
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -6- Blanking surface
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -7- depth
    V = TK.StringVar(win); V.set('2'); VARS.append(V)
    # -8- Blanking type
    V = TK.StringVar(win); V.set('cell_intersect'); VARS.append(V)
    # -9- Priorite
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    
    r = 0 # row
    # - Number of layers (depth) -
    B = TTK.Label(Frame, text="Depth")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of interpolated layers.')
    B = TTK.Entry(Frame, textvariable=VARS[7], background='White', width=2)
    B.grid(row=r, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of interpolated layers.')

    # - Blanking type -
    B = TTK.OptionMenu(Frame, VARS[8], 'cell_intersect', 'cell_intersect_opt',
                       'center_in', 'node_in')
    B.grid(row=r, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Blanking type.')
    r += 1
    
    # - initCellN -
    B = TTK.Button(Frame, text="Init cellN", command=initCellN)
    B.grid(row=r, column=0, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                      text='Init the cellN to 1. in the tree or selection.')
    r += 1
    
    # - applyBCOverlap -
    B = TTK.Button(Frame, text="ApplyBCOverlaps", command=applyOverlap)
    B.grid(row=r, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Modify the cellN (=2) near BCOverlaps and doubly defined BCs.')
    
    # - setHoleInterpolatedPoints -
    B = TTK.Button(Frame, text="SetHoleInterpPoints", command=setHoleInterpolatedPoints)
    B.grid(row=r, column=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Modify the cellN (=2) near holes.')
    r += 1

    # - XRay delta  -
    B = TTK.Label(Frame, text="XRay delta")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='The created holes will expand of this value.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=r, column=1, sticky=TK.EW)

    # - Xray tol -
    B = TTK.Label(Frame, text="Tol")
    B.grid(row=r, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Two surface points separated of this value\nwill be considered as identical.')
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=r, column=3, sticky=TK.EW)
    r += 1

    # - Surface -
    B = TTK.Button(Frame, text="Surf", command=setSurface,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=r, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set blanking surfaces.\nThey must define ONE closed body.')
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White')
    BB = CTK.infoBulle(parent=B, text='Blanking surfaces.')
    B.grid(row=r, column=0, columnspan=3, sticky=TK.EW)
    r += 1
    
    # - blanking -
    B = TTK.Button(Frame, text="Blank cells", command=blank)
    B.grid(row=r, column=0, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Blank cells with with surface.')
    r += 1

    # - set priority -
    B = TTK.Button(Frame, text="Set Priority", command=setPriority)
    B.grid(row=r, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set a priority to a base.')
    B = TTK.OptionMenu(Frame, VARS[9], '0', '1', '2', '3')
    B.grid(row=r, column=3, sticky=TK.EW)
    r += 1

    # - double wall active ? for optimizeOverlap
    B = TTK.Label(Frame, text="DoubleWall")
    B.grid(row=r, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[0], 'inactive', 'active')
    B.grid(row=r, column=3, sticky=TK.EW)

    # - optimizeOverlap -
    B = TTK.Button(Frame, text="Optimize overlap", command= optimize)
    B.grid(row=r, column=0, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Optimize the overlapping')
    r += 1
    
    # - createOversetHoles -
    B = TTK.Button(Frame, text="Create OversetHoles nodes",
                   command= createOversetHoles)
    B.grid(row=r, column=0, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create the nodes OversetHoles in pyTree.')
    r += 1
    
#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW); updateApp()

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(event=None): return

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
    (win, menu, file, tools) = CTK.minimal('tkChimera '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
