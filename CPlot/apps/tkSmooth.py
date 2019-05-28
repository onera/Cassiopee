# - mesh smoother -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Transform.PyTree as T
import Generator.PyTree as G
import Post.PyTree as P
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
    VARS[1].set(selected)
    
#==============================================================================
def smooth():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    
    smooth = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(smooth) != 1:
        CTK.TXT.insert('START', 'Smooth iter is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    smooth = smooth[0]

    eps = CTK.varsFromWidget(VARS[4].get(), type=1)
    if len(eps) != 1:
        CTK.TXT.insert('START', 'Eps is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    eps = eps[0]

    ntype = VARS[8].get()
    if ntype == 'Volume': ntype = 0
    elif ntype == 'Scale': ntype = 1
    elif ntype == 'Taubin': ntype = 2
    else: ntype = 0

    # Constraint strength
    strength = CTK.varsFromWidget(VARS[2].get(), type=1)
    if len(strength) != 1:
        CTK.TXT.insert('START', 'Constraint strength is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    strength = strength[0]

    # Constraints
    fixedConstraints = []; projConstraints = []

    name = VARS[1].get()
    names = name.split(';')    
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        base = Internal.getNodeFromName1(CTK.t, sname[0])
        if base is not None:
            nodes = Internal.getNodesFromType1(base, 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): fixedConstraints.append(z)

    CTK.saveTree()
    
    # Get all zones
    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        zones.append(z)

    # Mesh unique
    try:
        A = C.convertArray2Tetra(zones)
        A = T.join(A); A = G.close(A)
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: smooth')
        CTK.TXT.insert('START', 'Some zones are invalid for smoothing.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # pb surfacique ou volumique?
    dims = Internal.getZoneDim(A)
    pbDim = 3
    if dims[3] == 'TRI': pbDim = 2
        
    # Keep external faces
    if VARS[3].get() == 1 or pbDim == 3: 
        try: ext = P.exteriorFaces(A)
        except: ext = []
    if VARS[3].get() == 1 and ext != []: fixedConstraints.append(ext)

    # Keep sharp edges
    if VARS[5].get() == 1:
        angle = CTK.varsFromWidget(VARS[6].get(), type=1)
        if len(angle) != 1:
            CTK.TXT.insert('START', 'Sharp edges angle is incorrect.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        angle = angle[0]
        if pbDim == 2:
            try:
                sharp = P.sharpEdges(A, angle)
                fixedConstraints += sharp
            except: pass
        else:
            try:
                sharp = P.sharpEdges(ext, angle)
                fixedConstraints += sharp
            except: pass

    # Project on surface
    Pj = VARS[7].get()
    if pbDim == 2: projSurf = A
    else:
        # projSurf = ext
        # Pour l'instant, on ne sait pas projeter un maillage volumique
        # sur une surface
        Pj = 0
    
    # Smooth
    fail = False
    try:
        if Pj == 0:
            zones = T.smooth(zones, eps=eps, niter=smooth, type=ntype,
                             fixedConstraints=fixedConstraints, 
                             projConstraints=projConstraints, delta=strength)
        else:
            for s in xrange(smooth):
                zones = T.smooth(zones, eps=eps, niter=2, type=ntype,
                                 fixedConstraints=fixedConstraints, 
                                 projConstraints=projConstraints, 
                                 delta=strength)
                zones = T.projectOrtho(zones, [projSurf])
    except Exception as e:
        fail = True
        Panels.displayErrors([0,str(e)], header='Error: smooth')

    if fail == False:
        c = 0
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            a = zones[c]; c += 1
            CTK.replace(CTK.t, nob, noz, a)
            CTK.TXT.insert('START', 'Mesh smoothed.\n')
    else:
        CTK.TXT.insert('START', 'Smooth fails.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkSmooth', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Smooth meshes.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=0)
    Frame.columnconfigure(2, weight=0)
    Frame.columnconfigure(3, weight=0)
    Frame.columnconfigure(4, weight=0)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkSmooth')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Smoother niter -
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    if 'tkSmoothIter' in CTK.PREFS: V.set(CTK.PREFS['tkSmoothIter'])
    # -1- Constraint
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -2- Constraint strength
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)
    if 'tkSmoothConsStrength' in CTK.PREFS: 
        V.set(CTK.PREFS['tkSmoothConsStrength'])
    # -3- Constraint external faces
    V = TK.IntVar(win); V.set(1); VARS.append(V)
    # -4- smooth eps
    V = TK.StringVar(win); V.set('0.5'); VARS.append(V)
    if 'tkSmoothEps' in CTK.PREFS: 
        V.set(CTK.PREFS['tkSmoothEps'])
    # -5- Constraint sharp edges
    V = TK.IntVar(win); V.set(0); VARS.append(V)
    # -6- Sharp edges detection angle
    V = TK.StringVar(win); V.set('30.'); VARS.append(V)
    if 'tkSmoothSharpAngle' in CTK.PREFS: 
        V.set(CTK.PREFS['tkSmoothSharpAngle'])
    # -7- Project on surface
    V = TK.IntVar(win); V.set(0); VARS.append(V)
    # -8- Type de smoothing
    V = TK.StringVar(win); V.set('Volume'); VARS.append(V)
    if 'tkSmoothType' in CTK.PREFS: 
        V.set(CTK.PREFS['tkSmoothType'])

    # - Smoother power -
    #B = TTK.Label(Frame, text="Eps")
    #B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[8], 'Volume', 'Scale', 'Taubin')
    BB = CTK.infoBulle(parent=B, text='Smoothing algorithm.')
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=4)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Smoother power.')
    B = TTK.Label(Frame, text="Iter")
    B.grid(row=0, column=3, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=4)
    B.grid(row=0, column=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of smoother iterations.')

    # - Constraint strength -
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=4)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Strength of constraints.')

    # - Keep external faces -
    B = TTK.Checkbutton(Frame, text='EF', variable=VARS[3])
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Keep external faces.')

    # - Keep sharp edges faces -
    B = TTK.Checkbutton(Frame, text='SE', variable=VARS[5])
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Keep sharp edges.')

    # - Project on surface
    B = TTK.Checkbutton(Frame, text='PJ', variable=VARS[7])
    B.grid(row=1, column=3, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Project on surface.')

    # - Sharp edges detection
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White', width=4)
    B.grid(row=1, column=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Sharp edges detection angle.')

    # - Smoother constraint -
    B = TTK.Button(Frame, command=setConstraint,
                   image=iconics.PHOTO[8], padx=0, pady=0)
    BB = CTK.infoBulle(parent=B, text='Set constraint curves for smoother (points dont move).')
    B.grid(row=2, column=4, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=8)
    B.grid(row=2, column=0, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Constraint curves for smoother (points dont move).')
    
    # - Smooth -
    B = TTK.Button(Frame, text="Smooth", command=smooth)
    BB = CTK.infoBulle(parent=B, text='Smooth mesh.')
    B.grid(row=3, column=0, columnspan=5, sticky=TK.EW)
    
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
    CTK.PREFS['tkSmoothIter'] = VARS[0].get()
    CTK.PREFS['tkSmoothConsStrength'] = VARS[2].get()
    CTK.PREFS['tkSmoothEps'] = VARS[4].get()
    CTK.PREFS['tkSmoothSharpAngle'] = VARS[6].get()
    CTK.PREFS['tkSmoothType'] = VARS[8].get()
    CTK.savePrefFile()
    
#==============================================================================
def resetApp():
    VARS[0].set('10')
    VARS[2].set('0.1')
    VARS[4].set('0.5')
    VARS[6].set('30.')
    VARS[8].set('Volume')
    CTK.PREFS['tkSmoothIter'] = VARS[0].get()
    CTK.PREFS['tkSmoothConsStrength'] = VARS[2].get()
    CTK.PREFS['tkSmoothEps'] = VARS[4].get()
    CTK.PREFS['tkSmoothSharpAngle'] = VARS[6].get()
    CTK.PREFS['tkSmoothType'] = VARS[8].get()
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
    (win, menu, file, tools) = CTK.minimal('tkSmooth '+C.__version__)

    createApp(win); activateApp()

    # - Main loop -
    win.mainloop()
