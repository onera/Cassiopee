# - tkDist2Walls -
"""Compute distance to walls."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Dist2Walls.PyTree as DTW
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setSurfaces():
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
def compute():
    if CTK.t == []: return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    surf = 0; walls = []
    name = VARS[1].get()
    if name == '': names = []
    else: names = name.split(';')
    for v in names:
        surf = 1
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: walls.append(z)

    if surf == 0:
        walls = C.extractBCOfType(CTK.t, 'BCWall')
        walls += C.extractBCOfType(CTK.t, 'BCWallViscous')
        walls += C.extractBCOfType(CTK.t, 'BCWallInviscid')

    tp = C.newPyTree(['Base'])
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        tp[2][1][2].append(z)

    try:
        if VARS[2].get() == 'absolute': signed = 0
        else: signed = 1
        tp = DTW.distance2Walls(tp, walls, type=VARS[0].get(), loc=VARS[3].get(),
                                signed=signed)
        c = 0
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            CTK.t[2][nob][2][noz] = tp[2][1][2][c]
            c += 1
        #C._fillMissingVariables(CTK.t)
        CTK.TKTREE.updateApp()
        CTK.display(CTK.t)
        CTK.TXT.insert('START', 'Distance to walls computed.\n')
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: dist2Walls')
        CTK.TXT.insert('START', 'Distance to walls failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkDist2Walls  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Compute wall distance.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=4)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkDist2Walls')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Type de distance -
    V = TK.StringVar(win); V.set('ortho'); VARS.append(V)
    if 'tkDist2WallsType' in CTK.PREFS:
        V.set(CTK.PREFS['tkDist2WallsType'])
    # -1- Surfaces -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -2- Signed ou absolute -
    V = TK.StringVar(win); V.set('absolute'); VARS.append(V)
    if 'tkDist2WallsSigned' in CTK.PREFS:
        V.set(CTK.PREFS['tkDist2WallsSigned'])
    # -3- Vars location -
    V = TK.StringVar(win); V.set('nodes'); VARS.append(V)
    if 'tkVariablesLoc' in CTK.PREFS:
        V.set(CTK.PREFS['tkVariablesLoc'])

    # - Surfaces -
    B = TTK.Button(Frame, text="Surfaces", command=setSurfaces)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set surfaces for distances computation.\nIf not set, BCWalls are used.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White')
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)

    # - Algorithm -
    B = TTK.OptionMenu(Frame, VARS[0], 'ortho', 'mininterf')
    B.grid(row=1, column=1, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[2], 'absolute', 'signed')
    B.grid(row=1, column=0, sticky=TK.EW)

    # - Compute -
    B = TTK.Button(Frame, text="Compute", command=compute)
    B.grid(row=2, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[3], 'nodes', 'centers')
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Compute the wall distance.\nTree is modified.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SolverNoteBook'].add(WIDGETS['frame'], text='tkDist2Walls')
    except: pass
    CTK.WIDGETS['SolverNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['SolverNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkDist2WallsType'] = VARS[0].get()
    CTK.PREFS['tkDist2WallsSigned'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('ortho')
    VARS[2].set('absolute')
    CTK.PREFS['tkDist2WallsType'] = VARS[0].get()
    CTK.PREFS['tkDist2WallsSigned'] = VARS[2].get()
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
    (win, menu, file, tools) = CTK.minimal('tkDist2Walls '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
