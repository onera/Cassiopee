# - tkInit -
"""Initialize flow fields."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Dist2Walls.PyTree as DW

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def initWallDistance():
    if CTK.t == []: return
    bodies = C.extractBCOfType(CTK.t, 'BCWall')
    for c in range(len(bodies)):
        try: bodies[c] = C.node2ExtCenter(bodies[c])
        except: pass
    tb = C.newPyTree(['Base', 2]); tb[2][1][2] += bodies
    CTK.saveTree()
    try:
        CTK.t = DW.distance2Walls(CTK.t, tb, loc='centers', type='ortho')
        CTK.TXT.insert('START', 'Distance to wall computed.\n')
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: wallDistance')
        CTK.TXT.insert('START', 'Wall distance fails.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def initSolution():
    if CTK.t == []: return
    CTK.saveTree()
    state = Internal.getNodeFromType(CTK.t, 'ReferenceState_t')
    if state is None:
        CTK.TXT.insert('START', 'state is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Check for GoverningEquations
    eqs = Internal.getNodeFromType(CTK.t, 'GoverningEquations_t')
    Model = 'NSTurbulent'
    if eqs is not None: Model = Internal.getValue(eqs)

    zvars = ['Density', 'MomentumX', 'MomentumY', 'MomentumZ',
             'EnergyStagnationDensity']
    for v in zvars:
        node = Internal.getNodeFromName(state, v)
        if node is not None:
            val = float(node[1][0])
            C._initVars(CTK.t, 'centers:'+v, val)
        else:
            CTK.TXT.insert('START', v + ' is missing (tkState).\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return
    if Model == 'NSTurbulent':
        zvars = ['TurbulentSANuTildeDensity', 'TurbulentEnergyKineticDensity',
                 'TurbulentDissipationDensity']
        for v in zvars:
            node = Internal.getNodeFromName(state, v)
            if node is not None:
                val = float(node[1][0])
                C._initVars(CTK.t, 'centers:'+v, val)

    CTK.TXT.insert('START', 'Solution initialized.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkInit  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Init solution fields.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkInit')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Zone filter regexp -
    #V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Init solution -
    B = TTK.Button(Frame, text="Initialize solution from state",
                   command=initSolution)
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Initialize the solution in centers (require state).')

    # - Init wall distances -
    B = TTK.Button(Frame, text="Initialize wall distances",
                   command=initWallDistance)
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Compute wall distances.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SolverNoteBook'].add(WIDGETS['frame'], text='tkInit')
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
    (win, menu, file, tools) = CTK.minimal('tkInit '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
