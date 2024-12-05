# - elsA solver app for elsAxdt -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def export(event=None):
    if CTK.t == []: return
    import Converter.elsAProfile
    filename = VARS[0].get()
    tp = Converter.elsAProfile.convert2elsAxdt(CTK.t)
    C.convertPyTree2File(tp, filename)
    CTK.TXT.insert('START', 'File '+filename+' exported.\n')

#==============================================================================
def adapt(event=None):
    if CTK.t == []: return
    import Converter.elsAProfile
    CTK.t = Converter.elsAProfile.convert2elsAxdt(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.TXT.insert('START', 'Tree adapted for elsA.CGNS.\n')

#==============================================================================
def elsAHybridNode(event=None):
    if CTK.t == []: return
    import Converter.Internal as Internal
    Internal._createElsaHybrid(CTK.t, VARS[1].get())
    CTK.TKTREE.updateApp()
    CTK.TXT.insert('START', 'elsAHybrid node created.\n')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkElsaSolver  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Export to Onera\nelsA fluid solver.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkElsaSolver')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- file name for export -
    V = TK.StringVar(win); V.set('elsA.cgns'); VARS.append(V)
    if 'tkElsaSolverFile' in CTK.PREFS: 
        V.set(CTK.PREFS['tkElsaSolverFile'])
    # -1- Method for createElsaHybrid
    V = TK.IntVar(win); V.set(0); VARS.append(V)

    # - Export -
    B = TTK.Button(Frame, text="Export to", command=export)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Export tree to a CGNS file suitable for elsA.CGNS.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=0, column=1, sticky=TK.EW)
    B.bind('<Return>', export)

    # - Adapt in memory tree - 
    B = TTK.Button(Frame, text="Adapt current tree", command=adapt)
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Adapt in memory tree to make it suitable for elsA.CGNS.')

    # - Create elsA Hybrid node - 
    B = TTK.Button(Frame, text="Create elsA Hybrid Node", command=elsAHybridNode)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create elsAHybrid node for unstructured mesh')
    B = TK.Entry(Frame, textvariable=VARS[1], background='White')
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Method: 0 (elsA < 3.8.01), 1 (elsA >= 3.8.1)')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SolverNoteBook'].add(WIDGETS['frame'], text='tkElsaSolver')
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
    CTK.PREFS['tkElsaSolverFile'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('elsA.cgns')
    CTK.PREFS['tkElsaSolverFile'] = VARS[0].get()
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
    (win, menu, file, tools) = CTK.minimal('tkElsaSolver '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
