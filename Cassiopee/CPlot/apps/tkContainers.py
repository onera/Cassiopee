# - tkContainers -
"""Data container setup (FlowSolution, GridCoordinates)."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.PyTree as C
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Cree une liste des containers GridCoordinates_t
#==============================================================================
def updateGridCoordinates(event=None):
    f = Internal.getNodesFromType3(CTK.t, 'GridCoordinates_t')
    vars = ['GridCoordinates']
    for b in f: vars.append(b[0])
    seen = set()
    seen_add = seen.add
    vars = [x for x in vars if x not in seen and not seen_add(x)]
    if 'GridCoordinates' in WIDGETS:
        WIDGETS['GridCoordinates']['values'] = vars

#==============================================================================
# Cree une liste des containers FlowSolution_t + Vertex
#==============================================================================
def updateFlowSolution(event=None):
    f = Internal.getNodesFromType3(CTK.t, 'FlowSolution_t')
    zvars = ['FlowSolution']
    for b in f:
        loc = Internal.getNodesFromType1(b, 'GridLocation_t')
        if loc == []: zvars.append(b[0])
        else:
            loc = loc[0]; v = Internal.getValue(loc)
            if v == 'Vertex': zvars.append(b[0])
    seen = set()
    seen_add = seen.add
    zvars = [x for x in zvars if x not in seen and not seen_add(x)]
    if 'FlowSolution' in WIDGETS:
        WIDGETS['FlowSolution']['values'] = zvars

#==============================================================================
# Cree une liste des containers FlowSolution_t + CellCenter
#==============================================================================
def updateFlowSolutionCenters(event=None):
    f = Internal.getNodesFromType3(CTK.t, 'FlowSolution_t')
    zvars = ['FlowSolution#Centers']
    for b in f:
        loc = Internal.getNodesFromType1(b, 'GridLocation_t')
        if loc != []:
            loc = loc[0]; v = Internal.getValue(loc)
            if v == 'CellCenter': zvars.append(b[0])
    seen = set()
    seen_add = seen.add
    zvars = [x for x in zvars if x not in seen and not seen_add(x)]
    if 'FlowSolutionCenters' in WIDGETS:
        WIDGETS['FlowSolutionCenters']['values'] = zvars

#==============================================================================
def setNames(event=None):
    Internal.__GridCoordinates__ = VARS[0].get()
    Internal.__FlowSolutionNodes__ = VARS[1].get()
    Internal.__FlowSolutionCenters__ = VARS[2].get()
    CTK.TXT.insert('START', 'Container names set.\n')
    CTK.display(CTK.t)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkContainers  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage container names.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkContainers')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- GridCoordinates container -
    V = TK.StringVar(win); V.set('GridCoordinates'); VARS.append(V)
    if 'GridCoordinatesContainer' in CTK.PREFS:
        V.set(CTK.PREFS['GridCoordinatesContainer'])

    # -1- FlowSolutionNodes container -
    V = TK.StringVar(win); V.set('FlowSolution'); VARS.append(V)
    if 'FlowSolutionNodesContainer' in CTK.PREFS:
        V.set(CTK.PREFS['FlowSolutionNodesContainer'])
    # -2- FlowSolutionCenters container -
    V = TK.StringVar(win); V.set('FlowSolution#Centers'); VARS.append(V)
    if 'FlowSolutionCentersContainer' in CTK.PREFS:
        V.set(CTK.PREFS['FlowSolutionCentersContainer'])

    # - GridCoordinates -
    B = TTK.Label(Frame, text="GridCoordinates")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='This name will be used to find coordinates node in zones.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.Entry(F, textvariable=VARS[0], background='White')
        B.grid(sticky=TK.EW)
        F.bind('<Return>', setNames)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['GridCoordinates'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[0],
                         values=[], state='normal')
        B.bind("<<ComboboxSelected>>", setNames)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateGridCoordinates)
        B.bind('<Return>', setNames)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['GridCoordinates'] = B

    # - FlowSolutionNodes -
    B = TTK.Label(Frame, text="FlowSolutionNodes")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='This name will be used to find solution node in zones.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.Entry(F, textvariable=VARS[1], background='White')
        B.grid(sticky=TK.EW)
        F.bind('<Return>', setNames)
        F.grid(row=1, column=1, sticky=TK.EW)
        WIDGETS['FlowSolution'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[1],
                         values=[], state='normal')
        B.bind("<<ComboboxSelected>>", setNames)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFlowSolution)
        B.bind('<Return>', setNames)
        F.grid(row=1, column=1, sticky=TK.EW)
        WIDGETS['FlowSolution'] = B

    # - FlowSolutionCenters -
    B = TTK.Label(Frame, text="FlowSolutionCenters")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='This name will be used to find centers solution node in zones.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.Entry(F, textvariable=VARS[2], background='White')
        B.grid(sticky=TK.EW)
        F.bind('<Return>', setNames)
        F.grid(row=2, column=1, sticky=TK.EW)
        WIDGETS['FlowSolutionCenters'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[2],
                         values=[], state='normal')
        B.bind("<<ComboboxSelected>>", setNames)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFlowSolutionCenters)
        B.bind('<Return>', setNames)
        F.grid(row=2, column=1, sticky=TK.EW)
        WIDGETS['FlowSolutionCenters'] = B

    # - set -
    B = TTK.Button(Frame, text="Set", command=setNames)
    BB = CTK.infoBulle(parent=B, text='Change the container names used by Cassiopee functions.')
    B.grid(row=3, column=0, columnspan=2, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['StateNoteBook'].add(WIDGETS['frame'], text='tkContainers')
    except: pass
    CTK.WIDGETS['StateNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['StateNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp():
    VARS[0].set(Internal.__GridCoordinates__)
    VARS[1].set(Internal.__FlowSolutionNodes__)
    VARS[2].set(Internal.__FlowSolutionCenters__)

#==============================================================================
def saveApp():
    CTK.PREFS['GridCoordinatesContainer'] = VARS[0].get()
    CTK.PREFS['FlowSolutionNodesContainer'] = VARS[1].get()
    CTK.PREFS['FlowSolutionCentersContainer'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('GridCoordinates')
    VARS[1].set('FlowSolution')
    VARS[2].set('FlowSolution#Centers')
    CTK.PREFS['GridCoordinatesContainer'] = VARS[0].get()
    CTK.PREFS['FlowSolutionNodesContainer'] = VARS[1].get()
    CTK.PREFS['FlowSolutionCentersContainer'] = VARS[2].get()
    CTK.savePrefFile()
    setNames()

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
    (win, menu, file, tools) = CTK.minimal('tkContainers '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
