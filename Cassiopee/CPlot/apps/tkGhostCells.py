# - tkGhostCells -
"""Add or remove ghost cells."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def addGhostCells(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    N = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(N) != 1:
        CTK.TXT.insert('START', 'Number of ghost cell layers is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    N = N[0]

    CTK.saveTree()

    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        zones += [z]

    zones = Internal.addGhostCells(CTK.t, zones, N, adaptBCs=1)
    c = 0
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        CTK.replace(CTK.t, nob, noz, zones[c])
        c += 1

    CTK.TXT.insert('START', 'Ghost cells added.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def rmGhostCells(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    N = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(N) != 1:
        CTK.TXT.insert('START', 'Number of ghost cell layers is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    N = N[0]
    CTK.saveTree()

    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        zones += [z]

    zones = Internal.rmGhostCells(CTK.t, zones, N, adaptBCs=1)
    c = 0
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        CTK.replace(CTK.t, nob, noz, zones[c])
        c += 1

    CTK.TXT.insert('START', 'Ghost cells removed.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkGhostCells  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='My personal applet.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkGhostCells')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- number of layers of ghost cells -
    V = TK.StringVar(win); V.set('2'); VARS.append(V)
    # -1- methode pour les coins -
    V = TK.StringVar(win); V.set('None'); VARS.append(V)

    # - Number of ghost cell layers -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of ghost cell layers.')

    # - methode pour les coins -
    B = TTK.OptionMenu(Frame, VARS[1], 'None')
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)

    # - add ghost cells -
    B = TTK.Button(Frame, text="Add", command=addGhostCells)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Add layers of ghost cells.')

    # - rm ghost cells -
    B = TTK.Button(Frame, text="Rm", command=rmGhostCells)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Rm layers of ghost cells.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['BlockNoteBook'].add(WIDGETS['frame'], text='tkGhostCells')
    except: pass
    CTK.WIDGETS['BlockNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['BlockNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkGhostCells '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
