# - tkTetraMesher -
"""Generate TETRA meshes."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Generator.PyTree as G
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def tetraMesher():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    if VARS[0].get() == 'netgen': algo = 0
    else: algo = 1

    CTK.setCursor(2, WIDGETS['tetraMesher'])
    CTK.saveTree()
    out = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        out.append(CTK.t[2][nob][2][noz])

    try:
        mesh = G.tetraMesher(out, algo=algo, grading=0.3, recoverBC=True)
        CTK.t = C.addBase2PyTree(CTK.t, 'MESHES')
        bases = Internal.getNodesFromName1(CTK.t, 'MESHES')
        nob = C.getNobOfBase(bases[0], CTK.t)
        CTK.add(CTK.t, nob, -1, mesh)
        CTK.TXT.insert('START', 'Tetra mesh created.\n')
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: TetraMesher')
        CTK.TXT.insert('START', 'Tetra mesh failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    CTK.setCursor(0, WIDGETS['tetraMesher'])

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkTetraMesher  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    #Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkTetraMesher')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Mesher type -
    V = TK.StringVar(win); V.set('tetgen'); VARS.append(V)
    if 'tkTetraMesherType' in CTK.PREFS: V.set(CTK.PREFS['tkTetraMesherType'])

    # - mesher menu -
    #B = TTK.OptionMenu(Frame, VARS[0], 'netgen', 'tetgen')
    #B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)

    # - Run -
    B = TTK.Button(Frame, text="tetraMesher", command=tetraMesher)
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    WIDGETS['tetraMesher'] = B
    BB = CTK.infoBulle(parent=B, text='Mesh with TETRAs or TRIs.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MeshNoteBook'].add(WIDGETS['frame'], text='tkTetraMesher')
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
    CTK.PREFS['tkTetraMesherType'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('netgen')
    CTK.PREFS['tkTetraMesherType'] = VARS[0].get()
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
    (win, menu, file, tools) = CTK.minimal('tkTetraMesher '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
