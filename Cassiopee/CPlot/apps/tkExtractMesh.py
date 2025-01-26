# - tkExtractMesh -
"""Extract mesh by interpolation."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Post.PyTree as P
import Converter.Internal as Internal
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setSources():
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
    VARS[0].set(selected)

#==============================================================================
def extractMesh():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Source meshes
    name = VARS[0].get()
    names = name.split(';')
    sources = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): sources.append(z)

    # Maillage a interpoler
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    meshes = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        meshes.append(z)

    fail = False
    try:
        meshes = P.extractMesh(sources, meshes, order=2, tol=1.e-6)
    except Exception as e:
        fail = True
        Panels.displayErrors([0,str(e)], header='Error: extractMesh')
    # Remise dans l'arbre
    if (fail == False):
        c = 0
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            CTK.t[2][nob][2][noz] = meshes[c]; c += 1
        CTK.TXT.insert('START', 'Interpolation performed.\n')
    else:
        CTK.TXT.insert('START', 'Interpolation fails.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkExtractMesh  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Interpolate solution on\nanother mesh.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=4)
    Frame.columnconfigure(1, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkExtractMesh')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Source grids -
    V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Source meshes -
    B = TTK.Button(Frame, text="Source", command=setSources,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    BB = CTK.infoBulle(parent=B, text='Set source meshes for interpolation.')
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    BB = CTK.infoBulle(parent=B, text='Source meshes for interpolation.')
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)

    # - Extract mesh -
    B = TTK.Button(Frame, text="Extract mesh", command=extractMesh)
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract solution from given source meshes.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['PostNoteBook'].add(WIDGETS['frame'], text='tkExtractMesh')
    except: pass
    CTK.WIDGETS['PostNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['PostNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkExtractMesh '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
