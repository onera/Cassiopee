# - PolyLineMesher app -
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Generator.PyTree as G

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def generatePLM():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) == 0:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    hf = CTK.varsFromWidget(VARS[1].get(), type=1)
    if len(hf) != 1:
        CTK.TXT.insert('START', 'First cell height is incorrect.\n'); return
    hf = hf[0]
    h = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(h) != 1:
        CTK.TXT.insert('START', 'Mesh height is incorrect.\n'); return
    h = h[0]
    density = CTK.varsFromWidget(VARS[2].get(), type=1)
    if len(density) != 1:
        CTK.TXT.insert('START', 'Grid point density is incorrect.\n'); return
    density = density[0]

    try:
        CTK.saveTree()
        CTK.t = C.addBase2PyTree(CTK.t, 'MESHES')
        bases = Internal.getNodesFromName1(CTK.t, 'MESHES')
        gnob = C.getNobOfBase(bases[0], CTK.t)
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            B = G.polyLineMesher(z, h, hf, density)
            for i in B[0]: CTK.add(CTK.t, gnob, -1, i)
        CTK.TXT.insert('START', 'PLM mesh created.\n')
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: PLM')
        CTK.TXT.insert('START', 'PLM mesh failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                          text='tkPLM  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='PLM 2D automatic mesher.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkPLM')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- height of mesh -
    V = TK.StringVar(win); V.set('50.'); VARS.append(V)
    # -1- height of first cell of mesh -
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)
    # -2- grid point density -
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)

    # - hf -
    B = TK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Height of first cell.')

    # - h -
    B = TK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Height of generated mesh.')

    # - density -
    B = TK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Grid point density.')

    # - Generate -
    B = TK.Button(Frame, text="Generate", command=generatePLM)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Generate PLM grid.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.NSEW)

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
    (win, menu, file, tools) = CTK.minimal('tkPLM '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
