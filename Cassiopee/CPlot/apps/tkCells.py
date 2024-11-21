# - tkCells -
"""Cell mesh manipulations."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Post.PyTree as P
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import time

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def selectCells(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    formula = VARS[0].get()
    strict = VARS[1].get()
    strict = int(strict)
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = P.selectCells(z, formula, strict=strict)
        CTK.replace(CTK.t, nob, noz, z)

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TXT.insert('START', 'Cells selected.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Suppress cells
#==============================================================================
def suppressCells():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    W = WIDGETS['suppress']
    if not CTK.__BUSY__:
        CPlot.unselectAllZones()
        CTK.__BUSY__ = True
        TTK.sunkButton(W)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            l = []
            while l == []:
                nz = CPlot.getSelectedZone()
                l = CPlot.getActivePointIndex()
                CPlot.unselectAllZones()
                time.sleep(CPlot.__timeStep__)
                W.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                nob = CTK.Nb[nz]+1
                noz = CTK.Nz[nz]
                CTK.saveTree()
                z = CTK.t[2][nob][2][noz]
                C._initVars(z, 'centers:__tag__', 1)
                C.setValue(z, 'centers:__tag__', l[1], 0)
                try:
                    z = P.selectCells2(z, 'centers:__tag__')
                    CTK.replace(CTK.t, nob, noz, z)
                except: pass
                CTK.TKTREE.updateApp()
                CPlot.render()
        CTK.__BUSY__ = False
        TTK.raiseButton(W)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(W)
        CPlot.setState(cursor=0)

#==============================================================================
# Refine cells
#==============================================================================
def refineCells():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    W = WIDGETS['refine']
    if not CTK.__BUSY__:
        CPlot.unselectAllZones()
        CTK.__BUSY__ = True
        TTK.sunkButton(W)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            l = []
            while l == []:
                nz = CPlot.getSelectedZone()
                l = CPlot.getActivePointIndex()
                CPlot.unselectAllZones()
                time.sleep(CPlot.__timeStep__)
                W.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                nob = CTK.Nb[nz]+1
                noz = CTK.Nz[nz]
                CTK.saveTree()
                z = CTK.t[2][nob][2][noz]
                C._initVars(z, 'centers:__tag__', 0)
                C.setValue(z, 'centers:__tag__', l[1], 1)
                try:
                    z = P.refine(z, '__tag__')
                    CTK.replace(CTK.t, nob, noz, z)
                except: pass
                CTK.TKTREE.updateApp()
                CPlot.render()
        CTK.__BUSY__ = False
        TTK.raiseButton(W)
        CPlot.setState(cursor=0)
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(W)
        CPlot.setState(cursor=0)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCells  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Per cell operation\non meshes.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=4)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkCells')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- selection formula -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- strict -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)

    # - Cells suppress/refine -
    B = TTK.Button(Frame, text="Suppress cell mode", command=suppressCells)
    B.grid(row=0, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Suppress selected cell.')
    WIDGETS['suppress'] = B
    B = TTK.Button(Frame, text="Refine cell mode", command=refineCells)
    B.grid(row=1, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refine selected cell (TRI mesh only).')
    WIDGETS['refine'] = B

    # - Select cells -
    B = TTK.Button(Frame, text="SelectCells", command=selectCells)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Select cells that matches formula.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.bind('<Return>', selectCells)
    BB = CTK.infoBulle(parent=B, text='Selection formula.')
    B.grid(row=2, column=1, sticky=TK.EW)
    B = TTK.Checkbutton(Frame, text='', variable=VARS[1])
    BB = CTK.infoBulle(parent=B, text='Strict/not strict selection.')
    B.grid(row=2, column=2, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MeshNoteBook'].add(WIDGETS['frame'], text='tkCells')
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
    (win, menu, file, tools) = CTK.minimal('tkCells '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
