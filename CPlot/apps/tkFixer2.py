# - Fix holes in meshes -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T
import Intersector.PyTree as XOR

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Fix one gap
#==============================================================================
def fixGap():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # get bump factor
    factor = 2*WIDGETS['bump'].get() / 100.
    
    # Contours
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
    contours = []
    for nz in nzs:
      nob = CTK.Nb[nz]+1
      noz = CTK.Nz[nz]
      contours.append(CTK.t[2][nob][2][noz])

    contours = C.convertArray2Tetra(contours)
    contours = T.join(contours)
    contours = G.close(contours)
    contours = T.splitManifold(contours)

    fail = False; out = []
    for contour in contours:
        try:
          p = G.fittingPlaster(contour, bumpFactor=factor)
          b = G.gapfixer(contour, p)
          out.append(b)
        except Exception as e:
          fail = True
          Panels.displayErrors([0,str(e)], header='Error: gapfixer on %s.'%contour[0])
  
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'FIXED', 2)
    base = Internal.getNodesFromName1(CTK.t, 'FIXED')[0]
    nob = C.getNobOfBase(base, CTK.t)
    for b in out:
       CTK.add(CTK.t, nob, -1, b)
       #CTK.add(CTK.t, nob, -1, p)
    if not fail:
       #C._fillMissingVariables(CTK.t)
       CTK.TXT.insert('START', 'Gap fixed.\n')
    else:
      CTK.TXT.insert('START', 'Gap fixing fails.\n')
      CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Fix multiple gaps, eliminate overlap in patches
#==============================================================================
def fixGaps():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    
    # get mode
    modeString = VARS[0].get()
    mode = 2; coplanar = 0 # unknown
    if modeString == 'Centers': mode = 1; coplanar = 0
    elif modeString == 'Nodes': mode = 0; coplanar = 0
    elif modeString == 'Slice': mode = 0; coplanar = 1
    
    # Patches
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    patches = []
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        patches.append(CTK.t[2][nob][2][noz])
    nob0 = CTK.Nb[nzs[0]]+1
    patches = C.convertArray2Tetra(patches)
    patches = G.close(patches)
    b = G.gapsmanager(patches, mode=mode, coplanar=coplanar)
    CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
    CTK.t[2][nob0][2] += b
    
    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Gaps fixed.\n')    
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def conformUnstr():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    tol = CTK.varsFromWidget(VARS[1].get(), type=1)
    if len(tol) != 1:
        CTK.TXT.insert('START', 'Tolerance is incorrect.\n') 
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = tol[0]

    patches = []
    CTK.saveTree()

    if len(nzs) == 2: # two by two
        nz = nzs[0]
        nob1 = CTK.Nb[nz]+1
        noz1 = CTK.Nz[nz]
        z1 = CTK.t[2][nob1][2][noz1]
        nz = nzs[1]
        nob2 = CTK.Nb[nz]+1
        noz2 = CTK.Nz[nz]
        z2 = CTK.t[2][nob2][2][noz2]
        zo1 = XOR.conformUnstr(z1, z2, tol=tol, itermax=1)
        zo2 = XOR.conformUnstr(z2, z1, tol=tol, itermax=1)
        CTK.replace(CTK.t, nob1, noz1, zo1)
        CTK.replace(CTK.t, nob2, noz2, zo2)
    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            z = XOR.conformUnstr(z, tol=tol, itermax=1)
            CTK.replace(CTK.t, nob, noz, z)
    
    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Surface conformized.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFixer2', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Fix holes in surfaces.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkFixer2')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- mode pour gapsmanager -
    V = TK.StringVar(win); V.set('Nodes'); VARS.append(V)
    
    # -1- tol pour conformUnstr -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)

    # - Slider -
    B = TTK.Scale(Frame, from_=-50, to=50, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, value=0)
    WIDGETS['bump'] = B
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Bump.')
    
    # - Fix gap in contour -
    B = TTK.Button(Frame, text="Fix gap in contour", command=fixGap)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fix gap in a contour.')

    # - Fix gaps in a set of patches (with overlap) -
    B = TTK.Button(Frame, text="Fix gaps in patches", command=fixGaps)
    B.grid(row=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fix gaps in patches (even with overlap).')
    B = TTK.OptionMenu(Frame, VARS[0], 'Nodes', 'Centers', 'Unknown', 'Slice')
    B.grid(row=1, column=1, sticky=TK.EW)
    
    # - conformUnstr -
    B = TTK.Button(Frame, text="conformUnstr", command=conformUnstr)
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Conformize a TRI surface.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Tolerance used in conformUnstr operations.\n0. means automatic setting.')

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
    (win, menu, file, tools) = CTK.minimal('tkFixer '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
