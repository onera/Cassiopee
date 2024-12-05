# - tkFixer -
"""Fix holes in mesh."""
try: import tkinter as TK
except: import Tkinter as TK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Generator.PyTree as G
import Transform.PyTree as T

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setContour():
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
def setSurface():
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
def fixGap():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Contours
    name = VARS[0].get()
    names = name.split(';')
    contours = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if (bases != []):
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): contours.append(z)
    # Surfaces
    name = VARS[1].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if (bases != []):
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if (z[0] == sname[1]): surfaces.append(z)

    p = G.plaster(contours, surfaces)
    contours = C.convertArray2Tetra(contours)
    contours = T.join(contours)
    contours = G.close(contours)
    b = G.gapfixer(contours, p)
    CTK.saveTree()
    CTK.t[2][1][2].append(b)
    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Gap fixed.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                          text='tkFixe  [ + ]  r', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Fix holes in surfaces.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - VARS -
    # -0- Contour to fix -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- Surface of contour to fix
    V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Contour -
    B = TK.Button(Frame, text="Contour", command=setContour)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set contour to fix.')
    B = TK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=0, column=1, sticky=TK.EW)

    # - Surface -
    B = TK.Button(Frame, text="Surface", command=setSurface)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set surface to fix.')
    B = TK.Entry(Frame, textvariable=VARS[1], background='White')
    B.grid(row=1, column=1, sticky=TK.EW)

    # - Fix -
    B = TK.Button(Frame, text="Fix gap", command=fixGap)
    B.grid(row=2, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fix gap.')

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
    (win, menu, file, tools) = CTK.minimal('tkFixer '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
