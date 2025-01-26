# - Extension of a surface starting from an edge of the surface following normals -
try: import tkinter as TK
except: import Tkinter as TK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Generator.PyTree as G

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def extendSurf():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # surfaces
    name = VARS[0].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: surfaces.append(z)

    # - Hauteur de chaque maille -
    dhloc = CTK.varsFromWidget(VARS[1].get(), type=1); dhloc = dhloc[0]
    N = CTK.varsFromWidget(VARS[2].get(), type=2); N = N[0]
    dh = G.cart((0.,0.,0.), (dhloc,1.,1.), (N+1,1,1))

    # - nb d'iterations de lissage -
    nit = CTK.varsFromWidget(VARS[3].get(), type=2); nit = nit[0]

    # - contour -
    nzs = CPlot.getSelectedZones()
    if (nzs == []):
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    contours = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        contour = C.convertBAR2Struct(z)
        contours.append(contour)

    # Extension
    zlist = G.buildExtension(contours, surfaces, dh, niter=nit)
    c = 0
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        CTK.t[2][nob][2][noz] = zlist[c]; c += 1

    CTK.TXT.insert('START', 'Surface extension done.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def setSurface():
    if (CTK.t == []): return
    if (CTK.__MAINTREE__ <= 0):
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if (nzs == []):
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
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                          text='tkExtension  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Surface extension.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - VARS -
    # -0- Surface to extend -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- Hauteur de chaque maille -
    V = TK.StringVar(win); V.set('1.e-1'); VARS.append(V)
    # -2- Nombre de layers a ajouter
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -3- Nombre d'iterations de lissage
    V = TK.StringVar(win); V.set('0'); VARS.append(V)

    # - Surface -
    B = TK.Button(Frame, text="Ref surf", command=setSurface)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reference surface.')

    B = TK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)

    # - Extension data
    B = TK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Height of each layer.')

    B = TK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of layers.')

    B = TK.Entry(Frame, textvariable=VARS[3], background='White', width=5)
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of smoothing iterations.')

    # Extend surface starting from contour
    B = TK.Button(Frame, text="Extend", command=extendSurf)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extend.')

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
    if (len(sys.argv) == 2):
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkExtendSurf '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
