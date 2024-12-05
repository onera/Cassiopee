# - tkBody -
"""Creates closed and watertight bodies."""
try: import tkinter as TK
except: import Tkinter as TK
import Converter.Internal as Internal
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Generator.PyTree as G
import Post.PyTree as P

# local widgets list
WIDGETS = []; VARS = []

#==============================================================================
# Extraction des BCWall,BCWallInviscid et BCWallViscous des bases basename
# dans une base BODY#basename
#==============================================================================
def extractBodies():
    pref = 'BODY#'
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    for nz in nzs:
        nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        baseName = CTK.t[2][nob][0]; bodyName = pref+baseName
        CTK.t = C.addBase2PyTree(CTK.t, bodyName, cellDim=2)
        nodes = Internal.getNodesFromName(CTK.t, bodyName)
        p = Internal.getParentOfNode(CTK.t, nodes[0])
        walls = C.extractBCOfType(z, 'BCWall')
        walls += C.extractBCOfType(z, 'BCWallInviscid')
        walls += C.extractBCOfType(z, 'BCWallViscous')
        CTK.t[2][p[1]][2] = CTK.t[2][p[1]][2] + walls

    CTK.TXT.insert('START', 'Walls extracted.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)        

#=========================================================================
def extractExternalContours():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        contours = P.exteriorFacesStructured(z)
        CTK.t[2][nob][2] = CTK.t[2][nob][2] + contours
    CTK.TXT.insert('START', 'External contours extracted.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#=========================================================================
def stitchedHat():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    args = VARS[1].get(); args = args.split(';')
    if len(args) != 2: return
    eps = float(args[0]); eps2 = float(args[1])

    args = VARS[2].get(); args = args.split(';')
    if len(args) != 3: return
    offx = float(args[0]); offy = float(args[1]); offz = float(args[2])

    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            CTK.t[2][nob][2][noz] = G.stitchedHat(CTK.t[2][nob][2][noz],
                                                  (offx,offy,offz), \
                                                  eps, eps2)

    CTK.TXT.insert('START', 'Stitched hat created.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)      
    return

#==============================================================================
def pointedHat():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    args = VARS[3].get(); args = args.split(';')
    if len(args) != 3: return
    x0 = float(args[0]); y0 = float(args[1]); z0 = float(args[2])

    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    else :
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            CTK.t[2][nob][2][noz] = G.pointedHat(CTK.t[2][nob][2][noz],
                                                 (x0,y0,z0))

    CTK.TXT.insert('START', 'Pointed hat created.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)      
    return

#==============================================================================
def closeBody():
    type = VARS[0].get()
    if type == 'stitchedHat': return stitchedHat()
    else: return pointedHat()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE, 
                          text='tkBody  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Close surface meshes.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkBody')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- close bodies -
    V = TK.StringVar(win); V.set('stitchedHat'); VARS.append(V)
    # -1- stitchedHat : eps; eps2
    V = TK.StringVar(win); V.set('1.e-5;1.e-5'); VARS.append(V)
    # -2- stitchedHat : offx,offy,offz
    V = TK.StringVar(win); V.set('0.;0.;0.'); VARS.append(V)
    # -3- pointedHat : x,y,z 
    V = TK.StringVar(win); V.set('0.;0.;0.'); VARS.append(V)

    # - extractBodies -
    B = TK.Button(Frame, text="Extract BCWall", command=extractBodies)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract bodies for zones.')

    # - Select exterior faces -
    B = TK.Button(Frame, text="Extract external contours", \
                  command=extractExternalContours)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract external contours of body.')

    # - closeBodies
    B = TK.Button(Frame, text="Close body", command=closeBody)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Close selected contour with hat.')
    B = TK.OptionMenu(Frame, VARS[0], 'stitchedHat', 'pointedHat')
    B.grid(row=1, column=1, sticky=TK.EW)

    B = TK.Label(Frame, text="eps;eps2")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Tolerances for stitchedHat.')
    B = TK.Entry(Frame, textvariable=VARS[1])
    B.grid(row=2, column=1, sticky=TK.EW)

    B = TK.Label(Frame, text="offx;offy;offz")
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Offsets for stitchedHat.')
    B = TK.Entry(Frame, textvariable=VARS[2])
    B.grid(row=3, column=1, sticky=TK.EW)

    B = TK.Label(Frame, text="x0;y0;z0")
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Coordinates of top of pointed hat.')
    B = TK.Entry(Frame, textvariable=VARS[3])
    B.grid(row=4, column=1, sticky=TK.EW)

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
    (win, menu, file, tools) = CTK.minimal('tkBody '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
