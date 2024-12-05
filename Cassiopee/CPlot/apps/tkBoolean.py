# - tkBoolean -
"""Boolean operators (union, intersection, ...)."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Intersector.PyTree as XOR

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# IN: CPlot.getSelectedZones
#==============================================================================
def union():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) < 2:
        CTK.TXT.insert('START', 'Please, select two or more surfaces.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    tol = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(tol) != 1:
        CTK.TXT.insert('START', 'Tolerance is incorrect.\n') 
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = tol[0]

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['union'])
    zlist = []
    deletedZoneNames = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        deletedZoneNames.append(CTK.t[2][nob][0]+Internal.SEP1+CTK.t[2][nob][2][noz][0])
        z = CTK.t[2][nob][2][noz]
        zlist.append(z)

    try: 
        j = XOR.booleanUnion(zlist[0], zlist[1], tol=tol)
    except Exception as e:
        CTK.setCursor(0, WIDGETS['union'])
        Panels.displayErrors([0,str(e)], header='Error: union')
        CTK.TXT.insert('START', 'Union failed\n'); return

    for nz in range(len(zlist)-2):
        try: 
            j = XOR.booleanUnion(j, zlist[nz+2], tol=tol)
        except Exception as e:
            CTK.setCursor(0, WIDGETS['union'])
            Panels.displayErrors([0,str(e)], header='Error: union')
            CTK.TXT.insert('START', 'Union failed.\n'); return

    CTK.setCursor(0, WIDGETS['union'])    
    CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
    CPlot.delete(deletedZoneNames)
    CTK.add(CTK.t, CTK.Nb[0]+1, -1, j)

    CTK.TXT.insert('START', 'Union performed.\n')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def difference():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) != 2:
        CTK.TXT.insert('START', 'Please, select two surfaces.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    tol = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(tol) != 1:
        CTK.TXT.insert('START', 'Tolerance is incorrect.\n') 
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = tol[0]

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['difference'])
    deletedZoneNames = []
    nz = nzs[0]
    nob1 = CTK.Nb[nz]+1
    noz1 = CTK.Nz[nz]
    deletedZoneNames.append(CTK.t[2][nob1][0]+Internal.SEP1+CTK.t[2][nob1][2][noz1][0])
    z1 = CTK.t[2][nob1][2][noz1]
    nz = nzs[1]
    nob2 = CTK.Nb[nz]+1
    noz2 = CTK.Nz[nz]
    deletedZoneNames.append(CTK.t[2][nob2][0]+Internal.SEP1+CTK.t[2][nob2][2][noz2][0])
    z2 = CTK.t[2][nob2][2][noz2]

    try:
        j = XOR.booleanMinus(z1, z2, tol=tol)
        CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
        CPlot.delete(deletedZoneNames)
        CTK.add(CTK.t, nob1, -1, j)
        CTK.TXT.insert('START', 'Difference performed.\n')
    except Exception as e:
        CTK.setCursor(0, WIDGETS['difference'])
        Panels.displayErrors([0,str(e)], header='Error: difference')
        CTK.TXT.insert('START', 'Difference failed.\n')
    #C._fillMissingVariables(CTK.t)
    CTK.setCursor(0, WIDGETS['difference'])

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def difference2():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) != 2:
        CTK.TXT.insert('START', 'Please, select two surfaces.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    tol = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(tol) != 1:
        CTK.TXT.insert('START', 'Tolerance is incorrect.\n') 
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = tol[0]

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['revdiff'])
    deletedZoneNames = []
    nz = nzs[0]
    nob1 = CTK.Nb[nz]+1
    noz1 = CTK.Nz[nz]
    deletedZoneNames.append(CTK.t[2][nob1][0]+Internal.SEP1+CTK.t[2][nob1][2][noz1][0])
    z1 = CTK.t[2][nob1][2][noz1]
    nz = nzs[1]
    nob2 = CTK.Nb[nz]+1
    noz2 = CTK.Nz[nz]
    deletedZoneNames.append(CTK.t[2][nob2][0]+Internal.SEP1+CTK.t[2][nob2][2][noz2][0])
    z2 = CTK.t[2][nob2][2][noz2]

    try:
        j = XOR.booleanMinus(z2, z1, tol=tol)
        CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
        CPlot.delete(deletedZoneNames)
        CTK.add(CTK.t, nob1, -1, j)
        CTK.TXT.insert('START', 'Difference performed.\n')
    except Exception as e:
        CTK.setCursor(0, WIDGETS['revdiff'])
        Panels.displayErrors([0,str(e)], header='Error: difference')
        CTK.TXT.insert('START', 'Difference failed.\n')

    CTK.setCursor(0, WIDGETS['revdiff'])    
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def intersection():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) < 2:
        CTK.TXT.insert('START', 'Please, select two or more surfaces.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    tol = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(tol) != 1:
        CTK.TXT.insert('START', 'Tolerance is incorrect.\n') 
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    tol = tol[0]

    CTK.saveTree()
    CTK.setCursor(2, WIDGETS['intersection'])

    zlist = []
    deletedZoneNames = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        deletedZoneNames.append(CTK.t[2][nob][0]+Internal.SEP1+CTK.t[2][nob][2][noz][0])
        z = CTK.t[2][nob][2][noz]
        zlist.append(z)

    try: j = XOR.booleanIntersection(zlist[0], zlist[1], tol=tol)
    except Exception as e:
        CTK.setCursor(0, WIDGETS['intersection'])
        Panels.displayErrors([0,str(e)], header='Error: intersection')
        CTK.TXT.insert('START', 'Intersection failed.\n'); return

    for nz in range(len(zlist)-2):
        try: j = XOR.booleanIntersection(j, zlist[nz+2], tol=tol)
        except Exception as e:
            CTK.setCursor(0, WIDGETS['intersection'])
            Panels.displayErrors([0,str(e)], header='Error: intersection')
            CTK.TXT.insert('START', 'Intersection failed.\n'); return

    CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
    CPlot.delete(deletedZoneNames)
    CTK.setCursor(0, WIDGETS['intersection'])
    CTK.add(CTK.t, CTK.Nb[0]+1, -1, j)
    CTK.TXT.insert('START', 'Intersection performed.\n')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkBoolean  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Boolean operations\on surfaces.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkBoolean')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- tolerance -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)

    # - Tolerance -
    B = TTK.Label(Frame, text='Tolerance')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Tolerance used in boolean operations.\n0. means automatic setting.')

    # - Union -
    B = TTK.Button(Frame, text="Union", command=union)
    WIDGETS['union'] = B
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Union of two surfaces.')

    # - Difference -
    B = TTK.Button(Frame, text="Difference", command=difference)
    WIDGETS['difference'] = B
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Difference of two surfaces.')

    # - Intersection -
    B = TTK.Button(Frame, text="Intersection", command=intersection)
    WIDGETS['intersection'] = B
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Intersect two surfaces.')

    # - Rev. Difference -
    B = TTK.Button(Frame, text="Rev. Diff", command=difference2)
    WIDGETS['revdiff'] = B
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reversed difference of two surfaces.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkBoolean')
    except: pass
    CTK.WIDGETS['SurfNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['SurfNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkBoolean '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
