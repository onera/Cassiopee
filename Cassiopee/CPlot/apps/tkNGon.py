# - tkNGon -
"""Perform NGON operations."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Transform.PyTree as T
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Converti des zones en NGon
# IN: t, cplot.selectedZones
# OUT: t avec zones remplacees et affiche
#==============================================================================
def convert2NGon():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            a = C.convertArray2NGon(z)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = False; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Zones converted to NGon.\n')
    else:
        Panels.displayErrors(errors, header='Error: convert2NGon')
        CTK.TXT.insert('START', 'NGon conversion fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Break elements
# IN: t, cplot.selectedZones
# OUT: t avec zones remplacees et affiche
#==============================================================================
def breakElts():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            zones = T.breakElements(z)
            if (len(zones) > 0): CTK.replace(CTK.t, nob, noz, zones[0])
            for zz in zones[1:]: CTK.add(CTK.t, nob, -1, zz)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Zones converted to basic elements.\n')
    else:
        Panels.displayErrors(errors, header='Error: breakElts')
        CTK.TXT.insert('START', 'Break elts fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Calcul le dual des zones
# IN: t, cplot.selectedZones
# OUT: t avec zones remplacees et affiche
#==============================================================================
def dual():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            a = T.dual(z)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Dual of zones computed.\n')
    else:
        Panels.displayErrors(errors, header='Error: dual')
        CTK.TXT.insert('START', 'Dual fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Conformise un NGON topologique (hanging nodes)
# IN: t, cplot.selectedZones
# OUT: t avec zones remplacees et affiche
#==============================================================================
def conformize():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            zp = C.conformizeNGon(z)
            CTK.replace(CTK.t, nob, noz, zp)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Zones conformized.\n')
    else:
        Panels.displayErrors(errors, header='Error: conformize')
        CTK.TXT.insert('START', 'Conformize fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Converti des NGons v3 en v4
# IN: t, cplot.selectedZones
# OUT: t avec zones remplacees et affiche
#==============================================================================
def adaptNGon32NGon4():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            a = Internal.adaptNGon32NGon4(z)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = False; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Zones converted to NGon v4.\n')
    else:
        Panels.displayErrors(errors, header='Error: adaptNGon32NGon4')
        CTK.TXT.insert('START', 'NGon conversion fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Converti des NGons v4 en v3
# IN: t, cplot.selectedZones
# OUT: t avec zones remplacees et affiche
#==============================================================================
def adaptNGon42NGon3():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        try:
            a = Internal.adaptNGon42NGon3(z)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = False; errors += [0,str(e)]

    if not fail:
        CTK.TXT.insert('START', 'Zones converted to NGon v3.\n')
    else:
        Panels.displayErrors(errors, header='Error: adaptNGon42NGon3')
        CTK.TXT.insert('START', 'NGon conversion fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkNGon  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage Polyhedral (NGON) meshes.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkNGon')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -

    # - Convert2NGon -
    B = TTK.Button(Frame, text="Convert2NGon", command=convert2NGon)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Convert selection to NGon.')

    # - Break Elements -
    B = TTK.Button(Frame, text="BreakElts", command=breakElts)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Convert a NGon zone to basic elements.')

    # - Dual -
    B = TTK.Button(Frame, text="Dual", command=dual)
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Compute the dual mesh of a NGon.')

    # - Conformize -
    B = TTK.Button(Frame, text="Conformize", command=conformize)
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Conformize a NGon (split faces to match hanging nodes).')

    # - adaptNGONv3 to NGONv4
    B = TTK.Button(Frame, text="NGonv3->v4", command=adaptNGon32NGon4)
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Adapt a NGon with CGNSv3 format to NGon with CGNSv4 format.')

    # - adaptNGONv4 to NGONv3
    B = TTK.Button(Frame, text="NGonv4->v3", command=adaptNGon42NGon3)
    B.grid(row=2, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Adapt a NGon with CGNSv4 format to NGon with CGNSv3 format.')

#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['BlockNoteBook'].add(WIDGETS['frame'], text='tkNGon')
    except: pass
    CTK.WIDGETS['BlockNoteBook'].select(WIDGETS['frame'])

#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['BlockNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkNGon '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
