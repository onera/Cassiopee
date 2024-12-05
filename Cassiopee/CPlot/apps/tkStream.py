# - tkStream -
"""Compute streamlines."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Post.PyTree as P

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def updateVarNameList__(no):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        vars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz])
    m = WIDGETS['variable'+str(no)].children['menu']
    m.delete(0, TK.END)
    allvars = []
    if len(vars) > 0:
        for v in vars[0]: allvars.append(v)
    for i in allvars:
        m.add_command(label=i, command=lambda v=VARS[no], l=i:v.set(l))

def updateVarNameList2__(no):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        vars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz])

    allvars = []
    if len(vars) > 0:
        for v in vars[0]: allvars.append(v)
    if 'variable'+str(no) in WIDGETS:
        WIDGETS['variable'+str(no)]['values'] = allvars

#==============================================================================
def updateVarNameList1(event=None):
    updateVarNameList__(1)
def updateVarNameList2(event=None):
    updateVarNameList__(2)
def updateVarNameList3(event=None):
    updateVarNameList__(3)
def updateVarNameList1_2(event=None):
    updateVarNameList2__(1)
def updateVarNameList2_2(event=None):
    updateVarNameList2__(2)
def updateVarNameList3_2(event=None):
    updateVarNameList2__(3)

#==============================================================================
def streamSurface():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    npts = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(npts) != 1:
        CTK.TXT.insert('START', 'Number of points in stream incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    npts = npts[0]
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    v1 = VARS[1].get(); v2 = VARS[2].get(); v3 = VARS[3].get()
    streams = []
    fail = False; errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = C.convertArray2Tetra(z)
        try:
            stream = P.streamSurf(CTK.t, z, [v1, v2, v3], N=npts)
            streams.append(stream)
        except Exception as e:
            fail = True; errors += [0,str(e)]
    CTK.setCursor(2, WIDGETS['streamSurface'])
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'STREAMS', 2)
    b = Internal.getNodesFromName1(CTK.t, 'STREAMS')
    nob = C.getNobOfBase(b[0], CTK.t)
    for i in streams: CTK.add(CTK.t, nob, -1, i)
    if not fail:
        CTK.TXT.insert('START', 'Stream surface created.\n')
    else:
        Panels.displayErrors(errors, header='Error: streamSurf')
        CTK.TXT.insert('START', 'Sream surface fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['streamSurface'])


#==============================================================================
def streamLine():
    if CTK.t == []: return
    npts = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(npts) != 1:
        CTK.TXT.insert('START', 'Number of points in stream incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    npts = npts[0]
    v1 = VARS[1].get(); v2 = VARS[2].get(); v3 = VARS[3].get()

    l = CPlot.getActivePoint()
    if l == []:
        CTK.TXT.insert('START', 'Click to select starting point...\n')
        return
    CTK.setCursor(2, WIDGETS['streamLine'])
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'STREAMS', 2)
    b = Internal.getNodesFromName1(CTK.t, 'STREAMS')
    nob = C.getNobOfBase(b[0], CTK.t)

    # Arbre source (on enleve les Bases CANVAS, CONTOURS et STREAMS)
    print((l[0], l[1], l[2]))
    source = C.newPyTree()
    bases = Internal.getBases(CTK.t)
    for b in bases:
        if b[0] != 'CANVAS' and b[0] != 'CONTOURS' and b[0] != 'STREAMS':
            source[2].append(b)
    try:
        #stream = P.streamLine(source, (l[0], l[1], l[2]), [v1, v2, v3], N=npts)
        #CTK.add(CTK.t, nob, -1, stream)
        stream = P.streamLine2(source, (l[0], l[1], l[2]), [v1, v2, v3], N=npts)
        for s in stream: CTK.add(CTK.t, nob, -1, s)
        CTK.TXT.insert('START', 'Stream line created.\n')
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: streamLine')
        CTK.TXT.insert('START', 'Stream line fails.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['streamLine'])

#==============================================================================
def streamRibbon():
    if CTK.t == []: return
    npts = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(npts) != 1:
        CTK.TXT.insert('START', 'Number of points in stream incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error') ; return
    npts = npts[0]
    v1 = VARS[1].get(); v2 = VARS[2].get(); v3 = VARS[3].get()

    l = CPlot.getActivePoint()
    if l == []:
        CTK.TXT.insert('START', 'Click to select starting point...\n')
        return

    CTK.setCursor(2, WIDGETS['streamRibbon'])
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'STREAMS', 2)
    b = Internal.getNodesFromName1(CTK.t, 'STREAMS')
    nob = C.getNobOfBase(b[0], CTK.t)

    # Arbre source (on enleve les Bases CANVAS, CONTOURS et STREAMS)
    source = C.newPyTree()
    bases = Internal.getBases(CTK.t)
    for b in bases:
        if b[0] != 'CANVAS' and b[0] != 'CONTOURS' and b[0] != 'STREAMS':
            source[2].append(b)

    try:
        stream = P.streamRibbon(source, (l[0], l[1], l[2]), (0,0,0.01),
                                [v1, v2, v3], N=npts)
        CTK.add(CTK.t, nob, -1, stream)
        CTK.TXT.insert('START', 'Stream ribbon created.\n')
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: streamRibbon')
        CTK.TXT.insert('START', 'Stream ribbon fails.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    CTK.setCursor(0, WIDGETS['streamRibbon'])

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkStream  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Compute stream lines/surfaces.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkStream')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- nptsmax -
    V = TK.StringVar(win); V.set('2000'); VARS.append(V)
    if 'tkStreamNpts' in CTK.PREFS:
        V.set(CTK.PREFS['tkStreamNpts'])
    # -1- Var0 for vector -
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -2- Var1 for vector -
    V = TK.StringVar(win); V.set('CoordinateY'); VARS.append(V)
    # -3- Var2 for vector -
    V = TK.StringVar(win); V.set('CoordinateZ'); VARS.append(V)

    # - Menu des variables -
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[1], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList1)
        F.grid(row=0, column=0, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 1.')
        WIDGETS['variable1'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[1],
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList1_2)
        F.grid(row=0, column=0, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 1.')
        WIDGETS['variable1'] = B

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[2], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=0, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 2.')
        WIDGETS['variable2'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[2],
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2_2)
        F.grid(row=0, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 2.')
        WIDGETS['variable2'] = B

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[3], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList3)
        F.grid(row=0, column=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 3.')
        WIDGETS['variable3'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[3],
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList3_2)
        F.grid(row=0, column=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 3.')
        WIDGETS['variable3'] = B

    # - nptsmax -
    B = TTK.Label(Frame, text="nptsmax")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Max number of points of streams.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='Max number of points of streams.')
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)

    # - Stream line -
    B = TTK.Button(Frame, text="Line", command=streamLine)
    WIDGETS['streamLine'] = B
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Draw a stream line.')

    # - Stream ribbon -
    B = TTK.Button(Frame, text="Ribbon", command=streamRibbon)
    WIDGETS['streamRibbon'] = B
    B.grid(row=2, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Draw a stream ribbon.')

    # - Stream surface -
    B = TTK.Button(Frame, text="Surface", command=streamSurface)
    WIDGETS['streamSurface'] = B
    B.grid(row=2, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Draw a stream surface from a BAR.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['PostNoteBook'].add(WIDGETS['frame'], text='tkStream')
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
def saveApp():
    CTK.PREFS['tkStreamNpts'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('2000')
    CTK.PREFS['tkStreamNpts'] = VARS[0].get()
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
    (win, menu, file, tools) = CTK.minimal('tkStream '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
