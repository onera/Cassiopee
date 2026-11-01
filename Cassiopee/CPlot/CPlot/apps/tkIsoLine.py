# - tkIsoLine -
"""Draw isolines."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Converter.Internal as Internal
import Post.PyTree as P
import Transform.PyTree as T

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def updateVarNameList(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        vars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz])
    m = WIDGETS['field'].children['menu']
    m.delete(0, TK.END)
    if len(vars) == 0: return
    for i in vars[0]:
        m.add_command(label=i, command=lambda v=VARS[0],l=i:v.set(l))

#==============================================================================
def updateVarNameList2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        vars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz])
    if len(vars) == 0: return
    if 'field' in WIDGETS:
        WIDGETS['field']['values'] = vars[0]

#==============================================================================
def drawIsoLines():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    field = VARS[0].get()
    nlevels = VARS[1].get()

    try: nlevels = int(nlevels)
    except: levels = 30
    if nlevels < 2: nlevels = 2

    nzs = CPlot.getSelectedZones()
    CTK.saveTree()

    fmin = CTK.varsFromWidget(VARS[3].get(), type=1)
    if fmin == []: fmin = C.getMinValue(CTK.t, field)
    else: fmin = fmin[0]
    fmax = CTK.varsFromWidget(VARS[4].get(), type=1)
    if fmax == []: fmax = C.getMaxValue(CTK.t, field)
    else: fmax = fmax[0]
    if nzs == []:
        z = Internal.getZones(CTK.t)
    else:
        z = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z.append(CTK.t[2][nob][2][noz])

    isos = []
    nlevels += 1 # pour etre coeherent avec les niveaux d'iso solides
    fail = False; errors = []
    for v in range(nlevels):
        value = fmin + (fmax-fmin)/(nlevels-1)*v
        for zone in z:
            try:
                i = P.isoLine(zone, field, value)
                isos.append(i)
            except Exception as e:
                fail = True; errors += [0,str(e)]

    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    bases = Internal.getNodesFromName1(CTK.t, 'CONTOURS')
    nob = C.getNobOfBase(bases[0], CTK.t)
    if isos != []:
        isos = T.join(isos)
        CTK.add(CTK.t, nob, -1, isos)
    if not fail:
        CTK.TXT.insert('START', 'Isolines extracted.\n')
    else:
        Panels.displayErrors(errors, header='Error: Isolines')
        CTK.TXT.insert('START', 'Isolines fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def extractIsoLine():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    field = VARS[0].get()
    value = VARS[2].get()

    try: value = float(value)
    except: value = 1.

    nzs = CPlot.getSelectedZones()
    CTK.saveTree()

    if nzs == []:
        z = Internal.getZones(CTK.t)
    else:
        z = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z.append(CTK.t[2][nob][2][noz])

    isos = []
    fail = False; errors = []
    for zone in z:
        try:
            i = P.isoLine(zone, field, value)
            isos.append(i)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    CTK.t = C.addBase2PyTree(CTK.t, 'CONTOURS', 1)
    bases = Internal.getNodesFromName1(CTK.t, 'CONTOURS')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in isos: CTK.add(CTK.t, nob, -1, i)
    if not fail:
        CTK.TXT.insert('START', 'Isolines extracted.\n')
    else:
        Panels.displayErrors(errors, header='Error: Isolines')
        CTK.TXT.insert('START', 'Isolines fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkIsoLine  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Compute isolines.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkIsoLine')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Field name -
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -1- nlevels -
    V = TK.StringVar(win); V.set('25'); VARS.append(V)
    if 'tkIsoLineLevels' in CTK.PREFS:
        V.set(CTK.PREFS['tkIsoLineLevels'])
    # -2- value -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkIsoLineValue' in CTK.PREFS:
        V.set(CTK.PREFS['tkIsoLineValue'])
    # -3- min iso
    V = TK.StringVar(win); V.set('MIN'); VARS.append(V)
    if 'tkIsoLineMin' in CTK.PREFS:
        V.set(CTK.PREFS['tkIsoLineMin'])
    # -4- max iso
    V = TK.StringVar(win); V.set('MAX'); VARS.append(V)
    V = TK.StringVar(win); V.set('MIN'); VARS.append(V)
    if 'tkIsoLineMax' in CTK.PREFS:
        V.set(CTK.PREFS['tkIsoLineMax'])

    # - field name -
    B = TTK.Label(Frame, text="Field:")
    B.grid(row=0, column=0, sticky=TK.EW)
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[0], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList)
        F.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Extracted field.')
        WIDGETS['field'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[0],
                         values=[], state='readonly')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Extracted field.')
        WIDGETS['field'] = B

    if CTK.t != []:
        vars = C.getVarNames(CTK.t)
        if len(vars) > 0:
            if (len(vars[0])>0): VARS[0].set(vars[0][0])

    # - nlevels -
    #B = TK.Label(Frame, text="Nlevels:")
    #B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=7)
    BB = CTK.infoBulle(parent=B, text='Number of levels.')
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=7)
    BB = CTK.infoBulle(parent=B, text='Min value.')
    B.grid(row=1, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=7)
    BB = CTK.infoBulle(parent=B, text='Max value.')
    B.grid(row=1, column=2, sticky=TK.EW)

    # - Draw all isolines -
    B = TTK.Button(Frame, text="Extract isolines", command=drawIsoLines)
    BB = CTK.infoBulle(parent=B, text='Extract isolines.')
    B.grid(row=2, column=0, columnspan=3, sticky=TK.EW)

    # - Value -
    B = TTK.Label(Frame, text="Value:")
    B.grid(row=3, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White')
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)

    # - Extract one isoline -
    B = TTK.Button(Frame, text="Extract isoline", command=extractIsoLine)
    B.grid(row=4, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract one isoline to CONTOURS.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['PostNoteBook'].add(WIDGETS['frame'], text='tkIsoLine')
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
    CTK.PREFS['tkIsoLineLevels'] = VARS[1].get()
    CTK.PREFS['tkIsoLineValue'] = VARS[2].get()
    CTK.PREFS['tkIsoLineMin'] = VARS[3].get()
    CTK.PREFS['tkIsoLineMax'] = VARS[4].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[1].set('25')
    VARS[2].set('1.')
    VARS[3].set('MIN')
    VARS[4].set('MAX')
    CTK.PREFS['tkIsoLineLevels'] = VARS[1].get()
    CTK.PREFS['tkIsoLineValue'] = VARS[2].get()
    CTK.PREFS['tkIsoLineMin'] = VARS[3].get()
    CTK.PREFS['tkIsoLineMax'] = VARS[4].get()
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
    (win, menu, file, tools) = CTK.minimal('tkIsoLine '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
