# - plot 1d data (plots) -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Geom.PyTree as D
import Post.PyTree as P
import Transform.PyTree as T

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def updateVar1NameList(event=None):
    if CTK.t == []: return
    vars = C.getVarNames(CTK.t)
    vars[0] += ['s']
    m = WIDGETS['var1'].children['menu']
    m.delete(0, TK.END)
    if len(vars) == 0: return
    for i in vars[0]:
        m.add_command(label=i, command=lambda v=VARS[3],l=i:v.set(l))

def updateVar1NameList2(event=None):
    if CTK.t == []: return
    vars = C.getVarNames(CTK.t)
    vars[0] += ['s']
    if len(vars) == 0: return
    if 'var1' in WIDGETS:
        WIDGETS['var1']['values'] = vars[0]

#==============================================================================
def updateVar2NameList(event=None):
    if CTK.t == []: return
    vars = C.getVarNames(CTK.t)
    vars[0] += ['s']
    m = WIDGETS['var2'].children['menu']
    m.delete(0, TK.END)
    if len(vars) == 0: return
    for i in vars[0]:
        m.add_command(label=i, command=lambda v=VARS[4],l=i:v.set(l))

def updateVar2NameList2(event=None):
    if CTK.t == []: return
    vars = C.getVarNames(CTK.t)
    vars[0] += ['s']
    if len(vars) == 0: return
    if 'var2' in WIDGETS:
        WIDGETS['var2']['values'] = vars[0]

#==============================================================================
# Fonction display 1D
#==============================================================================
def display1D(event=None):
    if CTK.t == []: return

    # Get slot
    try: slot = int(VARS[5].get())
    except: slot = 0
    # Get grid size
    try:
        gridSize = VARS[1].get()
        grids = gridSize.split(';')
        if (len(grids) == 1): gridSize = (int(grids[0]),1)
        else: gridSize = (int(grids[0]), int(grids[1]))
    except: gridSize = (1,1)
    CPlot.setState(gridSize=gridSize)
    # Get grid pos
    try:
        gridPos = VARS[2].get()
        grids = gridPos.split(';')
        if len(grids) == 1: gridPos = (int(grids[0]),1)
        else: gridPos = (int(grids[0]), int(grids[1]))
    except: gridPos = (0,0)

    # Recupere la direction pour la coupe ou 'Elements'
    dir = VARS[0].get()
    if dir == 'None': CPlot.display1D([], slot=slot); return # clear

    # Recupere le pt pour la coupe ou les elements 1D
    if dir == 'Elements': # elements -> recupere les elements
        if CTK.__MAINTREE__ <= 0:
            CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        nzs = CPlot.getSelectedZones()
        if nzs == []:
            CTK.TXT.insert('START', 'Selection is empty.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        points = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            selected = CTK.t[2][nob][0]+'/'+z[0]
            points.append(selected)
    elif dir == 'I' or dir == 'J' or dir == 'K': # indice -> recupere les indices + la zone
        if CTK.__MAINTREE__ <= 0:
            CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        nz = CPlot.getSelectedZone()
        if nz == -1:
            CTK.TXT.insert('START', 'Selection is empty.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        points = []
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        selected = CTK.t[2][nob][0]+'/'+z[0]
        index = CPlot.getActivePointIndex()
        points = (selected,index)
    else: # les coupes -> recupere les coord du pt 
        point = CPlot.getActivePoint()
        if point == []: point = (0.,0.,0.)

    # Recupere les variables a afficher
    var1 = VARS[3].get()
    var1 = var1.replace('centers:', '')
    var2 = VARS[4].get()
    var2 = var2.replace('centers:', '')

    # Recupere les zones actives
    actives = []
    zones = Internal.getZones(CTK.t)
    if CTK.__MAINTREE__ == 1: 
        nzs = CPlot.getActiveZones()
        for nz in nzs: actives.append(zones[nz]) 
    else: actives = zones
    if actives == []: return

    if (dir == 'X (Y)'):
        elts = P.isoSurfMC(actives, 'CoordinateY', point[1])
        if elts != []:
            elts2 = P.isoSurfMC(elts, 'CoordinateZ', point[2])
            if (elts2 != []): elts = elts2
    elif (dir == 'Y (X)'):
        elts = P.isoSurfMC(actives, 'CoordinateX', point[0])
        if elts != []:
            elts2 = P.isoSurfMC(elts, 'CoordinateZ', point[2])
            if (elts2 != []): elts = elts2
    elif (dir == 'Z (X)'):
        elts = P.isoSurfMC(actives, 'CoordinateX', point[0])
        if (elts != []):
            elts2 = P.isoSurfMC(elts, 'CoordinateY', point[1])
            if (elts2 != []): elts = elts2
    elif (dir == 'X (Z)'):
        elts = P.isoSurfMC(actives, 'CoordinateZ', point[2])
        if elts != []:
            elts2 = P.isoSurfMC(elts, 'CoordinateY', point[1])
            if (elts2 != []): elts = elts2
    elif (dir == 'Y (Z)'):
        elts = P.isoSurfMC(actives, 'CoordinateZ', point[2])
        if elts != []:
            elts2 = P.isoSurfMC(elts, 'CoordinateX', point[0])
            if (elts2 != []): elts = elts2
    elif (dir == 'Z (Y)'):
        elts = P.isoSurfMC(actives, 'CoordinateY', point[1])
        if (elts != []):
            elts2 = P.isoSurfMC(elts, 'CoordinateX', point[0])
            if (elts2 != []): elts = elts2
    elif (dir == 'I'):
        v = points[0]; ind = points[1]
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        elts = []
        if bases != []:
            zones = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in zones:
                if (z[0] == sname[1]):
                    try:
                        zp = C.center2Node(z, Internal.__FlowSolutionCenters__)
                        zp = T.subzone(zp, (1,ind[3],ind[4]), (-1,ind[3],ind[4]))
                        elts.append(zp)
                    except: pass
    elif (dir == 'J'):
        v = points[0]; ind = points[1]
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        elts = []
        if bases != []:
            zones = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in zones:
                if (z[0] == sname[1]):
                    try:
                        zp = C.center2Node(z, Internal.__FlowSolutionCenters__)
                        zp = T.subzone(zp, (ind[2],1,ind[4]), (ind[2],-1,ind[4]))
                        elts.append(zp)
                    except: pass
    elif (dir == 'K'):
        v = points[0]; ind = points[1]
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        elts = []
        if bases != []:
            zones = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in zones:
                if (z[0] == sname[1]):
                    try:
                        zp = C.center2Node(z, Internal.__FlowSolutionCenters__)
                        zp = T.subzone(zp, (ind[2],ind[3],1), (ind[2],ind[3],-1))
                        elts.append(zp)
                    except: pass
    elif (dir == 'Elements'):
        elts = []
        for v in points:
            v = v.lstrip(); v = v.rstrip()
            sname = v.split('/', 1)
            bases = Internal.getNodesFromName1(CTK.t, sname[0])
            if (bases != []):
                zones = Internal.getNodesFromType1(bases[0], 'Zone_t')
                for z in zones:
                    if (z[0] == sname[1]): elts.append(z)
    if elts == []:
        CTK.TXT.insert('START', 'Nothing to display.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    try: elts = D.getCurvilinearAbscissa(elts)
    except: pass

    # Fit first axis
    pos = WIDGETS['rangePos'].get() / 50.-1.
    zoom = WIDGETS['rangeZoom'].get() / 120.
    minv1 = C.getMinValue(elts, var1)
    maxv1 = C.getMaxValue(elts, var1)
    if (maxv1-minv1 < 1.e-6): maxv1 += 5.e-7; minv1 -= 5.e-7

    # active point localisation
    nz = CPlot.getSelectedZone()
    if (nz != -1):
        ind = CPlot.getActivePointIndex()
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        f1 = C.getValue(z, var1, ind[0])
        try:
            r1min = (f1-minv1)*zoom+minv1 + pos*(1.-zoom)*(maxv1-minv1)
            r1max = (f1-maxv1)*zoom+maxv1 + pos*(1.-zoom)*(maxv1-minv1)
        except: # var1 not found in z, le cherche dans elts
            xf1 = C.getValue(z, 'CoordinateX', ind[0])
            yf1 = C.getValue(z, 'CoordinateY', ind[0])
            zf1 = C.getValue(z, 'CoordinateZ', ind[0])
            f1 = minv1+0.5*(maxv1-minv1)
            r1min = 0.5*(maxv1-minv1)*zoom+minv1 +pos*(1.-zoom)*(maxv1-minv1)
            r1max = -0.5*(maxv1-minv1)*zoom+maxv1 +pos*(1.-zoom)*(maxv1-minv1)
    else:
        f1 = minv1+0.5*(maxv1-minv1)
        r1min = 0.5*(maxv1-minv1)*zoom+minv1 +pos*(1.-zoom)*(maxv1-minv1)
        r1max = -0.5*(maxv1-minv1)*zoom+maxv1 +pos*(1.-zoom)*(maxv1-minv1)

    # Fit second axis
    p = P.selectCells(elts, '({%s} < %20.16g) & ({%s} > %20.16g)'%(var1,r1max,var1,r1min))
    minv2 = C.getMinValue(p, var2)
    maxv2 = C.getMaxValue(p, var2)

    # display
    CPlot.display1D(p, slot=slot, bgBlend=0., gridPos=gridPos, 
                    var1=var1, var2=var2, 
                    r1=(r1min,r1max), r2=(minv2,maxv2))
    CTK.TXT.insert('START', 'Plot displayed.\n')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkPlot  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='My personal applet.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkPlot')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Direction -
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    if 'tkPlotDirection' in CTK.PREFS: V.set(CTK.PREFS['tkPlotDirection'])
    # -1- grid size -
    V = TK.StringVar(win); V.set('3;3'); VARS.append(V)
    if 'tkPlotGridSize' in CTK.PREFS: V.set(CTK.PREFS['tkPlotGridSize'])
    # -2- grid pos -
    V = TK.StringVar(win); V.set('0;0'); VARS.append(V)
    if 'tkPlotGridPos' in CTK.PREFS: V.set(CTK.PREFS['tkPlotGridPos'])
    # -3- Var1
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -4- Var2
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -5- Slot
    V = TK.StringVar(win); V.set('0'); VARS.append(V)

    # - Slot -
    B = TTK.Entry(Frame, textvariable=VARS[5], width=2)
    BB = CTK.infoBulle(parent=B, text='Slot.')
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)

    # - grid size -
    B = TTK.Entry(Frame, textvariable=VARS[1], width=2)
    BB = CTK.infoBulle(parent=B, text='Grid size (1;1).')
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)

    # - grid pos -
    B = TTK.Entry(Frame, textvariable=VARS[2], width=2)
    BB = CTK.infoBulle(parent=B, text='Position of slot in grid (0;0).')
    B.grid(row=0, column=2, columnspan=1, sticky=TK.EW)

    # - Element 1D -
    B = TTK.OptionMenu(Frame, VARS[0], 'None', 'X (Y)', 'Y (X)', 'Z (X)', 'X (Z)', 'Y (Z)', 'Z (Y)', 'I', 'J', 'K', 'Elements')
    B.grid(row=1, column=0, columnspan=3, sticky=TK.EW)

    # Var1
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[3], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVar1NameList)
        F.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 1 (abcsiss).')
        WIDGETS['var1'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[3], 
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVar1NameList2)
        F.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 1 (absciss).')
        WIDGETS['var1'] = B

    # Var2
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[4], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVar2NameList)
        F.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 2 (ordinates).')
        WIDGETS['var2'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[4], 
                         values=[], state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVar2NameList2)
        F.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable 2 (ordinates).')
        WIDGETS['var2'] = B

    # - Ranges -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=display1D, value=50)
    WIDGETS['rangePos'] = B
    B.grid(row=3, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Move range.')
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=display1D, value=0)
    WIDGETS['rangeZoom'] = B
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Zoom range.')

    # - Set in slot -
    B = TTK.Button(Frame, text="Set", command=display1D)
    B.grid(row=4, column=0, columnspan=3, sticky=TK.EW)

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
def saveApp():
    CTK.PREFS['tkPlotDirection'] = VARS[0].get()
    CTK.PREFS['tkPlotGridSize'] = VARS[1].get()
    CTK.PREFS['tkPlotGridPos'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('None')
    VARS[1].set('1;1')
    VARS[2].set('0;0')
    CTK.PREFS['tkPlotDirection'] = VARS[0].get()
    CTK.PREFS['tkPlotGridSize'] = VARS[1].get()
    CTK.PREFS['tkPlotGridPos'] = VARS[2].get()
    CTK.savePrefFile()

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
    (win, menu, file, tools) = CTK.minimal('tkPlot '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
