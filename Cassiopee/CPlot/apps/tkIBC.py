# - IBC app -
"""Interface to set IBCs."""

try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Geom.PyTree as D
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Geom.IBM as D_IBM
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

# Set IBM data in zone
# if zone is 1D STRUCT or BAR: remesh
def _setDataInZone(z, snear, ibctype, dfar, inv):
    # set IBM data in .Solver#define
    n = Internal.createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
    Internal.createUniqueChild(n, 'snear', 'DataArray_t', value=snear)
    Internal.createUniqueChild(n, 'ibctype', 'DataArray_t', value=ibctype)
    Internal.createUniqueChild(n, 'dfar', 'DataArray_t', value=dfar)
    Internal.createUniqueChild(n, 'inv', 'DataArray_t', value=inv)
    # Set Extractions triggers in .Solver#define
    if VARS[4].get() == "1": value = 1
    else: value = 0
    Internal.createUniqueChild(n, 'extractWalls', 'DataArray_t', value=value)
    if VARS[5].get() == "1": value = 1
    else: value = 0
    Internal.createUniqueChild(n, 'extractLoads', 'DataArray_t', value=value)

    # remesh surface eventually
    dim = Internal.getZoneDim(z)
    remesh = False
    if dim[0] == 'Structured' and dim[2] == 1 and dim[3] == 1: remesh = True
    if dim[0] == 'Unstructured' and dim[3] == 'BAR': remesh = True
    if remesh: D._uniformize(z, h=float(snear))
    return None


#==============================================================================
# Creates a symmetry plane
#==============================================================================
def symmetrize():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    axis = VARS[7].get()

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    bodySymName=None
    snear_sym = 0

    nob = CTK.Nb[0]+1
    noz = CTK.Nz[0]
    bodySymName = CTK.t[2][nob][0]
    z = CTK.t[2][nob][2][noz]
    snear_sym = Internal.getValue(Internal.getNodeFromName(z,'snear'))

    if axis == 'Around XY-': dir_sym = 1
    elif axis == 'Around XZ-': dir_sym = 2
    elif axis == 'Around YZ-': dir_sym = 3

    D_IBM._symetrizePb(CTK.t, bodySymName, snear_sym, dir_sym=dir_sym)
    CTK.TXT.insert('START', 'Symmetry plane has been created with snear=%f.\n'%snear_sym)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    CPlot.render()

#==============================================================================
# Set data in selected zones
#==============================================================================
def setData():
    if CTK.t == []: return

    snear = VARS[0].get()
    ibctype = VARS[1].get()
    dfar = VARS[2].get()
    if VARS[3].get() == 'out': inv = 0
    else: inv = 1

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        _setDataInZone(z, snear, ibctype, dfar, inv)
        CTK.replace(CTK.t, nob, noz, z)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.TXT.insert('START', 'IBC data set on surface.\n')

#==============================================================================
# View undefined IBC
#==============================================================================
def ViewUndefinedIBC():
    CTK.TXT.insert('START', 'Display undefined IBC zones.\n')
    if CTK.t == []: return

    nzs             = Internal.getNodesFromType2(CTK.t, 'Zone_t')
    VARSlocal       = 10e10
    ZoneLocalString = ''
    ZoneLocal       = []
    bases           = CTK.t[2][1:]
    for b in bases:
        for zone in Internal.getZones(b):
            n = Internal.getNodeFromPath(zone, '.Solver#define/snear')
            if n is not None:
                val = Internal.getValue(n)
                i   = CPlot.getCPlotNumber(CTK.t, b[0], zone[0])
                ZoneLocal.append((i,0))
            else:
                val = -1000
                VARSlocal = min(val,VARSlocal)
                ZoneLocalString += zone[0]+','

    if VARSlocal < 0:
        #VARS[6].set(ZoneLocalString[:-1])
        TTK.setButtonRed(WIDGETS['ViewUndefinedIBC'])

    if VARSlocal > 0:
        #VARS[6].set(ZoneLocalString[:])
        TTK.setButtonGreen(WIDGETS['ViewUndefinedIBC'])

    CPlot.setActiveZones(ZoneLocal)
    WIDGETS['ViewUndefinedIBC'].update()
    CPlot.setState(ghostifyDeactivatedZones=1)

#==============================================================================
# View all defined IBC
#==============================================================================
def ViewAllDefinedIBC(t):
    natives = set()
    zones = Internal.getZones(t)
    for zone in zones:
        n = Internal.getNodeFromPath(zone, '.Solver#define/ibctype')
        if n is not None:
            natives.add(Internal.getValue(n))
    natives = list(natives)
    natives.sort(key=str.lower)
    return natives

#==============================================================================
# Automatic update of all defined IBC
#==============================================================================
# Pour list box
def updateIBCNameList(event=None):
    if CTK.t == []: return

    lb = WIDGETS['IBCLB']
    lb.focus_set()
    varsbc = ['-All IBC-']+ViewAllDefinedIBC(CTK.t)

    # Remplace tous les elements
    lb.delete(0, TK.END)
    for i, value in enumerate(['-All IBC-']+ViewAllDefinedIBC(CTK.t)): lb.insert(i, value)
    return lb

#==============================================================================
# View specific IBC
#==============================================================================
def ViewIBC(event=None):
    if CTK.t == []: return

    nz  = len(Internal.getZones(CTK.t))
    nzs = CPlot.getActiveZones()
    active = [(i,1) for i in range(nz)]
    CPlot.setActiveZones(active)

    IBCTypes = []
    selection = WIDGETS['IBCLB'].curselection()
    for s in selection:
        t = WIDGETS['IBCLB'].get(s)
        IBCTypes.append(t)

    nzs = Internal.getNodesFromType2(CTK.t, 'Zone_t')
    ZoneLocal = []
    bases = CTK.t[2][1:]
    for b in bases:
        for zone in Internal.getZones(b):
            n = Internal.getNodeFromPath(zone, '.Solver#define/ibctype')

            if n is not None:
                val = Internal.getValue(n)
                if val not in IBCTypes and '-All IBC-' not in IBCTypes:
                    i = CPlot.getCPlotNumber(CTK.t, b[0], zone[0])
                    ZoneLocal.append((i,0))
            else:
                i = CPlot.getCPlotNumber(CTK.t, b[0], zone[0])
                ZoneLocal.append((i,0))
    CPlot.setActiveZones(ZoneLocal)
    CPlot.setState(ghostifyDeactivatedZones=1)
    return

#==============================================================================
# Get data from selected zone
#==============================================================================
def getData():
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        zone = Internal.getNodeFromType2(CTK.t, 'Zone_t')
    else: # get first of selection
        nz = nzs[0]
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        zone = CTK.t[2][nob][2][noz]
    if zone is not None:
        n = Internal.getNodeFromPath(zone, '.Solver#define/snear')
        if n is not None:
            val = Internal.getValue(n)
            VARS[0].set(val)
        else:
            VARS[0].set(-1000)
        n = Internal.getNodeFromPath(zone, '.Solver#define/ibctype')
        if n is not None:
            val = Internal.getValue(n)
            VARS[1].set(val)
        else:
            VARS[1].set('None')
        n = Internal.getNodeFromPath(zone, '.Solver#define/dfar')
        if n is not None:
            val = Internal.getValue(n)
            VARS[2].set(val)
        else:
            VARS[2].set(-1000)
        n = Internal.getNodeFromPath(zone, '.Solver#define/inv')
        if n is not None:
            val = Internal.getValue(n)
            if val == 0: VARS[3].set('out')
            else: VARS[3].set('in')
        n = Internal.getNodeFromPath(zone, '.Solver#define/extractWalls')
        if n is not None:
            val = Internal.getValue(n)
            if val == 0: VARS[4].set('0')
            else: VARS[4].set('1')
        n = Internal.getNodeFromPath(zone, '.Solver#define/extractLoads')
        if n is not None:
            val = Internal.getValue(n)
            if val == 0: VARS[5].set('0')
            else: VARS[5].set('1')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkIBC  [ + ]  ', font=CTK.FRAMEFONT,
                           takefocus=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkIBC')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Snear -
    V = TK.DoubleVar(win); V.set(0.01); VARS.append(V)
    # -1- IBC type -
    V = TK.StringVar(win); V.set('WallLaw'); VARS.append(V)
    # -2- dfar local -
    V = TK.DoubleVar(win); V.set(20.); VARS.append(V)
    # -3- mask inv or not -
    V = TK.StringVar(win); V.set('out'); VARS.append(V)
    # -4- extract fields on surface
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -5- extract loads for each component
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -6- Zones missing IBC info -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -7- Symmetry plane -
    V = TK.StringVar(win); V.set('Around XZ-'); VARS.append(V)

    # - Snear settings -
    B = TTK.Label(Frame, text="snear")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='The generated grid spacing for selected curve.')
    B = TTK.Entry(Frame, textvariable=VARS[0], width=4, background="White")
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)

    # - dfar settings  -
    B = TTK.Label(Frame, text="dfar")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='The distance from the center of object to the far boundary.\nIf set to -1, not taken into account.')
    B = TTK.Entry(Frame, textvariable=VARS[2], width=4, background="White")
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)

    # - IBC type -
    B = TTK.Label(Frame, text="IBC type")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of Immersed Boundary Condition.')
    B = TTK.OptionMenu(Frame, VARS[1], 'slip', 'noslip', 'Log', 'Musker', 'SA', 'outpress', 'inj', 'TBLE', 'slip_cr', 'overlap','wiremodel','MuskerLinear','SALinear','None')
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)

    # - Mask settings (in or out) -
    B = TTK.Label(Frame, text="Fluid")
    B.grid(row=3, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[3], 'out', 'in')
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)

    # - Symmetry plane -
    B = TTK.Button(Frame, text="Set symmetry plane", command=symmetrize)
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Create a symmetry plane.')
    B = TTK.OptionMenu(Frame, VARS[7], 'Around YZ-', 'Around XZ-', 'Around XY-')
    B.grid(row=4, column=1, columnspan=2, sticky=TK.EW)

    # - Extract fields on surface
    B = TTK.Label(Frame, text="Extract")
    B.grid(row=5, column=0, sticky=TK.EW)
    B = TTK.Checkbutton(Frame, text='Wall Fields', variable=VARS[4])
    BB = CTK.infoBulle(parent=B, text='Extract various fields on surface.')
    B.grid(row=5, column=1, columnspan=2, sticky=TK.EW)

    # - Extract loads on components
    B = TTK.Label(Frame, text="Extract")
    B.grid(row=6, column=0, sticky=TK.EW)
    B = TTK.Checkbutton(Frame, text='Loads', variable=VARS[5])
    BB = CTK.infoBulle(parent=B, text='Extract loads for each component.')
    B.grid(row=6, column=1, columnspan=2, sticky=TK.EW)

    # - Set data -
    B = TTK.Button(Frame, text="Set data", command=setData)
    BB = CTK.infoBulle(parent=B, text='Set data into selected zone.')
    B.grid(row=7, column=0, columnspan=1, sticky=TK.EW)

    # - Get data -
    B = TTK.Button(Frame, text="Get data", command=getData,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    BB = CTK.infoBulle(parent=B, text='Get data from selected zone.')
    B.grid(row=7, column=1, columnspan=2, sticky=TK.EW)

    # - View type de IBC -
    B = TTK.Button(Frame, text="View IBC", command=ViewIBC)
    B.grid(row=8, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='View specified IBC.\nTree is NOT modified.')

    # - Type of IBC - ListBox Frame -
    LBFrame = TTK.Frame(Frame)
    LBFrame.grid(row=8, column=1, rowspan=4, columnspan=2, sticky=TK.EW)
    LBFrame.rowconfigure(0, weight=1)
    LBFrame.columnconfigure(0, weight=1)
    LBFrame.columnconfigure(1, weight=0)
    SB  = TTK.Scrollbar(LBFrame)
    LB  = TTK.Listbox(LBFrame, selectmode=TK.EXTENDED, height=5)
    LB.bind('<Double-1>', ViewIBC)
    LB.bind('<Enter>', updateIBCNameList)
    for i, value in enumerate(['-All IBC-']+ViewAllDefinedIBC(CTK.t)): LB.insert(i, value)
    SB.config(command=LB.yview)
    LB.config(yscrollcommand=SB.set)
    LB.grid(row=0, column=0, sticky=TK.NSEW)
    SB.grid(row=0, column=1, sticky=TK.NSEW)
    WIDGETS['IBCLB'] = LB

    # - View Undefined IBC data -
    B = TTK.Button(Frame, text="View Undefined\n     IBC", command=ViewUndefinedIBC)
    WIDGETS['ViewUndefinedIBC'] = B
    BB = CTK.infoBulle(parent=B, text='View Undefined IBC.')
    B.grid(row=9, column=0, sticky=TK.EW)

    ## - Zones that are missing IBC info  -
    #B = TTK.Label(Frame, text="No IBC (Zones)")
    #B.grid(row=8, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Zones that are missing IBC info.')
    #B = TTK.Entry(Frame, textvariable=VARS[6], width=4, background="White")
    #B.grid(row=8, column=1, columnspan=2, sticky=TK.EW)


#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['BCNoteBook'].add(WIDGETS['frame'], text='tkIBC')
    except: pass
    CTK.WIDGETS['BCNoteBook'].select(WIDGETS['frame'])
    getData()

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['BCNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkIBC '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
