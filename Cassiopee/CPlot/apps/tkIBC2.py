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
def _setDataInZone(z, bLocal, snear, ibctype, dfar, inv):
    # set IBM data in .Solver#define
    
    if VARS[1].get()!='rectilinear':
        if Internal.getNodeFromName1(z, '.Solver#define'):
            Internal._rmNode(z,Internal.getNodeFromName1(z, '.Solver#define'))
        n = Internal.createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        Internal.createUniqueChild(n, 'snear', 'DataArray_t', value=snear)
        Internal.createUniqueChild(n, 'ibctype', 'DataArray_t', value=ibctype)
        Internal.createUniqueChild(n, 'dfar', 'DataArray_t', value=dfar)
        Internal.createUniqueChild(n, 'inv', 'DataArray_t', value=inv)
    else:
        if Internal.getNodeFromName1(bLocal, '.Solver#define'):
            Internal._rmNode(z,Internal.getNodeFromName1(bLocal, '.Solver#define'))
        n = Internal.createUniqueChild(bLocal, '.Solver#define', 'UserDefinedData_t')
        Internal.createUniqueChild(n, 'dirx'       , 'DataArray_t', value=int(VARS[12].get()))
        Internal.createUniqueChild(n, 'diry'       , 'DataArray_t', value=int(VARS[13].get()))
        Internal.createUniqueChild(n, 'dirz'       , 'DataArray_t', value=int(VARS[14].get()))
        if VARS[15].get()=='coarse(0)': granLocal = 0
        elif VARS[15].get()=='fine(1)': granLocal = 1
        Internal.createUniqueChild(n, 'granularity', 'DataArray_t', value=granLocal)

    if VARS[1].get() in ['noslip', 'Log', 'Musker', 'SA', 'TBLE', 'MuskerLinear','SALinear']:
        # Set Extractions triggers in .Solver#define
        if VARS[4].get() == "1": value = 1
        else: value = 0
        Internal.createUniqueChild(n, 'extractWalls', 'DataArray_t', value=value)
        if VARS[5].get() == "1": value = 1
        else: value = 0
        Internal.createUniqueChild(n, 'extractLoads', 'DataArray_t', value=value)
    if VARS[1].get()=='outpress':
        Internal.createUniqueChild(n, 'pStatic', 'DataArray_t'          , value=float(VARS[8].get()))
        Internal.createUniqueChild(n, 'isDensityConstant', 'DataArray_t', value=float(VARS[9].get()))
    if VARS[1].get()=='inj':
        Internal.createUniqueChild(n, 'StagnationPressure', 'DataArray_t', value=float(VARS[10].get()))
        Internal.createUniqueChild(n, 'StagnationEnthalpy', 'DataArray_t', value=float(VARS[11].get()))
        Internal.createUniqueChild(n, 'dirx', 'DataArray_t', value=float(VARS[12].get()))
        Internal.createUniqueChild(n, 'diry', 'DataArray_t', value=float(VARS[13].get()))
        Internal.createUniqueChild(n, 'dirz', 'DataArray_t', value=float(VARS[14].get()))
    if VARS[1].get()=='wmm':
        Internal.createUniqueChild(n, 'diameter', 'DataArray_t', value=float(VARS[12].get()))
        Internal.createUniqueChild(n, 'ct'      , 'DataArray_t', value=float(VARS[13].get()))
        Internal.createUniqueChild(n, 'k'       , 'DataArray_t', value=float(VARS[14].get()))
        
    # remesh surface eventually
    dim = Internal.getZoneDim(z)
    remesh = False
    if dim[0] == 'Structured' and dim[2] == 1 and dim[3] == 1: remesh = True
    if dim[0] == 'Unstructured' and dim[3] == 'BAR': remesh = True
    if remesh: D._uniformize(z, h=float(snear))
    return None



# Change the mode 
def setMode(event=None):
    mode = VARS[1].get()
    if mode == 'wall': VARS[1].set('noslip')
    if mode in ['noslip', 'Log', 'Musker', 'SA', 'TBLE', 'MuskerLinear','SALinear']: mode = 'wall'
    imode = 0
    WIDGETS['slip'].grid_forget()
    WIDGETS['wall'].grid_forget()
    WIDGETS['outpress'].grid_forget()
    WIDGETS['inj'].grid_forget()
    WIDGETS['rec'].grid_forget()
    WIDGETS['wmm'].grid_forget()

    if mode == 'slip':
        imode = 0; WIDGETS['slip'].grid(row=10, column=0, columnspan=2, sticky=TK.EW)
    elif mode == 'wall':
        imode = 0; WIDGETS['wall'].grid(row=10, column=0, columnspan=4, sticky=TK.EW)
    elif mode == 'outpress':
        imode = 0; WIDGETS['outpress'].grid(row=10, column=0, columnspan=3, sticky=TK.EW)
    elif mode == 'inj':
        imode = 0; WIDGETS['inj'].grid(row=10, column=0, columnspan=3, sticky=TK.EW)
    elif mode == 'rectilinear':
        imode = 0; WIDGETS['rec'].grid(row=10, column=0, columnspan=3, sticky=TK.EW)
    elif mode == 'wiremodel':
        imode = 0; WIDGETS['wmm'].grid(row=10, column=0, columnspan=3, sticky=TK.EW)  
    CPlot.setState(mode=imode)
    CTK.TXT.insert('START', 'Mode %s displayed.\n'%mode)

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
        bLocal = CTK.t[2][nob]
        z      = bLocal[2][noz]
        _setDataInZone(z, bLocal, snear, ibctype, dfar, inv)
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
            setMode()
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
                           text='tkIBC2  [ + ]  ', font=CTK.FRAMEFONT, 
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
    CTK.addPinMenu(FrameMenu, 'tkIBC2')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Snear -
    V = TK.DoubleVar(win); V.set(0.01); VARS.append(V)
    # -1- IBC type -
    V = TK.StringVar(win); V.set('slip'); VARS.append(V)
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
    # -8- OutletPressure
    V = TK.DoubleVar(win); V.set(101325); VARS.append(V)
    # -9- Constant outlet density
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -10- Injection total Pressure
    V = TK.DoubleVar(win); V.set(101325); VARS.append(V)
    # -11- Injection total enthalpy
    V = TK.DoubleVar(win); V.set(101325); VARS.append(V)
    # -12- Injection dirx
    V = TK.DoubleVar(win); V.set(1); VARS.append(V)
    # -13- Injection diry
    V = TK.DoubleVar(win); V.set(0); VARS.append(V)
    # -14- Injection dirz
    V = TK.DoubleVar(win); V.set(0); VARS.append(V)
    # -15- Rectilinear granularity
    V = TK.StringVar(win); V.set('coarse'); VARS.append(V)

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
    B = TTK.OptionMenu(Frame, VARS[1], 'slip', 'wall', 'outpress', 'inj', 'slip_cr', 'overlap','wiremodel','rectilinear','None', command=setMode)
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)

    # - Mask settings (in or out) -
    B = TTK.Label(Frame, text="Fluid")
    B.grid(row=3, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[3], 'out', 'in')
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)

    # - Set data -
    B = TTK.Button(Frame, text="Set data", command=setData)
    BB = CTK.infoBulle(parent=B, text='Set data into selected zone.')
    B.grid(row=4, column=0, columnspan=1, sticky=TK.EW)

    # - Get data -
    B = TTK.Button(Frame, text="Get data", command=getData,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    BB = CTK.infoBulle(parent=B, text='Get data from selected zone.')
    B.grid(row=4, column=1, columnspan=2, sticky=TK.EW)

    # - View type de IBC -
    B = TTK.Button(Frame, text="View IBC", command=ViewIBC)
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='View specified IBC.\nTree is NOT modified.')

    # - Type of IBC - ListBox Frame -
    LBFrame = TTK.Frame(Frame)
    LBFrame.grid(row=5, column=1, rowspan=4, columnspan=2, sticky=TK.EW)
    LBFrame.rowconfigure(0, weight=1)
    LBFrame.columnconfigure(0, weight=1)
    LBFrame.columnconfigure(1, weight=0)
    SB  = TTK.Scrollbar(LBFrame)
    LB  = TTK.Listbox(LBFrame, selectmode=TK.EXTENDED, height=5)
    LB.bind('<Double-1>', ViewIBC)
    LB.bind('<Enter>', updateIBCNameList)
    for i, value in enumerate(['-All IBC-']+ViewAllDefinedIBC(CTK.t)): LB.insert(i, value)
    SB.config(command = LB.yview)
    LB.config(yscrollcommand = SB.set)
    LB.grid(row=0, column=0, sticky=TK.NSEW)
    SB.grid(row=0, column=1, sticky=TK.NSEW)
    WIDGETS['IBCLB'] = LB

    # - View Undefined IBC data -
    B = TTK.Button(Frame, text="View Undefined\n     IBC", command=ViewUndefinedIBC)
    WIDGETS['ViewUndefinedIBC'] = B
    BB = CTK.infoBulle(parent=B, text='View Undefined IBC.')
    B.grid(row=6, column=0, sticky=TK.EW)


    ## WIDGETS THAT APPEAR & DISAPPEAR

    # - Symmetry plane -
    slip = TTK.LabelFrame(Frame, borderwidth=2, relief="solid", text="Slip Parameters:")
    slip.columnconfigure(0, weight=1)
    slip.columnconfigure(1, weight=1)
    slip.grid(row=10, column=0, columnspan=2)
    WIDGETS['slip'] = slip
    
    B = TTK.Button(slip, text="Set symmetry plane", command=symmetrize)
    BB = CTK.infoBulle(parent=B, text='Create a symmetry plane.')
    B.grid(row=10, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(slip, VARS[7], 'Around YZ-', 'Around XZ-', 'Around XY-')
    B.grid(row=10, column=1, sticky=TK.EW)


    ## Wall Conditions
    wall = TTK.LabelFrame(Frame, borderwidth=2, relief="solid", text="Wall Parameters:")
    wall.columnconfigure(0, weight=0)
    wall.columnconfigure(1, weight=0)
    wall.columnconfigure(2, weight=0)
    wall.columnconfigure(3, weight=0)
    wall.grid(row=10, column=0, columnspan=4)
    WIDGETS['wall'] = wall
    
    # - Extract fields on surface
    B = TTK.Label(wall, text="Extract")
    B.grid(row=10, column=0, sticky=TK.EW)
    B = TTK.Checkbutton(wall, text='Wall Fields', variable=VARS[4])
    B.grid(row=10, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract various fields on surface.')
    
    # - Extract loads on components
    B = TTK.Checkbutton(wall, text='Loads', variable=VARS[5])
    B.grid(row=10, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract loads for each component.')

    B = TTK.Label(wall, text="Wall Type")
    B.grid(row=11, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of Immersed Boundary Condition Wall Conditions.')
    B = TTK.OptionMenu(wall, VARS[1], 'noslip', 'Log', 'Musker', 'SA', 'TBLE', 'MuskerLinear','SALinear')
    B.grid(row=11, column=1, columnspan=2, sticky=TK.EW)


    ## Outlet Pressure
    outpress = TTK.LabelFrame(Frame, borderwidth=2, relief="solid", text="Outlet Pressure Parameters:")
    outpress.columnconfigure(0, weight=0)
    outpress.columnconfigure(1, weight=1)
    outpress.columnconfigure(2, weight=0)
    outpress.grid(row=10, column=0, columnspan=3)
    WIDGETS['outpress'] = outpress
    
    B = TTK.Label(outpress, text="Static Pressure")
    B.grid(row=10, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Static Pressure at outlet')
    B = TTK.Entry(outpress, textvariable=VARS[8], width=4, background="White")
    B.grid(row=10, column=1, columnspan=2, sticky=TK.EW)

    B = TTK.Checkbutton(outpress, text='Constant density at outlet', variable=VARS[9])
    B.grid(row=11, column=0, columnspan=2,sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set a fixed outlet density')

    ## Injection
    inj = TTK.LabelFrame(Frame, borderwidth=2, relief="solid", text="Injection Parameters:")
    inj.columnconfigure(0, weight=0)
    inj.columnconfigure(1, weight=1)
    inj.columnconfigure(2, weight=0)
    inj.grid(row=10, column=0, columnspan=3)
    WIDGETS['inj'] = inj
    
    B = TTK.Label(inj, text="Total Pressure")
    B.grid(row=10, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Total Pressure at injection inlet')
    B = TTK.Entry(inj, textvariable=VARS[10], width=4, background="White")
    B.grid(row=10, column=1, columnspan=2, sticky=TK.EW)

    B = TTK.Label(inj, text="Total Enthalpy")
    B.grid(row=11, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Total Enthalpy at injection inlet')
    B = TTK.Entry(inj, textvariable=VARS[11], width=4, background="White")
    B.grid(row=11, column=1, columnspan=2, sticky=TK.EW)

    B = TTK.Label(inj, text="Unit Normal: x")
    B.grid(row=12, column=0, sticky=TK.EW)
    B = TTK.Entry(inj, textvariable=VARS[12], width=4, background="White")
    B.grid(row=12, column=1, columnspan=2, sticky=TK.EW)

    B = TTK.Label(inj, text="Unit Normal: y")
    B.grid(row=13, column=0, sticky=TK.EW)
    B = TTK.Entry(inj, textvariable=VARS[13], width=4, background="White")
    B.grid(row=13, column=1, columnspan=2, sticky=TK.EW)

    B = TTK.Label(inj, text="Unit Normal: z")
    B.grid(row=14, column=0, sticky=TK.EW)
    B = TTK.Entry(inj, textvariable=VARS[14], width=4, background="White")
    B.grid(row=14, column=1, columnspan=2, sticky=TK.EW)

    ## Rectilienar
    rec = TTK.LabelFrame(Frame, borderwidth=2, relief="solid", text="Rectilinear Parameters:")
    rec.columnconfigure(0, weight=1)
    rec.columnconfigure(1, weight=1)
    rec.columnconfigure(2, weight=1)
    rec.grid(row=10, column=0, columnspan=3)
    WIDGETS['rec'] = rec

    B = TTK.Label(rec, text="x-dir")
    B.grid(row=10, column=0, sticky=TK.EW)
    B = TTK.Label(rec, text="y-dir")
    B.grid(row=10, column=1, sticky=TK.EW)
    B = TTK.Label(rec, text="z-dir")
    B.grid(row=10, column=2, sticky=TK.EW)

    B = TTK.Entry(rec, textvariable=VARS[12], width=4, background="White")
    B.grid(row=11, column=0, sticky=TK.EW)
    B = TTK.Entry(rec, textvariable=VARS[13], width=4, background="White")
    B.grid(row=11, column=1, sticky=TK.EW)
    B = TTK.Entry(rec, textvariable=VARS[14], width=4, background="White")
    B.grid(row=11, column=2, sticky=TK.EW)

    B = TTK.Label(rec, text="Rectilinear granularity")
    B.grid(row=12, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Specify the granularity of the rectilinear approach.')
    B = TTK.OptionMenu(rec, VARS[15], 'coarse(0)', 'fine(1)')
    B.grid(row=12, column=1, sticky=TK.EW)

    ## Wire Mesh Model
    wmm = TTK.LabelFrame(Frame, borderwidth=2, relief="solid", text="Wire Mesh Model Parameters:")
    wmm.columnconfigure(0, weight=1)
    wmm.columnconfigure(1, weight=1)
    wmm.columnconfigure(2, weight=1)
    wmm.grid(row=10, column=0, columnspan=3)
    WIDGETS['wmm'] = wmm

    B = TTK.Label(wmm, text="DiameterWire")
    B.grid(row=10, column=0, sticky=TK.EW)
    B = TTK.Label(wmm, text="CtWire")
    B.grid(row=10, column=1, sticky=TK.EW)
    B = TTK.Label(wmm, text="KWire")
    B.grid(row=10, column=2, sticky=TK.EW)

    B = TTK.Entry(wmm, textvariable=VARS[12], width=4, background="White")
    B.grid(row=11, column=0, sticky=TK.EW)
    B = TTK.Entry(wmm, textvariable=VARS[13], width=4, background="White")
    B.grid(row=11, column=1, sticky=TK.EW)
    B = TTK.Entry(wmm, textvariable=VARS[14], width=4, background="White")
    B.grid(row=11, column=2, sticky=TK.EW)

    if 'tkViewMode' in CTK.PREFS: setMode()
    
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
    try: CTK.WIDGETS['BCNoteBook'].add(WIDGETS['frame'], text='tkIBC2')
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
    (win, menu, file, tools) = CTK.minimal('tkIBC2 '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
