# - Fast solvers app -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Fast.PyTree as Fast
import CPlot.iconics as iconics
import os, math, os.path
# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setData():
    if CTK.t == []: return
    
    temporal_scheme = VARS[0].get()
    ss_iteration = int(VARS[1].get())
    scheme = VARS[4].get()
    time_step = VARS[5].get()
    snear = VARS[6].get()
    ibctype = VARS[7].get()
    dfar = VARS[8].get()

    numb = {'temporal_scheme':temporal_scheme,
            'ss_iteration':ss_iteration}
    numz = {'scheme':scheme, 'time_step':time_step}

    nzs = CPlot.getSelectedZones()
    CTK.saveTree()
    if nzs == []:
        Fast._setNum2Base(CTK.t, numb)
        Fast._setNum2Zones(CTK.t, numz)
        zones = Internal.getZones(CTK.t)
        for z in zones:
            n = Internal.createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
            Internal.createUniqueChild(n, 'snear', 'DataArray_t', value=snear)
            Internal.createUniqueChild(n, 'ibctype', 'DataArray_t', value=ibctype)
            Internal.createUniqueChild(n, 'dfar', 'DataArray_t', value=dfar)

    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            b, c = Internal.getParentOfNode(CTK.t, z)
            Fast._setNum2Base(b, numb)
            Fast._setNum2Zones(z, numz)
            n = Internal.createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
            Internal.createUniqueChild(n, 'snear', 'DataArray_t', snear)
            Internal.createUniqueChild(n, 'ibctype', 'DataArray_t', ibctype)
            Internal.createUniqueChild(n, 'dfar', 'DataArray_t', value=dfar)

    CTK.TXT.insert('START', 'Solver data set.\n')
    
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
        n = Internal.getNodeFromPath(zone, '.Solver#define/time_step')
        if n is not None:
            val = Internal.getValue(n)
            VARS[5].set(val)
        n = Internal.getNodeFromPath(zone, '.Solver#define/snear')
        if n is not None:
            val = Internal.getValue(n)
            VARS[6].set(val)
        n = Internal.getNodeFromPath(zone, '.Solver#define/ibctype')
        if n is not None:
            val = Internal.getValue(n)
            VARS[7].set(val)
        n = Internal.getNodeFromPath(zone, '.Solver#define/dfar')
        if n is not None:
            val = Internal.getValue(n)
            VARS[8].set(val)

        d, c = Internal.getParentOfNode(CTK.t, zone)
        n = Internal.getNodeFromPath(d, '.Solver#define/temporal_scheme')
        if n is not None:
            val = Internal.getValue(n)
            VARS[0].set(val)
        n = Internal.getNodeFromPath(zone, '.Solver#define/scheme')
        if n is not None:
            val = Internal.getValue(n)
            VARS[4].set(val)

#==============================================================================
# Write setup file
#==============================================================================
def writeSetupFile():
    if CTK.t == []: return
    
    # EquationDimension
    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: dim = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); dim = 3

    # GoverningEquations
    nodes = Internal.getNodesFromName(CTK.t, 'GoverningEquations')
    if nodes != []:
        equations = Internal.getValue(nodes[0])
        if equations == 'Euler': model = 'euler'
        if equations == 'NSLaminar': model = 'nslam'
        if equations == 'NSTurbulent': model = 'nstur'
    else:
        CTK.TXT.insert('START', 'GoverningEquations is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    # Turbulence model
    if model == 'nstur':
        nodes = Internal.getNodesFromName(CTK.t, 'TurbulenceModel')
        if nodes != []:
            tm = Internal.getValue(nodes[0])
            if tm == 'OneEquation_SpalartAllmaras': model += '_sa'
            elif tm == 'TwoEquation_Wilcox': model += '_kw'
            elif tm == 'TwoEquation_MenterSST': model += '_kw'
            else:
                CTK.TXT.insert('START', 'This turbulence model is not accepted by Cassiopee solver.\n')
                CTK.TXT.insert('START', 'Error: ', 'Error')
                return
    
    # ReferenceState
    nodes = Internal.getNodesFromName(CTK.t, 'ReferenceState')
    if nodes == []:
        CTK.TXT.insert('START', 'Reference state is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    state = nodes[0]

    # Mach
    nodes = Internal.getNodesFromName(state, 'Mach')
    if nodes != []: Mach = Internal.getValue(nodes[0])
    else:
        CTK.TXT.insert('START', 'Mach is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    
    # Reynolds
    nodes = Internal.getNodesFromName(state, 'Reynolds')
    if nodes != []: Reynolds = Internal.getValue(nodes[0])
    elif equations == 'NSLaminar' or equations == 'NSTurbulent':
        CTK.TXT.insert('START', 'Reynolds is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    else: Reynolds = 1.
    if Reynolds <= 0.: Reynolds = 1.
    
    # Incidences
    node = Internal.getNodeFromName(state, 'VelocityX')
    if node is not None: UInf = Internal.getValue(node)
    else: UInf = 0.
    node = Internal.getNodeFromName(state, 'VelocityY')
    if node is not None: VInf = Internal.getValue(node)
    else: VInf = 0.
    node = Internal.getNodeFromName(state, 'VelocityZ')
    if node is not None: WInf = Internal.getValue(node)
    else: WInf = 0.
    if (UInf != 0.):
        aly = math.atan(WInf/UInf)
        alz = math.atan( math.cos(aly)*VInf/UInf )
    else:
        aly = 0.; alz = 0.
    alphaZ = alz*180./math.pi
    alphaY = aly*180./math.pi
    
    # tree file
    treeFile = os.path.splitext(CTK.FILE)[0]+'.cgns'
    
    f = open('setup.py', 'w')
    
    text = [
'import Converter.PyTree as C\n',
'import Fast.PyTree as Fast\n',
'import FastS.PyTree as FastS\n',
]

    f.write(text)
    CTK.TXT.insert('START', 'File setup.py written.\n')
    f.close()

#==============================================================================
# IN: le maillage
# IN: les BCs (match + autres)
# IN: Distribution (option)
# IN: Condition initiale (option)
# Genere les fichiers :  
# - setup.py (script de calcul)
# - t.cgns : maillage + solution initiale + 
# - tD.cgns: donneurs
#==============================================================================
def generate():
    writeSetupFile()
    #generateT()
    #generateTD()
    return

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFastSolver', font=CTK.FRAMEFONT, 
                           takefocus=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkFastSolver')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- temporal_scheme -
    V = TK.StringVar(win); V.set('explicit'); VARS.append(V)
    # -1- ss_iteration -
    V = TK.StringVar(win); V.set('20'); VARS.append(V)
    # -2- modulo_verif -
    V = TK.StringVar(win); V.set('200'); VARS.append(V)
    # -3- restart_fields -
    V = TK.StringVar(win); V.set('1');VARS.append(V)
    # -4- scheme -
    V = TK.StringVar(win); V.set('ausmpred'); VARS.append(V)
    # -5- Time step -
    V = TK.DoubleVar(win); V.set(0.002); VARS.append(V)
    # -6- Snear -
    V = TK.DoubleVar(win); V.set(0.01); VARS.append(V)
    # -7- IBC type -
    V = TK.StringVar(win); V.set('Musker'); VARS.append(V)
    # -8- dfar local -
    V = TK.DoubleVar(win); V.set(-1.); VARS.append(V)


    # - temporal scheme -
    B = TTK.Label(Frame, text="temporal_scheme")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Time integration.')
    B = TTK.OptionMenu(Frame, VARS[0], 'explicit', 'implicit', 'implicit_local')
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)

    # - ss_iteration -
    #B = TTK.Label(Frame, text="ss_iteration")
    #B.grid(row=1, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Nbre de sous iterations max.')
    #B = TTK.Entry(Frame, textvariable=VARS[1])
    #B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)

    # - scheme -
    B = TTK.Label(Frame, text="scheme")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Numerical scheme.')
    B = TTK.OptionMenu(Frame, VARS[4], 'ausmpred', 'roe_min')
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)

    # - time_step -
    B = TTK.Label(Frame, text="time_step")
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Time step.')
    B = TTK.Entry(Frame, textvariable=VARS[5])
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)
    
    #- Snear settings  -
    B = TTK.Label(Frame, text="snear")
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[6])
    B.grid(row=4, column=1, columnspan=2, sticky=TK.EW)
    
    #- dfar settings  -
    B = TTK.Label(Frame, text="dfar")
    B.grid(row=5, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[8])
    B.grid(row=5, column=1, columnspan=2, sticky=TK.EW)

    # - IBC type -
    B = TTK.Label(Frame, text="IBC type")
    B.grid(row=6, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of IBC.')
    B = TTK.OptionMenu(Frame, VARS[7], 'slip', 'noslip', 'Log', 'Musker', 'outpress', 'inj', 'TBLE')
    B.grid(row=6, column=1, columnspan=2, sticky=TK.EW)

    # - Set data -
    B = TTK.Button(Frame, text="Set data", command=setData)
    BB = CTK.infoBulle(parent=B, text='Set data into selected zone.')
    B.grid(row=7, column=0, columnspan=2, sticky=TK.EW)
    B = TTK.Button(Frame, command=getData,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    BB = CTK.infoBulle(parent=B, text='Get data from selected zone.')
    B.grid(row=7, column=2, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW)

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
    (win, menu, file, tools) = CTK.minimal('tkFastSolver '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
