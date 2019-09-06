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

BODY = None

#==============================================================================
def changeMode(event=None):
    global BODY
    mode = VARS[10].get()
    if mode == 'Body':
        if BODY is not None: 
            CTK.t = BODY
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.TKTREE.updateApp()
            CTK.display(CTK.t)
        CTK.TXT.insert('START', 'Revert to body tree.\n')
    elif mode == 'PrevStep':
        # Reload from restart.cgns
        try: 
            CTK.t = C.convertFile2PyTree('restart.cgns')
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.TKTREE.updateApp()
            CTK.display(CTK.t)
            CTK.TXT.insert('START', 'Revert to previous solution step.\n')
        except: pass
        VARS[10].set('Main')
    else:
        CTK.TXT.insert('START', 'Display main computation tree.\n')

#==============================================================================
# Set data in selected zones
#==============================================================================
def setData():
    global BODY
    if CTK.t == []: return
    mode = VARS[10].get()
    if mode == 'Body': BODY = Internal.copyRef(CTK.t)

    temporal_scheme = VARS[0].get()
    ss_iteration = int(VARS[1].get())
    scheme = VARS[4].get()
    time_step = VARS[5].get()
    snear = VARS[6].get()
    ibctype = VARS[7].get()
    dfar = VARS[8].get()
    if VARS[12].get() == 'out': inv = 0
    else: inv = 1
    
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
            Internal.createUniqueChild(n, 'inv', 'DataArray_t', value=inv)

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
            Internal.createUniqueChild(n, 'inv', 'DataArray_t', value=inv)

    CTK.TXT.insert('START', 'Solver data set.\n')
    
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
        n = Internal.getNodeFromPath(zone, '.Solver#define/inv')
        if n is not None:
            val = Internal.getValue(n)
            if val == 0: VARS[12].set('out')
            else: VARS[12].set('in')

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
def run(event=None):
    mode = VARS[10].get()
    if mode == 'Body':
        global BODY
        BODY = Internal.copyRef(CTK.t)
        prepare() # save t, tc
        VARS[10].set('Main')
        CTK.t = CTK.upgradeTree(CTK.t)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CTK.display(CTK.t)
    
    # Set CPlot to scalar mode to monitor solution
    if CPlot.getState('mode') == 0: # mesh
        CPlot.setState(mode=3, scalarField='Density')

    compute()
    CTK.display(CTK.t)

#==============================================================================
# A partir de CTK.t considere comme les bodies
def prepare():
    if CTK.t == []: return

    # Save preventif
    C.convertPyTree2File(CTK.t, 'body.cgns')

    # Recupere la base REFINE
    b = Internal.getNodeFromName1(CTK.t, 'REFINE')
    if b is not None:
        tbox = C.newPyTree()
        tbox[2].append(b)
        Internal._rmNodesFromName1(CTK.t, 'REFINE')
    else: tbox = None

    import Apps.Fast.IBM as App
    myApp = App.IBM(format='single')
    CTK.t, tc = myApp.prepare(CTK.t, t_out='t.cgns', tc_out='tc.cgns', vmin=21, tbox=tbox, check=False)
    return None

#==============================================================================
# Lance des iterations
def compute():
    if CTK.t == []: return

    import Apps.Fast.IBM as App
    import Converter.Internal as Internal
    import Converter.PyTree as C

    # Save preventif
    C.convertPyTree2File(CTK.t, 'restart.cgns')

    temporal_scheme = VARS[0].get()
    scheme = VARS[4].get()
    a = VARS[11].get()
    if a == 'cfl': time_step_nature = 'local'
    else: time_step_nature = 'global'
    val = float(VARS[5].get())
    if time_step_nature == 'local': ss_iteration = 1
    else: ss_iteration = 20

    myApp = App.IBM(format='single')
    myApp.set(numb={
    "temporal_scheme": temporal_scheme,
    "ss_iteration": ss_iteration,
    "omp_mode": 1,
    "modulo_verif": 50
    })

    myApp.set(numz={
    "time_step": val,
    "scheme": scheme,
    "time_step_nature": time_step_nature,
    "cfl": val
    })

    nit = VARS[9].get()
    moduloVerif = 50

    # Compute
    #CTK.t, tc = myApp.compute('restart.cgns', 'tc.cgns', t_out=None, tc_out='tc_restart.cgns', nit=nit)
    
    # open compute
    CTK.t, tc, ts, metrics, graph = myApp.setup('restart.cgns', 'tc.cgns')

    import FastS.PyTree as FastS
    
    it0 = 0; time0 = 0.
    first = Internal.getNodeFromName1(CTK.t, 'Iteration')
    if first is not None: it0 = Internal.getValue(first)
    first = Internal.getNodeFromName1(CTK.t, 'Time')
    if first is not None: time0 = Internal.getValue(first)
    time_step = Internal.getNodeFromName(CTK.t, 'time_step')
    time_step = Internal.getValue(time_step)

    for it in range(1,nit+1):
        FastS._compute(CTK.t, metrics, it, tc, graph)
        time0 += time_step
        if it%50 == 0:
            CTK.TXT.insert('START', '%d / %d - %f\n'%(it+it0,it0+nit,time0))
            CTK.TXT.update()
        if it%moduloVerif == 0:
            FastS.display_temporal_criteria(CTK.t, metrics, it, format='single')
            #CTK.display(CTK.t)
    Internal.createUniqueChild(CTK.t, 'Iteration', 'DataArray_t', value=it0+nit)
    Internal.createUniqueChild(CTK.t, 'Time', 'DataArray_t', value=time0)
    
    return None

#===============================================================
def writeFiles():
    writePrepFile()
    writeComputeFile()
    CTK.TXT.insert('START', 'Write prep.py and compute.py.\n')
    return None

#==============================================================================
# Write prep file
#==============================================================================
def writePrepFile():
    if CTK.t == []: return

    mode = VARS[10].get()
    if mode == 'Body': tbody = CTK.t
    elif BODY is not None: tbody = BODY

    # Recupere la base REFINE si presente
    b = Internal.getNodeFromName1(tbody, 'REFINE')
    if b is not None:
        tbox = C.newPyTree()
        tbox[2].append(b)
        tbody = Internal.rmNodesFromName1(tbody, 'REFINE')
    else: tbox = None
    
    # Save preventif
    C.convertPyTree2File(tbody, 'body.cgns')
    C.convertPyTree2File(tbox, 'tbox.cgns')
    
    f = open('prep.py', 'w')
    
    text= """
import Apps.Fast.IBM as App
myApp = App.IBM(format='single')
"""

    if tbox is None: text +="myApp.prepare('body.cgns', t_out='t.cgns', tc_out='tc.cgns', check=False)"
    else: text += "myApp.prepare('body.cgns', t_out='t.cgns', tbox='tbox.cgns', tc_out='tc.cgns', check=False)"

    f.write(text)
    CTK.TXT.insert('START', 'File prep.py written.\n')
    f.close()

#==============================================================================
# Write compute file
#==============================================================================
def writeComputeFile():
    if CTK.t == []: return
            
    temporal_scheme = VARS[0].get()
    scheme = VARS[4].get()
    a = VARS[11].get()
    if a == 'cfl': time_step_nature = 'local'
    else: time_step_nature = 'global'
    val = float(VARS[5].get())
    if time_step_nature == 'global': time_step = val; cfl = 4.
    else: time_step = 0.1; cfl = val
    if time_step_nature == 'local': ss_iteration = 3
    else: ss_iteration = 30
    nit = VARS[9].get()

    f = open('compute.py', 'w')
    
    text= """
import Apps.Fast.IBM as App

myApp = App.IBM(format='single')
myApp.set(numb={
    "temporal_scheme": "%s",
    "ss_iteration": %d,
    "omp_mode": 1,
    "modulo_verif": 50
})
myApp.set(numz={
    "time_step": %g,
    "scheme": "%s",
    "time_step_nature": "%s",
    "cfl": %g,
})

# Compute
myApp.compute('t.cgns', 'tc.cgns', t_out='restart.cgns', tc_out='tc_restart.cgns', nit=%d)
"""%(temporal_scheme, ss_iteration, time_step, scheme, time_step_nature, cfl, nit)

    f.write(text)
    CTK.TXT.insert('START', 'File compute.py written.\n')
    f.close()

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
    V = TK.StringVar(win); V.set('senseur'); VARS.append(V)
    # -5- Time step -
    V = TK.DoubleVar(win); V.set(0.004); VARS.append(V)
    # -6- Snear -
    V = TK.DoubleVar(win); V.set(0.01); VARS.append(V)
    # -7- IBC type -
    V = TK.StringVar(win); V.set('Musker'); VARS.append(V)
    # -8- dfar local -
    V = TK.DoubleVar(win); V.set(20.); VARS.append(V)
    # -9- nbre d'iterations -
    V = TK.IntVar(win); V.set(100); VARS.append(V)
    # -10 - body or main mode
    V = TK.StringVar(win); V.set('Body'); VARS.append(V)
    # -11- time_step or cfl
    V = TK.StringVar(win); V.set('time_step'); VARS.append(V)
    # -12- mask inv or not -
    V = TK.StringVar(win); V.set('out'); VARS.append(V)

    #- Mode -
    B = TTK.Label(Frame, text="Mode")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Currently viewed tree.\nIn body mode:\nOne Base per body\nBase REFINE for refinement zones.')
    B = TTK.OptionMenu(Frame, VARS[10], 'Body', 'Main', 'PrevStep', command=changeMode)
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)

    #- Snear settings  -
    B = TTK.Label(Frame, text="snear")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='The generated grid spacing for selected curve.')
    B = TTK.Entry(Frame, textvariable=VARS[6], width=4, background="White")
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    
    #- dfar settings  -
    B = TTK.Label(Frame, text="dfar")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='The distance from the center of object to the far boundary.\nIf set to -1, not taken into account.')
    B = TTK.Entry(Frame, textvariable=VARS[8], width=4, background="White")
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)

    # - IBC type -
    B = TTK.Label(Frame, text="IBC type")
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Type of Immersed boundary condition.')
    B = TTK.OptionMenu(Frame, VARS[7], 'slip', 'noslip', 'Log', 'Musker', 'outpress', 'inj', 'TBLE')
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)

    #- Mask settings (in or out)  -
    B = TTK.Label(Frame, text="Fluid")
    B.grid(row=4, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[12], 'out', 'in')
    B.grid(row=4, column=1, columnspan=2, sticky=TK.EW)

    # - Set data -
    B = TTK.Button(Frame, text="Set data", command=setData)
    BB = CTK.infoBulle(parent=B, text='Set data into selected zone.')
    B.grid(row=5, column=0, columnspan=2, sticky=TK.EW)
    B = TTK.Button(Frame, command=getData,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    BB = CTK.infoBulle(parent=B, text='Get data from selected zone.')
    B.grid(row=5, column=2, sticky=TK.EW)

    # - temporal scheme -
    B = TTK.Label(Frame, text="time_scheme")
    B.grid(row=6, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Time integration.')
    B = TTK.OptionMenu(Frame, VARS[0], 'explicit', 'implicit', 'implicit_local')
    B.grid(row=6, column=1, columnspan=2, sticky=TK.EW)

    # - ss_iteration -
    #B = TTK.Label(Frame, text="ss_iteration")
    #B.grid(row=1, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Nbre de sous iterations max.')
    #B = TTK.Entry(Frame, textvariable=VARS[1])
    #B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)

    # - scheme -
    B = TTK.Label(Frame, text="scheme")
    B.grid(row=7, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Numerical scheme.')
    B = TTK.OptionMenu(Frame, VARS[4], 'roe_min', 'ausmpred', 'senseur')
    B.grid(row=7, column=1, columnspan=2, sticky=TK.EW)

    # - time_step -
    B = TTK.OptionMenu(Frame, VARS[11], "time_step", "cfl")
    B.grid(row=8, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Time step.')
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White')
    B.grid(row=8, column=1, columnspan=2, sticky=TK.EW)

    # - compute -
    B = TTK.Button(Frame, text="Compute", command=run)
    BB = CTK.infoBulle(parent=B, text='Launch computation.')
    B.grid(row=9, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[9], width=5, background='White')
    B.grid(row=9, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Button(Frame, text="Files", command=writeFiles)
    BB = CTK.infoBulle(parent=B, text='Write files to run elsewhere.')
    B.grid(row=9, column=2, sticky=TK.EW)
    

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
    (win, menu, file, tools) = CTK.minimal('tkFastSolver '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
