# - tkDriver -
"""GUI for Roms/Driver."""
import tkinter as TK
import CPlot.Ttk as TTK
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.PyTree as C

# local widgets list
WIDGETS = {}; VARS = []

# DRIVER
DRIVER = None
ENTITY = None
MODE = 0 # 0: mesh, 1: dxdmu

#==============================================================================
# set driver and current entity
def setDriver(driver, entity):
    global DRIVER, ENTITY
    DRIVER = driver
    ENTITY = entity
    return None

#==============================================================================
# get the free params as set in driver
def getParamFromDriver():
    params = {}
    for f in DRIVER.freeParams:
        params[f.name] = DRIVER.scalars[f.name].v
    return params

#==============================================================================
# get the entities as set in driver
def getEntitiesFromDriver():
    entities = []
    #for f in DRIVER.edges:
    #    entities.append(f)
    for f in DRIVER.sketches:
        entities.append(f)
    for f in DRIVER.surfaces:
        entities.append(f)
    return entities

#==============================================================================
def setParameterName(event=None):
    # set parameter value
    paramName = VARS[0].get()
    params = getParamFromDriver()
    r = DRIVER.scalars[paramName].range
    value = params[paramName] # current value
    VARS[1].set(str(value))
    delta = max(r[1]-r[0], 1.e-6)
    alpha = 100./delta
    beta = 100. - alpha*r[1]
    scaleValue = beta + alpha*value
    WIDGETS['valueSlider'].set(scaleValue)
    update()
    return None

#==============================================================================
def setParameterValue(event=None):
    paramName = VARS[0].get()
    r = DRIVER.scalars[paramName].range
    value = CTK.varsFromWidget(VARS[1].get(), type=1)[0]
    delta = max(r[1]-r[0], 1.e-6)
    alpha = 100./delta
    beta = 100. - alpha*r[1]
    scaleValue = beta + alpha*value
    WIDGETS['valueSlider'].set(scaleValue)
    update()
    return None

#==============================================================================
def setParameterValueWithScale(event=None):
    val = WIDGETS['valueSlider'].get()
    paramName = VARS[0].get()
    r = DRIVER.scalars[paramName].range
    vmax = r[1]
    vmin = r[0]
    value = vmin + val*(vmax-vmin)/100.
    VARS[1].set(str(value))
    update()
    return None

#==============================================================================
def setEntityName(event=None):
    global ENTITY
    name = VARS[2].get()
    if name in DRIVER.surfaces:
        ENTITY = DRIVER.surfaces[name]
    elif name in DRIVER.sketches:
        ENTITY = DRIVER.sketches[name]
    elif name in DRIVER.edges:
        ENTITY = DRIVER.edges[name]
    update()

#==============================================================================
# update driver, hook, mesh from current parameter
def update(event=None):
    paramName = VARS[0].get()
    paramValue = CTK.varsFromWidget(VARS[1].get(), type=1)[0]
    params = getParamFromDriver()
    params[paramName] = paramValue
    DRIVER.instantiate(params)
    if MODE == 0: # mesh
        if ENTITY.h is None:
            ENTITY.h = (0.1, 0.1, 0.01) # a mieux regler automatiquement
        m = ENTITY.Mesh(method=0)
        CPlot.display(m, mode='Mesh', meshStyle=4, bgColor=1)
    elif MODE == 1: # dXdmu
        if ENTITY.h is None:
            ENTITY.h = (0.1, 0.1, 0.01) # a mieux regler automatiquement
        m = ENTITY.MeshAsReference()
        DRIVER._dXdmu(ENTITY, Mesh=m, freeParams=[VARS[0].get()], deps=1.e-4)
        CPlot.display(m, mode='Vector',
                      vectorField1='dXd0', vectorField2='dYd0', vectorField3='dZd0',
                      vectorStyle=1, bgColor=1)

#==============================================================================
# switch mode:
# 0: display mesh
# 1: display dXdmu
def switchMode(event=None):
    global MODE
    MODE += 1
    if MODE == 2: MODE = 0
    if MODE == 0: WIDGETS['Mode'].config(text='dXdmu')
    elif MODE == 1: WIDGETS['Mode'].config(text='mesh')
    update()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkDriver  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage container names.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=4)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkDriver')
    WIDGETS['frameMenu'] = FrameMenu

    # Get default values from driver
    entities = getEntitiesFromDriver()
    freeParams = DRIVER.freeParams
    params = getParamFromDriver()
    paramName = list(params.keys())[0]
    r = DRIVER.scalars[paramName].range
    paramValue = params[paramName]
    delta = max(r[1]-r[0], 1.e-6)
    alpha = 100./delta
    beta = 100. - alpha*r[1]
    scaleValue = beta + alpha*paramValue

    # - VARS -
    # -0- Current parameter name -
    V = TK.StringVar(win); V.set(paramName); VARS.append(V)

    # -1- Current parameter value -
    V = TK.StringVar(win); V.set(paramValue); VARS.append(V)

    # -2- Current entity
    V = TK.StringVar(win); V.set(ENTITY.name); VARS.append(V)

    # - Entity chooser -
    B = TTK.Label(Frame, text="Entity")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Current entity name.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.Entry(F, textvariable=VARS[2], background='White')
        B.grid(sticky=TK.EW)
        F.bind('<Return>', setEntityName)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['EntityName'] = B
    else:
        B = TTK.Combobox(F, textvariable=VARS[2],
                         values=entities, state='readonly')
        B.bind("<<ComboboxSelected>>", setEntityName)
        B.grid(sticky=TK.EW)
        B.bind('<Return>', setEntityName)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['ParameterName'] = B

    # - Parameter chooser -
    B = TTK.Label(Frame, text="Parameter")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Current parameter name.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.Entry(F, textvariable=VARS[0], background='White')
        B.grid(sticky=TK.EW)
        F.bind('<Return>', setParameterName)
        F.grid(row=1, column=1, sticky=TK.EW)
        WIDGETS['ParameterName'] = B
    else:
        B = TTK.Combobox(F, textvariable=VARS[0],
                         values=freeParams, state='readonly')
        B.bind("<<ComboboxSelected>>", setParameterName)
        B.grid(sticky=TK.EW)
        B.bind('<Return>', setParameterName)
        F.grid(row=1, column=1, sticky=TK.EW)
        WIDGETS['ParameterName'] = B

    # - current parameter value -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=10)
    B.grid(row=2, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Parameter value.')
    B.bind('<Return>', setParameterValue)

    # - slider -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  command=setParameterValueWithScale, borderwidth=1, value=scaleValue)
    WIDGETS['valueSlider'] = B
    B.grid(row=3, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Parameter value.')

    # - trigger mode -
    B = TTK.Button(Frame, text="dXdmu", command=switchMode)
    B.grid(row=4, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Switch mode.')
    WIDGETS['Mode'] = B

    # update for first time
    update()

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.NSEW)
    #try: CTK.WIDGETS['StateNoteBook'].add(WIDGETS['frame'], text='tkDriver')
    #except: pass
    #CTK.WIDGETS['StateNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()
    #CTK.WIDGETS['StateNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp():
    return None

#==============================================================================
def saveApp():
    CTK.PREFS['parameterName'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('None')
    CTK.PREFS['parameterName'] = VARS[0].get()
    CTK.savePrefFile()
    setParameterName()

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
    (win, menu, file, tools) = CTK.minimal('tkDriver '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
