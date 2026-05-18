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

#==============================================================================
# set driver and current entity
def setDriver(driver, entity):
    global DRIVER, ENTITY
    DRIVER = driver
    ENTITY = entity
    return None

#==============================================================================
# get the current free params as set in driver
def getParamFromDriver():
    params = {}
    for f in DRIVER.freeParams:
        params[f.name] = DRIVER.scalars[f.name].v
    return params 

#==============================================================================
def setParameterName(event=None):
    # set parameter value
    paramName = VARS[0].get()
    params = getParamFromDriver()
    r = DRIVER.scalars[paramName].range
    value = params[paramName] # current value
    VARS[1].set(str(value))
    alpha = 100./(r[1]-r[0])
    beta = 100. - alpha*r[1]
    scaleValue = beta + alpha*value
    WIDGETS['valueSlider'].set(scaleValue)
    update()
    return None

#==============================================================================
def setParameterValue(event=None):
    paramName = VARS[0].get()
    params = getParamFromDriver()
    r = DRIVER.scalars[paramName].range
    value = CTK.varsFromWidget(VARS[1].get(), type=1)[0]
    alpha = 100./(r[1]-r[0])
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
# update driver, hook, mesh from current parameter
def update(event=None):
    paramName = VARS[0].get()
    paramValue = CTK.varsFromWidget(VARS[1].get(), type=1)[0]
    params = getParamFromDriver()
    params[paramName] = paramValue
    DRIVER.instantiate(params)
    m = ENTITY.Mesh()
    CPlot.display(m)

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
    freeParams = DRIVER.freeParams
    params = getParamFromDriver()
    paramName = list(params.keys())[0]
    r = DRIVER.scalars[paramName].range
    paramValue = params[paramName]
    alpha = 100./(r[1]-r[0])
    beta = 100. - alpha*r[1]
    scaleValue = beta + alpha*paramValue

    # - VARS -
    # -0- Current parameter name -
    V = TK.StringVar(win); V.set(paramName); VARS.append(V)

    # -1- Current parameter value -
    V = TK.StringVar(win); V.set(paramValue); VARS.append(V)

    # - Parameter chooser -
    B = TTK.Label(Frame, text="Parameters")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Current parameter name.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.Entry(F, textvariable=VARS[0], background='White')
        B.grid(sticky=TK.EW)
        F.bind('<Return>', setParameterName)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['ParameterName'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[0],
                         values=freeParams, state='normal')
        B.bind("<<ComboboxSelected>>", setParameterName)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', setParameterName)
        B.bind('<Return>', setParameterName)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['ParameterName'] = B

    # - current parameter value -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=10)
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Parameter value.')
    B.bind('<Return>', setParameterValue)

    # - slider -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  command=setParameterValueWithScale, borderwidth=1, value=scaleValue)
    WIDGETS['valueSlider'] = B
    B.grid(row=2, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Parameter value.')

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
