# - tkDb -
"""GUI for Roms/DB."""
import tkinter as TK
import CPlot.Ttk as TTK
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.PyTree as C
import numpy

# local widgets list
WIDGETS = {}; VARS = []

# DB
DB = None
# dict of current parameter/values
PARAMS = {}
# dict of min values
RANGEMIN = {}
# dict of max values
RANGEMAX = {}
# sorted points of parameters
POINTS = {}

#==============================================================================
# set data base
def setDb(db):
    global DB
    DB = db
    return None

#==============================================================================
# Find parameters and minmax from db
def initParametersFromDb():
    global PARAMS, RANGEMIN, RANGEMAX
    q = DB.query() # query all
    points = DB.fetchPoints(q) # may be huge in memory
    npoints = len(points)
    if npoints == 0: return
    p = points[0]
    for k in p: RANGEMIN[k] = +1.e6
    for k in p: RANGEMAX[k] = -1.e6
    for k in p: POINTS[k] = numpy.zeros( (npoints), dtype=numpy.float64 )
    for c, p in enumerate(points):
        for k in p: RANGEMIN[k] = min(RANGEMIN[k], p[k])
        for k in p: RANGEMAX[k] = max(RANGEMAX[k], p[k])
        for k in p: POINTS[k][c] = p[k]
    for k in points[0]: PARAMS[k] = RANGEMIN[k]
    for k in points[0]: POINTS[k].sort()
    print("RANGEMIN=", RANGEMIN)
    print("RANGEMAX=", RANGEMAX)
    print("PARAMS=", PARAMS)
    return None

#==============================================================================
# update t from current parameters, display
def update(event=None):
    paramName = VARS[0].get()
    paramValue = CTK.varsFromWidget(VARS[1].get(), type=1)[0]
    PARAMS[paramName] = paramValue
    print("PARAMS=", PARAMS)
    q = DB.query(PARAMS)
    t = DB.fetchTree(q)
    CPlot.display(t, mode='Mesh', bgColor=1)

#==============================================================================
def setParameterValueWithScale(event=None):
    val = WIDGETS['valueSlider'].get()
    paramName = VARS[0].get()
    vmax = RANGEMAX[paramName]
    vmin = RANGEMIN[paramName]
    value = vmin + val*(vmax-vmin)/100.
    idx = numpy.searchsorted(POINTS[paramName], value)
    value = POINTS[paramName][idx] # nearest
    VARS[1].set(str(value))
    update()
    return None

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkDb  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
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

    # Get default values
    initParametersFromDb()
    params = list(PARAMS.keys())

    # - VARS -
    # -0- Current parameter name -
    V = TK.StringVar(win); V.set(params[0]); VARS.append(V)

    # -1- Current parameter value -
    V = TK.StringVar(win); V.set(PARAMS[params[0]]); VARS.append(V)

    # - Parameter chooser -
    B = TTK.Label(Frame, text="Parameter")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Current parameter name.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.Entry(F, textvariable=VARS[0], background='White')
        B.grid(sticky=TK.EW)
        F.bind('<Return>', update)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['ParameterName'] = B
    else:
        B = TTK.Combobox(F, textvariable=VARS[0],
                         values=params, state='readonly')
        B.bind("<<ComboboxSelected>>", update)
        B.grid(sticky=TK.EW)
        B.bind('<Return>', update)
        F.grid(row=0, column=1, sticky=TK.EW)
        WIDGETS['ParameterName'] = B

    # - current parameter value -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=10)
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Parameter value.')
    B.bind('<Return>', update)

    # - slider -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1) #, value=scaleValue)
    B.bind("<ButtonRelease-1>", setParameterValueWithScale)
    WIDGETS['valueSlider'] = B
    B.grid(row=2, column=0, columnspan=2, sticky=TK.EW)
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
    #setParameterName()

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
