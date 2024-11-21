# - tkTime -
"""Time machine."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import RigidMotion.PyTree as RM

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def buildWalls(zones):
    Z = C.extractBCOfType(zones, 'BCWall', topTree=CTK.t)
    return Z

#==============================================================================
def stop(event=None):
    CTK.__BUSY__ = False
    CPlot.setState(cursor=0)

#==============================================================================
def setSlider(event=None):
    t0 = CTK.varsFromWidget(VARS[0].get(), type=1)[0]
    tf = CTK.varsFromWidget(VARS[2].get(), type=1)[0]
    st = WIDGETS['slider'].get()
    step = (tf - t0)/100.
    time = st*step+t0
    VARS[1].set(str(time))
    VARS[4].set('Time [%g].'%time)
    #setTime()

#==============================================================================
def setTime(event=None):
    if CTK.t == []: return
    walls = VARS[3].get()
    time = CTK.varsFromWidget(VARS[1].get(), type=1)
    if len(time) != 1:
        CTK.TXT.insert('START', 'Invalid time.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    time = time[0]
    t0 = CTK.varsFromWidget(VARS[0].get(), type=1)[0]
    tf = CTK.varsFromWidget(VARS[2].get(), type=1)[0]
    step = (tf - t0)/100.

    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()
    if walls == '1' and CTK.dt == []:
        zones = Internal.getZones(CTK.t)
        Z = buildWalls(zones)
        CTK.dt = C.newPyTree(['Base']); CTK.dt[2][1][2] += Z

    if walls == '1': temp = RM.evalPosition(CTK.dt, time)
    else: temp = RM.evalPosition(CTK.t, time)

    WIDGETS['slider'].set((time-t0)/step); WIDGETS['slider'].update()
    if time == 0.: CTK.display(CTK.t, mainTree=CTK.MAIN)
    else: CTK.display(temp, mainTree=CTK.TIME)

#==============================================================================
def playForward(event=None):
    if CTK.t == []: return
    walls = VARS[3].get()
    t0 = CTK.varsFromWidget(VARS[0].get(), type=1)[0]
    time = CTK.varsFromWidget(VARS[1].get(), type=1)[0]
    tf = CTK.varsFromWidget(VARS[2].get(), type=1)[0]
    step = (tf - t0)/100.

    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()

    if walls == '1' and CTK.dt == []:
        zones = Internal.getNodesFromType(CTK.t, 'Zone_t')
        Z = buildWalls(zones)
        CTK.dt = C.newPyTree(['Base']); CTK.dt[2][1][2] += Z

    CTK.__BUSY__ = True
    CPlot.setState(cursor=2)
    while time < tf and CTK.__BUSY__:
        if walls == '1': temp = RM.evalPosition(CTK.dt, time)
        else: temp = RM.evalPosition(CTK.t, time)
        CTK.display(temp, mainTree=CTK.TIME)
        time += step; VARS[1].set(str(time))
        WIDGETS['slider'].set((time-t0)/step)
        WIDGETS['time'].update(); WIDGETS['slider'].update()
    CTK.__BUSY__ = False
    CPlot.setState(cursor=0)

#==============================================================================
def playBackward(event=None):
    if CTK.t == []: return
    walls = VARS[3].get()
    t0 = CTK.varsFromWidget(VARS[0].get(), type=1)[0]
    time = CTK.varsFromWidget(VARS[1].get(), type=1)[0]
    tf = CTK.varsFromWidget(VARS[2].get(), type=1)[0]
    step = (tf - t0)/100.

    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()

    if walls == '1' and CTK.dt == []:
        zones = Internal.getNodesFromType(CTK.t, 'Zone_t')
        Z = buildWalls(zones)
        CTK.dt = C.newPyTree(['Base']); CTK.dt[2][1][2] += Z

    CTK.__BUSY__ = True
    CPlot.setState(cursor=2)
    while time > t0 and CTK.__BUSY__:
        if walls == '1':  temp = RM.evalPosition(CTK.dt, time)
        else: temp = RM.evalPosition(CTK.t, time)
        CTK.display(temp, mainTree=CTK.TIME)
        time -= step; VARS[1].set(str(time))
        WIDGETS['slider'].set((time-t0)/step)
        WIDGETS['time'].update(); WIDGETS['slider'].update()
    CTK.__BUSY__ = False
    CPlot.setState(cursor=0)

#==============================================================================
def setWalls(event=None):
    walls = VARS[3].get()
    if walls == '1':
        zones = Internal.getNodesFromType(CTK.t, 'Zone_t')
        Z = buildWalls(zones)
        CTK.dt = C.newPyTree(['Base']); CTK.dt[2][1][2] += Z
    else: CTK.dt = []

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkTime  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage time.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkTime')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- t0 -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    if 'tkTimeT0' in CTK.PREFS: V.set(CTK.PREFS['tkTimeT0'])
    # -1- t -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -2- tf -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkTimeTF' in CTK.PREFS: V.set(CTK.PREFS['tkTimeTF'])
    # -3- Show only BCWall
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkTimeWall' in CTK.PREFS: V.set(CTK.PREFS['tkTimeWall'])
    # -4- time info bulle
    V = TK.StringVar(win); V.set('Time.'); VARS.append(V)

    # - Times -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Starting time.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=0, column=1, sticky=TK.EW)
    B.bind('<Return>', setTime)
    BB = CTK.infoBulle(parent=B, text='Current time.')
    WIDGETS['time'] = B
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Final time.')

    # - Slider -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, command=setSlider, value=0)
    WIDGETS['slider'] = B
    B.grid(row=1, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[4])

    # - Options -
    B = TTK.Checkbutton(Frame, text='Walls', variable=VARS[3], command=setWalls)
    BB = CTK.infoBulle(parent=B, text='Show only BCWalls.')
    B.grid(row=2, column=0, sticky=TK.EW)

    # - Play, stop -
    B = TTK.Button(Frame, text=">", command=playForward)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Play forward.')

    B = TTK.Button(Frame, text="||", command=stop)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Stop playing.')

    B = TTK.Button(Frame, text="<", command=playBackward)
    B.grid(row=3, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Play backward.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MotionNoteBook'].add(WIDGETS['frame'], text='tkTime')
    except: pass
    CTK.WIDGETS['MotionNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['MotionNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkTimeT0'] = VARS[0].get()
    CTK.PREFS['tkTimeTF'] = VARS[2].get()
    CTK.PREFS['tkTimeWall'] = VARS[3].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('0.')
    VARS[2].set('1.')
    VARS[3].set('0')
    CTK.PREFS['tkTimeT0'] = VARS[0].get()
    CTK.PREFS['tkTimeTF'] = VARS[2].get()
    CTK.PREFS['tkTimeWall'] = VARS[3].get()
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
    (win, menu, file, tools) = CTK.minimal('tkTime '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
