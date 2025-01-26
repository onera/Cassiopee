# - CPlot control app  -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import KCore

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def oneovern(event=None):
    type = VARS[1].get()
    if type == 'All points': CTK.__ONEOVERN__ = 1
    elif type == 'One over 2': CTK.__ONEOVERN__ = 2
    elif type == 'One over 3': CTK.__ONEOVERN__ = 3
    elif type == 'One over 4': CTK.__ONEOVERN__ = 4
    elif type == 'Exclusive': CTK.__ONEOVERN__ = 0
    if CTK.t == []: return
    else: CTK.display(CTK.t)

#==============================================================================
def setViewMode(v, l):
    v.set(l)
    type = v.get()
    if type == 'All fields':
        CTK.__FIELD__ = '__all__'
        if CTK.t != []: CTK.display(CTK.t)
    elif type == 'No field':
        CTK.__FIELD__ = '__none__'
        if CTK.t != []: CPlot.setMode(0); CTK.display(CTK.t)
    else:
        CTK.__FIELD__ = type
        if CTK.t != []: CTK.display(CTK.t)

def setViewMode2(event=None):
    type = VARS[0].get()
    if type == 'All fields':
        CTK.__FIELD__ = '__all__'
        if CTK.t != []: CTK.display(CTK.t)
    elif type == 'No field':
        CTK.__FIELD__ = '__none__'
        if CTK.t != []: CPlot.setMode(0); CTK.display(CTK.t)
    else:
        CTK.__FIELD__ = type
        if CTK.t != []: CTK.display(CTK.t)

#==============================================================================
def updateVarNameList(event=None):
    m = WIDGETS['fields'].children['menu']
    m.delete(0, TK.END)
    for i in ['All fields', 'No field']:
        m.add_command(label=i, command=lambda v=VARS[0],l=i:setViewMode(v,l))
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        vars = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
    if len(vars) == 0: return
    for i in vars[0]:
        m.add_command(label=i, command=lambda v=VARS[0],l=i:setViewMode(v,l))

#==============================================================================
def updateVarNameList2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        vars = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
    if len(vars) == 0: return
    vars = ['All fields', 'No field'] + vars[0]
    if 'fields' in WIDGETS:
        WIDGETS['fields']['values'] = vars

#==============================================================================
def getThreads():
    nt = KCore.kcore.getOmpMaxThreads()
    VARS[2].set(str(nt))
    WIDGETS['threads'].update()
    return

#==============================================================================
def getMaxProcs():
    try:
        import multiprocessing
        return str(multiprocessing.cpu_count())
    except: return '?'

#==============================================================================
def setThreads(event=None):
    nt = VARS[2].get()
    try:
        nti = int(nt)
        KCore.kcore.setOmpMaxThreads(nti)
        CTK.TXT.insert('START', 'Num threads set to %d.\n'%nti)
    except:
        CTK.TXT.insert('START', 'Bad thread number.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
    return

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkPerfo  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Improve performance of Cassiopee.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=4)
    Frame.columnconfigure(2, weight=1)

    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkPerfo')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Displayed field -
    V = TK.StringVar(win); V.set('All fields'); VARS.append(V)
    # -1- Oneovern -
    V = TK.StringVar(win); V.set('All points'); VARS.append(V)
    if 'tkPerfoPoints' in CTK.PREFS: V.set(CTK.PREFS['tkPerfoPoints'])
    # -2- Threads
    V = TK.StringVar(win); V.set('1'); VARS.append(V)

    # - Field transmission -
    B = TTK.Label(Frame, text="Display")
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Describe what part of the tree is displayed.\nCan fasten the display. In memory tree is not modified.')

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[0], 'All fields', 'No field')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList)
        F.grid(row=0, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Display only certain fields.')
        WIDGETS['fields'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[0],
                         values=['All fields', 'No field'],
                         state='readonly', width=10)
        B.grid(sticky=TK.EW)
        B.bind('<<ComboboxSelected>>', setViewMode2)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=0, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Display only certain fields.')
        WIDGETS['fields'] = B

    # - Point transmission -
    B = TTK.OptionMenu(Frame, VARS[1], 'All points', 'One over 2', 'One over 3',
                       'One over 4', 'Exclusive', command=oneovern)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Display less points.')

    # - OMP threads -
    B = TTK.Label(Frame, text="Threads")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='The number of CPU cores used by Cassiopee.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[2], background='White', width=3)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of threads\n(init value is OMP_NUM_THREADS).')
    B.bind('<Return>', setThreads)
    WIDGETS['threads'] = B
    mxProcs = getMaxProcs()
    B = TK.Label(F, text='/'+mxProcs)
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)

    B = TTK.Button(Frame, text="Set", command=setThreads)
    B.grid(row=1, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the number of threads.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    getThreads()
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['StateNoteBook'].add(WIDGETS['frame'], text='tkPerfo')
    except: pass
    CTK.WIDGETS['StateNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['StateNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkPerfoPoints'] = VARS[1].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[1].set('All points')
    CTK.PREFS['tkPerfoPoints'] = VARS[1].get()
    CTK.savePrefFile()
    oneovern()

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
    (win, menu, file, tools) = CTK.minimal('tkPerfo '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
