# - global preferences -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Met a jour PREFS a partir des donnees de l'interface
#==============================================================================
def buildPrefs():
    v = VARS[1].get(); va = '0'
    if v == 'Active': va = '1'
    elif v == 'Inactive': va = '0'
    CTK.PREFS['undo'] = va
    v = VARS[2].get(); va = '0'
    if v == 'Active': va = '1'
    elif v == 'Inactive': va = '0'
    CTK.PREFS['displayInfo'] = va
    v = VARS[3].get(); va = '0'
    if v == 'Black': va = '0'
    elif v == 'White': va = '1'
    elif v == 'Grey': va = '2'
    elif v == 'Yellow': va = '3'
    elif v == 'Blue': va = '4'
    elif v == 'Red': va = '5'
    elif v == 'Paper1': va = '6'
    elif v == 'Paper2': va = '7'
    elif v == 'Arch': va = '8'
    elif v == 'Jarvis': va = '9'
    elif v == 'Onera': va = '10'
    elif v == 'Clouds': va = '11'
    elif v == 'Blueprint': va = '12'
    elif v == 'Custom': va = '13'
    CTK.PREFS['bgColor'] = va
    v = VARS[6].get()
    CTK.PREFS['envmap'] = v
    v = VARS[7].get()
    CTK.PREFS['selectionStyle'] = v
    v = VARS[8].get()
    CTK.PREFS['exportResolution'] = v
    v = VARS[9].get()
    CTK.PREFS['guitheme'] = v
    autoapps = VARS[4].get()
    CTK.PREFS['auto'] = autoapps
    return

#==============================================================================
def setPrefs(event=None):
    buildPrefs()
    CTK.setPrefs()

def setUndo(event=None):
    v = VARS[1].get()
    if v == 'Active': CTK.__UNDO__ = True
    else: CTK.__UNDO__ = False
    va = '0'
    if v == 'Active': va = '1'
    elif v == 'Inactive': va = '0'
    CTK.PREFS['undo'] = va

def setBgColor(event=None):
    v = VARS[3].get(); va = 0
    if v == 'Black': va = 0
    elif v == 'White': va = 1
    elif v == 'Grey': va = 2
    elif v == 'Yellow': va = 3
    elif v == 'Blue': va = 4
    elif v == 'Red': va = 5
    elif v == 'Paper1': va = 6
    elif v == 'Paper2': va = 7
    elif v == 'Arch': va = 8
    elif v == 'Jarvis': va = 9
    elif v == 'Onera': va = 10
    elif v == 'Clouds': va = 11
    elif v == 'Blueprint': va = 12
    elif v == 'Custom':
        try: import tkinter.filedialog as tkFileDialog
        except: import tkFileDialog
        if 'backgroundFile' in CTK.PREFS: initFile = CTK.PREFS['backgroundFile']
        else: initFile = 'paperBackground1.png'
        files = tkFileDialog.askopenfilenames(
            filetypes=[('png image file', '*.png')], initialfile=initFile, multiple=0)
        if files == '' or files is None or files == (): # user cancel
            return
        files = CTK.fixFileString__(files, initFile)
        CPlot.setState(backgroundFile=files[0])
        CTK.PREFS['backgroundFile'] = files[0]
        va = 13

    CPlot.setState(bgColor=va)
    CTK.PREFS['bgColor'] = str(va)

def setDisplayInfo(event=None):
    v = VARS[2].get(); va = 0
    if v == 'Active': va = 1
    elif v == 'Inactive': va = 0
    CPlot.setState(displayInfo=va)
    CPlot.setState(displayBB=va)
    CTK.PREFS['displayInfo'] = str(va)

def setSelectionStyle(event=None):
    v = VARS[7].get(); style = 0
    if v == 'Blue': style = 0
    elif v == 'Alpha': style = 1
    CPlot.setState(selectionStyle=style)
    CTK.PREFS['selectionStyle'] = v

def setOnDragStyle(event=None):
    v = VARS[10].get(); ondrag = 0
    if v == 'None': ondrag = 0
    elif v == 'BBox': ondrag = 1
    CPlot.setState(simplifyOnDrag=ondrag)
    CTK.PREFS['simplifyOnDrag'] = v

def setEnvmap(event=None):
    val = VARS[6].get()
    CPlot.setState(envmap=val)
    CTK.PREFS['envmap'] = val

def setExportRes(event=None):
    CTK.PREFS['exportResolution'] = VARS[8].get()

def setTheme(event=None):
    theme = VARS[9].get()
    TTK.setTheme(theme)
    CTK.PREFS['guitheme'] = theme

def updateThemeList(event=None):
    m = WIDGETS['guitheme'].children['menu']
    m.delete(0, TK.END)
    allvars = TTK.getAvailableThemes()
    for i in allvars:
        m.add_command(label=i, command=lambda v=VARS[9],l=i:TTK.setTheme(l))

def updateThemeList2(event=None):
    allvars = TTK.getAvailableThemes()
    WIDGETS['guitheme']['values'] = allvars

#==============================================================================
def savePrefs():
    buildPrefs()
    CTK.savePrefFile()
    CTK.TXT.insert('START', 'Prefs file saved.\n')

#==============================================================================
def getOpenedApps():
    import tkCassiopee
    allapps = tkCassiopee.ALLAPPS
    opened = []
    for i in allapps:
        if (i != '---'):
            app = __import__(i)
            info = app.WIDGETS['frame'].grid_info()
            if info != {}: opened.append(i)

    st = ''
    for i in opened: st += i+';'
    st = st[:-1]
    VARS[4].set(st)

def getClassicApps():
    st = 'tkMeshInfo;tkBlock;tkBC;tkView;tkSlice'
    VARS[4].set(st)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkPrefs  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Set Cassiopee preferences.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    #FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkPrefs')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Default mesh display style -
    V = TK.StringVar(win); V.set('Monocolor wires+solid'); VARS.append(V)
    # -1- Undo
    V = TK.StringVar(win); V.set('Active'); VARS.append(V)
    # -2- DisplayInfo
    V = TK.StringVar(win); V.set('Active'); VARS.append(V)
    # -3- bgColor
    V = TK.StringVar(win); V.set('Black'); VARS.append(V)
    # -4- auto open apps
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -5- Default solid display style -
    V = TK.StringVar(win); V.set('Monocolor/1-side'); VARS.append(V)
    # -6- Default envmap
    V = TK.StringVar(win); V.set('windtunnel.png'); VARS.append(V)
    # -7- Selection style
    V = TK.StringVar(win); V.set('Blue'); VARS.append(V)
    # -8- Export resolution
    V = TK.StringVar(win); V.set('1920x1080'); VARS.append(V)
    # -9- GUI Theme
    V = TK.StringVar(win); V.set('default'); VARS.append(V)
    # -10- simplify on drag
    V = TK.StringVar(win); V.set('None'); VARS.append(V)

    # Init VARS par le fichier de preferences
    CTK.loadPrefFile()
    for i in CTK.PREFS:
        k1 = CTK.PREFS[i]
        if i == 'undo':
            if k1 == '1': VARS[1].set('Active')
            elif k1 == '0': VARS[1].set('Inactive')
        elif i == 'displayInfo':
            if k1 == '1': VARS[2].set('Active')
            elif k1 == '0': VARS[2].set('Inactive')
        elif i == 'bgColor':
            if k1 == '0': VARS[3].set('Black')
            elif k1 == '1': VARS[3].set('White')
            elif k1 == '2': VARS[3].set('Grey')
            elif k1 == '3': VARS[3].set('Yellow')
            elif k1 == '4': VARS[3].set('Blue')
            elif k1 == '5': VARS[3].set('Red')
            elif k1 == '6': VARS[3].set('Paper1')
            elif k1 == '7': VARS[3].set('Paper2')
            elif k1 == '8': VARS[3].set('Arch')
            elif k1 == '9': VARS[3].set('Jarvis')
            elif k1 == '10': VARS[3].set('Onera')
            elif k1 == '11': VARS[3].set('Clouds')
            elif k1 == '12': VARS[3].set('Blueprint')
            elif k1 == '13': VARS[3].set('Custom')

        elif i == 'auto': VARS[4].set(k1)
        elif i == 'envmap': VARS[6].set(k1)
        elif i == 'selectionStyle': VARS[7].set(k1)
        elif i == 'simplifyOnDrag': VARS[10].set(k1)
        elif i == 'exportResolution': VARS[8].set(k1)
        elif i == 'guitheme': VARS[9].set(k1)

    r = 0
    # - gui theme -
    B = TTK.Label(Frame, text="GUI Theme")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Style for GUI (buttons,...).')

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TTK.OptionMenu(F, VARS[9], 'default')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateThemeList)
        F.grid(row=r, column=1, sticky=TK.EW)
        WIDGETS['guitheme'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[9],
                         values=[], state='readonly', width=10)
        B.bind('<<ComboboxSelected>>', setTheme)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateThemeList2)
        F.grid(row=r, column=1, sticky=TK.EW)
        WIDGETS['guitheme'] = B
    r += 1

    # - Default display info style -
    B = TTK.Label(Frame, text="Display Info")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Default display info style.')
    B = TTK.OptionMenu(Frame, VARS[2], 'Active', 'Inactive', command=setDisplayInfo)
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1

    # - Default background color -
    B = TTK.Label(Frame, text="Background")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Default background color or image.')
    B = TTK.OptionMenu(Frame, VARS[3], 'Black', 'White', 'Grey', 'Yellow',
                       'Blue', 'Red', 'Paper1', 'Paper2',
                       'Arch', 'Jarvis', 'Onera', 'Clouds',
                       'Blueprint', 'Custom', command=setBgColor)
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1

    # - Undo -
    B = TTK.Label(Frame, text="Undo")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Undo is active.')
    B = TTK.OptionMenu(Frame, VARS[1], 'Active', 'Inactive', command=setUndo)
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1

    # - selection style -
    B = TTK.Label(Frame, text="Selection")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Selection style.')
    B = TTK.OptionMenu(Frame, VARS[7], 'Blue', 'Alpha', command=setSelectionStyle)
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1

    # - simplify on drag mode -
    B = TTK.Label(Frame, text="On drag")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Kind of simplification when dragging.')
    B = TTK.OptionMenu(Frame, VARS[10], 'None', 'BBox', command=setOnDragStyle)
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1

    # - Envmap -
    B = TTK.Label(Frame, text="Envmap")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Environment map (chrome/glass render).')
    B = TTK.OptionMenu(Frame, VARS[6], 'windtunnel.png', 'sky.png', 'city.png', 'forest.png', 'house.png', command=setEnvmap)
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1

    # - export resolution -
    B = TTK.Label(Frame, text="Export")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Resolution for image export.')
    B = TTK.Entry(Frame, textvariable=VARS[8], background='White', width=5)
    B.bind('<Return>', setExportRes)
    BB = CTK.infoBulle(parent=B, text='Resolution for image export.')
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1


    # - Auto open apps -
    #B = TK.Button(Frame, text='Opened apps', command=getOpenedApps)
    #BB = CTK.infoBulle(parent=B, text='Keep currently opened apps for next restart.')
    #B.grid(row=r, column=0, sticky=TK.EW)
    #B = TK.Button(Frame, text='Classic apps', command=getClassicApps)
    #BB = CTK.infoBulle(parent=B, text='Keep a classic set of apps for next restart.')
    #B.grid(row=r, column=1, sticky=TK.EW)
    #r += 1
    #B = TK.Label(Frame, text="Auto open")
    #B.grid(row=r, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Open automatically these apps for next restart (tkTree;tkBC;...)')
    #B = TK.Entry(Frame, textvariable=VARS[4], background='White')
    #BB = CTK.infoBulle(parent=B, text='Apps opened for next restart (tkTree;tkBC;...)')
    #B.grid(row=r, column=1, sticky=TK.EW)
    #r += 1

    # - Set prefs -
    #B = TK.Button(Frame, text="Set prefs", command=setPrefs)
    #B.grid(row=r, column=0, columnspan=2, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Set prefs in Cassiopee.')
    #WIDGETS.append(B); r += 1

    # - Save prefs -
    B = TTK.Button(Frame, text="Save prefs", command=savePrefs)
    B.grid(row=r, column=0, columnspan=2, sticky=TK.EW)
    B = CTK.infoBulle(parent=B, text='Save pref file (.cassiopee).')
    r += 1

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['StateNoteBook'].add(WIDGETS['frame'], text='tkPrefs')
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
    savePrefs()

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
    (win, menu, file, tools) = CTK.minimal('tkPrefs '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
