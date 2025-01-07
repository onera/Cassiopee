# - tkExtractBC -
"""Extract Bcs in a pyTree."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

# Known BCs
from tkBC import getAllDefinedBC

#==============================================================================
def updateFamilyBCNameList(event=None):
    if CTK.t == []: return
    m = WIDGETS['BC'].children['menu']
    m.delete(0, TK.END)
    varsp = ['-All BC-']+getAllDefinedBC(CTK.t)
    for i in varsp:
        m.add_command(label=i, command=lambda v=VARS[0],l=i:v.set(l))

def updateFamilyBCNameList2(event=None):
    if CTK.t == []: return
    varsp = ['-All BC-']+getAllDefinedBC(CTK.t)
    if 'BC' in WIDGETS: WIDGETS['BC']['values'] = varsp

#==============================================================================
def setBC2Recover():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    selected = ''
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        selected += CTK.t[2][nob][0]+'/'+z[0]+';'
    selected = selected[0:-1]
    VARS[1].set(selected)

#==============================================================================
# Extrait en fonction du type de BC demandee
# IN: t
# IN: VARS[0] state (type de BC demandee)
# OUT: base EXTRACTED_BC ajoute
#==============================================================================
def extract(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    BCtype = VARS[0].get()
    if BCtype not in Internal.KNOWNBCS: BCtype = 'FamilySpecified:'+BCtype

    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'EXTRACTED_BC', 2)
    base = Internal.getNodeFromName1(CTK.t, 'EXTRACTED_BC')
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []: zones = CTK.t
    else:
        zones = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            zones.append(CTK.t[2][nob][2][noz])

    if BCtype != 'FamilySpecified:-All BC-':
        Z = C.extractBCOfType(zones, BCtype, topTree=CTK.t)
        for i in Z: Internal._createChild(i, 'BCType', 'UserDefined_t', BCtype)
    else:
        (Zp, BCNames, BCTypes) = C.getBCs(zones)
        Z = []
        for i in Zp: Z += i
        for i, z in enumerate(Z): Internal._createChild(z, 'BCType', 'UserDefined_t', BCTypes[i])

    nob = C.getNobOfBase(base, CTK.t)
    for i in Z:
        i[0] = C.getZoneName(i[0])
        CTK.add(CTK.t, nob, -1, i)
    #C._fillMissingVariables(CTK.t) # a cause du BC data set
    CTK.TXT.insert('START', 'BCs of type %s extracted.\n'%BCtype)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
# Recover the BCs on zones
#==============================================================================
def recover(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Zone BC to recover
    name = VARS[1].get()
    names = name.split(';')
    BCs = []; BCNames = []; BCTypes = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]:
                    BCs.append(z)
                    r = sname[1].split('/')
                    BCNames.append(r[1])
                    n = Internal.getNodeFromName1(z, 'BCType')
                    if n is not None: BCTypes.append(Internal.getValue(n))
                    else: BCTypes.append('UserDefined')
    CTK.saveTree()
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        zones = Internal.getZones(CTK.t)
    else:
        zones = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
            zones.append(CTK.t[2][nob][2][noz])

    fail = False
    for z in zones:
        zdim = Internal.getZoneDim(z)
        if zdim[3] == 'NGON': C._recoverBCs(z, (BCs, BCNames, BCTypes))
        else: fail = True

    if not fail: CTK.TXT.insert('START', 'BCs recovered.\n')
    else:
        CTK.TXT.insert('START', 'BC recover fails. Need NGON zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CTK.TKTREE.updateApp()

#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkExtractBC  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Extract boundary conditions.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=2)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkExtractBC')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # - 0 - Type de BC -
    V = TK.StringVar(win); V.set('-All BC-'); VARS.append(V)
    if 'tkExtractBCType' in CTK.PREFS:
        V.set(CTK.PREFS['tkExtractBCType'])
    # - 1 - List of zoneBCs to recover
    V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Extract type de BC -
    B = TTK.Button(Frame, text="Extract", command=extract)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Extract specified BC.\nTree is modified.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[0], *(Internal.KNOWNBCS))
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFamilyBCNameList)
        F.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['BC'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[0],
                         values=Internal.KNOWNBCS,
                         state='readonly', width=10)
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateFamilyBCNameList2)
        F.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
        WIDGETS['BC'] = B

    # - Recover BCs -
    B = TTK.Button(Frame, text="Recover", command=recover)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Recover BCs.\nTree is modified.')
    B = TTK.Button(Frame, command=setBC2Recover,
                   image=iconics.PHOTO[8], padx=0, pady=0)
    B.grid(row=1, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set the BC (zones) to recover.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=8)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='BCs to recover.')

#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['BCNoteBook'].add(WIDGETS['frame'], text='tkExtractBC')
    except: pass
    CTK.WIDGETS['BCNoteBook'].select(WIDGETS['frame'])

#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['BCNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
def saveApp():
    CTK.PREFS['tkExtractBCType'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('-All BC-')
    CTK.PREFS['tkExtractBCType'] = VARS[0].get()
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
    (win, menu, file, tools) = CTK.minimal('tkExtractBC'+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
