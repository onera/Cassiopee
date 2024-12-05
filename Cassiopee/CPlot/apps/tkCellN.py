# - tkCellN -
"""CellNatureField visualisation and extraction."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Post.PyTree as P
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Connector.PyTree as X

# local widgets list
WIDGETS = {}; VARS = []
eps = 1.e-12

#==============================================================================
# set chimera info
#==============================================================================
def chimeraInfo():
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    typename = VARS[1].get()
    CTK.saveTree()
    X._chimeraInfo(CTK.t,typename)
    CTK.TXT.insert('START', 'Field %s added.\n'%typename)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    # try:
    #     CTK.t = X.chimeraInfo(CTK.t,typename)
    #     CTK.TXT.insert('START', 'Field %s added.\n'%typename)
    #     CTK.TKTREE.updateApp()
    #     CTK.display(CTK.t)
    # except Exception, e:
    #     Panels.displayErrors([0,str(e)], header='Error: chimeraInfo')
    #     CTK.TXT.insert('START', 'Fail to add chimera variable %s.\n'%typename)
    #     CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# Filters
#==============================================================================
def Filter1(c1, c2, c3, c4, c5, c6, c7, c8):
    if c1 > 1+eps: return 1
    if c1 < -eps: return 1
    if c2 > 1+eps: return 1
    if c2 < -eps: return 1
    if c3 > 1+eps: return 1
    if c3 < -eps: return 1
    if c4 > 1+eps: return 1
    if c4 < -eps: return 1
    if c5 > 1+eps: return 1
    if c5 < -eps: return 1
    if c6 > 1+eps: return 1
    if c6 < -eps: return 1
    if c7 > 1+eps: return 1
    if c7 < -eps: return 1
    if c8 > 1+eps: return 1
    if c8 < -eps: return 1
    return 0

#==============================================================================
# Filter cellN
# Pour les zones actives de t
#==============================================================================
def view():
    if CTK.t == []: return
    type = VARS[0].get()
    if type == 'Mesh': CTK.display(CTK.t); return

    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()

    tp = Internal.appendBaseName2ZoneName(CTK.t, separator=Internal.SEP1)

    active = []
    zones = Internal.getZones(tp)
    for z in CTK.__MAINACTIVEZONES__: active.append(zones[z])

    Z = None
    if type == 'cf>1':
        Z = P.selectCells(active, Filter1, ['interpCoefs1', 'interpCoefs2',
                                            'interpCoefs3', 'interpCoefs4',
                                            'interpCoefs5', 'interpCoefs6',
                                            'interpCoefs7', 'interpCoefs8'])
    elif type == 'cellN=-99999':
        Z = selectWithFormula(active, '{cellN} == -99999')
    elif type == 'cellN=1':
        Z = selectWithFormula(active, '{cellN} == 1')
    elif type == 'cellN=0':
        Z = selectWithFormula(active, '{cellN} == 0')
    elif type == 'cellN=2':
        Z = selectWithFormula(active, '{cellN} == 2')
    elif type == 'cellN<0':
        Z = selectWithFormula(active, '{cellN} < 0')
    elif type == '0<cellN<1':
        Z = selectWithFormula(active, '({cellN}>0) & ({cellN}<1)')
    elif type == 'Orphan points':
        Z = X.extractChimeraInfo(active,'orphan')
    elif type == 'Extrapolated points':
        Z = X.extractChimeraInfo(active,'extrapolated')

    if Z is not None:
        CTK.TXT.insert('START', 'Filter '+type+' displayed.\n')
        CTK.dt = C.newPyTree(['Base'])
        CTK.dt[2][1][2] += Z
        CTK.display(CTK.dt, mainTree=CTK.CELLN)
    else:
        CTK.TXT.insert('START', 'Nothing to display.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# La variable var existe t elle dans la premiere zone de l'arbre?
# Retourne 0: non
# Retourne 1: oui, en noeuds
# Retourne 2: oui, en centres
#==============================================================================
def findVar(var):
    v = C.getVarNames(CTK.t)
    if len(v) > 0: vars = v[0]
    else: vars = []
    for i in vars:
        if var == i: return 1
        if 'centers:'+var == i: return 2
    return 0

#==============================================================================
# Cree un champ __tag__ dans chaque zone, selectionne, efface __tag__.
#==============================================================================
def selectTag(zones, Filter, vars):
    if len(vars) == 0: return None
    res = findVar(vars[0])
    if res == 0: return None
    if res == 1: # tag en noeuds
        Z = C.initVars(zones, '__tag__', Filter, vars)
        Z = P.selectCells2(Z, '__tag__')
        C._rmVars(Z, '__tag__')
    else: # tag en centres
        Z = C.initVars(zones, 'centers:__tag__', Filter, vars)
        Z = P.selectCells2(Z, 'centers:__tag__')
        C._rmVars(Z, 'centers:__tag__')
    for z in Z: z[0] = C.getZoneName(z[0])
    return Z

#==============================================================================
# Selectionne avec une formule en cellN
#==============================================================================
def selectWithFormula(zones, formula):
    res = findVar('cellN')
    if res == 0:
        res = findVar('cellnf')
        if res == 0: return None
        if res == 1: formula = formula.replace('cellN', 'cellnf')
        if res == 2: formula = formula.replace('cellN', 'centers:cellnf')
    elif res == 2:
        formula = formula.replace('cellN', 'centers:cellN')

    if res == 1: # cellN en noeuds
        Z = P.selectCells(zones, formula)
    else: # cellN en centres
        Z = P.selectCells(zones, formula)
    for z in Z: z[0] = C.getZoneName(z[0])
    return Z

#==============================================================================
def extract():
    if CTK.t == []: return
    type = VARS[0].get()

    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()

    active = []
    zones = Internal.getZones(CTK.t)
    for z in CTK.__MAINACTIVEZONES__: active.append(CTK.t[2][CTK.Nb[z]+1][2][CTK.Nz[z]])

    Z = None
    if type == 'cellN=-99999':
        Z = selectWithFormula(active, '{cellN} == -99999')
    elif type == 'cellN=1':
        Z = selectWithFormula(active, '{cellN} == 1')
    elif type == 'cellN=0':
        Z = selectWithFormula(active, '{cellN} == 0')
    elif type == 'cellN=2':
        Z = selectWithFormula(active, '{cellN} == 2')
    elif type == 'cellN<0':
        Z = selectWithFormula(active, '{cellN}<0')
    elif type == '0<cellN<1':
        Z = selectWithFormula(active, '({cellN}>0) & ({cellN}<1)')
    elif type == 'Interpolated points':
        Z = X.extractChimeraInfo(zones,type='interpolated',loc='centers')
        if Z == []: Z = None
    elif type == 'Extrapolated points':
        Z = X.extractChimeraInfo(zones,type='extrapolated',loc='centers')
        if Z == []: Z = None
    elif type == 'Orphan points':
        Z = X.extractChimeraInfo(zones,type='orphan',loc='centers')
        if Z == []: Z = None
    elif type == 'cf>1':
        Z = X.extractChimeraInfo(zones,type='cf>1',loc='centers')
        if Z == []: Z = None

    if Z is not None:
        CTK.TXT.insert('START', 'Filter '+type+' extracted.\n')
        C._addBase2PyTree(CTK.t, 'EXTRACT')
        b = Internal.getNodesFromName1(CTK.t, 'EXTRACT')
        base = b[0]
        base[2] += Z
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        #C._fillMissingVariables(CTK.t)
        CTK.TKTREE.updateApp()
        CTK.display(CTK.t)
    else:
        CTK.TXT.insert('START', 'Nothing extracted.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCellN  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Chimera cellN analysis.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkCellN')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- cellN filter -
    norow = 0
    V = TK.StringVar(win); V.set('cellN=0'); VARS.append(V)
    if 'tkCellNFilter' in CTK.PREFS: V.set(CTK.PREFS['tkCellNFilter'])

    # Filter
    B = TTK.OptionMenu(Frame, VARS[0], 'Mesh', 'cellN=0', 'cellN=-99999',
                       'cellN=2', 'cellN<0', '0<cellN<1', 'cellN=1', 'cf>1', 'Orphan points',
                       'Extrapolated points','Interpolated points')
    B.grid(row=norow, column=0, columnspan=2, sticky=TK.EW)

    # - View cellN -
    norow += 1
    B = TTK.Button(Frame, text="View", command=view)
    B.grid(row=norow, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='View the location of specified cellN.\nTree is NOT modified.')

    # - Extract cellN -
    B = TTK.Button(Frame, text="Extract", command=extract)
    B.grid(row=norow, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract the location of specified cellN.\nTree is modified.')

    # - chimeraInfo
    # -1- chimeraInfo type
    norow += 1
    V = TK.StringVar(win); V.set('interpolated'); VARS.append(V)
    if 'tkChimeraInfoType' in CTK.PREFS: V.set(CTK.PREFS['tkChimeraInfoType'])

    # Filter
    B = TTK.Button(Frame, text="Chimera info", command=chimeraInfo)
    B.grid(row=norow, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[1], 'interpolated', 'extrapolated', 'orphan',
                       'cellRatio','donorAspect')
    B.grid(row=norow, column=1, sticky=TK.EW)

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['VisuNoteBook'].add(WIDGETS['frame'], text='tkCellN')
    except: pass
    CTK.WIDGETS['VisuNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['VisuNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkCellNFilter'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('cellN=0')
    CTK.PREFS['tkCellNFilter'] = VARS[0].get()
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
    (win, menu, file, tools) = CTK.minimal('tkCellN '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
