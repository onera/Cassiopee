# - tkDistributor -
"""Block distribution over processors."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Distributor2.PyTree as D
import Transform.PyTree as T
import Connector.PyTree as X
import Converter.Internal as Internal
import numpy, math
# local widgets list
WIDGETS = {}; VARS = []
STATS = {}

#==============================================================================
# Distribute
#==============================================================================
def distribute(event=None):
    global STATS
    if CTK.t == []: return
    try: NProc = int(VARS[0].get())
    except:
        CTK.TXT.insert('START', 'distribute: NProc is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    try: comSpeed = float(VARS[1].get())
    except:
        CTK.TXT.insert('START', 'distribute: ComSpeed is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    algo = VARS[2].get()
    useCom =  VARS[3].get()

    CTK.saveTree()
    CTK.t, STATS = D.distribute(CTK.t, NProc, perfo=(1.,0.,comSpeed),
                                useCom=useCom, algorithm=algo)
    CTK.TXT.insert('START', 'Blocks distributed.\n')
    CTK.TKTREE.updateApp()
    updateStats()

#==============================================================================
# Split and distribute
#==============================================================================
def splitAndDistribute(event=None):
    global STATS
    if CTK.t == []: return
    try: NProc = int(VARS[0].get())
    except:
        CTK.TXT.insert('START', 'distribute: NProc is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    try: comSpeed = float(VARS[1].get())
    except:
        CTK.TXT.insert('START', 'distribute: ComSpeed is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    algo = VARS[2].get()
    useCom = VARS[3].get()
    level = int(VARS[5].get())
    CTK.saveTree()

    try:
        #CTK.t = T.splitNParts(CTK.t, NProc, multigrid=level)
        CTK.t = T.splitSize(CTK.t, N=0, multigrid=level, R=NProc, type=2)
        # no need to check inconsistant match (they have been deleted)
        node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
        if node is not None:
            ndim = Internal.getValue(node)
            # Manque le reglage de la tol
        else:
            CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
            CTK.TXT.insert('START', 'Warning: ', 'Warning')
            ndim = 3
        CTK.t = X.connectMatch(CTK.t, dim=ndim)

        CTK.display(CTK.t)
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: distribute/split')
        CTK.TXT.insert('START', 'splitSize fails for at least one zone.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)

    CTK.t, STATS = D.distribute(CTK.t, NProc, perfo=(1.,0.,comSpeed),
                                useCom=useCom, algorithm=algo)
    CTK.TXT.insert('START', 'Blocks split and distributed.\n')
    CTK.TKTREE.updateApp()
    updateStats()

#==============================================================================
# setProc to selection
#==============================================================================
def setProc(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    proc = int(VARS[4].get())
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        nodes = Internal.getNodesFromName1(z, '.Solver#Param')
        if nodes != []: param = nodes[0]
        else:
            param = ['.Solver#Param', None, [], 'UserDefinedData_t']
            z[2].append(param)
        v = numpy.zeros((1,1), dtype=Internal.E_NpyInt); v[0,0] = proc
        nodes = Internal.getNodesFromName(param, 'proc')
        if nodes != []:
            a = nodes[0]; a[1] = v
        else:
            a = ['proc', v, [], 'DataArray_t']
            param[2].append(a)

    CTK.TXT.insert('START', 'Proc set in selection.\n')
    CTK.TKTREE.updateApp()
    computeStats(); updateStats()

#==============================================================================
# setProcField
# Note: je pense que cette fonctionalite est obsolete et avantageusement
# remplacee par tkFilter
#==============================================================================
def setProcField():
    if CTK.t == []: return
    CTK.saveTree()
    zones = Internal.getZones(CTK.t)
    for z in zones:
        param = Internal.getNodesFromName(z, '.Solver#Param')
        if param != []:
            nodes = Internal.getNodesFromName1(param[0], 'proc')
            (r,c) = Internal.getParentOfNode(CTK.t, z)
            if nodes != []:
                value = nodes[0][1][0,0]
                C._initVars(z, 'proc', value)
            else: C._initVars(z, 'proc', -1.)
        else: C._initVars(z, 'proc', -1.)
    CTK.TXT.insert('START', 'Field proc set.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# Compute stats
# IN: NProc + arbre (proc)
# OUT: STATS
#==============================================================================
def computeStats():
    if CTK.t == []: return
    zones = Internal.getZones(CTK.t)

    try: NProc = int(VARS[0].get())
    except:
        CTK.TXT.insert('START', 'stats: NProc is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    a = numpy.zeros((NProc), dtype=Internal.E_NpyInt)
    m = 0; ntot = 0
    for z in zones:
        param = Internal.getNodesFromName1(z, '.Solver#Param')
        if param != []:
            proc = Internal.getNodesFromName1(param[0], 'proc')
            if proc != []:
                value = proc[0][1][0,0]
                dim = Internal.getZoneDim(z)
                if dim[0] == 'Structured': size = dim[1]*dim[2]*dim[3]
                else: size = dim[1]
                a[value] += size
                ntot += size
    m = ntot*1. / NProc
    varRMS = 0.; varMin = 1.e6; varMax = 0.
    for i in a:
        v = abs(i-m)
        varMin = min(varMin, v)
        varMax = max(varMax, v)
        varRMS += v*v

    varMin = varMin / m
    varMax = varMax / m
    varRMS = math.sqrt(varRMS) / (NProc*m)

    STATS['meanPtsPerProc'] = m
    STATS['varRMS'] = varRMS
    STATS['varMin'] = varMin
    STATS['varMax'] = varMax
    STATS['nptsCom'] = 0.
    STATS['comRatio'] = 0.
    return

#==============================================================================
# Update le canvas des stats
# IN: NProc (VARS[0])
# IN: STATS
#==============================================================================
def updateStats():
    if CTK.t == []: return
    # Update canvas
    c = WIDGETS['canvas']
    c.delete(TK.ALL)
    width = int(c.cget("width"))
    height = int(c.cget("height"))
    c.create_line(0, height/2., width, height/2., fill='red')

    try: NProc = int(VARS[0].get())
    except:
        CTK.TXT.insert('START', 'stats: NProc is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    # Calcul du nombre de pts par proc
    zones = Internal.getZones(CTK.t)
    a = numpy.zeros((NProc), dtype=Internal.E_NpyInt)
    for z in zones:
        param = Internal.getNodesFromName1(z, '.Solver#Param')
        if param != []:
            proc = Internal.getNodesFromName1(param[0], 'proc')
            if proc != []:
                value = proc[0][1][0,0]
                dim = Internal.getZoneDim(z)
                if dim[0] == 'Structured':
                    a[value] += dim[1]*dim[2]*dim[3]
                else:
                    a[value] += dim[1]

    m = STATS['meanPtsPerProc']
    fmin = 1.e10; fmax = 0
    for i in a:
        fmin = min(i, fmin); fmax = max(i, fmax)

    alpha = min(1./abs(fmax-m+1.e-6), 1./abs(fmin-m+1.e-6))
    alpha = min(alpha, 1./m)

    barWidth = width*1. / (NProc*1.)
    for i in range(NProc):
        v = -alpha*(a[i] - m*1.)
        if a[i] == 0: fillColor = 'yellow'
        elif i%2 == 0: fillColor = 'blue'
        else: fillColor = 'red'

        c.create_rectangle(i*barWidth, height/2., (i+1)*barWidth,
                           v*height/2.+height/2.,
                           fill=fillColor)

    varRMS = int(STATS['varRMS']*10000)/100.
    varMin = int(STATS['varMin']*10000)/100.
    varMax = int(STATS['varMax']*10000)/100.
    nptsCom = STATS['nptsCom']
    comRatio = int(STATS['comRatio']*10000)/100.
    CTK.TXT.insert('START', 'stats: mean='+str(m)+
                   ' pts, varMax='+str(varMax)+'%, nptsCom='+str(nptsCom)+
                   ' npts, ratio='+str(comRatio)+'%.\n')

    # Update info bulle
    b = WIDGETS['bulle']
    stats = 'meanPtsPerProc='+str(STATS['meanPtsPerProc'])+'.\n'
    stats += 'varMax='+str(varMax)+'%, varMin='+\
             str(varMin)+'%, varRMS='+str(varRMS)+'%.\n'
    stats += 'nptsCom='+str(STATS['nptsCom'])+', ratio='+str(comRatio)+'%.'
    b.label.configure(text=stats)

#==============================================================================
# Essai d'ajuster le nbre de procs
# (essai les classes de sator)
# Note: je ne suis pas sur que ca soit tres utile
#==============================================================================
def adjustNProc():
    global STATS
    if CTK.t == []: return
    try: NProc = int(VARS[0].get())
    except: NProc = 1
    try: comSpeed = float(VARS[1].get())
    except:
        CTK.TXT.insert('START', 'distribute: ComSpeed is invalid.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return
    algo = VARS[2].get()
    useCom = VARS[3].get()

    classes = [8,16,32,64,128,256,512]
    CTK.saveTree()

    for c in range(len(classes)):
        if NProc > classes[c]: break;

    CTK.t, stat1 = D.distribute(CTK.t, classes[c], perfo=(1.,0.,comSpeed),
                                useCom=useCom, algorithm=algo)
    if c > 0:
        CTK.t, stat2 = D.distribute(CTK.t, classes[c-1],
                                    perfo=(1.,0.,comSpeed),
                                    useCom=useCom, algorithm=algo)
    if c < len(classes)-1:
        CTK.t, stat3 = D.distribute(CTK.t, classes[c+1],
                                    perfo=(1.,0.,comSpeed),
                                    useCom=useCom, algorithm=algo)

    # Best
    VARS[0].set(str(classes[c]))
    CTK.t, STATS = D.distribute(CTK.t, classes[c], perfo=(1.,0.,comSpeed),
                                useCom=useCom, algorithm=algo)

    CTK.TXT.insert('START', 'Blocks distributed.\n')
    CTK.TKTREE.updateApp()
    updateStats()
    return

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkDistributor  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Distribute blocks\nover processors.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkDistributor')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- NProc -
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    if 'tkDistributorNProc' in CTK.PREFS:
        V.set(CTK.PREFS['tkDistributorNProc'])
    # -1- ComSpeed -
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)
    if 'tkDistributorComSpeed' in CTK.PREFS:
        V.set(CTK.PREFS['tkDistributorComSpeed'])
    # -2- Algorithm
    V = TK.StringVar(win); V.set('graph'); VARS.append(V)
    if 'tkDistributorAlgorithm' in CTK.PREFS:
        V.set(CTK.PREFS['tkDistributorAlgorithm'])
    # -3- Communication types
    V = TK.StringVar(win); V.set('all'); VARS.append(V)
    if 'tkDistributorComType' in CTK.PREFS:
        V.set(CTK.PREFS['tkDistributorComType'])
    # -4- Manual proc setting
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -5- Multigrid level
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkDistributorMultigrid' in CTK.PREFS:
        V.set(CTK.PREFS['tkDistributorMultigrid'])

    # - NProc -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=3)
    BB = CTK.infoBulle(parent=B, text='Number of processors.')
    B.grid(row=0, column=0, sticky=TK.EW)
    #B.bind('<Return>', distribute)

    # - ComSpeed -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=3)
    BB = CTK.infoBulle(parent=B, text='Weight of communication\ncompared to solver.')
    B.grid(row=0, column=1, sticky=TK.EW)

    # - Multigrid level for split -
    B = TTK.OptionMenu(Frame, VARS[5], '0', '1', '2' )
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Multigrid level.')

    # - Algorithms -
    B = TTK.OptionMenu(Frame, VARS[2], 'graph', 'gradient0', 'gradient1', 'genetic',
                       'fast')
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Distribution algorithm.')

    # - Com types -
    B = TTK.OptionMenu(Frame, VARS[3], 'all', 'match',
                       'overlap', 'bbox', 'none')
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Communication types taken into account while distributing.')

    # - Distribute -
    B = TTK.Button(Frame, text="Distribute tree", command=distribute)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Distribute tree over processors.\nAll zones are reassigned.')

    # - Split and distribute -
    B = TTK.Button(Frame, text="Split and distribute", command=splitAndDistribute)
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Split and distribute tree over processors.\nAll zones are reassigned.')

    # - ProcField -
    #B = TK.Button(Frame, text="Set proc field", command=setProcField)
    #B.grid(row=1, column=2, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Create a field with the attributed proc\nfor each block.')

    # - set proc manually -
    B = TTK.Button(Frame, text="Set proc", command=setProc)
    B.grid(row=3, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set proc manually in selection.')

    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=5)
    BB = CTK.infoBulle(parent=B, text='Processor number.')
    B.grid(row=3, column=1, columnspan=2, sticky=TK.EW)
    B.bind('<Return>', setProc)

    # -  Canvas -
    B = TK.Canvas(Frame, width=250, height=60, background='White')
    B.grid(row=4, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Distribution stats.')
    WIDGETS['canvas'] = B
    WIDGETS['bulle'] = BB

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SolverNoteBook'].add(WIDGETS['frame'], text='tkDistributor')
    except: pass
    CTK.WIDGETS['SolverNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['SolverNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkDistributorNProc'] = VARS[0].get()
    CTK.PREFS['tkDistributorComSpeed'] = VARS[1].get()
    CTK.PREFS['tkDistributorAlgorithm'] = VARS[2].get()
    CTK.PREFS['tkDistributorComType'] = VARS[3].get()
    CTK.PREFS['tkDistributorMultigrid'] = VARS[5].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('10')
    VARS[1].set('0.1')
    VARS[2].set('graph')
    VARS[3].set('all')
    VARS[5].set('0')
    CTK.PREFS['tkDistributorNProc'] = VARS[0].get()
    CTK.PREFS['tkDistributorComSpeed'] = VARS[1].get()
    CTK.PREFS['tkDistributorAlgorithm'] = VARS[2].get()
    CTK.PREFS['tkDistributorComType'] = VARS[3].get()
    CTK.PREFS['tkDistributorMultigrid'] = VARS[5].get()
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
    (win, menu, file, tools) = CTK.minimal('tkDistributor '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
