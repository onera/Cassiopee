# - tkVariables -
"""Variable manager."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Post.PyTree as P
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Pour remove (optionMenu)
def updateVarNameList1(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        varsl = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        varsl = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
    m = WIDGETS['var1'].children['menu']
    m.delete(0, TK.END)
    zvars = ['All', 'FlowSolutionNodes', 'FlowSolutionCenters']
    if len(varsl) != 0: zvars += varsl[0]
    for i in zvars:
        m.add_command(label=i, command=lambda v=VARS[5],l=i:v.set(l))

# Pour remove (combobox)
def updateVarNameList1_2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        varsl = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        varsl = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)

    zvars = ['All', 'FlowSolutionNodes', 'FlowSolutionCenters']
    if len(varsl) != 0: zvars += varsl[0]
    if 'var1' in WIDGETS:
        WIDGETS['var1']['values'] = zvars

#==============================================================================
# Pour le gradient (optionMenu)
def updateVarNameList2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        zvars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        zvars = C.getVarNames(CTK.t[2][nob][2][noz])
    m = WIDGETS['var2'].children['menu']
    m.delete(0, TK.END)
    if len(zvars) == 0: return
    for i in zvars[0]:
        m.add_command(label=i, command=lambda v=VARS[2],l=i:v.set(l))

# Pour le gradient (combobox)
def updateVarNameList2_2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        zvars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        zvars = C.getVarNames(CTK.t[2][nob][2][noz])    
    if len(zvars) == 0: return
    if 'var2' in WIDGETS:
        WIDGETS['var2']['values'] = zvars[0]

#==============================================================================
# Pour center2Node - seult les variables en centres (optionMenu)
def updateVarNameList3(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        varsl = C.getVarNames(CTK.t, excludeXYZ=True, loc='centers')
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        varsl = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True,
                              loc='centers')
    m = WIDGETS['var3'].children['menu']
    m.delete(0, TK.END)
    zvars = ['FlowSolutionCenters']
    if len(varsl) != 0: zvars += varsl[0] 
    for i in zvars:
        m.add_command(label=i, command=lambda v=VARS[8],l=i:v.set(l))

# Pour center2Node - seult les variables en centres (combobox)
def updateVarNameList3_2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        varsl = C.getVarNames(CTK.t, excludeXYZ=True, loc='centers')
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        varsl = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True,
                              loc='centers')

    zvars = ['FlowSolutionCenters']
    if len(varsl) != 0: zvars += varsl[0]

    if 'var3' in WIDGETS:
        WIDGETS['var3']['values'] = zvars

#==============================================================================
# Pour node2Center - seult les variables en noeuds (optionMenu)
def updateVarNameList4(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        varsl = C.getVarNames(CTK.t, excludeXYZ=True, loc='nodes')
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        varsl = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True,
                              loc='nodes')
    m = WIDGETS['var4'].children['menu']
    m.delete(0, TK.END)
    zvars = ['FlowSolutionNodes']
    if len(varsl) != 0: zvars += varsl[0]
    for i in zvars:
        m.add_command(label=i, command=lambda v=VARS[9],l=i:v.set(l))

# Pour node2Center - seult les variables en noeuds (combobox)
def updateVarNameList4_2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        varsl = C.getVarNames(CTK.t, excludeXYZ=True, loc='nodes')
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        varsl = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True,
                              loc='nodes')
    zvars = ['FlowSolutionNodes']
    if len(varsl) != 0: zvars += varsl[0]
    if 'var4' in WIDGETS:
        WIDGETS['var4']['values'] = zvars

#==============================================================================
# Pour rename (optionMenu)
def updateVarNameList5(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        varsl = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        varsl = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)
    m = WIDGETS['var5'].children['menu']
    m.delete(0, TK.END)
    for i in varsl[0]:
        m.add_command(label=i, command=lambda v=VARS[11],l=i:v.set(l))

# Pour rename (combobox)
def updateVarNameList5_2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        varsl = C.getVarNames(CTK.t, excludeXYZ=True)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        varsl = C.getVarNames(CTK.t[2][nob][2][noz], excludeXYZ=True)

    if 'var5' in WIDGETS:
        WIDGETS['var5']['values'] = varsl[0]

#==============================================================================
def rmVar():
    if CTK.t == []: return
    CTK.saveTree()
    var = VARS[5].get()
    if var == 'All':
        C._rmVars(CTK.t, Internal.__FlowSolutionNodes__)
        C._rmVars(CTK.t, Internal.__FlowSolutionCenters__)
    elif var == 'FlowSolutionNodes':
        C._rmVars(CTK.t, Internal.__FlowSolutionNodes__)
    elif var == 'FlowSolutionCenters':
        C._rmVars(CTK.t, Internal.__FlowSolutionCenters__)
    else: C._rmVars(CTK.t, var)
    CTK.TXT.insert('START', 'Variable %s removed.\n'%var)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
def renameVar():
    if CTK.t == []: return
    CTK.saveTree()
    varp = VARS[11].get()
    varn = VARS[10].get()
    CTK.t = P.renameVars(CTK.t, [varp], [varn])
    CTK.TXT.insert('START', 'Variable %s replaced in all tree.\n'%varp)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def center2NodeVar():
    if CTK.t == []: return
    CTK.saveTree()
    var = VARS[8].get()
    if var == 'FlowSolutionCenters':
        CTK.t = C.center2Node(CTK.t, Internal.__FlowSolutionCenters__)
    else: CTK.t = C.center2Node(CTK.t, var)
    CTK.TXT.insert('START', 'Variable %s put to nodes in all tree.\n'%var)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def node2CenterVar():
    if CTK.t == []: return
    CTK.saveTree()
    var = VARS[9].get()
    if var == 'FlowSolutionNodes':
        CTK.t = C.node2Center(CTK.t, Internal.__FlowSolutionNodes__)
    else: CTK.t = C.node2Center(CTK.t, var)
    CTK.TXT.insert('START', 'Variable %s put to centers in all tree.\n'%var)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def chooseImportFile(event=None):
    try: import tkFileDialog
    except: import tkinter.filedialog as tkFileDialog 
    init = VARS[4].get()
    init = init.split(';')[0]
    files = tkFileDialog.askopenfilenames(
        filetypes=CTK.fileTypes, initialfile=init, multiple=1)
    if files == '' or files == None or files == (): # user cancel
        return
    # strangely, initfile is part of the return
    files = CTK.fixFileString__(files, init)
    s = ''
    for f in files: s += f+';'
    VARS[4].set(s)
    importFile()

#==============================================================================
def importFile(event=None):
    if CTK.t == []: return
    s = VARS[4].get(); s = s.split(';')
    try:
        t1 = []
        for filename in s:
            if filename != '':
                t2 = C.convertFile2PyTree(filename)
                # Fusion des bases de t et t2
                if t1 == []: t1 = t2
                else: t1 = C.mergeTrees(t1, t2)
    except:
        CTK.TXT.insert('START', 'Import failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if t1 == []:
        CTK.TXT.insert('START', 'Import failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()
    nzs = CPlot.getSelectedZones()

    # Essaie de trouver une methode adaptee
    method = 1 # match by geom
    import sets
    zoneNames = sets.Set(C.getZoneNames(CTK.t, prefixByBase=False))
    zoneNames1 = sets.Set(C.getZoneNames(t1, prefixByBase=False))
    inter = zoneNames & zoneNames1
    linter = len(inter)*1.
    comp = min(len(zoneNames), len(zoneNames1))*1.
    if linter / comp > 0.9: method = 0 # try match by name (safer)

    if CTK.__MAINTREE__ <= 0 or nzs == []:
        CTK.t = P.importVariables(t1, CTK.t, method=method)
    else:
        zones = C.newPyTree(['Base'])
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            zone = CTK.t[2][nob][2][noz]
            zones[2][1][2].append(zone)
        zones = P.importVariables(t1, zones, method=method)
        c = 0
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            CTK.t[2][nob][2][noz] = zones[2][1][2][c]; c += 1
    CTK.TXT.insert('START', 'Variable file %s imported.\n'%filename)
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
def computeVariables():
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    varname = VARS[0].get()

    # Adimensionnement - RefState must exist (full)
    state = Internal.getNodeFromName(CTK.t, 'ReferenceState')
    if state is None:
        CTK.TXT.insert('START', 'ReferenceState is missing (tkState).\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # If gamma exists -> suppose full RefState, we take all from it
    gamma = Internal.getNodeFromName(state, 'Gamma')
    fail = False
    if gamma is not None: gamma = Internal.getValue(gamma)
    else: fail = True

    Cv = Internal.getNodeFromName(state, 'Cv')
    if Cv is not None: Cv = Internal.getValue(Cv)
    else: fail = True

    TInf = Internal.getNodeFromName(state, 'Temperature')
    if TInf is not None: TInf = Internal.getValue(TInf)
    else: fail = True 

    Ts = Internal.getNodeFromName(state, 'Ts')
    if Ts is not None: Ts = Internal.getValue(Ts)
    else: fail = True

    Cs = Internal.getNodeFromName(state, 'Cs')
    if Cs is not None: Cs = Internal.getValue(Cs)
    else: fail = True

    Mus = Internal.getNodeFromName(state, 'Mus')
    if Mus is not None: Mus = Internal.getValue(Mus)
    else: fail = True

    if fail:
        CTK.TXT.insert('START', 'Reference state is not full. Use tkState.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    rgp = (gamma-1)*Cv
    loc = VARS[6].get()    
    CTK.saveTree()
    if (varname == 'Vorticity' or varname == 'VorticityMagnitude' or
        varname == 'QCriterion' or varname == 'ShearStress' or
        varname == 'SkinFriction' or varname == 'SkinFrictionTangential'): # extra variables
        varloc = loc+':'+varname
        if CTK.__MAINTREE__ <= 0 or nzs == []:
            try:
                CTK.t = P.computeExtraVariable(CTK.t, varloc, gamma=gamma,
                                               rgp=rgp, Cs=Cs, mus=Mus,
                                               Ts=Ts)
                CTK.TXT.insert('START', 'Variable %s computed.\n'%varloc)
            except Exception as e:
                Panels.displayErrors([0,str(e)], header='Error: computeExtraVariables')
                CTK.TXT.insert('START', 'Computation of variable %s failed.\n'%varloc)
                CTK.TXT.insert('START', 'Error: ', 'Error')

        else:
            fail = False; errors = []
            for nz in nzs:
                nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
                try:
                    CTK.t[2][nob][2][noz] = \
                        P.computeExtraVariable(CTK.t[2][nob][2][noz], varloc,
                                               gamma=gamma, rgp=rgp, Cs=Cs,
                                               mus=Mus, Ts=Ts)
                except Exception as e:
                    fail = True; errors += [0,str(e)]

            if not fail:
                CTK.TXT.insert('START', 'Variable %s computed.\n'%varloc)
            else:
                Panels.displayErrors(errors, header='Error: computeExtraVariables')
                CTK.TXT.insert('START', 'Computation of variable %s failed.\n'%varloc)
                CTK.TXT.insert('START', 'Error: ', 'Error')

    else: # std variables 
        varloc = loc+':'+varname
        if CTK.__MAINTREE__ <= 0 or nzs == []:
            try:
                CTK.t = P.computeVariables(CTK.t, [varloc],
                                           gamma=gamma, rgp=rgp, Cs=Cs,
                                           mus=Mus, Ts=Ts)
                CTK.TXT.insert('START', 'Variable %s computed.\n'%varloc)
            except Exception as e:
                Panels.displayErrors([0,str(e)], header='Error: computeVariables')
                CTK.TXT.insert('START', 'Computation of variable %s failed.\n'%varloc)
                CTK.TXT.insert('START', 'Error: ', 'Error')   
        else:
            fail = False; errors = []
            for nz in nzs:
                nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
                try:
                    CTK.t[2][nob][2][noz] = \
                        P.computeVariables(CTK.t[2][nob][2][noz], [varloc],
                                           gamma=gamma, rgp=rgp, Cs=Cs,
                                           mus=Mus, Ts=Ts)
                except Exception as e:
                    fail = True; errors += [0,str(e)]

            if not fail: 
                CTK.TXT.insert('START', 'Variable %s computed.\n'%varloc)
            else:
                Panels.displayErrors(errors, header='Error: computeVariables')
                CTK.TXT.insert('START', 'Computation of variable %s failed.\n'%varloc)
                CTK.TXT.insert('START', 'Error: ', 'Error')
    #C._fillMissingVariables(CTK.t)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
def addVar(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    varname = VARS[1].get()
    CTK.saveTree()
    s = varname.split('=')
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        if len(s) > 1: C._initVars(CTK.t, varname)
        else: C._addVars(CTK.t, varname)
    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            if len(s) > 1: C._initVars(z, varname)
            else: C._addVars(z, varname)
    CTK.TXT.insert('START', 'Variable %s added.\n'%varname)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
def computeGrad():
    if CTK.t == []: return
    varname = VARS[2].get()
    CTK.saveTree()
    try: CTK.t = P.computeGrad(CTK.t, varname)
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: computeGrad')
        CTK.TXT.insert('START', 'Gradient computation failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.TXT.insert('START', 'Gradient of %s computed on all tree.\n'%varname)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
def computeNormGrad():
    if CTK.t == []: return
    varname = VARS[2].get()
    CTK.saveTree()
    try: CTK.t = P.computeNormGrad(CTK.t, varname)
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: computeNormGrad')
        CTK.TXT.insert('START', 'Gradient\'s norm computation failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.TXT.insert('START', 'Gradient\'s norm of %s computed on all tree.\n'%varname)
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
def computeNormCurl():
    if CTK.t == []: return
    zvars = VARS[3].get()
    zvars = zvars.replace(' ', '')
    zvars = zvars.split(';')
    CTK.saveTree()
    try: CTK.t = P.computeNormCurl(CTK.t, zvars)
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: computeNormCurl')
        CTK.TXT.insert('START', 'Curl\'s norm computation failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.TXT.insert('START', 'Curl\'s norm computed on all tree.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
def computeCurl():
    if CTK.t == []: return
    zvars = VARS[3].get()
    zvars = zvars.replace(' ', '')
    zvars = zvars.split(';')
    CTK.saveTree()
    try: CTK.t = P.computeCurl(CTK.t, zvars)
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: computeCurl')
        CTK.TXT.insert('START', 'Curl computation failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.TXT.insert('START', 'Curl computed on all tree.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()

#==============================================================================
# Fill missing variables
#==============================================================================
def fillMissingVariables():
    if CTK.t == []: return
    CTK.saveTree()
    C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Missing variables filled.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkVariable  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage field variables.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkVariables')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- computeVariable name
    V = TK.StringVar(win); V.set('Pressure'); VARS.append(V)
    if 'tkVariablesName' in CTK.PREFS: 
        V.set(CTK.PREFS['tkVariablesName'])
    # -1- addVar
    V = TK.StringVar(win); V.set('Density'); VARS.append(V)
    if 'tkVariablesAddVar' in CTK.PREFS: 
        V.set(CTK.PREFS['tkVariablesAddVar'])
    # -2- computeGrad -
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -3- computeCurl -
    V = TK.StringVar(win); V.set('MomentumX;MomentumY;MomentumZ')
    VARS.append(V)
    # -4- importFile -
    V = TK.StringVar(win); V.set('output.plt'); VARS.append(V)
    if 'tkVariablesImportFile' in CTK.PREFS: 
        V.set(CTK.PREFS['tkVariablesImportFile'])
    # -5- Rm variable
    V = TK.StringVar(win); V.set('All'); VARS.append(V)
    # -6- Var location
    V = TK.StringVar(win); V.set('centers'); VARS.append(V)
    if 'tkVariablesLoc' in CTK.PREFS: 
        V.set(CTK.PREFS['tkVariablesLoc'])
    # -7- adim type
    V = TK.StringVar(win); V.set('Adim1 (ro,a,T)'); VARS.append(V)
    if 'tkVariablesAdim' in CTK.PREFS: 
        V.set(CTK.PREFS['tkVariablesAdim'])
    # -8- center2Node variable
    V = TK.StringVar(win); V.set('FlowSolutionCenters'); VARS.append(V)
    # -9- node2Center variable
    V = TK.StringVar(win); V.set('FlowSolutionNodes'); VARS.append(V)
    # -10- renameVar variable - new
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -11- renameVar variable - prev
    V = TK.StringVar(win); V.set('x'); VARS.append(V)

    # - importFile -
    norow = 0
    B = TTK.Button(Frame, text="Import file", command=importFile)
    B.grid(row=norow, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Import solution file into existing pyTree.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=0)
    B = TTK.Entry(F, textvariable=VARS[4], background='White')
    B.bind('<Return>', importFile)
    B.grid(row=norow, column=0, sticky=TK.EW)
    B = TTK.Button(F, text="...", padx=0, command=chooseImportFile)
    BB = CTK.infoBulle(parent=B, text='Select solution file.')
    B.grid(row=norow, column=1, sticky=TK.EW)
    F.grid(row=norow, column=1, sticky=TK.EW)

    # - addVar -
    norow += 1
    B = TTK.Button(Frame, text="Add/modify variable", command=addVar)
    B.grid(row=norow, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Add a new variable into pyTree. A formula of type: Density=3*{CoordinateX} can be specified.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White')
    B.bind('<Return>', addVar)
    B.grid(row=norow, column=1, sticky=TK.EW)

    # - rmVar -
    norow += 1
    B = TTK.Button(Frame, text="Rm variable", command=rmVar)
    B.grid(row=norow, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Rm variable from pyTree.')

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[5], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList1)
        F.grid(row=norow, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Removed variable.')
        WIDGETS['var1'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[5], 
                         values=[], state='readonly')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList1_2)
        F.grid(row=norow, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Removed variable.')
        WIDGETS['var1'] = B

    # - renameVar - 
    #norow+= 1
    #F = TTK.Frame(Frame, borderwidth=0)
    #F.columnconfigure(0, weight=1)
    #if ttk is None:
    #    B = TK.OptionMenu(F, VARS[11], '')
    #    B.grid(sticky=TK.EW)
    #    F.bind('<Enter>', updateVarNameList5)
    #    F.grid(row=norow, column=0, sticky=TK.EW)
    #    BB = CTK.infoBulle(parent=B, text='Renamed variable.')
    #    WIDGETS['var5'] = B
    #else:
    #    B = ttk.Combobox(F, textvariable=VARS[11], 
    #                     values=[], state='readonly')
    #    B.grid(sticky=TK.EW)
    #    F.bind('<Enter>', updateVarNameList5_2)
    #    F.grid(row=norow, column=0, sticky=TK.EW)
    #    BB = CTK.infoBulle(parent=B, text='Renamed variable.')
    #    WIDGETS['var5'] = B
    #B = TK.Entry(Frame, textvariable=VARS[10], background='White')
    #B.bind('<Return>', renameVar)
    #B.grid(row=norow, column=1, sticky=TK.EW)
    #norow+=1
    #B = TK.Button(Frame, text="Rename variable", command=renameVar)
    #B.grid(row=norow, column=0, columnspan=2, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Rename variable from pyTree.') 

    # - center2Node var -
    norow += 1
    B = TTK.Button(Frame, text="Center2Node", command=center2NodeVar)
    B.grid(row=norow, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Put a center variable to nodes in pyTree.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[8], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList3)
        F.grid(row=norow, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Center variable to be set in nodes.')
        WIDGETS['var3'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[8], 
                         values=[], state='readonly')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList3_2)
        F.grid(row=norow, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Center variable to be set in nodes.')
        WIDGETS['var3'] = B

    # - node2Center var -
    norow+=1
    B = TTK.Button(Frame, text="Node2Center", command=node2CenterVar)
    B.grid(row=norow, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Put a node variable to centers in pyTree.')

    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[9], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList4)
        F.grid(row=norow, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Node variable to be set in centers.')
        WIDGETS['var4'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[9], 
                         values=[], state='readonly')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList4_2)
        F.grid(row=norow, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Node variable to be set in centers.')
        WIDGETS['var4'] = B

    # - computeGrad -
    norow+=1
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=1)
    B = TTK.Button(F, text="Grad", command=computeGrad)
    BB = CTK.infoBulle(parent=B, text='Compute gradient of variables.')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text="Norm", command=computeNormGrad)
    BB = CTK.infoBulle(parent=B, text='Compute gradient\' norm of variables.')
    B.grid(row=0, column=1, sticky=TK.EW)
    F.grid(row=norow, column=0, sticky=TK.EW)
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    if ttk is None:
        B = TK.OptionMenu(F, VARS[2], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=norow, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable for gradient.')
        WIDGETS['var2'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[2], 
                         values=[], state='readonly')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2_2)
        F.grid(row=norow, column=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Variable for gradient.')
        WIDGETS['var2'] = B

    # - computeCurl -
    norow += 1
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)
    F.columnconfigure(1, weight=1)
    B = TTK.Button(F, text="Curl ", command=computeCurl)
    BB = CTK.infoBulle(parent=B, text='Compute curl of variables.')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text="Norm", command=computeNormCurl)
    BB = CTK.infoBulle(parent=B, text='Compute gradient\' norm of variables.')
    B.grid(row=0, column=1, sticky=TK.EW)
    F.grid(row=norow, column=0, sticky=TK.EW)

    B = TTK.Entry(Frame, textvariable=VARS[3], background='White')
    BB = CTK.infoBulle(parent=B, text='Variables for curl.')
    B.grid(row=norow, column=1, sticky=TK.EW)

    # - computeVariables -
    norow += 1
    #B = TTK.OptionMenu(Frame, VARS[7], 'Adim1 (ro,a,T)', 'Adim2 (ro,u,T)', 'dim')
    #B.grid(row=norow, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Use this adimensioning for variable computation.')
    B = TTK.OptionMenu(Frame, VARS[6], 'nodes', 'centers')
    B.grid(row=norow, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Computed variable will be localized here.')

    norow += 1
    B = TTK.Button(Frame, text="Compute variable", command=computeVariables)
    B.grid(row=norow, column=0, sticky=TK.EW)
    B = TTK.OptionMenu(Frame, VARS[0],'Pressure','Temperature',
                       'VelocityMagnitude','VelocityX',\
                       'VelocityY','VelocityZ','Enthalpy','Entropy','Mach',\
                       'ViscosityMolecular','PressureStagnation',\
                       'TemperatureStagnation', 'PressureDynamic',
                       'Vorticity', 'VorticityMagnitude', 'QCriterion',
                       'ShearStress', 'SkinFriction', 'SkinFrictionTangential')
    B.grid(row=norow, column=1, sticky=TK.EW)

    # fill missing variables
    norow+=1
    B = TTK.Button(Frame, text="Fill missing variables",
                   command=fillMissingVariables)
    B.grid(row=norow, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='All zones will have the same variables.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['PostNoteBook'].add(WIDGETS['frame'], text='tkVariables')
    except: pass
    CTK.WIDGETS['PostNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['PostNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkVariablesName'] = VARS[0].get()
    CTK.PREFS['tkVariablesAddVar'] = VARS[1].get()
    CTK.PREFS['tkVariablesImportFile'] = VARS[4].get()
    CTK.PREFS['tkVariablesLoc'] = VARS[6].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('Pressure')
    VARS[1].set('Density')
    VARS[4].set('output.plt')
    VARS[6].set('centers')
    CTK.PREFS['tkVariablesName'] = VARS[0].get()
    CTK.PREFS['tkVariablesAddVar'] = VARS[1].get()
    CTK.PREFS['tkVariablesImportFile'] = VARS[4].get()
    CTK.PREFS['tkVariablesLoc'] = VARS[6].get()
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
    (win, menu, file, tools) = CTK.minimal('tkVariable '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
