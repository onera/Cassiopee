# -- tkFilter --
"""Filters for viewing certain zones."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Distributor2.PyTree as D2
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Set filter zones depending on regexp
# IN: t, regexp
# OUT: t modifie et affiche
# Activate : active les matchs et desactive le reste
# Deactivate : desactive les matchs
# Select : selectionne les matchs et deselectionne le reste
#==============================================================================
def setFilter(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # Get filter type
    filterType = VARS[1].get()
    actionType = VARS[2].get()
    # Filter by name
    if filterType == 'By Zone name':
        rexp = VARS[0].get()
        bases = CTK.t[2][1:]
        active = []
        for b in bases:
            baseName = b[0]
            for z in b[2]:
                if z[3] == 'Zone_t':
                    zoneName = baseName + '/' + z[0]
                    i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                    if CTK.matchString(rexp, zoneName): active.append((i,1))
                    else: active.append((i,0))

        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by name.\n')

    # Filter by Zone family
    if filterType == 'By Zone family':
        rexp = VARS[0].get()
        bases = CTK.t[2][1:]
        active = []
        familySelect = C.getFamilyZones(bases, rexp)
        listFmlySlct=[]
        for i in familySelect: listFmlySlct.append(i[0])
        del familySelect
        for b in bases:
            baseName = b[0]
            for z in b[2]:
                if z[3] == 'Zone_t':
                    zoneName = baseName + '/' + z[0]
                    i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                    if z[0] in listFmlySlct: active.append((i,1))
                    else: active.append((i,0))

        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by Zone Family Name.\n')

    # Filter by number
    elif filterType == 'By number':
        no = VARS[0].get()
        try: no = int(no)
        except:
            CTK.TXT.insert('START', 'Filter value must be int.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return

        bases = CTK.t[2][1:]
        active = []
        c = 0
        for b in bases:
            baseName = b[0]
            for z in b[2]:
                if z[3] == 'Zone_t':
                    i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                    if no == c: active.append((i,1))
                    else: active.append((i,0))
                    c += 1
        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by number.\n')

    # Filter by size
    elif filterType  == 'By size >':
        size = VARS[0].get()
        try: size = int(size)
        except:
            CTK.TXT.insert('START', 'Filter value must be int.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return
        bases = CTK.t[2][1:]
        active = []
        for b in bases:
            for z in b[2]:
                if z[3] == 'Zone_t':
                    dim = Internal.getZoneDim(z)
                    if dim[0] == 'Structured': np = dim[1]*dim[2]*dim[3]
                    else: np = dim[1]
                    i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                    if np > size: active.append((i,1))
                    else: active.append((i,0))

        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by size.\n')

    elif filterType  == 'By size <':
        size = VARS[0].get()
        try: size = int(size)
        except:
            CTK.TXT.insert('START', 'Filter value must be int.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return
        bases = CTK.t[2][1:]
        active = []
        for b in bases:
            for z in b[2]:
                if z[3] == 'Zone_t':
                    dim = Internal.getZoneDim(z)
                    if dim[0] == 'Structured':
                        np = dim[1]*dim[2]*dim[3]
                    else:
                        np = dim[1]
                    i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                    if np < size: active.append((i,1))
                    else: active.append((i,0))

        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by size.\n')

    elif filterType  == 'By MG lvl =':
        lvl = VARS[0].get()
        try: lvl = int(lvl)
        except:
            CTK.TXT.insert('START', 'Filter value must be int.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return
        fac = 2**lvl
        bases = CTK.t[2][1:]
        active = []
        for b in bases:
            for z in b[2]:
                if z[3] == 'Zone_t':
                    dim = Internal.getZoneDim(z)
                    i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                    if dim[0] == 'Structured':
                        celldim = dim[4]
                        ni = dim[1]; nj = dim[2]; nk = dim[3]
                        if celldim == 2:
                            if ((ni-1)%fac == 0 and (nj-1)%fac == 0):
                                active.append((i,1))
                            else: active.append((i,0))
                        elif celldim == 3:
                            if ((ni-1)%fac == 0 and (nj-1)%fac == 0 and (nk-1)%fac == 0):
                                active.append((i,1))
                            else: active.append((i,0))
                        else: active.append((i,0))
                    else: active.append((i,0))
        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by multigrid level.\n')

    elif filterType  == 'By MG lvl !=':
        lvl = VARS[0].get()
        try: lvl = int(lvl)
        except:
            CTK.TXT.insert('START', 'Filter value must be int.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return
        fac = 2**lvl
        bases = CTK.t[2][1:]
        active = []
        for b in bases:
            for z in b[2]:
                if z[3] == 'Zone_t':
                    dim = Internal.getZoneDim(z)
                    if dim[0] == 'Structured':
                        celldim = dim[4]
                        ni = dim[1]; nj = dim[2]; nk = dim[3]
                        i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                        if celldim == 2:
                            if ((ni-1)%fac == 0 and (nj-1)%fac == 0):
                                active.append((i,0))
                            else: active.append((i,1))
                        elif celldim == 3:
                            if ((ni-1)%fac == 0 and (nj-1)%fac == 0 and (nk-1)%fac == 0):
                                active.append((i,0))
                            else: active.append((i,1))
                        else: active.append((i,0))
                    else: active.append((i,0))

        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by multigrid level.\n')

    elif filterType == 'By proc':
        myProc = VARS[0].get()
        try: myProc = int(myProc)
        except:
            CTK.TXT.insert('START', 'Filter value must be int.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return
        bases = CTK.t[2][1:]
        active = []
        for b in bases:
            for z in b[2]:
                if z[3] == 'Zone_t':
                    i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                    proc = Internal.getNodeFromName(z, 'proc')
                    if proc is not None and D2.getProc(z) == myProc:
                        active.append((i,1))
                    else:
                        active.append((i,0))
        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by proc.\n')

    elif filterType == 'By priority':
        prio = VARS[0].get()
        try: prio = int(prio)
        except:
            CTK.TXT.insert('START', 'Filter value must be int.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return
        bases = CTK.t[2][1:]
        active = []
        for b in bases:
            prios = Internal.getNodesFromName2(b, 'Priority')
            if prios != [] and prios[0][1][0,0] == prio:
                for z in b[2]:
                    if z[3] == 'Zone_t':
                        i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                        active.append((i,1))
            elif prios == [] and prio == 0:
                for z in b[2]:
                    if z[3] == 'Zone_t':
                        i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                        active.append((i,1))
            else:
                for z in b[2]:
                    if z[3] == 'Zone_t':
                        i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                        active.append((i,0))
        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by priority.\n')

    elif filterType  == 'By formula (and)' or filterType == 'By formula (or)':
        if filterType  == 'By formula (and)': testType = 0
        else: testType = 1
        formula = VARS[0].get()
        b1 = formula.find('{')
        b2 = formula.find('}')
        if b1 == -1 or b2 == -1:
            CTK.TXT.insert('START', 'Filter value must be int.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return

        field = formula[b1+1:b2]; # print field
        try: # try an eval
            ff = formula.replace('{'+field+'}', '12.')
            res = eval(ff)
        except:
            CTK.TXT.insert('START', 'Filter supports only simple formulas.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
            return
        b1 = formula.find('==')
        formula1 = formula; formula2 = formula
        if b1 != -1:
            formula1 = formula.replace('==', '<=')
            formula2 = formula.replace('==', '>=')
        else:
            b1 = formula.find('=')
            if b1 != -1:
                formula1 = formula.replace('=', '<=')
                formula2 = formula.replace('=', '>=')

        bases = CTK.t[2][1:]
        active = []
        for b in bases:
            for z in b[2]:
                if z[3] == 'Zone_t':
                    fmin = C.getMinValue(z, field)
                    fmax = C.getMaxValue(z, field)
                    ff = formula1.replace('{'+field+'}', str(fmin))
                    res1 = eval(ff)
                    ff = formula2.replace('{'+field+'}', str(fmax))
                    res2 = eval(ff)
                    i = CPlot.getCPlotNumber(CTK.t, b[0], z[0])
                    if testType == 0:
                        if res1 == True and res2 == True: active.append((i,1))
                        else: active.append((i,0))
                    else:
                        if res1 == True or res2 == True: active.append((i,1))
                        else: active.append((i,0))
        if actionType == 'Activate': CPlot.setActiveZones(active)
        elif actionType == 'Deactivate':
            inactive = []
            for i in active:
                if i[1] == 1: inactive.append((i[0],0))
            CPlot.setActiveZones(inactive)
        else: CPlot.setSelectedZones(active)
        CTK.TXT.insert('START', 'Filtered by formula.\n')

#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFilter  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Display/select certain blocks.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkFilter')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- filter value -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    if 'tkFilterValue' in CTK.PREFS: V.set(CTK.PREFS['tkFilterValue'])
    # -1- filter type
    V = TK.StringVar(win); V.set('By Zone name'); VARS.append(V)
    if 'tkFilterName' in CTK.PREFS: V.set(CTK.PREFS['tkFilterName'])
    # -2- action
    V = TK.StringVar(win); V.set('Activate'); VARS.append(V)
    if 'tkFilterAction' in CTK.PREFS: V.set(CTK.PREFS['tkFilterAction'])

    # - Options -
    B = TTK.OptionMenu(Frame, VARS[2], 'Activate', 'Deactivate', 'Select')
    BB = CTK.infoBulle(parent=B, text='Action to be performed on filtered zones.')
    B.grid(row=0, column=0, sticky=TK.EW)

    B = TTK.OptionMenu(Frame, VARS[1], 'By Zone name', 'By Zone family' ,'By size >',
                       'By size <', 'By MG lvl =', 'By MG lvl !=',
                       'By proc', 'By priority', 'By number',
                       'By formula (or)', 'By formula (and)')
    BB = CTK.infoBulle(parent=B, text='Filter criteria.')
    B.grid(row=0, column=1, sticky=TK.EW)

    # - Set filter -
    B = TTK.Button(Frame, text="Filter", command=setFilter)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Apply filter.')

    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=20)
    BB = CTK.infoBulle(parent=B, text='Matching string.')
    B.grid(row=1, column=1, sticky=TK.EW)
    B.bind('<Return>', setFilter)

#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['StateNoteBook'].add(WIDGETS['frame'], text='tkFilter')
    except: pass
    CTK.WIDGETS['StateNoteBook'].select(WIDGETS['frame'])

#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['StateNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkFilterValue'] = VARS[0].get()
    CTK.PREFS['tkFilterName'] = VARS[1].get()
    CTK.PREFS['tkFilterAction'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('')
    VARS[1].set('By Name')
    VARS[2].set('Activate')
    CTK.PREFS['tkFilterValue'] = VARS[0].get()
    CTK.PREFS['tkFilterName'] = VARS[1].get()
    CTK.PREFS['tkFilterAction'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def displayFrameMenu(event=None):
    WIDGETS['frameMenu'].tk_popup(event.x_root+50, event.y_root, 0)

#==============================================================================
if __name__ == "__main__":
    import sys
    if (len(sys.argv) == 2):
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkFilter '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
