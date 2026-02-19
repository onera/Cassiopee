# - display probe and view info -
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Cree la liste des variables dans l'arbre (optionMenu)
#==============================================================================
def updateVarNameList(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        zvars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        zvars = C.getVarNames(CTK.t[2][nob][2][noz])
    m = WIDGETS['variable'].children['menu']
    m.delete(0, TK.END)
    if len(zvars) == 0: return
    for i in zvars[0]:
        m.add_command(label=i, command=lambda v=VARS[5],l=i:v.set(l))

#==============================================================================
# Cree la liste des variables dans l'arbre (combobox)
#==============================================================================
def updateVarNameList2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if CTK.__MAINTREE__ <= 0 or nzs == []:
        zvars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        zvars = C.getVarNames(CTK.t[2][nob][2][noz])

    if 'variable' in WIDGETS:
        WIDGETS['variable']['values'] = zvars[0]

#==============================================================================
# Recupere la variable de la souris
#==============================================================================
def getVariableValue(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nob = CTK.Nb[0]+1
    noz = CTK.Nz[0]
    zvars = C.getVarNames(CTK.t[2][nob][2][noz])[0]

    index = CPlot.getActivePointIndex()
    if index != []:
        indv = index[0] # vertex
        inde = index[1] # element

    point = CPlot.getActivePoint()
    field = CPlot.getActivePointF()
    if point == []:
        CTK.TXT.insert('START', 'Please, select only a point.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    cvar = VARS[0].get()
    if cvar == 'CoordinateX':
        VARS[1].set(str(point[0]))
    elif cvar == 'CoordinateY':
        VARS[1].set(str(point[1]))
    elif cvar == 'CoordinateZ':
        VARS[1].set(str(point[2]))
    else:
        ivar = 0
        for v in zvars:
            if v == cvar: break
            if v[0:10] != 'Coordinate': ivar += 1
        VARS[1].set(str(field[ivar]))

#==============================================================================
# Set variable in tree
#==============================================================================
def setVariableValue(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if len(nzs) > 1:
        CTK.TXT.insert('START', 'Please, select only one zone.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    val = CTK.varsFromWidget(VARS[1].get(), type=1)

    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    ind = CPlot.getActivePointIndex()
    if ind == []:
        CTK.TXT.insert('START', 'No selected point.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    indv = ind[0]; inde = ind[1]

    CTK.saveTree()
    zp = Internal.copyTree(z)
    cvar = VARS[0].get()
    svar = cvar.split(':')
    if len(svar) == 2 and svar[0] == 'centers':
        C.setValue(zp, cvar, inde, val[0])
    else:
        C.setValue(zp, cvar, indv, val[0])
    CPlot.replace(CTK.t, nob, noz, zp)
    CTK.TXT.insert('START', 'Point modified.\n')
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkProbe  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Display mesh informations.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=3)
    Frame.columnconfigure(1, weight=0)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkProbe')
    WIDGETS['frameMenu'] = FrameMenu

    #- VARS -
    # -0- Field name -
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -1- Field value -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)

    # - Field variable -
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[0], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList)
        F.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Selected var name.')
        WIDGETS['variable'] = B
    else:
        B = ttk.Combobox(F, textvariable=VARS[0],
                         values=[], state='readonly')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Selected var name.')
        WIDGETS['variable'] = B

    B = TTK.Button(Frame, image=iconics.PHOTO[8],
                   command=getVariableValue, padx=0)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get variable value from mouse.')

    B = TTK.Entry(Frame, textvariable=VARS[1], background='White',
                  width=10)
    B.bind('<Return>', setVariableValue)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Variable value.')

    B = TTK.Button(Frame, text="Set", command=setVariableValue)
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Modify value.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['StateNoteBook'].add(WIDGETS['frame'], text='tkProbe')
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
    (win, menu, file, tools) = CTK.minimal('tkProbe '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
