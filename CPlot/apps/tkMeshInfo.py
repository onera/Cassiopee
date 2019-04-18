# - display mesh info -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# String with groups
# IN: un entier
# OUT: une chaine ou l'on separe par groupe de 3 les nombres
#==============================================================================
def strwg__(i):
    s = '%d' % i
    groups = []
    while s and s[-1].isdigit():
        groups.append(s[-3:])
        s = s[:-3]
    return s + ' '.join(reversed(groups))

#================================================================================
# Compute mesh information : nnodes, ncells, nfaces, ...
# IN: dim = getZoneDim
#================================================================================
def computeMeshInfo(z, dim):
    if dim[0] == 'Structured':
        # Compute nodes number
        np = dim[1]*dim[2]*dim[3]
        # Compute cell number
        if dim[4] == 1:
            ncells = dim[1]-1
        elif dim[4] == 2:
            ncells = (dim[1]-1)*(dim[2]-1)
        else:
            ncells = (dim[1]-1)*(dim[2]-1)*(dim[3]-1)
        # Compute face number
        if dim[4] == 1:
            nfaces = dim[1]
        elif dim[4] == 2:
            nfaces = dim[1]*(dim[2]-1)+dim[2]*(dim[1]-1)
        elif dim[4] == 3:
            nfaces = dim[1]*(dim[2]-1)*(dim[3]-1)+ \
                     dim[2]*(dim[3]-1)*(dim[1]-1)+ \
                     dim[3]*(dim[1]-1)*(dim[2]-1)
    else:
        # Compute nodes number
        np = dim[1]
        # Compute cell number
        ncells = dim[2]
        # Compute face number
        nfaces = 0
        if dim[3] == 'BAR': nfaces = np
        elif dim[3] == 'NGON':
            # Find NFACE
            elts = Internal.getNodesFromType2(z, 'Elements_t')
            nfaces = 0
            for e in elts:
                if Internal.getValue(e)[0] == 22:
                    erange = Internal.getNodeFromName1(e, 'ElementRange')
                    nfaces += erange[1][1]-erange[1][0]+1
        elif dim[3] == 'NODE':
            ncells = 0
            nfaces = 0
        else: # BE, very expensive!!
            if dim[1] < 100000: # Trouver un autre moyen
                zp = C.convertArray2NGon(z)
                elts = Internal.getNodesFromType2(zp, 'Elements_t')
                nfaces = 0
                for e in elts:
                    if Internal.getValue(e)[0] == 22:
                        erange = Internal.getNodeFromName1(e, 'ElementRange')
                        nfaces += erange[1][1]-erange[1][0]+1
            else: nfaces = 0
    return np, ncells, nfaces

#==============================================================================
# Cree la liste des variables dans l'arbre (optionMenu)
#==============================================================================
def updateVarNameList(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if (CTK.__MAINTREE__ <= 0 or nzs == []):
        vars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz])
    m = WIDGETS['variable'].children['menu']
    m.delete(0, TK.END)
    if len(vars) == 0: return
    for i in vars[0]:
        m.add_command(label=i, command=lambda v=VARS[5],l=i:v.set(l))

#==============================================================================
# Cree la liste des variables dans l'arbre (combobox)
#==============================================================================
def updateVarNameList2(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if (CTK.__MAINTREE__ <= 0 or nzs == []):
        vars = C.getVarNames(CTK.t)
    else:
        nob = CTK.Nb[0]+1
        noz = CTK.Nz[0]
        vars = C.getVarNames(CTK.t[2][nob][2][noz])

    if 'variable' in WIDGETS:
        WIDGETS['variable']['values'] = vars[0]

#==============================================================================
def updateInfo(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    
    npTot = 0; ncellsTot = 0; nfacesTot = 0
    nzones = 0
    minv = 1e6; maxv = -1e6
    var = VARS[5].get()
    failed = False
    if nzs == []:
        zones = Internal.getZones(CTK.t)
        for z in zones:
            dim = Internal.getZoneDim(z)
            try: minv = min(minv, C.getMinValue(z, var))
            except: failed = True
            try: maxv = max(maxv, C.getMaxValue(z, var))
            except: failed = True
            np, ncells, nfaces = computeMeshInfo(z, dim)
            npTot += np
            ncellsTot += ncells
            nfacesTot += nfaces
            nzones += 1
    else:
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            dim = Internal.getZoneDim(z)
            minv = min(minv, C.getMinValue(z, var))
            maxv = max(maxv, C.getMaxValue(z, var))
            np, ncells, nfaces = computeMeshInfo(z, dim)
            npTot += np
            ncellsTot += ncells
            nfacesTot += nfaces
            nzones += 1

    VARS[0].set(strwg__(npTot))
    VARS[6].set(strwg__(ncellsTot))
    VARS[7].set(strwg__(nfacesTot))
    VARS[1].set(str(minv))
    VARS[2].set(str(maxv))
    VARS[4].set('')
    
    if nzs == []:
        VARS[3].set('Tree')
        if nzones == 0 or nzones == 1:
            VARS[4].set(strwg__(nzones)+' zone')
        else: VARS[4].set(strwg__(nzones)+' zones')
    elif len(nzs) > 1:
        VARS[3].set('Multiple')
        VARS[4].set(str(nzones)+' zones')
    else:
        nz = nzs[0]
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        VARS[3].set(CTK.t[2][nob][0]+'/'+z[0])
        dim = Internal.getZoneDim(z)
        if dim[0] == 'Structured':
            VARS[4].set('Structured')
            VARS[0].set(str(npTot) + ' ('+str(dim[1])+'x'+str(dim[2])+'x'+str(dim[3])+')')
        else:
            VARS[4].set(dim[3])
    if not failed:
        CTK.TXT.insert('START', 'Info updated.\n')
    else:
        CTK.TXT.insert('START', 'Variable min-max was not computed.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkMeshInfo', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Display mesh informations.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkMeshInfo')
    WIDGETS['frameMenu'] = FrameMenu

    #- VARS -
    # -0- Npts -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- Var min -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -2- Var max -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -3- Selected block name
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -4- Selected block type
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -5- Selected variable
    V = TK.StringVar(win); V.set('CoordinateX'); VARS.append(V)
    # -6- NCells -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -7- NFaces -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    

    # - selected block name -
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=15)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Currently selected zone name.')

    # - selected block type -
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White')
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Currently selected zone type.')
    
    # - Npts -
    B = TTK.Label(Frame, text="Npts")
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of points in selection.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White',
                  width=10)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of points in selection.')

    # - Ncells -
    B = TTK.Label(Frame, text="NCells")
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of cells in selection.')
    B = TTK.Entry(Frame, textvariable=VARS[6], background='White',
                  width=10)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of cells in selection.')

    # - NFaces -
    B = TTK.Label(Frame, text="NFaces")
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of faces in selection.')
    B = TTK.Entry(Frame, textvariable=VARS[7], background='White',
                  width=10)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of faces in selection.')

    # - Variable -
    B = TTK.Label(Frame, text="Variable: ")
    B.grid(row=4, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Selected var name.')
    F = TTK.Frame(Frame, borderwidth=0)
    F.columnconfigure(0, weight=1)

    if ttk is None:
        B = TK.OptionMenu(F, VARS[5], '')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList)
        F.grid(row=4, column=1, columnspan=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Selected var name.')
        WIDGETS['variable'] = B
    else: 
        B = ttk.Combobox(F, textvariable=VARS[5], 
                         values=[], state='readonly')
        B.grid(sticky=TK.EW)
        F.bind('<Enter>', updateVarNameList2)
        F.grid(row=4, column=1, columnspan=1, sticky=TK.EW)
        BB = CTK.infoBulle(parent=B, text='Selected var name.')
        WIDGETS['variable'] = B
    
    # - Min de la variable - 
    B = TTK.Label(Frame, text="Min:")
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Min of variable.')
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White',
                  width=10)
    B.grid(row=5, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Min of variable.')
    
    # - Max de la variable -    
    B = TTK.Label(Frame, text="Max:")
    B.grid(row=6, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Max of variable.')
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White',
                  width=10)
    B.grid(row=6, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Max of variable.')
    
    # - update -
    B = TTK.Button(Frame, text="Update info", command=updateInfo)
    B.grid(row=7, column=0, columnspan=2, sticky=TK.EW)
    
#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.EW)

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

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
    (win, menu, file, tools) = CTK.minimal('tkMeshInfo '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
