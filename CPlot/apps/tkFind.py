# - find nodes and elements of given number -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Converter.Internal as Internal
import Transform.PyTree as T
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def find(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if len(nzs) != 1:
        CTK.TXT.insert('START', 'Only one block must be selected.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    type = VARS[0].get()
    if type == 'Node index':
        inds = CTK.varsFromWidget(VARS[1].get(), type=3)
        if len(inds) == 0:
            CTK.TXT.insert('START', 'Invalid index.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        [px,py,pz] = C.getValue(z, Internal.__GridCoordinates__, inds[0])
        CTK.TXT.insert('START', 'Node %s found.\n'%VARS[1].get())
    elif type == 'Coordinates':
        vars = CTK.varsFromWidget(VARS[1].get(), type=1)
        if len(vars) != 3:
            CTK.TXT.insert('START', 'Invalid coordinates.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        [px,py,pz] = vars
        CTK.TXT.insert('START', 'Coordinates %s found.\n'%VARS[1].get())
    else: # type = element index
        inds = CTK.varsFromWidget(VARS[1].get(), type=3)
        if len(inds) == 0:
            CTK.TXT.insert('START', 'Invalid index.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        zp = Internal.copyRef(z)
        C._deleteAllBCAndSolutions__(zp)
        zp = C.node2Center(z)
        [px,py,pz] = C.getValue(zp, Internal.__GridCoordinates__, inds[0])
        CTK.TXT.insert('START', 'Element %s found.\n'%VARS[1].get())

    (xeye,yeye,zeye) = CPlot.getState('posEye')
    (xcam,ycam,zcam) = CPlot.getState('posCam')
    dx = xcam-xeye; dy = ycam-yeye; dz = zcam-zeye
    CPlot.setState(posEye=(px,py,pz))
    CPlot.setState(posCam=(px+dx,py+dy,pz+dz))
    CPlot.setActivePoint(px,py,pz)

#==============================================================================
def extract(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if len(nzs) != 1:
        CTK.TXT.insert('START', 'Only one block must be selected.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    type = VARS[0].get()
    if type == 'Node index':
        CTK.TXT.insert('START', 'Nodes are invalid for extraction.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        inds = CTK.varsFromWidget(VARS[1].get(), type=3)
        if len(inds) == 0:
            CTK.TXT.insert('START', 'Invalid index.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        ind = inds[0]
        [px,py,pz] = C.getValue(z, Internal.__GridCoordinates__, ind)
        #try:
        a = T.subzone(z, inds)
        b = Internal.createUniqueChild(CTK.t, 'EXTRACT', 'CGNSBase_t')
        nob = C.getNobOfBase(b, CTK.t)
        CTK.add(CTK.t, nob, -1, a)
        #except: pass
        #CTK.TXT.insert('START', 'Extraction not implemented for Nodes.\n')
        #CTK.TXT.insert('START', 'Error: ', 'Error')
    elif type == 'Coordinates':
        CTK.TXT.insert('START', 'Coordinates are invalid for extraction.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        [px,py,pz] = vars
    else: # type = element index
        inds = CTK.varsFromWidget(VARS[1].get(), type=3)
        if len(inds) == 0:
            CTK.TXT.insert('START', 'Invalid index.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error'); return
        ind = inds[0]
        try:
            if Internal.getZoneType(z) == 1: # struct
                dims = Internal.getZoneDim(z)
                indices = []
                for ind in inds: # convert to monoindex
                    if isinstance(ind, tuple):
                        ind = ind[0]-1+(ind[1]-1)*dims[1]+(ind[2]-1)*dims[1]*dims[2]
                        indices.append(ind)
                    else: indices.append(ind)
                a = T.subzone(C.convertArray2Hexa(z), indices, type='elements')
            else:  # unstruct
                a = T.subzone(z, inds, type='elements')
            b = Internal.createUniqueChild(CTK.t, 'EXTRACT', 'CGNSBase_t')
            nob = C.getNobOfBase(b, CTK.t) 
            CTK.add(CTK.t, nob, -1, a)
        except: pass
        #zp = C.node2Center(z)
        #[px,py,pz] = C.getValue(zp, Internal.__GridCoordinates__, ind)

    #(xeye,yeye,zeye) = CPlot.getState('posEye')
    #(xcam,ycam,zcam) = CPlot.getState('posCam')
    #dx = xcam-xeye; dy = ycam-yeye; dz = zcam-zeye
    #CPlot.setState(posEye=(px,py,pz))
    #CPlot.setState(posCam=(px-0.1*dx,py-0.1*dy,pz-0.1*dz))
    #CPlot.setActivePoint(px,py,pz)
    CTK.TXT.insert('START', 'Cell extracted.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFind', font=CTK.FRAMEFONT, takefocus=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=4)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkFind')
    WIDGETS['frameMenu'] = FrameMenu
    
    # - VARS -
    # -0- type (node/elements/coordinates) -
    V = TK.StringVar(win); V.set('Node index'); VARS.append(V)
    if 'tkFindDataType' in CTK.PREFS: V.set(CTK.PREFS['tkFindDataType'])
    # -1- value -
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    
    # - Value -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White')
    B.bind('<Return>', find)
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Index to find/extract in selected mesh.\nCan be a global index (0 starts)\nor (i,j,k) (1 starts) for structured grids.\nFor extraction, you can define a list of indices.')

    # - type -
    B = TTK.OptionMenu(Frame, VARS[0], 'Node index', 'Element index', 'Coordinates')
    BB = CTK.infoBulle(parent=B, text='Type of data to be searched.')
    B.grid(row=0, column=0, sticky=TK.EW)

    # - Find -
    B = TTK.Button(Frame, text="Find", command=find)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Find the given index/coordinates in selected mesh.')

    # - Extract -
    B = TTK.Button(Frame, text="Extract", command=extract)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract the given index/coordinates cells.')

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
def saveApp():
    CTK.PREFS['tkFindDataType'] = VARS[0].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('Node index')
    CTK.PREFS['tkFindDataType'] = VARS[0].get()
    CTK.savePrefFile()
    
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
    (win, menu, file, tools) = CTK.minimal('tkPersonal '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
