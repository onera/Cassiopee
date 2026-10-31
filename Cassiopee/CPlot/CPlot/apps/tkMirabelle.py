# - tkMirabelle -
"""Remesh structured meshes."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Transform.PyTree as T
import Apps.Mesh.Mirabelle3 as Mirabelle
import Geom.PyTree as D

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def extractEdge(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if len(nzs) > 1:
        CTK.TXT.insert('START', 'Please select only one block.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # selected block
    nob = CTK.Nb[nzs[0]]+1
    noz = CTK.Nz[nzs[0]]
    z = CTK.t[2][nob][2][noz]
    dimz = Internal.getZoneDim(z)
    if dimz[0] == 'Unstructured':
        CTK.TXT.insert('START', 'Selected block must be structured.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    ni = dimz[1]; nj = dimz[2]; nk = dimz[3]

    # selected range
    ind = CPlot.getActivePointIndex()
    diri = (ind[2] == 1) or (ind[2] == ni)
    dirj = (ind[3] == 1) or (ind[3] == nj)
    dirk = (ind[4] == 1) or (ind[4] == nk)
    if diri and dirj and not dirk:
        rdir = 5
        erange = [ind[2],ind[2],ind[3],ind[3],1,nk]
    elif diri and not dirj and dirk:
        rdir = 3
        erange = [ind[2],ind[2],1,nj,ind[4],ind[4]]
    elif not diri and dirj and dirk:
        rdir = 1
        erange = [1,ni,ind[3],ind[3],ind[4],ind[4]]
    else:
        CTK.TXT.insert('START', 'Not a valid edge.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.saveTree()

    print('indmin=',erange[0],erange[2],erange[4])
    print('indmax=',erange[1],erange[3],erange[5])
    zt = T.subzone(z, (erange[0],erange[2],erange[4]), (erange[1],erange[3],erange[5]))
    Internal.newUserDefinedData('block', value=z[0], parent=zt)
    Internal.newUserDefinedData('rdir', value=rdir, parent=zt)

    b = Internal.createUniqueChild(CTK.t, 'EXTRACTED', 'CGNSBase_t')
    nob = C.getNobOfBase(b, CTK.t)
    CTK.add(CTK.t, nob, -1, zt)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def propagateEdge(event=None):

    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if len(nzs) > 1:
        CTK.TXT.insert('START', 'Please select only one block.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # selected linelet
    nob = CTK.Nb[nzs[0]]+1
    noz = CTK.Nz[nzs[0]]
    linelet = CTK.t[2][nob][2][noz]
    # remove linelet from CTK.t
    CTK.t[2][nob][2].pop(noz)

    block = Internal.getNodeFromName1(linelet, 'block')
    rdir = Internal.getNodeFromName1(linelet, 'rdir')
    if block is None or rdir is None:
        CTK.TXT.insert('START', 'Edge must have been extracted with extract edge.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    block = Internal.getValue(block)
    rdir = Internal.getValue(rdir)

    # passe la linelet en abscisse curviligne
    linelet = D.getCurvilinearAbscissa(linelet)
    C._initVars(linelet, '{CoordinateX}={s}')
    C._initVars(linelet, '{CoordinateY}=0.')
    C._initVars(linelet, '{CoordinateZ}=0.')

    # getting match graph
    g = Mirabelle.getGraph(CTK.t)

    stack = [(block, rdir)]

    print("== Running propagate...")
    treated = []
    Mirabelle._propagate(CTK.t, g, stack, treated, linelet)

    # Adapte les donneurs a la fin
    Mirabelle._adaptDonorRanges(CTK.t)

    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.display(CTK.t)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkMirabelle  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Structured remeshing.\nCtrl+w to close applet.', temps=0, btype=1)
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
    CTK.addPinMenu(FrameMenu, 'tkMirabelle')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Zone filter regexp -
    #V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Extract edge -
    B = TTK.Button(Frame, text="Extract edge", command=extractEdge)
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Extract mesh edge for remeshing.')

    # - Propagate edge -
    B = TTK.Button(Frame, text="Propagate edge", command=propagateEdge)
    B.grid(row=1, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Propagate remeshed edge in structured mesh.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.NSEW)

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
    (win, menu, file, tools) = CTK.minimal('tkMirabelle '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
