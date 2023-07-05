# - tkFilterSurfs -
"""Filter or offset a surface."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Geom.PyTree 
# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def remap():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    a = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        a.append(z)

    # density
    density = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(density) != 1:
        CTK.TXT.insert('START', 'Density is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    pointsPerUnitLength = density[0]
    
    # offset
    offset = CTK.varsFromWidget(VARS[1].get(), type=1)
    if len(offset) != 1:
        CTK.TXT.insert('START', 'Offset is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    offset = offset[0]

    CTK.saveTree()
    CPlot.setState(cursor=2)
    if VARS[2].get() == '0': algo = 0
    else: algo = 1
    iso = Geom.PyTree.offsetSurface(a, offset, pointsPerUnitLength, algo)
    CPlot.setState(cursor=0)

    if iso != []:
        nob = CTK.Nb[nzs[0]]+1
        for i in iso: CTK.add(CTK.t, nob, -1, i)
    
        #C._fillMissingVariables(CTK.t)
        CTK.TXT.insert('START', 'Surface filtered and offset (offset=%g).\n'%offset)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
    else:
         CTK.TXT.insert('START', 'Surface filter failed.\n')
         CTK.TXT.insert('START', 'Error: ', 'Error')     
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFilterSurfs  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Filter or offset a surface.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=2)
    Frame.columnconfigure(3, weight=0)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkFilterSurfs')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Point density -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkFilterSurfsDensity' in CTK.PREFS: 
        V.set(CTK.PREFS['tkFilterSurfsDensity'])
    # -1- Offset -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    if 'tkFilterSurfsOffset' in CTK.PREFS: 
        V.set(CTK.PREFS['tkFilterSurfsOffset'])
    # -2- Algorithm cart/octree
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkFilterSurfsType' in CTK.PREFS: 
        V.set(CTK.PREFS['tkFilterSurfsType'])

    # - Point density -
    B = TTK.Label(Frame, text="density/offset")
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=8)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Point density (npts/length unit).')

    # - Offset -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=8)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Distance of surface offset.')

    # - Algorithm -
    B = TTK.Checkbutton(Frame, text='', variable=VARS[2])
    B.grid(row=0, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Toggle octree algorithm.')
    
    # - Remap -
    B = TTK.Button(Frame, text="Filter/offset", command=remap)
    BB = CTK.infoBulle(parent=B, text='Remap/offset a surface with a given spacing.')
    B.grid(row=2, column=0, columnspan=4, sticky=TK.EW)
    
#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkFilterSurfs')
    CTK.WIDGETS['SurfNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['SurfNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkFilterSurfsDensity'] = VARS[0].get()
    CTK.PREFS['tkFilterSurfsOffset'] = VARS[1].get()
    CTK.PREFS['tkFilterSurfsType'] = VARS[2].get()
    CTK.savePrefFile()
    
#==============================================================================
def resetApp():
    VARS[0].set('1.')
    VARS[1].set('0.')
    VARS[2].set('0')
    CTK.PREFS['tkFilterSurfsDensity'] = VARS[0].get()
    CTK.PREFS['tkFilterSurfsOffset'] = VARS[1].get()
    CTK.PREFS['tkFilterSurfsType'] = VARS[1].get()
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
    (win, menu, file, tools) = CTK.minimal('tkFilterSurfs '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
