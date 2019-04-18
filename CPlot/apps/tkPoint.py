# - an applet for creating points -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Geom.PyTree as D
import Converter.Internal as Internal
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def getPointCoordinates():
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    point = CPlot.getActivePoint()
    if point != []:
        VARS[0].set(str(point[0])+';'+str(point[1])+';'+str(point[2]))

#==============================================================================
def createPoint(event=None):
    point = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(point) != 3:
        CTK.TXT.insert('START', 'Point coords are incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'POINTS', 1)
    base = Internal.getNodesFromName1(CTK.t, 'POINTS')
    base = base[0]
    nob = C.getNobOfBase(base, CTK.t)
    a = D.point((point[0], point[1], point[2]))
    CTK.add(CTK.t, nob, -1, a)
    CTK.TXT.insert('START', 'Point '+VARS[0].get()+' created.\n')
    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def modifyPoint(event=None):
    if CTK.t == []: return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
       CTK.TXT.insert('START', 'Selection is empty.\n')
       CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if len(nzs) > 1:
       CTK.TXT.insert('START', 'Please, select only one zone.\n')
       CTK.TXT.insert('START', 'Error: ', 'Error'); return
    point = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(point) != 3:
        CTK.TXT.insert('START', 'Point coords are incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
       
    nz = nzs[0]
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    ind = CPlot.getActivePointIndex()
    if ind == []:
       CTK.TXT.insert('START', 'No selected point.\n')    
       CTK.TXT.insert('START', 'Error: ', 'Error'); return
    ind = ind[0]
    
    CTK.saveTree()
    zp = Internal.copyTree(z)
    C.setValue(zp, 'CoordinateX', ind, point[0])
    C.setValue(zp, 'CoordinateY', ind, point[1])
    C.setValue(zp, 'CoordinateZ', ind, point[2])
    
    CPlot.replace(CTK.t, nob, noz, zp)
    CTK.TXT.insert('START', 'Point modified.\n')
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkPoint', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Generate points.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=4)
    Frame.columnconfigure(1, weight=4)
    Frame.columnconfigure(2, weight=0)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkPoint')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- point coordinates -
    V = TK.StringVar(win); V.set('0;0;0'); VARS.append(V)
    
    # - Buttons -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=15)
    B.bind('<Return>', createPoint)
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Point coordinates (x;y;z).')
    B = TTK.Button(Frame, image=iconics.PHOTO[8],
                   command=getPointCoordinates, padx=0)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get point coordinates from mouse.')
    
    B = TTK.Button(Frame, text="Create", command=createPoint)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create a point.')

    B = TTK.Button(Frame, text="Modify", command=modifyPoint)
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Modify point coordinate.')

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
    (win, menu, file, tools) = CTK.minimal('tkPoint '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
