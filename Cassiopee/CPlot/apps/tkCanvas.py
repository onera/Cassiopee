# - tkCanvas -
"""Create drawing canvas."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot as CP
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

# Globals
CANVASSIZE = -1.
XC = 0; YC = 0; ZC = 0
UX = 0; UY = 0; UZ = 0
LX = 0; LY = 0; LZ = 0
DIRX = 0; DIRY = 0; DIRZ = 0

#==============================================================================
def initCanvas(event=None):
    zdir = VARS[0].get()
    if zdir == 'None' and CTK.t == []: return
    global CANVASSIZE, XC, YC, ZC, UX, UY, UZ, LX, LY, LZ, DIRX, DIRY, DIRZ
    if CANVASSIZE == -1:
        if CTK.t == []:
            CANVASSIZE = 1
            dx = 1
            VARS[2].set(str(dx))
            VARS[3].set(str(dx))
            XC = 0.
            YC = 0.
            ZC = 0.
        else:
            bb = G.bbox(CTK.t)
            CANVASSIZE = max(bb[3]-bb[0], bb[4]-bb[1], bb[5]-bb[2])
            dx = (bb[3]-bb[0] + bb[4]-bb[1] + bb[5]-bb[2])/30.
            VARS[2].set(str(dx))
            VARS[3].set(str(dx))
            XC = 0.5*(bb[3]+bb[0])
            YC = 0.5*(bb[4]+bb[1])
            ZC = 0.5*(bb[5]+bb[2])

    nzs = CPlot.getSelectedZones()
    if nzs != [] and dir != 'None':
        point = CPlot.getActivePoint()
        if point != []:
            if zdir == 'YZ': XC = point[0]
            elif zdir == 'XZ': YC = point[1]
            elif zdir == 'XY': ZC = point[2]
            else: XC = point[0]; YC = point[1]; ZC = point[2]

    if zdir == 'None':
        deleteCanvasBase()
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render(); return
    elif zdir == 'View':
        posCam = CPlot.getState('posCam')
        posEye = CPlot.getState('posEye')
        dirCam = CPlot.getState('dirCam')
        if nzs == []:
            XC = posEye[0]; YC = posEye[1]; ZC = posEye[2]
        LX = posEye[0]-posCam[0]
        LY = posEye[1]-posCam[1]
        LZ = posEye[2]-posCam[2]
        if LX*LX + LY*LY + LZ*LZ < 1.e-10: LX = -1
        if dirCam[0]*dirCam[0] + dirCam[1]*dirCam[1] + dirCam[2]*dirCam[2] == 0.:
            dirCam = (0,0,1)
        DIRX = dirCam[0]; DIRY = dirCam[1]; DIRZ = dirCam[2]

        UX = dirCam[1]*LZ - dirCam[2]*LY
        UY = dirCam[2]*LX - dirCam[0]*LZ
        UZ = dirCam[0]*LY - dirCam[1]*LX

    elif zdir == 'YZ':
        DIRX = 1; DIRY = 0; DIRZ = 0
        UX = 0; UY = 1; UZ = 0
        LX = 0; LY = 0; LZ = 1
    elif zdir == 'XZ':
        DIRX = 0; DIRY = 1; DIRZ = 0
        UX = 1; UY = 0; UZ = 0
        LX = 0; LY = 0; LZ = 1
    else:
        DIRX = 0; DIRY = 0; DIRZ = 1
        UX = 1; UY = 0; UZ = 0
        LX = 0; LY = 0; LZ = 1
    CTK.TXT.insert('START', 'Set a canvas.\n')
    setCanvas()

#==============================================================================
def deleteCanvasBase():
    nodes = Internal.getNodesFromName1(CTK.t, 'CANVAS')
    if nodes == []: return
    canvas = nodes[0]
    # delete from plotter
    zones = Internal.getNodesFromType(canvas, 'Zone_t')
    dels = []
    for z in zones: dels.append(canvas[0]+Internal.SEP1+z[0])
    CPlot.delete(dels)
    ret = Internal.getParentOfNode(CTK.t, canvas)
    del ret[0][2][ret[1]]

#==============================================================================
def setCanvas(event=None):
    zdir = VARS[0].get()
    CTK.saveTree()
    deleteCanvasBase()
    CTK.t = C.addBase2PyTree(CTK.t, 'CANVAS', 2)

    size = CANVASSIZE
    if zdir == 'View':
        a = G.cart( (-size, 0, -size), (2*size, 2*size, 2*size), (2,1,2) )
        a = T.translate(a, (XC, YC, ZC) )
        a = T.rotate(a, (XC, YC, ZC),
                     ((1,0,0),(0,1,0),(0,0,1)),
                     ((-DIRX,-DIRY,-DIRZ), (LX,LY,LZ), (UX,UY,UZ)))
        VARS[1].set(str(XC)+';'+str(YC)+';'+str(ZC))
    elif zdir == 'YZ':
        a = G.cart( (0, -size, -size), (2*size, 2*size, 2*size), (1,2,2) )
        a = T.translate(a, (XC, YC, ZC) )
        VARS[1].set(str(XC))
    elif zdir == 'XZ':
        a = G.cart( (-size, 0, -size), (2*size, 2*size, 2*size), (2,1,2) )
        a = T.translate(a, (XC, YC, ZC) )
        VARS[1].set(str(YC))
    else:
        a = G.cart( (-size, -size, 0), (2*size, 2*size, 2*size), (2,2,1) )
        a = T.translate(a, (XC, YC, ZC) )
        VARS[1].set(str(ZC))
    nodes = Internal.getNodesFromName1(CTK.t, 'CANVAS')
    base = nodes[0]
    nob = C.getNobOfBase(base, CTK.t)
    if CP.__slot__ is None:
        CTK.t[2][nob][2].append(a); CTK.display(CTK.t)
    else: CTK.add(CTK.t, nob, -1, a)

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def moveUp():
    if CTK.t == []: return
    zdir = VARS[0].get()
    step = VARS[2].get()
    try: step = float(step)
    except: step = 1.

    global XC, YC, ZC
    if zdir == 'View':
        XC = XC + step*LX
        YC = YC + step*LY
        ZC = ZC + step*LZ
    elif zdir == 'YZ':
        XC = XC + step
    elif zdir == 'XZ':
        YC = YC + step
    else:
        ZC = ZC + step
    CTK.TXT.insert('START', 'Move canvas up.\n')
    setCanvas()

#==============================================================================
def moveDown():
    if CTK.t == []: return
    zdir = VARS[0].get()
    step = VARS[2].get()
    try: step = float(step)
    except: step = 1.

    global XC, YC, ZC
    if zdir == 'View':
        XC = XC - step*LX
        YC = YC - step*LY
        ZC = ZC - step*LZ
    elif zdir == 'YZ':
        XC = XC - step
    elif zdir == 'XZ':
        YC = YC - step
    else:
        ZC = ZC - step
    CTK.TXT.insert('START', 'Move canvas down.\n')
    setCanvas()

#==============================================================================
def enlarge():
    global CANVASSIZE
    if CTK.t == []: return
    step = VARS[3].get()
    try: step = float(step)
    except: step = 1.
    CANVASSIZE = CANVASSIZE + step
    CTK.TXT.insert('START', 'Enlarge canvas.\n')
    setCanvas()

#==============================================================================
def reduce():
    global CANVASSIZE
    if CTK.t == []: return
    step = VARS[3].get()
    try: step = float(step)
    except: step = 1.
    CANVASSIZE = CANVASSIZE - step
    if CANVASSIZE < 1.e-3: CANVASSIZE = 1.e-3
    CTK.TXT.insert('START', 'Reduce canvas.\n')
    setCanvas()

#==============================================================================
def setCanvasPos(event=None):
    if CTK.t == []: return
    global XC, YC, ZC
    zdir = VARS[0].get()
    res = CTK.varsFromWidget(VARS[1].get(), type=1)
    if zdir == 'View':
        if len(res) != 3:
            CTK.TXT.insert('START', 'xc;yc;zc is incorrect.\n'); return
        else:
            XC = res[0]; YC = res[1]; ZC = res[2]
    else:
        if len(res) != 1:
            CTK.TXT.insert('START', 'pos is incorrect.\n'); return
        else:
            if zdir == 'YZ': XC = res[0]
            elif zdir == 'XZ': YC = res[0]
            else: ZC = res[0]
    setCanvas()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCanvas  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Manage a canvas\nfor drawing..\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=4)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkCanvas')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Canvas dir
    V = TK.StringVar(win); V.set('None'); VARS.append(V)
    # -1- Canvas pos
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -2- Canvas dir step
    V = TK.StringVar(win); V.set('0.1'); VARS.append(V)
    if 'tkCanvasDirStep' in CTK.PREFS:
        V.set(CTK.PREFS['tkCanvasDirStep'])
    # -3- Canvas enlarge step
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkCanvasEnlargeStep' in CTK.PREFS:
        V.set(CTK.PREFS['tkCanvasEnlargeStep'])

    # - Direction -
    B = TTK.OptionMenu(Frame, VARS[0],
                       'None', 'View', 'XY', 'XZ', 'YZ',
                       command=initCanvas)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create a canvas.')

    # - Refresh -
    B = TTK.Button(Frame, text="", command=initCanvas,
                   image=iconics.PHOTO[8], padx=0, pady=2)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Refresh canvas.\nUpdate it to make it pass through selection point.')

    # - Position -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=7)
    B.bind('<Return>', setCanvasPos)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Canvas position.')

    # - Move canvas -
    B = TTK.Button(Frame, text="+", command=moveUp)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Move canvas right.')
    B = TTK.Button(Frame, text="-", command=moveDown)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Move canvas left.')
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=7)
    B.grid(row=1, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Canvas move step.')

    # - Enlarge / reduce -
    B = TTK.Button(Frame, text="Enlarge", command=enlarge)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Enlarge canvas.')
    B = TTK.Button(Frame, text="Reduce", command=reduce)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Reduce canvas.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=7)
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Canvas enlarge step.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['EdgeNoteBook'].add(WIDGETS['frame'], text='tkCanvas')
    except: pass
    CTK.WIDGETS['EdgeNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['EdgeNoteBook'].hide(WIDGETS['frame'])

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(): return

#==============================================================================
def saveApp():
    CTK.PREFS['tkCanvasDirStep'] = VARS[2].get()
    CTK.PREFS['tkCanvasEnlargeStep'] = VARS[3].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[2].set('0.1')
    VARS[3].set('1.')
    CTK.PREFS['tkCanvasDirStep'] = VARS[2].get()
    CTK.PREFS['tkCanvasEnlargeStep'] = VARS[3].get()
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
    (win, menu, file, tools) = CTK.minimal('tkCanvas '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
