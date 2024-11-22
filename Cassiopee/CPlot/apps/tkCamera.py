# - tkCamera -
"""Display Camera information."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Get camera info from state
#==============================================================================
def getInfo(event=None):
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')
    VARS[0].set('(%f,% f, %f)'%(posCam[0],posCam[1],posCam[2]))
    VARS[1].set('(%f, %f, %f)'%(posEye[0],posEye[1],posEye[2]))
    VARS[2].set('(%f, %f, %f)'%(dirCam[0],dirCam[1],dirCam[2]))

#==============================================================================
# Set camera position
#==============================================================================
def setInfo(event=None):
    posCam = VARS[0].get()
    posCam = CTK.varsFromWidget(posCam, 1)
    if posCam != [] and len(posCam) == 3:
        CPlot.setState(posCam=posCam)
    posEye = VARS[1].get()
    posEye = CTK.varsFromWidget(posEye, 1)
    if posEye != [] and len(posEye) == 3:
        CPlot.setState(posEye=posEye)
    dirCam = VARS[2].get()
    dirCam = CTK.varsFromWidget(dirCam, 1)
    if dirCam != [] and len(dirCam) == 3:
        CPlot.setState(dirCam=dirCam)
    CTK.TXT.insert('START', 'Set camera position.\n')

#==============================================================================
# Write CPlot command
#==============================================================================
def exportInfo(event=None):
    com = 'CPlot.display( t, '
    posCam = CPlot.getState('posCam'); com += 'posCam=%s'%str(posCam)
    posEye = CPlot.getState('posEye'); com += ', posEye=%s'%str(posEye)
    dirCam = CPlot.getState('dirCam'); com += ', dirCam=%s'%str(dirCam)
    mode = CPlot.getState('mode')
    if mode == 0: com += ', mode="Mesh"'
    elif mode == 1: com += ', mode="Solid"'
    elif mode == 2: com += ', mode="Render"'
    elif mode == 3: com += ', mode="Scalar"'
    if mode == 0: # mesh
        meshStyle = CPlot.getState('meshStyle'); com += ', meshStyle=%d'%meshStyle
    elif mode == 1: # solid
        solidStyle = CPlot.getState('solidStyle'); com += ', solidStyle=%d'%solidStyle
    elif mode == 3: # scalar
        scalarField = CPlot.getState('scalarField')
        varnames = C.getVarNames(CTK.t, excludeXYZ=True)
        field = str(scalarField)
        if len(varnames) > 0:
            varnames = varnames[0]
            if scalarField < len(varnames): field = varnames[scalarField]
        com += ', scalarField="%s"'%field

        scalarStyle = CPlot.getState('scalarStyle'); com += ', scalarStyle=%d'%scalarStyle
        colormap = CPlot.getState('colormap'); com += ', colormap=%d'%colormap
        isoEdges = CPlot.getState('isoEdges'); com += ', isoEdges=%g'%isoEdges
        isoScale = CPlot.getState('isoScale'); 
        isoScale = [field]+isoScale
        com += ', isoScales=%s'%str(isoScale)

    #bgColor = CPlot.getState('bgColor'); com += ', bgColor=%d'%bgColor
    com += ' )'
    # Ecriture a l'ecran
    print('\n')
    print(com)
    print('\n')

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCamera  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Display mesh informations.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkCamera')
    WIDGETS['frameMenu'] = FrameMenu

    #- VARS -
    # -0- posCam -
    V = TK.StringVar(win); V.set('(0., 0., 0.)'); VARS.append(V)
    # -1- posEye -
    V = TK.StringVar(win); V.set('(0., 0., 0.)'); VARS.append(V)
    # -2- dirCam -
    V = TK.StringVar(win); V.set('(0., 0., 0.)'); VARS.append(V)

    # - posCam -
    B = TTK.Label(Frame, text="posCam: ")
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=15)
    B.grid(row=0, column=1, columnspan=2, sticky=TK.EW)
    B.bind('<Return>', setInfo)
    BB = CTK.infoBulle(parent=B, text='Camera position.')

    # - posEye -
    B = TTK.Label(Frame, text="posEye: ")
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=15)
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    B.bind('<Return>', setInfo)
    BB = CTK.infoBulle(parent=B, text='Eye position.')

    # - dirCam -
    B = TTK.Label(Frame, text="dirCam: ")
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=15)
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
    B.bind('<Return>', setInfo)
    BB = CTK.infoBulle(parent=B, text='Camera direction.')

    # - get -
    B = TTK.Button(Frame, text="Get", command=getInfo)
    B.grid(row=3, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get camera position from view.')

    # - set -
    B = TTK.Button(Frame, text="Set", command=setInfo)
    B.grid(row=3, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set Camera position in view.')

    # - export -
    B = TTK.Button(Frame, text="Export", command=exportInfo)
    B.grid(row=3, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Write display command line to terminal.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['StateNoteBook'].add(WIDGETS['frame'], text='tkCamera')
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
    (win, menu, file, tools) = CTK.minimal('tkCamera '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
