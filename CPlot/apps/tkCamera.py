# - display Camera information -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

def getInfo(event=None):
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')
    VARS[0].set('(%f,% f, %f)'%(posCam[0],posCam[1],posCam[2]))
    VARS[1].set('(%f, %f, %f)'%(posEye[0],posEye[1],posEye[2]))
    VARS[2].set('(%f, %f, %f)'%(dirCam[0],dirCam[1],dirCam[2]))

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

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkCamera', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Display mesh informations.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=3)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkCamera')
    WIDGETS['frameMenu'] = FrameMenu

    #- VARS -
    # -0- posCam -
    V = TK.StringVar(win); V.set('(0., 0., 0.)'); VARS.append(V)
    # -1- posEye -
    V = TK.StringVar(win); V.set('(0., 0.; 0.)'); VARS.append(V)
    # -2- dirCam -
    V = TK.StringVar(win); V.set('(0., 0., 0.)'); VARS.append(V)
    
    # - posCam -
    B = TTK.Label(Frame, text="posCam: ")
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=15)
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)
    B.bind('<Return>', setInfo)
    BB = CTK.infoBulle(parent=B, text='Camera position.')

    # - posEye -
    B = TTK.Label(Frame, text="posEye: ")
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=15)
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    B.bind('<Return>', setInfo)
    BB = CTK.infoBulle(parent=B, text='Eye position.')

    # - dirCam -
    B = TTK.Label(Frame, text="dirCam: ")
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=15)
    B.grid(row=2, column=1, columnspan=1, sticky=TK.EW)
    B.bind('<Return>', setInfo)
    BB = CTK.infoBulle(parent=B, text='Camera direction.')
    
    # - get -
    B = TTK.Button(Frame, text="Get", command=getInfo)
    B.grid(row=3, column=1, columnspan=1, sticky=TK.EW)
    
    # - set -
    B = TTK.Button(Frame, text="Set", command=setInfo)
    B.grid(row=3, column=0, columnspan=1, sticky=TK.EW)
    
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
    (win, menu, file, tools) = CTK.minimal('tkCamera '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
