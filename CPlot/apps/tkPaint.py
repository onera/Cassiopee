# - paint app -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def paint():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    # get patined var
    field = CPlot.getState('scalarField')
    if field == -1:
        CTK.TXT.insert('START', 'Scalar field is not set.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CPlot.unselectAllZones()

    CTK.saveTree()
    w = WIDGETS['paint']
    if not CTK.__BUSY__:
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)

        while CTK.__BUSY__ == True:
            l = []
            while l == []:
                nz = CPlot.getSelectedZone()
                l = CPlot.getActivePointIndex()
                #time.sleep(CPlot.__timeStep__)
                w.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                nob = CTK.Nb[nz]+1
                noz = CTK.Nz[nz]
                CTK.saveTree()
                value = float(WIDGETS['value'].get())
                #width = float(WIDGETS['width'].get())
                #brushType = VARS[2].get()
                z = CTK.t[2][nob][2][noz]
                posCam = CPlot.getState('posCam')
                posEye = CPlot.getState('posEye')
                vect = (posEye[0]-posCam[0], posEye[1]-posCam[1],
                        posEye[2]-posCam[2])

                click = CPlot.getActivePoint()
                point = (click[0], click[1], click[2])
                ind = CPlot.getActivePointIndex()
                #hook = C.createHook
                varNames = C.getVarNames(z)[0]
                if len(varNames) > field+3:
                    var = varNames[field+3]
                    C.setValue(z, var, ind[0], value)
                else: 
                    CTK.TXT.insert('START', 'Field %d not found. Use scalar mode.\n'%field)
                    CTK.TXT.insert('START', 'Error: ', 'Error')
                CTK.replace(CTK.t, nob, noz, z)

                CTK.TKTREE.updateApp()
                CPlot.unselectAllZones()
                CPlot.render()
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)
    else:
       CTK.__BUSY__ = False
       TTK.raiseButton(w)
       CPlot.setState(cursor=0)

#==============================================================================
# Get value from mouse pointer
#==============================================================================
def getValue():
    nz = CPlot.getSelectedZone()
    l = CPlot.getActivePointIndex()
    if l == []: return
    field = CPlot.getState('scalarField')
    if field == -1: return
    ind = CPlot.getActivePointIndex()
    nob = CTK.Nb[nz]+1
    noz = CTK.Nz[nz]
    z = CTK.t[2][nob][2][noz]
    varNames = C.getVarNames(z)[0]
    if len(varNames) > field+3:
        var = varNames[field+3]
        val = C.getValue(z, var, ind[0])
        VARS[0].set(str(val))
    else:
        CTK.TXT.insert('START', 'Field %d not found. Use scalar mode.\n'%field)
        CTK.TXT.insert('START', 'Error: ', 'Error')
    #WIDGETS['value'].set(str(val))

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkPaint', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Paint fields.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=2)
    Frame.columnconfigure(1, weight=0)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkPaint')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- value -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -1- width -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -2- Brush
    V = TK.StringVar(win); V.set('Sphere'); VARS.append(V)
    
    # - Value -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Value.')
    WIDGETS['value'] = B
    B = TTK.Button(Frame, command=getValue,
                   image=iconics.PHOTO[8], padx=0, pady=0)
    B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Get value from mouse pointer.')

    # - Width -
    #B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    #B.grid(row=0, column=1, columnspan=1, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Width.')
    #WIDGETS['width'] = B

    # - Brush -
    #B = TTK.OptionMenu(Frame, VARS[2], 'Sphere')
    #BB = CTK.infoBulle(parent=B, text='Brush shape.')
    #B.grid(row=1, column=0, sticky=TK.EW)

    # - Paint mode -
    B = TTK.Button(Frame, text="Paint mode", command=paint)
    BB = CTK.infoBulle(parent=B, text='Enter paint mode.\nValue of field (scalar mode) is modified.')
    B.grid(row=0, column=2, columnspan=1, sticky=TK.EW)
    WIDGETS['paint'] = B
    
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
    (win, menu, file, tools) = CTK.minimal('tkPaint '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
