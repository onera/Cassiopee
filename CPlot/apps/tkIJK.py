# - show only subzones of structured grids -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Transform.PyTree as T
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

def setPlusI(event=None):
    set__(VARS[0], VARS[1], +1)

def setMoinsI(event=None):
    set__(VARS[0], VARS[1], -1)
    
def setFullI(event=None):
    setFull__(VARS[0], VARS[1])
    
def setPlusJ(event=None):
    set__(VARS[2], VARS[3], +1)

def setMoinsJ(event=None):
    set__(VARS[2], VARS[3], -1)

def setFullJ(event=None):
    setFull__(VARS[2], VARS[3])

def setPlusK(event=None):
    set__(VARS[4], VARS[5], +1)

def setMoinsK(event=None):
    set__(VARS[4], VARS[5], -1)

def setFullK(event=None):
    setFull__(VARS[4], VARS[5])

def set__(V0, V1, off):
    v = CTK.varsFromWidget(V0.get(), type=2)
    if len(v) > 0: imin = v[0]
    v = CTK.varsFromWidget(V1.get(), type=2)
    if len(v) > 0: imax = v[0]
    m = min(imin, imax)+off
    m = max(m,1)
    V0.set(str(m))
    V1.set(str(m))
    show()

def setFull__(V0, V1):
    V0.set(str(1))
    V1.set(str(-1))
    show()
       
#==============================================================================
# show
#==============================================================================
def show(event=None):
    if CTK.t == []: return
    # Get indices
    v = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(v) > 0: imin = v[0]
    v = CTK.varsFromWidget(VARS[1].get(), type=2)
    if len(v) > 0: imax = v[0]
    v = CTK.varsFromWidget(VARS[2].get(), type=2)
    if len(v) > 0: jmin = v[0]
    v = CTK.varsFromWidget(VARS[3].get(), type=2)
    if len(v) > 0: jmax = v[0]
    v = CTK.varsFromWidget(VARS[4].get(), type=2)
    if len(v) > 0: kmin = v[0]
    v = CTK.varsFromWidget(VARS[5].get(), type=2)
    if len(v) > 0: kmax = v[0]
    
    if CTK.__MAINTREE__ == 1:
        CTK.__MAINACTIVEZONES__ = CPlot.getActiveZones()
    active = []
    tp = Internal.appendBaseName2ZoneName(CTK.t, updateRef=False, 
                                          separator=Internal.SEP1)
    for z in CTK.__MAINACTIVEZONES__: active.append(tp[2][CTK.Nb[z]+1][2][CTK.Nz[z]])

    temp = C.newPyTree(['Base']); temp[2][1][2] += active
            
    CTK.dt = C.newPyTree(['Base'])
    for z in Internal.getZones(temp):
        try:
            zp = T.subzone(z, (imin,jmin,kmin), (imax,jmax,kmax)); zp[0] = z[0]
            CTK.dt[2][1][2].append(zp)
        except: pass

    CTK.display(CTK.dt, mainTree=CTK.SLICE)
    if CTK.TKPLOTXY is not None: CTK.TKPLOTXY.updateApp()
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkIJK', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Display mesh informations.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=0)
    Frame.columnconfigure(2, weight=0)
    
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkIJK')
    WIDGETS['frameMenu'] = FrameMenu

    #- VARS -
    # -0- Imin -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -1- Imax -
    V = TK.StringVar(win); V.set('-1'); VARS.append(V)
    # -2- Jmin -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -3- Jmax -
    V = TK.StringVar(win); V.set('-1'); VARS.append(V)
    # -4- Kmin -
    V = TK.StringVar(win); V.set('1'); VARS.append(V)
    # -5- Kmax -
    V = TK.StringVar(win); V.set('-1'); VARS.append(V)

    # - Subzone indices -
    F = TTK.LabelFrame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE,
                        text='I', font=CTK.FRAMEFONT, takefocus=1)
    F.grid(row=0, column=0)
    
    B = TTK.Button(F, text='+', width=5, command=setPlusI)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text='-', width=5, command=setMoinsI)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(F, text='X', width=1, command=setFullI)
    B.grid(row=0, column=2, sticky=TK.EW)
    
    B = TTK.Entry(F, textvariable=VARS[0], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[1], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    
    F = TTK.LabelFrame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE,
                       text='J', font=CTK.FRAMEFONT, takefocus=1)
    F.grid(row=0, column=1)
    
    B = TTK.Button(F, text='+', width=5, command=setPlusJ)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text='-', width=5, command=setMoinsJ)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(F, text='X', width=1, command=setFullJ)
    B.grid(row=0, column=2, sticky=TK.EW)
    
    B = TTK.Entry(F, textvariable=VARS[2], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[3], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)
    
    F = TTK.LabelFrame(Frame, borderwidth=2, relief=CTK.FRAMESTYLE,
                        text='K', font=CTK.FRAMEFONT, takefocus=1)
    F.grid(row=0, column=2)
    
    B = TTK.Button(F, text='+', width=5, command=setPlusK)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Button(F, text='-', width=5, command=setMoinsK)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Button(F, text='X', width=1, command=setFullK)
    B.grid(row=0, column=2, sticky=TK.EW)
    
    B = TTK.Entry(F, textvariable=VARS[4], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(F, textvariable=VARS[5], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=1, column=1, columnspan=2, sticky=TK.EW)

    #B = TTK.Button(Frame, text="Show", command=show)
    #B.grid(row=1, column=0, columnspan=6, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Show subzones.')

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
    (win, menu, file, tools) = CTK.minimal('tkIJJK '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
