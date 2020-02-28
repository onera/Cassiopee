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
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    Frame.columnconfigure(3, weight=1)
    Frame.columnconfigure(4, weight=1)
    Frame.columnconfigure(5, weight=1)
    
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
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=0, column=2, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=0, column=3, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[4], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=0, column=4, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White', width=5)
    B.bind('<Return>', show)
    B.grid(row=0, column=5, sticky=TK.EW)

    B = TTK.Button(Frame, text="Show", command=show)
    B.grid(row=1, column=0, columnspan=6, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Show subzones.')

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
