# - MMGs -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import Generator.PyTree as G
import Converter.Internal as Internal

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def remesh():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    ridgeAngle = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(ridgeAngle) != 1:
            CTK.TXT.insert('START', 'Invalid ridge angle.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    ridgeAngle = ridgeAngle[0]
    hausd = CTK.varsFromWidget(VARS[1].get(), type=1)
    if len(hausd) != 1:
            CTK.TXT.insert('START', 'Invalid hausdorff distance.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    hausd = hausd[0]
    hmin = CTK.varsFromWidget(VARS[2].get(), type=1)
    if len(hmin) != 1:
            CTK.TXT.insert('START', 'Invalid hmin.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    hmin = hmin[0]
    hmax = CTK.varsFromWidget(VARS[3].get(), type=1)
    if len(hmax) != 1:
            CTK.TXT.insert('START', 'Invalid hmax.\n')
            CTK.TXT.insert('START', 'Error: ', 'Error')
    hmax = hmax[0]
    
    CTK.saveTree()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        zp = G.mmgs(z, ridgeAngle=ridgeAngle, hausd=hausd, 
                    hmin=hmin, hmax=hmax)
        CTK.replace(CTK.t, nob, noz, zp)
        CTK.TXT.insert('START', 'Surface remeshed.\n')
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkMMGs', font=CTK.FRAMEFONT, takefocus=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkMMGs')
    WIDGETS['frameMenu'] = FrameMenu
    
    # - VARS -
    # -0- ridgeAngle -
    V = TK.StringVar(win); V.set('45.'); VARS.append(V)
    if 'tkMMGsRidgeAngle' in CTK.PREFS: V.set(CTK.PREFS['tkMMGsRidgeAngle'])
    # -1- hausd -
    V = TK.StringVar(win); V.set('0.01'); VARS.append(V)
    if 'tkMMGsHausd' in CTK.PREFS: V.set(CTK.PREFS['tkMMGsHausd'])
    # -2- hmin -
    V = TK.StringVar(win); V.set('0.0'); VARS.append(V)
    # -3- hmax -
    V = TK.StringVar(win); V.set('0.0'); VARS.append(V)

    # RidgeAngle
    B = TTK.Label(Frame, text="ridgeAngle")
    BB = CTK.infoBulle(parent=B, text='Ridge angle (degree).')
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=10)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Ridge angle (degree).')

    # hausd
    B = TTK.Label(Frame, text="hausd")
    BB = CTK.infoBulle(parent=B, text='Chordal error (hausdorff distance).')
    B.grid(row=1, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=10)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Chordal error (hausdorff distance).')

    # hmin
    B = TTK.Label(Frame, text="hmin")
    BB = CTK.infoBulle(parent=B, text='Minimum step in remeshed surface.')
    B.grid(row=2, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=10)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Minimum step in remeshed surface.')

    # hmax
    B = TTK.Label(Frame, text="hmax")
    BB = CTK.infoBulle(parent=B, text='Maximum step in remeshed surface.')
    B.grid(row=3, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=10)
    B.grid(row=3, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Maximum step in remeshed surface.')

    # - Remesh -
    B = TTK.Button(Frame, text="Remesh", command=remesh)
    B.grid(row=4, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Remesh a TRI surface.')
    
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
    CTK.PREFS['tkTetraMesherType'] = VARS[0].get()
    CTK.savePrefFile()
    
#==============================================================================
def resetApp():
    VARS[0].set('netgen')
    CTK.PREFS['tkTetraMesherType'] = VARS[0].get()
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
    (win, menu, file, tools) = CTK.minimal('tkTetraMesher '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
