# - mapuv app -
try: import tkinter as TK
except: import Tkinter as TK
import Geom.PyTree as D
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setNormalDeviation(event=None):
    VARS[1].set('Normal deviation [%.2f]'%(WIDGETS['deviation'].get() / 10.))

#==============================================================================
def generateUVMap():
    if CTK.t == []: return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.setCursor(2, WIDGETS['generate'])
    CTK.saveTree()

    ndeviation = WIDGETS['deviation'].get()/10.
    print("ndeviation=", ndeviation)

    for nz in nzs:
        nob = CTK.Nb[nz]+1; noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        (z, color, normal) = D.getUV(z, ndeviation, 0.)
        CPlot.replace(CTK.t, nob, noz, z)
        C.convertPyTree2File(color, 'color#%s.png'%z[0])
        C.convertPyTree2File(normal, 'bump#%s.png'%z[0])
    CTK.TXT.insert('START', 'UV map generated.\n')
    CPlot.render()
    CTK.setCursor(0, WIDGETS['generate'])

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkMapUV  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Map uv.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkMapUV')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- setting (not used for now) -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -1- normal deviation info bulle
    V = TK.StringVar(win); V.set('Normal deviation.'); VARS.append(V)

    # - set Normal deviation -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL,
                  command=setNormalDeviation, showvalue=0, borderwidth=1, value=100)
    WIDGETS['deviation'] = B
    WIDGETS['deviation'].set(20.)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[1])

    # - Generate UV map -
    B = TTK.Button(Frame, text='Generate UV map', command=generateUVMap)
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    WIDGETS['generate'] = B
    BB = CTK.infoBulle(parent=B, text='Generate UV map for selection.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkMapUV')
    except: pass
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
    (win, menu, file, tools) = CTK.minimal('tkPaint '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
