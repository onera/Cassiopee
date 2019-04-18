# - mesure d'un modele -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import math, time

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def measure():
    if CTK.t == []: return
    prev = []
    w = WIDGETS['button']
    if CTK.__BUSY__ == False:
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            CPlot.unselectAllZones()
            l = []
            while (l == []):
                l = CPlot.getActivePoint()
                time.sleep(CPlot.__timeStep__)
                w.update()
                if CTK.__BUSY__ == False: break
            if CTK.__BUSY__ == True:
                if prev == []:
                    prev = l
                    CTK.TXT.insert('START', 'Click second point...\n')
                elif prev != l:
                    dist = (l[0]-prev[0])*(l[0]-prev[0])+\
                    (l[1]-prev[1])*(l[1]-prev[1])+\
                    (l[2]-prev[2])*(l[2]-prev[2])
                    dist = math.sqrt(dist)
                    CTK.TXT.insert('START', 'dist= '+str(dist)+'\n')
                    time.sleep(CPlot.__timeStep__)
                    prev = []
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
        CPlot.setState(cursor=0)
    else:
       CTK.__BUSY__ = False
       TTK.raiseButton(w)
       CPlot.setState(cursor=0)
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                          text='tkRuler', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Take measures by clicking.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkRuler')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Zone filter regexp -
    #V = TK.StringVar(win); V.set(''); VARS.append(V)
    
    # - Buttons -
    B = TTK.Button(Frame, text="Measure mode", command=measure)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Click on two points to obtain the distance.')
    WIDGETS['button'] = B
    
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
    (win, menu, file, tools) = CTK.minimal('tkRuler '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
