# - sculpt app -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Transform.PyTree as T
import Generator.PyTree as G
import Geom.PyTree as D
import Intersector.PyTree as XOR
import time

# local widgets list
WIDGETS = {}; VARS = []
TOOLS = []

#==============================================================================
# Cree les outils, les mets dans TOOLS
#==============================================================================
def createTools():
    global TOOLS
    if TOOLS != []: return # deja crees
    # square
    P0 = (0,0,0)
    # pointe
    P0 = (0,0,0)
    P1 = (-1,-1,1)
    P2 = (-1,1,1)
    P3 = (-1,1,-1)
    P4 = (-1,-1,-1)
    t1 = D.triangle(P0,P1,P2)
    t2 = D.triangle(P0,P2,P3)
    t3 = D.triangle(P0,P3,P4)
    t4 = D.triangle(P0,P4,P1)
    t = T.join([t1,t2,t3,t4])
    t = G.close(t)
    TOOLS.append(t)
    # biseau horizontal
    P0 = (0,-1,0)
    P1 = (0,1,0)
    P1 = (-1,-1,1)
    P2 = (-1,1,1)
    P3 = (-1,1,-1)
    P4 = (-1,-1,-1)
    t1 = D.triangle(P0,P1,P2)
    t2 = D.triangle(P0,P2,P3)
    t3 = D.triangle(P0,P3,P4)
    t4 = D.triangle(P0,P4,P1)
    t = T.join([t1,t2,t3,t4])
    t = G.close(t)
    # biseau vertical
    # pointe spherique
    return TOOLS

#==============================================================================
def sculpt():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    TOOLS = createTools()
    #CTK.display(TOOLS)

    bbox = G.bbox(CTK.t)
    size = max(bbox[3]-bbox[0], bbox[4]-bbox[1], bbox[5]-bbox[2])
    CPlot.unselectAllZones()

    w = WIDGETS['sculpt']
    if not CTK.__BUSY__:
        CTK.__BUSY__ = True
        TTK.sunkButton(w)
        CPlot.setState(cursor=1)
        while CTK.__BUSY__:
            l = []
            while l == []:
                nz = CPlot.getSelectedZone()
                l = CPlot.getActivePointIndex()
                time.sleep(CPlot.__timeStep__)
                w.update()
                if not CTK.__BUSY__: break
            if CTK.__BUSY__:
                nob = CTK.Nb[nz]+1
                noz = CTK.Nz[nz]
                CTK.saveTree()
                depth = 0.5*WIDGETS['depth'].get()/100.
                depth = size*depth
                width = 0.5*WIDGETS['width'].get()/100.
                width = size*width
                brushType = VARS[2].get()
                z = CTK.t[2][nob][2][noz]
                posCam = CPlot.getState('posCam')
                posEye = CPlot.getState('posEye')
                vect = (posEye[0]-posCam[0], posEye[1]-posCam[1],
                        posEye[2]-posCam[2])
                if brushType == 'Deform':
                    click = CPlot.getActivePoint()
                    point = (click[0], click[1], click[2])
                    z = T.deformPoint(z, point, vect, depth, width)
                    CTK.replace(CTK.t, nob, noz, z)
                elif brushType == 'Sphere':
                    click = CPlot.getActivePoint()
                    center = (click[0], click[1], click[2])
                    s = D.sphere(center, depth, N=10)
                    s = C.convertArray2Tetra(s)
                    s = G.close(s)
                    z = C.convertArray2Tetra(z)
                    z = G.close(z)
                    z = XOR.booleanMinus(z, s)
                    CTK.replace(CTK.t, nob, noz, z)
                elif brushType == 'Cube':
                    click = CPlot.getActivePoint()
                    center = (click[0], click[1], click[2])
                    s = D.sphere(center, depth, N=20)
                    s = C.convertArray2Tetra(s)
                    s = G.close(s)
                    z = C.convertArray2Tetra(z)
                    z = G.close(z)
                    z = XOR.booleanMinus(z, s)
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
def setDepth(event=None):
    depth = 0.5*WIDGETS['depth'].get()/100.
    VARS[4].set('Depth [%.2f %%]'%depth)

def setWidth(event=None):
    width = 0.5*WIDGETS['width'].get()/100.
    VARS[5].set('Width [%.2f %%]'%width)

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkSculpt  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Sculpt surfaces by deformations.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkSculpt')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- depth -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -1- width -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    # -2- Brush
    V = TK.StringVar(win); V.set('Deform'); VARS.append(V)
    # -3- Negative/positive depth
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -4- Depth info bulle
    V = TK.StringVar(win); V.set('Depth.'); VARS.append(V)
    # -5- Width info bulle
    V = TK.StringVar(win); V.set('Width.'); VARS.append(V)

    # - Depth -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  command=setDepth, borderwidth=1, value=50)
    WIDGETS['depth'] = B
    B.grid(row=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[4])
    #B = TK.Checkbutton(Frame, text='', variable=VARS[3])
    #BB = CTK.infoBulle(parent=B, text='Check for additive depth.')
    #B.grid(row=0, column=1, sticky=TK.EW)

    # - Width -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  command=setWidth, borderwidth=1, value=50)
    WIDGETS['width'] = B
    B.grid(row=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, textVariable=VARS[5])

    # - Brush -
    B = TTK.OptionMenu(Frame, VARS[2], 'Deform', 'Sphere')
    B.grid(row=2, column=0, sticky=TK.EW)

    # - Sculpt mode -
    B = TTK.Button(Frame, text="Sculpt mode", command=sculpt)
    BB = CTK.infoBulle(parent=B, text='Enter sculpt mode.')
    B.grid(row=2, column=1, columnspan=1, sticky=TK.EW)
    WIDGETS['sculpt'] = B

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkSculpt')
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
    (win, menu, file, tools) = CTK.minimal('tkSculpt '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
