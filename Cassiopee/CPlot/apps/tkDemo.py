# - tkDemo -
"""Demo like app."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Generator.PyTree as G
import Converter.Internal as Internal
import Transform.PyTree as T
import math, time
import CPlot.iconics as iconics

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setPath():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    selected = ''
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        selected += CTK.t[2][nob][0]+'/'+z[0]+';'
    selected = selected[0:-1]
    VARS[0].set(selected)

#==============================================================================
def rotate():
    if CTK.t == []: return
    if not CTK.__BUSY__:
        bb = G.bbox(CTK.t)
        xc = 0.5*(bb[3]+bb[0])
        yc = 0.5*(bb[4]+bb[1])
        zc = 0.5*(bb[5]+bb[2])
        pos = CPlot.getState('posCam')
        posCam = [pos[0], pos[1], pos[2]]
        #posEye = [xc, yc, zc]
        #dirCam = [0,0,1]
        CTK.__BUSY__ = True
        TTK.sunkButton(WIDGETS['rotate'])
        CPlot.setState(cursor=2)
        CTK.setCursor(2, WIDGETS['rotate'])
        i = 0
        while CTK.__BUSY__:
            speed = WIDGETS['speed'].get() * 0.0006 / 100.
            cs = math.cos(speed*i * math.pi/180)
            ss = math.sin(speed*i * math.pi/180)
            px = cs * (posCam[0]-xc) + ss * (posCam[1]-yc) + xc
            py = -ss * (posCam[0]-xc) + cs * (posCam[1]-yc) + yc
            posCam[0] = px; posCam[1] = py
            CPlot.setState(posCam=posCam)
            time.sleep(CPlot.__timeStep__)
            WIDGETS['rotate'].update()
            i += 1
        CTK.__BUSY__ = False
        TTK.raiseButton(WIDGETS['rotate'])
        CPlot.setState(cursor=0)
        CTK.setCursor(0, WIDGETS['rotate'])
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(WIDGETS['rotate'])
        CPlot.setState(cursor=0)
        CTK.setCursor(0, WIDGETS['rotate'])

#==============================================================================
def fly():
    if CTK.t == []: return
    if not CTK.__BUSY__:
        bb = G.bbox(CTK.t)
        xc = 0.5*(bb[3]+bb[0])
        yc = 0.5*(bb[4]+bb[1])
        zc = 0.5*(bb[5]+bb[2])
        pos = CPlot.getState('posCam')
        posCam = [pos[0], pos[1], pos[2]]
        step = 0.0001
        sigma = 10; beta = 8/3; ro = 28
        CTK.__BUSY__ = True
        TTK.sunkButton(WIDGETS['fly'])
        CPlot.setState(cursor=2)
        CTK.setCursor(2, WIDGETS['fly'])
        i = 0
        while CTK.__BUSY__:
            speed = WIDGETS['speed'].get() / 50.
            CPlot.setState(posCam=posCam)
            time.sleep(CPlot.__timeStep__)
            x = posCam[0]-xc; y = posCam[2]-yc; z = posCam[1]-zc
            x = 10*x /max((bb[3]-bb[0]), 1.e-10)
            y = 10*y /max((bb[4]-bb[1]), 1.e-10)
            z = 30*z /max((bb[5]-bb[2]), 1.e-10)
            xp = x + step*speed*sigma*(y - x)
            yp = y + step*speed*(ro*x - y - x*z)
            zp = z + step*speed*(x*y - beta*z)
            xp = xp * (bb[3]-bb[0])*0.1 + xc
            yp = yp * (bb[4]-bb[1])*0.1 + yc
            zp = zp * (bb[5]-bb[2])/30. + zc
            posCam[0] = xp; posCam[1] = zp; posCam[2] = yp
            WIDGETS['fly'].update()
            i += 1
        CTK.__BUSY__ = False
        TTK.raiseButton(WIDGETS['fly'])
        CPlot.setState(cursor=0)
        CTK.setCursor(0, WIDGETS['fly'])
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(WIDGETS['fly'])
        CPlot.setState(cursor=0)
        CTK.setCursor(0, WIDGETS['fly'])

#==============================================================================
def orbite():
    if CTK.t == []: return
    if not CTK.__BUSY__:
        name = VARS[0].get()
        names = name.split(';')
        # Get paths
        paths = []
        for v in names:
            v = v.lstrip(); v = v.rstrip()
            sname = v.split('/', 1)
            bases = Internal.getNodesFromName1(CTK.t, sname[0])
            if bases != []:
                nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
                for z in nodes:
                    if z[0] == sname[1]: paths.append(z)
        # Keep only 1D arrays
        path = []
        for p in paths:
            dim = Internal.getZoneDim(p)
            if dim[0] == 'Unstructured' and dim[3] == 'BAR': path.append(p)
            if dim[0] == 'Structured' and dim[2] == 1 and dim[3] == 1:
                path.append(C.convertArray2Tetra(p))
        if path == []: return
        path = T.join(path)
        path = G.close(path)
        path = C.convertBAR2Struct(path)
        dim = Internal.getZoneDim(path)
        N = dim[1]

        CTK.__BUSY__ = True
        TTK.sunkButton(WIDGETS['orbite'])
        CPlot.setState(cursor=2)
        CTK.setCursor(2, WIDGETS['orbite'])
        i = 0
        while CTK.__BUSY__:
            speed = 100. - WIDGETS['speed'].get()
            time.sleep(CPlot.__timeStep__*speed*0.06)
            if i > N-1: i = 0
            if i+N/10 > N-1: inc = 1
            else: inc = N/10
            posCam = C.getValue(path, Internal.__GridCoordinates__, i)
            posEye = C.getValue(path, Internal.__GridCoordinates__, i+inc)
            CPlot.setState(posCam=posCam, posEye=posEye)
            WIDGETS['orbite'].update()
            i += 1
        CTK.__BUSY__ = False
        TTK.raiseButton(WIDGETS['orbite'])
        CPlot.setState(cursor=0)
        CTK.setCursor(0, WIDGETS['orbite'])
    else:
        CTK.__BUSY__ = False
        TTK.raiseButton(WIDGETS['orbite'])
        CPlot.setState(cursor=0)
        CTK.setCursor(0, WIDGETS['orbite'])

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkDemo  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Demo mode.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=2)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkDemo')
    WIDGETS['frameMenu'] = FrameMenu

    # -0- Orbiter path -
    V = TK.StringVar(win); V.set(''); VARS.append(V)

    # - Slider -
    B = TTK.Scale(Frame, from_=0, to=100, orient=TK.HORIZONTAL, showvalue=0,
                  borderwidth=1, value=50)
    WIDGETS['speed'] = B
    B.grid(row=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Speed.')

    # - Rotate -
    B = TTK.Button(Frame, text="Rotate", command=rotate)
    B.grid(row=1, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Rotate your model automatically.')
    WIDGETS['rotate'] = B

    # - Fly -
    B = TTK.Button(Frame, text="Fly", command=fly)
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Fly around your model.')
    WIDGETS['fly'] = B

    # - Orbite -
    B = TTK.Button(Frame, text="Path", command=setPath,
                   image=iconics.PHOTO[8], compound=TK.RIGHT, padx=0, pady=0)
    B.grid(row=2, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Set orbiter path.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Orbiter path (curve).')

    B = TTK.Button(Frame, text="Orbite", command=orbite)
    B.grid(row=3, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Orbite following a path.')
    WIDGETS['orbite'] = B

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['RenderNoteBook'].add(WIDGETS['frame'], text='tkDemo')
    except: pass
    CTK.WIDGETS['RenderNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['RenderNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkDemo '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
