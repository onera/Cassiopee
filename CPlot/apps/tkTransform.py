# - simple transformations of mesh -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import KCore.Vector as Vector
import time

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def symetrize():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    axis = VARS[4].get()
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
    sel = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        sel.append(z)
    bb = G.bbox(sel)

    xmin = bb[0]; ymin = bb[1]; zmin = bb[2]
    xmax = bb[3]; ymax = bb[4]; zmax = bb[5]
    if axis == 'around XY-':
        X = ((xmin+xmax)*0.5, (ymin+ymax)*0.5, zmin)
        axe1 = (1.,0.,0.); axe2 = (0.,1.,0.)
    elif axis == 'XY+':
        X = ((xmin+xmax)*0.5, (ymin+ymax)*0.5, zmax)
        axe1 = (1.,0.,0.); axe2 = (0.,1.,0.)
    elif axis == 'around XZ-':
        X = ((xmin+xmax)*0.5, ymin, (zmin+zmax)*0.5)
        axe1 = (1.,0.,0.); axe2 = (0.,0.,1.)
    elif axis == 'around XZ+':
        X = ((xmin+xmax)*0.5, ymax, (zmin+zmax)*0.5)
        axe1 = (1.,0.,0.); axe2 = (0.,0.,1.)
    elif axis == 'around YZ-':
        X = (xmin, (ymin+ymax)*0.5, (zmin+zmax)*0.5)
        axe1 = (0.,1.,0.); axe2 = (0.,0.,1.)
    elif axis == 'around YZ+':
        X = (xmax, (ymin+ymax)*0.5, (zmin+zmax)*0.5)
        axe1 = (0.,1.,0.); axe2 = (0.,0.,1.)
    elif axis == 'around view':
        X = CPlot.getState('posEye')
        Y = CPlot.getState('posCam')
        axe1 = (X[0]-Y[0], X[1]-Y[1], X[2]-Y[2])
        axe2 = CPlot.getState('dirCam')
    else: X=(0.,0.,0.); axe1 = (1.,0.,0.); axe2 = (0.,1.,0.)
    
    CTK.saveTree()
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = T.symetrize(CTK.t[2][nob][2][noz], (X[0],X[1],X[2]), axe1, axe2)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TXT.insert('START', 'Zones have been symmetrized.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def rotate(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    axis = VARS[2].get()
    angle = CTK.varsFromWidget(VARS[3].get(), type=1)
    if len(angle) == 1: angle = angle[0]; X=None
    elif len(angle) == 4: X=(angle[1],angle[2],angle[3]); angle = angle[0] 
    else: 
        CTK.TXT.insert('START', 'Invalid angle or angle+rotation center.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if axis == 'around X': axe = (1.,0.,0.)
    elif axis == 'around Y': axe = (0.,1.,0.)
    elif axis == 'around Z': axe = (0.,0.,1.)
    elif axis == 'around view':
        pos = CPlot.getState('posCam')
        eye = CPlot.getState('posEye')
        axe = (eye[0]-pos[0], eye[1]-pos[1], eye[2]-pos[2])
    else: axe = (0.,0.,1.)
    try: angle = float(angle)
    except: angle = 0.
    
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
        
    CTK.saveTree()
    if X is None:
        sel = []
        for nz in nzs:
            nob = CTK.Nb[nz]+1
            noz = CTK.Nz[nz]
            z = CTK.t[2][nob][2][noz]
            sel.append(z)
        X = G.barycenter(sel)

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = T.rotate(CTK.t[2][nob][2][noz], (X[0],X[1],X[2]), axe, angle)
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TXT.insert('START', 'Zones have been rotated.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()
    
#==============================================================================
def translate():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    v = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(v) != 3:
        CTK.TXT.insert('START', 'Translation vector is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    axis = VARS[6].get()
    if axis == 'along view':
        posCam = CPlot.getState('posCam')
        posEye = CPlot.getState('posEye')
        dirCam = CPlot.getState('dirCam')
        axe1 = (posEye[0]-posCam[0], posEye[1]-posCam[1], posEye[2]-posCam[2])
        axe2 = dirCam
        axe3 = (axe1[1]*axe2[2]-axe1[2]*axe2[1],
                axe1[2]*axe2[0]-axe1[0]*axe2[2],
                axe1[0]*axe2[1]-axe1[1]*axe2[0])
        axe1 = Vector.normalize(axe1)
        axe2 = Vector.normalize(axe2)
        axe3 = Vector.normalize(axe3)
        ax = v[0]*axe1[0]+v[1]*axe2[0]+v[2]*axe3[0]
        ay = v[0]*axe1[1]+v[1]*axe2[1]+v[2]*axe3[1]
        az = v[0]*axe1[2]+v[1]*axe2[2]+v[2]*axe3[2]
        v[0] = ax; v[1] = ay; v[2] = az

    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        a = T.translate(CTK.t[2][nob][2][noz], (v[0], v[1], v[2]))
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TXT.insert('START', 'Zones have been translated.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def translateClick():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    prev = []
    w = WIDGETS['translate']
    if CTK.__BUSY__ == False:
        CTK.__BUSY__ = True
        nzs = CPlot.getSelectedZones()
        CTK.TXT.insert('START', 'Click start point...\n')
        TTK.sunkButton(w)
        while CTK.__BUSY__:
            CPlot.unselectAllZones()
            l = []
            while l == []:
                l = CPlot.getActivePoint()
                time.sleep(CPlot.__timeStep__)
                w.update()
                if CTK.__BUSY__ == False: break
            if CTK.__BUSY__:
                if prev == []:
                    prev = l
                    if nzs == []: nzs = CPlot.getSelectedZones()
                    CTK.TXT.insert('START', 'Click end point...\n')
                elif prev != l:
                    CTK.saveTree()
                    vx = l[0]-prev[0]
                    vy = l[1]-prev[1]
                    vz = l[2]-prev[2]
                    for nz in nzs:
                        nob = CTK.Nb[nz]+1
                        noz = CTK.Nz[nz]
                        z = CTK.t[2][nob][2][noz]
                        a = T.translate(z, (vx, vy, vz))
                        CTK.replace(CTK.t, nob, noz, a)
                    CTK.TKTREE.updateApp()
                    CTK.TXT.insert('START', 'Zones translated.\n')
                    CPlot.render()
                    prev = []
                    break
        CTK.__BUSY__ = False
        TTK.raiseButton(w)
    else:
       CTK.__BUSY__ = False
       TTK.raiseButton(w)

#==============================================================================
def scale():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    v = CTK.varsFromWidget(VARS[1].get(), type=1)
    if len(v) != 1 and len(v) != 3 and len(v) != 6:
        CTK.TXT.insert('START', 'Scale factor is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    axis = VARS[5].get()
    if axis == 'along XYZ':
        axe1 = (1,0,0); axe2 = (0,1,0); axe3 = (0,0,1)
    else: # view
         posCam = CPlot.getState('posCam')
         posEye = CPlot.getState('posEye')
         dirCam = CPlot.getState('dirCam')
         axe1 = (posEye[0]-posCam[0], posEye[1]-posCam[1], posEye[2]-posCam[2])
         axe2 = dirCam
         axe3 = (axe1[1]*axe2[2]-axe1[2]*axe2[1],
                 axe1[2]*axe2[0]-axe1[0]*axe2[2],
                 axe1[0]*axe2[1]-axe1[1]*axe2[0])

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()

    selection = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        selection.append(z)
    if len(v) == 6: # center is given
        X = [v[3],v[4],v[5]]
    else:
        X = G.barycenter(selection)
    if len(v) == 1 and v[0] == 0.: # scale unitaire 
        bbox = G.bbox(selection)
        dx = bbox[3]-bbox[0]
        dy = bbox[4]-bbox[1]
        dz = bbox[5]-bbox[2]
        if dx >= dy and dx >= dz: v[0] = 1./dx
        if dy >= dx and dy >= dz: v[0] = 1./dy
        if dz >= dy and dz >= dx: v[0] = 1./dz
    
    list = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        list.append( (nob,noz,nz) )
        z = CTK.t[2][nob][2][noz]
        if len(v) == 1:
            a = T.homothety(z, (X[0],X[1],X[2]), v[0])
        else:
            z = T.contract(z, (X[0],X[1],X[2]), axe2, axe3, v[0])
            z = T.contract(z, (X[0],X[1],X[2]), axe1, axe3, v[1])
            a = T.contract(z, (X[0],X[1],X[2]), axe1, axe2, v[2])
        CTK.replace(CTK.t, nob, noz, a)
    CTK.TXT.insert('START', 'Zones have been scaled.\n')
    CTK.TKTREE.updateApp()
    CPlot.render()

#=========================================================================
def changeFrame():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    mode = VARS[7].get()
    assert(mode in dir(T))

    args = CTK.varsFromWidget(VARS[8].get(), type=1)
    if len(args) != 6:
        CTK.TXT.insert('START', '{} requires 6 values.\nOrigin: x0;y0;z0\n Axis tx;ty;tz.\n'.format(mode))
        CTK.TXT.insert('START', 'Error: ', 'Error')
        return

    origin = (args[0], args[1], args[2])
    axis   = (args[3], args[4], args[5])
    CTK.saveTree()

    fail = False;
    errors = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        try:
            func = getattr(T, mode)
            a = func(CTK.t[2][nob][2][noz], origin, axis)
            CTK.replace(CTK.t, nob, noz, a)
        except Exception as e:
            fail = True; errors += [0,str(e)]

    if not fail: CTK.TXT.insert('START', '{} done.\n'.format(mode))
    else:
        Panels.displayErrors(errors, header='Error: {}'.format(mode))
        CTK.TXT.insert('START', '{} fails for at least one zone.\n'.format(mode))
        CTK.TXT.insert('START', 'Warning: ', 'Warning')
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkTransform', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='General block transformations.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=4)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkTransform')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Translation vector -
    V = TK.StringVar(win); V.set('0;0;0'); VARS.append(V)
    # -1- Scale factors -
    V = TK.StringVar(win); V.set('1;1;1'); VARS.append(V)
    # -2- Rotate axis
    V = TK.StringVar(win); V.set('around X'); VARS.append(V)
    # -3- Rotation angle
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    # -4- Symetry axis
    V = TK.StringVar(win); V.set('around XY-'); VARS.append(V)
    # -5- Scale axis
    V = TK.StringVar(win); V.set('along XYZ'); VARS.append(V)
    # -6- Translation axis
    V = TK.StringVar(win); V.set('along XYZ'); VARS.append(V)
    # -7- cart2Cyl or cyl2Cart
    V = TK.StringVar(win); V.set('cart2Cyl'); VARS.append(V)
    # -8- Origin and axis for cart2Cyl or cyl2Cart
    V = TK.StringVar(win); V.set('(0.,0.,0.); (1.,0.,0.)'); VARS.append(V)

    # - Translate -
    B = TTK.Button(Frame, text="Translate", command=translate)
    B.grid(row=0, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Translate zone(s) of pyTree.')
    B = TTK.OptionMenu(Frame, VARS[6], 'along XYZ', 'along view')
    B.grid(row=0, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=5)
    B.grid(row=0, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Translation vector.')

    # - Translate from here to here -
    B = TTK.Button(Frame, text="Translate by clicking", command=translateClick)
    B.grid(row=1, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Translate zone(s) of pyTree \nby clicking.')
    WIDGETS['translate'] = B

    # - Scale -
    B = TTK.Button(Frame, text="Scale", command=scale)
    B.grid(row=2, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Scale zone(s) of pyTree.')
    B = TTK.OptionMenu(Frame, VARS[5], 'along XYZ', 'along view')
    B.grid(row=2, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    B.grid(row=2, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Scale factors (+ optionally center for scaling): uniform f value or fx;fy;fz (+ optionally cx;cy;cz)\nFor automatic adimensioning, use 0.')

    # - Rotate -
    B = TTK.Button(Frame, text="Rotate", command=rotate)
    B.grid(row=3, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Rotate zone(s) of pyTree.')
    B = TTK.OptionMenu(Frame, VARS[2], 'around X', 'around Y',
                       'around Z', 'around view')
    B.grid(row=3, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=5)
    B.grid(row=3, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='angle (degrees) or \nangle; Xc;Yc;Zc (angle+rotation center)\nIf center is not specified, rotate around barycenter of zones.')

    # - Symetrize -
    B = TTK.Button(Frame, text="Mirror", command=symetrize)
    B.grid(row=4, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Symetrize zone(s) of pyTree.')
    B = TTK.OptionMenu(Frame, VARS[4], 'around XY-', 'around XY+', 'around YZ-',
                       'around YZ+', 'around XZ-', 'around XZ+', 'around view')
    B.grid(row=4, column=1, columnspan=2, sticky=TK.EW)
    
    # - cart2Cyl and cyl2Cart -
    B = TTK.Button(Frame, text="Apply", command=changeFrame)
    B.grid(row=5, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Convert grid coordinates from cartesian->cylindrical or cylindrical->cartesian frame.\nTree is modified.')
    B = TTK.OptionMenu(Frame, VARS[7], 'cart2Cyl', 'cyl2Cart')
    B.grid(row=5, column=1, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[8], background='White')
    B.grid(row=5, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Origin: x0;y0;z0\nAxis: tx;ty;tz.')

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
    if len(sys.argv) == 2:
        CTK.FILE = sys.argv[1]
        try:
            CTK.t = C.convertFile2PyTree(CTK.FILE)
            (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
            CTK.display(CTK.t)
        except: pass

    # Main window
    (win, menu, file, tools) = CTK.minimal('tkTransform '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
