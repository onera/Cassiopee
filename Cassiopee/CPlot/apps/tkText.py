# - tkText -
"""Create Texts."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Geom.PyTree as D
import Transform.PyTree as T
import Generator.PyTree as G
import KCore.Vector as Vector
import math

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def createText(event=None):
    CTK.saveTree()
    text = VARS[0].get()
    type = VARS[1].get()
    font = VARS[3].get()

    smoothness = VARS[2].get()
    smooth = 0
    if smoothness == 'Regular': smooth = 0
    elif smoothness == 'Smooth': smooth = 2
    elif smoothness == 'Very smooth': smooth = 4
    CTK.t = C.addBase2PyTree(CTK.t, 'TEXT', 3)
    nodes = Internal.getNodesFromName1(CTK.t, 'TEXT')
    base = nodes[0]
    if type == '3D': a = D.text3D(text, smooth=smooth, font=font)
    elif type == '2D': a = D.text2D(text, smooth=smooth, font=font)
    elif type == '1D':
        a = D.text1D(text, smooth=smooth, font=font)
        a = C.convertArray2Tetra(a); a = T.join(a)

    # Modification de l'angle, de la position et de la taille du texte
    # en fonction du point de vue
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')
    BB = G.bbox(a)
    xc = 0.5*(BB[3]+BB[0]); yc = 0.5*(BB[4]+BB[1]); zc = 0.5*(BB[5]+BB[2])
    a = T.translate(a, (posEye[0]-xc, posEye[1]-yc, posEye[2]-zc) )
    lx = posEye[0]-posCam[0]
    ly = posEye[1]-posCam[1]
    lz = posEye[2]-posCam[2]
    if lx*lx + ly*ly + lz*lz == 0.: lx = -1
    if dirCam[0]*dirCam[0] + dirCam[1]*dirCam[1] + dirCam[2]*dirCam[2] == 0.:
        dirCam = (0,0,1)
    ll = math.sqrt(lx*lx + ly*ly + lz*lz)
    a = T.homothety(a, (posEye[0], posEye[1], posEye[2]), 0.01*ll)
    ux = dirCam[1]*lz - dirCam[2]*ly
    uy = dirCam[2]*lx - dirCam[0]*lz
    uz = dirCam[0]*ly - dirCam[1]*lx
    a = T.rotate(a, (posEye[0], posEye[1], posEye[2]),
                 ((1,0,0),(0,1,0),(0,0,1)),
                 ((-ux,-uy,-uz),dirCam,(lx,ly,lz)))

    nob = C.getNobOfBase(base, CTK.t)
    CTK.add(CTK.t, nob, -1, a)
    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Text created.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Replace text
#==============================================================================
def replaceText(event=None):
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    CTK.saveTree()
    # Recupere l'OBB de la selection
    Z = []; dels = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        Z.append(CTK.t[2][nob][2][noz])
        dels.append(CTK.t[2][nob][0]+Internal.SEP1+CTK.t[2][nob][2][noz][0])
    nob0 = CTK.Nb[nzs[0]]+1
    try: a = T.join(Z)
    except: a = Z[0]
    OBB = G.BB(a, method='OBB')

    P0 = C.getValue(OBB, 'GridCoordinates', 0)
    P1 = C.getValue(OBB, 'GridCoordinates', 1)
    P2 = C.getValue(OBB, 'GridCoordinates', 2)
    P3 = C.getValue(OBB, 'GridCoordinates', 4)
    v1 = Vector.sub(P1, P0)
    n1 = Vector.norm(v1)
    v2 = Vector.sub(P2, P0)
    n2 = Vector.norm(v2)
    v3 = Vector.sub(P3, P0)
    n3 = Vector.norm(v3)
    v1p = v1; v2p = v2; v3p = v3
    if abs(n1) < 1.e-12: v1 = Vector.cross(v2,v3); v1p = (0.,0.,0.)
    elif abs(n2) < 1.e-12: v2 = Vector.cross(v1,v3); v2p = (0.,0.,0.)
    elif abs(n3) < 1.e-12: v3 = Vector.cross(v2,v3); v3p = (0.,0.,0.)

    # Essaie de matcher les vecteur sur la vue p1,p2,p3
    # On suppose que dirCam doit etre e2, ...
    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')
    e2 = dirCam
    e3 = Vector.sub(posCam, posEye)
    e1 = Vector.cross(e2, e3)

    f1 = None; f2 = None; f3 = None; Pt = P0
    s1 = Vector.dot(e1, v1)
    s2 = Vector.dot(e1, v2)
    s3 = Vector.dot(e1, v3)
    if abs(s1) > abs(s2) and abs(s1) > abs(s3):
        if s1 > 0: f1 = v1
        else: f1 = Vector.mul(-1.,v1); Pt = Vector.add(Pt,v1p)
    elif abs(s2) > abs(s1) and abs(s2) > abs(s3):
        if s2 > 0: f1 = v2
        else: f1 = Vector.mul(-1.,v2); Pt = Vector.add(Pt,v2p)
    elif abs(s3) > abs(s1) and abs(s3) > abs(s2):
        if s3 > 0: f1 = v3
        else: f1 = Vector.mul(-1.,v3); Pt = Vector.add(Pt,v3p)
    s1 = Vector.dot(e2, v1)
    s2 = Vector.dot(e2, v2)
    s3 = Vector.dot(e2, v3)
    if abs(s1) > abs(s2) and abs(s1) > abs(s3):
        if s1 > 0: f2 = v1
        else: f2 = Vector.mul(-1.,v1); Pt = Vector.add(Pt,v1p)
    elif abs(s2) > abs(s1) and abs(s2) > abs(s3):
        if s2 > 0: f2 = v2
        else: f2 = Vector.mul(-1.,v2); Pt = Vector.add(Pt,v2p)
    elif abs(s3) > abs(s1) and abs(s3) > abs(s2):
        if s3 > 0: f2 = v3
        else: f2 = Vector.mul(-1.,v3); Pt = Vector.add(Pt,v3p)
    s1 = Vector.dot(e3, v1)
    s2 = Vector.dot(e3, v2)
    s3 = Vector.dot(e3, v3)
    if abs(s1) > abs(s2) and abs(s1) > abs(s3):
        if s1 > 0: f3 = v1
        else: f3 = Vector.mul(-1.,v1); Pt = Vector.add(Pt,v1p)
    elif abs(s2) > abs(s1) and abs(s2) > abs(s3):
        if s2 > 0: f3 = v2
        else: f3 = Vector.mul(-1.,v2); Pt = Vector.add(Pt,v2p)
    elif abs(s3) > abs(s1) and abs(s3) > abs(s2):
        if s3 > 0: f3 = v3
        else: f3 = Vector.mul(-1.,v3); Pt = Vector.add(Pt,v3p)
    (x0,y0,z0) = Pt
    n2 = Vector.norm(f2)

    # Cree le texte
    text = VARS[0].get()
    type = VARS[1].get()
    font = VARS[3].get()

    smoothness = VARS[2].get()
    smooth = 0
    if smoothness == 'Regular': smooth = 0
    elif smoothness == 'Smooth': smooth = 2
    elif smoothness == 'Very smooth': smooth = 4
    if type == '3D': a = D.text3D(text, smooth=smooth, font=font)
    elif type == '2D': a = D.text2D(text, smooth=smooth, font=font)
    elif type == '1D':
        a = D.text1D(text, smooth=smooth, font=font)
        a = C.convertArray2Tetra(a); a = T.join(a)
    BB = G.bbox(a)
    h2 = BB[4]-BB[1]

    # Scale, positionne le texte
    factor = n2/h2
    a = T.homothety(a, (BB[0],BB[1],BB[2]), factor)
    a = T.translate(a, (x0-BB[0],y0-BB[1],z0-BB[2]))
    a = T.rotate(a, (x0,y0,z0), ((1,0,0),(0,1,0), (0,0,1)), ((f1,f2,f3)))

    CTK.t = CPlot.deleteSelection(CTK.t, CTK.Nb, CTK.Nz, nzs)
    CPlot.delete(dels)
    CTK.add(CTK.t, nob0, -1, a)
    #C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Text replaced.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkText  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=1)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkText')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- corps du texte -
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- 1D/2D/3D -
    V = TK.StringVar(win); V.set('3D'); VARS.append(V)
    if 'tkTextDim' in CTK.PREFS: V.set(CTK.PREFS['tkTextDim'])
    # -2-  Smoothness -
    V = TK.StringVar(win); V.set('Regular'); VARS.append(V)
    if 'tkTextSmoothness' in CTK.PREFS:
        V.set(CTK.PREFS['tkTextSmoothness'])
    # -3- Font -
    V = TK.StringVar(win); V.set('vera'); VARS.append(V)
    if 'tkTextFont' in CTK.PREFS: V.set(CTK.PREFS['tkTextFont'])

    # - 1D/2D/3D -
    B = TTK.OptionMenu(Frame, VARS[1], '3D', '2D', '1D')
    BB = CTK.infoBulle(parent=B, text='Text type.')
    B.grid(row=0, column=0, sticky=TK.EW)

    # - Font -
    B = TTK.OptionMenu(Frame, VARS[3], 'vera', 'chancery',
                       'courier', 'text1', 'nimbus')
    BB = CTK.infoBulle(parent=B, text='Font type.')
    B.grid(row=0, column=1, sticky=TK.EW)

    # - Smoothness -
    B = TTK.OptionMenu(Frame, VARS[2], 'Regular', 'Smooth', 'Very smooth')
    BB = CTK.infoBulle(parent=B, text='Font smoothness.')
    B.grid(row=0, column=2, sticky=TK.EW)

    # - Text -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.bind('<Return>', createText)
    B.grid(row=1, column=0, columnspan=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Your text.')

    # - Create text -
    B = TTK.Button(Frame, text="Create text", command=createText)
    B.grid(row=2, column=0, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create the text in TEXT basis.')

    # - Replace text -
    B = TTK.Button(Frame, text="Replace", command=replaceText)
    B.grid(row=2, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Replace selected text with new text.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkText')
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
def saveApp():
    CTK.PREFS['tkTextDim'] = VARS[1].get()
    CTK.PREFS['tkTextSmoothness'] = VARS[2].get()
    CTK.PREFS['tkTextFont'] = VARS[3].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[1].set('3D')
    VARS[2].set('Regular')
    VARS[3].set('vera')
    CTK.PREFS['tkTextDim'] = VARS[1].get()
    CTK.PREFS['tkTextSmoothness'] = VARS[2].get()
    CTK.PREFS['tkTextFont'] = VARS[3].get()
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
    (win, menu, file, tools) = CTK.minimal('tkText '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
