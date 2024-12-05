# - tkBasicSurfs -
"""Create basic surfaces."""
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.CPlot as CP
import Geom.PyTree as D 
import Transform.PyTree as T
import Generator.PyTree as G
import Converter.Internal as Internal
import math
from Geom.Parametrics import base

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Maille un triangle en structure (decoupe en 3 quads)
#==============================================================================
def meshTri(P0, P1, P2, N):
    C01 = (0.5*(P0[0]+P1[0]), 0.5*(P0[1]+P1[1]), 0.5*(P0[2]+P1[2]))
    C12 = (0.5*(P1[0]+P2[0]), 0.5*(P1[1]+P2[1]), 0.5*(P1[2]+P2[2]))
    C02 = (0.5*(P0[0]+P2[0]), 0.5*(P0[1]+P2[1]), 0.5*(P0[2]+P2[2]))
    C = (1./3.*(P0[0]+P1[0]+P2[0]),
         1./3.*(P0[1]+P1[1]+P2[1]),
         1./3.*(P0[2]+P1[2]+P2[2]))

    l1 = D.line(P0, C01, N)
    l2 = D.line(C01, C, N)
    l3 = D.line(C, C02, N)
    l4 = D.line(C02, P0, N)
    m1 = G.TFI([l1, l2, l3, l4])
    m1 = T.reorder(m1, (-1,2,3))

    l1 = D.line(C01, P1, N)
    l2 = D.line(P1, C12, N)
    l3 = D.line(C12, C, N)
    l4 = D.line(C, C01, N)
    m2 = G.TFI([l1, l2, l3, l4])
    m2 = T.reorder(m2, (-1,2,3))

    l1 = D.line(C, C12, N)
    l2 = D.line(C12, P2, N)
    l3 = D.line(P2, C02, N)
    l4 = D.line(C02, C, N)
    m3 = G.TFI([l1, l2, l3, l4])
    m3 = T.reorder(m3, (-1,2,3))    
    return [m1, m2, m3]

#==============================================================================
# Maille un cercle en structure (decoupe en 5 quads)
#==============================================================================
def meshCircle(center, R, N):
    coeff = R*math.sqrt(2.)*0.25
    x = center[0]; y = center[1]; z = center[2]
    c = D.circle( center, R, tetas=-45., tetae=45., N=N)
    l1 = D.line((x+coeff,y-coeff,z), (x+coeff,y+coeff,z), N=N)
    l2 = D.line((x+coeff,y-coeff,z), (x+2*coeff,y-2*coeff,z), N=N)
    l3 = D.line((x+coeff,y+coeff,z), (x+2*coeff,y+2*coeff,z), N=N)
    m1 = G.TFI([c, l1, l2, l3])

    c = D.circle( center, R, tetas=45., tetae=45.+90., N=N)
    l1 = D.line((x+coeff,y+coeff,z), (x-coeff,y+coeff,z), N=N)
    l2 = D.line((x+coeff,y+coeff,z), (x+2*coeff,y+2*coeff,z), N=N)
    l3 = D.line((x-coeff,y+coeff,z), (x-2*coeff,y+2*coeff,z), N=N)
    m2 = G.TFI([c, l1, l2, l3])

    c = D.circle( center, R, tetas=45.+90, tetae=45.+180., N=N)
    l1 = D.line((x-coeff,y+coeff,z), (x-coeff,y-coeff,z), N=N)
    l2 = D.line((x-coeff,y+coeff,z), (x-2*coeff,y+2*coeff,z), N=N)
    l3 = D.line((x-coeff,y-coeff,z), (x-2*coeff,y-2*coeff,z), N=N)
    m3 = G.TFI([c, l1, l2, l3])

    c = D.circle( center, R, tetas=45.+180, tetae=45.+270., N=N)
    l1 = D.line((x-coeff,y-coeff,z), (x+coeff,y-coeff,z), N=N)
    l2 = D.line((x-coeff,y-coeff,z), (x-2*coeff,y-2*coeff,z), N=N)
    l3 = D.line((x+coeff,y-coeff,z), (x+2*coeff,y-2*coeff,z), N=N)
    m4 = G.TFI([c, l1, l2, l3])

    h = 2*coeff/(N-1)
    m5 = G.cart((x-coeff,y-coeff,z), (h,h,h), (N, N, 1))
    m5 = T.reorder(m5, (-1,2,3))
    return [m1,m2,m3,m4,m5]

#==============================================================================
def generate(event=None):
    CTK.setCursor(2, WIDGETS['generate'], WIDGETS['Npts'])
    CTK.saveTree()
    N = CTK.varsFromWidget(VARS[0].get(), type=2)
    if len(N) != 1:
        CTK.TXT.insert('START', 'NPts is incorrect.\n'); return
    N = N[0]
    eltType = VARS[1].get()
    surfType = VARS[2].get()
    if surfType == 'Sphere':
        s = D.sphere6((0,0,0), 0.5, N=N)
        xc = 0; yc = 0; zc = 0
    elif surfType == 'Plane':
        h = 1./(N-1)
        s1 = G.cart( (-0.5,-0.5,-0.5), (h,h,h), (N, 1, N) )
        s = [s1]
        xc = 0; yc = 0; zc = 0
    elif surfType == 'Cube':
        h = 1./(N-1)
        s1 = G.cart( (-0.5,-0.5,-0.5), (h,h,h), (N, N, 1) )
        s1 = T.reorder(s1, (-1,2,3))
        s2 = G.cart( (-0.5,-0.5,0.5), (h,h,h), (N, N, 1) )
        s3 = G.cart( (-0.5,-0.5,-0.5), (h,h,h), (N, 1, N) )
        s4 = G.cart( (-0.5,0.5,-0.5), (h,h,h), (N, 1, N) )
        s4 = T.reorder(s4, (-1,2,3))
        s5 = G.cart( (-0.5,-0.5,-0.5), (h,h,h), (1, N, N) )
        s5 = T.reorder(s5, (1,-2,3))
        s6 = G.cart( (0.5,-0.5,-0.5), (h,h,h), (1, N, N) )
        s = [s1, s2, s3, s4, s5, s6]
        xc = 0; yc = 0; zc = 0
    elif surfType == 'Tetra':
        m1 = meshTri([0,0,0], [1,0,0], [0,1,0], N=N)
        m1 = T.reorder(m1, (-1,2,3))
        m2 = meshTri([0,0,0], [1,0,0], [0,0,1], N=N)
        m3 = meshTri([0,0,0], [0,1,0], [0,0,1], N=N)
        m3 = T.reorder(m3, (-1,2,3))
        m4 = meshTri([1,0,0], [0,1,0], [0,0,1], N=N)
        s = m1 + m2 + m3 + m4
        xc = 0.5; yc = 0.5; zc = 0.5
    elif surfType == 'Pyramid':
        h = 1./(2*N-2)
        m0 = G.cart( (-0.5,-0.5,-0.5), (h,h,h), (2*N-1, 2*N-1, 1) )
        m0 = T.reorder(m0, (-1,2,3))
        m1 = meshTri([-0.5,-0.5,-0.5], [0.5,-0.5,-0.5], [0,0,0.5], N=N)
        m2 = meshTri([-0.5,-0.5,-0.5], [-0.5,0.5,-0.5], [0,0,0.5], N=N)
        m2 = T.reorder(m2, (-1,2,3))
        m3 = meshTri([-0.5,0.5,-0.5], [0.5,0.5,-0.5], [0,0,0.5], N=N)
        m3 = T.reorder(m3, (-1,2,3))
        m4 = meshTri([0.5,-0.5,-0.5], [0.5,0.5,-0.5], [0,0,0.5], N=N)
        s = [m0] + m1 + m2 + m3 + m4
        xc = 0.; yc = 0.; zc = 0.
    elif surfType == 'Cylinder':
        m0 = meshCircle((0,0,-0.5), 0.5, N)
        m1 = meshCircle((0,0,0.5), 0.5, N)
        m1 = T.reorder(m1, (-1,2,3))
        m2 = D.circle( (0,0,-0.5), 0.5, tetas=-45, tetae=-45+360, N=4*N-3)
        l = D.line((0,0,-0.5), (0,0,0.5), N=N)
        m2 = D.lineDrive(m2, l)
        s = m0 + m1 + [m2]
        xc = 0.; yc = 0.; zc = 0.
    elif surfType == 'Cone':
        s = [D.cone((0.,0,0), 1, 0.1, 1, N=N)]
        (xc, yc, zc) = G.barycenter(s)
    else: # Geom parametrics surfaces
        formula = base[surfType]
        if formula.replace('{u}', '') == formula: # curve
            s = D.curve(base[surfType], N)
        else:
            s = D.surface(base[surfType], N)
        (xc, yc, zc) = G.barycenter(s)
        s = [s]

    if eltType == 'TRI':
        s = C.convertArray2Tetra(s)
        s = T.join(s); s = G.close(s)
    elif eltType == 'QUAD':
        s = C.convertArray2Hexa(s)
        s = T.join(s); s = G.close(s)

    posCam = CPlot.getState('posCam')
    posEye = CPlot.getState('posEye')
    dirCam = CPlot.getState('dirCam')

    s = T.translate(s, (posEye[0]-xc, posEye[1]-yc, posEye[2]-zc) )
    lx = posEye[0]-posCam[0]
    ly = posEye[1]-posCam[1]
    lz = posEye[2]-posCam[2]
    if lx*lx + ly*ly + lz*lz < 1.e-10: lx = -1
    if dirCam[0]*dirCam[0] + dirCam[1]*dirCam[1] + dirCam[2]*dirCam[2] == 0.:
        dirCam = (0,0,1)
    ll = math.sqrt(lx*lx + ly*ly + lz*lz)
    s = T.homothety(s, (posEye[0], posEye[1], posEye[2]), 0.5*ll)

    ux = dirCam[1]*lz - dirCam[2]*ly
    uy = dirCam[2]*lx - dirCam[0]*lz
    uz = dirCam[0]*ly - dirCam[1]*lx
    s = T.rotate(s, (posEye[0], posEye[1], posEye[2]),
                 ((1,0,0),(0,1,0),(0,0,1)),
                 ((-ux,-uy,-uz), (lx,ly,lz), dirCam))

    CTK.t = C.addBase2PyTree(CTK.t, 'SURFACES', 2)
    b = Internal.getNodeFromName1(CTK.t, 'SURFACES')

    if eltType == 'TRI' or eltType == 'QUAD':
        nob = C.getNobOfBase(b, CTK.t)
        CTK.add(CTK.t, nob, -1, s)
    else:
        nob = C.getNobOfBase(b, CTK.t)
        if CP.__slot__ is None: 
            CTK.t[2][nob][2] += s; CTK.display(CTK.t)
        else:
            for i in s: CTK.add(CTK.t, nob, -1, i)

    CTK.TXT.insert('START', 'Surface created.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()
    CTK.setCursor(0, WIDGETS['generate'], WIDGETS['Npts'])

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    ttk = CTK.importTtk()

    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkBasicSurfs  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Create basic surfaces.\nCtrl+w to close applet.', btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(0, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkBasicSurfs')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- NPts -
    V = TK.StringVar(win); V.set('10'); VARS.append(V)
    if 'tkBasicSurfsNpts' in CTK.PREFS: 
        V.set(CTK.PREFS['tkBasicSurfsNpts'])
    # -1- Type d'elements
    V = TK.StringVar(win); V.set('TRI'); VARS.append(V)
    if 'tkBasicSurfsElts' in CTK.PREFS: 
        V.set(CTK.PREFS['tkBasicSurfsElts'])
    # -2- Type de surface
    V = TK.StringVar(win); V.set('Sphere'); VARS.append(V)
    if 'tkBasicSurfsType' in CTK.PREFS: 
        V.set(CTK.PREFS['tkBasicSurfsType'])

    # - Npts -
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    B.grid(row=0, column=0, sticky=TK.EW)
    WIDGETS['Npts'] = B
    BB = CTK.infoBulle(parent=B, text='Number of generated points.')
    B.bind('<Return>', generate)

    # - Type d'elements: TRI ou STRUCT -
    B = TTK.OptionMenu(Frame, VARS[1], 'TRI', 'QUAD', 'STRUCT')
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Surface element type.')

    # - Type de surface -
    SURFTYPES = ['Sphere', 'Cube', 'Tetra', 'Pyramid', 'Cylinder', 'Plane', 'Cone']
    baseKeys = list(base.keys())
    baseKeys.sort(key=str.lower)
    SURFTYPES += baseKeys
    if ttk is None:
        B = TK.OptionMenu(Frame, VARS[2], *SURFTYPES)
    else:
        B = TTK.Combobox(Frame, textvariable=VARS[2], values=SURFTYPES, 
                         state='readonly', width=10)

    B.grid(row=1, column=0, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Type of generated surface.')
    B = TTK.Button(Frame, text="Generate", command=generate)
    WIDGETS['generate'] = B
    B.grid(row=1, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Generate surface.')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['SurfNoteBook'].add(WIDGETS['frame'], text='tkBasicSurfs')
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
    CTK.PREFS['tkBasicSurfsNpts'] = VARS[0].get()
    CTK.PREFS['tkBasicSurfsElts'] = VARS[1].get()
    CTK.PREFS['tkBasicSurfsType'] = VARS[2].get()
    CTK.savePrefFile()

#==============================================================================
def resetApp():
    VARS[0].set('10')
    VARS[1].set('TRI')
    VARS[2].set('Sphere')
    CTK.PREFS['tkBasicSurfsNpts'] = VARS[0].get()
    CTK.PREFS['tkBasicSurfsElts'] = VARS[1].get()
    CTK.PREFS['tkBasicSurfsType'] = VARS[2].get()
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
    (win, menu, file, tools) = CTK.minimal('tkBasicSurfs '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
