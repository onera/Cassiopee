# - surface filter and offset -
try: import Tkinter as TK
except: import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.Internal as Internal
import Generator.PyTree as G
import Post.PyTree as P
import Dist2Walls.PyTree

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Calcul la distance a la paroi par differents algos
# IN: b: maillage ou l'on calcule la distance
# IN: a: maillage du corps
# IN: loc='centers', 'nodes': localisation du champ de distance
#==============================================================================
def compDistance(b, a, loc):
    try: 
        b = Dist2Walls.PyTree.distance2Walls(b, a, type='ortho', 
                                             loc=loc, signed=1)
        fail = False
    except: fail = True
    if fail:
        try:
            b = Dist2Walls.PyTree.distance2Walls(b, a, type='ortho', 
                                                 loc=loc, signed=0)
            fail = False
        except: fail = True
    if fail:
        try:
            b = Dist2Walls.PyTree.distance2Walls(b, a, type='mininterf', 
                                                 loc=loc, signed=0)
            fail = False
        except: raise
    return b

#==============================================================================
def withCart(a, offset, density):
    # Calcul la grille cartesienne
    BB = G.bbox(a)
    xmin = BB[0]; ymin = BB[1]; zmin = BB[2]
    xmax = BB[3]; ymax = BB[4]; zmax = BB[5]
    ni = density*(xmax-xmin); nj = density*(ymax-ymin);
    nk = density*(zmax-zmin)
    if ni < 2: ni = 2
    if nj < 2: nj = 2
    if nk < 2: nk = 2
    hi = (xmax-xmin)/(ni-1); hj = (ymax-ymin)/(nj-1); hk = (zmax-zmin)/(nk-1)
    h = min(hi, hj); h = min(h, hk); h = max(h, 1.e-6)
    ni = int((xmax-xmin)/h)+7; nj = int((ymax-ymin)/h)+7
    nk = int((zmax-zmin)/h)+7
    ni += int(2*offset/h); nj += int(2*offset/h); nk += int(2*offset/h)
    b = G.cart( (xmin-3*h-offset, ymin-3*h-offset, zmin-3*h-offset), (h, h, h), (ni,nj,nk) )

    # Calcul la distance a la paroi
    b = compDistance(b, a, loc='nodes')
    #C.convertPyTree2File(b, 'out.cgns')

    # Extraction isoSurf
    iso = P.isoSurfMC([b], 'TurbulentDistance', value=offset)
    return iso

#==============================================================================
def withOctree(a, offset, density):
    # step
    tol = 1./density
    
    # octree
    snears = []; sec = 0
    for z in a:
        bb = G.bbox(z)
        rx = bb[3]-bb[0]; ry = bb[4]-bb[1]; rz = bb[5]-bb[2]
        snear = min(rx, ry); snear = min(snear, rz)
        snear = 0.1*snear
        sec = max(sec, snear)
        snears.append(snear)
    o = G.octree(a, snears, dfar=offset+sec)
    o = compDistance(o, a, loc='nodes')

    # iteration d'adaptation
    nit = 0
    while nit < 10:
        print('iterating: %d...'%nit)

        o = C.node2Center(o, 'TurbulentDistance')
        o = G.getVolumeMap(o)

        # adapt
        C._initVars(o, '{centers:vol}={centers:vol}**0.33333')
        # was 2.1 factor
        C._initVars(o, '{centers:indicator}=logical_and({centers:vol} > %20.16g , abs({centers:TurbulentDistance}-%20.16g) < 1.*{centers:vol})'%(tol,offset))
        o1 = G.adaptOctree(o, 'centers:indicator')

        #C.convertPyTree2File([o1]+a, 'out%d.cgns'%nit)

        # check convergence
        dim1 = Internal.getZoneDim(o1); dim = Internal.getZoneDim(o)
        if dim1 == dim: break

        #if (nit%2 == 0): o1 = P.extractMesh([o], o1)
        o1 = compDistance(o1, a, loc='nodes')
        o = o1
        nit += 1
    
    o = C.rmVars(o, 'centers:TurbulentDistance')
    o = C.rmVars(o, 'centers:vol')
    o = C.rmVars(o, 'centers:indicator')
    #C.convertPyTree2File(o, 'out.cgns')

    # Iso surface
    iso = P.isoSurfMC([o], 'TurbulentDistance', value=offset)
    return iso

#==============================================================================
def remap():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    a = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        a.append(z)

    # density
    density = CTK.varsFromWidget(VARS[0].get(), type=1)
    if len(density) != 1:
        CTK.TXT.insert('START', 'Density is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    density = density[0]
    
    # offset
    offset = CTK.varsFromWidget(VARS[1].get(), type=1)
    if len(offset) != 1:
        CTK.TXT.insert('START', 'Offset is incorrect.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    offset = offset[0]

    CTK.saveTree()

    if VARS[2].get() == '0': iso = withCart(a, offset, density)
    else: iso = withOctree(a, offset, density)
        
    if iso != []:
        nob = CTK.Nb[nzs[0]]+1
        for i in iso: CTK.add(CTK.t, nob, -1, i)
    
        #C._fillMissingVariables(CTK.t)
        CTK.TXT.insert('START', 'Surface filtered and offset (offset=%g).\n'%offset)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
    else:
         CTK.TXT.insert('START', 'Surface filter failed.\n')
         CTK.TXT.insert('START', 'Error: ', 'Error')     
    
#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkFilterSurfs', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Filter or offset a surface.\nCtrl+c to close applet.', temps=0, btype=1)
    Frame.bind('<Control-c>', hideApp)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=1)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=2)
    Frame.columnconfigure(3, weight=0)
    WIDGETS['frame'] = Frame
    
    # - Frame menu -
    FrameMenu = TK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+c', command=hideApp)
    FrameMenu.add_command(label='Save', command=saveApp)
    FrameMenu.add_command(label='Reset', command=resetApp)
    CTK.addPinMenu(FrameMenu, 'tkFilterSurfs')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Point density -
    V = TK.StringVar(win); V.set('1.'); VARS.append(V)
    if 'tkFilterSurfsDensity' in CTK.PREFS: 
        V.set(CTK.PREFS['tkFilterSurfsDensity'])
    # -1- Offset -
    V = TK.StringVar(win); V.set('0.'); VARS.append(V)
    if 'tkFilterSurfsOffset' in CTK.PREFS: 
        V.set(CTK.PREFS['tkFilterSurfsOffset'])
    # -2- Algorithm cart/octree
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    if 'tkFilterSurfsType' in CTK.PREFS: 
        V.set(CTK.PREFS['tkFilterSurfsType'])

    # - Point density -
    B = TTK.Label(Frame, text="density/offset")
    B.grid(row=0, column=0, sticky=TK.EW)
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White', width=8)
    B.grid(row=0, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Point density (npts/length unit).')

    # - Offset -
    B = TTK.Entry(Frame, textvariable=VARS[1], background='White', width=8)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Distance of surface offset.')

    # - Algorithm -
    B = TTK.Checkbutton(Frame, text='', variable=VARS[2])
    B.grid(row=0, column=3, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Toggle octree algorithm.')
    
    # - Remap -
    B = TTK.Button(Frame, text="Filter/offset", command=remap)
    BB = CTK.infoBulle(parent=B, text='Remap/offset a surface with a given spacing.')
    B.grid(row=2, column=0, columnspan=4, sticky=TK.EW)
    
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
    CTK.PREFS['tkFilterSurfsDensity'] = VARS[0].get()
    CTK.PREFS['tkFilterSurfsOffset'] = VARS[1].get()
    CTK.PREFS['tkFilterSurfsType'] = VARS[2].get()
    CTK.savePrefFile()
    
#==============================================================================
def resetApp():
    VARS[0].set('1.')
    VARS[1].set('0.')
    VARS[2].set('0')
    CTK.PREFS['tkFilterSurfsDensity'] = VARS[0].get()
    CTK.PREFS['tkFilterSurfsOffset'] = VARS[1].get()
    CTK.PREFS['tkFilterSurfsType'] = VARS[1].get()
    CTK.savePrefFile()

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
    (win, menu, file, tools) = CTK.minimal('tkFilterSurfs '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
