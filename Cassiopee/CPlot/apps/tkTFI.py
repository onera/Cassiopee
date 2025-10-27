# - tkTFI -
"""Transfinite interpolation mesher."""
import tkinter as TK
import CPlot.Ttk as TTK
import Converter.PyTree as C
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import CPlot.Panels as Panels
import CPlot.iconics as iconics
import Converter.Internal as Internal
import Generator.PyTree as G
import Generator.TFIs as TFIs
import Transform.PyTree as T
import Converter
import Generator

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
def setSurface():
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

#=============================================================================
def getSurfaces():
    name = VARS[0].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/')
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            zones = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in zones:
                if z[0] == sname[1]: surfaces.append(z)
    return surfaces

#==============================================================================
def TFI():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) == 0:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    surf = getSurfaces()

    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        dim = Internal.getZoneDim(z)
        if dim[3] == 'BAR':
            zp = C.convertBAR2Struct(z); zones.append(zp)
        else: zones.append(z)

    try:
        CTK.saveTree()
        mesh = G.TFI(zones)
        if surf != []: mesh = T.projectOrthoSmooth(mesh, surf)
        CTK.t = C.addBase2PyTree(CTK.t, 'MESHES')
        bases = Internal.getNodesFromName1(CTK.t, 'MESHES')
        nob = C.getNobOfBase(bases[0], CTK.t)
        CTK.add(CTK.t, nob, -1, mesh)
        CTK.TXT.insert('START', 'TFI mesh created.\n')
        #C._fillMissingVariables(CTK.t)
        (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
        CTK.TKTREE.updateApp()
        CPlot.render()
    except Exception as e:
        Panels.displayErrors([0,str(e)], header='Error: TFI')
        CTK.TXT.insert('START', 'TFI mesh failed.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')

#==============================================================================
# IN: a1,a2,a3: les 3 cotes du triangle
#==============================================================================
def trimesh(a1, a2, a3):
    N1 = a1[2]; N2 = a2[2]; N3 = a3[2]

    # Verif de N
    Nt = N3-N2+N1+1
    if Nt//2-Nt*0.5 != 0: return [0, 'N3-N2+N1 must be odd.', 0]
    N = Nt//2
    if N < 2: return [0, 'invalid number of points for this operation.', 0]
    if N > N1-1: return [0, 'invalid number of points for this operation.',0]
    if N > N2-1: return [0, 'invalid number of points for this operation.',0]

    return TFIs.TFITri(a1, a2, a3)

#==============================================================================
# Cree un seul maillage a partir de 2 courbes
# a1: straight, a2: round
# Il faut N1-N2 pair.
# Celle qui a le plus de points est prise pour round.
#==============================================================================
def mono2mesh(a1, a2):
    N1 = a1[2]; N2 = a2[2]
    diff = N2-N1
    if diff//2 != diff*0.5: return ['N1-N2 must be even.']
    return TFIs.TFIMono(a1, a2)

#==============================================================================
# Cree un seul maillage a partir de 3 courbes
# a1: straight, a2, a3: round
# Celle qui a le plus de points est prise pour round
#==============================================================================
def mono1mesh(a1, a2, a3):
    return []

#==============================================================================
# Evalue la qualite du maillage m
# Plus le score est elevee pour le maillage est mauvais
#==============================================================================
def quality(meshes):
    score = 0.
    for m in meshes:
        ortho = Generator.getOrthogonalityMap(m)
        vol = Generator.getVolumeMap(m)
        min1 = Converter.getMinValue(ortho, 'orthogonality')
        max1 = Converter.getMaxValue(ortho, 'orthogonality')
        min2 = Converter.getMinValue(vol, 'vol')
        max2 = Converter.getMaxValue(vol, 'vol')
        score = max(score, min1); score = max(score, max1)
        if min2 < 1.e-12 and max2 > 0: score += 1000.
        elif max2 < 1.e-12 and min2 > 0: score += 1000.
    return score

#==============================================================================
def OTFI():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) == 0:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    surf = getSurfaces()

    CTK.setCursor(2, WIDGETS['o'])

    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = C.convertArray2Hexa(z)
        zones.append(z)
    zones = T.join(zones); zones = G.close(zones)
    a = C.convertBAR2Struct(z)

    weight = CTK.varsFromWidget(VARS[1].get(), type=1); weight = weight[0]

    # Nombre de pts
    Nt = Internal.getZoneDim(a)[1]
    if Nt//2 - Nt*0.5 == 0:
        CTK.TXT.insert('START', 'Number of points of countour must be odd.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error')
        CTK.setCursor(0, WIDGETS['o'])
        return

    coords = C.getFields(Internal.__GridCoordinates__, a, api=1)[0]

    optWeight = 0; optOffset = 0; optScore = 1.e6
    Nt = coords[2]
    for j in range(-Nt//4,Nt//4+1):
        for i in range(3,10):
            try:
                [m,m1,m2,m3,m4] = TFIs.TFIO__(coords, i, j)
                score = quality([m,m1,m2,m3])
                if score < optScore:
                    optWeight = i; optOffset = j; optScore = score
            except: pass
    print('resulting weight=%g, offset=%g.'%(optWeight,optOffset))
    print('resulting score=%g.'%optScore)
    [m,m1,m2,m3,m4] = TFIs.TFIO__(coords, optWeight, optOffset)

    m = C.convertArrays2ZoneNode('TFI1', [m])
    m1 = C.convertArrays2ZoneNode('TFI2', [m1])
    m2 = C.convertArrays2ZoneNode('TFI3', [m2])
    m3 = C.convertArrays2ZoneNode('TFI4', [m3])
    m4 = C.convertArrays2ZoneNode('TFI5', [m4])

    if surf != []:
        m = T.projectOrtho(m, surf)
        m1 = T.projectOrthoSmooth(m1, surf)
        m2 = T.projectOrthoSmooth(m2, surf)
        m3 = T.projectOrthoSmooth(m3, surf)
        m4 = T.projectOrthoSmooth(m4, surf)

    CTK.saveTree()
    CTK.setCursor(0, WIDGETS['o'])

    CTK.t = C.addBase2PyTree(CTK.t, 'MESHES')
    bases = Internal.getNodesFromName1(CTK.t, 'MESHES')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in [m,m1,m2,m3,m4]: CTK.add(CTK.t, nob, -1, i)
    CTK.TXT.insert('START', 'O-TFI mesh created.\n')

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def HOTFI():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) == 0:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    surf = getSurfaces()

    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = C.convertBAR2Struct(z)
        zones.append(z)

    if len(zones) != 2:
        CTK.TXT.insert('START', 'HO TFI takes 2 contours.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    weight = CTK.varsFromWidget(VARS[1].get(), type=1); weight = weight[0]

    # Nombre de pts (tous les 2 pairs ou tous les 2 impairs)
    Nt1 = Internal.getZoneDim(zones[0])[1]
    Nt2 = Internal.getZoneDim(zones[1])[1]
    if Nt1//2 - Nt1*0.5 == 0 and Nt2//2 - Nt2*0.5 != 0:
        CTK.TXT.insert('START', 'Number of points of countours must be all odd or all even.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    if Nt1//2 - Nt1*0.5 != 0 and Nt2//2 - Nt2*0.5 == 0:
        CTK.TXT.insert('START', 'Number of points of countours must be all odd or all even.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    CTK.setCursor(2, WIDGETS['ho'])

    coords1 = C.getFields(Internal.__GridCoordinates__, zones[0], api=1)[0]
    coords2 = C.getFields(Internal.__GridCoordinates__, zones[1], api=1)[0]

    optWeight = 0; optOffset = 0; optScore = 1.e6
    Nt2 = coords2[2]
    for j in range(-Nt2//8,Nt2//8):
        for i in range(2,10):
            try:
                [m,m1,m2,m3] = TFIs.TFIHalfO__(coords1, coords2, i, j)
                score = quality([m,m1,m2,m3])
                if score < optScore:
                    optWeight = i; optScore = score; optOffset = j
            except: pass
    print('Resulting score=%g'%optScore)
    [m,m1,m2,m3] = TFIs.TFIHalfO__(coords1, coords2, optWeight, optOffset)

    m = C.convertArrays2ZoneNode('TFI1', [m])
    m1 = C.convertArrays2ZoneNode('TFI2', [m1])
    m2 = C.convertArrays2ZoneNode('TFI3', [m2])
    m3 = C.convertArrays2ZoneNode('TFI4', [m3])

    if surf != []:
        m = T.projectOrthoSmooth(m, surf)
        m1 = T.projectOrthoSmooth(m1, surf)
        m2 = T.projectOrthoSmooth(m2, surf)
        m3 = T.projectOrthoSmooth(m3, surf)

    CTK.saveTree()
    CTK.setCursor(0, WIDGETS['ho'])

    CTK.t = C.addBase2PyTree(CTK.t, 'MESHES')
    bases = Internal.getNodesFromName1(CTK.t, 'MESHES')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in [m,m1,m2,m3]: CTK.add(CTK.t, nob, -1, i)
    CTK.TXT.insert('START', 'HO-TFI mesh created.\n')

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def TRITFI():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) == 0:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    surf = getSurfaces()

    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = C.convertBAR2Struct(z)
        zones.append(z)

    if len(zones) != 3:
        CTK.TXT.insert('START', 'TRI TFI takes 3 contours.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    coords1 = C.getFields(Internal.__GridCoordinates__, zones[0], api=1)[0]
    coords2 = C.getFields(Internal.__GridCoordinates__, zones[1], api=1)[0]
    coords3 = C.getFields(Internal.__GridCoordinates__, zones[2], api=1)[0]

    [m1,m2,m3] = trimesh(coords1, coords2, coords3)
    if m1 == 0:
        CTK.TXT.insert('START', m2+'\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    m1 = C.convertArrays2ZoneNode('TFI1', [m1])
    m2 = C.convertArrays2ZoneNode('TFI2', [m2])
    m3 = C.convertArrays2ZoneNode('TFI3', [m3])

    if surf != []:
        m1 = T.projectOrthoSmooth(m1, surf)
        m2 = T.projectOrthoSmooth(m2, surf)
        m3 = T.projectOrthoSmooth(m3, surf)

    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'MESHES')
    bases = Internal.getNodesFromName1(CTK.t, 'MESHES')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in [m1,m2,m3]: CTK.add(CTK.t, nob, -1, i)
    CTK.TXT.insert('START', 'TRI-TFI mesh created.\n')

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
def MONO2TFI():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if len(nzs) == 0:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    surf = getSurfaces()

    zones = []
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        z = C.convertBAR2Struct(z)
        zones.append(z)

    if len(zones) != 2:
        CTK.TXT.insert('START', 'MONO2 TFI takes 2 contours.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    coords1 = C.getFields(Internal.__GridCoordinates__, zones[0], api=1)[0]
    coords2 = C.getFields(Internal.__GridCoordinates__, zones[1], api=1)[0]

    [m] = mono2mesh(coords1, coords2)
    if isinstance(m, str):
        CTK.TXT.insert('START', m+'\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    m = C.convertArrays2ZoneNode('TFI', [m])
    if surf != []: m = T.projectOrthoSmooth(m, surf)

    CTK.saveTree()
    CTK.t = C.addBase2PyTree(CTK.t, 'MESHES')
    bases = Internal.getNodesFromName1(CTK.t, 'MESHES')
    nob = C.getNobOfBase(bases[0], CTK.t)
    CTK.add(CTK.t, nob, -1, m)
    CTK.TXT.insert('START', 'HO-TFI mesh created.\n')

    #C._fillMissingVariables(CTK.t)
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkTFI  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Transfinite interpolation mesh generation.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=2)
    Frame.columnconfigure(1, weight=2)
    Frame.columnconfigure(2, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkTFI')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- underlaying surface
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -1- weight
    V = TK.StringVar(win); V.set('10.'); VARS.append(V)

    # - Model surface -
    B = TTK.Button(Frame, text="Surf", command=setSurface,
                   image=iconics.PHOTO[8], padx=0, pady=0, compound=TK.RIGHT)
    B.grid(row=0, column=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Underlaying model surfaces.')
    B = TTK.Entry(Frame, textvariable=VARS[0], background='White')
    BB = CTK.infoBulle(parent=B, text='Underlaying model surfaces.')
    B.grid(row=0, column=0, columnspan=2, sticky=TK.EW)

    # - O - TFI -
    B = TTK.Button(Frame, text="O", command=OTFI)
    WIDGETS['o'] = B
    B.grid(row=1, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Perform a O-TFI from a close contour (npts must be odd).')
    # - HO - TFI
    B = TTK.Button(Frame, text="Half-O", command=HOTFI)
    WIDGETS['ho'] = B
    B.grid(row=1, column=1, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Perform a Half O-TFI from 2 contours (longer contour is taken to be the round part).')
    # - TRI - TFI -
    B = TTK.Button(Frame, text="TRI", command=TRITFI)
    B.grid(row=1, column=2, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Perform a TRI-TFI from 3 contours.')
    #B = TK.Entry(Frame, textvariable=VARS[1], background='White', width=5)
    #B.grid(row=1, column=2, sticky=TK.EW)
    #BB = CTK.infoBulle(parent=B, text='Weight of curve.')

    # - MONO2 - TFI -
    B = TTK.Button(Frame, text="MONO2", command=MONO2TFI)
    B.grid(row=2, column=0, columnspan=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Perform a TFI from two contours (longer contour is split).')

    # - TFI des cotes -
    B = TTK.Button(Frame, text="TFI", command=TFI)
    B.grid(row=2, column=1, columnspan=2, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Perform TFI from 4 edges (2D) or 6 faces (3D).')

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    #WIDGETS['frame'].grid(sticky=TK.NSEW)
    try: CTK.WIDGETS['MeshNoteBook'].add(WIDGETS['frame'], text='tkTFI')
    except: pass
    CTK.WIDGETS['MeshNoteBook'].select(WIDGETS['frame'])

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    #WIDGETS['frame'].grid_forget()
    CTK.WIDGETS['MeshNoteBook'].hide(WIDGETS['frame'])

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
    (win, menu, file, tools) = CTK.minimal('tkTFI '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
