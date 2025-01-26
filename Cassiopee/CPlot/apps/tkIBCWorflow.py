# - IBC app -
try: import tkinter as TK
except: import Tkinter as TK
import CPlot.Ttk as TTK
import CPlot.PyTree as CPlot
import CPlot.Tk as CTK
import Converter.PyTree as C
import Connector.PyTree as X
import Converter.Internal as Internal
import Dist2Walls.PyTree as DTW
import Generator.PyTree as G
import Post.PyTree as P
import Transform.PyTree as T
import Converter
import numpy

# local widgets list
WIDGETS = {}; VARS = []

#==============================================================================
# Creation des premiers pts IBC
# IN: a: zone avec cellN=2 pour les pts IBC
# OUT: retourne les pts IBC sous forme de 'NODE'
#==============================================================================
def getIBCFrontForZone__(a):
    f0 =  P.selectCells(a, '{centers:cellN} == 2.')
    f0 = T.join(f0); f0 = G.close(f0)

    # recuperation des champs en centres perdus par selectCells
    a2 = C.initVars(a, 'centers:cellN', 1.)
    ta = C.newPyTree(['Base', a2])
    tf0 = C.newPyTree(['Base', f0])
    tf0 = P.extractMesh(ta, tf0)
    tf0 = C.rmVars(tf0,['centers:cellN','cellN'])
    coords = C.getFields(Internal.__GridCoordinates__,tf0)
    solc =  C.getFields(Internal.__FlowSolutionCenters__,tf0)
    coords = Converter.node2Center(coords) # passage en centres pour recuperer les coordonnees des centres
    coords = Converter.addVars([coords,solc])
    # passage en 'NODE'
    coords = Converter.convertArray2Node(coords)[0]
    # Sortie : conversion en zone node CGNS
    return C.convertArrays2ZoneNode('Front_'+a[0], [coords])

#==============================================================================
# Calcul de la normale et de delta pour les pts du front
# IN: fc1: front IBC sous forme de zone
# IN: parentZone: zone dont provient la zone IBC, localisee comme
# souhaitee pour les interpolations (ex en centres avec cell fict)
# OUT: front fc1 avec delta
#          la normale = delta * grad TurbulentDistance
#          l'indice global dans le maillage parent correspondant au pt IBC
#==============================================================================
def getIBCFrontInfo__(fc1, parentZone, dhloc, toldist=1.e-10):
    listPts = [] # on ne modifie pas localement le maillage
    eps = 0. # decalage de distance pour listpts

    # Determination de l indice du pt dans le maillage donneur
    hook = C.createHook(parentZone, function='nodes')
    indices,distances = C.nearestNodes(hook,fc1)
    C.freeHook(hook)

    coords1 = C.getFields(Internal.__GridCoordinates__,fc1)[0]
    Converter._initVars(coords1,'ind',-1.)
    Converter._initVars(coords1,'indI',-1.)
    Converter._initVars(coords1,'indJ',-1.)
    Converter._initVars(coords1,'indK',-1.)
    coords1 = Converter.extractVars(coords1,['ind','indI','indJ','indK'])

    npts = coords1[1].shape[1]
    dimGC = Internal.getZoneDim(parentZone)
    nigc = dimGC[1]; njgc = dimGC[2]; nkgc = dimGC[3]; nigcnjgc = nigc*njgc
    dhLoc = [dhloc]*npts # tableau dhLoc
    xt = C.getField('CoordinateX',parentZone)[0][1]
    yt = C.getField('CoordinateY',parentZone)[0][1]
    zt = C.getField('CoordinateZ',parentZone)[0][1]
    for ind in range(npts):
        index = indices[ind]-1
        dist  = distances[ind]
        if dist < toldist:
            indk = index/nigcnjgc
            indj = (index-indk*nigcnjgc)/nigc
            indi = index - indj*nigc - indk*nigcnjgc
            coords1[1][0,ind] = index
            coords1[1][1,ind] = indi
            coords1[1][2,ind] = indj
            coords1[1][3,ind] = indk
            if indi == nigc-1: dxloc = abs(xt[0,index]-xt[0,index-1])
            else: dxloc = abs(xt[0,index]-xt[0,index+1])
            if indj == njgc-1: dyloc = abs(yt[0,index]-yt[0,index-nigc])
            else: dyloc = abs(yt[0,index]-yt[0,index+nigc])
            if indk == nkgc-1: dzloc = abs(zt[0,index]-zt[0,index-nigcnjgc])
            else: dzloc = abs(zt[0,index]-zt[0,index+nigcnjgc])

            if dxloc < toldist: dxloc = 1.e10
            if dyloc < toldist: dyloc = 1.e10
            if dzloc < toldist: dzloc = 1.e10
            dhLoc[ind] = min(dxloc,dyloc,dzloc)

    C.setFields([coords1], fc1, loc='nodes')

    # delta * normale
    varnx = 'gradxTurbulentDistance'
    varny = 'gradyTurbulentDistance'
    varnz = 'gradzTurbulentDistance'
    fc1 = C.normalize(fc1, [varnx, varny, varnz])

    if listPts == []:
        C._initVars(fc1,'delta',0.)
        deltaa = C.getField('delta',fc1)[0]
        distance = C.getField('TurbulentDistance',fc1)[0][1]
        # formule d obtention du frontC2
        for ind in range(deltaa[1].shape[1]):
            dist = distance[0,ind]
            # NOUVELLE VERSION
            # # cas 1 : le centre est proche paroi, le point interpole est alors positionne a dhloc+eps de la paroi
        # if abs(dist) < dhLoc[ind]: deltaa[1][0,ind] =  2*dhLoc[ind] + eps
            # # cas 2 : le centre est loin de la paroi, le point interpole est alors positionne a dist+eps de la paroi
        # else: deltaa[1][0,ind] = 2.*abs(dist) + eps
            # FIN NOUVELLE VERSION

            # cas 1 : le centre est proche paroi, le point interpole est alors positionne a dhloc+eps de la paroi
            # cas 2 : le centre est loin de la paroi, le point interpole est alors positionne a dist+eps de la paroi
            if abs(dist) < dhloc: deltaa[1][0,ind] = abs(dist) + dhloc + eps
            else: deltaa[1][0,ind] = 2.*abs(dist) + eps
        C.setFields([deltaa], fc1, loc='nodes')

    else:
        # modification locale de delta : a modifier ?
        deltaa = C.getField('delta',fc1)[0]
        for ind in listPts: deltaa[1][0,ind] += eps
        C.setFields([deltaa], fc1, loc='nodes')

    C._initVars(fc1, '{nx} = {gradxTurbulentDistance} * {delta}')
    C._initVars(fc1, '{ny} = {gradyTurbulentDistance} * {delta}')
    C._initVars(fc1, '{nz} = {gradzTurbulentDistance} * {delta}')
    return fc1

#==============================================================================
def setSurface():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    selected = ''
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        selected += CTK.t[2][nob][0]+'/'+z[0]+';'
    selected = selected[0:-1]
    VARS[5].set(selected)

#==============================================================================
# blanking de la selection avec la surface fournie
#==============================================================================
def blank():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    # type de blanking
    type0 = VARS[0].get()
    # EquationDimension
    node = Internal.getNodeFromName(CTK.t, 'EquationDimension')
    if node is not None: dimPb = Internal.getValue(node)
    else:
        CTK.TXT.insert('START', 'EquationDimension not found (tkState). Using 3D.\n')
        CTK.TXT.insert('START', 'Warning: ', 'Warning'); dimPb = 3

    # Blanking surfaces
    name = VARS[5].get()
    names = name.split(';')
    surfaces = []
    for v in names:
        v = v.lstrip(); v = v.rstrip()
        sname = v.split('/', 1)
        bases = Internal.getNodesFromName1(CTK.t, sname[0])
        if bases != []:
            nodes = Internal.getNodesFromType1(bases[0], 'Zone_t')
            for z in nodes:
                if z[0] == sname[1]: surfaces.append(z)
    # Reglages XRay
    delta = float(VARS[2].get())
    tol = float(VARS[3].get())

    # Blank in/out ?
    # Create blanking Matrix
    BM = numpy.zeros((1, 1), dtype=Internal.E_NpyInt)
    isIn = VARS[4].get()
    if isIn == 'inside': BM[0,0] = 1
    else: BM[0,0] = -1

    # Creation de l'arbre temporaire
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    t = C.newPyTree(['Base'])
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        t[2][1][2].append(z)
    # BlankCells
    CTK.saveTree()
    depth = VARS[6].get(); depth = int(depth)
    t = X.blankCells(t, [surfaces], blankingMatrix=BM,
                     blankingType=type0,
                     delta=delta, dim=dimPb, tol=tol)
    C._initVars(t, '{centers:cellN}=minimum(1.,{centers:cellN})')
    t = X.cellN2OversetHoles(t)
    t = X.setHoleInterpolatedPoints(t, depth=-depth)

    tp = C.newPyTree(['Base']); donorNoz=[]
    for noz in range(len(t[2][1][2])):
        z = t[2][1][2][noz]
        valmax = C.getMaxValue(z, 'centers:cellN')
        if valmax == 2.: tp[2][1][2].append(z); donorNoz.append(noz)

    tp = DTW.distance2Walls(tp, surfaces, type='ortho', loc='centers', signed=1,dim=dimPb)
    tp = Internal.correctPyTree(tp, level=6)
    tp = C.center2Node(tp,'centers:TurbulentDistance')
    tp = P.computeGrad(tp, 'TurbulentDistance')
    for noz2 in range(len(donorNoz)):
        noz = donorNoz[noz2]
        t[2][1][2][noz] = tp[2][1][2][noz2]

    # Back to tree
    c = 0
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = t[2][1][2][c]
        CTK.t[2][nob][2][noz] = z
        c += 1

    C._fillMissingVariables(CTK.t)
    CTK.TXT.insert('START', 'Blanking done.\n')
    CTK.TKTREE.updateApp()
    CTK.display(CTK.t)

#==============================================================================
def getIBCFront():
    if CTK.t == []: return
    if CTK.__MAINTREE__ <= 0:
        CTK.TXT.insert('START', 'Fail on a temporary tree.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return
    nzs = CPlot.getSelectedZones()
    if nzs == []:
        CTK.TXT.insert('START', 'Selection is empty.\n')
        CTK.TXT.insert('START', 'Error: ', 'Error'); return

    front1 = []; nozDonors=[]# numerotation ds td
    td = C.newPyTree(['Donors'])
    nozloc = 0
    for nz in nzs:
        nob = CTK.Nb[nz]+1
        noz = CTK.Nz[nz]
        z = CTK.t[2][nob][2][noz]
        td[2][1][2].append(z)
        z2 = getIBCFrontForZone__(z)
        dims = Internal.getZoneDim(z2)
        if dims[1] != 0: front1.append(z2); nozDonors.append(nozloc)
        nozloc+= 1 # numerotation ds td

    # donneurs : add ghostcells ? localisation ?
    nghostcells = int(VARS[7].get())
    locD = VARS[8].get()

    if nghostcells > 0: td = Internal.addGhostCells(td, td, nghostcells, adaptBCs=0)
    if locD == 'centers': td = C.node2Center(td)
    elif locD == 'ext_centers': td = C.node2ExtCenter(td)
    CTK.saveTree()

    # dhloc : distmin a partir de laquelle on peut symetriser fc1
    dhloc = float(VARS[9].get())
    # get ind,indi,indj,indk,delta,nx,ny,nz pour les pts du front1
    for nof1 in range(len(front1)):
        fc1 = front1[nof1]
        fc1 = getIBCFrontInfo__(fc1, td[2][1][2][nozDonors[nof1]], dhloc)
        fc1 = C.rmNodes(fc1,'gradxTurbulentDistance')
        fc1 = C.rmNodes(fc1,'gradyTurbulentDistance')
        fc1 = C.rmNodes(fc1,'gradzTurbulentDistance')
        front1[nof1] = fc1

    # Front2
    front2 = T.deform(front1, ['nx','ny','nz'])

    # Interpolation
    interpType = VARS[10].get()
    interpOrder = 2
    if interpType == '2nd order': interpOrder = 2
    elif interpType == '3rd order Lagrangian': interpOrder = 3
    elif interpType == '5th order Lagrangian': interpOrder = 5
    front2 = X.setInterpData(front2, td, order=interpOrder,loc='nodes',penalty=1,nature=0)
    front2 = C.rmNodes(front2,Internal.__FlowSolutionCenters__)

    CTK.t = C.addBase2PyTree(CTK.t, 'IBCFront')
    bases = Internal.getNodesFromName1(CTK.t, 'IBCFront')
    nob = C.getNobOfBase(bases[0], CTK.t)
    for i in front2: CTK.add(CTK.t, nob, -1, i)

    C._fillMissingVariables(CTK.t)
    CTK.t = C.rmVars(CTK.t,'centers:cellN')
    CTK.TXT.insert('START', 'IBC front zones added.\n')
    (CTK.Nb, CTK.Nz) = CPlot.updateCPlotNumbering(CTK.t)
    CTK.TKTREE.updateApp()
    CPlot.render()

#==============================================================================
# Create app widgets
#==============================================================================
def createApp(win):
    # - Frame -
    Frame = TTK.LabelFrame(win, borderwidth=2, relief=CTK.FRAMESTYLE,
                           text='tkIBC  [ + ]  ', font=CTK.FRAMEFONT, takefocus=1)
    #BB = CTK.infoBulle(parent=Frame, text='Compute OBC points.\nCtrl+w to close applet.', temps=0, btype=1)
    Frame.bind('<Control-w>', hideApp)
    Frame.bind('<ButtonRelease-1>', displayFrameMenu)
    Frame.bind('<ButtonRelease-3>', displayFrameMenu)
    Frame.bind('<Enter>', lambda event : Frame.focus_set())
    Frame.columnconfigure(0, weight=0)
    Frame.columnconfigure(1, weight=1)
    WIDGETS['frame'] = Frame

    # - Frame menu -
    FrameMenu = TTK.Menu(Frame, tearoff=0)
    FrameMenu.add_command(label='Close', accelerator='Ctrl+w', command=hideApp)
    CTK.addPinMenu(FrameMenu, 'tkIBC')
    WIDGETS['frameMenu'] = FrameMenu

    # - VARS -
    # -0- Blanking type
    V = TK.StringVar(win); V.set('cell_intersect_opt'); VARS.append(V)
    # -1- Blanking surface
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -2- delta XRay -
    V = TK.StringVar(win); V.set('1.e-10'); VARS.append(V)
    # -3- tolerance XRay -
    V = TK.StringVar(win); V.set('1.e-8'); VARS.append(V)
    # -4- blank inside bodies - outside
    V = TK.StringVar(win); V.set('inside'); VARS.append(V)
    # -5- Blanking surface
    V = TK.StringVar(win); V.set(''); VARS.append(V)
    # -6- Nb de rangees de pts du front1
    V = TK.StringVar(win); V.set('2'); VARS.append(V)
    # -7- Nb de rangees de ghost cells du maillage donneur
    V = TK.StringVar(win); V.set('0'); VARS.append(V)
    # -8- Localisation du maillage donneur
    V = TK.StringVar(win); V.set('centers'); VARS.append(V)
    # -9- distance minimale a la paroi des pts du front interieur pour pouvoir construire le symetrique
    V = TK.StringVar(win); V.set('0.00973457'); VARS.append(V)
    # -10- Interpolation type
    V = TK.StringVar(win); V.set('2nd order'); VARS.append(V)


    r = 0 # row
    # - inside/outside -
    B = TTK.Label(Frame,text='Blanking region')
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Blank region: inside or outside surfaces.')
    B = TTK.OptionMenu(Frame, VARS[4], 'inside', 'outside')
    B.grid(row=r, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Blank region: inside or outside surfaces.')
    r += 1
    # - Blanking type -
    B = TTK.Label(Frame,text='Blanking type')
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Blanking type.')
    B = TTK.OptionMenu(Frame, VARS[0], 'cell_intersect', 'cell_intersect_opt',
                       'center_in', 'node_in')
    B.grid(row=r, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Blanking type.')
    r += 1

    # - XRay delta  -
    B = TTK.Label(Frame, text="XRay delta")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='The created holes will expand of this value.')
    B = TTK.Entry(Frame, textvariable=VARS[2], background='White', width=5)
    B.grid(row=r, column=1, sticky=TK.EW)
    r+=1

    # - Xray tol -
    B = TTK.Label(Frame, text="Tol")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Two surface points distant from this value\nwill be considered as identical.')
    B = TTK.Entry(Frame, textvariable=VARS[3], background='White', width=5)
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1

    # - Surface -
    B = TTK.Button(Frame, text="Bodies", command=setSurface)
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Blanking bodies.')
    B = TTK.Entry(Frame, textvariable=VARS[5], background='White')
    B.grid(row=r, column=1, columnspan=3, sticky=TK.EW)
    r += 1

    # - depth -
    B = TTK.Label(Frame, text="Depth")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of IBC front layers.')
    #B = TK.Entry(Frame, textvariable=VARS[6], background='White', width=5)
    B = TTK.OptionMenu(Frame, VARS[6], '1', '2', '3', '4', '5')
    B.grid(row=r, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Number of IBC front layers.')
    r += 1

    # - blanking -
    B = TTK.Button(Frame, text="Blank by bodies", command=blank)
    B.grid(row=r, column=0, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Blank cells with with surface.')
    r += 1

    # - Donor zones : add ghost cells -
    B = TTK.Label(Frame, text="Add ghost cells")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Add ghost cells to donor zones.')
    B = TTK.OptionMenu(Frame, VARS[7], '0','1', '2', '3', '4', '5')
    B.grid(row=r, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Add ghost cells to donor zones.')
    r += 1

    # - Donor zones : location
    B = TTK.Label(Frame,text='Donors location')
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Donor zones location.')
    B = TTK.OptionMenu(Frame, VARS[8], 'nodes', 'centers','ext_centers')
    B.grid(row=r, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Donor zones location.')
    r += 1

    # - dhloc : modification de delta selon la formule ecrite en dur dans getIBCFrontInfo
    B = TTK.Label(Frame, text="distMin")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B,
                       text='Minimum distance to symmetrize inside and resulting fronts.')
    B = TTK.Entry(Frame, textvariable=VARS[9], background='White', width=5)
    B.grid(row=r, column=1, sticky=TK.EW)
    r += 1

    # - interpolation type :
    B = TTK.Label(Frame, text="Interpolation type")
    B.grid(row=r, column=0, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Interpolation type of IBC points.')
    B = TTK.OptionMenu(Frame, VARS[10], '2nd order','3rd order Lagrangian', '5th order Lagrangian')
    B.grid(row=r, column=1, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Interpolation type of IBC points.')
    r += 1


    # - create IBC front -
    B = TTK.Button(Frame, text="Create IBC front", command=getIBCFront)
    B.grid(row=r, column=0, columnspan=4, sticky=TK.EW)
    BB = CTK.infoBulle(parent=B, text='Create inside front.')
    r += 1

#==============================================================================
# Called to display widgets
#==============================================================================
def showApp():
    WIDGETS['frame'].grid(sticky=TK.NSEW); updateApp()

#==============================================================================
# Called to hide widgets
#==============================================================================
def hideApp(event=None):
    WIDGETS['frame'].grid_forget()

#==============================================================================
# Update widgets when global pyTree t changes
#==============================================================================
def updateApp(event=None): return

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
    (win, menu, file, tools) = CTK.minimal('tkIBC '+C.__version__)

    createApp(win); showApp()

    # - Main loop -
    win.mainloop()
