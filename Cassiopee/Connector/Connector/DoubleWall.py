# - Double Wall functions -
from . import Connector
__version__ = Connector.__version__

try: range = xrange
except: pass

try:
    import Geom
    import Geom.PyTree as D
    import Transform
    import Transform.PyTree as T
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
    import Converter.Mpi as Cmpi
except:
    raise ImportError("Connector.DoubleWall requires Transform, Converter, Geom modules.")

#------------------------------------------------------------------------------
# Retourne les ranges (en noeuds) des frontieres paroi de la zone z
#------------------------------------------------------------------------------
def getBCWallRanges__(z, familyNames=[]):
    walls = []
    bnds = Internal.getNodesFromType2(z, 'BC_t')
    for bc in bnds:
        v = Internal.getValue(bc)
        if v == 'BCWallViscous' or v == 'BCWall' or v == 'BCWallInviscid':
            range0 = Internal.getNodeFromName1(bc, 'PointRange')
            r = Internal.range2Window(range0[1])
            i1 = int(r[0]); j1 = int(r[2]); k1 = int(r[4])
            i2 = int(r[1]); j2 = int(r[3]); k2 = int(r[5])
            walls.append([i1,i2,j1,j2,k1,k2])
        elif v == 'FamilySpecified':
            for fname in familyNames:
                nodes = C.getFamilyBCs(z, fname)
                for fambc in nodes:
                    range0 = Internal.getNodeFromName1(fambc, 'PointRange')
                    r = Internal.range2Window(range0[1])
                    i1 = int(r[0]); j1 = int(r[2]); k1 = int(r[4])
                    i2 = int(r[1]); j2 = int(r[3]); k2 = int(r[5])
                    walls.append([i1,i2,j1,j2,k1,k2])
    return walls

#===============================================================================
# CI-DESSOUS: SPECIFIQUE DONNEURS CENTRES ETENDUS
# validite : setInterpolations, optimizeOverlap
#===============================================================================
#------------------------------------------------------------------------------
# Retourne les premiers pts (centres ou noeuds selon loc) pres des fenetres
# de range wallRanges de la zone z sous forme d'une seule zone HEXA
# Coordonnees + indcellw + dir1 + dir2 + dir3
#------------------------------------------------------------------------------
def getFirstPointsInfo0__(z, wallRanges, loc='nodes', ghostCells=False):
    if loc == 'nodes': shift = 0
    elif loc == 'centers':
        shift = 1 #decalage d'indices de -1 si centres
        z = Converter.node2Center(z)
    else: raise ValueError("getFirstPointsInfo0__: loc must be nodes or centers.")

    # CB: ghost cells in t
    if ghostCells: s2 = 2
    else: s2 = 0

    ni = z[2]; nj = z[3]; nk = z[4]; ninj = ni*nj
    indwt = Converter.array('indcellw',ni,nj,nk); indca = indwt[1]
    for ind in range(indca.shape[1]): indca[0,ind] = float(ind)

    wallsc = []
    dirz = Converter.array('dir1,dir2,dir3',ni,nj,nk)
    # Premiere passe : on ajoute dans z2 les bonnes valeurs de dir1,dir2,dir3
    for w in wallRanges:
        i1 = w[0]; i2 = w[1]; j1 = w[2]; j2 = w[3]; k1 = w[4]; k2 = w[5]

        if i1 == i2:
            if i1 == 1: diri=1; ic = 0
            else: diri=-1; ic = i1-shift-1
            for k in range(k1-1-s2,k2-shift+s2):
                for j in range(j1-1-s2,j2-shift+s2):
                    ind = ic + j*ni + k*ninj
                    dirz[1][0,ind] = diri

        elif j1 == j2:
            if j1 == 1: diri=1; jc = 0
            else: diri=-1; jc = j1-shift-1
            for k in range(k1-1-s2,k2-shift+s2):
                for i in range(i1-1-s2,i2-shift+s2):
                    ind = i + jc*ni + k*ninj
                    dirz[1][1,ind] = diri

        else: # k1 == k2
            if k1 == 1: diri=1; kc = 0
            else: diri=-1; kc = k1-shift-1
            for j in range(j1-1-s2,j2-shift+s2):
                for i in range(i1-1-s2,i2-shift+s2):
                    ind = i + j*ni + kc*ninj
                    dirz[1][2,ind] = diri

    z = Converter.addVars([z,indwt,dirz])
    # 2e passe : on subzone pour avoir les bons dir1,dir2,dir3
    for w in wallRanges:
        i1 = w[0]; i2 = w[1]; j1 = w[2]; j2 = w[3]; k1 = w[4]; k2 = w[5]
        if i1 == i2:
            if i1 == 1: diri=1; ic = 1
            else: diri=-1; ic = i1-shift
            z2w = Transform.subzone(z,(ic,j1-s2,k1-s2),(ic,j2-shift+s2,k2-shift+s2))

        elif j1 == j2:
            if j1 == 1: diri=1; jc = 1
            else: diri=-1; jc = j1-shift
            #print('range', i1-s2, i2-shift+s2, k1-s2, k2-shift+s2, flush=True)
            z2w = Transform.subzone(z,(i1-s2,jc,k1-s2),(i2-shift+s2,jc,k2-shift+s2))

        else: # k1 == k2
            if k1 == 1: diri=1; kc = 1
            else: diri=-1; kc = k1-shift
            z2w = Transform.subzone(z,(i1-s2,j1-s2,kc),(i2-shift+s2,j2-shift+s2,kc))
        wallsc.append(z2w)

    # conversion en 'QUAD'
    # wallsc = T.merge(wallsc)# NE SERT A RIEN
    wallsc = Converter.convertArray2Hexa(wallsc)
    wallsc = Transform.join(wallsc)
    return wallsc

#------------------------------------------------------------------------------
# Retourne les premiers pts (centres ou noeuds selon loc) pres des fenetres
# de range wallRanges de la zone z
# Coordonnees + dir1 + dir2 + dir3 + curvature height
#------------------------------------------------------------------------------
def getFirstPointsInfo__(z, wallRanges, loc='nodes', ghostCells=False):
    coords = C.getFields(Internal.__GridCoordinates__, z)[0]
    wallsc = getFirstPointsInfo0__(coords, wallRanges, loc, ghostCells)
    # calcul de la hauteur de courbure
    hmax = Geom.getCurvatureHeight(wallsc)
    wallsc = Converter.addVars([wallsc, hmax])
    # conversion en 'NODE'
    wallsc = Converter.convertArray2Node(wallsc)
    return wallsc

#----------------------------------------------------------------------------------------
# Met a jour les coordonnees des points frontieres d'un maillage surfacique
# Pour traiter les surfaces de projection en centres etendus issues de BCWall splittees
# ----------------------------------------------------------------------------------------
def modifyBorders__(a, iminL, imaxL, jminL, jmaxL):
    return C.TZGC1(a, 'nodes', True, Connector.modifyBorders__, iminL, imaxL, jminL, jmaxL)

# ----------------------------------------------------------------------------------------
# determination de toutes les sous fenetres definissant une frontiere i=imax par ex
# ----------------------------------------------------------------------------------------
def getExtCenterSubWin(wloc, allIndices, dirW=0):
    wloc = C.node2ExtCenter(wloc)
    dimsW = Internal.getZoneDim(wloc)
    niW=dimsW[1]; njW=dimsW[2]; nkW=dimsW[3]
    if   dirW == 1: wLoc = T.subzone(wloc,(2,1,1),(2,njW,nkW))
    elif dirW == 2: wLoc = T.subzone(wloc,(1,2,1),(niW,2,nkW))
    elif dirW == 3: wLoc = T.subzone(wloc,(1,1,2),(niW,njW,2))
    surfs=[]

    for indices in allIndices:
        imin = indices[0]; imax = indices[3]
        jmin = indices[1]; jmax = indices[4]
        kmin = indices[2]; kmax = indices[5]
        dimS = Internal.getZoneDim(wLoc)
        niS = dimS[1]; njS = dimS[2]; nkS = dimS[3]

        if dirW == 1:
            iminL = jmin ; jminL = kmin
            imaxL = jmax+1; jmaxL = kmax+1
        elif dirW == 2:
            iminL = imin; jminL = kmin
            imaxL = imax+1; jmaxL = kmax+1
        else:
            iminL = imin; jminL = jmin
            imaxL = imax+1; jmaxL = jmax+1

        sLoc = modifyBorders__(wLoc,iminL,imaxL,jminL,jmaxL)
        if dirW == 1: sLoc = T.subzone(sLoc,(jmin,kmin,1),(jmax+1,kmax+1,1))
        elif dirW == 2: sLoc = T.subzone(sLoc,(imin,kmin,1),(imax+1,kmax+1,1))
        elif dirW == 3: sLoc = T.subzone(sLoc,(imin,jmin,1),(imax+1,jmax+1,1))
        surfs.append(sLoc)

    return surfs

#-------------------------------------------------------------------------------
def getExtCenterSurfaces__(a, walls):
    dims0 = Internal.getZoneDim(a)
    ni0 = dims0[1]; nj0 = dims0[2]; nk0 = dims0[3]
    indicesI1=[]; indicesI2=[];indicesJ1=[];indicesJ2=[];indicesK1=[];indicesK2=[]
    for w in walls:
        i1 = int(w[0]); i2=int(w[1]); j1 = int(w[2]); j2 = int(w[3]); k1 = int(w[4]); k2 = int(w[5])
        if i1 == i2:
            if i1==1: indicesI1.append([i1,j1,k1,i2,j2,k2])
            else: indicesI2.append([i1,j1,k1,i2,j2,k2])
        elif j1 == j2:
            if j1 == 1: indicesJ1.append([i1,j1,k1,i2,j2,k2])
            else: indicesJ2.append([i1,j1,k1,i2,j2,k2])
        elif k1 == k2:
            if k1 == 1: indicesK1.append([i1,j1,k1,i2,j2,k2])
            else: indicesK2.append([i1,j1,k1,i2,j2,k2])

    surfs = []
    if indicesI1 != []:
        wloc = T.subzone(a,(1,1,1),(2,nj0,nk0))
        surfs += getExtCenterSubWin(wloc,indicesI1,1)
    if indicesI2 != []:
        wloc = T.subzone(a,(ni0-1,1,1),(ni0,nj0,nk0))
        surfs += getExtCenterSubWin(wloc,indicesI2,1)
    if indicesJ1 != []:
        wloc = T.subzone(a,(1,1,1),(ni0,2,nk0))
        surfs += getExtCenterSubWin(wloc,indicesJ1,2)
    if indicesJ2 != []:
        wloc = T.subzone(a,(1,nj0-1,1),(ni0,nj0,nk0))
        surfs += getExtCenterSubWin(wloc,indicesJ2,2)
    if indicesK1 != []:
        wloc = T.subzone(a,(1,1,1),(ni0,nj0,2))
        surfs += getExtCenterSubWin(wloc,indicesK1,3)
    if indicesK2 != []:
        wloc = T.subzone(a,(1,1,nk0-1),(ni0,nj0,nk0))
        surfs += getExtCenterSubWin(wloc,indicesK2,3)

    return surfs
#-------------------------------------------------------------------------------
# Information sur les parois de a (pour le double wall)
# utilise pour setInterpolations
# Retourne: res sous forme d'arrays pour chaque base, pour chaque zone
# res[0]: indices,dir1,dir2,dir3,hmax des premiers points pres de la paroi
# res[1]: parois decrites par des maillages TRI avec point milieu. Contient hmax.
#--------------------------------------------------------------------------------
def extractDoubleWallInfo__(t):
    firstCenters = [] # indices des premiers centres pres des parois, leur direction et leur hauteur de courbure
    wallSurfacesEC = [] # surfaces en centres etendus, TRI
    bases = Internal.getNodesFromType1(t, 'CGNSBase_t')
    wallBndIndicesN = [] # indices des parois
    # 1. Extraction des surfaces paroi et indices correspondants
    # liste des familles "paroi"
    famwalls = C.getFamilyBCNamesOfType(t, 'BCWall')
    famwalls += C.getFamilyBCNamesOfType(t, 'BCWallViscous')
    famwalls += C.getFamilyBCNamesOfType(t, 'BCWallInviscid')
    for base1 in bases:
        firstCentersb = [] # liste des surfaces de la base
        wallsb = []
        zones1 = Internal.getNodesFromType1(base1, 'Zone_t')
        for z1 in zones1:
            walls = getBCWallRanges__(z1, famwalls)
            if walls == []: firstCenters1 = []
            else:
                # firstCenters : CoordinateX,CoordinateY,CoordinateZ,indcellw,dir1,dir2,dir3,hmax
                firstCenters1 = getFirstPointsInfo__(z1, walls, loc='centers')
            wallsb.append(walls); firstCentersb.append(firstCenters1)
        # par base : on recupere les infos
        firstCenters.append(firstCentersb); wallBndIndicesN.append(wallsb)
    # fin parcours de toutes les bases
    # 2. on passe en centres etendus les parois de projection
    # Il faut assurer la continuite de hmax au passage des raccords : passage en tri de la base
    #
    nob=0; wallSurfacesEC=[]
    for base1 in bases:
        zones1 = Internal.getNodesFromType1(base1, "Zone_t")
        surfacesExt=[]; surfacesPerBase=[]
        noz=0
        for z1 in zones1:
            surfs = getExtCenterSurfaces__(z1, wallBndIndicesN[nob][noz])
            if surfs != []:
                surfs = C.getAllFields(surfs, loc='nodes')
                #surfs = T.merge(surfs)
                surfs = Converter.convertArray2Hexa(surfs)
                surfs = Transform.join(surfs)
                surfacesPerBase.append(surfs)

            surfacesExt.append(surfs)
            noz += 1
        # Join des surfaces de la meme base pour assurer la continuite de hmax
        surfacesPerBase = Transform.join(surfacesPerBase)
        # Calcul de la hauteur de courbure
        if surfacesPerBase != []:
            hsRef = Geom.getCurvatureHeight(surfacesPerBase)
            hsRef = Converter.node2Center(hsRef)
            hsRef = Converter.center2Node(hsRef) # lissage en prenant le premier voisinage
            #
            # Identify nodes
            hsRef = hsRef[1]
            hook = Converter.createHook(surfacesPerBase,'nodes')
            for nosz in range(len(surfacesExt)): # identification par zone
                surfECZ = surfacesExt[nosz]
                if surfECZ != []:
                    surfECZ = Converter.addVars(surfECZ,'hmax')
                    hsLoc = Converter.extractVars(surfECZ,['hmax'])
                    indices = Converter.identifyNodes(hook, surfECZ)
                    i = 0
                    for indh in indices:
                        hsLoc[1][0,i] = hsRef[0,indh-1]
                        i += 1
                    surfECZ = Converter.addVars([surfECZ,hsLoc])
                    surfECZ = Converter.convertArray2Tetra(surfECZ, split='withBarycenters')
                    surfacesExt[nosz]=[surfECZ]
            for nozOfNob in range(len(firstCenters[nob])):
                wallZ = firstCenters[nob][nozOfNob]
                if wallZ != []:
                    indices2 = Converter.identifyNodes(hook,wallZ)
                    hsLoc = Converter.addVars(wallZ,'hmax')
                    hsLoc = Converter.extractVars(hsLoc,['hmax'])
                    i = 0
                    for indh in indices2:
                        hsLoc[1][0,i] = hsRef[0,indh-1]
                        i += 1
                    firstCenters[nob][nozOfNob] = Converter.addVars([wallZ,hsLoc])
            Converter.freeHook(hook)
        wallSurfacesEC.append(surfacesExt)
        nob += 1
    return [firstCenters, wallSurfacesEC]

#=========================================================================
# Change wall with explicit wall windows for Fast t,tc
# listOfMismatch : liste des BCs en mismatch
# listOfMismatch1 = ['Base1/cart1/wall1']
# Prend les fenetres de mismatch2, les projetent sur les fenetres de
# mismatch1 et modifie tc
#=========================================================================
def _changeWall2(t, tc, listOfMismatch1, listOfMismatch2, familyBC1, familyBC2, ghostCells=False, check=False, surfaceCenters1=None):

    # STEP0 : gather lists
    listOfMismatch1 = Cmpi.allreduce(listOfMismatch1)
    listOfMismatch2 = Cmpi.allreduce(listOfMismatch2)

    if Cmpi.rank==0:
        print("Info: changeWall: listOfMismatch1 (%s)"%(listOfMismatch1))
        print("Info: changeWall: listOfMismatch2 (%s)"%(listOfMismatch2))

    # STEP1 : gather proj surfaces of mismatch1 (surfaceCenters1)
    if surfaceCenters1 is None:
        walls1 = []
        for w1 in listOfMismatch1:
            name = w1.rsplit('/', 1)
            z1 = Internal.getNodeFromPath(t, name[0])
            if z1 is not None:
                walls = C.extractBCOfType(z1, 'BCWall')
                walls1 += walls
                walls = C.extractBCOfType(z1, 'BCWallExchange')
                walls1 += walls
                walls = C.extractBCOfType(z1, 'BCWallmodel')
                walls1 += walls

        #C.convertPyTree2File(walls1, 'walls%d.cgns'%Cmpi.rank)
        walls1 = Cmpi.allgatherZones(walls1, coord=True, variables=['centers:cellN'])
        #if Cmpi.rank == 0: C.convertPyTree2File(walls1, 'walls.cgns')

        if walls1 == []: return None

        walls1 = C.convertArray2Hexa(walls1)
        walls1 = T.join(walls1)

        walls1 = C.node2Center(walls1)
        walls1 = C.convertArray2Tetra(walls1, split='withBarycenters')
        #walls1 = T.reorderAll(walls1)
        walls1 = T.reorder(walls1, (1,))

        D._getCurvatureHeight(walls1)
        surfaceCenters1 = C.getAllFields(walls1, 'nodes') # array with __GridCoordinates__ + __FlowSolutionNodes_
        if check and Cmpi.rank==0: Converter.convertArrays2File(surfaceCenters1, 'surfaceCenters1_%s_to_%s.plt'%(familyBC2, familyBC1))

    # STEP2 : project surfaces of mismatch2 (a2c (domain) and firstWallCenters2 (wall))
    for w2 in listOfMismatch2:
        name = w2.rsplit('/', 1)
        z2 = Internal.getNodeFromPath(t, name[0])
        if z2 is not None:
            pr = Internal.getNodeFromPath(z2, 'ZoneBC/'+name[1]+'/PointRange')
            if pr is not None:
                wr = Internal.range2Window(pr[1])
                firstWallCenters2 = getFirstPointsInfo__(z2, [wr], loc='centers', ghostCells=ghostCells)
                if check: Converter.convertArrays2File(firstWallCenters2, 'firstWallCenters2_%s_to_%s_rank%d.plt'%(familyBC2, familyBC1, Cmpi.rank))

                z2c = Internal.getNodeFromPath(tc, name[0])
                if z2c is not None:
                    a2c = C.getFields(Internal.__GridCoordinates__, z2c)[0]
                    cellN = C.getField('cellN', z2c)[0]
                    a2c = Converter.addVars([a2c,cellN]) # array at centers with cellN
                    if check: Converter.convertArrays2File(a2c, 'surfaceCenters2_%s_to_%s_rank%d.plt'%(familyBC2, familyBC1, Cmpi.rank))

                    a2cp = Connector.changeWall__(a2c, firstWallCenters2, surfaceCenters1, planarTol=0.)

                    C.setFields([a2cp], z2c, 'nodes', False) # modify tc
                    if check: Converter.convertArrays2File(a2cp, 'surfaceCenters2p_%s_to_%s_rank%d.plt'%(familyBC2, familyBC1, Cmpi.rank))

    return None

#=========================================================================
# Generate the projection surface for _changeWall2
#=========================================================================
def getProjSurfaceForDoubleWall(t, listOfMismatch1, familyBC1, check=False):

    listOfMismatch1 = Cmpi.allreduce(listOfMismatch1)

    walls1 = []
    for w1 in listOfMismatch1:
        name = w1.rsplit('/', 1)
        z1 = Internal.getNodeFromPath(t, name[0])
        if z1 is not None:
            walls = C.extractBCOfType(z1, 'BCWall')
            walls1 += walls

    walls1 = Cmpi.allgatherZones(walls1, coord=True, variables=['centers:cellN'])

    if walls1 == []: return None

    walls1 = C.convertArray2Hexa(walls1)
    walls1 = T.join(walls1)

    walls1 = C.node2Center(walls1)
    walls1 = C.convertArray2Tetra(walls1, split='withBarycenters')
    walls1 = T.reorder(walls1, (1,))

    D._getCurvatureHeight(walls1)
    surfaceCenters1 = C.getAllFields(walls1, 'nodes') # array with __GridCoordinates__ + __FlowSolutionNodes_
    if check and Cmpi.rank==0: Converter.convertArrays2File(surfaceCenters1, 'surfaceCenters1_%s.plt'%(familyBC1))

    return surfaceCenters1
