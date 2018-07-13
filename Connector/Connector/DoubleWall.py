# - Double Wall functions - 
import Connector
__version__ = Connector.__version__

try:
    import Geom
    import Transform
    import Transform.PyTree as T
    import Converter
    import Converter.PyTree as C
    import Converter.Internal as Internal
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
        if (v == 'BCWallViscous' or v == 'BCWall' or v == 'BCWallInviscid'):
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
# de range wallRanges de la zone z sous forme d une seule zone HEXA
# Coordonnees + indcellw + dir1 + dir2 + dir3
#------------------------------------------------------------------------------
def getFirstPointsInfo0__(z, wallRanges, loc='nodes'):
    if loc == 'nodes':  shift = 0
    elif loc == 'centers': 
        shift = 1 #decalage d indices de -1 si centres
        z = Converter.node2Center(z)
    else: raise ValueError("getFirstPointsInfo0__: loc must be nodes or centers.")

    ni = z[2]; nj = z[3]; nk = z[4]; ninj = ni*nj
    indwt = Converter.array('indcellw',ni,nj,nk); indca = indwt[1]
    for ind in xrange(indca.shape[1]): indca[0,ind] = float(ind)

    wallsc = []
    dirz = Converter.array('dir1,dir2,dir3',ni,nj,nk)
    # Premiere passe : on ajoute dans z2 les bonnes valeurs de dir1,dir2,dir3
    for w in wallRanges:
        i1 = w[0]; i2 = w[1]; j1 = w[2]; j2 = w[3]; k1 = w[4]; k2 = w[5]
        
        if i1 == i2:
            if i1 == 1: diri=1; ic = 0
            else: diri=-1; ic = i1-shift-1
            for k in range(k1-1,k2-shift):
                for j in range(j1-1,j2-shift):
                    ind = ic + j*ni + k*ninj
                    dirz[1][0,ind] = diri   
            
        elif j1 == j2:
            if j1 == 1: diri=1; jc = 0
            else: diri=-1; jc = j1-shift-1
            for k in range(k1-1,k2-shift):
                for i in range(i1-1,i2-shift):
                    ind = i + jc*ni + k*ninj
                    dirz[1][1,ind] = diri   
        
        else: # k1 == k2
            if k1 == 1: diri=1; kc = 0
            else: diri=-1; kc = k1-shift-1
            for j in range(j1-1,j2-shift):
                for i in range(i1-1,i2-shift):
                    ind = i + j*ni + kc*ninj
                    dirz[1][2,ind] = diri  

    z = Converter.addVars([z,indwt,dirz])
    # 2e passe : on subzone pour avoir les bons dir1,dir2,dir3
    for w in wallRanges:
        i1 = w[0]; i2 = w[1]; j1 = w[2]; j2 = w[3]; k1 = w[4]; k2 = w[5]
        if i1 == i2:
            if i1 == 1: diri=1; ic = 1
            else: diri=-1; ic = i1-shift            
            z2w = Transform.subzone(z,(ic,j1,k1),(ic,j2-shift,k2-shift))
            
        elif j1 == j2:
            if j1 == 1: diri=1; jc = 1
            else: diri=-1; jc = j1-shift
            z2w = Transform.subzone(z,(i1,jc,k1),(i2-shift,jc,k2-shift))

        else: # k1 == k2
            if k1 == 1: diri=1; kc = 1
            else: diri=-1; kc = k1-shift
            z2w = Transform.subzone(z,(i1,j1,kc),(i2-shift,j2-shift,kc))
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
def getFirstPointsInfo__(z, wallRanges, loc='nodes'):
    coords = C.getFields(Internal.__GridCoordinates__,z)[0]
    wallsc = getFirstPointsInfo0__(coords, wallRanges, loc)
    # calcul de la hauteur de courbure
    hmax = Geom.getCurvatureHeight(wallsc)
    wallsc = Converter.addVars([wallsc,hmax])
    # conversion en 'NODE'
    wallsc = Converter.convertArray2Node(wallsc)
    return wallsc

#----------------------------------------------------------------------------------------
# Met a jour les coordonnees des points frontieres d un maillage surfacique 
# Pour traiter les surfaces de projection en centres etendus issues de BCWall splittees
# ----------------------------------------------------------------------------------------
def modifyBorders__(a, iminL,imaxL,jminL,jmaxL):  
    return C.TZGC(a, 'nodes', Connector.modifyBorders__, iminL,imaxL,jminL,jmaxL)
# ----------------------------------------------------------------------------------------
# determination de toutes les sous fenetres definissant une frontiere i=imax par ex 
# ----------------------------------------------------------------------------------------
def getExtCenterSubWin(wloc,allIndices,dirW=0):
    wloc = C.node2ExtCenter(wloc)
    dimsW = Internal.getZoneDim(wloc)        
    niW=dimsW[1]; njW=dimsW[2];nkW=dimsW[3]
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
            iminL = jmin
            jminL = kmin
            imaxL = jmax+1
            jmaxL = kmax+1
        elif dirW == 2:
            iminL = imin
            jminL = kmin
            imaxL = imax+1
            jmaxL = kmax+1
        else: 
            iminL = imin
            jminL = jmin
            imaxL = imax+1
            jmaxL = jmax+1

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
        surfs+=getExtCenterSubWin(wloc,indicesI1,1)
    if indicesI2 != []:    
        wloc = T.subzone(a,(ni0-1,1,1),(ni0,nj0,nk0))
        surfs+=getExtCenterSubWin(wloc,indicesI2,1)
    if indicesJ1 != []:
        wloc = T.subzone(a,(1,1,1),(ni0,2,nk0))
        surfs+=getExtCenterSubWin(wloc,indicesJ1,2)
    if indicesJ2 != []:
        wloc = T.subzone(a,(1,nj0-1,1),(ni0,nj0,nk0)) 
        surfs+=getExtCenterSubWin(wloc,indicesJ2,2)
    if indicesK1 != []:
        wloc = T.subzone(a,(1,1,1),(ni0,nj0,2)) 
        surfs+=getExtCenterSubWin(wloc,indicesK1,3)
    if indicesK2 != []:
        wloc = T.subzone(a,(1,1,nk0-1),(ni0,nj0,nk0))         
        surfs+=getExtCenterSubWin(wloc,indicesK2,3)

    return surfs     
#-------------------------------------------------------------------------------
# Information sur les parois de a (pour le double wall) 
# utilise pour setInterpolations 
# Retourne: res sous forme d'arrays pour chaque base, pour chaque zone
# res[0]: indices,dir1,dir2,dir3,hmax des premiers points pres de la paroi
# res[1]: parois decrites par des maillages TRI avec point milieu.Contient hmax
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
    nob=0; wallSurfacesEC=[]; 
    for base1 in bases:
        zones1 = Internal.getNodesFromType1(base1,"Zone_t")
        surfacesExt=[]; surfacesPerBase=[]
        noz=0
        for z1 in zones1:
            surfs = getExtCenterSurfaces__(z1, wallBndIndicesN[nob][noz])
            if surfs != []:
                surfs = C.getAllFields(surfs,loc='nodes')
                #surfs = T.merge(surfs)
                surfs = Converter.convertArray2Hexa(surfs)
                surfs = Transform.join(surfs)
                surfacesPerBase.append(surfs)                                
            
            surfacesExt.append(surfs)
            noz+=1
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
            for nosz in xrange(len(surfacesExt)): # identification par zone
                surfECZ = surfacesExt[nosz]
                if surfECZ != []:
                    surfECZ = Converter.addVars(surfECZ,'hmax')
                    hsLoc = Converter.extractVars(surfECZ,['hmax'])
                    indices = Converter.identifyNodes(hook, surfECZ)
                    i = 0
                    for indh in indices: 
                        hsLoc[1][0,i]=hsRef[0,indh-1]
                        i+=1                
                    surfECZ = Converter.addVars([surfECZ,hsLoc])
                    surfECZ = Converter.convertArray2Tetra(surfECZ,split='withBarycenters')
                    surfacesExt[nosz]=[surfECZ]
            for nozOfNob in xrange(len(firstCenters[nob])):
                wallZ = firstCenters[nob][nozOfNob]
                if wallZ !=[]:
                    indices2 = Converter.identifyNodes(hook,wallZ)
                    hsLoc = Converter.addVars(wallZ,'hmax')
                    hsLoc = Converter.extractVars(hsLoc,['hmax'])
                    i = 0
                    for indh in indices2: 
                        hsLoc[1][0,i]=hsRef[0,indh-1]
                        i+=1                
                    firstCenters[nob][nozOfNob] = Converter.addVars([wallZ,hsLoc])
            Converter.freeHook(hook)
        wallSurfacesEC.append(surfacesExt)
        nob += 1
    return [firstCenters, wallSurfacesEC]
