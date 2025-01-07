# - OversetDataDtlocal -
# Module for internal functions used for overset info and for local time stepping
from . import Connector
from . import connector
from . import OversetData
import numpy
__version__ = Connector.__version__

try:
    import Converter.Internal as Internal
    import Converter.PyTree as C
    import Converter
    import KCore
except:
    raise ImportError("Connector.OversetDataDtlocal requires Converter module.")

#============================================================================================================
#============================================================================================================
#============================================================================================================

def setInterpData3(tR, tD, double_wall=0, order=2, penalty=1, nature=0,
                   method='lagrangian', loc='nodes', storage='direct',
                   hook=None,
                   topTreeRcv=None, topTreeDnr=None, sameName=1, dim=3, itype='both'):

    locR = loc
    aR = Internal.copyRef(tR)
    aD = Internal.copyRef(tD)

    bases  = Internal.getNodesFromType1(aR     , 'CGNSBase_t')
    dimmm  = Internal.getNodeFromName2(bases[0], 'EquationDimension')
    dimPb   = Internal.getValue(dimmm)

    zones = Internal.getZones(aR)
    niveaux_temps={}
    for z in zones:
        nodes = Internal.getNodesFromName(z, 'niveaux_temps')[0]
        local_time_level = int(nodes[1][0][0][0])
        niveaux_temps[z[0]] = local_time_level

    """
    dxmax = 0.0 

    zones = Internal.getZones(aR)
    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]
  
        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0])
        
        dx = min(dxx,dyy,dzz)
        if (dx > dxmax):dxmax=dx

    niveaux_temps = {} 
    for z in zones:
        nodes = Internal.getNodesFromName(z, 'GridCoordinates')
        coordx = nodes[0][2][0][1]
        coordy = nodes[0][2][1][1]
        coordz = nodes[0][2][2][1]

        dxx  = abs(coordx[1,0,0]   - coordx[0,0,0])
        dyy  = abs(coordy[0,1,0]   - coordy[0,0,0])
        dzz  = abs(coordz[0,0,1]   - coordz[0,0,0]) 
        dx = min(dxx,dyy,dzz)
        niveaux_temps[z[0]] = round(dxmax/dx)

        print round(dxmax/dx)

    """

    # Recherche pour les pts coincidents (base sur les GridConnectivity)
    if itype != 'chimera':
        if storage == 'direct': aR = setInterpDataForGhostCells2__(aR,aD,storage,loc)
        else:
            aD = setInterpDataForGhostCells2__(aR,aD,storage,loc)

     # Si pas de cellN receveur, on retourne
    if loc == 'nodes': cellNPresent = C.isNamePresent(aR, 'cellN')
    elif loc == 'centers': cellNPresent = C.isNamePresent(aR, 'centers:cellN')
    else:
        raise ValueError("setInterpData: invalid loc provided.")
    if cellNPresent == -1 or itype == 'abutting':
        if storage == 'direct': return aR
        else: return aD

    locCellND = 'nodes'
    aD = OversetData.addCellN__(aD, loc=locCellND)

    if method == 'conservative' and itype != 'abutting':
        if loc != 'centers':
            raise ValueError("setInterpData: conservative type is available only for loc='centers'.")
        else: return setInterpDataConservative__(aR, aD, storage=storage)

    zonesRcv = Internal.getZones(aR); nzonesRcv = len(zonesRcv)
    zonesDnr = Internal.getZones(aD); nzonesDnr = len(zonesDnr)

    #print  'zonesRcv = ' , zonesRcv[1]
    #print  'zonesDnr = ' , zonesDnr[1]
    #---------------------------------------------------------------------------
    # CAS DOUBLE WALL
    # Extraction des parois de projection issues des zones donneuses
    # Extraction des premiers centres ou premiers noeuds (selon locR) des zones receveuses
    #---------------------------------------------------------------------------
    donorSurfs = []; interpWallPts = []
    noWallsInDnr = 1 # parametre qui vaut 1 si pas de parois trouvees dans les zones donneuses, 0 si au moins une
    noWallsInRcv = 1 # idem mais pour les receveurs
    if double_wall == 1:
        import DoubleWall
        try: import Geom.PyTree as D
        except: raise ImportError("setInterpData+double wall requires Geom.PyTree module.")

        ttreeR = []
        if topTreeRcv is not None: ttreeR = topTreeRcv
        else:
            if Internal.isTopTree(aR): ttreeR = aR
            else:
                if Internal.getBases(aR) != []: ttreeR = aR# c est une base on peut recuperer les familles
                else: print('Warning: setInterpData+double wall: receptor zones may require a top tree.')
        ttreeD = []
        if topTreeDnr is not None: ttreeD = topTreeDnr
        else:
            if Internal.isTopTree(aD): ttreeD = aD
            else:
                if Internal.getBases(aD) != []: ttreeD = aD # c'est une base on peut recuperer les familles
                else: print('Warning: setInterpData+double wall: donors zones may require a top tree.')

        # Zones donneuses : on recupere les surfaces de projection double wall
        for zd in zonesDnr:
            walls = C.extractBCOfType(zd,'BCWall',topTree=ttreeD)
            if walls != []:
                noWallsInDnr = 0
                walls = D.getCurvatureHeight(walls)
                walls = C.convertArray2Tetra(walls, split="withBarycenters")
                walls = C.getAllFields(walls,loc='nodes')
            donorSurfs.append(walls)

        # Zones receveuses : on determine les 1ers points paroi (centres ou noeuds)
        # recup des familles de type paroi dans les zones receveuses
        famwallsR = C.getFamilyBCNamesOfType(ttreeR, 'BCWall')
        famwallsR += C.getFamilyBCNamesOfType(ttreeR, 'BCWallViscous')
        famwallsR += C.getFamilyBCNamesOfType(ttreeR, 'BCWallInviscid')

        # interpWallPts : selon la loc, on recupere les premiers centres ou noeuds
        for zr in zonesRcv:
            wallRanges = DoubleWall.getBCWallRanges__(zr,famwallsR)
            if wallRanges != []:
                noWallsInRcv = 0
                if locR == 'nodes': interpWallPts.append(DoubleWall.getFirstPointsInfo__(zr, wallRanges,loc='nodes'))
                else: interpWallPts.append(DoubleWall.getFirstPointsInfo__(zr, wallRanges,loc='centers'))
            else: interpWallPts.append([])

    if noWallsInDnr == 1 or noWallsInRcv == 1: double_wall = 0 # on desactive le double wall
    arraysD = C.getFields(Internal.__GridCoordinates__, zonesDnr)
    cellND = C.getField('cellN', zonesDnr)
    arraysD = Converter.addVars([arraysD,cellND])
    cellND = []
    #---------------------------------------------------------------------------
    # 1. Extraction des points a interpoler
    #    interpPts : un array si pas de double wall
    #              : une liste d arrays avec les coordonnees des pts a interpoler modifiees
    # 2. Calcul des coefs et donneurs
    #---------------------------------------------------------------------------
    if hook is not None:
        allHooks = hook[:]
        allHooksL = allHooks[:]
    else: allHooks = None

    arraysDL = arraysD[:]
    donorSurfsL = donorSurfs[:]
    zonesDnrL = zonesDnr[:]
    nzonesDnrL = nzonesDnr
    nozr = -1

    #dico={}
    #for z in zonesRcv:
    # listofjoins = Internal.getNodesFromType2(z, 'GridConnectivity_t')
    # if listofjoins is not None:
    #    prange_list=[]
    #    for join in listofjoins:
    #        prange = Internal.getNodeFromName1(join,'PointRange')[1]
    #        for i in range(3):
    #            if prange[i,1] == prange[i,0] and prange[i,1] != 1:
    #                  prange[i,1] =  prange[i,1]-1
    #                  prange[i,0] =  prange[i,0]-1
    #            elif prange[i,1] != prange[i,0] and prange[i,1] != 1 :
    #                  prange[i,1] =  prange[i,1]-1
    #        prange=numpy.reshape(prange,6)
    #        #print prange
    #        prange_list.append(prange)
    #    dico[z[0]] = prange_list
    #print dico.get('cart.0')

    for z in zonesRcv:
        dim_ = Internal.getZoneDim(z)
        #if nodes is not None:
        #nodesname = nodes[0][0][0:6]
        #if nodesname == 'nmatch':
        #print 'ok'
        #print nodes[0]


        nozr += 1
        if loc == 'nodes': cellNPresent = C.isNamePresent(z, 'cellN')
        else: cellNPresent = C.isNamePresent(z, 'centers:cellN')
        if cellNPresent != -1:
            if sameName == 1:
                arraysD = arraysDL[:]
                if hook is not None: allHooks = allHooksL[:]
                donorSurfs = donorSurfsL[:]
                zonesDnr = zonesDnrL[:]
                nzonesDnr = nzonesDnrL
                cL = 0; found = 0
                for zo in zonesDnr:
                    if zo[0] == z[0]: found = 1; break;
                    cL += 1
                if found == 1:
                    arraysD.pop(cL)
                    if hook is not None: allHooks.pop(cL)
                    zonesDnr.pop(cL); nzonesDnr = nzonesDnr-1
                    if double_wall == 1: donorSurfs.pop(cL)

            #-------------------------------------------
            # Etape 1: on recupere les pts a interpoler
            #-------------------------------------------
            interpPts = [];
            isdw = 0 # si double wall effectivement active, interpPts est une liste d'arrays
            if locR == 'nodes':
                an = C.getFields(Internal.__GridCoordinates__, z)[0]
                cellN = C.getField('cellN',z)[0]
                an = Converter.addVars([an,cellN])
                if double_wall == 1 and interpWallPts[nozr] != []: # dw: liste d arrays
                    isdw = 1
                    for nozd in range(nzonesDnr):
                        an2 = Connector.changeWall__(an, interpWallPts[nozr], donorSurfs[nozd], planarTol=0.)
                        interpPts.append(Connector.getInterpolatedPoints__(an2))
                else: # pas de dw: un seul array
                    interpPts = Connector.getInterpolatedPoints__(an)

            elif locR == 'centers':
                an = C.getFields(Internal.__GridCoordinates__, z)[0]
                ac = Converter.node2Center(an)
                cellN = C.getField('centers:cellN',z)[0]
                ac = Converter.addVars([ac,cellN])
                if double_wall == 1 and interpWallPts[nozr] != []:# dw : liste d arrays
                    isdw = 1
                    for nozd in range(nzonesDnr):
                        ac2 = Connector.changeWall__(ac, interpWallPts[nozr], donorSurfs[nozd], planarTol=0.)
                        interpPts.append(Connector.getInterpolatedPoints__(ac2))
                else:  # pas de dw : un seul array
                    interpPts = Connector.getInterpolatedPoints__(ac)
                    #print interpPts[1].shape
            #---------------------------------------------
            # Etape 2 : calcul des donnees d'interpolation
            #---------------------------------------------
            # resInterp = [rcvInd1D,donorInd1D,donorType,coefs,extrap,orphan, EXdirs]
            resInterp = Connector.setInterpData__(interpPts, arraysD, order=order, penalty=penalty, nature=nature, method=method, hook=allHooks, dim=dim)
            if resInterp is not None:
                # Bilan
                nborphan = 0; nbextrapolated = 0; nbinterpolated = 0
                if double_wall==0: nbinterpolated = interpPts[1].shape[1]
                else: nbinterpolated = interpPts[0][1].shape[1]
                if len(resInterp[4])>0:
                    for r in resInterp[4]: nbextrapolated += r.shape[0]
                nborphan = resInterp[5].size
                nbinterpolated = nbinterpolated-nbextrapolated-nborphan

                print('Zone %s: interpolated=%d ; extrapolated=%d ; orphan=%d'%(z[0], nbinterpolated, nbextrapolated, nborphan))
                if  nborphan>0: print('Warning: zone %s has %d orphan points !'%(z[0], nborphan))
                # on remet a une seule zone, attention si x,y,z sont necessaires ensuite
                # les coordonnees peuvent etre fausses a cause du double walls
                indcells=[]
                if loc == 'faces': vari = 'indcell1'
                else: vari = 'indcell'
                if isdw == 1:
                    posindcell = KCore.isNamePresent(interpPts[0],vari)
                    if posindcell != -1: indcells = interpPts[0][1][posindcell,:]
                else:
                    posindcell = KCore.isNamePresent(interpPts,vari)
                    if posindcell != -1: indcells = interpPts[1][posindcell,:]

                # on recupere les bons indices de pts interpoles (indcell ds interpPts), de EXdir
                # pour ceux obtenus en sortie de setInterpData
                # Orphelins
                if nborphan > 0:
                    for noi in range(nborphan):
                        noind = resInterp[5][noi]
                        resInterp[5][noi] = indcells[noind]
                # Interpoles/Extrapoles
                for noz in range(nzonesDnr):
                    ninterploc = resInterp[0][noz].size
                    if ninterploc > 0:# domaine d'interpolation
                        # Indices des receveurs
                        for noi in range(ninterploc):
                            index = resInterp[0][noz][noi]
                            resInterp[0][noz][noi] = indcells[index]
                if len(resInterp[4])>0: # pts extrapoles par zone donneuse
                    for noz in range(nzonesDnr):
                        nextraploc = resInterp[4][noz].size
                        if nextraploc > 0:# Extrapoles
                            for noi in range(nextraploc):
                                index = resInterp[4][noz][noi]
                                resInterp[4][noz][noi] = indcells[index]

                #----------------------------------
                # Etape 3: Stockage dans l'arbre
                # direct: on stocke dans aR
                # inverse: on stocke dans aD
                #----------------------------------

                levelrcv = niveaux_temps[z[0]]

                for noz in range(nzonesDnr):
                    # ajout des donnees d interpolation
                    ninterploc = resInterp[0][noz].size
                    if ninterploc>0: # domaine d'interpolation
                        if storage == 'direct':
                            if resInterp[4][noz].size == 0: extrapPts = numpy.array([], Internal.E_NpyInt)
                            else: extrapPts = resInterp[4][noz]
                            if resInterp[6] == []: EXdir = numpy.array([], Internal.E_NpyInt)
                            else: EXdir = resInterp[6][noz]
                            OversetData._createInterpRegion__(z, zonesDnr[noz][0], resInterp[0][noz], resInterp[1][noz], resInterp[3][noz], \
                                                              resInterp[2][noz], numpy.array([],numpy.float64), \
                                                              extrapPts, resInterp[5], tag='Receiver', loc=locR, EXDir=EXdir)

                        else: # inverse
                            if resInterp[4][noz].size == 0: extrapPts = numpy.array([], Internal.E_NpyInt)
                            else: extrapPts = resInterp[4][noz]
                            if resInterp[6] == []: EXdir = numpy.array([], Internal.E_NpyInt)
                            else: EXdir = resInterp[6][noz]
                            OversetData._createInterpRegion__(zonesDnr[noz], z[0], resInterp[1][noz], resInterp[0][noz], resInterp[3][noz], \
                                                              resInterp[2][noz], numpy.array([],numpy.float64), \
                                                              extrapPts, resInterp[5], tag='Donor', loc=locR, EXDir=EXdir)



                            dim__ = Internal.getZoneDim(zonesDnr[noz])
                            prange = numpy.zeros(6,dtype=Internal.E_NpyInt)
                            dirR=numpy.zeros(1,dtype=Internal.E_NpyInt)
                            leveldnr = niveaux_temps[zonesDnr[noz][0]]

                            #print('donneur= ', zonesDnr[noz][0],noz)
                            #print 'donneur= ', zonesRcv[noz][0],noz
                            #print('receveur= ', z[0])
                            #connector.indiceToCoordbis(resInterp[0][noz],prange,dirR,resInterp[0][noz].size,resInterp[2][noz][0],dim_[1],dim_[2],dim_[3])

                            prange=numpy.reshape(prange,(3,2))
                            #print prange
                            #dirR = GhostCells.getDirBorderStruct__(prange,dimPb)
                            #print dirR
                            prangebis=numpy.reshape(prange,6)
                            #print(prange)
                            info = zonesDnr[noz][2][len(zonesDnr[noz][2])-1]
                            info[2].append(['PointRange', prangebis , [], 'IndexArray_t'])

                            transfo=getTransfo(zonesDnr[noz],z)


                            prangedonor = numpy.zeros(6,dtype=Internal.E_NpyInt)
                            profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                            dirD=numpy.zeros(1,dtype=Internal.E_NpyInt)
                            connector.indiceToCoord2(resInterp[1][noz],prangedonor,transfo,profondeur,dirD,resInterp[2][noz],dirR[0],resInterp[2][noz].size,dim__[1]+1,dim__[2]+1,dim__[3]+1)

                            #print('dirR= ',dirR)
                            #print 'dimPb= ',dimPb
                            #print resInterp[3][noz]

                            connector.correctCoeffList(resInterp[1][noz],resInterp[3][noz],resInterp[2][noz],resInterp[2][noz].size,dim__[1]+1,dim__[2]+1,dim__[3]+1)

                            ### Determination du point pivot symetrique dans zoneR du point (imin,jmin,kmin) de la zoneD
                            pt_pivot=numpy.array(prange[0:3,0]) # Point (imin,jmin,kmin) de la zone R
                            if (abs(dirR)==1): # Le point ne peut varier que dans son plan (j,k) pour etre en face du point de la zone D
                                if (transfo[1] < 0) : pt_pivot[1] = prange[1,1] #Le vecteur j et son transforme sont opposes jmin -> jmax
                                if (transfo[2] < 0) : pt_pivot[2] = prange[2,1]
                            elif (abs(dirR)==2):  # Le point ne peut varier que dans son plan (i,k) pour etre en face du point de la zone D
                                if (transfo[0] < 0) : pt_pivot[0] = prange[0,1]
                                if (transfo[2] < 0) : pt_pivot[2] = prange[2,1]
                            else :  # Le point ne peut varier que dans son plan (i,j) pour etre en face du point de la zone D
                                if (transfo[0] < 0) : pt_pivot[0] = prange[0,1]
                                if (transfo[1] < 0) : pt_pivot[1] = prange[1,1]

                            info[2].append(['PointRangeDonor', prangedonor , [], 'IndexArray_t'])

                            info[2].append(['DirDonneur', dirD , [], 'IndexArray_t'])
                            info[2].append(['DirReceveur', dirR , [], 'IndexArray_t'])
                            info[2].append(['Transform', transfo , [], 'IndexArray_t'])
                            info[2].append(['PointPivot', pt_pivot , [], 'IndexArray_t'])
                            info[2].append(['Profondeur', profondeur , [], 'IndexArray_t'])
                            info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                            info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])

                            NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                            NMratio[0] = int(round(float(dim_[abs(transfo[0])])/float(dim__[1]+1)))
                            NMratio[1] = int(round(float(dim_[abs(transfo[1])])/float(dim__[2]+1)))
                            NMratio[2] = int(round(float(dim_[abs(transfo[2])])/float(dim__[3]+1)))
                            #print dim_[abs(transfo[0])], dim__[1]+1
                            #print dim_[abs(transfo[1])], dim__[2]+1
                            #print dim_[abs(transfo[2])], dim__[3]+1
                            #print profondeur
                            NMratio[abs(dirD)-1]=1
                            #print('NMratio= ',NMratio)

                            info[2].append(['NMratio', NMratio , [], 'IndexArray_t'])
                            info[2].append(['DnrZoneName', zonesDnr[noz][0] , [], 'IndexArray_t'])





    if storage != 'direct':
        OversetData._adaptForRANSLES__(tR, aD)

    # fin parcours des zones receveuses
    if storage == 'direct': return aR
    else:
        ztD = Internal.getZones(tD)
        zaD = Internal.getZones(aD)
        for i in range(len(ztD)): # enleve le cellN is interpData l'a ajoute
            ret = C.isNamePresent(ztD[i], 'cellN')
            if ret == -1: C._rmVars(zaD[i],['cellN'])
        return aD


#==============================================================================
# Cette fonction retourne la transformation vectorielle d'une zone a l'autre
#==============================================================================
def getTransfo(zdonor,zrcv):
    import KCore.Vector as Vector
    transfo = numpy.zeros(3,dtype=Internal.E_NpyInt)

    a = C.getFields(Internal.__GridCoordinates__, zdonor)[0]
    ni = a[2]; nj=a[3]; nk=a[4]

    if (nk == 1): #2D
        i = ni/2; j = nj/2; k = nk/2
        ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
        i=int(i);ip1=int(ip1)
        j=int(j);jp1=int(jp1)
        k=int(k);kp1=int(kp1)
        ind = i + j*ni + k*ni*nj
        P0 = [ a[1][0][ind], a[1][1][ind],0.0 ]
        ind = ip1 + j*ni + k*ni*nj
        P1 = [ a[1][0][ind], a[1][1][ind],0.0 ]
        ind = i + jp1*ni + k*ni*nj
        P2 = [ a[1][0][ind], a[1][1][ind],0.0 ]
        l1 = Vector.sub(P1,P0)
        l2 = Vector.sub(P2,P0)
        l1=Vector.normalize(l1)
        l2=Vector.normalize(l2)
        x=[1.0,0.0,0.0]
        y=[0.0,1.0,0.0]
        a=Vector.dot(l1,x);b=Vector.dot(l1,y)
        c=Vector.dot(l2,x);d=Vector.dot(l2,y)
        mat_ = numpy.array([[a,b],
                            [c,d]])

        a = C.getFields(Internal.__GridCoordinates__, zrcv)[0]
        ni = a[2]; nj=a[3]; nk=a[4]
        i = ni/2; j = nj/2; k = nk/2
        ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
        i=int(i);ip1=int(ip1)
        j=int(j);jp1=int(jp1)
        k=int(k);kp1=int(kp1)
        #ip1=i+1;jp1=j+1
        ind = i + j*ni + k*ni*nj
        #print ind
        P0_  = [ a[1][0][ind], a[1][1][ind],0.0 ]
        ind = ip1 + j*ni + k*ni*nj
        #print P0_
        #print ind
        P1_ = [ a[1][0][ind], a[1][1][ind],0.0 ]
        ind = i + jp1*ni + k*ni*nj
        #print P1_
        #print ind
        P2_ = [ a[1][0][ind], a[1][1][ind],0.0 ]
        l1rcv = Vector.sub(P1_,P0_)
        l2rcv = Vector.sub(P2_,P0_)
        l1rcv=Vector.normalize(l1rcv)
        l2rcv=Vector.normalize(l2rcv)
        #print l1rcv
        #print l2rcv
        x=[1.0,0.0,0.0]
        y=[0.0,1.0,0.0]
        a=Vector.dot(l1rcv,x);b=Vector.dot(l1rcv,y)
        c=Vector.dot(l2rcv,x);d=Vector.dot(l2rcv,y)
        mat = numpy.array([[a,b],
                          [c,d]])
        mat=mat.T
        a=mat_[0,0]*mat[0,0] + mat_[0,1]*mat[1,0]
        b=mat_[0,0]*mat[0,1] + mat_[0,1]*mat[1,1]
        c=mat_[1,0]*mat[0,0] + mat_[1,1]*mat[1,0]
        d=mat_[1,0]*mat[0,1] + mat_[1,1]*mat[1,1]
        mat__ = numpy.array([[a,b],
                            [c,d]])
        #print mat__[0,:]
        #print mat__[1,:]
        transfo[0]=(numpy.nonzero(mat__[0,:])[0][0]+1)*numpy.sign(mat__[0,numpy.nonzero(mat__[0,:])[0][0]])
        transfo[1]=(numpy.nonzero(mat__[1,:])[0][0]+1)*numpy.sign(mat__[1,numpy.nonzero(mat__[1,:])[0][0]])
        transfo[2]=3

        return transfo


    else : #3D

        i = ni/2; j = nj/2; k = nk/2
        ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
        i=int(i);ip1=int(ip1)
        j=int(j);jp1=int(jp1)
        k=int(k);kp1=int(kp1)
        ind = i + j*ni + k*ni*nj
        P0 = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = ip1 + j*ni + k*ni*nj
        P1 = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = i + jp1*ni + k*ni*nj
        P2 = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = i + j*ni + kp1*ni*nj
        P3 = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        l1 = Vector.sub(P1,P0)
        l2 = Vector.sub(P2,P0)
        l3 = Vector.sub(P3,P0)
        l1=Vector.normalize(l1)
        l2=Vector.normalize(l2)
        l3=Vector.normalize(l3)
        #print l1
        #print l2
        x=[1.0,0.0,0.0]
        y=[0.0,1.0,0.0]
        z=[0.0,0.0,1.0]
        a=Vector.dot(l1,x);b=Vector.dot(l1,y);c=Vector.dot(l1,z)
        d=Vector.dot(l2,x);e=Vector.dot(l2,y);f=Vector.dot(l2,z)
        g=Vector.dot(l3,x);h=Vector.dot(l3,y);i=Vector.dot(l3,z)
        mat_ = numpy.array([[a,b,c],
                            [d,e,f],
                            [g,h,i]])

        a = C.getFields(Internal.__GridCoordinates__, zrcv)[0]
        ni = a[2]; nj=a[3]; nk=a[4]
        i = ni/2; j = nj/2; k = nk/2
        ip1 = max(i+1,ni-1); jp1 = max(j+1,nj-1); kp1 = max(k+1,nk-1)
        i=int(i);ip1=int(ip1)
        j=int(j);jp1=int(jp1)
        k=int(k);kp1=int(kp1)
        ind = i + j*ni + k*ni*nj
        P0_  = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = ip1 + j*ni + k*ni*nj
        P1_ = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = i + jp1*ni + k*ni*nj
        P2_ = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        ind = i + j*ni + kp1*ni*nj
        P3_ = [ a[1][0][ind], a[1][1][ind],a[1][2][ind] ]
        l1rcv = Vector.sub(P1_,P0_)
        l2rcv = Vector.sub(P2_,P0_)
        l3rcv = Vector.sub(P3_,P0_)
        l1rcv=Vector.normalize(l1rcv)
        l2rcv=Vector.normalize(l2rcv)
        l3rcv=Vector.normalize(l3rcv)
        #print l1rcv
        #print l2rcv
        x=[1.0,0.0,0.0]
        y=[0.0,1.0,0.0]
        z=[0.0,0.0,1.0]
        a=Vector.dot(l1rcv,x);b=Vector.dot(l1rcv,y);c=Vector.dot(l1rcv,z)
        d=Vector.dot(l2rcv,x);e=Vector.dot(l2rcv,y);f=Vector.dot(l2rcv,z)
        g=Vector.dot(l3rcv,x);h=Vector.dot(l3rcv,y);i=Vector.dot(l3rcv,z)
        mat = numpy.array([[a,b,c],
                           [d,e,f],
                           [g,h,i]])
        #print 'mat= ',mat
        mat=mat.T
        a=mat_[0,0]*mat[0,0] + mat_[0,1]*mat[1,0] + mat_[0,2]*mat[2,0]
        b=mat_[0,0]*mat[0,1] + mat_[0,1]*mat[1,1] + mat_[0,2]*mat[2,1]
        c=mat_[0,0]*mat[0,2] + mat_[0,1]*mat[1,2] + mat_[0,2]*mat[2,2]
        d=mat_[1,0]*mat[0,0] + mat_[1,1]*mat[1,0] + mat_[1,2]*mat[2,0]
        e=mat_[1,0]*mat[0,1] + mat_[1,1]*mat[1,1] + mat_[1,2]*mat[2,1]
        f=mat_[1,0]*mat[0,2] + mat_[1,1]*mat[1,2] + mat_[1,2]*mat[2,2]
        g=mat_[2,0]*mat[0,0] + mat_[2,1]*mat[1,0] + mat_[2,2]*mat[2,0]
        h=mat_[2,0]*mat[0,1] + mat_[2,1]*mat[1,1] + mat_[2,2]*mat[2,1]
        i=mat_[2,0]*mat[0,2] + mat_[2,1]*mat[1,2] + mat_[2,2]*mat[2,2]

        mat__ = numpy.array([[a,b,c],
                             [d,e,f],
                             [g,h,i]])
        #print mat__[0,:]
        #print mat__[1,:]
        transfo[0]=(numpy.nonzero(mat__[0,:])[0][0]+1)*numpy.sign(mat__[0,numpy.nonzero(mat__[0,:])[0][0]])
        transfo[1]=(numpy.nonzero(mat__[1,:])[0][0]+1)*numpy.sign(mat__[1,numpy.nonzero(mat__[1,:])[0][0]])
        transfo[2]=(numpy.nonzero(mat__[2,:])[0][0]+1)*numpy.sign(mat__[2,numpy.nonzero(mat__[2,:])[0][0]])
        #print transfo

        return transfo
        #print 'mat= ',mat__


#------------------------------------------------------------------------------
# Creation or completion (if it already exists) of a node ZoneSubRegion_t for interpolation datas
# IN/OUT: z zone to be modified
# IN: zname: nom de la zone (donneuse si direct ou receveuse si inverse)
# IN: pointlist: liste des points ds PointList (points donneurs si stockage inverse, receveurs sinon)
# IN: pointlistdonor: liste des pts (points donneurs si stockage direct, receveurs sinon)
# IN: interpCoef: coefs d interpolation
# IN: interpType: type d interpolation
# IN: indicesExtrap: indices des pts extrapoles
# IN: indicesOrphan: indices de ts les pts orphelins (dupliques pour ttes les zones donneuses)
# IN: tag:'Receiver' or 'Donor'
# IN: loc='centers', 'nodes' ou 'faces'
# IN: EXDir: direction pour les pts EX (Chimere depth=1)
# IN: itype: 'abutting' pour transferts multidomaine, 'chimera' pour les interp chimere,
# 'ibc' pour les IBC
#------------------------------------------------------------------------------

#==============================================================================================================================
#==============================================================================================================================
#==============================================================================================================================
def setInterpDataForGhostCells2__(tR, tD, storage='direct', loc='nodes'):
    try: import Converter.GhostCells as GhostCells
    except: raise ImportError("setInterpDataForGhostCells__ requires Converter.GhostCells module.")
    # empty numpy arrays for zonesubregion nodes
    indicesExtrap = numpy.array([],Internal.E_NpyInt)
    indicesOrphan = numpy.array([],Internal.E_NpyInt)
    vols =  numpy.array([],numpy.float64)
    EXdir = numpy.array([],Internal.E_NpyInt)

    bases  = Internal.getNodesFromType1(tR     , 'CGNSBase_t')       # noeud
    dimmm  = Internal.getNodeFromName2(bases[0], 'EquationDimension')
    dimm   = Internal.getValue(dimmm)

    if loc == 'nodes': locR = 0; locS = 'Vertex'
    else: locR = 1; locS = 'CellCenter'

    aR = Internal.copyRef(tR)
    aD = Internal.copyRef(tD)
    for zp in Internal.getZones(aR):
        zname = zp[0]
        zoneDimR = Internal.getZoneDim(zp)
        if zoneDimR[0] == 'Unstructured':
            print('Warning: setInterpDataForGC not yet implemented for unstructured zones.')
        else: # Structured
            dimPb = zoneDimR[4]
            rindnode = Internal.getNodeFromType1(zp, 'Rind_t')
            if rindnode is not None: # rind indices exist : ghost cell data to be computed
                rindnode = rindnode[1]
                rindimin = rindnode[0]; rindimax = rindnode[1]
                rindjmin = 0; rindjmax = 0; rindkmin = 0; rindkmax = 0
                if dimPb > 1:
                    rindjmin = rindnode[2]; rindjmax = rindnode[3]
                    if dimPb == 3: rindkmin = rindnode[4]; rindkmax = rindnode[5]

                rindrcv = [rindimin,rindimax,rindjmin,rindjmax,rindkmin,rindkmax]
                imr = zoneDimR[1]; jmr = zoneDimR[2]; kmr = zoneDimR[3]
                if locR == 1: imr=imr-1; jmr=jmr-1; kmr=kmr-1
                listofjoins = Internal.getNodesFromType2(zp, 'GridConnectivity1to1_t')
                for join in listofjoins:
                    # receptor window ranges
                    prange = Internal.getNodeFromName1(join,'PointRange')[1]
                    # donor window ranges
                    prangedonor = Internal.getNodeFromName1(join,'PointRangeDonor')[1]

                    # trirac
                    if dimPb == 3: trirac = [1,2,3]
                    elif dimPb == 2: trirac = [1,2]
                    else: trirac = [1]
                    transfo = Internal.getNodeFromName1(join, 'Transform')
                    if transfo is not None:
                        trirac[0] = transfo[1][0]
                        if dimPb != 1: trirac[1] = transfo[1][1]
                        if dimPb == 3: trirac[2] = transfo[1][2]
                    Periodic = Internal.getNodeFromType2(join,'Periodic_t')
                    RotationAngle=None; RotationCenter=None
                    if Periodic is not None:
                        RotationAngle = Internal.getNodeFromName1(Periodic,'RotationAngle')
                        RotationCenter = Internal.getNodeFromName1(Periodic,'RotationCenter')
                        if RotationAngle is not None:
                            RotationAngle[1][0]=-RotationAngle[1][0]
                            RotationAngle[1][1]=-RotationAngle[1][1]
                            RotationAngle[1][2]=-RotationAngle[1][2]

                    # donor zone name
                    zdonorname = Internal.getValue(join)
                    zdonor = Internal.getZones(aR)
                    zdonor = Internal.getNodesFromName(zdonor, zdonorname)
                    zdonor = zdonor[0] # Warning: if a zone name is not unique, the first zone is kept
                    # check if the donor zone is defined in aD
                    zdonorp = Internal.getNodesFromName2(aD, zdonorname)
                    if zdonorp == []:
                        raise ValueError("setInterpDataForGhostCells: donor zone not found in donor pyTree.")
                    zdonorp = Internal.getNodesFromType2(zdonorp, 'Zone_t')
                    if zdonorp == []:
                        raise ValueError("setInterpDataForGhostCells: donor zone not found in donor pyTree.")
                    if len(zdonorp)  > 1 :
                        print('Warning: setInterpDataForGhostCells: zone name %s defined several times in donor pyTree.'%zdonorname)
                    zdonorp = zdonorp[0]
                    # donor zone dimensions
                    zoneDimD = Internal.getZoneDim(zdonor)
                    imd = zoneDimD[1]; jmd = zoneDimD[2]; kmd = zoneDimD[3]
                    if locR == 1: imd=imd-1; jmd=jmd-1; kmd=kmd-1

                    # rind for donor
                    rindnode = Internal.getNodeFromType1(zdonor, 'Rind_t')
                    if rindnode is not None: # rind indices exist: ghost cell data to be computed
                        rindnode = rindnode[1]
                        rindimin = rindnode[0]; rindimax = rindnode[1]
                        rindjmin = 0; rindjmax = 0; rindkmin = 0; rindkmax = 0
                        if dimPb > 1:
                            rindjmin = rindnode[2]; rindjmax = rindnode[3]
                            if dimPb == 3: rindkmin = rindnode[4]; rindkmax = rindnode[5]
                        #
                        rinddnr = [rindimin,rindimax,rindjmin,rindjmax,rindkmin,rindkmax]
                        # get directions of receptor and donor borders
                        dirR = GhostCells.getDirBorderStruct__(prange,dimPb)
                        dirD = GhostCells.getDirBorderStruct__(prangedonor,dimPb)

                        # get list of border nodes and direction of border
                        if dimPb != 1:
                            # arrayOfIndicesR : indices globaux des 1ers points pres de la frontiere
                            # definie par le prange dimensionnes dans le maillage ghost cells
                            shift = 0 # les indices sont ceux des pts en frontiere max GC
                            #arrayOfIndicesR = GhostCells.getBorderIndicesStruct__(prange,zoneDimR,dirR,0,locS,dimPb,shift)
                            [ arrayOfIndicesR, dim1, dim2]  = GhostCells.getBorderIndicesStruct__(prange,zoneDimR,dirR,0,locS,dimPb,shift)


                            # listOfIndicesD : indices globaux des 1ers pts donneurs associes a ceux definis par
                            # arrayOfIndicesR
                            shift = 0
                            if dirD  ==-1: shift += rinddnr[0]
                            elif dirD== 1: shift += rinddnr[1]
                            elif dirD==-2: shift += rinddnr[2]
                            elif dirD== 2: shift += rinddnr[3]
                            elif dirD==-3: shift += rinddnr[4]
                            else: shift+= rinddnr[5]
                            if dirR  ==-1: shift += rindimin
                            elif dirR== 1: shift += rindimax
                            elif dirR==-2: shift += rindjmin
                            elif dirR== 2: shift += rindjmax
                            elif dirR==-3: shift += rindkmin
                            else: shift+= rindkmax
                            if locR == 1: shift -= 1
                            listOfIndicesD = GhostCells.getJoinDonorIndicesStruct__(prange,prangedonor,zoneDimD,dirD, trirac,0,locS,dimPb, shift)

                    # Increments
                    absdirR=abs(dirR); absdirD=abs(dirD)
                    # list of directions of receptor zone
                    dir_rcv = [0]*dimPb; dir_rcv[absdirR-1] = 1
                    # list of direction of donor zone
                    dir_dnr = [0]*dimPb; dir_dnr[absdirD-1] = 1
                    if dimPb == 3:
                        incrR = GhostCells.increment__(dir_rcv[0],dir_rcv[1],dir_rcv[2], imr, jmr, kmr, 0)
                        incrD = GhostCells.increment__(dir_dnr[0],dir_dnr[1],dir_dnr[2], imd, jmd, kmd, 0)
                    elif dimPb == 2:
                        incrR = GhostCells.increment__(dir_rcv[0],dir_rcv[1], imr, jmr, 0)
                        incrD = GhostCells.increment__(dir_dnr[0],dir_dnr[1], imd, jmd, 0)
                    if dirR > 0: incrR = -incrR
                    if dirD < 0: incrD = -incrD

                    if dirR  ==-1: depth = rinddnr[0]
                    elif dirR== 1: depth = rinddnr[1]
                    elif dirR==-2: depth = rinddnr[2]
                    elif dirR== 2:  depth = rinddnr[3]
                    elif dirR==-3: depth = rinddnr[4]
                    else: depth = rinddnr[5]

                    res = connector.setInterpDataForGC(arrayOfIndicesR, listOfIndicesD,
                                                       dimPb, locR, depth, incrR, incrD)
                    # Stockage
                    vol = numpy.array([],numpy.float64)
                    prefix = 'ID_'
                    if RotationAngle is not None:
                        val = Internal.getValue(RotationAngle)
                        if val[0]>0. or val[1]>0. or val[2]>0.: prefix='IDPERP'
                        else: prefix = 'IDPERM'
                    if storage == 'direct':
                        OversetData._createInterpRegion__(zp, zdonorname,res[0],res[1],res[3],res[2],vols,indicesExtrap,\
                                                          indicesOrphan,tag='Receiver',loc=loc,EXDir=EXdir,itype='abutting',\
                                                          prefix=prefix, RotationAngle=RotationAngle, RotationCenter=RotationCenter)

                    else:
                        OversetData._createInterpRegion__(zdonorp, zname, res[1], res[0], res[3], res[2], vols, indicesExtrap,\
                                                          indicesOrphan, tag='Donor',loc=loc,EXDir=EXdir,itype='abutting',\
                                                          prefix=prefix,RotationAngle=RotationAngle, RotationCenter=RotationCenter)

                    zoneDimR2 = list(zoneDimR)
                    zoneDimR2[4] = dimm

                    zoneDimD2 = list(zoneDimD)
                    zoneDimD2[4] = dimm

                    #print  'zoneDimD2= ', zoneDimD2

                    GhostCells._adaptBCStruct3__(join,zoneDimR2,zoneDimD2, -2)

                    transfo = Internal.getNodeFromName1(join, 'Transform')

                    #print 'transfo= ', transfo[1]

                    prange = Internal.getNodeFromName1(join,'PointRange')[1]
                    #direction = GhostCells.getDirBorderStruct__(prange,zoneDimR2[4])
                    #prange=numpy.append(prange,[direction,0])

                    #print 'prangeav= ',prange
                    ### passage en indices "cell centered"
                    #for j in range(2):
                    for i in range(3):
                        if prange[i,1] == prange[i,0] and prange[i,1] != 1:
                            prange[i,1] =  prange[i,1]-1
                            prange[i,0] =  prange[i,0]-1
                        elif prange[i,1] != prange[i,0] and prange[i,1] != 1 :
                            prange[i,1] =  prange[i,1]-1


                    #print 'prangeap= ', prange


                    prangedonor = Internal.getNodeFromName1(join,'PointRangeDonor')[1]

                    #directiond = GhostCells.getDirBorderStruct__(prangedonor,zoneDimD2[4])
                    #prangedonor=numpy.append(prangedonor,[directiond,0])

                    ### passage en indices "cell centered"
                    for i in range(3):
                        if prangedonor[i,1] == prangedonor[i,0] and prangedonor[i,1] != 1:
                            prangedonor[i,1] =  prangedonor[i,1]-1
                            prangedonor[i,0] =  prangedonor[i,0]-1
                        elif prangedonor[i,1] != prangedonor[i,0] and prangedonor[i,1] != 1 :
                            prangedonor[i,1] =  prangedonor[i,1]-1


                    ### Determination du point pivot symetrique dans zoneR du point (imin,jmin,kmin) de la zoneD
                    pt_pivot=numpy.array(prange[0:3,0]) # Point (imin,jmin,kmin) de la zone R
                    if abs(dirR)==1: # Le point ne peut varier que dans son plan (j,k) pour etre en face du point de la zone D
                        if (transfo[1][1] < 0) : pt_pivot[1] = prange[1,1] #Le vecteur j et son transforme sont opposes jmin -> jmax
                        if (transfo[1][2] < 0) : pt_pivot[2] = prange[2,1]
                    elif (abs(dirR)==2):  # Le point ne peut varier que dans son plan (i,k) pour etre en face du point de la zone D
                        if (transfo[1][0] < 0) : pt_pivot[0] = prange[0,1]
                        if (transfo[1][2] < 0) : pt_pivot[2] = prange[2,1]
                    else :  # Le point ne peut varier que dans son plan (i,j) pour etre en face du point de la zone D
                        if (transfo[1][0] < 0) : pt_pivot[0] = prange[0,1]
                        if (transfo[1][1] < 0) : pt_pivot[1] = prange[1,1]



                    #print 'pt_pivot_ap= ',pt_pivot
                    # Inversion du vector transform pour l explicite local
                    transfo_inv=numpy.array(transfo[1][0:3])
                    for i in range(0,3):
                        transfo_inv[abs(transfo[1][i])-1]=i+1
                        if transfo[1][i] < 0: transfo_inv[abs(transfo[1][i])-1]=-(i+1)

                    dirR = GhostCells.getDirBorderStruct__(prange,dimm)
                    dirD = GhostCells.getDirBorderStruct__(prangedonor,dimm)

                    #print 'dirR= ', dirR
                    #print 'dirD= ', dirD

                    prange = numpy.reshape(prange,6)
                    prangedonor= numpy.reshape(prangedonor,6)

                    profondeur=numpy.zeros(1,dtype=Internal.E_NpyInt)
                    profondeur[0]=1
                    #transfo = Internal.getNodeFromName1(join, 'Transform')
                    #print 'transfo= ', transfo[1]
                    #print 'pt_pivot_ap= ',pt_pivot
                    NMratio = numpy.zeros(3,dtype=Internal.E_NpyInt)
                    NMratio[0]=1
                    NMratio[1]=1
                    NMratio[2]=1

                    #print 'zname= ', zname
                    #print zdonor
                    z_     = Internal.getNodeFromName2(aR, zname)
                    node   = Internal.getNodeFromName1(z_, 'FlowSolution#Centers')
                    sol    = Internal.getNodeFromName1(node,  'niveaux_temps')
                    levelrcv  = Internal.getValue(sol)[0][0][0]

                    z_     = Internal.getNodeFromName2(aR, zdonorname)
                    node   = Internal.getNodeFromName1(z_, 'FlowSolution#Centers')
                    sol    = Internal.getNodeFromName1(node,  'niveaux_temps')
                    leveldnr  = Internal.getValue(sol)[0][0][0]


                    if storage == 'direct':
                        info = zp[2][len(zp[2])-1]
                        info[2].append(['PointRange',      prange , [], 'IndexArray_t'])
                        info[2].append(['PointRangeDonor',      prangedonor , [], 'IndexArray_t'])
                        info[2].append(['DirReceveur',      dirR , [], 'IndexArray_t'])
                        info[2].append(['DirDonneur',      dirD , [], 'IndexArray_t'])
                        info[2].append(['Transform',      transfo_inv , [], 'IndexArray_t'])
                        info[2].append(['PointPivot',      pt_pivot , [], 'IndexArray_t'])
                        info[2].append(['Profondeur',      profondeur , [], 'IndexArray_t'])
                        info[2].append(['NMratio',      NMratio , [], 'IndexArray_t'])
                        info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                        info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])

                    else:
                        info = zdonorp[2][len(zdonorp[2])-1]
                        info[2].append(['PointRange',      prange , [], 'IndexArray_t'])
                        info[2].append(['PointRangeDonor',     prangedonor , [], 'IndexArray_t'])
                        info[2].append(['DirReceveur',      dirR , [], 'IndexArray_t'])
                        info[2].append(['DirDonneur',      dirD , [], 'IndexArray_t'])
                        info[2].append(['Transform',      transfo_inv , [], 'IndexArray_t'])
                        info[2].append(['PointPivot',      pt_pivot , [], 'IndexArray_t'])
                        info[2].append(['Profondeur',      profondeur , [], 'IndexArray_t'])
                        info[2].append(['NMratio',      NMratio , [], 'IndexArray_t'])
                        info[2].append(['LevelZRcv', levelrcv , [], 'IndexArray_t'])
                        info[2].append(['LevelZDnr', leveldnr , [], 'IndexArray_t'])

    if storage == 'direct': return aR
    else: return aD

#===============================================================================================================
#===============================================================================================================
#===============================================================================================================
