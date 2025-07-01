# - OversetDataElsA -
# Module for internal functions used for overset info
from . import Connector
import numpy

try:
    import Converter.Internal as Internal
    import Converter.PyTree as C
    import Converter
    from . import OversetData as XOD
except:
    raise ImportError("Connector.OversetDataElsA requires Converter and Connector.OversetData modules.")

#=============================================================================
#=============================================================================
# 1. CALCUL DES DONNEES D INTERPOLATIONS
#=============================================================================
#=============================================================================


#=============================================================================
# Calcul des coefficients d'interpolation et stockage dans l'arbre CGNS/Python
# ----------------------------------------------------------------------------
# setSeqInterpolation: calcul en sequentiel (stockage direct ou inverse
# setDistInterpolation: calcul dans un environnement distribue (stockage de type inverse)
# -------------------
# le celln doit etre defini dans l arbre et valoir :
# 0 pour les pts masques, 2 pour les pts interpoles
# -------------------
# IN: t: arbre sur lequel les interpolations sont calculees
# IN: loc: localisation des interpolations ('cell':cellules ou 'face':faces)
# IN: double_wall: activation de la technique double wall
# IN: prefixFile: prefixe pour le nom du fichier elsA
# IN: sameBase: si sameBase=1, les zones d'une meme base peuvent etre des zones d'interpolation
# IN: parallelDatas: liste de donnees a fournir dans un contexte parallele (=[] en sequentiel)
#     la liste contient: la liste de voisinage, le rang et la liste des cellules ou points EX interpoles.
# -------------------
# prefixFile == '': stockage dans l'arbre CGNS/Python uniquement
# prefixFile != '': stockage dans l'arbre CGNS/Python et ecriture dans un
#                   fichier relisible par elsA
#==============================================================================
def setInterpolations(t, loc='cell', double_wall=0, storage='inverse', prefixFile='',
                      sameBase=0, solver='elsA', nGhostCells=2, parallelDatas=[], cfMax=30., planarTol=0., check=True):
    """Compute interpolation data for chimera and for elsA solver."""
    a = Internal.copyRef(t)
    _setInterpolations(a, loc=loc, double_wall=double_wall, storage=storage, prefixFile=prefixFile,
                       sameBase=sameBase, solver=solver, nGhostCells=nGhostCells, parallelDatas=parallelDatas,
                       cfMax=cfMax, planarTol=planarTol, check=check)
    return a

def _setInterpolations(t, loc='cell', double_wall=0, storage='inverse', prefixFile='', sameBase=0,
                       solver='elsA', nGhostCells=2, parallelDatas=[], cfMax=30., planarTol=0., check=True):
    if storage != 'direct' and storage != 'inverse':
        raise TypeError("Warning: setInterpolations: bad value for attribute 'storage'.")

    if storage == 'direct':
        if prefixFile == '':
            print('Warning: _setInterpolations: inverse storage is mandatory if no Chimera connectivity files are written.')
        elif parallelDatas != []:
            print('Warning: _setInterpolations: inverse storage is activated (mandatory in a distributed mode).')
            storage='inverse'

    # Solveur :
    if solver == 'elsA': Solver = 1
    elif solver == 'Cassiopee': Solver = 2

    # Localisation
    if loc == 'cell': depth=2
    elif loc == 'face': depth=1
    else: raise TypeError("setInterpolations: bad value for attribute 'loc'.")
    if parallelDatas == []: # mode sequentiel
        _setSeqInterpolations(t, depth=depth, double_wall=double_wall, storage=storage, prefixFile=prefixFile,
                              sameBase=sameBase, solver=Solver, nGhostCells=nGhostCells, cfMax=cfMax, planarTol=planarTol,
                              check=check)
    else: # mode distribue
        if nGhostCells != 2: print('Warning: _setInterpolations: nGhostCells must be 2 in distributed mode.')
        _setDistInterpolations(t, parallelDatas, depth, double_wall, sameBase, Solver, cfMax, check=check)
    return None

#==============================================================================
# Calcul des coefficients d'interpolation et stockage dans l'arbre CGNS/Python
# dans un environnement distribue
#==============================================================================
def _setDistInterpolations(a, parallelDatas=[], depth=2, double_wall=0,
                           sameBase=0, solver=1, cfMax=30., check=True):

    try: import Generator.PyTree as G
    except: raise ImportError("_setDistInterpolations: requires Generator module.")

    if double_wall == 1: from . import DoubleWall
    print("Warning: _setDistInterpolations: periodic Chimera not yet implemented.")

    if parallelDatas == []: return None
    else:
        if len(parallelDatas) != 3:
            raise TypeError("setDistInterpolations: missing datas in parallelDatas.")
        else:
            graph = parallelDatas[0]
            rank = parallelDatas[1]
            interpDatas = parallelDatas[2]

    if rank not in graph:
        localGraph = {}
    else:
        localGraph = graph[rank]
    # ----------------------------------------
    # Initialisation
    # ----------------------------------------
    # creation et initialisation a 1 du champ cellN s'il n'est pas deja present dans l'arbre
    XOD._addCellN__(a); XOD._addCellN__(a, loc='centers')

    # ----------------------------------------------------------------------------------------
    # Calcul des cellules et coefficients d'interpolation et ecriture dans un arbre CGNS/Python
    # ----------------------------------------------------------------------------------------
    # Calcul des donnees sur les domaines d interpolations (coordonnees, noms des zones, cellN, surfaces)
    # initialisation de liste de parois et de surfaces pour le traitement 'double_wall'
    firstWallCenters = []; surfacesExtC = []; walls1 = []

    if double_wall == 1:
        dwInfo = DoubleWall.extractDoubleWallInfo__(a)
        firstWallCenters = dwInfo[0]; surfacesExtC = dwInfo[1]
    interpolationZonesD={}; interpolationZonesNameD={}; interpolationCellnD={}; surfiD={}

    # initialisation des cles
    for oppNode in localGraph:
        oppZones = graph[oppNode][rank]
        for s in range(len(interpDatas[oppNode])):
            oppZone = oppZones[s]
            interpolationZonesD[oppZone] = []
            interpolationZonesNameD[oppZone] = []
            interpolationCellnD[oppZone] = []
            surfiD[oppZone] = []

    bases = Internal.getBases(a); nbases = len(bases)
    # remplissage des dictionnaires
    for oppNode in localGraph:
        oppZones = graph[oppNode][rank]
        for s in range(len(interpDatas[oppNode])):
            oppZone = oppZones[s]
            for nob in range(nbases):
                zones = Internal.getZones(bases[nob])
                for noz in range(len(zones)):
                    z = zones[noz]
                    zn = C.getFields(Internal.__GridCoordinates__,z)[0]
                    cellN = C.getField('centers:cellN', z)[0]
                    interpolationZonesNameD[oppZone].append(z[0])
                    interpolationZonesD[oppZone].append(zn)
                    interpolationCellnD[oppZone].append(cellN)
                    if (surfacesExtC != []): surfiD[oppZone].append(surfacesExtC[nob][noz])
                    del zn

    # 1. deals with depth=2
    if depth == 2:
        for oppNode in localGraph:
            oppZones = graph[oppNode][rank]
            for s in range(len(interpDatas[oppNode])):
                interpCells=interpDatas[oppNode][s]
                oppZone = oppZones[s]
                interpolationZonesName = interpolationZonesNameD[oppZone]
                interpolationZones=interpolationZonesD[oppZone]
                interpolationCelln=interpolationCellnD[oppZone]
                surfi = surfiD[oppZone]
                # recuperation des domaines d interpolations
                # Calcul des cellules d interpolations et des coefficients associes
                # -----------------------------------------------------------------
                nBlksI = len(interpolationZones)
                if nBlksI > 0:
                    listOfInterpCells = []
                    # traitement double_wall
                    for i in range(nBlksI):
                        # changewall des centres et EX si paroi existe pour le domaine d interpolation
                        if surfi ==[] or walls1 == []:# attention: walls1 est tjs [] ici, pourquoi ?
                            listOfInterpCells.append(interpCells)
                        else:
                            if surfi[i] == []:
                                listOfInterpCells.append(interpCells)
                            else: # surface existe
                                ### DBG TO DO ...
                                #zc2 = Connector.changeWall__(zc, firstWallCenters1, surfi[i])
                                #interpCells = Connector.getInterpolatedPoints__(zc2)
                                #del zc2
                                listOfInterpCells.append(interpCells)
                    resInterp = Connector.setInterpolations__(oppZone, 1, 1, listOfInterpCells,
                                                              interpolationZones, interpolationCelln, isEX=0, cfMax=cfMax,check=check)
                    del listOfInterpCells

                    for nozd in range(nBlksI):
                        zdonorname = interpolationZonesName[nozd]
                        zdonor = Internal.getNodeFromName(a, zdonorname)
                        _interpInverseStorage(oppZone,zdonor,nozd,resInterp,depth)
    else:
        #2. deals with depth =1
        for oppNode in localGraph:
            oppZones = graph[oppNode][rank]
            for s in range(len(interpDatas[oppNode])):
                EXPts=interpDatas[oppNode][s]
                oppZone = oppZones[s]
                interpolationZones=interpolationZonesD[oppZone]
                interpolationCelln=interpolationCellnD[oppZone]
                interpolationZonesName = interpolationZonesNameD[oppZone]
                surfi = surfiD[oppZone]

                # recuperation des domaines d interpolations
                # Calcul des cellules d interpolations et des coefficients associes
                # -----------------------------------------------------------------
                nBlksI = len(interpolationZones)
                if nBlksI > 0:
                    listOfEXPts = [];
                    for i in range(nBlksI):
                        if surfi == [] or walls1 == []:
                            listOfEXPts.append(EXPts);
                        else:
                            if surfi[i] == []:
                                listOfEXPts.append(EXPts)
                            else:
                                ### DBG TO DO ...
                                #expts = Connector.changeWallEX__(EXPts, zc, zn, walls1, surfi[i],\
                                #                                 doubleWallTol)
                                #listOfEXPts.append(expts);
                                listOfEXPts.append(EXPts);
                    resInterp = Connector.setInterpolations__(oppZone, 1, 1, listOfEXPts, interpolationZones, \
                                                              interpolationCelln, isEX=1, cfMax=cfMax, check=check)
                    del listOfEXPts

                    for nozd in range(nBlksI):
                        zdonorname = interpolationZonesName[nozd]
                        zdonor = Internal.getNodeFromName(a, zdonorname)
                        _interpInverseStorage(oppZone,zdonor,nozd,resInterp,depth)
    return None

#=============================================================================
# Calcul des coefficients d'interpolation et stockage dans l'arbre CGNS/Python
# dans un environnement sequentiel
#==============================================================================
def _setSeqInterpolations(a, depth=2, double_wall=0, storage='inverse', prefixFile='',
                          sameBase=0, solver=1, nGhostCells=2,
                          cfMax=30., planarTol=0., check=True):
    import Generator.PyTree as G
    if double_wall == 1: from . import DoubleWall

    # Ecriture ou non dans des fichiers elsA
    writeInFile=0
    if prefixFile != '': writeInFile=1
    # ----------------------------------------
    # Initialisation
    # ----------------------------------------
    # creation et initialisation a 1 du champ cellN s'il n'est pas deja present dans l'arbre
    XOD._addCellN__(a, loc='centers') # aux centres
    bases = Internal.getBases(a); nbases = len(bases)

    # Correspondance entre une zone et son Id globale pour l'ecriture dans les fichiers elsA
    ZonesId={}; Id = 0
    if writeInFile == 1:
        # declaration de dictionnaire pour le stockage de donnees par bloc donneur pour l'ecriture en fichier elsA
        listRcvId = {}; nbInterpCellsForDonor = {}; listIndicesRcv = {}; listEXdir={}; listIndicesDonor = {}; listInterpolantsDonor = {}; listCellN={}; listZoneDim={}; listInterpTypes = {}
    for nob in range(nbases):
        zones = Internal.getNodesFromType1(bases[nob], 'Zone_t')
        for noz in range(len(zones)):
            ZonesId[zones[noz][0]]=Id; Id = Id+1
    ntotZones = len(ZonesId)

    #=======================================================================
    # 2 - Recherche des periodicites :
    #     duplication des blocs periodiques dans les bases associees
    #     creation d un noeud fils au niveau de la zone dupliquee de nom 'TemporaryPeriodicZone'
    #=======================================================================
    for nob in range(len(a[2])):
        if a[2][nob][3]=='CGNSBase_t':
            C._addPeriodicZones__(a[2][nob])
    # initialisation de liste de parois et de surfaces pour le traitement 'double_wall'
    wallBndIndicesN=[]; surfacesExtC=[]; surfacesN=[]

    # -------------------------------------------------------------------------
    # Calcul des cellules et coefficients d interpolation et ecriture dans un
    # arbre CGNS/Python
    # -------------------------------------------------------------------------
    if writeInFile == 1:
        nbInterpCellsForDonor = [0]*ntotZones
    if double_wall == 1:
        dwInfo = DoubleWall.extractDoubleWallInfo__(a)
        firstWallCenters = dwInfo[0]; surfacesExtC = dwInfo[1]

    intersectionDict = XOD.getIntersectingDomains(a, method='AABB')

    for nob in range(nbases):
        if bases[nob][0] == 'CARTESIAN': sameBase = 1 # peut etre a enlever ??
        zones = Internal.getNodesFromType1(bases[nob], 'Zone_t')
        nzones = len(zones)
        for noz in range(nzones):
            z = zones[noz]
            zname = z[0]
            if Internal.getNodeFromName1(z,'TempPeriodicZone') is not None: break

            # calcul d'arrays contenant les coefficients d'interpolation et les indices des cellules donneuses
            zn = C.getFields(Internal.__GridCoordinates__, z)[0]
            nirc = Internal.getZoneDim(z)[1]-1; njrc = Internal.getZoneDim(z)[2]-1; nkrc = Internal.getZoneDim(z)[3]-1
            if writeInFile == 1: listZoneDim[ZonesId[z[0]]]=[nirc,njrc,nkrc]
            # calcul des cellules d interpolations et des coefficients d interpolations
            zc = Converter.node2Center(zn)
            cellN = C.getField('centers:cellN',z)[0] # celln aux centres des cellules
            if writeInFile == 1:
                listCellN[ZonesId[z[0]]]=C.getField('centers:cellN',z, api=2)[0]#cellN[1]
            zc = Converter.addVars([zc,cellN]); del cellN
            firstWallCenters1 = []
            if double_wall == 1: firstWallCenters1 = firstWallCenters[nob][noz]

            # traitement des centres
            interpCells = Connector.getInterpolatedPoints__(zc)
            interpCells = Converter.extractVars(interpCells,['CoordinateX','CoordinateY','CoordinateZ','indcell'])
            # traitement des pts EX
            if depth == 1: EXPts = Connector.getEXPoints__(zn, zc)
            if interpCells[1].size != 0: # pts interpoles
                # recuperation des domaines d interpolations
                res = getInterpolationDomains__(bases,nob,noz,zn,surfacest=surfacesExtC,sameBase=sameBase, intersectionDict=intersectionDict[zname])
                interpolationZonesName = res[0]; interpolationZones = res[1]; interpolationCelln = res[2]; surfi = res[3]; periodicZones = res[4]
                nBlksI = len(interpolationZones)
                # Calcul des cellules d'interpolations et des coefficients associes
                # -----------------------------------------------------------------
                if nBlksI > 0:
                    if depth == 2:
                        listOfInterpCells = []
                        # traitement double_wall
                        for i in range(nBlksI):
                            # changewall des centres et EX si paroi existe pour le domaine d interpolation
                            if (firstWallCenters1 == [] or surfi == []):
                                listOfInterpCells.append(interpCells)
                            else:
                                if surfi[i] == []:
                                    listOfInterpCells.append(interpCells)
                                else: # surface existe
                                    zc2 = Connector.changeWall__(zc, firstWallCenters1, surfi[i], planarTol)
                                    interpCells = Connector.getInterpolatedPoints__(zc2)
                                    interpCells = Converter.extractVars(interpCells,['CoordinateX','CoordinateY','CoordinateZ','indcell'])
                                    del zc2
                                    listOfInterpCells.append(interpCells)
                        resInterp = Connector.setInterpolations__(z[0],nirc, njrc, listOfInterpCells, interpolationZones, \
                                                                  interpolationCelln, isEX=0, cfMax=cfMax, zoneId=ZonesId[z[0]]+1, check=check)
                        del listOfInterpCells

                    # traitement des pts EX
                    elif depth == 1:
                        listOfEXPts = []
                        for i in range(nBlksI):
                            zdonorId = ZonesId[interpolationZonesName[i]]
                            if (surfi == [] or firstWallCenters1 == []): listOfEXPts.append(EXPts)
                            else:
                                if surfi[i] == []: listOfEXPts.append(EXPts)
                                else:
                                    expts = Connector.changeWallEX__(EXPts, zc, zn, firstWallCenters1, surfi[i], planarTol)
                                    listOfEXPts.append(expts)
                        del zn
                        resInterp = Connector.setInterpolations__(z[0],nirc, njrc, listOfEXPts, interpolationZones, \
                                                                  interpolationCelln, isEX=1, cfMax=cfMax, zoneId=ZonesId[z[0]]+1, check=check)
                        del listOfEXPts
                    for nozd in range(nBlksI):
                        indicesRcv = resInterp[0][nozd]
                        indicesDonor=resInterp[1][nozd]
                        indicesDonorExtC = resInterp[2][nozd]
                        interpCoef=resInterp[3][nozd];
                        interpVol=resInterp[4][nozd]; interpType=resInterp[5][nozd]
                        indicesExtrap = resInterp[6][nozd]; indicesOrphan = resInterp[7][nozd]
                        if depth == 1: EXdir = resInterp[8][nozd]
                        zdonorname = interpolationZonesName[nozd]
                        isperiodic = periodicZones[nozd]
                        if isperiodic == 2:
                            print('Periodic interpolation from +theta: %d.'%interpType.shape[0])
                            interpType = 102*numpy.ones((interpType.shape[0]), Internal.E_NpyInt)
                        elif isperiodic == 3:
                            print('Periodic interpolation from -theta: %d.'%interpType.shape[0])
                            interpType = 103*numpy.ones((interpType.shape[0]), Internal.E_NpyInt)
                        resInterp[5][nozd] = interpType
                        zdonor = Internal.getNodesFromName(a,zdonorname)[0]
                        # ----------------------------------------
                        # stockage direct (sur le bloc interpole)
                        # ----------------------------------------
                        if storage == 'direct':
                            _interpDirectStorage(z,zdonorname,nozd,resInterp,depth)
                        # -------------------------------------
                        # stockage inverse (sur le bloc donneur)
                        # -------------------------------------
                        else:
                            _interpInverseStorage(z[0],zdonor,nozd,resInterp,depth)
                        # -----------------------------------------------------
                        # Stockage de donnees pour l ecriture des fichiers elsA
                        # -----------------------------------------------------
                        # - nombre total de points interpoles par le bloc donneur zdonorname
                        # - liste des donnees a ecrire par bloc donneur (indices cell. interpolees, indices cell. donneuses, coefficients d interpolation)
                        if writeInFile == 1:
                            donorId = ZonesId[zdonorname]
                            # dimensions de la grille d interpolation
                            ni = Internal.getZoneDim(zdonor)[1]; nj = Internal.getZoneDim(zdonor)[2]; nk = Internal.getZoneDim(zdonor)[3]
                            nbInterpCellsForDonor[donorId] = nbInterpCellsForDonor[donorId] + len(indicesDonor)
                            if donorId in listRcvId:
                                listRcvId[donorId].append(ZonesId[z[0]])
                                listIndicesRcv[donorId].append(indicesRcv)
                                listIndicesDonor[donorId].append(indicesDonorExtC)
                                listInterpolantsDonor[donorId].append(interpCoef)
                                listInterpTypes[donorId].append(interpType)
                            else:
                                listRcvId[donorId]=[ZonesId[z[0]]]
                                listIndicesRcv[donorId]=[indicesRcv]
                                listIndicesDonor[donorId]=[indicesDonorExtC]
                                listInterpolantsDonor[donorId]=[interpCoef]
                                listInterpTypes[donorId]=[interpType]

                            if depth ==1:
                                if donorId in listEXdir: listEXdir[donorId].append(EXdir)
                                else: listEXdir[donorId]=[EXdir]

                else: # no interpolation domain found
                    if depth == 2: indicesOrphan = numpy.array(interpCells[1][3], dtype=Internal.E_NpyInt)
                    elif depth == 1: indicesOrphan = numpy.array(EXPts[1][3], dtype=Internal.E_NpyInt)
                    # on cree une zone subregion avec les pts orphelins
                    nameSubRegion='Orphan_'+z[0]
                    z[2].append([nameSubRegion, None, [],'ZoneSubRegion_t'])
                    info = z[2][len(z[2])-1]
                    Internal._createChild(info, 'ZoneRole', 'DataArray_t', value='Receiver')
                    if depth == 1: info[2].append(['OrphanFaceList',numpy.unique(indicesOrphan), [], 'DataArray_t'])
                    else: info[2].append(['OrphanPointList', numpy.unique(indicesOrphan), [], 'DataArray_t'])
    # Delete duplicated periodic zones
    C._removeDuplicatedPeriodicZones__(a)

    # -------------------------------
    # Ecriture dans un fichier "elsA"
    # -------------------------------
    # Parcours des zones. On regarde si ce sont des zones donneuses. Si oui, on ecrit le(s) fichier(s) elsA correspondant(s).
    if writeInFile == 1:
        if depth == 2: EX=0
        elif depth == 1: EX=1
        Connector.writeCoefs(ntotZones, listRcvId,listIndicesRcv,listEXdir,listIndicesDonor,
                             listInterpolantsDonor,
                             listInterpTypes,listCellN,listZoneDim,nbInterpCellsForDonor,
                             prefixFile, isEX=EX, solver=solver, nGhostCells=nGhostCells)
    return None

#-----------------------------------------------------------------------------
# IN: zname: name of interpolated zone
# IN: zdonor: donor zone
# OUT: modified donor zone
#-----------------------------------------------------------------------------
def _interpInverseStorage(zname, zdonor, nozd, resInterp, depth):
    if depth == 1: loc = 'faces'
    else: loc = 'centers'
    indicesDonor=resInterp[1][nozd]
    indicesDonorExtC=resInterp[2][nozd]
    indicesOrphan = resInterp[7][nozd]
    if indicesDonor.size > 0 or indicesOrphan.size > 0:
        indRcv = resInterp[0][nozd]; interpCoef=resInterp[3][nozd];
        interpVol=resInterp[4][nozd]; interpType=resInterp[5][nozd]
        indicesExtrap = resInterp[6][nozd];
        if depth == 2: EXdir=numpy.array([], Internal.E_NpyInt)
        elif depth == 1: EXdir = resInterp[8][nozd]
        cellIndExtrap = indicesExtrap; cellIndOrphan = indicesOrphan
        coef=interpCoef; vol=interpVol; interptype=interpType
        # creation des noeuds ZoneSubRegion_t contenant les coefficients d'interpolation et les indices des cellules donneuses
        # Il y a un noeud ZoneSubRegion_t par paire de blocs interpole / bloc donneur
        XOD._createInterpRegion__(zdonor,zname, indicesDonor, indRcv, coef,
                                  interptype, vol, cellIndExtrap,cellIndOrphan,
                                  EXdir,indicesDonorExtC, 'Donor',loc=loc)
    return None

#-----------------------------------------------------------------------------
# IN: zdonorname: name of donor zone
# IN: z: interpolated zone
# OUT: modified interpolated zone
#-----------------------------------------------------------------------------
def _interpDirectStorage(z, zdonorname, nozd, resInterp, depth):
    if depth == 1: loc = 'faces'
    else: loc = 'centers'
    indicesDonor=resInterp[1][nozd]
    indicesDonorExtC=resInterp[2][nozd]
    indicesOrphan = resInterp[7][nozd]
    if indicesDonor.size > 0 or indicesOrphan.size >0:
        indRcv = resInterp[0][nozd]; interpCoef=resInterp[3][nozd]
        interpVol=resInterp[4][nozd]; interpType=resInterp[5][nozd]
        indicesExtrap = resInterp[6][nozd];
        if depth == 2: EXdir=numpy.array([], Internal.E_NpyInt)
        elif depth == 1: EXdir = resInterp[8][nozd]
        cellIndExtrap = indicesExtrap; cellIndOrphan = indicesOrphan
        coef = interpCoef; vol = interpVol; interptype = interpType
        # creation des noeuds ZoneSubRegion_t contenant les coefficients d'interpolation et
        # les indices des cellules donneuses
        # Il y a un noeud ZoneSubRegion_t par paire de blocs interpole / bloc donneur
        XOD._createInterpRegion__(z, zdonorname, indRcv, indicesDonor, coef,
                                  interptype, vol, cellIndExtrap, cellIndOrphan,
                                  EXdir,indicesDonorExtC, 'Receiver', loc=loc)
    return None



#=============================================================================
#=============================================================================
# 2. TRANSFERTS
#=============================================================================
#=============================================================================

#------------------------------------------------------------------------------
# Calcul des transferts Chimere
# IN: vars: liste des noms des champs sur lesquels le transfert est applique
#------------------------------------------------------------------------------
def chimeraTransfer(t, storage='inverse', variables=[], loc='cell',mesh='extended'):
    vars2 = []
    for v in variables:
        v2 = v.split(':',1)
        if len(v2) == 2 and v2[0] == 'centers': vars2.append(v)
        # DBG else: print('Warning: chimeraTransfer: only variables located at centers taken into account.')
    if vars2 == []:
        raise ValueError("chimeraTransfer: no variable to transfer.")

    if storage == 'direct': return directChimeraTransfer__(t, vars2, loc)
    else: return inverseChimeraTransfer__(t, vars2, locinterp=loc, mesh=mesh)

#------------------------------------------------------------------------------
# Calcul des transferts Chimere a partir d'un stockage direct
# Les champs sur lesquels on applique les transferts Chimere sont
# entres par l'utilisateur
# IN: vars: liste des noms des champs sur lesquels le transfert est applique
# locinterp = 'cell' ou 'face'
# Pour l'instant on suppose que les interpolations sont realisees aux centres
# donc les transferts aussi
# de STEPHANIE : cette fonction semble obsolete...
#------------------------------------------------------------------------------
def directChimeraTransfer__(t, variables, locinterp):
    print('Warning: directChimeraTransfer__: donor mesh is located at cell centers.')
    print('For elsA simulations: please use inverse storage.')

    loc = 'centers'; depth = 2
    if locinterp == 'face': depth = 1

    tp = Internal.copyRef(t)
    zones = Internal.getZones(tp); nzones = len(zones)
    for noz in range(nzones):
        z = zones[noz]
        parent1,d1 = Internal.getParentOfNode(tp,z)
        lRcvArrays = []
        for i in variables:
            f = C.getField(i,z)[0]
            if f != []: lRcvArrays.append(f)
        if lRcvArrays != []: rcvArray = Converter.addVars(lRcvArrays)
        # Les interpolations sont stockees dans les noeuds ZoneSubRegion_t
        subRegions2 = Internal.getNodesFromType1(z,'ZoneSubRegion_t')
        subRegions=[]
        for s in subRegions2:
            sname = s[0][0:2]
            if sname=='ID': subRegions.append(s)
        for s in subRegions:
            idn = Internal.getNodeFromName1(s,'InterpolantsDonor')
            if idn is not None: # la subRegion decrit des interpolations Chimere
                interpDonor = idn[1]
                if depth == 1:
                    cellDonorEX = Internal.getNodeFromName1(s,'FaceListDonor')[1]
                    interpDonorEX = Internal.getNodeFromName1(s,'FaceInterpolantsDonor')[1]
                    interpDonorEXType = Internal.getNodeFromName1(s,'FaceInterpolantsType')[1]
                interpDonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                cellRcv   = Internal.getNodeFromName1(s,'PointList')[1]
                cellDonor = Internal.getNodeFromName1(s,'PointListDonor')[1]
                zdonorname = Internal.getValue(s)
                zdonors = Internal.getNodesFromName2(tp, zdonorname)
                zdonor = Internal.getNodesFromType1(zdonors,'Zone_t')[0]
                lDonorArrays = []
                for i in variables:
                    f = C.getField(i,zdonor)[0]
                    if f !=[]: lDonorArrays.append(f)
                if lDonorArrays != []: donorArray = Converter.addVars(lDonorArrays)
                # dimensions de la grille d interpolation
                ni = Internal.getZoneDim(zdonor)[1]; nj = Internal.getZoneDim(zdonor)[2]; nk = Internal.getZoneDim(zdonor)[3]
                # dimensions de la grille interpolee
                nir = Internal.getZoneDim(z)[1]; njr = Internal.getZoneDim(z)[2]
                if loc == 'centers':
                    ni = ni-1; nj = nj-1; nk = nk-1; nir = nir-1; njr = njr-1
                rcvArray = Connector.chimeraTransfer(cellRcv, cellDonor, interpDonorType, interpDonor,
                                                     rcvArray, donorArray)
##                 # traitement des pts EX
##                 if depth == 1:
##                     rcvArray=Connector.chimeraTransfer(cellRcv, cellDonorEX, interpDonorEXType,interpDonorEX,
##                rcvArray,donorArray)
                C.setFields([rcvArray], z, loc)
                parent1[2][d1] = z
    return tp

#=============================================================================
# Retourne la liste des domaines d'interpolation et la liste des celln associes
# en centres pour une zone noz0 donnee d une base nob0 donnee a partir des
# autres bases. Si sameBase=1: recherche aussi sur sa base
# Si periodicite, il faut avoir duplique les zones periodiques avant
# car la recherche des zones periodiques se fait par 'TempPeriodicZone'
# Fonction appelee par setInterpolations
#=============================================================================
def getInterpolationDomains__(bases, nob0, noz0, zn0, surfacest=[],
                              sameBase=0, intersectionDict=None):
    surf = 0; surfaces = []; periodiczones = []
    if surfacest != []: surf = 1
    nbases = len(bases)
    names = []; zones = []; cellns = []
    for nob in range(nbases):
        if nob != nob0:
            nodes = Internal.getNodesFromType1(bases[nob], 'Zone_t')
            for noz in range(len(nodes)):
                zname = nodes[noz][0]; intersect = 1
                zn = C.getFields(Internal.__GridCoordinates__,nodes[noz])[0]
                if intersectionDict is not None:
                    if zname not in intersectionDict: intersect = 0

                if intersect == 1:
                    cellN = C.getField('centers:cellN', nodes[noz])[0]
                    donorName = nodes[noz][0]
                    dupzoneinfo = Internal.getNodeFromName1(nodes[noz],'TempPeriodicZone')
                    isper = 0
                    if dupzoneinfo is not None:
                        # periodicite par rotation ? On regarde le signe de la rotation
                        rotinfo = Internal.getNodeFromName(dupzoneinfo,'SignOfRotationAngle')
                        if rotinfo is not None:
                            sign = Internal.getValue(rotinfo)
                            if sign == 1: isper=2
                            else: isper=3
                        donorName = Internal.getValue(dupzoneinfo)
                    periodiczones.append(isper)
                    names.append(donorName); zones.append(zn); cellns.append(cellN)
                    if (surf == 1): surfaces.append(surfacest[nob][noz])
        else:
            if sameBase == 1:
                nodes = Internal.getNodesFromType1(bases[nob], 'Zone_t')
                for noz in range(len(nodes)):
                    if (noz != noz0):
                        zname = nodes[noz][0]; intersect = 1
                        zn = C.getFields(Internal.__GridCoordinates__,nodes[noz])[0]
                        if intersectionDict is not None:
                            if zname not in intersectionDict: intersect = 0

                        if intersect == 1:
                            cellN = C.getField('centers:cellN', nodes[noz])[0]
                            donorName = nodes[noz][0]
                            dupzoneinfo = Internal.getNodeFromName1(nodes[noz],'TempPeriodicZone')
                            if dupzoneinfo is not None:
                                donorName = Internal.getValue(dupzoneinfo)
                                periodiczones.append(1)
                            else: periodiczones.append(0)
                            names.append(donorName); zones.append(zn); cellns.append(cellN)
                            if surf == 1: surfaces.append(surfacest[nob][noz])
    return [names, zones, cellns, surfaces, periodiczones]

#------------------------------------------------------------------------------
# Calcul des transferts Chimere a partir d'un stockage inverse
# Les noms des champs sur lesquels on applique les transferts Chimere sont
# entres par l'utilisateur
# IN: variables: liste des noms des champs sur lesquels le transfert est applique
# IN: locinterp='cell' ou 'face', ne marche que pour cell actuellement
# IN: mesh='standard' or 'extended'
#------------------------------------------------------------------------------
def inverseChimeraTransfer__(t, variables, locinterp='centers', mesh='extended'):

    depth = 2 # valeur par defaut pour des interpolations 'centers'
    if locinterp == 'face': depth = 1

    tp = Internal.copyRef(t)
    # On regarde si t est une zone ou une liste de zones
    isZoneOrListZone = False
    if Internal.isStdNode(tp) == -1: # noeud
        if tp[3] == 'Zone_t': isZoneOrListZone = True; tp = [tp]
    elif Internal.isStdNode(tp) == 0: #liste
        for zt in tp:
            if zt[3] == 'Zone_t': isZoneOrListZone = True
            else: isZoneOrListZone = False; break

    if isZoneOrListZone: zones = tp
    else: zones = Internal.getZones(tp)

    # liste des zones locales (sur le processeur courant)
    localZonesName = []
    for z in zones: localZonesName.append(z[0])

    # Definition de dictionnaires pour les champs interpoles :
    # Cle : nom zone receveuse. Valeurs : champs interpole + indices pts receveurs + volume cellule donneuse
    rcvFields = {}

    if depth == 1: rcvLocalFields = {}

    # Parcours des zones locales (qui sont zones donneuses)
    nzones = len(zones)
    for noz in range(nzones):
        zdonor = zones[noz]
        # Recuperation des champs des zones donneuses pour les grandeurs a transferer
        # Ces champs sont stockes dans un array unique donorArray
        lDonorArrays = []
        for i in variables:
            f = C.getField(i,zdonor)[0]
            if f !=[]:
                # Passage du maillage donneur en centres etendus
                if mesh == 'extended': lDonorArrays.append(Converter.center2ExtCenter(f))
                else: lDonorArrays.append(f)
        if lDonorArrays != []:
            donorArray = Converter.addVars(lDonorArrays)
           # Parcours des ZoneSubRegion pour effectuer les transferts
           # (Les interpolations sont stockees dans les noeuds ZoneSubRegion_t)
            subRegions = Internal.getNodesFromType1(zdonor,'ZoneSubRegion_t')
            for s in subRegions:
                sname = Internal.getName(s)[0:2]
                if sname == 'ID': # la ZoneSubRegion decrit des interpolations Chimere
                    idn = Internal.getNodeFromName1(s,'InterpolantsDonor')
                    idEX = Internal.getNodeFromName1(s,'FaceInterpolantsDonor')
                    if idn is not None or idEX is not None: # la ZoneSubRegion decrit des interpolations Chimere
                        rcvZoneName = Internal.getValue(s) # nom de la zone interpolee
                        #A.  traitement locinterp = centers
                        if idn is not None and depth == 2:
                            # A.1 recuperation des champs pour l interpolation
                            if mesh == 'extended': cellDonorList = Internal.getNodeFromName1(s,'PointListExtC')
                            else: cellDonorList = Internal.getNodeFromName1(s,'PointList')
                            #if cellDonorList is not None and cellDonorList[1] != []: # la region traite des interpolations depth = 2
                            if cellDonorList is not None and cellDonorList[1].size != 0: # la region traite des interpolations depth = 2
                                cellDonor = cellDonorList[1]
                                cellRcv   = Internal.getNodeFromName1(s,'PointListDonor')[1]
                                interpDonor = idn[1]
                                cellVol   = Internal.getNodeFromName1(s,'InterpolantsVol')[1]
                                interpDonorType = Internal.getNodeFromName1(s,'InterpolantsType')[1]
                                # A.2 creation du array pour stocker la solution interpolee
                                if rcvZoneName not in localZonesName: # zone receveuse pas sur le processeur locale
                                    vars2 = ''
                                    for v in variables:
                                        v2 = v.split(':',1)[1]
                                        if variables.index(v) == len(variables)-1: vars2=vars2+v2
                                        else: vars2=vars2+v2+','
                                    cellRcvPara = numpy.arange(cellRcv.shape[0], dtype=Internal.E_NpyInt)
                                    rcvField = numpy.empty((donorArray[1].shape[0],cellRcv.shape[0]), order='F')
                                    rcvArray = [vars2,rcvField,cellRcv.shape[0],1,1] # tableau monodimensionnel
                                else: # zone receveuse sur le processeur locale
                                    cellRcvPara = cellRcv
                                    zs = Internal.getNodesFromName2(tp,rcvZoneName)
                                    for zr in zs:
                                        if zr[3] == 'Zone_t':
                                            z = zr; break
                                    lRcvArrays = []
                                    for i in variables:
                                        f = C.getField(i,z)[0]
                                        if f !=[]: lRcvArrays.append(f)
                                    if lRcvArrays != []: rcvArray = Converter.addVars(lRcvArrays)
                                    else:
                                        raise ValueError("inverseChimeraTransfer: field to interpolate not in zone %d."%rcvZoneName)
                                # A.3 interpolation Chimere
                                rcvField = Connector.chimeraTransfer(cellRcvPara, cellDonor, interpDonorType,
                                                                     interpDonor, rcvArray, donorArray)[1]
                                if rcvZoneName not in localZonesName:
                                    # redimensionnement du tableau au nombre de cellules interpolees
                                    rcvField2 = numpy.empty((donorArray[1].shape[0]+2,cellRcv.shape[0]), order='F')
                                    rcvField2[:rcvField.shape[0]]=rcvField
                                    rcvField2[donorArray[1].shape[0]]=cellRcv
                                    rcvField2[donorArray[1].shape[0]+1]=cellVol
                                    if rcvZoneName in rcvFields:
                                        rcvFields[rcvZoneName].append(rcvField2)
                                    else:
                                        rcvFields[rcvZoneName]=[rcvField2]
                                else: # solution mise directement dans l arbre local
                                    rcvArray[1]=rcvField
                                    C.setFields([rcvArray], z, locinterp)
                        # traitement locinterp = face
                        elif idEX != None and depth == 1:
                            if mesh == 'extended': cellDonorListEX = Internal.getNodeFromName1(s,'FaceListExtC')
                            else : cellDonorListEX = Internal.getNodeFromName1(s,'FaceList')
                            if cellDonorListEX is not None: # la region traite des interpolations depth = 1
                                cellDonorEX = cellDonorListEX[1]
                                if cellDonorEX.size != 0:
                                    cellRcvEX = Internal.getNodeFromName1(s,'FaceListDonor')[1]
                                    cellRcvEXdummy=numpy.arange(cellRcvEX.shape[0], dtype=Internal.E_NpyInt)
                                    EXdir = Internal.getNodeFromName1(s,'FaceDirection')[1]
                                    interpDonorEX = idEX[1]
                                    cellVolEX   = Internal.getNodeFromName1(s,'FaceInterpolantsVol')[1]
                                    interpDonorEXType = Internal.getNodeFromName1(s,'FaceInterpolantsType')[1]
                                    rcvArrayEX = Converter.array(donorArray[0],cellRcvEX.shape[0], 1, 1)
                                    rcvFieldEX = Connector.chimeraTransfer(cellRcvEXdummy,cellDonorEX,interpDonorEXType,
                                                                           interpDonorEX, rcvArrayEX, donorArray)[1]
                                    rcvField2 = numpy.empty((donorArray[1].shape[0]+3,cellRcvEX.shape[0]))
                                    rcvField2[:rcvFieldEX.shape[0]]=rcvFieldEX
                                    rcvField2[donorArray[1].shape[0]]=cellRcvEX
                                    rcvField2[donorArray[1].shape[0]+1]=cellVolEX
                                    rcvField2[donorArray[1].shape[0]+2]=EXdir
                                    if rcvZoneName not in localZonesName:
                                        if rcvZoneName in rcvFields:
                                            rcvFields[rcvZoneName].append(rcvField2)
                                        else:
                                            rcvFields[rcvZoneName]=[rcvField2]
                                    else: # solution locale
                                        if rcvZoneName in rcvLocalFields:
                                            rcvLocalFields[rcvZoneName].append(rcvField2)
                                        else:
                                            rcvLocalFields[rcvZoneName]=[rcvField2]
    if depth==1: return rcvLocalFields, rcvFields
    return tp, rcvFields

#=============================================================================
#=============================================================================
# 3. INFO CHIMERE EXTRAITES
#=============================================================================
#=============================================================================

#=============================================================================
# Get overset info : if a cell is interpolated, orphan, extrapolated
#                    or the quality of the donor cell
# 'orphan' (0,1), 'extrapolated'(sum |cf|), 'interpolated'(0,1)
# 'cellRatio' : max(volD/volR,volR/volD)
# 'donorAspect' : edgeRatio of donor cell
# Compliant with setInterpolations storage
#=============================================================================
def chimeraInfo(a, type='interpolated'):
    """Extract Overset information when computed with setInterpolations."""
    t = Internal.copyRef(a)
    _chimeraInfo(t,type=type)
    return t

def _chimeraInfo(t,type='interpolated'):
    if type == 'interpolated':
        C._initVars(t,'centers:interpolated', 0.)
        _chimeraNatureOfCells__(t,'interpolated')
    elif type == 'extrapolated':
        C._initVars(t,'centers:extrapolated', 0.)
        _chimeraNatureOfCells__(t,'extrapolated')
    elif type =='orphan':
        C._initVars(t,'centers:orphan', 0.)
        _chimeraNatureOfCells__(t,'orphan')
    elif type =='cellRatio':
        C._initVars(t,'centers:cellRatio', 0.)
        _chimeraCellRatio__(t)
    elif type == 'donorAspect':
        C._initVars(t,'centers:donorAspect', 0.)
        _chimeraDonorAspect__(t)
    else: raise ValueError("chimeraInfo: type is not valid.")
    return None

#===============================================================================
# Information on overset quality : ratio of volumes between donors and receivers
# Compliant with setInterpolations
#===============================================================================
def _chimeraCellRatio__(t):
    try: import Generator.PyTree as G
    except: raise ImportError("chimeraInfo: requires Generator module.")
    G._getVolumeMap(t)
    zones = Internal.getZones(t)

    # First pass : direct storage
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname=='ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        volR = C.getField('centers:vol',zr)[0][1]
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                ListVolD = Internal.getNodesFromName2(s,"InterpolantsVol")
                if ListVolD != []:
                    ListVolD = ListVolD[0][1]
                    ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                    nindI = ListRcv.size
                    if nindI > 0:
                        field = Converter.array('cellRatio',nindI,1,1)
                        ListIndex = []
                        for noind in range(nindI):
                            index = ListRcv[noind]
                            volr = volR[0,index]
                            vold = ListVolD[noind]
                            cr = max(volr/vold,vold/volr)
                            field[1][0,noind] = cr
                        C._setPartialFields(zr, [field], [ListRcv], loc='centers')

    # 2nd pass: inverse storage
    zones = Internal.getZones(t)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname=='ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(t,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if zr != []:
                    zr = zr[0]
                    volR = C.getField('centers:vol',zr)[0][1]
                    ListVolD = Internal.getNodesFromName2(s,"InterpolantsVol")
                    if ListVolD != []:
                        ListVolD = ListVolD[0][1]
                        ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                        nindI = ListRcv.size
                        if nindI > 0:
                            field = Converter.array('cellRatio',nindI,1,1)
                            for noind in range(nindI):
                                index = ListRcv[noind]
                                volr = volR[0,index]
                                vold = ListVolD[noind]
                                cr = max(volr/vold,vold/volr)
                                field[1][0,noind] = cr
                            C._setPartialFields(zr, [field], [ListRcv], loc='centers')

    C._rmVars(t,'centers:vol') # faut il la detruire ou non ? pas de test leger pour savoir si c'etait ds l arbre avant
    return None

#=============================================================================
# Information on overset quality : aspect ratio of donor cells
# Compliant with setInterpolations
#==============================================================================
def _chimeraDonorAspect__(t):
    try: import Generator.PyTree as G
    except: raise ImportError("chimeraInfo: requires Generator module.")
    G._getEdgeRatio(t)
    zones = Internal.getZones(t)

    # First pass : direct storage
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname=='ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                zdnrname = Internal.getValue(s)
                zdonors = Internal.getNodesFromName2(t,zdnrname)
                zd = Internal.getNodesFromType1(zdonors,'Zone_t')
                if zd == []: raise ValueError("chimeraInfo: donor zone %s not found."%zdnrname)
                else: zd = zd[0]
                #
                ListDnr = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                #
                ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                nindI = ListRcv.size
                #
                if nindI > 0:
                    edgeRatio = C.getField('centers:EdgeRatio',zd)[0][1]
                    field = Converter.array('donorAspect',nindI,1,1)
                    for noind in range(nindI):
                        indD = ListDnr[noind]
                        field[1][0,noind] = edgeRatio[0,indD]
                    C._setPartialFields(zr, [field], [ListRcv], loc='centers')

    # Inverse storage
    zones = Internal.getZones(t)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname=='ID':
                idn = Internal.getNodesFromName1(s, 'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        dimsD = Internal.getZoneDim(zd)
        nidc = max(dimsD[1]-1,1)
        njdc = max(dimsD[2]-1,1)
        nidcnjdc = nidc*njdc
        edgeRatio = C.getField('centers:EdgeRatio',zd)[0][1]
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(t,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if zr != []:
                    zr = zr[0]
                    ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                    nindI = ListRcv.size
                    if nindI > 0:
                        ListDnr = Internal.getNodesFromName1(s,'PointList')[0][1]
                        field = Converter.array('donorAspect',nindI,1,1)
                        for noind in range(nindI):
                            indD = ListDnr[noind]
                            field[1][0,noind] = edgeRatio[0,indD]
                        C._setPartialFields(zr, [field], [ListRcv],loc='centers')

    C._rmVars(t, 'centers:EdgeRatio')
    return None

# =============================================================================
# extract info for interpolated/extrapolated/orphan cells
# Compliant with setInterpolations
# =============================================================================
def _chimeraNatureOfCells__(t, nature):
    zones = Internal.getZones(t)
    # First pass : direct storage
    for zr in zones:
        subRegions2 = Internal.getNodesFromType1(zr,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname== 'ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Receiver': # direct storage ok
                if nature == 'interpolated':
                    ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                    if ListRcv is not None:
                        field = Converter.array('interpolated',ListRcv.size,1,1)
                        field = Converter.initVars(field,'interpolated',1.)
                        C._setPartialFields(zr, [field], [ListRcv], loc='centers')

                elif nature == 'extrapolated':
                    ListExtrap = Internal.getNodesFromName1(s,'ExtrapPointList')
                    if ListExtrap != []:
                        ListRcv = Internal.getNodesFromName1(s,'PointList')[0][1]
                        ListExtrap = ListExtrap[0][1]
                        DonorTypes = Internal.getNodesFromName1(s,'InterpolantsType')[0][1]
                        Coefs = Internal.getNodesFromName1(s,'InterpolantsDonor')[0][1]
                        # Somme des |coefs| : necessite ListRcv, ListExtrap, Coefs, DonorType
                        field = Connector.connector.getExtrapAbsCoefs(ListRcv, ListExtrap, DonorTypes, Coefs)
                        C._setPartialFields(zr, [field], [ListExtrap],loc='centers')
                elif nature == 'orphan':
                    orphans = Internal.getNodesFromName(s, 'OrphanPointList')
                    if orphans != []:
                        ListOrphan = orphans[0][1]
                        field = Converter.array('orphan',ListOrphan.size,1,1)
                        field = Converter.initVars(field,'orphan',1.)
                        C._setPartialFields(zr, [field], [ListOrphan],loc='centers')

    # 2nd pass: storage in donor zones
    zones = Internal.getZones(t)
    for zd in zones:
        subRegions2 = Internal.getNodesFromType1(zd,'ZoneSubRegion_t')
        subRegions = []
        for s in subRegions2:
            sname = s[0][0:2]
            if sname== 'ID':
                idn = Internal.getNodesFromName1(s,'InterpolantsDonor')
                if idn != []: # la subRegion decrit des interpolations
                    subRegions.append(s)
        subRegions2 = []
        for s in subRegions:
            zoneRole = Internal.getNodesFromName2(s,'ZoneRole')[0]
            zoneRole = Internal.getValue(zoneRole)
            if zoneRole == 'Donor': # inverse storage ok
                zrcvname = Internal.getValue(s)
                zreceivers = Internal.getNodesFromName2(t,zrcvname)
                zr = Internal.getNodesFromType1(zreceivers,'Zone_t')
                if zr != []:
                    zr = zr[0]
                    if nature == 'interpolated':
                        ListRcv = Internal.getNodesFromName1(s, 'PointListDonor')[0][1]
                        if ListRcv is not None:
                            field = Converter.array('interpolated',ListRcv.size,1,1)
                            field = Converter.initVars(field,'interpolated',1.)
                            C._setPartialFields(zr, [field], [ListRcv], loc='centers')
                    elif nature == 'extrapolated':
                        ListExtrap = Internal.getNodesFromName1(s,'ExtrapPointList')
                        if ListExtrap != []:
                            ListExtrap = ListExtrap[0][1]
                            DonorTypes = Internal.getNodesFromName1(s,'InterpolantsType')[0][1]
                            Coefs = Internal.getNodesFromName1(s,'InterpolantsDonor')[0][1]
                            ListRcv = Internal.getNodesFromName1(s,'PointListDonor')[0][1]
                            # Somme des |coefs| : necessite ListRcv, ListExtrap, Coefs, DonorType
                            field = Connector.connector.getExtrapAbsCoefs(ListRcv, ListExtrap, DonorTypes, Coefs)
                            C._setPartialFields(zr, [field], [ListExtrap], loc='centers')
                    elif nature == 'orphan':
                        orphans = Internal.getNodesFromName(s, 'OrphanPointList')
                        if orphans != []:
                            ListOrphan = orphans[0][1]
                            field = Converter.array('orphan',ListOrphan.size,1,1)
                            field = Converter.initVars(field,'orphan',1.)
                            C._setPartialFields(zr, [field], [ListOrphan], loc='centers')
    return None
