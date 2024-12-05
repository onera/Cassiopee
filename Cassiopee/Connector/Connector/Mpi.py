# Interface pour MPI
import Converter.Mpi as Cmpi
from . import PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import Converter.converter
import numpy
from . import connector
import RigidMotion.PyTree as RM

#==============================================================================
# optimizeOverlap
# IN: t: full/loaded skel/partial
# IN: graph: graph d'intersection si deja calcule
# IN: intersectionsDict: Dictionnaire d'intersections
# OUT: arbre partiel avec overlap optimise
#==============================================================================
def optimizeOverlap(t, double_wall=0, priorities=[], graph=None,
                    intersectionsDict=None):
    if graph is None:
        tb = Cmpi.createBBoxTree(t)
        graph = Cmpi.computeGraph(tb, type='bbox2',
                                  intersectionsDict=intersectionsDict)
    tl = Cmpi.addXZones(t, graph)
    tl = Cmpi.convert2PartialTree(tl)
    # print info
    # zones = Internal.getZones(tl)
    # print('Rank %d has %d zones.'%(Cmpi.rank, len(zones)))
    tl = X.optimizeOverlap(tl, double_wall, priorities, intersectionsDict)
    tl = Cmpi.rmXZones(tl)
    return tl

#==============================================================================
# connectMatch
#==============================================================================
def connectMatch(a, tol=1.e-6, dim=3):

    # Ajout des bandelettes
    Cmpi._addBXZones(a, depth=2)

    # Construction des raccords
    a = X.connectMatch(a, tol=tol, dim=dim)

    # Suppression des XZones et correction des matchs
    Cmpi._rmBXZones(a)

    # Fusion des fenetres des raccords
    a = mergeWindows(a)

    return a

#==============================================================================
def connectNearMatch(a, ratio=2, tol=1.e-6, dim=3):
    """Find boundaries that matches with a given ratio."""
    if not isinstance(ratio, list):
        iratio = ratio
    else:
        iratio = 1
        for r in ratio: iratio=max(iratio,r)

    # Ajout des bandelettes
    Cmpi._addBXZones(a, depth=iratio+1)

    # Construction des raccords
    a = X.connectNearMatch(a, ratio=2, tol=tol, dim=dim)

    for z in Internal.getZones(a):
        gcs = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
        for g in gcs:
            nodes = Internal.getNodesFromType1(g, 'GridConnectivity_t')
            for n in nodes:
                gctype = Internal.getNodeFromType(n,'GridConnectivityType_t')
                if Internal.getValue(gctype)=='Abutting':
                    nmratio = Internal.getNodeFromName(n,'NMRatio')
                    nmratio = Internal.getValue(nmratio)
                    fratio = 1.
                    for i in nmratio: fratio *= i

                    if fratio==1.:
                        Internal._rmNodesByName(z,n[0])

    # Suppression des XZones et correction des matchs
    Cmpi._rmBXZones(a)

    # Fusion des fenetres des raccords
    a = mergeWindows(a)

    return a

#==============================================================================
# connectMatchPeriodic
#==============================================================================
def connectMatchPeriodic(a, rotationCenter=[0.,0.,0.],
                         rotationAngle=[0.,0.,0.],
                         translation=[0.,0.,0.], tol=1.e-6, dim=3,
                         unitAngle=None):

    # Ajout des bandelettes
    Cmpi._addBXZones(a, depth=2,allB=True)

    # Construction des raccords
    a = X.connectMatchPeriodic(a,rotationCenter,rotationAngle,translation,tol,dim,unitAngle)

    # Suppression des XZones et correction des matchs
    Cmpi._rmBXZones(a)

    # Fusion des fenetres des raccords
    a = mergeWindows(a)

    return a

#==============================================================================
def giveName2Window(p, zname, zopp):
    if p[0] == p[1]:
        if p[0] == 1: pos = zname+'_imin_'+zopp
        else: pos = zname+'_imax_'+zopp

    elif p[2] == p[3]:
        if p[2] == 1: pos = zname+'_jmin_'+zopp
        else: pos = zname+'_jmax_'+zopp

    elif p[4] == p[5]:
        if p[4] == 1: pos = zname+'_kmin_'+zopp
        else: pos = zname+'_kmax_'+zopp

    return pos

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def mergeWindows(t):
    # Merge grid connectivities created after addBXZones
    zones = Internal.getZones(t)
    for z in zones:
        xz = Internal.getNodeFromName1(z, 'XZone')
        if xz is None:
            # Construction du dictionnaire des matchs
            dico = {}
            gcs   = Internal.getNodesFromType1(z, 'ZoneGridConnectivity_t')
            for g in gcs:
                nodes = Internal.getNodesFromType1(g, 'GridConnectivity1to1_t')

                for n in Internal.getNodesFromType1(g,'GridConnectivity_t'):
                    gctype = Internal.getNodeFromType(n,'GridConnectivityType_t')
                    if Internal.getValue(gctype)=='Abutting': nodes.append(n)

                for n in nodes:
                    pr    = Internal.getNodeFromName1(n, 'PointRange')
                    p     = Internal.range2Window(pr[1])
                    zopp  = Internal.getValue(n)
                    pos   = giveName2Window(p,z[0],zopp)

                    if pos not in dico.keys(): dico[pos] = [n[0]]
                    else: dico[pos].append(n[0])

            # Test si match peuvent etre fusionnes
            for match in dico.keys():
                if len(dico[match]) > 1:
                    sumSurf = 0
                    pglob   = [None]*6
                    for name in dico[match]:
                        node    = Internal.getNodeFromName(z,name)
                        pr      = Internal.getNodeFromName1(node, 'PointRange')
                        p       = Internal.range2Window(pr[1])
                        surf    = max(1, p[1]-p[0])*max(1, p[3]-p[2])*max(1, p[5]-p[4])
                        sumSurf = sumSurf + surf

                        if pglob[0] is None:
                            pglob[0] = p[0] ; pglob[1] = p[1]
                            pglob[2] = p[2] ; pglob[3] = p[3]
                            pglob[4] = p[4] ; pglob[5] = p[5]
                        else:
                            if pglob[0] > p[0] : pglob[0] = p[0]
                            if pglob[1] < p[1] : pglob[1] = p[1]
                            if pglob[2] > p[2] : pglob[2] = p[2]
                            if pglob[3] < p[3] : pglob[3] = p[3]
                            if pglob[4] > p[4] : pglob[4] = p[4]
                            if pglob[5] < p[5] : pglob[5] = p[5]

                    surfMatch = max(1,(pglob[1]-pglob[0]))*max(1,(pglob[3]-pglob[2]))*max(1,(pglob[5]-pglob[4]))

                    # Fusion des matchs
                    if surfMatch == sumSurf:
                        # Fenetre du match donneur
                        pglobD   = [None]*6
                        for name in dico[match]:
                            node    = Internal.getNodeFromName(z,name)
                            prd     = Internal.getNodeFromName2(node, 'PointRangeDonor')
                            pd      = Internal.range2Window(prd[1])
                            #print(" name = %s, p = "%(name),p, prd)

                            if pglobD[0] is None:
                                pglobD[0] = pd[0] ; pglobD[1] = pd[1]
                                pglobD[2] = pd[2] ; pglobD[3] = pd[3]
                                pglobD[4] = pd[4] ; pglobD[5] = pd[5]
                            else:
                                if pglobD[0] > pd[0] : pglobD[0] = pd[0]
                                if pglobD[1] < pd[1] : pglobD[1] = pd[1]
                                if pglobD[2] > pd[2] : pglobD[2] = pd[2]
                                if pglobD[3] < pd[3] : pglobD[3] = pd[3]
                                if pglobD[4] > pd[4] : pglobD[4] = pd[4]
                                if pglobD[5] < pd[5] : pglobD[5] = pd[5]

                        # Modif du 1er match et suppression des autres
                        first = True
                        for name in dico[match]:
                            if first:
                                first = False
                                modifMatch = dico[match][0]
                                node    = Internal.getNodeFromName(z,modifMatch)
                                pr      = Internal.getNodeFromName1(node, 'PointRange')
                                prd     = Internal.getNodeFromName2(node, 'PointRangeDonor')
                                pglob   = Internal.window2Range(pglob)
                                pglobD  = Internal.window2Range(pglobD)
                                Internal.setValue(pr,  pglob)
                                Internal.setValue(prd, pglobD)
                            else:
                                Internal._rmNodesByName(z, name)

                    else:
                        print("Warning: mergeWindows: in zone ",z[0], " fail to merge matches: ", dico[match])

    return t

#===============================================================================
# setInterpTransfers
# Warning: inverse storage!
# IN: aR: arbre des receveurs
# IN: aD: arbre des donneurs
# IN: type: ID: interpolation, IBCD: IBCs, ALLD: interp+IBCs
# IN: bcType  0: glissement
#             1: adherence
#             2: loi de paroi log
#             3: loi de paroi Musker,4: outpress, 5 inj, 6 TBLE-SA
# IN: varType=1,2,3: variablesIBC define (ro,rou,rov,row,roE(,ronutilde)),(ro,u,v,w,t(,nutilde)),(ro,u,v,w,p(,nutilde))
# Adim: KCore.adim1 for Minf=0.1
#===============================================================================
def setInterpTransfers(aR, aD, variables=[], cellNVariable='',
                       variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                       bcType=0, varType=1, graph=None,
                       procDict=None, type='ALLD',
                       Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                       Cs=0.3831337844872463, Ts=1.0, alpha=1.):
    tp = Internal.copyRef(aR)
    compact = 0
    _setInterpTransfers(tp, aD, variables=variables, cellNVariable=cellNVariable, variablesIBC=variablesIBC,
                        bcType=bcType, varType=varType,  compact=compact, graph=graph,
                        procDict=procDict, type=type, Gamma=Gamma, Cv=Cv, MuS=MuS, Cs=Cs, Ts=Ts)
    return tp
#===============================================================================
def _setInterpTransfers(aR, aD, variables=[], cellNVariable='',
                        variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'],
                        bcType=0, varType=1, compact=0, graph=None,
                        procDict=None, type='ALLD',
                        Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08,
                        Cs=0.3831337844872463, Ts=1.0, alpha=1.):

    if procDict is None: procDict = Cmpi.getProcDict(aD)
    if graph is None: graph = Cmpi.computeGraph(aD, type=type)

    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    zonesD = Internal.getZones(aD)
    for zD in zonesD:
        infos = X.setInterpTransfersD(zD, variables=variables, cellNVariable=cellNVariable, variablesIBC=variablesIBC,
                                      bcType=bcType, varType=varType, compact=compact, Gamma=Gamma, Cv=Cv, MuS=MuS,
                                      Cs=Cs, Ts=Ts, alpha=alpha)
        for n in infos:
            rcvName = n[0]
            proc = procDict[rcvName]
            if proc == Cmpi.rank:
                field = n[1]
                #print('direct', Cmpi.rank, rcvName)
                if field != []:
                    listIndices = n[2]
                    z = Internal.getNodeFromName2(aR, rcvName)
                    C._setPartialFields(z, [field], [listIndices], loc=n[3])
            else:
                rcvNode = procDict[rcvName]
                #print(Cmpi.rank, 'envoie a ',rcvNode)
                if rcvNode not in datas: datas[rcvNode] = [n]
                else: datas[rcvNode] += [n]
                #print datas
    # Envoie des numpys suivant le graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas:
        #print(Cmpi.rank, 'recoit de',i, '->', len(rcvDatas[i]))
        for n in rcvDatas[i]:
            rcvName = n[0]
            #print('reception', Cmpi.rank, rcvName)
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = Internal.getNodeFromName2(aR, rcvName)
                C._setPartialFields(z, [field], [listIndices], loc=n[3])
    return None

#===============================================================================
# setInterpTransfers for pressure gradients information (compact = 0)
#===============================================================================
def _setInterpTransfersForPressureGradients(aR, aD, ibctypes=[], secondOrder=False, procDict=None, graph=None):

    if procDict is None: procDict = Cmpi.getProcDict(aD)
    if graph is None: graph = Cmpi.computeGraph(aD, type=type)

    datas = {}
    zonesD = Internal.getZones(aD)
    for zD in zonesD:
        infos = X._setIBCTransfersDForPressureGradients(aD, ibctypes=ibctypes, secondOrder=secondOrder)
        for n in infos:
            rcvName = n[0]
            proc = procDict[rcvName]
            if proc == Cmpi.rank:
                field = n[1]
                if field != []:
                    listIndices = n[2]
                    z = Internal.getNodeFromName2(aR, rcvName)
                    C._setPartialFields(z, [field], [listIndices], loc=n[3])
            else:
                rcvNode = procDict[rcvName]
                if rcvNode not in datas: datas[rcvNode] = [n]
                else: datas[rcvNode] += [n]
    # Envoie des numpys suivant le graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas:
        for n in rcvDatas[i]:
            rcvName = n[0]
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = Internal.getNodeFromName2(aR, rcvName)
                C._setPartialFields(z, [field], [listIndices], loc=n[3])
    return None

#===============================================================================
# __setInterpTransfers - version optimisee de _setInterpTransfers: arbre t et tc compact,
# moins de python + de C
#
# Warning: inverse storage!
# IN: zones: list zones receveurs
# IN: zoneD: list zones donneurs
# IN: type: ID: interpolation, IBCD: IBCs, ALLD: interp+IBCs
# IN: bcType  0: glissement
#             1: adherence
#             2: loi de paroi log
#             3: loi de paroi Musker
# IN: varType=1,2,3: variablesIBC define (ro,rou,rov,row,roE(,ronutilde)),(ro,u,v,w,t(,nutilde)),(ro,u,v,w,p(,nutilde))
# Adim: KCore.adim1 for Minf=0.1
#===============================================================================
def __setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, type_transfert, it_target,
                         nstep, nitmax, rk, exploc, num_passage, varType=1, compact=1,
                         graph=None, procDict=None, isWireModel_int=0):

    ##isWireModel_int: either -1, 0 , or 1
    ## 1:: YES :wire model treatment
    ## 0:: NO  :wire model treatment
    ##-1:: NO  :wire model treatment BUT a treatment on locks for IBC for ___setInterpTransfers
    isWireModel_intv2 = max(isWireModel_int,0)
    isSetPartialFieldsCheck = max(abs(isWireModel_int),0)

    ##for moving IBMs
    isIbmMoving_int  = 0
    #modif Ivan: souci si rigidext et ibm fixe: empeche le comportememnt normal de impli-local
    #motionType = int(Internal.getNodeFromName(zones, "Parameter_real")[1][64])
    #motionType = Cmpi.allreduce(motionType, op=Cmpi.MAX)
    #if motionType==3: isIbmMoving_int=1

    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    nbcomIBC    = param_int[2]
    shift_graph = nbcomIBC + param_int[3+nbcomIBC] + 3

    for comm_P2P in range(1,param_int[1]+1):
        pt_ech = param_int[comm_P2P + shift_graph]
        dest   = param_int[pt_ech]

        no_transfert = comm_P2P
        if dest == Cmpi.rank: #transfert intra_processus
            connector.___setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, it_target, varType,
                                            type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage,
                                            isWireModel_intv2)

        else:
            rank  = Cmpi.rank
            infos = connector.__setInterpTransfersD(zones, zonesD, vars, dtloc, param_int, param_real, it_target, varType,
                                                    type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage, rank,
                                                    isWireModel_int, isIbmMoving_int)
            if infos != []:
                for n in infos:
                    rcvNode = dest
                    if rcvNode not in datas: datas[rcvNode] = [n]
                    else: datas[rcvNode] += [n]

##[AJ] Keep for Now
##if isWireModel_int==0 and isSetPartialFieldsCheck==1:
##for i in datas:
##    for n in datas[i]:
##        rcvName   = n[0]
##        fieldname = n[1][0]
##        field = n[1]
##        #print(fieldname,flush=True)
##        count=0
##        dd=fieldname.split(',')
##        for i in dd:
##            print(i,field[1][count],'rank=',Cmpi.rank,flush=True)
##            count+=1

    # Envoie des numpys suivant le graph
    if graph is not None:
        rcvDatas = Cmpi.sendRecvC(datas, graph)
        #rcvDatas = Cmpi.sendRecv(datas, graph)
    else: rcvDatas = {}

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas:
        #if Cmpi.rank==0: print(Cmpi.rank, 'recoit de',i, '->', len(rcvDatas[i]), 'nstep=',nstep,flush=True)
        for n in rcvDatas[i]:
            rcvName = n[0]
            field = n[1]

            isSetPartialFields = True
            if isSetPartialFieldsCheck==1 and field != []:
                minfld = numpy.ndarray.min(field[1][0])
                maxfld = numpy.ndarray.max(field[1][0])
                if (maxfld == minfld and maxfld < -1e05): isSetPartialFields=False

            if isSetPartialFields:
                listIndices = n[2]
                z = zones[rcvName]
                C._setPartialFields(z, [field], [listIndices], loc='centers')

    return None

#===============================================================================
# Copie de __setInterpTransfers() pour la maj des info de gradient de pression
# transfert de ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature', 'nuSA'] -> 6 vars
# transfert de ['gradx/y/zDensity', 'gradx/y/zTemperature'] -> 6 varsGrad
# varType 22 : tc2/tc -> RCV ZONES
# varType 23 : RCV ZONES -> tc
#===============================================================================
def __setInterpTransfers4GradP(zones, zonesD, vars, param_int, param_real, type_transfert, it_target,
                               nstep, nitmax, rk, exploc, num_passage, varType=1, compact=1,
                               graph=None, procDict=None):
    isWireModelPrep= False
    isWireModel    = False
    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    datasGradP = {}
    nbcomIBC    = param_int[2]
    shift_graph = nbcomIBC + param_int[3+nbcomIBC] + 3

    for comm_P2P in range(1,param_int[1]+1):
        pt_ech = param_int[comm_P2P + shift_graph]
        dest   = param_int[pt_ech]

        no_transfert = comm_P2P
        if dest == Cmpi.rank: #transfert intra_processus
            connector.___setInterpTransfers4GradP(zones, zonesD, vars, param_int, param_real, it_target, varType,
                                                  type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage,
                                                  isWireModel)

        else:
            if varType != 24:
                allInfos = connector.__setInterpTransfersD4GradP(zones, zonesD, vars, param_int, param_real, it_target, varType,
                                                                 type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage)
                infos = allInfos[0]
                infosGrad = allInfos[1]

                if infos != []:
                    for n in infos:
                        rcvNode = dest
                        if rcvNode not in datas: datas[rcvNode] = [n]
                        else: datas[rcvNode] += [n]


                if infosGrad != []:
                    for n in infosGrad:
                        rcvNode = dest
                        if rcvNode not in datasGradP: datasGradP[rcvNode] = [n]
                        else: datasGradP[rcvNode] += [n]

    # Envoie des numpys suivant le graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas:
        for n in rcvDatas[i]:
            rcvName = n[0]
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = zones[rcvName]
                C._setPartialFields(z, [field], [listIndices], loc='centers')

    # Envoie des numpys suivant le graph
    rcvDatas = Cmpi.sendRecv(datasGradP, graph)

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas:
        for n in rcvDatas[i]:
            rcvName = n[0]
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = zones[rcvName]
                C._setPartialFields(z, [field], [listIndices], loc='centers')
    return None

#---------------------------------------------------------------------------------------------------------
# Transferts instationnaires en parallele
# avec prise en compte du mouvement
# absFrame=True: les coordonnees de t sont deja dans le repere absolu en entree
# interpInDnrFrame=True: interpolation avec les coordonnees des pts a interpoler dans le repere relatif au donneur
# applicable en mouvement rigide; en mvt avec deformation: mettre False
# verbose: 0 (rien), 1 (bilan interpolation), 2 (ecrit les indices de pts orphelins),
# 3 (met cellN=-1 pour les points orphelins)
#---------------------------------------------------------------------------------------------------------
def _transfer2(t, tc, variables, graph, intersectionDict, dictOfADT,
               dictOfNobOfRcvZones, dictOfNozOfRcvZones,
               dictOfNobOfDnrZones, dictOfNozOfDnrZones,
               dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC,
               time=0., absFrame=True, procDict=None, cellNName='cellN',
               interpInDnrFrame=True, order=2,
               hook=None, verbose=1):

    if procDict is None: procDict = Cmpi.getProcDict(tc)
    if Cmpi.size == 1:
        for name in procDict: procDict[name]=0
    # dictionnaire des matrices de mouvement pour passer du repere relatif d'une zone au repere absolu
    dictOfMotionMatR2A={}; dictOfMotionMatA2R={}
    coordsD=[0.,0.,0.]; coordsC=[0.,0.,0.] # XAbs = coordsD + coordsC + Mat*(XRel-coordsC)
    dictOfFields={}; dictOfIndices={}

    # 1. Formation data interpolation globale
    #Cmpi.trace("1. transfer2 - start")
    datas={}; listOfLocalData = []; interpDatas={}

    if hook is not None and len(hook) == 2:
        listOfLocalData = hook[0]; interpDatas = hook[1]
    else: # empty hook
        for z in Internal.getZones(t):
            zname = Internal.getName(z)
            if zname not in dictOfNobOfRcvZones: continue

            # coordonnees dans le repere absolu de la zone receptrice
            # on les recupere de zc pour eviter un node2center des coordonnees de z
            nobc = dictOfNobOfRcvZonesC[zname]  #no base
            nozc = dictOfNozOfRcvZonesC[zname]  #no zone
            zc = tc[2][nobc][2][nozc]
            if zc[0] != zname: # check
                raise ValueError("_transfer: t and tc skeletons must be identical.")

            C._cpVars(z, 'centers:'+cellNName, zc, cellNName)
            res = X.getInterpolatedPoints(zc, loc='nodes', cellNName=cellNName)
            if res is not None:
                indicesI, XI, YI, ZI = res  #indiceI des pts cellN=2 et coord des pts
                # passage des coordonnees du recepteur dans le repere absolu
                # si mouvement gere par FastS -> les coordonnees dans z sont deja les coordonnees en absolu
                if not absFrame:
                    if zname in dictOfMotionMatR2A:
                        MatRel2AbsR = RM.getMotionMatrixForZone(z, time=time, F=None)
                        dictOfMotionMatR2A[zname]=MatRel2AbsR
                    else:
                        MatRel2AbsR = dictOfMotionMatR2A[zname]
                    RM._moveN([XI,YI,ZI],coordsD,coordsC,MatRel2AbsR)

                procR = procDict[zname]
                for znamed in intersectionDict[zname]:
                    procD = procDict[znamed]
                    if procD == Cmpi.rank: # local
                        # local delayed
                        listOfLocalData.append([zname, znamed, indicesI, XI, YI, ZI])
                    else:
                        if procD not in datas: datas[procD] = [[zname, znamed, indicesI, XI, YI, ZI]]
                        else: datas[procD].append([zname, znamed, indicesI, XI, YI, ZI])

    # 2. envoie data interpolation globale en asynchrone
        #Cmpi.trace("2. transfer2")
        reqs = []
        if graph != {}:
            if Cmpi.rank in graph:
                g = graph[Cmpi.rank] # graph du proc courant
                for oppNode in g:
                    # Envoie les datas necessaires au noeud oppose
                    #print('%d: On envoie a %d: %s'%(rank,oppNode,g[oppNode]))
                    if oppNode in datas:
                        s = Converter.converter.iSend(datas[oppNode], oppNode, Cmpi.rank, Cmpi.KCOMM)
                    else:
                        s = Converter.converter.iSend(None, oppNode, Cmpi.rank, Cmpi.KCOMM)
                    reqs.append(s)

    # 3. interpolation locale
    #Cmpi.trace("3. transfer2")
    for z in listOfLocalData:
        zname   = z[0]
        znamed  = z[1]
        indicesI= z[2]
        XI      = z[3]
        YI      = z[4]
        ZI      = z[5]
        nobc = dictOfNobOfDnrZones[znamed]
        nozc = dictOfNozOfDnrZones[znamed]
        zdnr = tc[2][nobc][2][nozc]
        adt = dictOfADT[znamed]
        if adt is None: interpDataType = 0
        else: interpDataType = 1
        if interpInDnrFrame: [XIRel,YIRel,ZIRel] = RM.evalPositionM1([XI,YI,ZI], zdnr, time)
        else: [XIRel,YIRel,ZIRel] = [XI,YI,ZI]

        # transfers avec coordonnees dans le repere relatif
        if interpInDnrFrame and Internal.getNodeFromName1(zdnr, 'TimeMotion') is not None:
            # On suppose que le precond est dans init quand il y a un TimeMotion
            # Si il y en a pas, on suppose que le precond est construit dans courant
            GC1 = Internal.getNodeFromName1(zdnr, 'GridCoordinates')
            GC2 = Internal.getNodeFromName1(zdnr, 'GridCoordinates#Init')
            TEMP = GC1[2]; GC1[2] = GC2[2]; GC2[2] = TEMP

        fields = X.transferFields(zdnr, XIRel, YIRel, ZIRel, order=order, hook=adt, variables=variables, interpDataType=interpDataType)

        # hack par CB
        if interpInDnrFrame and Internal.getNodeFromName1(zdnr, 'TimeMotion') is not None:
            TEMP = GC1[2]; GC1[2] = GC2[2]; GC2[2] = TEMP

        if zname not in dictOfFields:
            dictOfFields[zname] = [fields]
            dictOfIndices[zname] = indicesI
        else:
            dictOfFields[zname].append(fields)

    # 4. reception des donnees d'interpolation globales
    #Cmpi.trace("4. transfer2")
    if hook is not None and len(hook) == 0:
        if graph != {}:
            for node in graph:
                if Cmpi.rank in graph[node]:
                    rec = Converter.converter.recv(node, Cmpi.rank, Cmpi.KCOMM)
                    if rec is not None: interpDatas[node] = rec
        a = Converter.converter.waitAll(reqs)

    # set hook
    if hook is not None and len(hook) == 0: hook += [listOfLocalData, interpDatas]

    # 5. interpolation globales
    #Cmpi.trace("5. transfer2")
    transferedDatas={}
    for i in interpDatas:
        for n in interpDatas[i]:
            zdnrname = n[1]
            zrcvname = n[0]
            indicesR = n[2]
            XI = n[3]; YI = n[4]; ZI = n[5]
            nobc = dictOfNobOfDnrZones[zdnrname]
            nozc = dictOfNozOfDnrZones[zdnrname]
            zdnr = tc[2][nobc][2][nozc]
            adt = dictOfADT[zdnrname]
            if adt is None: interpDataType = 0
            else: interpDataType = 1
            if interpInDnrFrame: [XIRel,YIRel,ZIRel] = RM.evalPositionM1([XI,YI,ZI], zdnr, time)
            else: [XIRel,YIRel,ZIRel] = [XI,YI,ZI]

            # [XIRel,YIRel,ZIRel] = RM.moveN([XI,YI,ZI],coordsC,coordsD,MatAbs2RelD)
            # transferts avec coordonnees dans le repere relatif
            if interpInDnrFrame and Internal.getNodeFromName1(zdnr, 'TimeMotion') is not None:
                GC1 = Internal.getNodeFromName1(zdnr, 'GridCoordinates')
                GC2 = Internal.getNodeFromName1(zdnr, 'GridCoordinates#Init')
                TEMP = GC1[2]; GC1[2] = GC2[2]; GC2[2] = TEMP

            fields = X.transferFields(zdnr, XIRel, YIRel, ZIRel, hook=adt, variables=variables, interpDataType=interpDataType)

            # hack par CB
            if interpInDnrFrame and Internal.getNodeFromName1(zdnr, 'TimeMotion') is not None:
                TEMP = GC1[2]; GC1[2] = GC2[2]; GC2[2] = TEMP

            procR = procDict[zrcvname]
            if procR not in transferedDatas:
                transferedDatas[procR]=[[zrcvname,indicesR,fields]]
            else:
                transferedDatas[procR].append([zrcvname,indicesR,fields])

    # 6. envoie des numpys des donnees interpolees suivant le graphe
    #Cmpi.trace("6. transfer2")
    rcvDatas = Cmpi.sendRecvC(transferedDatas, graph)
    #rcvDatas = Cmpi.sendRecv(transferedDatas, graph)

    # 7. remise des donnees interpolees chez les zones receveuses
    # une fois que tous les donneurs potentiels on calcule et envoye leurs donnees
    #Cmpi.trace("7. transfer2")
    for i in rcvDatas:
        for n in rcvDatas[i]:
            zrcvname = n[0]
            indicesI = n[1]
            fields = n[2]
            if zrcvname not in dictOfFields:
                dictOfFields[zrcvname] = [fields]
                dictOfIndices[zrcvname] = indicesI
            else:
                dictOfFields[zrcvname].append(fields)

    for zrcvname in dictOfIndices:
        nob = dictOfNobOfRcvZones[zrcvname]
        noz = dictOfNozOfRcvZones[zrcvname]
        z = t[2][nob][2][noz]
        allInterpFields = dictOfFields[zrcvname]
        indicesI = dictOfIndices[zrcvname]
        C._filterPartialFields(z, allInterpFields, indicesI, loc='centers', startFrom=0,
                               filterName='donorVol', verbose=verbose)

    #Cmpi.trace("8. transfer2 end")
    return None

#=========================================================================
# partie delicate :
# IN: t, tc: arbres partiels locaux
# IN: sameBase=1 (itype='chimera'): autorise l'interpolation dans la meme base
# memes arguments que setInterpData
#=========================================================================
def _setInterpData(aR, aD, order=2, penalty=1, nature=0, extrap=1,
                   method='lagrangian', loc='nodes', storage='direct',
                   interpDataType=1, hook=None, cartesian=False, sameBase=0,
                   topTreeRcv=None, topTreeDnr=None, sameName=1, verbose=2,
                   dim=3, itype='abutting'):
    """Compute interpolation data for abutting or chimera intergrid connectivity."""

    # create dictOfModels for adaptRANSLES in OversetData
    dictOfModels = {}
    for b in Internal.getBases(aR):
        model_b = Internal.getNodeFromName2(b, 'GoverningEquations')
        if model_b is not None: model_b = Internal.getValue(model_b)
        else: model_b = 'None'
        for z in Internal.getZones(b):
            model = Internal.getNodeFromName2(z, 'GoverningEquations')
            if model is None: model = model_b
            else: model = Internal.getValue(model)
            dictOfModels[z[0]] = [model]

    dictOfModels = Cmpi.allgatherDict(dictOfModels)
    dictOfModels = {key:value[0] for key,value in dictOfModels.items()}

    # Le graph doit correspondre au probleme
    if itype == 'abutting':
        graph = Cmpi.computeGraph(aR, type='match', reduction=True)
        Cmpi._addXZones(aR, graph, variables=[], noCoordinates=True,
                        cartesian=False, zoneGC=True, keepOldNodes=False)
        Cmpi._addXZones(aD, graph, variables=[], noCoordinates=True,
                        cartesian=False, zoneGC=True, keepOldNodes=False)
        X._setInterpData(aR, aD, order=order, penalty=penalty, nature=nature, extrap=extrap,
                         method=method, loc=loc, storage=storage, interpDataType=interpDataType, hook=hook,
                         topTreeRcv=topTreeRcv, topTreeDnr=topTreeDnr,
                         sameName=sameName, dim=dim, itype=itype, dictOfModels=dictOfModels)
        Cmpi._rmXZones(aR); Cmpi._rmXZones(aD)

    elif itype == 'chimeraOld': # ancienne version
        tbbc = Cmpi.createBBoxTree(aD)
        interDict = X.getIntersectingDomains(tbbc)
        # on ne conserve que les intersections inter base
        baseNames = {}
        for b in Internal.getBases(tbbc):
            for z in Internal.getZones(b): baseNames[z[0]] = b[0]
        for i in interDict:
            bi = baseNames[i]
            out = []
            for z in interDict[i]:
                if bi != baseNames[z]: out.append(z)
            interDict[i] = out

        graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
        Cmpi._addXZones(aR, graph, variables=['centers:cellN'], noCoordinates=False,
                        cartesian=False, zoneGC=False, keepOldNodes=False)
        Cmpi._addXZones(aD, graph, variables=['centers:cellN'], noCoordinates=False,
                        cartesian=False, zoneGC=False, keepOldNodes=False)

        X._setInterpData(aR, aD, order=order, penalty=penalty, nature=nature, extrap=extrap,
                         method=method, loc=loc, storage=storage, interpDataType=interpDataType, hook=hook,
                         topTreeRcv=topTreeRcv, topTreeDnr=topTreeDnr,
                         sameName=sameName, dim=dim, itype=itype, dictOfModels=dictOfModels)

        Cmpi._rmXZones(aR); Cmpi._rmXZones(aD)

    elif itype == 'chimera': # nouvelle version
        tbbc = Cmpi.createBBoxTree(aD)
        interDict = X.getIntersectingDomains(tbbc)
        procDict = Cmpi.getProcDict(aD)

        # Get baseName for each zone
        baseNames = {}
        for b in Internal.getBases(tbbc):
            for z in Internal.getZones(b): baseNames[z[0]] = b[0]

        # on ne conserve que les intersections inter bases
        if sameBase == 0:
            for i in interDict:
                bi = baseNames[i]
                out = []
                for z in interDict[i]:
                    if bi != baseNames[z]: out.append(z)
                interDict[i] = out

        # Perform addXZones on aD
        graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
        Cmpi._addXZones(aD, graph, variables=['cellN'], noCoordinates=False,
                        cartesian=cartesian, zoneGC=False, keepOldNodes=False)

        # serialisation eventuelle
        #graphs = Cmpi.splitGraph(graph)
        #for g in graphs:
        #    Cmpi._addXZones(aD, g, variables=['centers:cellN'], noCoordinates=False,
        #                    cartesian=False, zoneGC=False, keepOldNodes=False)

        # Build hook on local aD zones
        hooks = {};
        for b in Internal.getBases(aD):
            if b[0] == 'CARTESIAN':
                for z in Internal.getZones(b):
                    hooks[z[0]] = None # must be None for Cartesian
            else:
                for z in Internal.getZones(b):
                    hooks[z[0]] = C.createHook(z, 'adt')

        datas = {}
        # InterpData par zone
        for zr in Internal.getZones(aR):
            zrname = Internal.getName(zr)
            baseNameRcv = baseNames[zrname]
            dnrZones = []
            for zdname in interDict[zrname]:
                zd = Internal.getNodeFromName2(aD, zdname)
                baseNameDnr = baseNames[zd[0]]
                if sameBase == 0:
                    if baseNameDnr != baseNameRcv: dnrZones.append(zd)
                else:
                    dnrZones.append(zd)

            hookL = []; interpDataTypeL = []
            for z in dnrZones:
                h = hooks[z[0]]
                hookL.append(h)
                if h is None: interpDataTypeL.append(0)
                else: interpDataTypeL.append(1)

            if dnrZones != []:
                X._setInterpData(zr, dnrZones, order=order, penalty=penalty,
                                 nature=nature, extrap=extrap, verbose=verbose,
                                 method=method, loc=loc, storage=storage,
                                 interpDataType=interpDataTypeL, hook=hookL,
                                 topTreeRcv=topTreeRcv, topTreeDnr=topTreeDnr,
                                 sameName=sameName, dim=dim, itype="chimera", dictOfModels=dictOfModels)

            for zd in dnrZones:
                zdname = zd[0]
                destProc = procDict[zdname]

                IDs = []
                for i in zd[2]:
                    if i[0][0:2] == 'ID':
                        if Internal.getValue(i) == zrname: IDs.append(i)

                if IDs != []:
                    if destProc == Cmpi.rank:
                        zD = Internal.getNodeFromName2(aD, zdname)
                        zD[2] += IDs
                    else:
                        if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                        else: datas[destProc].append([zdname,IDs])
                else:
                    if destProc not in datas: datas[destProc] = []

        Cmpi._rmXZones(aD)
        for h in hooks: C.freeHook(hooks[h])

        destDatas = Cmpi.sendRecv(datas, graph)
        for i in destDatas:
            for n in destDatas[i]:
                zdname = n[0]
                IDs = n[1]
                if IDs != []:
                    zD = Internal.getNodeFromName2(aD, zdname)
                    zD[2] += IDs
        datas = {}; destDatas = None

    return None

def setInterpData2(tR, tD, order=2, loc='centers', cartesian=False):
    """Compute interpolation data for 2 different trees."""
    aD = Internal.copyRef(tD)
    aR = Internal.copyRef(tR)
    _setInterpData2(aR, aD, order=order, loc=loc, cartesian=cartesian)
    return aD

def _setInterpData2(tR, tD, order=2, loc='centers', cartesian=False):
    """Compute interpolation data for 2 different trees."""

    if loc == 'nodes': varcelln = 'cellN'
    else: varcelln = 'centers:cellN'

    # Clean previous IDs if necessary
    Internal._rmNodesFromType(tD, 'ZoneSubRegion_t')
    Internal._rmNodesFromName(tD, 'GridCoordinates#Init')

    if cartesian: interpDataType = 0 # 0 if tc is cartesian
    else: interpDataType = 1
    locR = loc
    # Compute BBoxTrees
    tsBB = Cmpi.createBBoxTree(tR)
    procDicts = Cmpi.getProcDict(tsBB)
    tDBB = Cmpi.createBBoxTree(tD)
    procDictD = Cmpi.getProcDict(tDBB)
    interDicts = X.getIntersectingDomains(tsBB, tDBB)
    interDictD2R = X.getIntersectingDomains(tDBB, tsBB)

    graph = Cmpi.computeGraph(tDBB, type='bbox3', intersectionsDict=interDictD2R,
                              procDict=procDictD, procDict2=procDicts, t2=tsBB, reduction=True)
    graph2 = Cmpi.computeGraph(tsBB, type='bbox3', intersectionsDict=interDicts,
                               procDict=procDicts, procDict2=procDictD, t2=tDBB, reduction=True)
    Cmpi._addXZones(tD, graph, variables=['cellN'], cartesian=cartesian, subr=False, keepOldNodes=False)

    datas = {}
    for zs in Internal.getZones(tR):
        zrname = Internal.getName(zs)
        dnrZones = []
        for zdname in interDicts[zrname]:
            zd = Internal.getNodeFromName2(tD, zdname)
            dnrZones.append(zd)

        cellNPresent = C.isNamePresent(zs, varcelln)
        if cellNPresent==-1: C._initVars(zs, varcelln, 2.) # interp all

        if dnrZones != []:
            X._setInterpData(zs, dnrZones, nature=1, penalty=1, order=order, loc=locR, storage='inverse',
                             sameName=0, interpDataType=interpDataType, itype='chimera')
        if cellNPresent == -1:
            C._rmVars(zs, [varcelln])
        for zd in dnrZones:
            zdname = zd[0]
            destProc = procDictD[zdname]

            IDs = []
            for i in zd[2]:
                if i[0][0:2] == 'ID':
                    if Internal.getValue(i) == zrname: IDs.append(i)

            if IDs != []:
                if destProc == Cmpi.rank:
                    zD = Internal.getNodeFromName2(tD, zdname)
                    zD[2] += IDs
                else:
                    if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                    else: datas[destProc].append([zdname,IDs])
            else:
                if destProc not in datas: datas[destProc] = []

    Cmpi._rmXZones(tD)
    destDatas = Cmpi.sendRecv(datas, graph2)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tD, zname)
                zD[2] += IDs
    datas = {}; destDatas = None

    return None


#==============================================================================
def __setInterpTransfers_WireModel(zones, zonesD, vars, dtloc, param_int, param_real, type_transfert, nitrun,
                                   nstep, nitmax, rk, exploc, num_passage, varType=1, compact=1,
                                   graph=None, procDict=None, graphIBCD=None, graphInvIBCD_WM=None, nvars=5):

    variablesIBC=['Density_WM', 'VelocityX_WM', 'VelocityY_WM', 'VelocityZ_WM', 'Temperature_WM', 'TurbulentSANuTilde_WM']
    if nvars == 5: variablesIBC=['Density_WM', 'VelocityX_WM', 'VelocityY_WM', 'VelocityZ_WM', 'Temperature_WM']

    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    datasGradP = {}
    nbcomIBC    = param_int[2]
    shift_graph = nbcomIBC + param_int[3+nbcomIBC] + 3

    for comm_P2P in range(1,param_int[1]+1):
        pt_ech = param_int[comm_P2P + shift_graph]
        dest   = param_int[pt_ech]

        no_transfert = comm_P2P
        if dest == Cmpi.rank: #transfert intra_processus
            isWireModel_int=2
            connector.___setInterpTransfers(zones, zonesD, vars, dtloc, param_int, param_real, nitrun, varType,
                                            type_transfert, no_transfert, nstep, nitmax, rk, exploc, num_passage,
                                            isWireModel_int)

    datas = {}
    znr   = {}
    for z in zones: znr[z[0]] = z

    for zd in zonesD:
        subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
        for s in subRegions:
            sname = s[0].split('_')[1]
            zname = s[0].split('_')[-1]
            dest = procDict[zname]
            if sname == '140' and dest != Cmpi.rank:
                ListDonor  = numpy.copy(Internal.getNodeFromName1(s, 'PointList')[1])
                ListRcv    = numpy.copy(Internal.getNodeFromName1(s, 'PointListDonor')[1])
                dens_wm    = numpy.copy(Internal.getNodeFromName1(s, 'Density_WM')[1])
                velx_wm    = numpy.copy(Internal.getNodeFromName1(s, 'VelocityX_WM')[1])
                vely_wm    = numpy.copy(Internal.getNodeFromName1(s, 'VelocityY_WM')[1])
                velz_wm    = numpy.copy(Internal.getNodeFromName1(s, 'VelocityZ_WM')[1])
                temp_wm    = numpy.copy(Internal.getNodeFromName1(s, 'Temperature_WM')[1])
                sanu_wm    = numpy.copy(Internal.getNodeFromName1(s, 'TurbulentSANuTilde_WM')[1])
                infos = [zd[0], zname, ListDonor, ListRcv, dens_wm, velx_wm, vely_wm, velz_wm, temp_wm, sanu_wm]
                rcvNode = dest
                if rcvNode not in datas: datas[rcvNode] = [infos]
                else: datas[rcvNode] += [infos]

    rcvDatas = Cmpi.sendRecv(datas, graphIBCD)
    datas    = {}
    for dest in rcvDatas:
        for [name, zname, ListDonor, ListRcv, dens_wm, velx_wm, vely_wm, velz_wm, temp_wm, sanu_wm] in rcvDatas[dest]:
            zr = znr[zname]
            connector._WM_getVal2tc(zr, variablesIBC, ListRcv,
                                    dens_wm, velx_wm, vely_wm, velz_wm, temp_wm, sanu_wm,
                                    1, nvars,
                                    Internal.__GridCoordinates__,
                                    Internal.__FlowSolutionNodes__,
                                    Internal.__FlowSolutionCenters__)
            rcvNode = dest
            infos = [name, zname, dens_wm, velx_wm, vely_wm, velz_wm, temp_wm, sanu_wm]
            if rcvNode not in datas: datas[rcvNode] = [infos]
            else: datas[rcvNode] += [infos]

    rcvDatas = Cmpi.sendRecv(datas, graphInvIBCD_WM)
    for dest in rcvDatas:
        for [name, zname, dens_wm_new, velx_wm_new, vely_wm_new, velz_wm_new, temp_wm_new, sanu_wm_new] in rcvDatas[dest]:
            for zd in zonesD:
                if zd[0] == name:
                    subRegions = Internal.getNodesFromType1(zd, 'ZoneSubRegion_t')
                    for s in subRegions:
                        sname = s[0].split('_')[1]
                        znameD = s[0].split('_')[-1]
                        if sname == '140' and znameD == zname:
                            ListRcv   = Internal.getNodeFromName1(s, 'PointListDonor')[1]
                            dens_wm    = Internal.getNodeFromName1(s, 'Density_WM')[1]
                            velx_wm    = Internal.getNodeFromName1(s, 'VelocityX_WM')[1]
                            vely_wm    = Internal.getNodeFromName1(s, 'VelocityY_WM')[1]
                            velz_wm    = Internal.getNodeFromName1(s, 'VelocityZ_WM')[1]
                            temp_wm    = Internal.getNodeFromName1(s, 'Temperature_WM')[1]
                            sanu_wm    = Internal.getNodeFromName1(s, 'TurbulentSANuTilde_WM')[1]
                            connector._WM_setVal2tc(dens_wm_new, velx_wm_new, vely_wm_new, velz_wm_new, temp_wm_new, sanu_wm_new,
                                                    dens_wm    , velx_wm    , vely_wm    , velz_wm    , temp_wm    , sanu_wm    )
    return None
