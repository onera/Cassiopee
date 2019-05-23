# Interface pour MPI
import Converter.Mpi as Cmpi
from . import PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
from . import connector
import RigidMotion.PyTree as RM
import numpy

try: range = xrange
except: pass

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
    zones = Internal.getZones(tl)
    #print 'Rank %d has %d zones.'%(Cmpi.rank, len(zones))
    tl = X.optimizeOverlap(tl, double_wall, priorities, intersectionsDict)
    tl = Cmpi.rmXZones(tl)
    return tl

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
                       Cs=0.3831337844872463, Ts=1.0):
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
                        Cs=0.3831337844872463, Ts=1.0):

    if procDict is None: procDict = Cmpi.getProcDict(aD)
    if graph is None: graph = Cmpi.computeGraph(aD, type=type)

    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    zonesD = Internal.getZones(aD)
    for zD in zonesD:
        infos = X.setInterpTransfersD(zD, variables=variables, cellNVariable=cellNVariable, variablesIBC=variablesIBC, 
                                      bcType=bcType, varType=varType, compact=compact, Gamma=Gamma, Cv=Cv, MuS=MuS, Cs=Cs, Ts=Ts)
        for n in infos:
            rcvName = n[0]
            proc = procDict[rcvName]
            if proc == Cmpi.rank:
                field = n[1]
                #print 'direct', Cmpi.rank, rcvName
                if field != []:
                    listIndices = n[2]
                    z = Internal.getNodeFromName2(aR, rcvName)
                    C._setPartialFields(z, [field], [listIndices], loc=n[3])
            else:
                rcvNode = procDict[rcvName]
                #print Cmpi.rank, 'envoie a ',rcvNode
                if rcvNode not in datas: datas[rcvNode] = [n]
                else: datas[rcvNode] += [n]
                #print datas
    # Envoie des numpys suivant le graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas:
        #print Cmpi.rank, 'recoit de',i, '->', len(rcvDatas[i])
        for n in rcvDatas[i]:
            rcvName = n[0]
            #print 'reception', Cmpi.rank, rcvName
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = Internal.getNodeFromName2(aR, rcvName)
                C._setPartialFields(z, [field], [listIndices], loc=n[3])
    return None

#===============================================================================
# __setInterpTransfers  version optimiser de _setInterpTransfers: arbre t et tc compact, moins de python + de C
#
# Warning: inverse storage!
# IN: zones: list zone receveurs
# IN: zoneD: list zone donneurs
# IN: type: ID: interpolation, IBCD: IBCs, ALLD: interp+IBCs
# IN: bcType  0: glissement
#             1: adherence
#             2: loi de paroi log
#             3: loi de paroi Musker
# IN: varType=1,2,3: variablesIBC define (ro,rou,rov,row,roE(,ronutilde)),(ro,u,v,w,t(,nutilde)),(ro,u,v,w,p(,nutilde))
# Adim: KCore.adim1 for Minf=0.1
#===============================================================================
def __setInterpTransfers(zones, zonesD, vars, param_int, param_real, type_transfert, nitrun,
                         nstep, nitmax, rk, exploc, num_passage, bcType=0, varType=1, compact=1,
                         graph=None, procDict=None,
                         Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08, Cs=0.3831337844872463, Ts=1.0):

    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    nbcomIBC    = param_int[1]
    shift_graph = nbcomIBC + param_int[2+nbcomIBC] + 2

    for comm_P2P in range(1,param_int[0]+1):
        pt_ech = param_int[comm_P2P + shift_graph]
        dest   = param_int[pt_ech]

        no_transfert = comm_P2P
        if dest == Cmpi.rank: #transfert intra_processus
            #print 'transfert local', type_transfert
            connector.___setInterpTransfers(zones, zonesD, vars, param_int, param_real, nitrun, varType, bcType, 
                                            type_transfert, no_transfert,  nstep, nitmax, rk, exploc, num_passage, Gamma,Cv,MuS,Cs,Ts)

        else:
            #print 'transfert global', type_transfert
            infos = connector.__setInterpTransfersD(zones, zonesD, vars, param_int, param_real, nitrun, varType, bcType, 
                                                    type_transfert, no_transfert,  nstep, nitmax, rk, exploc, num_passage, Gamma,Cv,MuS,Cs,Ts) 
            if infos != []:
               for n in infos:
                  rcvNode = dest
                  #print Cmpi.rank, 'envoie a ',rcvNode, ' le paquet : ', n
                  if rcvNode not in datas: datas[rcvNode] = [n]
                  else: datas[rcvNode] += [n]
                  #print datas
    
    # Envoie des numpys suivant le graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas:
        #if Cmpi.rank==0: print Cmpi.rank, 'recoit de',i, '->', len(rcvDatas[i])
        for n in rcvDatas[i]:
            rcvName = n[0]
            #if Cmpi.rank==0: print 'reception', Cmpi.rank, 'no zone', zones[ rcvName ][0]
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = zones[rcvName]
                C._setPartialFields(z, [field], [listIndices], loc='centers')
    return None

#---------------------------------------------------------------------------------------------------------
# Transferts instationnaires en parallele
# avec prise en compte du mouvement
# absFrame = True : les coordonnees de t sont deja dans le repere absolu en entree
#---------------------------------------------------------------------------------------------------------
def _transfer(t, tc, variables, graph, intersectionDict, dictOfADT, 
              dictOfNobOfRcvZones, dictOfNozOfRcvZones,
              dictOfNobOfDnrZones, dictOfNozOfDnrZones, 
              dictOfNobOfRcvZonesC, dictOfNozOfRcvZonesC, 
              time=0., absFrame=True, procDict=None, cellNName='cellN'):
    if procDict is None: procDict = Cmpi.getProcDict(tc)
    
    # dictionnaire des matrices de mouvement pour passer du repere relatif d une zone au repere absolu
    dictOfMotionMatR2A={}
    dictOfMotionMatA2R={}
    coordsD=[0.,0.,0.];  coordsC= [0.,0.,0.] # XAbs = coordsD + Mat*(XRel-coordsC)
    dictOfFields={}; dictOfIndices={}
    
    datas={}
    for z in Internal.getZones(t):
        zname = Internal.getName(z)
        if zname not in dictOfNobOfRcvZones: continue

        # coordonnees dans le repere absolu de la zone receptrice
        # on les recupere de zc pour eviter un node2center des coordonnees de z
        nobc = dictOfNobOfRcvZonesC[zname]
        nozc = dictOfNozOfRcvZonesC[zname]
        zc = tc[2][nobc][2][nozc]
        if zc[0] != zname:# check
            raise ValueError("_transfer: t and tc skeletons must be identical.")

        C._cpVars(z,'centers:'+cellNName, zc, cellNName)
        res = X.getInterpolatedPoints(zc,loc='nodes', cellNName=cellNName) 
        # print 'Zone %s du proc %d a interpoler'%(zname, Cmpi.rank)

        if res is not None: 
            # print 'Res not None : zone %s du proc %d a interpoler'%(zname, Cmpi.rank)

            indicesI, XI, YI, ZI = res
            # passage des coordonnees du recepteur dans le repere absolu
            # si mouvement gere par FastS -> les coordonnees dans z sont deja les coordonnees en absolu
            if not absFrame: 
                if zname in dictOfMotionMatR2A:
                    MatRel2AbsR=RM.getMotionMatrixForZone(z, time=time, F=None)
                    dictOfMotionMatR2A[zname]=MatRel2AbsR
                else:
                    MatRel2AbsR = dictOfMotionMatR2A[zname]
                RM._moveN([XI,YI,ZI],coordsD,coordsC,MatRel2AbsR)

            procR = procDict[zname]
            for znamed in intersectionDict[zname]:
                procD = procDict[znamed]
                if procD == Cmpi.rank:
                    nobc = dictOfNobOfDnrZones[znamed]
                    nozc = dictOfNozOfDnrZones[znamed]
                    zdnr = tc[2][nobc][2][nozc]
                    adt = dictOfADT[znamed]
                    if znamed in dictOfMotionMatA2R:
                        MatAbs2RelD=dictOfMotionMatA2R[znamed]
                    else:                        
                        if znamed in dictOfMotionMatR2A:
                            MatRel2AbsD = dictOfMotionMatR2A[znamed]
                            MatAbs2RelD = numpy.transpose(MatRel2AbsD)
                            dictOfMotionMatA2R[znamed] = MatAbs2RelD
                        else:
                            MatRel2AbsD=RM.getMotionMatrixForZone(zdnr, time=time, F=None)
                            dictOfMotionMatR2A[znamed]=MatRel2AbsD
                            MatAbs2RelD = numpy.transpose(MatRel2AbsD)
                            dictOfMotionMatA2R[znamed] = MatAbs2RelD
                    [XIRel, YIRel, ZIRel] = RM.moveN([XI,YI,ZI],coordsC,coordsD,MatAbs2RelD)

                    # transfers avec coordonnees dans le repere relatif 
                    fields = X.transferFields(zdnr, XIRel, YIRel, ZIRel, hook=adt, variables=variables)
                    if zname not in dictOfFields:
                        dictOfFields[zname]=[fields]
                        dictOfIndices[zname]=indicesI
                    else:
                        dictOfFields[zname].append(fields)

                else:                    
                    # print ' ECHANGE GLOBAL entre recepteur %s du proc %d et donneur %s du proc %d '%(zname, Cmpi.rank, znamed, procD)
                    if procD not in datas:
                        datas[procD] = [[zname, znamed, indicesI, XI, YI, ZI]]
                    else: datas[procD].append([zname, znamed, indicesI, XI, YI, ZI])

    # print 'Proc  : ', Cmpi.rank, ' envoie les donnees : ' ,datas.keys()
    # print ' a partir du graphe ', graph
    # 1er envoi : envoi des numpys des donnees a interpoler suivant le graphe
    interpDatas = Cmpi.sendRecv(datas,graph)

    # recuperation par le proc donneur des donnees pour faire les transferts    
    transferedDatas={}
    for i in interpDatas:
        #print Cmpi.rank, 'recoit de',i, '->', len(interpDatas[i])
        for n in interpDatas[i]:
            zdnrname = n[1]
            zrcvname = n[0]
            indicesR = n[2]
            XI = n[3]; YI = n[4]; ZI = n[5]
            nobc = dictOfNobOfDnrZones[zdnrname]
            nozc = dictOfNozOfDnrZones[zdnrname]
            zdnr = tc[2][nobc][2][nozc]
            adt = dictOfADT[zdnrname]
            if zdnrname in dictOfMotionMatA2R:
                MatAbs2RelD=dictOfMotionMatA2R[zdnrname]
            else:
                if zdnrname in dictOfMotionMatR2A:
                    MatRel2AbsD = dictOfMotionMatR2A[zdnrname]
                    MatAbs2RelD = numpy.transpose(MatRel2AbsD)
                    dictOfMotionMatA2R[zdnrname] = MatAbs2RelD
                else:
                    MatRel2AbsD=RM.getMotionMatrixForZone(zdnr, time=time, F=None)
                    dictOfMotionMatR2A[zdnrname]=MatRel2AbsD
                    MatAbs2RelD = numpy.transpose(MatRel2AbsD)
                    dictOfMotionMatA2R[zdnrname] = MatAbs2RelD
            
            [XIRel, YIRel, ZIRel] = RM.moveN([XI,YI,ZI],coordsC,coordsD,MatAbs2RelD)
            # transferts avec coordonnees dans le repere relatif 
            fields = X.transferFields(zdnr, XIRel, YIRel, ZIRel, hook=adt, variables=variables)
            procR = procDict[zrcvname]
            
            if procR not in transferedDatas:
                transferedDatas[procR]=[[zrcvname, indicesR, fields]]
            else:
                transferedDatas[procR].append([zrcvname,indicesR,fields])
            
    if transferedDatas != {}:
        # 2nd envoi : envoi des numpys des donnees  interpolees suivant le graphe
        rcvDatas = Cmpi.sendRecv(transferedDatas,graph)
                                                              
        # remise des donnees interpolees chez les zones receveuses
        # une fois que tous les donneurs potentiels ont calcule et envoye leurs donnees
        for i in rcvDatas:
            #print Cmpi.rank, 'recoit des donnees interpolees de',i, '->', len(rcvDatas[i])
            for n in rcvDatas[i]:
                zrcvname = n[0]
                indicesI = n[1]
                fields = n[2]
                if zrcvname not in dictOfFields:
                    dictOfFields[zrcvname]=[fields]
                    dictOfIndices[zrcvname]=indicesI
                else:
                    dictOfFields[zrcvname].append(fields)

    for zrcvname in dictOfIndices:
        nob = dictOfNobOfRcvZones[zrcvname]
        noz = dictOfNozOfRcvZones[zrcvname]
        z = t[2][nob][2][noz]
        allInterpFields = dictOfFields[zrcvname]
        indicesI = dictOfIndices[zrcvname]
        C._filterPartialFields(z, allInterpFields, indicesI, loc='centers', startFrom=0, filterName='donorVol')

    # SORTIE
    return None
