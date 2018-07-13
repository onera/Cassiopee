# Interface pour MPI
import Converter.Mpi as Cmpi
import PyTree as X
import Converter.Internal as Internal
import Converter.PyTree as C
import connector

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
#             3: loi de paroi Musker
# IN: varType=1,2,3: variablesIBC define (ro,rou,rov,row,roE(,ronutilde)),(ro,u,v,w,t(,nutilde)),(ro,u,v,w,p(,nutilde))
# Adim: KCore.adim1 for Minf=0.1
#===============================================================================
def setInterpTransfers(aR, aD, variables=[], 
                       variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'], 
                       bcType=0, varType=1, graph=None, 
                       procDict=None, type='ALLD', 
                       Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08, 
                       Cs=0.3831337844872463, Ts=1.0):
    tp = Internal.copyRef(aR)
    compact = 0
    _setInterpTransfers(tp, aD, variables, variablesIBC, 
                        bcType, varType,  compact, graph, 
                        procDict, type, Gamma, Cv, MuS, Cs, Ts)
    return tp
#===============================================================================
def _setInterpTransfers(aR, aD, variables=[], 
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
        infos = X.setInterpTransfersD(zD, variables, variablesIBC, bcType, varType, compact, Gamma, Cv, MuS, Cs, Ts)
        for n in infos:
            rcvName = n[0]
            proc = procDict[rcvName]
            if proc == Cmpi.rank:
                field = n[1]
                #print 'direct', Cmpi.rank, rcvName
                if field != []:
                    listIndices = n[2]
                    z = Internal.getNodesFromName2(aR, rcvName)[0]
                    C._setPartialFields(z, [field], [listIndices], loc=n[3])
            else:
                rcvNode = procDict[rcvName]
                #print Cmpi.rank, 'envoie a ',rcvNode
                if not datas.has_key(rcvNode): datas[rcvNode] = [n]
                else: datas[rcvNode] += [n]
                #print datas
    
    # Envoie des numpys suivant le graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas.keys():
        #print Cmpi.rank, 'recoit de',i, '->', len(rcvDatas[i])
        for n in rcvDatas[i]:
            rcvName = n[0]
            #print 'reception', Cmpi.rank, rcvName
            field = n[1]
            if field != []:
                listIndices = n[2]
                z = Internal.getNodesFromName2(aR, rcvName)[0]
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
                         bcType=0, varType=1, compact=1, graph=None, 
                         procDict=None,
                         Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08, Cs=0.3831337844872463, Ts=1.0):

    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    nbcomIBC    = param_int[1]
    shift_graph = nbcomIBC + param_int[2+nbcomIBC] + 2

    for comm_P2P in xrange(1,param_int[0]+1):
        pt_ech = param_int[comm_P2P + shift_graph]
        dest   = param_int[pt_ech]

        no_transfert = comm_P2P
        if dest == Cmpi.rank: #transfert intra_processus
            # print 'transfert local', no_transfert
            connector.___setInterpTransfers(zones, zonesD, vars, param_int, param_real, nitrun, varType, bcType, 
                                            type_transfert, no_transfert, Gamma,Cv,MuS,Cs,Ts)

        else:
            # print 'transfert global'
            infos = connector.__setInterpTransfersD(zones, zonesD, vars, param_int, param_real, nitrun, varType, bcType, 
                                                    type_transfert, no_transfert,Gamma,Cv,MuS,Cs,Ts) 
 
            for n in infos:
                rcvNode = dest
                #print Cmpi.rank, 'envoie a ',rcvNode, ' le paquet : ', n
                if not datas.has_key(rcvNode): datas[rcvNode] = [n]
                else: datas[rcvNode] += [n]
                #print datas
    
    # Envoie des numpys suivant le graph
    rcvDatas = Cmpi.sendRecv(datas, graph)

    # Remise des champs interpoles dans l'arbre receveur
    for i in rcvDatas.keys():
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

