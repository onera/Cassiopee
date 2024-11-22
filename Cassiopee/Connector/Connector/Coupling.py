# Fonctions basees sur Connector pour remplacer Pmpi.extractMesh
# pour faire du couplage CFD/externe ou du body force
# 2 fonctionnements possibles :
# update direct de la solution par interpolation
# ajout d'une contribution d'un donneur a un champ existant - plusieurs donneurs possibles
import Converter.Mpi as Cmpi 
import Converter.Internal as Internal 
import Converter.Filter as Filter
import Converter.PyTree as C
import Geom.PyTree as D
import Generator.PyTree as G
import Converter.Distributed as Distributed
import Distributor2.PyTree as D2
from . import PyTree as X
from . import Mpi as Xmpi
import numpy, os

# typeTransfer=0 : mise a jour du champ directement de l interpolant
# typeTransfer=1 : ajout des contributions des interpolants
def _interpolate(tR, tD, interpTree, graph, procDict, typeTransfer=0):
    FSC_SAV = Internal.__FlowSolutionCenters__
    FSN_SAV = Internal.__FlowSolutionNodes__

    ZSR = Internal.getNodeFromType(interpTree,"ZoneSubRegion_t")
    if ZSR is None: return None
    GL = Internal.getNodeFromType(ZSR,'GridLocation_t')
    if GL is not None:
        GL = Internal.getValue(GL)
        if GL == 'Vertex': loc='nodes'
        elif GL =='CellCenter': loc = 'centers'
    else:
        loc = 'nodes'

    zd = Internal.getNodeFromType(tD,"Zone_t")
    FSconts = Internal.getNodesFromType(zd,'FlowSolution_t')
    dictOfFSC={}; dictOfFSN={}
    for fs in FSconts:
        GL = Internal.getNodeFromType(fs,'GridLocation_t')
        if GL is not None:
            GL = Internal.getValue(GL)
            if GL == 'Vertex': locI='nodes'
            elif GL =='CellCenter': locI = 'centers'
        else:
            locI = 'nodes'

        varsI = []
        if locI==loc:
            if loc=='centers':
                fsname = Internal.getName(fs)
                for fnode in Internal.getNodesFromType(fs,'DataArray_t'):
                    fname = Internal.getName(fnode)
                    varsI.append(fname)
                dictOfFSC[fsname] = varsI

            else:
                fsname = Internal.getName(fs)
                for fnode in Internal.getNodesFromType(fs,'DataArray_t'):
                    fname = Internal.getName(fnode)            
                    varsI.append(fname)
                dictOfFSN[fsname] = varsI

    for fsname in dictOfFSC:
        Internal.__FlowSolutionCenters__=fsname
        varsI = dictOfFSC[fsname]
        for varl in varsI:
            C._cpVars(tD,'centers:%s'%varl, interpTree, varl)
            C._initVars(tR,'centers:%s'%varl,0)
        _setInterpTransfers(tR, interpTree, variables=varsI, cellNVariable='cellN',
                            graph=graph, procDict=procDict, type='ID', typeTransfer=typeTransfer)
    for fsname in dictOfFSN:
        Internal.__FlowSolutionNodes__=fsname
        varsI = dictOfFSN[fsname]
        for varl in varsI:
            C._cpVars(tD,'%s'%varl, interpTree, varl)
            C._initVars(tR,'%s'%varl,0)
        _setInterpTransfers(tR, interpTree, variables=varsI, cellNVariable='cellN',
                            graph=graph, procDict=procDict, type='ID', typeTransfer=typeTransfer)  

    Internal.__FlowSolutionCenters__ = FSC_SAV
    Internal.__FlowSolutionNodes__ = FSN_SAV
    return None

def prepareInterpData(tR, tD, order=2, loc='CellCenter', cartesian=False, cleanID=True, typeTransfer=0):
    if loc=='CellCenter':
        locR = 'centers'
        tc = C.node2Center(tD)
    else:
        locR = 'nodes'
        tc = Internal.copyRef(tD)

    [graphR, procDictR]=_setInterpData2__(tR, tc, order=2, loc=locR, cartesian=cartesian, cleanID=cleanID,
                                          typeTransfer=typeTransfer)
    Internal._rmNodesFromType(tc,'FlowSolution_t')
    Internal._rmNodesFromType(tc,'GridCoordinates_t')
    return [tc, graphR, procDictR]

def _setInterpData2__(tR, tD, order=2, loc='centers', cartesian=False, cleanID=True, typeTransfer=0):
    if loc == 'nodes': varcelln = 'cellN'
    else: varcelln = 'centers:cellN'    

    # Clean previous IDs if necessary
    if cleanID:
        Internal._rmNodesFromType(tD, 'ZoneSubRegion_t')
        Internal._rmNodesFromName(tD, 'GridCoordinates#Init')

    if cartesian: interpDataType = 0 # 0 if tc is cartesian
    else: interpDataType = 1
    locR = loc

    # Compute BBoxTrees
    tRBB = Cmpi.createBBoxTree(tR)
    procDictR = Cmpi.getProcDict(tRBB)
    tDBB = Cmpi.createBBoxTree(tD)
    procDictD = Cmpi.getProcDict(tDBB)
    interDictR2D = X.getIntersectingDomains(tRBB, tDBB)
    interDictD2R = X.getIntersectingDomains(tDBB, tRBB)

    graphR = Cmpi.computeGraph(tDBB, type='bbox3', intersectionsDict=interDictD2R,
                               procDict=procDictD, procDict2=procDictR, t2=tRBB, reduction=True)
    graphD = Cmpi.computeGraph(tRBB, type='bbox3', intersectionsDict=interDictR2D,
                               procDict=procDictR, procDict2=procDictD, t2=tDBB, reduction=True)
    Cmpi._addXZones(tD, graphR, variables=['cellN'], cartesian=cartesian, subr=False, keepOldNodes=False)

    datas = {}
    for zs in Internal.getZones(tR):
        zdim = Internal.getZoneDim(zs)
        zrname = Internal.getName(zs)
        dnrZones = []
        for zdname in interDictR2D[zrname]:
            zd = Internal.getNodeFromName2(tD, zdname)
            dnrZones.append(zd)

        cellNPresent = C.isNamePresent(zs, varcelln)
        if cellNPresent==-1: C._initVars(zs, varcelln, 2.) # interp all
        if dnrZones != []:
            if typeTransfer == 0:
                X._setInterpData(zs, dnrZones, nature=1, penalty=1, order=order, loc=loc, 
                                 storage='inverse', extrap=0, verbose=0,
                                 sameName=0, interpDataType=interpDataType, itype='chimera')
            else:
                for zd in dnrZones:
                    X._setInterpData(zs, zd, nature=1, penalty=1, order=order, loc=loc, 
                                     storage='inverse', extrap=0, verbose=0,
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
    destDatas = Cmpi.sendRecv(datas, graphD)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tD, zname)
                zD[2] += IDs
    datas = {}; destDatas = None

    # clean IDs
    for zD in Internal.getZones(tD):
        dicto={}
        for zsr in Internal.getNodesFromType(zD,"ZoneSubRegion_t"):
            idname = zsr[0]
            if idname not in dicto: dicto[idname] = 1
            else:
                Internal._rmNode(tD, zsr)
    return [graphR, procDictR]

#===============================================================================
# typeTransfer=0 : replace field by interpolated value
# typeTransfer=1 : sum to existing field the interpolated value 
def _setInterpTransfers(aR, aD, variables=[], cellNVariable='',
                        variablesIBC=['Density','MomentumX','MomentumY','MomentumZ','EnergyStagnationDensity'], 
                        bcType=0, varType=1, compact=0, graph=None, 
                        procDict=None, type='ALLD',
                        Gamma=1.4, Cv=1.7857142857142865, MuS=1.e-08, 
                        Cs=0.3831337844872463, Ts=1.0, alpha=1., typeTransfer=0):

    if procDict is None: procDict = Cmpi.getProcDict(aD)
    if graph is None: graph = Cmpi.computeGraph(aD, type=type)

    # Transferts locaux/globaux
    # Calcul des solutions interpolees par arbre donneur
    # On envoie aussi les indices receveurs pour l'instant
    datas = {}
    zonesD = Internal.getZones(aD)
    for zD in zonesD:
        infos = X.setInterpTransfersD(zD, variables=variables, cellNVariable=cellNVariable,
                                      variablesIBC=variablesIBC, 
                                      bcType=bcType, varType=varType, compact=compact,
                                      Gamma=Gamma, Cv=Cv, MuS=MuS, Cs=Cs, Ts=Ts, alpha=alpha)
        for n in infos:
            rcvName = n[0]
            proc = procDict[rcvName]
            if proc == Cmpi.rank:
                field = n[1]
                #print('direct', Cmpi.rank, rcvName)
                if field != []:
                    listIndices = n[2]
                    z = Internal.getNodeFromName2(aR, rcvName)
                    if typeTransfer==0: C._setPartialFields(z, [field], [listIndices], loc=n[3])
                    else:
                        C._updatePartialFields(z, [field], [listIndices], loc=n[3])

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
                if typeTransfer==0: C._setPartialFields(z, [field], [listIndices], loc=n[3])
                else:
                    C._updatePartialFields(z, [field], [listIndices], loc=n[3])

    return None

