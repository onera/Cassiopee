# Class for FastS "IBM"+"Overset" prepare and compute
import Fast.PyTree as Fast
import Converter.PyTree as C
import Generator.PyTree as G
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.PyTree as X
import Connector.ToolboxIBM as TIBM
import Dist2Walls.PyTree as DTW
import Distributor2.PyTree as D2
import Initiator.PyTree as I
import Converter.Mpi as Cmpi
import Converter.Filter as Filter
import Generator
import Post.PyTree as P
import Connector.connector as connector
import KCore.test as test
import Connector.Mpi as Xmpi
import Converter.Distributed as Distributed
import Connector.OversetData as XOD
from Apps.Fast.Common import Common
import Compressor.PyTree as Compressor 

# return True si z est une zone chimera
def isZoneChimera(z):
    if z[3] != 'Zone_t': print('Warning: isZoneChimera: not a Zone node.')
    n = Internal.getNodeFromType1(z, 'ZoneGridConnectivity_t')
    if n is None: return False
    n = Internal.getNodeFromType1(n, 'GridConnectivity_t')
    if n is None: return False
    n = Internal.getNodeFromType1(n, 'GridConnectivityType_t')
    if n is None: return False
    if Internal.getValue(n) == 'Overset': return True
    return False

#================================================================================
# IBMO prepare
#
#================================================================================
def prepare(t_case, t_out, tc_out,   
            vmin=21, check=False, NP=0,
            frontType=1, 
            expand=3, distrib=False, tinit=None, initWithBBox=-1.):
    DEPTH = 2
    IBCType = 1
    if frontType==2: 
        print("Warning: IBMO.prepare: not fully implemented yet for frontType=2 algorithm.")
        frontType=1

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD

    # Extraction de la liste des dfars de tb
    zones = Internal.getZones(tb)
    dfarList = [10.]*len(zones)
    for c, z in enumerate(zones): 
        n = Internal.getNodeFromName2(z, 'dfar')
        if n is not None: dfarList[c] = Internal.getValue(n)*1.

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree defined in %s.'%FILE)
    model = Internal.getValue(model)

    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced
    
    # tbchim : arbre des corps chimeres 
    tbchim = C.newPyTree()
    bases = Internal.getBases(tb)
    for b in bases:
        isChimera = False
        for z in b[2]:
            if isZoneChimera(z): isChimera = True; break
        if isChimera:
            base = Internal.newCGNSBase(b[0], parent=tbchim)
            zones = Internal.getZones(b)
            ext = C.extractBCOfType(zones, 'BCWall')
            base[2] += ext

    # tbibm : arbre des corps ibms en z = 0
    tbibm = C.newPyTree()
    bases = Internal.getBases(tb)
    for b in bases:
        zones = Internal.getZones(b)
        isChimera = False
        for z in zones:
            if isZoneChimera(z): isChimera = True; break
        if not isChimera: tbibm[2].append(b) 
    if dimPb == 2:
        C._initVars(tbchim,'CoordinateZ',0.)
        C._initVars(tbibm,'CoordinateZ',0.)

    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')

    # construction de l'arbre des corps pour l'octree : tbo
    tbo = Internal.copyRef(tb)
    bases = Internal.getBases(tbo)
    for b in bases:
        c = 0
        for z in b[2]:
            if z[3] == 'Zone_t':
                if isZoneChimera(z): 
                    ext = C.extractBCOfType(z, 'BCOverlap')
                    if len(ext)>0: b[2][c] = ext[0]
                    if len(ext)>1: b[2] += ext[1:]
            c += 1

    o = TIBM.buildOctree(tbo, snearFactor=1., dfarList=dfarList, 
                         dimPb=dimPb, vmin=vmin, rank=rank,
                         expand=expand)

    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = TIBM.buildParentOctrees__(o, tbo, snearFactor=4., dfarList=dfarList,
                                        dimPb=dimPb, vmin=vmin, rank=rank)
    test.printMem(">>> Octree unstruct [end]")

    # Split octree
    test.printMem(">>> Octree unstruct split [start]")
    bb = G.bbox(o)
    NPI = Cmpi.size
    if NPI == 1: p = Internal.copyRef(o) # keep reference
    else: p = T.splitNParts(o, N=NPI, recoverBC=False)[rank]
    del o
    test.printMem(">>> Octree unstruct split [end]")

    # fill vmin + merge in parallel
    test.printMem(">>> Octree struct [start]")
    res = TIBM.octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000)
    del p 
    if parento is not None:
        for po in parento: del po
    t = C.newPyTree(['CARTESIAN', res])
    zones = Internal.getZones(t)
    for z in zones: z[0] = z[0]+'X%d'%rank
    Cmpi._setProc(t, rank)
    C._addState(t, 'EquationDimension', dimPb)
    test.printMem(">>> Octree struct [end]")

    # extend cart grids
    test.printMem(">>> extended cart grids [start]")
    tbb = Cmpi.createBBoxTree(t)
    interDict = X.getIntersectingDomains(tbb)
    graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
    del tbb
    Cmpi._addXZones(t, graph, variables=[], cartesian=True)
    test.printMem(">>> extended cart grids [after add XZones]")
    zones = Internal.getZones(t)
    coords = C.getFields(Internal.__GridCoordinates__, zones, api=2)
    coords = Generator.generator.extendCartGrids(coords, 2+1, 1)
    C.setFields(coords, zones, 'nodes')
    Cmpi._rmXZones(t)
    coords = None; zones = None
    test.printMem(">>> extended cart grids (after rmXZones) [end]")

    TIBM._addBCOverlaps(t, bbox=bb)
    TIBM._addExternalBCs(t, bbox=bb, dimPb=dimPb)
    dz = 0.01
    if dimPb == 2:
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)

    # Construction de l'arbre des maillages chimeres
    t2 = C.newPyTree()
    bases = Internal.getBases(tb)
    for b in bases:
        zones = Internal.getZones(b)
        isChimera = False
        for z in zones:
            if isZoneChimera(z): isChimera = True; break
        if isChimera: t2[2].append(b) 
    # Ajout des ghost cells sur maillages chimeres
    C.addState(t2, 'EquationDimension', dimPb)
    Internal._addGhostCells(t2, t2, 2, adaptBCs=1, fillCorner=0)
    if dimPb == 2: 
        T._addkplane(t2)
        T._contract(t2, (0,0,0), (1,0,0), (0,1,0), dz)
        T._makeDirect(t2)

    # tball : arbre de tous les corps
    test.printMem(">>> Wall distance on Chimera grids [start]")
    tball = C.newPyTree()
    tball[2] += Internal.getBases(tbchim) + Internal.getBases(tbibm)
    if dimPb == 2:
        ZMEAN = dz*0.5
        C._initVars(tball, 'CoordinateZ', ZMEAN)
    # Distances a la paroi sur les grilles chimeres
    DTW._distance2Walls(t2, tball, loc='centers', type='ortho')
    test.printMem(">>> Wall distance on Chimera grids [end]")

    # ReferenceState 
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)

    # Distance a la paroi par rapport aux corps IBMs
    test.printMem(">>> Wall distance for Cartesian grids [start]")
    DTW._distance2Walls(t, tball, type='ortho', signed=0, dim=dimPb, loc='centers')
    test.printMem(">>> Wall distance for Cartesian grids [end]")
    if dimPb == 2: C._initVars(tball, 'CoordinateZ', 0.)

    # Calcul du cellNChim sur CART
    C._initVars(t, 'centers:cellNChim', 1.)
    X._applyBCOverlaps(t, depth=2, loc='centers', val=2, cellNName='cellNChim')
        
    # Calcul du cellN sur bases chimeres
    C._initVars(t2, 'centers:cellNChim', 1.)
    X._applyBCOverlaps(t2, depth=4, val=2, cellNName='cellNChim')
    X._applyBCOverlaps(t2, depth=2, val=0, cellNName='cellNChim')

    # Masquage de CART par chimere
    test.printMem(">>> Blanking [start]")
    if dimPb == 2:
        # Creation du corps 2D pour le preprocessing IBC
        T._addkplane(tbchim)
        T._contract(tbchim, (0,0,0), (1,0,0), (0,1,0), dz)
        T._addkplane(tbibm)
        T._contract(tbibm, (0,0,0), (1,0,0), (0,1,0), dz)

    # Masquage IBC
    C._initVars(t, 'centers:cellNIBC', 1.)
    t  = TIBM.blankByIBCBodies(t, tbibm, 'centers', dimPb, cellNName='cellNIBC')
    t  = TIBM.blankByIBCBodies(t, tbchim, 'centers', dimPb, cellNName='cellNChim')
    t2 = TIBM.blankByIBCBodies(t2, tbibm, 'centers', dimPb, cellNName='cellNChim')
    # a faire : masquer les bases chimeres entre elles

    # Signe la distance en fonction de cellNIBC et cellNChim
    TIBM._signDistance(t)

    # Points interpoles autour des points masques
    X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNChim',addGC=False)
    X._setHoleInterpolatedPoints(t2,depth=1,dir=1,loc='centers',cellNName='cellNChim',addGC=False)

    # determination des pts IBC
    X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNIBC',addGC=False) # pour les gradients
    if frontType < 2:
        X._setHoleInterpolatedPoints(t,depth=2,dir=0,loc='centers',cellNName='cellNIBC',addGC=False)
    else:
        X._setHoleInterpolatedPoints(t,depth=3,dir=0,loc='centers',cellNName='cellNIBC',addGC=False)         
    TIBM._removeBlankedGrids(t, loc='centers')
    test.printMem(">>> Blanking [end]")

    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
        dims = Internal.getZoneDim(i)
        npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))

    C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')
    for z in Internal.getZones(t):
        connector._updateNatureForIBM(z, 1,
                                      Internal.__GridCoordinates__,
                                      Internal.__FlowSolutionNodes__,
                                      Internal.__FlowSolutionCenters__)
        
    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # tp arbre complet (t cartesien + t2 chimere)
    tp1 = Internal.copyRef(t) # cartesien
    tp2 = Internal.copyRef(t2) # chimere
    C._initVars(tp2,'{centers:cellN}={centers:cellNChim}')
    # arbre complet tp
    tp = Internal.merge([tp1, tp2])
    test.printMem(">>> Interpdata [start]")
    tpc = C.node2Center(tp)

    # setInterpData parallel pour le chimere
    tbbc = Cmpi.createBBoxTree(tpc)
    interDict = X.getIntersectingDomains(tbbc)
    graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(tpc, graph, variables=['cellN'], cartesian=False)
    test.printMem(">>> Interpdata [after addXZones]")
    
    procDict = Cmpi.getProcDict(tpc)
    datas = {}

    # dictionnaire des bases pour pouvoir ensuite distinguer les donneurs Chimere inter base et abutting intra bases
    nobOfDnrBases = {}
    for nobd in range(len(tbbc[2])):
        if tbbc[2][nobd][3] == 'CGNSBase_t':
            for nozd in range(len(tbbc[2][nobd][2])):
                zdnr = tbbc[2][nobd][2][nozd]
                if zdnr[3] == 'Zone_t':
                    zdnrname = zdnr[0]
                    nobOfDnrBases[zdnrname]=nobd
    for nob in range(len(tp[2])):
        if Internal.getType(tp[2][nob])=='CGNSBase_t':
            basercv = tp[2][nob]
            basername = Internal.getName(basercv)

            # 1st step : Abutting IDs
            isCartBase=False
            if basername == 'CARTESIAN':
                isCartBase=True
            else:
                # on construit les donnees des raccords match
                X._setInterpData(tp[2][nob], tpc[2][nob], nature=1, penalty=1, loc='centers', storage='inverse', 
                                 sameName=1, interpDataType=1, itype='abutting')

            # 2nd step : Chimera
            for zrcv in Internal.getZones(basercv):
                zrname = Internal.getName(zrcv)
                dnrZonesChim = []
                if not isCartBase:
                    print(zrcv[0], ' : ' ,interDict[zrname])
                for zdname in interDict[zrname]:
                    zd = Internal.getNodeFromName2(tpc, zdname)
                    if isCartBase:
                        dnrZonesChim.append(zd)
                    else:                        
                        if nobOfDnrBases[zdname]!=nob:
                            dnrZonesChim.append(zd)

                # A OPTIMISER : des cartesiens et des ADT
                X._setInterpData(zrcv, dnrZonesChim, nature=1, penalty=1, loc='centers', storage='inverse', 
                                 sameName=1, interpDataType=1, itype='chimera')

                for zd in dnrZonesChim:
                    zdname = zd[0]
                    destProc = procDict[zdname]

                    IDs = []
                    for i in zd[2]:
                        if i[0][0:2] == 'ID':
                            if Internal.getValue(i)==zrname: IDs.append(i)

                    if IDs != []:
                        if destProc == rank:                
                            zD = Internal.getNodeFromName2(tpc, zdname)
                            zD[2] += IDs
                        else:
                            if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                            else: datas[destProc].append([zdname,IDs])
                    else:
                        if destProc not in datas: datas[destProc] = []
    nobOfDnrBases=None
    Cmpi._rmXZones(tpc)
    test.printMem(">>> Interpdata [after rmXZones]")
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tpc, zname)
                zD[2] += IDs
    datas = {}; destDatas = None; graph={}
    test.printMem(">>> Interpdata [after free]")
    test.printMem(">>> Interpdata [end]")

    # Traitement IBC : n est fait que sur le cartesien -> 
    tp_cart = C.newPyTree(); tp_cart[2].append(Internal.copyRef(tp[2][1]))
    tpc_cart = C.newPyTree(); tpc_cart[2].append(Internal.copyRef(tpc[2][1]))
    C._rmVars(tp,['centers:cellNFront','centers:cellNIBC'])

    # suppression des noeuds ID_* non cartesiens pour ne pas planter les transferts  
    for ID in Internal.getNodesFromType3(tpc_cart,'ZoneSubRegion_t'):
        if ID[0][0:7]!='ID_Cart': Internal._rmNode(tpc_cart,ID)

    C._initVars(tp_cart,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(tp_cart,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(tp_cart,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')    
    C._cpVars(tp_cart,'centers:cellNIBC',tp_cart,'centers:cellN')
    C._cpVars(tp_cart,'centers:cellN',tpc_cart,'cellN')
    
    # Transfert du cellNFront
    C._cpVars(tp_cart,'centers:cellNFront',tpc_cart,'cellNFront')
    # propager cellNFront
    Xmpi._setInterpTransfers(tp_cart,tpc_cart,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)
    # fin transfert front

    C._rmVars(tp_cart,['centers:cellNFront'])
    C._cpVars(tp_cart,'centers:TurbulentDistance',tpc_cart,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(tp[2][1],'centers:TurbulentDistance'))
    P._computeGrad2(tp_cart, 'centers:TurbulentDistance')
    
    test.printMem(">>> Building IBM front [start]")
    front = TIBM.getIBMFront(tpc_cart, 'cellNFront', dim=dimPb, frontType=frontType)
    front = TIBM.gatherFront(front)
    
    if check and rank == 0: C.convertPyTree2File(front, 'front.cgns')

    zonesRIBC = []
    for zrcv in Internal.getZones(tp_cart):
        if C.getMaxValue(zrcv, 'centers:cellNIBC')==2.:
            zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

    nbZonesIBC = len(zonesRIBC)
    if nbZonesIBC == 0: 
        res = [{},{},{}]
    else:
        res = TIBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tbibm, tfront=front, frontType=frontType,
                                   cellNName='cellNIBC', depth=2, IBCType=IBCType)

    front = None
    # cleaning
    C._rmVars(tpc_cart,['cellNChim','cellNIBC','TurbulentDistance','cellNFront'])
    # dans t, il faut cellNChim et cellNIBCDnr pour recalculer le cellN a la fin
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBC']
    C._rmVars(tp_cart, varsRM)
    test.printMem(">>> Building IBM front [end]")

    # graphe d'intersection des pts images de ce proc et des zones de tbbc
    tbbc_cart = C.newPyTree(); tbbc_cart[2].append(tbbc[2][1])
    zones = Internal.getZones(tbbc_cart)    
    allBBs = []
    dictOfCorrectedPtsByIBCType = res[0]
    dictOfWallPtsByIBCType = res[1] 
    dictOfInterpPtsByIBCType = res[2]
    interDictIBM={}
    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    zrname = zonesRIBC[nozr][0]
                    interpPtsBB = Generator.BB(allInterpPts[nozr])
                    for z in zones:
                        bba = C.getFields('GridCoordinates', z)[0]
                        if Generator.bboxIntersection(interpPtsBB,bba,isBB=True):
                            zname = z[0]
                            popp = Cmpi.getProc(z)
                            Distributed.updateGraph__(graph, popp, rank, zname)
                            if zrname not in interDictIBM: interDictIBM[zrname]=[zname]
                            else: 
                                if zname not in interDictIBM[zrname]: interDictIBM[zrname].append(zname)
    else: graph={}
    del tbbc_cart
    del tbbc
    allGraph = Cmpi.KCOMM.allgather(graph)
    #if rank == 0: print allGraph

    graph = {}
    for i in allGraph:
        for k in i:
            if not k in graph: graph[k] = {}
            for j in i[k]:
                if not j in graph[k]: graph[k][j] = []
                graph[k][j] += i[k][j]
                graph[k][j] = list(set(graph[k][j])) # pas utile?

    test.printMem(">>> Interpolating IBM [start]")
    # keyword subr=False to avoid memory overflow
    Cmpi._addXZones(tpc_cart, graph, variables=['cellN'], cartesian=True, subr=False)
    test.printMem(">>> Interpolating IBM [after addXZones]")

    ReferenceState = Internal.getNodeFromType2(tp, 'ReferenceState_t')
    nbZonesIBC = len(zonesRIBC)

    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in range(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    zrcv = zonesRIBC[nozr]
                    zrname = zrcv[0]
                    dnrZones = []
                    for zdname in interDictIBM[zrname]:
                        zd = Internal.getNodeFromName2(tpc_cart, zdname)
                        #if zd is not None: dnrZones.append(zd)
                        if zd is None: print('!!!Zone None', zrname, zdname)
                        else: dnrZones.append(zd)
                    XOD._setIBCDataForZone__(zrcv, dnrZones, allCorrectedPts[nozr], allWallPts[nozr], allInterpPts[nozr],
                                             nature=1, penalty=1, loc='centers', storage='inverse', dim=dimPb,
                                             interpDataType=0, ReferenceState=ReferenceState, bcType=ibcTypeL)

                    nozr += 1
                    for zd in dnrZones:       
                        zdname = zd[0]
                        destProc = procDict[zdname]
                        
                        #allIDs = Internal.getNodesFromName(zd, 'IBCD*')
                        #IDs = [] 
                        #for zsr in allIDs:
                        #    if Internal.getValue(zsr)==zrname: IDs.append(zsr)

                        IDs = []
                        for i in zd[2]:
                            if i[0][0:4] == 'IBCD':
                                if Internal.getValue(i)==zrname: IDs.append(i)

                        if IDs != []:
                            if destProc == rank:                
                                zD = Internal.getNodeFromName2(tpc_cart,zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []

    test.printMem(">>> Interpolating IBM [end]")
    Cmpi._rmXZones(tpc_cart)
    dictOfCorrectedPtsByIBCType = None
    dictOfWallPtsByIBCType = None
    dictOfInterpPtsByIBCType = None
    interDictIBM = None
    test.printMem(">>> Interpolating IBM [after rm XZones]")

    for zd_cart in Internal.getZones(tpc_cart):
        zd = Internal.getNodeFromName1(tpc[2][1],zd_cart[0])
        IBCDS = Internal.getNodesFromName(zd_cart,'IBCD_*')
        zd[2]+=IBCDS
    del tpc_cart

    Internal._rmNodesByName(tpc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(tpc, Internal.__GridCoordinates__)
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IBCDs = n[1]
            if IBCDs != []:
                zD = Internal.getNodeFromName2(tpc, zname)
                zD[2] += IBCDs

    datas = {}; graph = {}
    for i in range(Cmpi.size): datas[i] = [] # force

    # uniquement le cartesien
    C._initVars(tp_cart[2][1],'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:cellNChim','centers:cellNIBCDnr']
    if model == 'Euler': varsRM += ['centers:TurbulentDistance']
    C._rmVars(tp, varsRM)
    C._cpVars(tp_cart[2][1],'centers:cellN',tp[2][1],'centers:cellN')
    del tp_cart

    #-------------------------------------------
    # Recomputes distance field for Musker only
    #-------------------------------------------
    if model != 'Euler':
        ibctypes = set()
        for node in Internal.getNodesFromName(tball,'ibctype'):
            ibctypes.add(Internal.getValue(node))
        ibctypes = list(ibctypes)
        if model != 'Euler' and ('outpress' in ibctypes or 'inj' in ibctypes or 'slip' in ibctypes):
            test.printMem(">>> wall distance for viscous wall only [start]")
            for z in Internal.getZones(tball):
                ibc = Internal.getNodeFromName(z,'ibctype')
                if ibc is not None:
                    ibctypez = Internal.getValue(ibc)
                    if ibctypez=='outpress' or ibctypez=='inj' or ibctypez=='slip':
                        Internal._rmNode(tball,z)

            if dimPb == 2:
                tball2 = C.initVars(tball, 'CoordinateZ', dz*0.5)
                DTW._distance2Walls(tp,tball2,type='ortho', signed=0, dim=dimPb, loc='centers')
            else:
                DTW._distance2Walls(tp,tball,type='ortho', signed=0, dim=dimPb, loc='centers')
            test.printMem(">>> wall distance for Musker only [end]")

    test.printMem(">>> Saving [start]")
    # Sauvegarde des infos IBM
    if check:
        test.printMem(">>> Saving IBM infos [start]")
        tibm = TIBM.extractIBMInfo(tpc)

        # Avoid that two procs write the same information
        for z in Internal.getZones(tibm):
           if int(z[0][-1]) != rank:
              Internal._rmNodesByName(tibm, z[0])

        Cmpi.convertPyTree2File(tibm, 'IBMInfo.cgns')
        test.printMem(">>> Saving IBM infos [end]")
        del tibm

    # distribution par defaut (sur NP)
    # Perform the final distribution
    if distrib:
        tbbc = Cmpi.createBBoxTree(tpc)
        if NP == 0: NP = Cmpi.size
        stats = D2._distribute(tbbc, NP, algorithm='graph', useCom='ID')
        D2._copyDistribution(tpc, tbbc)
        D2._copyDistribution(tp, tbbc)
        del tbbc

    for zc in Internal.getZones(tpc[2][2:]):
        SDD = Internal.getNodesFromName1(zc,".Solver#define")
        Internal._rmNodesFromName(SDD, 'dfar')
        Internal._rmNodesFromName(SDD, 'snear')
    for zc in Internal.getZones(tp[2][2:]):
        SDD = Internal.getNodesFromName1(zc,".Solver#define")
        Internal._rmNodesFromName(SDD, 'dfar')
        Internal._rmNodesFromName(SDD, 'snear')

    # Save tc
    for zone in Internal.getZones(tpc[2][1]):
        zname=Internal.getName(zone)
        zname = zname.split('/')
        if len(zname)==2: Internal.setName(zone,zname[1])
    for zone in Internal.getZones(tpc[2][2]):
        Internal._rmNodesFromName(zone,'RANSLES')
    if isinstance(tc_out, str): Cmpi.convertPyTree2File(tpc, tc_out, ignoreProcNodes=True)
    
    # Initialisation
    if tinit is None: I._initConst(tp, loc='centers')
    else: 
        import Post.PyTree as Pmpi
        tp = Pmpi.extractMesh(tinit, tp, mode='accurate')
    if model != "Euler": C._initVars(tp, 'centers:ViscosityEddy', 0.)

    # Init with BBox : uniquement sur le cartesien autour des IBM - curviligne tel quel
    if initWithBBox>0.:
        print('initialisation par bounding box')
        import Geom.PyTree as D
        bodybb = C.newPyTree(['Base'])
        for base in Internal.getBases(tbibm):
            bbox = G.bbox(base)
            bodybbz = D.box(tuple(bbox[:3]),tuple(bbox[3:]), N=2, ntype='STRUCT')
            Internal._append(bodybb,bodybbz,'Base')
        T._scale(bodybb, factor=(initWithBBox,initWithBBox,initWithBBox))
        tbb = G.BB(tp)
        interDict = X.getIntersectingDomains(tbb,bodybb,taabb=tbb,taabb2=bodybb)
        for zone in Internal.getZones(tp[2][1]): 
            zname = Internal.getName(zone)
            if interDict[zname] != []:
                C._initVars(zone, 'centers:MomentumX', 0.)
                C._initVars(zone, 'centers:MomentumY', 0.)
                C._initVars(zone, 'centers:MomentumZ', 0.)

    # Save t
    if isinstance(t_out, str):
        #Compressor._compressCartesian(tp[2][1])
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    if Cmpi.size > 1: Cmpi.barrier()    
    return tp, tpc

#====================================================================================
class IBMO(Common):
    """Preparation et caculs avec le module FastS."""
    def __init__(self, format=None, numb=None, numz=None):
        Common.__init__(self, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["stephanie@onera.fr", "ash@onera.fr"]
        self.cartesian = True
        
    # Prepare 
    def prepare(self, t_case, t_out, tc_out, 
                vmin=21, check=False, frontType=1, NP=None, expand=3):
        if NP is None: NP = Cmpi.size
        if NP == 0: print('Preparing for a sequential computation.')
        else: print('Preparing for a computation on %d processors.'%NP)
        ret = prepare(t_case, t_out, tc_out, 
                      vmin=vmin, check=check, NP=NP,  
                      frontType=frontType, expand=expand)
        return ret