# Class for FastS "IBM"+"Overset" prepare and compute
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

def prepareMotion(t_case, t_out, tc_out, vmin=21, check=False, NP=0,
                  frontType=1, tbox=None, snearsf=None,
                  expand=3, distrib=False, tinit=None, initWithBBox=-1., dfarDir=0):
    import RigidMotion.PyTree as R
    
    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD
    tolMatch = 1.e-6
    DEPTH = 2
    IBCType = 1
    if frontType!=1: 
        print("Warning: IBMO.prepareMotion: currently implemented for frontType=1 algorithm only.")
        frontType=1

    if isinstance(t_case, str): 
        h = Filter.Handle(t_case)
        if distrib:
            tb = h.loadAndSplit(NParts=Cmpi.size)
        else:
            tb = h.loadAndSplit(NParts=max(NP, Cmpi.size))

        #tb = Cmpi.convertFile2SkeletonTree(t_case)
    else: tb = t_case

    # refinementSurfFile: surface meshes describing refinement zones
    if tbox is not None:
        if isinstance(tbox, str): tbox = C.convertFile2PyTree(tbox)
        else: tbox = tbox
        if snearsf is None:
            snearsf = []
            zones = Internal.getZones(tbox)
            for z in zones:
                sn = Internal.getNodeFromName2(z, 'snear')
                if sn is not None: snearsf.append(Internal.getValue(sn))
                else: snearsf.append(1.)

    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    # tbchim : arbre des corps chimeres 
    tbchim = C.newPyTree()
    tbov = C.newPyTree()
    tbblank= C.newPyTree()
    # tbibm : arbre des corps ibms en z = 0
    tbibm = C.newPyTree()    
    # chargement des zones Chimere par proc 
    bases = Internal.getBases(tb)
    baseNamesChim=[]

    for b in bases:
        isChimera = False
        for z in Internal.getZones(b):
            if isZoneChimera(z): 
                isChimera = True
                baseNamesChim.append(b[0])
                break

        if isChimera:
            zoneChimNames=[]
            for z in Internal.getZones(b): 
                proc = Internal.getNodeFromName2(z,'proc')
                proc = Internal.getValue(proc)
                if proc == rank:
                    zoneChimNames.append(z[0])
                else: 
                    Internal._rmNodesFromName1(b,z[0])
            if zoneChimNames != []:
                Cmpi._readZones(tb, t_case, rank=rank, zoneNames=zoneChimNames)
                basechim = Internal.newCGNSBase(b[0], parent=tbchim)
                baseov = Internal.newCGNSBase(b[0], parent=tbov)
                zones = Internal.getZones(b)
                for z in zones:
                    TM = Internal.getNodeFromName(z,'TimeMotion')
                    SD = Internal.getNodeFromName(z,'.Solver#define')
                    SP = Internal.getNodeFromName(z,'.Solver#Param')

                    zwalls = C.extractBCOfType(z, 'BCWall')
                    now=1
                    for zw in zwalls: 
                        zw[0] = b[0]+'_'+z[0]+'_W'+str(now)
                        if TM is not None: zw[2].append(TM)
                        if SD is not None: zw[2].append(SD)
                        if SP is not None: zw[2].append(SP)
                        now+=1
                    basechim[2] += zwalls

                    zovs = C.extractBCOfType(z, 'BCOverlap')
                    nov=1
                    for zo in zovs: 
                        zo[0] = b[0]+'_'+z[0]+'_OV'+str(nov)
                        if TM is not None: zo[2].append(TM)
                        if SD is not None: zo[2].append(SD)
                        if SP is not None: zo[2].append(SP)
                        nov+=1
                    baseov[2] += zovs

        else:
            zoneIBMNames = []
            for z in Internal.getZones(b): 
                proc = Internal.getNodeFromName2(z,'proc')
                proc = Internal.getValue(proc)
                if proc == rank:
                    zoneIBMNames.append(z[0])
                else: 
                    Internal._rmNodesFromName1(b,z[0])

            if zoneIBMNames!=[]:
                Cmpi._readZones(tb, t_case, rank=rank, zoneNames=zoneIBMNames)
                tbibm[2].append(b) 

    if dimPb == 2:
        C._initVars(tb, 'CoordinateZ', 0.)
        C._initVars(tbchim, 'CoordinateZ', 0.)
        C._initVars(tbibm, 'CoordinateZ', 0.)
        C._initVars(tbov, 'CoordinateZ', 0.)


    # allgather: everybody sees all the bodies and overlaps 
    tbchim = Cmpi.allgatherTree(tbchim)
    tbibm  = Cmpi.allgatherTree(tbibm)
    tbov   = Cmpi.allgatherTree(tbov)
    Internal._rmNodesFromName(tbov,'GridCoordinates#Init')
    Internal._rmNodesFromName(tbchim,'GridCoordinates#Init')
    Internal._rmNodesFromName(tb,'GridCoordinates#Init')

    R._copyGrid2GridInit(tbov)
    R._copyGrid2GridInit(tbchim)

    # Reference state
    refstate = C.getState(tb)

    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree.')
    model = Internal.getValue(model)

    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')

    # construction de l'arbre des corps pour l'octree : tbo
    tbo = C.newPyTree(['Surf_octree'])
    R._evalPosition(tbov, 0.)
    tbo[2][1][2] = Internal.getZones(tbibm)+Internal.getZones(tbov)
    # Extraction de la liste des dfars de tb
    zones = Internal.getZones(tbo)
    dfarList = [10.]*len(zones)
    for c, z in enumerate(zones): 
        n = Internal.getNodeFromName2(z, 'dfar')
        if n is not None: dfarList[c] = Internal.getValue(n)*1.

    o = TIBM.buildOctree(tbo, snearFactor=1., dfarList=dfarList, 
                         tbox=tbox, snearsf=snearsf,
                         dimPb=dimPb, vmin=vmin, rank=rank,
                         expand=expand, dfarDir=dfarDir)

    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = TIBM.buildParentOctrees__(o, tbo, snearFactor=4., dfarList=dfarList,
                                        tbox=tbox, snearsf=snearsf,
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
    coords, rinds = Generator.extendCartGrids(coords, ext=DEPTH+1, optimized=1, extBnd=0)
    C.setFields(coords, zones, 'nodes')
    for noz in range(len(zones)):
        Internal.newRind(value=rinds[noz], parent=zones[noz])
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
        if b[0] in baseNamesChim:
            t2[2].append(b)

    # Ajout des ghost cells sur maillages chimeres
    C.addState(t2, 'EquationDimension', dimPb)
    if Cmpi.size > 1:
        C._rmBCOfType(t2,'BCMatch')
        Cmpi._addBXZones(t2, depth=DEPTH*2, allB=False)
        t2 = X.connectMatch(t2, tol=tolMatch, dim=dimPb)
    Internal._addGhostCells(t2, t2, DEPTH, adaptBCs=1, fillCorner=0)

    # Suppression des XZones et correction des matchs 
    if Cmpi.size > 1:Cmpi._rmBXZones(t2)

    # Fusion des fenetres des raccords 
    if Cmpi.size > 1: t2 = Xmpi.mergeWindows(t2)

    if dimPb == 2: 
        T._addkplane(t2)
        T._contract(t2, (0,0,0), (1,0,0), (0,1,0), dz)
        T._makeDirect(t2)

    # tball : arbre de tous les corps
    test.printMem(">>> Wall distance on Chimera grids [start]")
    tball = C.newPyTree(['AllBodies'])
    tball[2][1][2] = Internal.getZones(tbchim) + Internal.getZones(tbibm)

    if dimPb == 2:
        ZMEAN = dz*0.5
        C._initVars(tball, 'CoordinateZ', ZMEAN)
    # Distances a la paroi sur les grilles curvilignes chimere
    DTW._distance2Walls(t2, tbchim, loc='centers', type='ortho')
    test.printMem(">>> Wall distance on Chimera grids [end]")

    # ReferenceState 
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)

    # Distance a la paroi pour les grilles cartesiennes
    test.printMem(">>> Wall distance for Cartesian grids [start]")
    DTW._distance2Walls(t, tbibm, type='ortho', signed=0, dim=dimPb, loc='centers')
    test.printMem(">>> Wall distance for Cartesian grids [end]")
    if dimPb == 2: C._initVars(tball, 'CoordinateZ', 0.)

    # Calcul du cellNChim sur CART
    C._initVars(t, 'centers:cellNChim', 1.)
    X._applyBCOverlaps(t, depth=DEPTH, loc='centers', val=2, cellNName='cellNChim')

    # Calcul du cellN sur bases Chimere - que du cellN Motion
    C._initVars(t2, 'centers:cellN', 1.)
    #===============================
    # donnees de raccords match
    #===============================
    test.printMem(">>> Abutting data for curvilinear grids [start]")
    tc2 = C.node2Center(t2)
    for nob in range(len(t2[2])):
        if Internal.getType(t2[2][nob])=='CGNSBase_t':
            graphMatch = Cmpi.computeGraph(t2[2][nob], type='match', reduction=True)
            Cmpi._addXZones(tc2[2][nob], graphMatch, noCoordinates=True, cartesian=False)
            Cmpi._addXZones(t2[2][nob], graphMatch, noCoordinates=True, cartesian=False)
            # on construit les donnees des raccords match
            X._setInterpData(t2[2][nob], tc2[2][nob], nature=1, penalty=1, loc='centers', 
                             storage='inverse', sameName=1, dim=dimPb, itype='abutting')
            Cmpi._rmXZones(tc2[2][nob])
            Cmpi._rmXZones(t2[2][nob])
               
    graphMatch={}
    test.printMem(">>> Abutting data [after free]")
    test.printMem(">>> Abutting data [end]")

    # C.convertPyTree2File(t2,'t_chim_%d.cgns'%rank)
    # C.convertPyTree2File(tc2,'tc_chim_%d.cgns'%rank)

    # Masquage de CART 
    test.printMem(">>> Blanking [start]")
    if dimPb == 2:
        # Creation du corps 2D pour le preprocessing IBC
        T._addkplane(tbibm)
        T._contract(tbibm, (0,0,0), (1,0,0), (0,1,0), dz)

    # Blanking CART base by IBC bodies
    C._initVars(t, 'centers:cellNIBC', 1.)
    t  = TIBM.blankByIBCBodies(t, tbibm, 'centers', dimPb, cellNName='cellNIBC')
    test.printMem(">>> Blanking [end]")

    # Signe la distance en fonction de cellNIBC et cellNChim
    TIBM._signDistance(t)

    # Points interpoles autour des points masques
    X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=1,loc='centers',cellNName='cellNChim',addGC=False)
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

    # setInterpData - Chimere entre grilles cartesiennes
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    test.printMem(">>> Interpdata [start]")
    tc = C.node2Center(t)
    # setInterpData parallel pour le chimere
    tbbc = Cmpi.createBBoxTree(tc)
    interDict = X.getIntersectingDomains(tbbc)
    graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=True)
    test.printMem(">>> Interpdata [after addXZones]")

    procDict = Cmpi.getProcDict(tc)
    datas = {}
    for zrcv in Internal.getZones(t):
        zrname = zrcv[0]
        dnrZones = []
        for zdname in interDict[zrname]:
            zd = Internal.getNodeFromName2(tc, zdname)
            dnrZones.append(zd)
        X._setInterpData(zrcv, dnrZones, nature=1, penalty=1, loc='centers', storage='inverse',
                         sameName=1, interpDataType=0, itype='chimera')
        for zd in dnrZones:
            zdname = zd[0]
            destProc = procDict[zdname]

            #allIDs = Internal.getNodesFromName(zd, 'ID*')
            #IDs = []
            #for zsr in allIDs:
            #    if Internal.getValue(zsr)==zrname: IDs.append(zsr)
            IDs = []
            for i in zd[2]:
                if i[0][0:2] == 'ID':
                    if Internal.getValue(i)==zrname: IDs.append(i)

            if IDs != []:
                if destProc == rank:
                    zD = Internal.getNodeFromName2(tc, zdname)
                    zD[2] += IDs
                else:
                    if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                    else: datas[destProc].append([zdname,IDs])
            else:
                if destProc not in datas: datas[destProc] = []
    Cmpi._rmXZones(tc)
    test.printMem(">>> Interpdata [after rmXZones]")
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IDs = n[1]
            if IDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IDs
    datas = {}; destDatas = None; graph={}
    test.printMem(">>> Interpdata [after free]")
    test.printMem(">>> Interpdata [end]")

    # fin interpData
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')

    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    # propager cellNVariable='cellNFront'
    Xmpi._setInterpTransfers(t,tc,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)

    C._rmVars(t,['centers:cellNFront'])
    C._cpVars(t,'centers:TurbulentDistance',tc,'TurbulentDistance')

    print('Minimum distance: %f.'%C.getMinValue(t,'centers:TurbulentDistance'))
    P._computeGrad2(t, 'centers:TurbulentDistance')

    test.printMem(">>> Building IBM front [start]")
    front = TIBM.getIBMFront(tc, 'cellNFront', dim=dimPb, frontType=frontType)
    front = TIBM.gatherFront(front)
    
    if check and rank == 0: C.convertPyTree2File(front, 'front.cgns')
    zonesRIBC = []
    for zrcv in Internal.getZones(t):
        if C.getMaxValue(zrcv, 'centers:cellNIBC')==2.:
            zrcvname = zrcv[0]; zonesRIBC.append(zrcv)

    nbZonesIBC = len(zonesRIBC)
    if nbZonesIBC == 0:
        res = [{},{},{}]
    else:
        res = TIBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tbibm, tfront=front, frontType=frontType,
                                   cellNName='cellNIBC', depth=DEPTH, IBCType=IBCType)
    # cleaning
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance','cellNFront'])
    # dans t, il faut cellNChim et cellNIBCDnr pour recalculer le cellN a la fin
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBC']
    C._rmVars(t, varsRM)
    front = None
    test.printMem(">>> Building IBM front [end]")

    # Interpolation IBC (front, tbbc)
    # graph d'intersection des pts images de ce proc et des zones de tbbc
    zones = Internal.getZones(tbbc)
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
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=True, subr=False)
    test.printMem(">>> Interpolating IBM [after addXZones]")

    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    nbZonesIBC = len(zonesRIBC)

    for i in range(Cmpi.size): datas[i] = [] # force

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
                        zd = Internal.getNodeFromName2(tc, zdname)
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
                                zD = Internal.getNodeFromName2(tc,zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc].append([zdname,IDs])
                        else:
                            if destProc not in datas: datas[destProc] = []

    test.printMem(">>> Interpolating IBM [end]")
    Cmpi._rmXZones(tc)
    dictOfCorrectedPtsByIBCType = None
    dictOfWallPtsByIBCType = None
    dictOfInterpPtsByIBCType = None
    interDictIBM = None
    test.printMem(">>> Interpolating IBM [after rm XZones]")

    Internal._rmNodesFromType(tc, 'FlowSolution_t')
    #Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IBCDs = n[1]
            if IBCDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IBCDs

    datas = {}; graph = {}
    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:cellNChim','centers:cellNIBCDnr']
    if model == 'Euler': varsRM += ['centers:TurbulentDistance']
    C._rmVars(t, varsRM)
    Internal._rmNodesFromName(tc,'FlowSolution')
    
    # Sauvegarde des infos IBM
    if check:
        test.printMem(">>> Saving IBM infos [start]")
        tibm = TIBM.extractIBMInfo(tc)

        # Avoid that two procs write the same information
        for z in Internal.getZones(tibm):
           if int(z[0][-1]) != rank:
              # Internal._rmNodesByName(tibm, z[0])
              z[0] = z[0]+"%{}".format(rank)

        Cmpi.convertPyTree2File(tibm, 'IBMInfo.cgns')
        test.printMem(">>> Saving IBM infos [end]")
        del tibm

    # merge trees
    tp1 = Internal.copyRef(t)
    tp2 = Internal.copyRef(t2)
    tp = Internal.merge([tp1,tp2])
    tpc1 = Internal.copyRef(tc)
    tpc2 = Internal.copyRef(tc2)
    tpc = Internal.merge([tpc1,tpc2])
    test.printMem(">>> Saving [start]")

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
    R._copyGrid2GridInit(tpc)
    Internal._rmNodesFromName(tpc,'ZoneRind')
    if isinstance(tc_out, str): 
        import Compressor.PyTree as Compressor
        Compressor._compressCartesian(tpc[2][1])
        Cmpi.convertPyTree2File(tpc, tc_out, ignoreProcNodes=True)
    
    # Initialisation
    if tinit is None: I._initConst(tp, loc='centers')
    else: 
        import Post.PyTree as Pmpi
        tp = Pmpi.extractMesh(tinit, tp, mode='accurate')
    if model != "Euler": C._initVars(tp, 'centers:ViscosityEddy', 0.)

    # Init with BBox : uniquement sur le cartesien autour des IBM - curviligne tel quel
    if initWithBBox>0.:
        print('Init momentum to 0 inside the bounding box')
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

    # clean t
    for z in Internal.getZones(tp):
        SDD = Internal.getNodeFromName1(z,".Solver#define")
        if SDD is not None:
            Internal._rmNodesFromName1(SDD,'inv')
            Internal._rmNodesFromName1(SDD,'ibctype')
            Internal._rmNodesFromName1(SDD,'dfar')
            Internal._rmNodesFromName1(SDD,'snear')
    Internal._rmNodesFromName(tp,'GridCoordinates#Init')
    R._copyGrid2GridInit(tp)

    # Save t
    if isinstance(t_out, str):
        #Compressor._compressCartesian(tp[2][1])
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    if Cmpi.size > 1: Cmpi.barrier()    
    return tp, tpc

#================================================================================
# IBMO prepare
#
#================================================================================
def prepare(t_case, t_out, tc_out,   
            vmin=21, check=False, NP=0,
            frontType=1, 
            expand=3, distrib=False, tinit=None, initWithBBox=-1., dfarDir=0):
    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD
    tolMatch = 1.e-6
    DEPTH = 2
    IBCType = 1
    if frontType!=1: 
        print("Warning: IBMO.prepare: currently implemented for frontType=1 algorithm only.")
        frontType=1

    if isinstance(t_case, str): 
        h = Filter.Handle(t_case)
        if distrib:
            tb = h.loadAndSplit(NParts=Cmpi.size)
        else:
            tb = h.loadAndSplit(NParts=max(NP, Cmpi.size))

        #tb = Cmpi.convertFile2SkeletonTree(t_case)
    else: tb = t_case

    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    # tbchim : arbre des corps chimeres 
    tbchim = C.newPyTree()
    tbov = C.newPyTree()
    # tbibm : arbre des corps ibms en z = 0
    tbibm = C.newPyTree()    
    # chargement des zones Chimere par proc 
    bases = Internal.getBases(tb)
    baseNamesChim=[]

    for b in bases:
        isChimera = False
        for z in Internal.getZones(b):
            if isZoneChimera(z): 
                isChimera = True
                baseNamesChim.append(b[0])
                break

        if isChimera:
            zoneChimNames=[]
            for z in Internal.getZones(b): 
                proc = Internal.getNodeFromName2(z,'proc')
                proc = Internal.getValue(proc)
                if proc == rank:
                    zoneChimNames.append(z[0])
                else: 
                    Internal._rmNodesFromName1(b,z[0])
            if zoneChimNames != []:
                Cmpi._readZones(tb, t_case, rank=rank, zoneNames=zoneChimNames)
                basechim = Internal.newCGNSBase(b[0], parent=tbchim)
                zones = Internal.getZones(b)
                basechim[2] += C.extractBCOfType(zones, 'BCWall')
                baseov = Internal.newCGNSBase(b[0], parent=tbov)
                baseov[2] += C.extractBCOfType(zones, 'BCOverlap')

        else:
            zoneIBMNames = []
            for z in Internal.getZones(b): 
                proc = Internal.getNodeFromName2(z,'proc')
                proc = Internal.getValue(proc)
                if proc == rank:
                    zoneIBMNames.append(z[0])
                else: 
                    Internal._rmNodesFromName1(b,z[0])

            if zoneIBMNames!=[]:
                Cmpi._readZones(tb, t_case, rank=rank, zoneNames=zoneIBMNames)
                tbibm[2].append(b) 


    if dimPb == 2:
        C._initVars(tb, 'CoordinateZ', 0.)
        C._initVars(tbchim, 'CoordinateZ', 0.)
        C._initVars(tbibm, 'CoordinateZ', 0.)
        C._initVars(tbov, 'CoordinateZ', 0.)

    # allgather: everybody sees all the bodies and overlaps 
    tbchim = Cmpi.allgatherTree(tbchim)
    tbibm  = Cmpi.allgatherTree(tbibm)
    tbov   = Cmpi.allgatherTree(tbov)

    # Reference state
    refstate = C.getState(tb)

    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree.')
    model = Internal.getValue(model)

    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')

    # construction de l'arbre des corps pour l'octree : tbo
    tbo = C.newPyTree(['Surf_octree'])
    tbo[2][1][2] = Internal.getZones(tbibm)+Internal.getZones(tbov)

    # Extraction de la liste des dfars de tb
    zones = Internal.getZones(tbo)
    dfarList = [10.]*len(zones)
    for c, z in enumerate(zones): 
        n = Internal.getNodeFromName2(z, 'dfar')
        if n is not None: dfarList[c] = Internal.getValue(n)*1.

    o = TIBM.buildOctree(tbo, snearFactor=1., dfarList=dfarList, 
                         dimPb=dimPb, vmin=vmin, rank=rank,
                         expand=expand, dfarDir=dfarDir)

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
    coords, rinds = Generator.extendCartGrids(coords, ext=DEPTH+1, optimized=1, extBnd=0)
    C.setFields(coords, zones, 'nodes')
    for noz in range(len(zones)):
        Internal.newRind(value=rinds[noz], parent=zones[noz])
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
        if b[0] in baseNamesChim:
            t2[2].append(b)

    # Ajout des ghost cells sur maillages chimeres
    C.addState(t2, 'EquationDimension', dimPb)
    C._rmBCOfType(t2,'BCMatch')
    Cmpi._addBXZones(t2, depth=DEPTH*2, allB=False)
    t2 = X.connectMatch(t2, tol=tolMatch, dim=dimPb)
    Internal._addGhostCells(t2, t2, DEPTH, adaptBCs=1, fillCorner=0)

    # Suppression des XZones et correction des matchs 
    Cmpi._rmBXZones(t2)

    # Fusion des fenetres des raccords 
    t2 = Xmpi.mergeWindows(t2)

    if dimPb == 2: 
        T._addkplane(t2)
        T._contract(t2, (0,0,0), (1,0,0), (0,1,0), dz)
        T._makeDirect(t2)

    # tball : arbre de tous les corps
    test.printMem(">>> Wall distance on Chimera grids [start]")
    tball = C.newPyTree(['AllBodies'])
    tball[2][1][2] = Internal.getZones(tbchim) + Internal.getZones(tbibm)
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

    # Distance a la paroi psur les grilles cartesiennes
    test.printMem(">>> Wall distance for Cartesian grids [start]")
    DTW._distance2Walls(t, tball, type='ortho', signed=0, dim=dimPb, loc='centers')
    test.printMem(">>> Wall distance for Cartesian grids [end]")
    if dimPb == 2: C._initVars(tball, 'CoordinateZ', 0.)

    # Calcul du cellNChim sur CART
    C._initVars(t, 'centers:cellNChim', 1.)
    X._applyBCOverlaps(t, depth=DEPTH, loc='centers', val=2, cellNName='cellNChim')
        
    # Calcul du cellN sur bases Chimere
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

    # Blanking IBC
    C._initVars(t, 'centers:cellNIBC', 1.)
    t  = TIBM.blankByIBCBodies(t, tbibm, 'centers', dimPb, cellNName='cellNIBC')
    t  = TIBM.blankByIBCBodies(t, tbchim, 'centers', dimPb, cellNName='cellNChim')
    t2 = TIBM.blankByIBCBodies(t2, tbibm, 'centers', dimPb, cellNName='cellNChim')
    test.printMem(">>> Blanking [end]")

    # 2-Blank between overset curvilinear grids
    test.printMem(">>> Blanking between curvilinear grids [start]")
    for bodyBase in Internal.getBases(tbchim):
        bodyBaseName = Internal.getName(bodyBase)
        for nob in range(len(t2[2])):
            if t2[2][nob][3]=='CGNSBase_t':
                baseName = Internal.getName(t2[2][nob])
                if baseName != bodyBaseName:
                    t2[2][nob] = TIBM.blankByIBCBodies(t2[2][nob], bodyBase, 'centers', dimPb, cellNName='cellNChim')
    test.printMem(">>> Blanking between curvilinear grids [end]")

    # Signe la distance en fonction de cellNIBC et cellNChim
    TIBM._signDistance(t)

    # Points interpoles autour des points masques
    X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=1,loc='centers',cellNName='cellNChim',addGC=False)
    X._setHoleInterpolatedPoints(t2,depth=DEPTH,dir=1,loc='centers',cellNName='cellNChim',addGC=False)

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

            isCartBaseR=False
            if basername == 'CARTESIAN':
               isCartBaseR=True

            # 2nd step : Chimera
            for zrcv in Internal.getZones(basercv):
                zrname = Internal.getName(zrcv)
                dnrZonesChim = []; listOfInterpDataTypes=[]
                for zdname in interDict[zrname]:
                    zd = Internal.getNodeFromName2(tpc, zdname)
                    isDnrCart=False
                    # base receveuse cartesienne : on prend les donneurs meme ds la meme base 
                    if isCartBaseR:
                        dnrZonesChim.append(zd)
                        prefix=zdname.split('.')
                        if len(prefix)>1:
                            prefix = prefix[0]
                            if prefix=='Cart': isDnrCart=True
                        if isDnrCart:
                            listOfInterpDataTypes.append(0)
                        else: 
                            listOfInterpDataTypes.append(1)
                    else:                        
                        if nobOfDnrBases[zdname]!=nob:
                            dnrZonesChim.append(zd)
                            prefix=zdname.split('.')
                            if len(prefix)>1:
                                prefix = prefix[0]
                                if prefix=='Cart': isDnrCart=True
                            if isDnrCart:
                                listOfInterpDataTypes.append(0)
                            else: 
                                listOfInterpDataTypes.append(1)
                #3- A OPTIMISER : des cartesiens et des ADT
                # for i in range(len(dnrZonesChim)):
                #     print(dnrZonesChim[i][0],listOfInterpDataTypes[i])
                X._setInterpData(zrcv, dnrZonesChim, nature=1, penalty=1, loc='centers', storage='inverse', 
                                 sameName=1, interpDataType=listOfInterpDataTypes, itype='chimera')

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

    #===============================
    # donnees de raccords match
    #===============================
    test.printMem(">>> Abutting data for curvilinear grids [start]")
    for nob in range(len(tp[2])):
        if Internal.getType(tp[2][nob])=='CGNSBase_t':
            basercv = tp[2][nob]
            basername = Internal.getName(basercv)

            # 1st step : Abutting IDs
            isCartBase=False
            if basername == 'CARTESIAN': isCartBase=True
            else:
                graphMatch = Cmpi.computeGraph(tp[2][nob], type='match', reduction=True)
                Cmpi._addXZones(tpc[2][nob], graphMatch, noCoordinates=True, cartesian=False)
                Cmpi._addXZones(tp[2][nob], graphMatch, noCoordinates=True, cartesian=False)
                # on construit les donnees des raccords match
                X._setInterpData(tp[2][nob], tpc[2][nob], nature=1, penalty=1, loc='centers', 
                                 storage='inverse', sameName=1, dim=dimPb, itype='abutting')
                Cmpi._rmXZones(tpc[2][nob])
                Cmpi._rmXZones(tp[2][nob])
                #C._rmBCOfType(tp[2][nob],'BCMatch')

    graphMatch={}
    test.printMem(">>> Abutting data [after free]")
    test.printMem(">>> Abutting data [end]")

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
    for i in range(Cmpi.size): datas[i] = [] # force

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
    #Internal._rmNodesByName(tpc, Internal.__GridCoordinates__)
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IBCDs = n[1]
            if IBCDs != []:
                zD = Internal.getNodeFromName2(tpc, zname)
                zD[2] += IBCDs

    datas = {}; graph = {}

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
    if isinstance(tc_out, str): 
        import Compressor.PyTree as Compressor
        Compressor._compressCartesian(tpc[2][1])
        Cmpi.convertPyTree2File(tpc, tc_out, ignoreProcNodes=True)
    
    # Initialisation
    if tinit is None: I._initConst(tp, loc='centers')
    else: 
        import Post.PyTree as Pmpi
        tp = Pmpi.extractMesh(tinit, tp, mode='accurate')
    if model != "Euler": C._initVars(tp, 'centers:ViscosityEddy', 0.)

    # Init with BBox : uniquement sur le cartesien autour des IBM - curviligne tel quel
    if initWithBBox>0.:
        print('Init momentum to 0 inside the bounding box')
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

    # clean t
    for z in Internal.getZones(tp):
        SDD = Internal.getNodeFromName1(z,".Solver#define")
        if SDD is not None:
            Internal._rmNodesFromName1(SDD,'inv')
            Internal._rmNodesFromName1(SDD,'ibctype')
            Internal._rmNodesFromName1(SDD,'dfar')
            Internal._rmNodesFromName1(SDD,'snear')

    # Save t
    if isinstance(t_out, str):
        #Compressor._compressCartesian(tp[2][1])
        Cmpi.convertPyTree2File(tp, t_out, ignoreProcNodes=True)

    if Cmpi.size > 1: Cmpi.barrier()    
    return tp, tpc

#====================================================================================
class IBMO(Common):
    """Preparation et calculs avec le module FastS."""
    def __init__(self, format=None, numb=None, numz=None):
        Common.__init__(self, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["stephanie@onera.fr", "ash@onera.fr"]
        self.cartesian = False
        
    # Prepare 
    def prepare(self, t_case, t_out, tc_out, distrib=False,
                vmin=21, check=False, frontType=1, NP=None, expand=3, dfarDir=0):
        if NP is None: NP = Cmpi.size
        if NP == 0: print('Preparing for a sequential computation.')
        else: print('Preparing for an IBMO computation on %d processors.'%NP)
        ret = prepare(t_case, t_out, tc_out, distrib=distrib,
                      vmin=vmin, check=check, NP=NP,  
                      frontType=frontType, expand=expand, dfarDir=dfarDir)
        return ret
