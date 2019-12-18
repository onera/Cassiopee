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
from Apps.Fast.Common import Common
import Generator
import KCore.test as test

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
            frontType=1, expand=3):
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
    
    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')

    # constuction de l'arbre des corps pour l'octree : tbo
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
        T._contract(t2, (0,0,0), (1,0,0), (0,1,0), 0.01)
        T._makeDirect(t2)

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

    # tbibm : arbre des corps ibms
    tbibm = C.newPyTree()
    bases = Internal.getBases(tb)
    for b in bases:
        zones = Internal.getZones(b)
        isChimera = False
        for z in zones:
            if isZoneChimera(z): isChimera = True; break
        if not isChimera: tbibm[2].append(b) 

    # tball : arbre de tous les corps
    tball = C.newPyTree()
    tball[2] += Internal.getBases(tbchim) + Internal.getBases(tbibm)

    # Distances a la paroi sur les grilles chimeres
    DTW._distance2Walls(t2, tball, loc='centers', type='ortho')

    # ReferenceState 
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)

    # Distance a la paroi par rapport aux corps IBMs
    test.printMem(">>> Wall distance [start]")
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]
        tball2 = C.initVars(tball,'CoordinateZ',dz*0.5)
        DTW._distance2Walls(t, tball2, type='ortho', signed=0, dim=dimPb, loc='centers')
    else:
        DTW._distance2Walls(t, tball, type='ortho', signed=0, dim=dimPb, loc='centers')

    # Calcul du cellNChim sur CART
    C._initVars(t, 'centers:cellNChim', 1.)
    X._applyBCOverlaps(t, depth=2, loc='centers', val=2, cellNName='cellNChim')
        
    # Calcul du cellN sur Chimere
    C._initVars(t2, 'centers:cellNChim', 1.)
    X._applyBCOverlaps(t2, depth=4, val=2, cellNName='cellNChim')
    X._applyBCOverlaps(t2, depth=2, val=0, cellNName='cellNChim')

    # Masquage de CART par chimere
    test.printMem(">>> Blanking [start]")
    t = TIBM.blankByIBCBodies(t, tbchim, 'centers', dimPb, cellNName='cellNChim')
    t2 = TIBM.blankByIBCBodies(t2, tbibm, 'centers', dimPb, cellNName='cellNChim')
    # a faire : masquer les bases chimeres entre elles

    # Masquage IBC
    C._initVars(t, 'centers:cellNIBC', 1.)
    t = TIBM.blankByIBCBodies(t, tbibm, 'centers', dimPb, cellNName='cellNIBC')

    # Signe la distance en fonction de cellNIBC et cellNChim
    TIBM._signDistance(t)

    # Point interpole autour des points masques
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
    tc = C.node2Center(tp1) # cartesien
    tc2 = C.node2Center(tp2) # chimere

    # arbre complet tp
    tp = Internal.merge(tp1, tp2)
    # tpc arbre complet des centres
    tpc = Internal.merge(tc, tc2)

    FSN = Internal.getNodesFromName3(tp, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(FSN, 'cellNFront')
    Internal._rmNodesByName(FSN, 'cellNIBC')
    Internal._rmNodesByName(FSN, 'TurbulentDistance')

    test.printMem(">>> Interpdata [start]")
    
    # setInterpData parallel pour le chimere
    tbbc = Cmpi.createBBoxTree(tpc)
    interDict = X.getIntersectingDomains(tbbc)
    graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(tpc, graph, variables=['cellN'])
    test.printMem(">>> Interpdata [after addXZones]")
    
    procDict = Cmpi.getProcDict(tpc)
    datas = {}
    for zrcv in Internal.getZones(tp):
        zrname = zrcv[0]
        dnrZones = []
        for zdname in interDict[zrname]:
            zd = Internal.getNodeFromName2(tpc, zdname)
            dnrZones.append(zd)
        X._setInterpData(zrcv, dnrZones, nature=1, penalty=1, loc='centers', storage='inverse', 
                         sameName=1, interpDataType=1, itype='all')
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
                    zD = Internal.getNodeFromName2(tpc, zdname)
                    zD[2] += IDs
                else:
                    if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                    else: datas[destProc].append([zdname,IDs])
            else:
                if destProc not in datas: datas[destProc] = []
    Cmpi._rmXZones(tpc)
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

    # Traitement IBC
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
        res = TIBM.getAllIBMPoints(zonesRIBC, loc='centers',tb=tb, tfront=front, frontType=frontType,
                                   cellNName='cellNIBC', depth=DEPTH, IBCType=IBCType)

    # cleaning
    C._rmVars(tc,['cellNChim','cellNIBC','TurbulentDistance','cellNFront'])
    # dans t, il faut cellNChim et cellNIBCDnr pour recalculer le cellN a la fin
    varsRM = ['centers:gradxTurbulentDistance','centers:gradyTurbulentDistance','centers:gradzTurbulentDistance','centers:cellNFront','centers:cellNIBC']
    C._rmVars(t, varsRM)
    front = None
    test.printMem(">>> Building IBM front [end]")

    # a finir


    return t

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
