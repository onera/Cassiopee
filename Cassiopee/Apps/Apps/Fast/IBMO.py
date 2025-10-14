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
import RigidMotion.PyTree as RM
from . import IBM as AppIBM
rank = Cmpi.rank
comm = Cmpi.COMM_WORLD
DEPTH = 2
dz = 0.01   # for 2D cases

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

def buildAllBodies__(t_case, distrib, motion=False):
    if isinstance(t_case, str):
        h = Filter.Handle(t_case)
        if distrib:
            tb = h.loadAndSplit(NParts=Cmpi.size)
        else:
            tb = h.loadAndSplit(NParts=max(1, Cmpi.size))
        #tb = Cmpi.convertFile2SkeletonTree(t_case)
    else: tb = t_case

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    # tbchim : arbre des corps chimeres
    tbchim = C.newPyTree()
    tbov = C.newPyTree()
    tbblank = C.newPyTree()
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
                    now = 1
                    for zw in zwalls:
                        zw[0] = b[0]+'_'+z[0]+'_W'+str(now)
                        if TM is not None: zw[2].append(TM)
                        if SD is not None: zw[2].append(SD)
                        if SP is not None: zw[2].append(SP)
                        now += 1
                    basechim[2] += zwalls

                    zovs = C.extractBCOfType(z, 'BCOverlap')
                    nov = 1
                    for zo in zovs:
                        zo[0] = b[0]+'_'+z[0]+'_OV'+str(nov)
                        if TM is not None: zo[2].append(TM)
                        if SD is not None: zo[2].append(SD)
                        if SP is not None: zo[2].append(SP)
                        nov += 1
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

    # for distance field at cell centers
    tball = C.newPyTree(['AllBodies'])
    tball[2][1][2] = Internal.getZones(tbchim) + Internal.getZones(tbibm)
    if dimPb == 2:
        ZMEAN = dz*0.5
        C._initVars(tball, 'CoordinateZ', ZMEAN)

    if motion:
        RM._copyGrid2GridInit(tbov)
        RM._copyGrid2GridInit(tbchim)

    return [tb, tbov, tbchim, tbibm, tball, baseNamesChim]

def createChimeraTree__(tb, tbchim, tball, baseNamesChim, dimPb=3, motion=False):
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
        t2 = X.connectMatch(t2, tol=1.e-6, dim=dimPb)
    Internal._addGhostCells(t2, t2, DEPTH, adaptBCs=1, fillCorner=0)

    # Suppression des XZones et correction des matchs
    if Cmpi.size > 1:Cmpi._rmBXZones(t2)
    # Suppression des XZones et correction des matchs
    if Cmpi.size > 1:Cmpi._rmBXZones(t2)

    # Fusion des fenetres des raccords
    if Cmpi.size > 1: t2 = Xmpi.mergeWindows(t2)

    if dimPb == 2:
        T._addkplane(t2)
        T._contract(t2, (0,0,0), (1,0,0), (0,1,0), dz)
        T._makeDirect(t2)

    # Distances a la paroi sur les grilles curvilignes chimere
    if motion:
        DTW._distance2Walls(t2, tbchim, loc='centers', type='ortho')
    else:
        DTW._distance2Walls(t2, tball, loc='centers', type='ortho')
    test.printMem(">>> Wall distance on Chimera grids [end]")

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

    return t2, tc2

def prepareMotion(t_case, t_out, tc_out, t_in=None, to=None, snears=0.01, dfars=10.,
                  tbox=None, snearsf=None, yplus=100., Lref=1.,
                  vmin=21, check=False, NP=0, format='single',
                  frontType=1, extrusion=False, smoothing=False, balancing=False,
                  distrib=True, expand=3, tinit=None, initWithBBox=-1., wallAdapt=None, yplusAdapt=100.,
                  dfarDir=0,dz_in=0.01,span_in=0.25,NPas_in=10,height_in=0.1,
                  correctionMultiCorpsF42=False, blankingF42=False, twoFronts=False,
                  redistribute=False,isDist2WallNearBodyOnly=False,isoverideheight=False,check2Donly=False,
                  dict_Nz={},isCartesianExtrude=False,isExtrudeByZone=False,directory_tmp_files='./'):
    motion = True
    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD
    tolMatch = 1.e-6
    DEPTH = 2
    IBCType = 1

    [tb, tbov, tbchim, tbibm, tball, baseNamesChim] = buildAllBodies__(t_case, distrib, motion)

    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    # Reference state
    refstate = C.getState(tb)

    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree.')
    model = Internal.getValue(model)

    # Construction de l'arbre des maillages chimeres
    t2,tc2 = createChimeraTree__(tb, tbchim, tball, baseNamesChim, dimPb=dimPb, motion=motion)
    # construction de l'arbre des corps pour l'octree : tbibm
    # if bodies in motion taken into account
    # either by a tbox or  input octree to defined
    t,tc = AppIBM.prepare1(tbibm, None, None, t_in=None, to=to,
                           snears=snears, dfars=dfars,
                           tbox=tbox, snearsf=snearsf, yplus=yplus, Lref=Lref,
                           vmin=vmin, check=check, NP=NP, format='single',
                           frontType=1, extrusion=extrusion, smoothing=smoothing, balancing=balancing,
                           distrib=distrib, expand=expand, tinit=tinit, initWithBBox=initWithBBox, wallAdapt=None, yplusAdapt=100.,
                           dfarDir=dfarDir,dz_in=dz_in,span_in=span_in,NPas_in=NPas_in,height_in=height_in,
                           correctionMultiCorpsF42=correctionMultiCorpsF42, blankingF42=blankingF42, twoFronts=twoFronts,
                           redistribute=redistribute,
                           isDist2WallNearBodyOnly=isDist2WallNearBodyOnly,
                           isoverideheight=isoverideheight,
                           check2Donly=check2Donly,
                           dict_Nz=dict_Nz,isCartesianExtrude=isCartesianExtrude,isExtrudeByZone=isExtrudeByZone,
                           directory_tmp_files=directory_tmp_files)

    # merge trees
    tp1 = Internal.copyRef(t)
    tp2 = Internal.copyRef(t2)
    tp = Internal.merge([tp1,tp2])
    tpc1 = Internal.copyRef(tc)
    tpc2 = Internal.copyRef(tc2)
    tpc = Internal.merge([tpc1,tpc2])
    test.printMem(">>> Saving [start]")

    # distribution on NP processors
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

    if motion: RM._copyGrid2GridInit(tpc)

    Internal._rmNodesFromName(tpc,'ZoneRind')
    if isinstance(tc_out, str):
        #import Compressor.PyTree as Compressor
        #Compressor._compressCartesian(tpc[2][1])
        Cmpi.convertPyTree2File(tpc, tc_out, ignoreProcNodes=True)

    # Initialization
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
    if motion: RM._copyGrid2GridInit(tp)

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
def prepare(t_case, t_out, tc_out, tblank=None, to=None,
            vmin=21, check=False, NP=0,
            frontType=1,  tbox=None, snearsf=None,
            expand=3, distrib=False, tinit=None, initWithBBox=-1.,
            dfarDir=0):
    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD
    tolMatch = 1.e-6
    DEPTH = 2
    IBCType = 1
    if frontType!=1:
        print("Warning: IBMO.prepare: currently implemented for frontType=1 algorithm only.")
        frontType=1

    #read blanking bodies - if None : blanking bodies for Chimera are bcwalls
    tbblank = None
    if tblank is not None:
        if isinstance(tblank,str):
            tbblank = C.convertFile2PyTree(tblank)
        else: tbblank = tblank

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
    if to is not None:
        if isinstance(to, str):
            o = C.convertFile2PyTree(to)
            o = Internal.getZones(o)[0]
        else:
            o = Internal.getZones(to)[0]
        parento = None
    else:
        # construction de l'arbre des corps pour l'octree : tbo
        tbo = C.newPyTree(['Surf_octree'])
        tbo[2][1][2] = Internal.getZones(tbibm)+Internal.getZones(tbov)

        # Extraction de la liste des dfars de tb
        zones = Internal.getZones(tbo)
        dfars = [10.]*len(zones)
        for c, z in enumerate(zones):
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None: dfars[c] = Internal.getValue(n)*1.

        o = TIBM.buildOctree(tbo, snearFactor=1., dfars=dfars,
                             dimPb=dimPb, vmin=vmin,
                             expand=expand, dfarDir=dfarDir)

        # build parent octree 3 levels higher
        # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
        parento = TIBM.buildParentOctrees__(o, tbo, snearFactor=4., dfars=dfars,
                                            dimPb=dimPb, vmin=vmin)
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
    coords = C.getFields(Internal.__GridCoordinates__, zones, api=3)
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
    if tbblank is None:
        t  = TIBM.blankByIBCBodies(t, tbchim, 'centers', dimPb, cellNName='cellNChim')
    else:
        t = TIBM.blankByIBCBodies(t, tbblank, 'centers', dimPb, cellNName='cellNChim')
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
    P._computeGrad2(tp_cart, 'centers:TurbulentDistance', withCellN=False)

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
                        bba = C.getFields('GridCoordinates', z, api=1)[0]
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
    if Cmpi.KCOMM is not None: allGraph = Cmpi.KCOMM.allgather(graph)
    else: allGraph = [graph]
    #if rank == 0: print(allGraph)

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
    model = Internal.getNodeFromName(t, 'GoverningEquations')
    if model is not None: model = Internal.getValue(model)
    else:                 model = "Euler"

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
                                             interpDataType=0, ReferenceState=ReferenceState, bcType=ibcTypeL,model=model)

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
