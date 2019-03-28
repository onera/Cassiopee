# Class for FastS "IBM" prepare and compute

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

# IN: maillage surfacique + reference State + snears

#================================================================================
# IBM prepare
# NP is the target number of processors
#================================================================================ 
def prepare(t_case, t_out, tc_out, snears=0.01, dfar=10., vmin=21, check=False, NP=0, format='single'):
    import Converter.Mpi as Cmpi
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    # sequential prep
    if rank == 0: ret = prepare0(t_case, t_out, tc_out, snears, dfar, vmin, check, NP, format)
    return ret

#================================================================================
# IBM prepare - seq
#================================================================================
def prepare0(t_case, t_out, tc_out, snears=0.01, dfar=10., vmin=21, check=False, NP=0, format='single'):

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    #-------------------------------------------------------
    # Refinement surfaces in the fluid
    #-------------------------------------------------------
    # snearsf: list of spacing required in the refinement surfaces
    snearsf = None    
    tbox = None
    # refinementSurfFile: surface meshes describing refinement zones
    #refinementSurfFile = 'refinementBody.cgns'
    #try: tbox = C.convertFile2PyTree(refinementSurfFile)
    #except: tbox=None # no refinement surface

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError, 'GoverningEquations is missing in input cgns.'
    # model: Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced

    #--------------------------------------------------------
    # Generates the full Cartesian mesh
    t = TIBM.generateIBMMesh(tb, vmin, snears, dfar, DEPTH=2,
                             tbox=tbox, snearsf=snearsf, check=check,
                             sizeMax=4000000)

    #------------------------------------------------------
    # distribute the mesh over NP processors
    if NP > 0:
        print('distribution over %d processors'%NP)
        stats = D2._distribute(t, NP)
        if check: print(stats)

    #------------------------------------------------
    # Add reference state to the pyTree and init fields
    # Add viscosity if model is not Euler
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)
    if check: C.convertPyTree2File(t, 'mesh1.cgns')

    #----------------------------------------
    # Computes distance field
    #----------------------------------------
    if dimPb == 2:
        z0 = Internal.getZones(t)
        bb = G.bbox(z0); dz = bb[5]-bb[2]
        tb2 = C.initVars(tb,'CoordinateZ',dz*0.5)
        DTW._distance2Walls(t,tb2,type='ortho',signed=0, dim=dimPb,loc='centers')
    else:
        DTW._distance2Walls(t,tb,type='ortho',signed=0, dim=dimPb,loc='centers')
    
    #----------------------------------------
    # Create IBM info
    #----------------------------------------
    t,tc = TIBM.prepareIBMData(t, tb, frontType=1, interpDataType=0)

    # arbre donneur
    D2._copyDistribution(tc, t)
    if isinstance(tc_out, str): Fast.save(tc, tc_out, split=format, NP=-NP)

    #----------------------------------------
    # Extraction des coordonnees des pts IBM
    #----------------------------------------
    if check:
        tibm = TIBM.extractIBMInfo(tc)
        C.convertPyTree2File(tibm, 'IBMInfo.cgns')
        del tibm

    # arbre de calcul
    I._initConst(t, loc='centers')
    if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)
    if isinstance(t_out, str): Fast.save(t, t_out, split=format, NP=-NP)
    return t, tc

#==================================================================================================
def prepare1(t_case, t_out, tc_out, snears=0.01, dfar=10., vmin=21, check=False, NP=0, format='single'):
    import Generator
    import Converter
    import Connector.connector as connector
    import Connector.Mpi as Xmpi
    import Post.PyTree as P
    import Converter.Distributed as Distributed
    import Connector.OversetData as XOD
    import KCore.test as test
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    DEPTH=2
    IBCType=1
    frontType=1

    rank = Cmpi.rank

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError, 'GoverningEquations is missing in input tree defined in %s.'%FILE
    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced

    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')
    surfaces = Internal.getZones(tb)
    # Build octree
    i = 0; surfaces=[]; snearso=[] # pas d'espace sur l'octree
    bodies = Internal.getZones(tb)
    if not isinstance(snears, list): snears = len(bodies)*[snears]
    if len(bodies) != len(snears):
        raise ValueError('generateIBMMesh: Number of bodies is not equal to the size of snears.')
    dxmin0 = 1.e10
    for s in bodies:
        sdd = Internal.getNodeFromName1(s, ".Solver#define")
        if sdd is not None:
            snearl = Internal.getNodeFromName1(sdd, "snear")
            if snearl is not None: 
                snearl = Internal.getValue(snearl)
                snears[i] = snearl
        dhloc = snears[i]*(vmin-1)
        surfaces += [s]; snearso += [dhloc]
        dxmin0 = min(dxmin0,dhloc)
        i += 1
    o = G.octree(surfaces, snearso, dfar=dfar, balancing=2)
    vmint = 31

    if vmin < vmint:
        #if rank==0: print('generateIBMMesh: octree finest level expanded (expandLayer activated).')
        to = C.newPyTree(['Base',o])
        to = TIBM.blankByIBCBodies(to, tb, 'centers', dimPb)
        C._initVars(o,"centers:indicator", 0.)
        cellN = C.getField("centers:cellN",to)[0]
        octreeA = C.getFields(Internal.__GridCoordinates__, o)[0]
        indic = C.getField("centers:indicator",o)[0]
        indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic,0,0)
        indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic,1,0) # CB
        indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic,2,0) # CB
        indic = Generator.generator.modifyIndicToExpandLayer(octreeA, indic,3,0) # CB
                                                                                      
        indic = Converter.addVars([indic,cellN])
        indic = Converter.initVars(indic,"{indicator}={indicator}*({cellN}>0.)")
        octreeA = Generator.adaptOctree(octreeA, indic, balancing=2)
        o = C.convertArrays2ZoneNode(o[0],[octreeA])

        to = C.newPyTree(['Base', o])
        G._getVolumeMap(to); volmin = C.getMinValue(to, 'centers:vol')

        dxmin = (volmin)**(1./dimPb)
        if rank==0: print('Minimum spacing of Cartesian mesh= %f (targeted %f)'%(dxmin/(vmin-1),dxmin0/(vmin-1)))
        C._rmVars(o,'centers:vol')

    test.printMem(">>> Octree unstruct [end]")

    # Split octree
    test.printMem(">>> Octree unstruct split [start]")
    bb = G.bbox(o)
    NPI = Cmpi.size
    if NPI == 1: p = Internal.copyRef(o) # keep reference
    else: p = T.splitNParts(o, N=NPI, recoverBC=False)[rank]
    del o
    if check: C.convertPyTree2File(p, 'octree_%d.cgns'%rank)
    test.printMem(">>> Octree unstruct split [end]")
    
    # fill vmin + merge in parallel
    test.printMem(">>> Octree struct [start]")
    res = TIBM.octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, merged=1, sizeMax=4000000); del p
    t = C.newPyTree(['CARTESIAN', res])
    zones = Internal.getZones(t)
    for z in zones: z[0] = z[0]+'_proc%d'%rank
    Cmpi._setProc(t,rank)
    C._addState(t, 'EquationDimension', dimPb)
    test.printMem(">>> Octree struct [end]")
    
    # Add xzones for ext
    test.printMem(">>> extended cart grids (XZones) [start]")
    tbb = Cmpi.createBBoxTree(t)
    interDict = X.getIntersectingDomains(tbb)
    graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(t, graph)
    zones = Internal.getZones(t)
    coords = C.getFields(Internal.__GridCoordinates__, zones, api=2)
    coords = Generator.generator.extendCartGrids(coords, DEPTH+1, 1)
    C.setFields(coords, zones, 'nodes')
    Cmpi._rmXZones(t)
    test.printMem(">>> extended cart grids (XZones) [end]")
    
    TIBM._addBCOverlaps(t, bbox=bb)
    TIBM._addExternalBCs(t, bbox=bb, dimPb=dimPb)

    dz = 0.01
    if dimPb == 2:
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)

    # ReferenceState 
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)

    # Distance a la paroi    
    test.printMem(">>> Wall distance [start]")
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = bb[5]-bb[2]
        dz = bb0[5]-bb0[2]
        tb2 = C.initVars(tb,'CoordinateZ',dz*0.5)
        DTW._distance2Walls(t,tb2,type='ortho',signed=0, dim=dimPb,loc='centers')
    else:
        DTW._distance2Walls(t,tb,type='ortho',signed=0, dim=dimPb,loc='centers')

    X._applyBCOverlaps(t, depth=DEPTH, loc='centers', val=2, cellNName='cellN')
    C._initVars(t,'{centers:cellNChim}={centers:cellN}')

    C._initVars(t,'centers:cellN',1.)
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, 'Zone_t')
        dims = Internal.getZoneDim(z0)
        npts = dims[1]*dims[2]*dims[3]
        zmin = C.getValue(z0,'CoordinateZ',0)
        zmax = C.getValue(z0,'CoordinateZ',npts-1)
        dz = zmax-zmin
        # Creation du corps 2D pour le preprocessing IBC
        T._addkplane(tb)
        T._contract(tb, (0,0,0), (1,0,0), (0,1,0), dz)
    test.printMem(">>> Wall distance [end]")
    
    test.printMem(">>> Blanking [start]")
    t = TIBM.blankByIBCBodies(t, tb, 'centers', dimPb)
    C._initVars(t, '{centers:cellNIBC}={centers:cellN}')
    TIBM._signDistance(t)

    C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
    # determination des pts IBC
    if IBCType == -1: X._setHoleInterpolatedPoints(t,depth=-DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
    elif IBCType == 1:
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellN',addGC=False) # pour les gradients
        if frontType < 2:
            X._setHoleInterpolatedPoints(t,depth=DEPTH,dir=0,loc='centers',cellNName='cellN',addGC=False)
        else:
            DEPTHL=DEPTH+1
            X._setHoleInterpolatedPoints(t,depth=DEPTHL,dir=0,loc='centers',cellNName='cellN',addGC=False)         
            #cree des pts extrapoles supplementaires
            TIBM._blankClosestTargetCells(t,cellNName='cellN',depth=DEPTHL)
    else:
        raise ValueError('prepareIBMData: not valid IBCType. Check model.')
    TIBM._removeBlankedGrids(t, loc='centers')
    test.printMem(">>> Blanking [end]")
    
    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
        dims = Internal.getZoneDim(i)
        npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))

    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')

    if IBCType==-1:
        #print('Points IBC interieurs: on repousse le front un peu plus loin.')
        C._initVars(t,'{centers:cellNDummy}=({centers:cellNIBC}>0.5)*({centers:cellNIBC}<1.5)')
        X._setHoleInterpolatedPoints(t,depth=1,dir=1,loc='centers',cellNName='cellNDummy',addGC=False)
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNDummy}>0.5, {centers:cellNDummy}<1.5)')
        C._rmVars(t, ['centers:cellNDummy'])
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
    else:
        C._initVars(t,'{centers:cellNFront}=logical_and({centers:cellNIBC}>0.5, {centers:cellNIBC}<1.5)')
        for z in Internal.getZones(t):
            connector._updateNatureForIBM(z, IBCType,
                                          Internal.__GridCoordinates__,
                                          Internal.__FlowSolutionNodes__,
                                          Internal.__FlowSolutionCenters__)
            
    # setInterpData - Chimere
    C._initVars(t,'{centers:cellN}=maximum(0.,{centers:cellNChim})')# vaut -3, 0, 1, 2 initialement

    # maillage donneur: on MET les pts IBC comme donneurs
    tp = Internal.copyRef(t)
    FSN = Internal.getNodesFromName3(tp, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(FSN, 'cellNFront')
    Internal._rmNodesByName(FSN, 'cellNIBC')
    Internal._rmNodesByName(FSN, 'TurbulentDistance')

    test.printMem(">>> Interpdata [start]")
    tc = C.node2Center(tp)
    
    # setInterpData parallel pour le chimere
    tbbc = Cmpi.createBBoxTree(tc)
    interDict = X.getIntersectingDomains(tbbc)
    graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(tc, graph)
    test.printMem(">>> Interpdata [after addDXzones]")
    
    procDict = Cmpi.getProcDict(tc)
    datas={}
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
        
            allIDs = Internal.getNodesFromName(zd, 'ID*')
            IDs = []              
            for zsr in allIDs:
                if Internal.getValue(zsr)==zrname: IDs.append(zsr)
            if IDs != []:
                if destProc == rank:                
                    zD = Internal.getNodeFromName2(tc,zdname)
                    zD[2] += IDs
                else:
                    if destProc not in datas: datas[destProc] = [[zdname,IDs]]
                    else: datas[destProc] += [[zdname,IDs]]
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

    # fin interpData
    C._initVars(t,'{centers:cellNIBCDnr}=minimum(2.,abs({centers:cellNIBC}))')
    C._initVars(t,'{centers:cellNIBC}=maximum(0.,{centers:cellNIBC})')# vaut -3, 0, 1, 2, 3 initialement
    C._initVars(t,'{centers:cellNIBC}={centers:cellNIBC}*({centers:cellNIBC}<2.5)')    
    C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
    C._cpVars(t,'centers:cellN',tc,'cellN')
    test.printMem(">>> Interpdata [end]")
    
    # Transfert du cellNFront
    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

    # propager cellNVariable='cellNFront'
    Xmpi._setInterpTransfers(t,tc,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)

    C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
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
    test.printMem(">>> Building IBM front [end]")

    # Interpolation IBC (front, tbbc)

    # graph d'intersection des pts images de ce proc et des zones de tbbc
    graph = {}
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
            for nozr in xrange(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    zrname = zonesRIBC[nozr][0]
                    interpPtsBB = Generator.BB(allInterpPts[nozr])
                    for z in zones:
                        bba = C.getFields('GridCoordinates', z)[0]
                        if Generator.bboxIntersection(interpPtsBB,bba,isBB=True):
                            popp = Cmpi.getProc(z)
                            Distributed.updateGraph__(graph, popp, rank, z[0])
                            if zrname not in interDictIBM: interDictIBM[zrname]=[z[0]]
                            else: interDictIBM[zrname].append(z[0])
    else: graph={}

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

    Cmpi._addXZones(tc, graph)
    test.printMem(">>> Interpolating IBM [start]")

    ReferenceState = Internal.getNodeFromType2(t, 'ReferenceState_t')
    nbZonesIBC = len(zonesRIBC)

    if dictOfCorrectedPtsByIBCType!={}:
        for ibcTypeL in dictOfCorrectedPtsByIBCType:
            allCorrectedPts = dictOfCorrectedPtsByIBCType[ibcTypeL]
            allWallPts = dictOfWallPtsByIBCType[ibcTypeL]
            allInterpPts = dictOfInterpPtsByIBCType[ibcTypeL]
            for nozr in xrange(nbZonesIBC):
                if allCorrectedPts[nozr] != []:
                    zrcv = zonesRIBC[nozr]
                    zrname = zrcv[0]
                    for zdname in interDictIBM[zrname]:
                        zd = Internal.getNodeFromName2(tc, zdname)
                        #if zd is not None: dnrZones.append(zd)
                        if zd is None: print '!!!Zone None', zrname, zdname
                        else: dnrZones.append(zd)

                    XOD._setIBCDataForZone__(zrcv,dnrZones,allCorrectedPts[nozr],allWallPts[nozr],allInterpPts[nozr],
                                             nature=1,penalty=1,loc='centers',storage='inverse', dim=dimPb,
                                             interpDataType=0,ReferenceState=ReferenceState,bcType=ibcTypeL)

                    nozr += 1
                    for zd in dnrZones:       
                        zdname = zd[0]
                        destProc = procDict[zdname]
                        allIDs = Internal.getNodesFromName(zd,'IBCD*')
                        IDs = []      
                        for zsr in allIDs:
                            if Internal.getValue(zsr)==zrname: IDs.append(zsr)
                        if IDs != []:
                            if destProc == rank:                
                                zD = Internal.getNodeFromName2(tc,zdname)
                                zD[2] += IDs
                            else:
                                if destProc not in datas: datas[destProc]=[[zdname,IDs]]
                                else: datas[destProc] += [[zdname,IDs]]

    Cmpi._rmXZones(tc)
    test.printMem(">>> Interpolating IBM [end]")

    Internal._rmNodesByName(tc, Internal.__FlowSolutionNodes__)
    Internal._rmNodesByName(tc, Internal.__GridCoordinates__)
    destDatas = Cmpi.sendRecv(datas, graph)
    for i in destDatas:
        for n in destDatas[i]:
            zname = n[0]
            IBCDs = n[1]
            if IBCDs != []:
                zD = Internal.getNodeFromName2(tc, zname)
                zD[2] += IBCDs

    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:cellNChim','centers:cellNIBCDnr']
    if model == 'Euler': varsRM += ['centers:TurbulentDistance']
    C._rmVars(t, varsRM)
    test.printMem(">>> Saving [start]")

    # ne marche pas pour l'instant en parallele
    if check: 
        tibm = TIBM.extractIBMInfo(tc)
        Cmpi.convertPyTree2File(tibm, 'IBMInfo.cgns')

    if isinstance(tc_out, str): Cmpi.convertPyTree2File(tc, tc_out)

    I._initConst(t, loc='centers')
    if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)
    if isinstance(t_out, str): Cmpi.convertPyTree2File(t, t_out)
    return t, tc

#=============================================================================
# Post
#==============================================================================
def post(t_case, t_in, tc_in, t_out, wall_out, NP=0, format='single'):
    import Post.PyTree as P
    
    if isinstance(t_in, str): t = C.convertFile2PyTree(t_in)
    else: t = t_in
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    #=============================
    # Supprime les champs inutiles
    #=============================
    vars = ['centers:Density_M1', 'centers:VelocityX_M1', 'centers:VelocityY_M1', 'centers:VelocityZ_M1', 'centers:Temperature_M1', 'centers:Density_P1', 'centers:VelocityX_P1', 'centers:VelocityY_P1', 'centers:VelocityZ_P1', 'centers:Temperature_P1','centers:TurbulentDistance']
    C._rmVars(t, vars)

    #=============================
    # Arbre de connectivite
    #=============================
    if isinstance(tc_in, str): tc = C.convertFile2PyTree(tc_in)
    else: tc = tc_in
    
    Internal._rmNodesByName(tc, 'GridCoordinates')

    #==========================================================
    # Extraction Cp, Cf, ... sur les surfaces par interpolation
    #==========================================================
    tb = C.convertArray2Tetra(tb)

    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError, 'GoverningEquations is missing in input cgns.'
    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    if model == 'Euler': bcType = 0
    elif model =='NSLaminar': bcType = 1
    else: bcType = 3 # Musker

    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
    ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
    Mus, Cs, Ts, Pr] = C.getState(tb)

    varType = 2 # IBM updated variables (rho,u,t)
    varsIBC = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
    vars = ['Density', 'VelocityX', 'VelocityY', 'VelocityZ', 'Temperature']
    if model != 'Euler': 
        vars += ['ViscosityEddy']
        if model == 'NSTurbulent': 
            vars += ['TurbulentSANuTilde']
            varsIBC += ['TurbulentSANuTilde']
            varType = 21

    for z in Internal.getNodesFromType2(t,"Zone_t"):
        zc = Internal.getNodeFromName(tc,z[0])
        for v in varsIBC: C._cpVars(z, 'centers:'+v, zc, v)

    X._setInterpTransfers(t, tc, variables=vars,
                          variablesIBC=varsIBC, bcType=bcType,
                          varType=varType, storage=1,
                          Gamma=Gamma, Cv=cvInf, MuS=Mus,
                          Cs=Cs, Ts=Ts)
    zw = TIBM.extractIBMWallFields(tc, tb=tb)

    RoUInf2I = 1./(RouInf*RouInf+RovInf*RovInf+RowInf*RowInf)
    C._initVars(zw,'{Cp}=2*%f*({Pressure}-%f)*%f'%(RoInf,PInf,RoUInf2I))
    if model != 'Euler':
        C._initVars(zw,'{Cf}=2*%f*{Density}*{utau}**2*%f'%(RoInf,RoUInf2I))

    Internal._rmNodesByName(zw, '.Solver#Param')
    Internal._rmNodesByName(zw, '.Solver#ownData')

    if isinstance(wall_out, str): C.convertPyTree2File(zw, wall_out)

    #===============================
    # En 2D, extrait un seul plan k
    #================================
    if dimPb == 2:
        t = T.subzone(t, (1,1,1), (-1,-1,1))
        C._initVars(tb, 'CoordinateZ', 0.) # forced

    #=================================
    # Calcul de mut/mu dans le volume
    #=================================
    if model != 'Euler':
        betas = Mus*(Ts+Cs)/(Ts**(3./2.))
        C._initVars(t,'{centers:ViscosityMolecular} = %20.16g*sqrt({centers:Temperature})/(1.+%20.16g/{centers:Temperature})'%(betas,Cs))
        C._initVars(t,'{centers:mutsmu}=({centers:ViscosityEddy})/({centers:ViscosityMolecular})-1.')

    #==============================
    # Sortie champs aux noeuds
    #==============================
    vars = ['centers:Density','centers:VelocityX', 'centers:Temperature','centers:ViscosityEddy', 
    'centers:TurbulentSANuTilde','centers:ViscosityMolecular', 'centers:mutsmu', 'centers:cellN']
    for v in vars: t = C.center2Node(t, v)
    Internal._rmNodesByName(t, 'FlowSolution#Centers')
    if isinstance(t_out, str): C.convertPyTree2File(t, t_out)

    return t, zw

#====================================================================================
# Redistrib on N processors
#====================================================================================
def distribute(t_in, tc_in, NP):
    if isinstance(tc_in, str):
        tcs = Cmpi.convertFile2SkeletonTree(tc_in, maxDepth=3)
    else: tcs = tc_in
    stats = D2._distribute(tcs, NP, algorithm='graph', useCom='ID') 
    print stats
    if isinstance(tc_in, str):
        paths = []; ns = []
        bases = Internal.getBases(tcs)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                nodes = Internal.getNodesFromName2(z, 'proc')
                for n in nodes:
                    p = 'CGNSTree/%s/%s/.Solver#Param/proc'%(b[0],z[0])
                    paths.append(p); ns.append(n)
        Filter.writeNodesFromPaths(tc, paths, ns, maxDepth=0, mode=1)

    if isinstance(t_in, str):
        ts = Cmpi.convertFile2SkeletonTree(t, maxDepth=3)
    else: ts = t_in
    D2._copyDistribution(ts, tcs)

    if isinstance(t_in, str):
        paths = []; ns = []
        bases = Internal.getBases(ts)
        for b in bases:
            zones = Internal.getZones(b)
            for z in zones:
                nodes = Internal.getNodesFromName2(z, 'proc')
                for n in nodes:
                    p = 'CGNSTree/%s/%s/.Solver#Param/proc'%(b[0],z[0])
                    paths.append(p); ns.append(n)
        Filter.writeNodesFromPaths(t, paths, ns, maxDepth=0, mode=1)

    # Affichage du nombre de points par proc - equilibrage ou pas
    NptsTot = 0
    for i in xrange(NP):
        NPTS = 0
        for z in Internal.getZones(ts):
            if Cmpi.getProc(z)==i: NPTS += C.getNPts(z)
        NptsTot += NPTS
        print 'Rank %d has %d points'%(i,NPTS) 
    print 'All points: %d million points'%(NptsTot/1.e6)
    return ts, tcs

#====================================================================================
class IBM(Common):
    """Preparation et caculs avec le module FastS."""
    def __init__(self, NP=None, format=None, numb=None, numz=None):
        Common.__init__(self, NP, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["ash@onera.fr"]
        
    # Prepare : n'utilise qu'un proc pour l'instant
    def prepare(self, t_case, t_out, tc_out, snears=0.01, dfar=10., vmin=21, check=False):
        NP = self.data['NP']
        if NP == 0: print('Preparing for a sequential computation.')
        else: print('Preparing for a computation on %d processors.'%NP)
        ret = prepare(t_case, t_out, tc_out, snears, dfar, vmin, check, NP, self.data['format'])
        return ret

    # post-processing: extrait la solution aux noeuds + le champs sur les surfaces
    def post(self, t_case, t_in, tc_in, t_out, wall_out):
        return post(t_case, t_in, tc_in, t_out, wall_out, self.data['NP'], self.data['format'])
