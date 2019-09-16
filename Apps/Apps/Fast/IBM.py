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

try: range = xrange
except: pass

# IN: maillage surfacique + reference State + snears

#================================================================================
# IBM prepare
# NP is the target number of processors for computation 
# (maybe different from the number of processors the prep is run on)
#================================================================================ 
def prepare(t_case, t_out, tc_out, snears=0.01, dfar=10., dfarList=[],
            tbox=None, snearsf=None,            
            vmin=21, check=False, NP=0, format='single',
            frontType=1, expand=3):
    import Converter.Mpi as Cmpi
    rank = Cmpi.rank; size = Cmpi.size
    ret = None
    # sequential prep
    if size == 1: ret = prepare0(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList, 
                                 tbox=tbox, snearsf=snearsf,
                                 vmin=vmin, check=check, NP=NP, format=format, frontType=frontType,
                                 expand=expand)
    # parallel prep
    else: ret = prepare1(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList, 
                         tbox=tbox, snearsf=snearsf,
                         vmin=vmin, check=check, NP=NP, format=format, frontType=frontType,
                         expand=expand)
    
    return ret

#================================================================================
# IBM prepare - seq
#================================================================================
def prepare0(t_case, t_out, tc_out, snears=0.01, dfar=10., dfarList=[],
             tbox=None, snearsf=None,
             vmin=21, check=False, NP=0, format='single',
             frontType=1, expand=3):
    import KCore.test as test
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case
    
    # list of dfars
    if dfarList == []:
        zones = Internal.getZones(tb)
        dfarList = [dfar*1.]*len(zones)
        for c, z in enumerate(zones): 
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None:
                dfarList[c] = Internal.getValue(n)*1.

    #-------------------------------------------------------
    # Refinement surfaces in the fluid
    #-------------------------------------------------------
    # snearsf: list of spacing required in the refinement surfaces
    if tbox is not None:
        if isinstance(tbox, str): tbox = C.convertFile2PyTree(tbox)
        else: tbox = tbox
        if snearsf is None:
            snearsf = []
            zones = Internal.getZones(tbox)
            for z in zones:
                sn = Internal.getNodeFromName2(z, 'snear')
                if sn is not None: snearsf.append(Internal.getValue(sn))
                else: snearf.append(1.)
    
    #--------------------------------------------------------
    # Get Reference State and model from body pyTree
    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input cgns.')
    # model: Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    # check Euler non consistant avec Musker
    if model == 'Euler': 
        for z in Internal.getZones(tb):
            ibctype = Internal.getNodeFromName2(z, 'ibctype')
            if ibctype is not None:
                ibctype = Internal.getValue(ibctype)
                if ibctype == 'Musker' or ibctype == 'Log': 
                    raise ValueError("In tb: governing equations (Euler) not consistent with ibc type (%s)"%(ibctype))

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)
    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced

    #--------------------------------------------------------
    # Generates the full Cartesian mesh
    t = TIBM.generateIBMMesh(tb, vmin=vmin, snears=snears, dfar=dfar, dfarList=dfarList, DEPTH=2,
                             tbox=tbox, snearsf=snearsf, check=check, sizeMax=1000000,
                             expand=expand)
    test.printMem(">>> Build octree full [end]")
    
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
    test.printMem(">>> wall distance [start]")
    if dimPb == 2:
        z0 = Internal.getZones(t)
        bb = G.bbox(z0); dz = bb[5]-bb[2]
        tb2 = C.initVars(tb, 'CoordinateZ', dz*0.5)
        DTW._distance2Walls(t,tb2,type='ortho', signed=0, dim=dimPb, loc='centers')
    else:
        DTW._distance2Walls(t,tb,type='ortho', signed=0, dim=dimPb, loc='centers')
    test.printMem(">>> wall distance [end]")
    
    #----------------------------------------
    # Create IBM info
    #----------------------------------------
    t,tc = TIBM.prepareIBMData(t, tb, frontType=frontType, interpDataType=0)
    test.printMem(">>> ibm data [end]")

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
    if isinstance(t_out, str): Fast.save(t, t_out, split=format, NP=-NP, cartesian=True)
    return t, tc

#==================================================================================================

#================================================================================
# IBM prepare - para
#
# extrusion: make an extrusion from a 2D profile. ATTENTION, each zone of the profile must be joined in one single zone
# smoothing : smooth the front during the front 2 specific treatment in the cases of local refinements
# balancing ; balance the entire distribution after the octree generation, useful for symetries
# distrib : new distribution at the end of prepare1
#================================================================================
def prepare1(t_case, t_out, tc_out, snears=0.01, dfar=10., dfarList=[],
             tbox=None, snearsf=None,  
             vmin=21, check=False, NP=0, format='single',
             frontType=1, extrusion=False, smoothing=False, balancing=False, 
             distrib=True, expand=3):
    import Generator
    import Connector.connector as connector
    import Connector.Mpi as Xmpi
    import Post.PyTree as P
    import Converter.Distributed as Distributed
    import Connector.OversetData as XOD
    import KCore.test as test
    from mpi4py import MPI
    import numpy
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD

    # list of dfars
    if dfarList == []:
        zones = Internal.getZones(tb)
        dfarList = [dfar*1.]*len(zones)
        for c, z in enumerate(zones): 
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None: dfarList[c] = Internal.getValue(n)*1.

    # a mettre dans la classe ou en parametre de prepare1 ??? 
    to = None

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
                else: snearf.append(1.)
    
    symmetry = 0
    fileout = None
    if check: fileout = 'octree.cgns'

    DEPTH=2
    IBCType=1

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree defined in %s.'%FILE)
    # model : Euler, NSLaminar, NSTurbulent
    model = Internal.getValue(model)

    # check Euler non consistant avec Musker
    if model == 'Euler': 
        for z in Internal.getZones(tb):
            ibctype = Internal.getNodeFromName2(z, 'ibctype')
            if ibctype is not None:
                ibctype = Internal.getValue(ibctype)
                if ibctype == 'Musker' or ibctype == 'Log': 
                    raise ValueError("In tb: governing equations (Euler) not consistent with ibc type (%s)"%(ibctype))

    if dimPb == 2: C._initVars(tb, 'CoordinateZ', 0.) # forced

    # Octree identical on all procs
    test.printMem('>>> Octree unstruct [start]')

    o = TIBM.buildOctree(tb, snears=snears, snearFactor=1., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                         dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=rank,
                         expand=expand)

    if rank==0 and check: C.convertPyTree2File(o, fileout)

    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = TIBM.buildParentOctrees__(o, tb, snears=snears, snearFactor=4., dfar=dfar, dfarList=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                                        dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=rank)
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
    
    # Add xzones for ext
    test.printMem(">>> extended cart grids [start]")
    tbb = Cmpi.createBBoxTree(t)
    interDict = X.getIntersectingDomains(tbb)
    graph = Cmpi.computeGraph(tbb, type='bbox', intersectionsDict=interDict, reduction=False)
    del tbb
    Cmpi._addXZones(t, graph, variables=[], cartesian=True)
    test.printMem(">>> extended cart grids [after add XZones]")
    zones = Internal.getZones(t)
    coords = C.getFields(Internal.__GridCoordinates__, zones, api=2)
    coords = Generator.generator.extendCartGrids(coords, DEPTH+1, 1)
    C.setFields(coords, zones, 'nodes')
    Cmpi._rmXZones(t)
    coords = None; zones = None
    test.printMem(">>> extended cart grids (after rmXZones) [end]")
    
    if not extrusion:
        TIBM._addBCOverlaps(t, bbox=bb)
        TIBM._addExternalBCs(t, bbox=bb, dimPb=dimPb)

    dz = 0.01
    if dimPb == 2:
        if not extrusion:
            T._addkplane(t)
            T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)
        if extrusion:
            chord = 1.
            NSplit = 1
            NPas = 200
            span = 0.25*chord
            dimPb = 3
            # Extrude 2D case
            T._addkplane(tb,N=NPas+4)
            for node in Internal.getNodesFromName(tb,'EquationDimension'): Internal.setValue(node,3)
            T._contract(tb, (0.,0.,0.), (1,0,0), (0,1,0), span/NPas)
            zmax = C.getMaxValue(tb,'CoordinateZ')
            T._translate(tb,(0.,0.,-0.5*zmax))
            # Close new 3D case
            for b in Internal.getBases(tb):
                name = Internal.getName(b)
                b = C.convertArray2Tetra(b)
                b = G.close(b)
                b = P.exteriorFaces(b)
                b = T.splitConnexity(b)
                for line in Internal.getZones(b):
                    closure = G.tetraMesher(line, algo=1)
                    tb = Internal.append(tb, closure, name)
            if rank==0: C.convertPyTree2File(tb, '3Dcase.cgns')
            # create new 3D tree
            t = T.subzone(t, (1,1,1), (-1,-1,1))
            bbox = G.bbox(t); bbox = [round(i,1) for i in bbox]
            bbox = numpy.array(bbox)
            # Share the boundaries of the entire mesh for BCFarfield
            comm.Barrier()
            minbox = numpy.zeros(3)
            maxbox = numpy.zeros(3)
            comm.Allreduce([bbox[0:3], MPI.DOUBLE], [minbox, MPI.DOUBLE], MPI.MIN)
            comm.Allreduce([bbox[3:], MPI.DOUBLE], [maxbox, MPI.DOUBLE], MPI.MAX)
            comm.Barrier()
            bbox[0:3] = minbox
            bbox[3:]  = maxbox
            C._rmBCOfType(t, 'BCFarfield')
            C._rmBCOfType(t, 'BCOverlap')
            Internal._rmNodesByType(t,'FlowSolution_t')
            for z in Internal.getZones(t):
                xmin = C.getValue( z, 'CoordinateX', (1,1,1))
                xmax = C.getValue( z, 'CoordinateX', (0,1,1))
                ymin = C.getValue( z, 'CoordinateY', (1,1,1))
                ymax = C.getValue( z, 'CoordinateY', (1,0,1))
                if abs(round(xmin-bbox[0]))==0.: C._addBC2Zone(z, 'external', 'BCFarfield', 'imin')
                if abs(round(xmax-bbox[3]))==0.: C._addBC2Zone(z, 'external', 'BCFarfield', 'imax')
                if abs(round(ymin-bbox[1]))==0.: C._addBC2Zone(z, 'external', 'BCFarfield', 'jmin')
                if abs(round(ymax-bbox[4]))==0.: C._addBC2Zone(z, 'external', 'BCFarfield', 'jmax')
            C._fillEmptyBCWith(t,'overlap','BCOverlap')
            T._addkplane(t,N=NPas+4)
            for node in Internal.getNodesFromName(t,'EquationDimension'): Internal.setValue(node,3)
            T._contract(t, (0.,0.,0.), (1,0,0), (0,1,0), span/NPas)
            T._translate(t,(0.,0.,-0.5*zmax))
            C._addBC2Zone(t, 'period', 'BCautoperiod', 'kmin')
            C._addBC2Zone(t, 'period', 'BCautoperiod', 'kmax')
            if check: Cmpi.convertPyTree2File(t, '3Dmesh.cgns')

    # ReferenceState 
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)

    # Balancing
    if balancing:
        test.printMem(">>> balancing [start]")
        Cmpi.convertPyTree2File(t, t_out)
        # Need to wait for all the procs to write their parts before the new distribution
        comm.Barrier()
        ts = Cmpi.convertFile2SkeletonTree(t_out)
        D2._distribute(ts, Cmpi.size, algorithm='graph')
        t = Cmpi.readZones(ts, t_out, rank=rank)
        t = Cmpi.convert2PartialTree(t)
        zones = Internal.getZones(t)
        for z in zones: z[0] = z[0] + 'X%d'%rank
        del ts
        test.printMem(">>> balancing [end]")

    # Distance a la paroi    
    test.printMem(">>> Wall distance [start]")
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]
        tb2 = C.initVars(tb,'CoordinateZ',dz*0.5)
        DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dimPb, loc='centers')
    else:
        DTW._distance2Walls(t, tb, type='ortho', signed=0, dim=dimPb, loc='centers')

    X._applyBCOverlaps(t, depth=DEPTH, loc='centers', val=2, cellNName='cellN')
    
    # Blank des corps chimere
    # applyBCOverlap des maillages de corps
    # SetHoleInterpolated points

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
    C._initVars(t,'{centers:cellNIBC}={centers:cellN}')
    
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
            #TIBM._blankClosestTargetCells(t,cellNName='cellN',depth=DEPTHL)
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
    tc = C.node2Center(tp); del tp
    
    # setInterpData parallel pour le chimere
    tbbc = Cmpi.createBBoxTree(tc)
    interDict = X.getIntersectingDomains(tbbc)
    graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)
    Cmpi._addXZones(tc, graph, variables=['cellN'], cartesian=True)
    test.printMem(">>> Interpdata [after addXZones]")
    
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

    ############################################################
    # Specific treatment for front 2
    ############################################################
    if frontType == 2:
        test.printMem(">>> pushBackImageFront2 [start]")

        # bboxDict needed for optimised AddXZones (i.e. "layers" not None)
        # Return a dict with the zones of t as keys and their specific bboxes as key values
        bboxDict = Cmpi.createBboxDict(t)
        tbbc = Cmpi.createBBoxTree(tc)
        interDict = X.getIntersectingDomains(tbbc)
        graph = Cmpi.computeGraph(tbbc, type='bbox', intersectionsDict=interDict, reduction=False)

        # if subr, the tree subregions are kept during the exchange 
        # if layers not None, only communicate the desired number of layers
        Cmpi._addXZones(tc, graph, variables=['cellNIBC','cellNChim'], cartesian=True, interDict=interDict, bboxDict=bboxDict, layers=4, subr=False)
        Cmpi._addXZones(t, graph,variables=['centers:cellNIBC', 'centers:cellNChim'], cartesian=True, interDict=interDict, bboxDict=bboxDict, layers=4, subr=False)

        # Zones of tc are modified after addXZones, new tbbc, interDict and intersectionDict
        tbbcx = Cmpi.createBBoxTree(tc)
        interDict = X.getIntersectingDomains(tbbcx)
        intersectionsDict = X.getIntersectingDomains(tbbcx, method='AABB', taabb=tbbcx)
        
        # Reconstruction of cellNFront and cellN from cellNIBC (reduce the communications)
        # cellNFront_origin and cellNIBC_origin are initialised to store the Data of cellNFront and cellNIBC before the transfers 
        C._initVars(t,'{centers:cellNFront}=({centers:cellNIBC}==1.)')
        C._initVars(t,'{centers:cellN}={centers:cellNIBC}')
        C._initVars(t,'{centers:cellNFront_origin}={centers:cellNFront}') 
        C._initVars(t,'{centers:cellNIBC_origin}={centers:cellNIBC}')
        C._initVars(t,'{centers:cellN_interp}=maximum(0.,{centers:cellNChim})')

        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')
        C._cpVars(t,'centers:cellNIBC',tc,'cellNIBC')
        C._cpVars(t,'centers:cellN',tc,'cellN')
        C._cpVars(t,'centers:cellN_interp',tc,'cellN_interp')
        C._cpVars(t,'centers:cellNFront_origin',tc,'cellNFront_origin')
        C._cpVars(t,'centers:cellNIBC_origin',tc,'cellNIBC_origin')

        # Find each zone that require the specific treatment 
        C._initVars(t,'{centers:cellNFront2}=1.-({centers:cellNFront}<1.)*(abs({centers:cellNChim})>1.)')
        # i.e. if cellNFront_origin == 2 and cellNFront == 1 ou -3 => cellNFront2 = 1

        # Transfers the information at each grid connection
        for z in Internal.getZones(t):
            cellNFront = Internal.getNodeFromName2(z,'cellNFront2')
            if cellNFront != []:
                cellNFront = cellNFront[1]
                sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
                sizeOne =  int(numpy.sum(cellNFront))
                if sizeOne < sizeTot:
                    X._setHoleInterpolatedPoints(z, depth=1, dir=0, loc='centers',cellNName='cellNFront2',addGC=False)
                    res = X.getInterpolatedPoints(z,loc='centers', cellNName='cellNFront2') # indices,X,Y,Z
                    if res is not None:
                        indicesI = res[0]
                        XI = res[1]; YI = res[2]; ZI = res[3]
                        allInterpFields=[]
                        for zc in Internal.getZones(tc):
                            if zc[0] in intersectionsDict[z[0]]:
                                C._cpVars(zc,'cellN_interp',zc,'cellN')
                                fields = X.transferFields(zc, XI, YI, ZI, hook=None, variables=['cellNFront_origin','cellNIBC_origin'], interpDataType=0,nature=1)
                                allInterpFields.append(fields)
                        if allInterpFields!=[]:
                            C._filterPartialFields(z, allInterpFields, indicesI, loc='centers', startFrom=0, filterName='donorVol',verbose=False)

        Cmpi._rmXZones(tc)
        Cmpi._rmXZones(t)

        # Update the cellNFront, cellNIBC and cellNIBCDnr fields
        for z in Internal.getZones(t):
            cellNFront = Internal.getNodeFromName2(z,'cellNFront2')
            if cellNFront != []:
                cellNFront = cellNFront[1]
                sizeTot = cellNFront.shape[0]*cellNFront.shape[1]*cellNFront.shape[2]
                sizeOne =  int(numpy.sum(cellNFront))
                if sizeOne < sizeTot:
                    C._initVars(z,'{centers:cellNFront}={centers:cellNFront}*({centers:cellNFront_origin}>0.5)') # Modification du Front uniquement lorsque celui-ci est repousse
                    # i.e. if cellNFront_origin == 0 and cellNFront == 1 => cellNfront = 0

                    C._initVars(z,'{centers:cellNIBC}={centers:cellNIBC}*(1.-({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)) \
                        + 2.*({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)')
                    # i.e. if cellNChim == 1 and cellNIBC_origin == 2 => cellNIBC = 2

                    C._initVars(z,'{centers:cellNIBCDnr}={centers:cellNIBCDnr}*(1.-({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)) \
                        + 2.*({centers:cellNChim}==1.)*({centers:cellNIBC_origin}>1.5)*({centers:cellNIBC_origin}<2.5)')

        C._cpVars(t,'centers:cellNIBC',tc,'cellNIBC')
        C._cpVars(t,'centers:cellNIBC',t,'centers:cellN')
        C._cpVars(t,'centers:cellN',tc,'cellN')

        C._rmVars(t,['centers:cellNFront2'])
        C._rmVars(t,['centers:cellNFront_origin'])
        C._rmVars(t,['centers:cellNIBC_origin'])
        C._rmVars(t,['centers:cellN_interp'])

        # Smooth the front in case of a local refinement - only work in 2D
        if smoothing and dimPb==2: 
            TIBM._smoothImageFront(t, tc)

        C._cpVars(t,'centers:cellNFront',tc,'cellNFront')

        Xmpi._setInterpTransfers(t,tc,variables=['cellNFront'], cellNVariable='cellNFront', compact=0)
        test.printMem(">>> pushBackImageFront2 [end]")
    ############################################################

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

    datas = {}; graph = {}

    C._initVars(t,'{centers:cellN}=minimum({centers:cellNChim}*{centers:cellNIBCDnr},2.)')
    varsRM = ['centers:cellNChim','centers:cellNIBCDnr']
    if model == 'Euler': varsRM += ['centers:TurbulentDistance']
    C._rmVars(t, varsRM)
    test.printMem(">>> Saving [start]")

    # Sauvegarde des infos IBM
    if check:
        test.printMem(">>> Saving IBM infos [start]")
        tibm = TIBM.extractIBMInfo(tc)

        # Avoid that two procs write the same information
        for z in Internal.getZones(tibm):
           if int(z[0][-1]) != rank:
              Internal._rmNodesByName(tibm, z[0])

        Cmpi.convertPyTree2File(tibm, 'IBMInfo.cgns')
        test.printMem(">>> Saving IBM infos [end]")
        del tibm

    # distribution par defaut (sur NP)
    tbbc = Cmpi.createBBoxTree(tc)
    
    # Perform the finale distribution
    if distrib:
        if NP == 0: NP = Cmpi.size
        stats = D2._distribute(tbbc, NP, algorithm='graph', useCom='ID')
        D2._copyDistribution(tc, tbbc)
        D2._copyDistribution(t, tbbc)

    del tbbc

    # Save tc
    if isinstance(tc_out, str): Cmpi.convertPyTree2File(tc, tc_out, ignoreProcNodes=True)
    I._initConst(t, loc='centers')
    if model != "Euler": C._initVars(t, 'centers:ViscosityEddy', 0.)
    # Save t
    if isinstance(t_out, str):
        import Compressor.PyTree as Compressor
        Compressor._compressCartesian(t)
        Cmpi.convertPyTree2File(t, t_out, ignoreProcNodes=True)

    if Cmpi.size > 1: Cmpi.barrier()
    return t, tc

#=============================================================================
# Post
# IN: t_case: fichier ou arbre du cas
# IN: t_in: fichier ou arbre de resultat
# IN: tc_in: fichier ou arbre de connectivite
# OUT: t_out ou None: fichier pour sortie du champ aux noeuds
# OUT: wall_out ou None: fichier pour sortie du champ sur la paroi  
#==============================================================================
def post(t_case, t_in, tc_in, t_out, wall_out):
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
    if model is None: raise ValueError('GoverningEquations is missing in input cgns.')
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

    for z in Internal.getNodesFromType2(t, "Zone_t"):
        zc = Internal.getNodeFromName(tc, z[0])
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
        C._initVars(t, 'CoordinateZ', 0.) # forced

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
    t = C.center2Node(t, vars)
    Internal._rmNodesByName(t, 'FlowSolution#Centers')
    if isinstance(t_out, str): C.convertPyTree2File(t, t_out)

    return t, zw

#=============================================================================
# Post Efforts
# IN: t_case: fichier ou arbre du cas
# IN: t_in: fichier ou arbre de resultat
# IN: tc_in: fichier ou arbre de connectivite
# OUT: wall_out ou None: fichier pour sortie des efforts sur la paroi aux centres
# IN: alpha: angle pour les efforts
# IN: beta: angle pour les efforts
#==============================================================================
def loads(t_case, tc_in, wall_out, alpha=0., beta=0., Sref=None):
    import Post.PyTree as P
    import Converter.Filter as Filter
    import math
    import numpy as numpy

    if isinstance(tc_in, str): tc = C.convertFile2PyTree(tc_in)
    else: tc = tc_in
    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    if Sref is None:
        C._initVars(tb, '__ONE__',1.)
        Sref = P.integ(tb, '__ONE__')[0]; print(Sref)
        C._rmVars(tb, ['__ONE__', 'centers:vol'])

    # Dans le cas du CRM version adimensionnee :
    # Sref = 2*191.8445/(7.00532**2)

    #==============================
    # Reference state
    #==============================
    [RoInf, RouInf, RovInf, RowInf, RoeInf, PInf, TInf, cvInf, MInf,
          ReInf, Cs, Gamma, RokInf, RoomegaInf, RonutildeInf,
          Mus, Cs, Ts, Pr] = C.getState(tb)

    alpha = math.radians(alpha)
    beta = math.radians(beta)

    dimPb = Internal.getValue(Internal.getNodeFromName(tb, 'EquationDimension'))

    q = 0.5*RoInf*(MInf*math.sqrt(Gamma*PInf/RoInf))**2

    #====================================
    # Extraction des grandeurs a la paroi
    #====================================
    zw = TIBM.extractIBMWallFields(tc, tb=tb)

    if dimPb == 2: zw = T.addkplane(zw)

    zw = C.convertArray2Tetra(zw)
    zw = T.reorderAll(zw, 1)
    C._initVars(zw, 'Cp=-({Pressure}-%f)/%f'%(PInf,q))

    #===========================
    # Calcul efforts de pression
    #===========================
    res = P.integNorm(zw, 'Cp')[0]
    res = [i/Sref for i in res]
    cd = res[0]*math.cos(alpha)*math.cos(beta) + res[2]*math.sin(alpha)*math.cos(beta)
    cl = res[2]*math.cos(alpha)*math.cos(beta) - res[0]*math.sin(alpha)*math.cos(beta)
    print("Pressure loads", cd, cl)

    #======================================
    # Calcul frottement et efforts visqueux
    #======================================
    if C.isNamePresent(zw, 'utau') != -1:
        C._initVars(zw, '{tau_wall}={Density}*{utau}**2')
    else:
        C._initVars(zw, '{tau_wall}=0.')

    G._getNormalMap(zw)
    zw = C.node2Center(zw, ['Cp', 'tau_wall', 'Pressure','VelocityX','VelocityY','VelocityZ'])
    if C.isNamePresent(zw, 'utau') != -1:
        zw = C.node2Center(zw, ['utau','yplus'])
    C._rmVars(zw, 'FlowSolution')
    C._normalize(zw, ['centers:sx','centers:sy','centers:sz'])

    # calcul du vecteur tangent
    C._initVars(zw, '{centers:tx}={centers:VelocityX}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sx}')
    C._initVars(zw, '{centers:ty}={centers:VelocityY}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sy}')
    C._initVars(zw, '{centers:tz}={centers:VelocityZ}-({centers:VelocityX}*{centers:sx}+{centers:VelocityY}*{centers:sy}+{centers:VelocityZ}*{centers:sz})*{centers:sz}')
    C._normalize(zw, ['centers:tx','centers:ty','centers:tz'])

    C._initVars(zw, '{centers:tauxx}=2*{centers:tau_wall}*{centers:tx}*{centers:sx}')
    C._initVars(zw, '{centers:tauyy}=2*{centers:tau_wall}*{centers:ty}*{centers:sy}')
    C._initVars(zw, '{centers:tauzz}=2*{centers:tau_wall}*{centers:tz}*{centers:sz}')
    C._initVars(zw, '{centers:tauxy}={centers:tau_wall}*({centers:tx}*{centers:sy}+{centers:ty}*{centers:sx})')
    C._initVars(zw, '{centers:tauxz}={centers:tau_wall}*({centers:tx}*{centers:sz}+{centers:tz}*{centers:sx})')
    C._initVars(zw, '{centers:tauyz}={centers:tau_wall}*({centers:ty}*{centers:sz}+{centers:tz}*{centers:sy})')

    # calcul frottement
    C._initVars(zw, '{centers:Fricx}={centers:tauxx}*{centers:sx}+{centers:tauxy}*{centers:sy}+{centers:tauxz}*{centers:sz}')
    C._initVars(zw, '{centers:Fricy}={centers:tauxy}*{centers:sx}+{centers:tauyy}*{centers:sy}+{centers:tauyz}*{centers:sz}')
    C._initVars(zw, '{centers:Fricz}={centers:tauxz}*{centers:sx}+{centers:tauyz}*{centers:sy}+{centers:tauzz}*{centers:sz}')

    # calcul effort complet
    C._initVars(zw, '{centers:Fx}={centers:Fricx}-({centers:Pressure}-%f)*{centers:sx}'%PInf)
    C._initVars(zw, '{centers:Fy}={centers:Fricy}-({centers:Pressure}-%f)*{centers:sy}'%PInf)
    C._initVars(zw, '{centers:Fz}={centers:Fricz}-({centers:Pressure}-%f)*{centers:sz}'%PInf)

    # calcul coefficient de frottement
    C._initVars(zw, '{centers:Cf}=(sqrt({centers:Fricx}**2+{centers:Fricy}**2+{centers:Fricz}**2))/%f'%q)

    zw = G.getVolumeMap(zw)
    effortX = P.integ(zw, 'centers:Fricx')[0]
    zw = G.getVolumeMap(zw)
    effortY = P.integ(zw, 'centers:Fricy')[0]
    zw = G.getVolumeMap(zw)
    effortZ = P.integ(zw, 'centers:Fricz')[0]

    cd = (effortX*math.cos(alpha)*math.cos(beta) + effortZ*math.sin(alpha)*math.cos(beta))/q
    cl = (effortZ*math.cos(alpha)*math.cos(beta) - effortX*math.sin(alpha)*math.cos(beta))/q
    print("Skin friction loads", cd, cl)

    vars = ['centers:sx','centers:sy','centers:sz','centers:tx','centers:ty','centers:tz','centers:tauxx','centers:tauyy','centers:tauzz','centers:tauxy','centers:tauxz',
'centers:tauyz','centers:Fricx','centers:Fricy','centers:Fricz','centers:Fx','centers:Fy','centers:Fz']
    zw = C.rmVars(zw, vars)
    if dimPb == 2: # reextrait en 2D
        zw = P.isoSurfMC(zw, "CoordinateZ", 0.)
        nodes = Internal.getNodesFromName(zw, 'CoordinateX')
        xmin = numpy.min(nodes[0][1])
        xmax = numpy.max(nodes[0][1])
        C._initVars(zw, 'xc=({CoordinateX}-%f)/(%f-%f)'%(xmin, xmax, xmin))

    if isinstance(wall_out, str): C.convertPyTree2File(zw, wall_out)
    return zw
    
#====================================================================================
# Redistrib on NP processors
#====================================================================================
def _distribute(t_in, tc_in, NP):
    if isinstance(tc_in, str):
        tcs = Cmpi.convertFile2SkeletonTree(tc_in, maxDepth=3)
    else: tcs = tc_in
    stats = D2._distribute(tcs, NP, algorithm='graph', useCom='ID')
    print(stats)
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
        Filter.writeNodesFromPaths(tc_in, paths, ns, maxDepth=0, mode=1)

    if isinstance(t_in, str):
        ts = Cmpi.convertFile2SkeletonTree(t_in, maxDepth=3)
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
        Filter.writeNodesFromPaths(t_in, paths, ns, maxDepth=0, mode=1)

    # Affichage du nombre de points par proc - equilibrage ou pas
    NptsTot = 0
    for i in range(NP):
        NPTS = 0
        for z in Internal.getZones(ts):
            if Cmpi.getProc(z) == i: NPTS += C.getNPts(z)
        NptsTot += NPTS
        print('Rank %d has %d points'%(i,NPTS))
    print('All points: %d million points'%(NptsTot/1.e6))
    return None

#====================================================================================
# Prend les snears dans t, les multiplie par factor
def snearFactor(t, factor=1.):
    tp = Internal.copyRef(t)
    _snearFactor(t, value)
    return tp

def _snearFactor(t, factor=1.):
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromName2(z, 'snear')
        for n in nodes:
            Internal._setValue(n, factor*Internal.getValue(n))
    return None

# Set IBC type in zones
def setIBCType(t, value):
    tp = Internal.copyRef(t)
    _setIBCType(t, value)
    return tp
    
def _setIBCType(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'ibctype', 'DataArray_t', value)
    return None

# Set snear in zones
def setSnear(t, value):
    tp = Internal.copyRef(t)
    _setSnear(t, value)
    return tp

def _setSnear(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'snear', 'DataArray_t', value)
    return None

# Set dfar in zones 
def setDfar(t, value):
    tp = Internal.copyRef(t)
    _setDfar(t, value)
    return tp

def _setDfar(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'dfar', 'DataArray_t', value)
    return None

def _modifIBCD(tc):
    # modif BC outpress
    nodes = Internal.getNodesFromName(tc, 'IBCD_4_*')
    for n in nodes:
        pres = Internal.getNodeFromName1(n, 'Pressure')
        Internal.setValue(pres, 74500.*numpy.ones(numpy.shape(pres[1])))
    # modif BC injection
    alpha = 0.*pi/180.
    nodes = I.getNodesFromName(tc, 'dirx')
    for n in nodes:
        d = numpy.ones(numpy.shape(n[1]))
        Internal.setValue(n, d)
    nodes = Internal.getNodesFromName(tc, 'StagnationEnthalpy')
    for n in nodes:
        Internal.setValue(n, stagEnt*numpy.ones(numpy.shape(n[1])))
    nodes = Internal.getNodesFromName(tc, 'StagnationPressure')
    for n in nodes:
        Internal.setValue(n, PiInj*numpy.ones(numpy.shape(n[1])))
    return None

#====================================================================================
class IBM(Common):
    """Preparation et caculs avec le module FastS."""
    def __init__(self, format=None, numb=None, numz=None):
        Common.__init__(self, format, numb, numz)
        self.__version__ = "0.0"
        self.authors = ["ash@onera.fr"]
        self.cartesian = True
        
    # Prepare 
    def prepare(self, t_case, t_out, tc_out, snears=0.01, dfar=10., dfarList=[], 
                tbox=None, snearsf=None,
                vmin=21, check=False, frontType=1, NP=None, expand=3):
        if NP is None: NP = Cmpi.size
        if NP == 0: print('Preparing for a sequential computation.')
        else: print('Preparing for a computation on %d processors.'%NP)
        ret = prepare(t_case, t_out, tc_out, snears=snears, dfar=dfar, dfarList=dfarList,
                      tbox=tbox, snearsf=snearsf,
                      vmin=vmin, check=check, NP=NP, format=self.data['format'], 
                      frontType=frontType, expand=expand)
        return ret

    # post-processing: extrait la solution aux noeuds + le champs sur les surfaces
    def post(self, t_case, t_in, tc_in, t_out, wall_out):
        return post(t_case, t_in, tc_in, t_out, wall_out)

    # post-processing: extrait les efforts sur les surfaces
    def loads(self, t_case, t_in, tc_in, wall_out, alpha=0., beta=0., Sref=None):
        return loads(t_case, t_in, tc_in, wall_out, alpha=0., beta=0., Sref=None)
