"""Toolbox for CODA."""
import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.ToolboxIBM as TIBM
import Transform.PyTree as T
import Converter
import Converter.Distributed as CD
import Generator.PyTree as G
import Post.PyTree as P
import Converter.Mpi as Cmpi
import numpy
import Dist2Walls.PyTree as DTW
import KCore.test as test
import Connector.PyTree as X
import Transform
import Geom.PyTree as D
import fnmatch

TOL = 1.e-6
rank = Cmpi.rank; size = Cmpi.size

# suppose 1 connect quad par zone
def extractBCOfType(t,bndType):
    res = []
    for z in Internal.getZones(t):
        zp = Internal.copyRef(z)
        Internal._rmNodesFromType(zp,"FlowSolution_t")
        connects = Internal.getNodesFromType1(z, "Elements_t")
        cnQUAD = None
        for cn in connects:
            val = Internal.getValue(cn)[0]
            if val == 7 and cnQUAD is None:  # QUAD
                cnQUAD = cn
                ER = Internal.getNodeFromName(cnQUAD,'ElementRange')
                ER[1][1] = ER[1][1]-ER[1][0]+1
                ER[1][0] = 1
                Internal.getValue(cn)[1] = 0
            else:
                Internal._rmNodesFromName(zp,cn[0])

        # physical BCs
        for bc in Internal.getNodesFromType2(z, 'BC_t'):
            if Internal.getValue(bc) == bndType:
                PL = Internal.getNodeFromName(bc,'PointList')
                GridLoc = Internal.getNodeFromType(bc,'GridLocation_t')
                GridLoc = Internal.getValue(GridLoc)
                if PL is not None and GridLoc=='FaceCenter':
                    PL = Internal.getValue(PL)[0].tolist() # list of elements
                    array = C.getFields(Internal.__GridCoordinates__, zp, api=1)[0]
                    array = Transform.transform.subzoneElements(array, PL)
                    res.append(C.convertArrays2ZoneNode(bc[0],[array]))
        # connect 1to1
        if bndType=='BCMatch':
            for gc in Internal.getNodesFromType2(z,'GridConnectivity1to1_t'):
                PL = Internal.getNodeFromName(gc,'PointList')
                GridLoc = Internal.getNodeFromType(gc,'GridLocation_t')
                GridLoc = Internal.getValue(GridLoc)
                if PL is not None and GridLoc == 'FaceCenter':
                    PL = Internal.getValue(PL)[0].tolist() # list of elements
                    array = C.getFields(Internal.__GridCoordinates__, zp, api=1)[0]
                    array = Transform.transform.subzoneElements(array, PL)
                    res.append(C.convertArrays2ZoneNode(gc[0],[array]))
    return res

def prepare(t_case, t, tskel, check=False):
    IBCType = 1 # Musker
    externalBCType = 'BCFarfield' # external boundaries
    prefzone ='CART_P' #zone names are eventually prefixed by this name

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    tBB = Cmpi.createBBoxTree(t)
    procDict = Cmpi.getProcDict(tBB)
    bbo = G.bbox(tBB)
    rank = Cmpi.rank
    comm = Cmpi.COMM_WORLD

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree.')
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
    if check:
        Cmpi.convertPyTree2File(tBB,'tBB.cgns')
    # only inter-proc connectivity are stored
    # PointList and PointListDonor are actually FaceList, start at 0, ordered as Cartesian !!!
    prefzone ='CART_P' #zone names are eventually prefixed by this name
    for z in Internal.getZones(t):
        gcnodes = Internal.getNodesFromType(z, 'GridConnectivity1to1_t')
        for gcnode in gcnodes:
            zdnrname = Internal.getValue(gcnode)
            dnrproc = procDict[zdnrname]
            if dnrproc == Cmpi.rank:
                Internal._rmNode(z, gcnode)
    graphM={}
    _create1To1Connectivity(t, tskel=tskel, dim=dimPb, convertOnly=True)
    dictOfAbuttingSurfaces={}
    for z in Internal.getZones(t):
        for gc in Internal.getNodesFromType2(z,'GridConnectivity_t'):
            indicesL = Internal.getNodeFromName2(gc,'PointList')
            indicesL = Internal.getValue(indicesL)
            zf=T.subzoneUnstruct__(z, indicesL, 'faces')
            Internal._rmNodesFromType(zf,"ZoneGridConnectivity_t")
            Internal._rmNodesFromType(zf,"ZoneBC_t")
            zdname = Internal.getValue(gc)
            PLD = Internal.getNodeFromName(gc,'PointListDonor')
            norankopp = procDict[zdname]
            if norankopp in dictOfAbuttingSurfaces:
                dictOfAbuttingSurfaces[norankopp] += [z[0], zf, zdname, PLD]
            else:
                dictOfAbuttingSurfaces[norankopp] = [[z[0], zf, zdname, PLD]]
                CD.updateGraph__(graphM, rank, norankopp, prefzone+'%d'%norankopp)

    # Ensures that the graph is the same for all the processors
    if Cmpi.KCOMM is not None: g = Cmpi.KCOMM.allgather(graphM)
    else: g = [graphM]
    graphM = {}
    for i in g:
        for k in i:
            if not k in graphM: graphM[k] = {}
            for j in i[k]:
                if not j in graphM[k]: graphM[k][j] = []
                graphM[k][j] += i[k][j]
                graphM[k][j] = list(set(graphM[k][j]))

    # Addkplane in 2D case
    dz = 0.01
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t,"Zone_t")
        dimZ = Internal.getZoneDim(z0)
        if dimZ[3]==1:
            T._addkplane(t)
            T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)
        elif dimZ[3]== 2:
            CoordZ = Internal.getNodeFromName(z0,'CoordinateZ')
            CoordZ = Internal.getValue(CoordZ)
            npts = dimZ[1]*dimZ[2]*dimZ[3]
            dz = CoordZ[-1,-1,1]-CoordZ[0,0,0]

    # ReferenceState
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)

    # Distance to IBCs (all IBCs)
    test.printMem(">>> Wall distance [start]")
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]
        tb2 = C.initVars(tb, 'CoordinateZ', dz*0.5)
        DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dimPb, loc='nodes')
    else:
        DTW._distance2Walls(t, tb, type='ortho', signed=0, dim=dimPb, loc='nodes')
    test.printMem(">>> Wall distance [end]")

    # blanking
    test.printMem(">>> Blanking [start]")
    C._initVars(t,'cellN',1.)
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, 'Zone_t')
        dims = Internal.getZoneDim(z0)
        npts = dims[1]*dims[2]*dims[3]
        zmin = C.getValue(z0,'CoordinateZ',0)
        zmax = C.getValue(z0,'CoordinateZ',npts-1)
        dz = zmax-zmin
        # Creation of the 2D body for IBM preprocessing
        T._addkplane(tb)
        T._contract(tb, (0,0,0), (1,0,0), (0,1,0), dz)
    t = TIBM.blankByIBCBodies(t, tb, 'nodes', dimPb)
    test.printMem(">>> Blanking [end]")
    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
        dims = Internal.getZoneDim(i)
        npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))
    C._initVars(t,'{TurbulentDistance}=-1.*({cellN}<1.)*{TurbulentDistance}+({cellN}>0.)*{TurbulentDistance}')
    print('Minimum distance: %f.'%C.getMinValue(t,'TurbulentDistance'))
    t = C.node2Center(t,['TurbulentDistance'])
    t = P.computeGrad(t, 'TurbulentDistance')
    graddvars=["centers:gradxTurbulentDistance",'centers:gradyTurbulentDistance','centers:gradzTurbulentDistance']
    t = C.center2Node(t, graddvars)
    #print("After computeGrad : Perform transfers of gradient correctly ????")
    #Internal._rmNodesFromName(t,Internal.__FlowSolutionCenters__)
    C._rmVars(t, graddvars)
    #
    # Extract front faces
    if IBCType == -1:
        X._setHoleInterpolatedPoints(t,loc='nodes',depth=-1)
        C._initVars(t,'{cellN}=minimum(1.,{cellN})')

    # Removal of fully blanked zones (must be done after gradient of distance for a correct gradient estimation near the obstacles)
    for z in Internal.getZones(t):
        if C.getMaxValue(z,'cellN') < 0.5:
            (parent,noz) = Internal.getParentOfNode(t, z)
            del parent[2][noz]

    print("Extract front faces of IBM target points...")
    # IBM target points
    frontType=1
    varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']

    front1 = []
    he = 0.
    frontDict = {}
    for z in Internal.getZones(t):
        f = P.frontFaces(z,'cellN')
        if Internal.getZoneDim(f)[1]>0:
            he = max(he, C.getMaxValue(f,'TurbulentDistance'))
            frontDict[z[0]]=f
        else:
            frontDict[z[0]]=[]
    if check:
        Cmpi.convertPyTree2File(front1,'targetFaces.cgns')

    print(" Compute IBM Wall points...")

    he = he*1.8 # distmax = sqrt(3)*dx => he min = distmax + dx + tol
    #varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
    allip_pts=[]; allwall_pts=[]; allimage_pts=[]

    zbcs=[]; bctypes=[]
    for bc in Internal.getNodesFromType(t,'BC_t'):
        bctype = Internal.getValue(bc)
        if bctype not in bctypes: bctypes.append(bctype)

    for bctype in bctypes:
        zbc = C.extractBCOfType(t,bctype)
        zbc = C.convertArray2Hexa(zbc)
        zbc = T.join(zbc)
        zbcs.append(zbc)

    for z in Internal.getZones(t):
        parentz,noz = Internal.getParentOfNode(t,z)
        ip_ptsZ = frontDict[z[0]]
        C._rmVars(z,varsn)
        if ip_ptsZ != []:
            ip_ptsZC = C.node2Center(ip_ptsZ)
            C._rmVars(ip_ptsZ,varsn)
            ip_ptsZC = C.convertArray2Node(ip_ptsZC)
            C._rmVars(z,['TurbulentDistance'])
            z = P.selectCells(z,"{cellN}==1.",strict=1)
            C._rmVars(z,['cellN'])
            #
            # IBM Wall points
            #
            zname = ip_ptsZ[0]
            ip_pts = C.getAllFields(ip_ptsZC,loc='nodes', api=1)[0]
            ip_pts = Converter.convertArray2Node(ip_pts)
            wallpts = T.projectAllDirs(ip_ptsZC, tb, varsn, oriented=0)
            C._normalize(wallpts,varsn)
            wallpts = C.convertArray2Node(wallpts)
            for var in varsn:
                C._initVars(wallpts,'{%s}=%g*{%s}'%(var,he,var))
            imagepts = T.deform(wallpts, vector=varsn)
            wallpts = C.getFields(Internal.__GridCoordinates__,wallpts, api=1)[0]
            imagepts = C.getFields(Internal.__GridCoordinates__,imagepts, api=1)[0]
            ip_pts = Converter.extractVars(ip_pts,['CoordinateX','CoordinateY','CoordinateZ'])
            wallpts = Converter.extractVars(wallpts,['CoordinateX','CoordinateY','CoordinateZ'])
            imagepts = Converter.extractVars(imagepts,['CoordinateX','CoordinateY','CoordinateZ'])
            if check:
                Converter.convertArrays2File(ip_pts,"targetPts_%s.plt"%zname)
                Converter.convertArrays2File(wallpts,"wallPts_%s.plt"%zname)
                Converter.convertArrays2File(imagepts,"imagePts_%s.plt"%zname)

            allip_pts.append(ip_pts)
            allimage_pts.append(imagepts)
            allwall_pts.append(wallpts)

        else:
            if C.getMinValue(z,'cellN')==0.:
                z = P.selectCells(z,"{cellN}==1.",strict=1)
                C._rmVars(z,['TurbulentDistance','cellN'])
            else:
                C._rmVars(z,['TurbulentDistance','cellN'])
                z = C.convertArray2Hexa(z)

        parentz[2][noz] = z

    front1=[]
    for zname in frontDict:
        if frontDict[zname] != []:
            front1.append(frontDict[zname])
    if front1 != []: front1 = T.join(front1)
    if check:
        C.convertPyTree2File(front1,'targetFaces_%d.cgns'%rank)

    zones = Internal.getZones(t)
    z = zones[0]
    for noz in range(1,len(zones)):
        z = T.join([z,zones[noz]])
    z[0] = prefzone+'%d'%rank

    # recoverBCs
    indicesF=[]; f = P.exteriorFaces(z, indices=indicesF)
    indicesF=indicesF[0]
    hook = C.createHook(f, 'elementCenters')
    for nobc, zbc in enumerate(zbcs):
        ids = C.identifyElements(hook, zbc, tol=TOL)
        ids = ids[ids[:] > -1]
        ids = ids.tolist()
        ids = [ids[i]-1 for i in range(len(ids))]
        zf = T.subzone(f,ids, type='elements')
        _addBC2ZoneLoc(z, zf[0], bctypes[nobc], zf)

    # recoverIBC familySpecified:IBMWall
    if front1 != []: # can be empty on a proc
        front1[0]='IBMWall'
        _addBC2ZoneLoc(z, 'IBMWall', 'FamilySpecified:IBMWall',front1)

        for bc in Internal.getNodesFromType(z,'BC_t'):
            famName = Internal.getNodeFromName(bc,'FamilyName')
            if famName is not None:
                if Internal.getValue(famName)=='IBMWall':
                    ibcdataset=Internal.createNode('BCDataSet','BCDataSet_t',parent=bc,value='Null')
                    dnrPts = Internal.createNode("DonorPointCoordinates",'BCData_t',parent=ibcdataset)
                    wallPts = Internal.createNode("WallPointCoordinates",'BCData_t',parent=ibcdataset)
                    allwall_pts = Transform.join(allwall_pts)
                    allimage_pts = Transform.join(allimage_pts)
                    coordsPI = Converter.extractVars(allimage_pts, ['CoordinateX','CoordinateY','CoordinateZ'])
                    dnrPts[2].append(['CoordinateX',coordsPI[1][0,:], [], 'DataArray_t'])
                    dnrPts[2].append(['CoordinateY',coordsPI[1][1,:], [], 'DataArray_t'])
                    dnrPts[2].append(['CoordinateZ',coordsPI[1][2,:], [], 'DataArray_t'])
                    coordsPW = Converter.extractVars(allwall_pts, ['CoordinateX','CoordinateY','CoordinateZ'])
                    wallPts[2].append(['CoordinateX',coordsPW[1][0,:], [], 'DataArray_t'])
                    wallPts[2].append(['CoordinateY',coordsPW[1][1,:], [], 'DataArray_t'])
                    wallPts[2].append(['CoordinateZ',coordsPW[1][2,:], [], 'DataArray_t'])

    # modify PL for abutting 1to1
    datas={}
    for norankopp in dictOfAbuttingSurfaces:
        infos = dictOfAbuttingSurfaces[norankopp]
        for info in infos:
            zrname = info[0]
            zf = info[1]
            zdname = info[2]
            PLD = info[3]

            ids = C.identifyElements(hook, zf, tol=TOL)
            ids = ids[ids[:] > -1]
            ids = ids.tolist()
            if ids != []:
                ids = [ids[i]-1 for i in range(len(ids))]
                zf = T.subzone(f,ids, type='elements')
                _addBC2ZoneLoc(z, zf[0], 'Abutting1to1', zf, zdnrName=zdname)
    C.freeHook(hook)
    mergeQuadConn(z)
    # convertBCAbutting to GC Abutting
    # Add abutting 1to1 on one side only-(smallest proc number chosen)
    for zbc in Internal.getNodesFromType(z,'BC_t'):
        if Internal.getValue(zbc)=='Abutting1to1':
            zdnrname = Internal.getNodeFromType1(zbc,'UserDefinedData_t')
            zdnrname = Internal.getValue(zdnrname)
            PL = Internal.getNodeFromName(zbc, 'PointList')
            norankopp = procDict[zdnrname]

            if rank < norankopp:
                ZGC = Internal.getNodeFromType1(z,'ZoneGridConnectivity_t')
                if ZGC is None:
                    ZGC = Internal.createNode('ZoneGridConnectivity', 'ZoneGridConnectivity_t', parent=z)

                matchname = 'match_%d_%d'%(rank,norankopp)
                gcnode = Internal.createNode(matchname,'GridConnectivity1to1_t', parent=ZGC, value='CART_P%d'%norankopp)
                gcnode[2].append(Internal.createNode('GridLocation','GridLocation_t',value='FaceCenter'))
                gcnode[2].append(PL)

            if rank > norankopp:
                matchname = 'match_%d_%d'%(norankopp, rank)
                if norankopp not in datas:
                    datas[norankopp]=[[matchname,PL[1]]]
                else:
                    datas[norankopp].append([matchname,PL[1]])

    #transfer PLD
    rcvDataM = Cmpi.sendRecv(datas, graphM)
    if rcvDataM != {}:
        for origProc in rcvDataM:
            dataorig = rcvDataM[origProc]
            for info in dataorig:
                matchname=info[0]
                PLD = info[1]
                gcnode = Internal.getNodeFromName(z,matchname)
                if gcnode is not None:
                    PLDN = Internal.createNode('PointListDonor','IndexArray_t',value=PLD,parent=gcnode)

    # remove all bc_t Abutting1to1
    for zbc in Internal.getNodesFromType(z,'BC_t'):
        if Internal.getValue(zbc)=='Abutting1to1': Internal._rmNode(z,zbc)

    t[2][1][2]=[z]

    Internal._renameNode(t,'FlowSolution#Centers','FlisWallDistance')
    return t

#==============================================
# IBM prepro for CODA with Octree mesh
#==============================================
def prepareOctree(t_case, t_out, vmin=5, dfarList=[], dfar=10., snears=0.01, NP=0,
                  tbox=None, snearsf=None, expand=3, check=False, fileout='octree.cgns'):
    recoverBC=False
    isHanging=False
    dfarDir=0
    IBCType=1
    externalBCType = 'BCFarfield'

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

    to = None
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

    # reference state
    refstate = C.getState(tb)
    # dimension du pb
    dimPb = Internal.getNodeFromName(tb, 'EquationDimension')
    dimPb = Internal.getValue(dimPb)

    model = Internal.getNodeFromName(tb, 'GoverningEquations')
    if model is None: raise ValueError('GoverningEquations is missing in input tree.')
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

    o = TIBM.buildOctree(tb, snears=snears, snearFactor=1., dfars=dfarList,
                         to=to, tbox=tbox, snearsf=snearsf,
                         dimPb=dimPb, vmin=vmin,
                         expand=expand, dfarDir=dfarDir)
    if rank==0 and check: C.convertPyTree2File(o, fileout)
    bbo = G.bbox(o)

    # build parent octree 3 levels higher
    # returns a list of 4 octants of the parent octree in 2D and 8 in 3D
    parento = TIBM.buildParentOctrees__(o, tb, snears=snears, snearFactor=4., dfars=dfarList, to=to, tbox=tbox, snearsf=snearsf,
                                        dimPb=dimPb, vmin=vmin)
    test.printMem(">>> Octree unstruct [end]")

    # Split octree
    test.printMem(">>> Octree unstruct split [start]")
    bb = G.bbox(o)
    NPI = Cmpi.size
    if NPI == 1: p = Internal.copyRef(o) # keep reference
    else: p = T.splitNParts(o, N=NPI, recoverBC=recoverBC)[rank]
    del o
    test.printMem(">>> Octree unstruct split [end]")
    # C.convertPyTree2File(p,"octree_%d.cgns"%rank)

    # fill vmin + merge in parallel
    test.printMem(">>> Octree struct [start]")
    res = TIBM.octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000)
    del p
    if parento is not None:
        for po in parento: del po

    t = C.newPyTree(['CARTESIAN', res])
    for z in Internal.getZones(t):
        zname = Internal.getName(z)+'_X%d'%rank
        Internal.setName(z,zname)
    Cmpi._setProc(t, rank)

    # only inter-proc connectivity are stored
    # PointList and PointListDonor are actually FaceList, start at 0, ordered as Cartesian !!!
    prefzone ='CART_P' #zone names are eventually prefixed by this name
    graphM={}
    _create1To1Connectivity(t, dim=dimPb)

    dictOfAbuttingSurfaces={}
    for z in Internal.getZones(t):
        for gc in Internal.getNodesFromType2(z,'GridConnectivity_t'):
            indicesL = Internal.getNodeFromName2(gc,'PointList')
            indicesL = Internal.getValue(indicesL)
            zf=T.subzoneUnstruct__(z, indicesL, 'faces')
            Internal._rmNodesFromType(zf,"ZoneGridConnectivity_t")
            Internal._rmNodesFromType(zf,"ZoneBC_t")
            zdname = Internal.getValue(gc)
            zdname = zdname.split("_X")
            norankopp = int(zdname[-1])
            if norankopp in dictOfAbuttingSurfaces:
                dictOfAbuttingSurfaces[norankopp]+=[zf]
            else:
                dictOfAbuttingSurfaces[norankopp]=[zf]
                CD.updateGraph__(graphM, rank, norankopp, prefzone+'%d'%norankopp)

    # Ensures that the graph is the same for all the processors

    if Cmpi.KCOMM is not None: g = Cmpi.KCOMM.allgather(graphM)
    else: g = [graphM]
    graphM = {}
    for i in g:
        for k in i:
            if not k in graphM: graphM[k] = {}
            for j in i[k]:
                if not j in graphM[k]: graphM[k][j] = []
                graphM[k][j] += i[k][j]
                graphM[k][j] = list(set(graphM[k][j]))
    #C.convertPyTree2File(t,"t_%d.cgns"%(Cmpi.rank))

    C._addState(t, 'EquationDimension', dimPb)
    test.printMem(">>> Octree struct [end]")

    dz = 0.01
    if dimPb == 2:
        T._addkplane(t)
        T._contract(t, (0,0,0), (1,0,0), (0,1,0), dz)

    # ReferenceState
    C._addState(t, state=refstate)
    C._addState(t, 'GoverningEquations', model)
    C._addState(t, 'EquationDimension', dimPb)

    # Distance to IBCs (all IBCs)
    test.printMem(">>> Wall distance [start]")
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]
        tb2 = C.initVars(tb, 'CoordinateZ', dz*0.5)
        DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dimPb, loc='nodes')
    else:
        DTW._distance2Walls(t, tb, type='ortho', signed=0, dim=dimPb, loc='nodes')
    test.printMem(">>> Wall distance [end]")

    # blanking
    test.printMem(">>> Blanking [start]")
    C._initVars(t,'cellN',1.)
    if dimPb == 2:
        z0 = Internal.getNodeFromType2(t, 'Zone_t')
        dims = Internal.getZoneDim(z0)
        npts = dims[1]*dims[2]*dims[3]
        zmin = C.getValue(z0,'CoordinateZ',0)
        zmax = C.getValue(z0,'CoordinateZ',npts-1)
        dz = zmax-zmin
        # Creation of the 2D body for IBM preprocessing
        T._addkplane(tb)
        T._contract(tb, (0,0,0), (1,0,0), (0,1,0), dz)
    t = TIBM.blankByIBCBodies(t, tb, 'nodes', dimPb)
    test.printMem(">>> Blanking [end]")
    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    for i in Internal.getZones(t):
        dims = Internal.getZoneDim(i)
        npts += dims[1]*dims[2]*dims[3]
    print('Final number of points=%5.4f millions.'%(npts/1000000.))
    C._initVars(t,'{TurbulentDistance}=-1.*({cellN}<1.)*{TurbulentDistance}+({cellN}>0.)*{TurbulentDistance}')
    print('Minimum distance: %f.'%C.getMinValue(t,'TurbulentDistance'))
    t = C.node2Center(t,['TurbulentDistance'])
    t = P.computeGrad(t, 'TurbulentDistance')
    graddvars=["centers:gradxTurbulentDistance",'centers:gradyTurbulentDistance','centers:gradzTurbulentDistance']
    t = C.center2Node(t,graddvars)
    #print("After computeGrad : Perform transfers of gradient correctly ????")
    #Internal._rmNodesFromName(t,Internal.__FlowSolutionCenters__)
    C._rmVars(t,graddvars)
    #
    # Extract front faces
    if IBCType == -1:
        X._setHoleInterpolatedPoints(t,loc='nodes',depth=-1)
        C._initVars(t,'{cellN}=minimum(1.,{cellN})')

    # Removal of fully blanked zones (must be done after gradient of distance for a correct gradient estimation near the obstacles)
    for z in Internal.getZones(t):
        if C.getMaxValue(z,'cellN') < 0.5:
            (parent,noz) = Internal.getParentOfNode(t, z)
            del parent[2][noz]

    print("Extract front faces of IBM target points...")
    # IBM target points
    frontType=1
    varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']

    front1 =[]
    he = 0.
    frontDict={}
    for z in Internal.getZones(t):
        f = P.frontFaces(z,'cellN')
        if Internal.getZoneDim(f)[1]>0:
            he = max(he, C.getMaxValue(f,'TurbulentDistance'))
            frontDict[z[0]]=f
        else:
            frontDict[z[0]]=[]
    if check:
        Cmpi.convertPyTree2File(front1,'targetFaces.cgns')

    print(" Compute IBM Wall points...")

    loc='FaceCenter'
    he = he*1.8 # distmax = sqrt(3)*dx => he min = distmax + dx + tol
    #varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
    allip_pts=[]; allwall_pts=[]; allimage_pts=[]
    _addExternalBCs(t, bbo, externalBCType, dimPb)

    extBCs=C.extractBCOfType(t,externalBCType)
    extBCs=C.convertArray2Hexa(extBCs)

    # join par niveau de resolution
    dhMin = 1e10
    dictOfExtBCs={}
    for extBC in Internal.getZones(extBCs):
        dxloc = abs(C.getValue(extBC,'CoordinateX',1)-C.getValue(extBC,'CoordinateX',0))
        dyloc = abs(C.getValue(extBC,'CoordinateY',1)-C.getValue(extBC,'CoordinateY',0))
        if dxloc > 0. and dhMin > dxloc: dhMin = dxloc
        if dyloc > 0. and dhMin > dyloc: dhMin = dyloc

    for extBC in Internal.getZones(extBCs):
        dxloc = abs(C.getValue(extBC,'CoordinateX',1)-C.getValue(extBC,'CoordinateX',0))
        dyloc = abs(C.getValue(extBC,'CoordinateY',1)-C.getValue(extBC,'CoordinateY',0))
        dh = max(dxloc,dyloc)
        lvl = int(dh/dhMin)
        if lvl not in dictOfExtBCs:
            dictOfExtBCs[lvl]=[extBC]
        else:
            dictOfExtBCs[lvl].append(extBC)
    extBCs=[]
    for lvl in dictOfExtBCs:
        extBCs.append(T.join(dictOfExtBCs[lvl]))

    for z in Internal.getZones(t):
        parentz,noz = Internal.getParentOfNode(t,z)
        ip_ptsZ = frontDict[z[0]]
        C._rmVars(z,varsn)
        if ip_ptsZ != []:
            ip_ptsZC = C.node2Center(ip_ptsZ)
            C._rmVars(ip_ptsZ,varsn)
            ip_ptsZC = C.convertArray2Node(ip_ptsZC)
            C._rmVars(z,['TurbulentDistance'])
            z = P.selectCells(z,"{cellN}==1.",strict=1)
            C._rmVars(z,['cellN'])
            #
            # IBM Wall points
            #
            zname = ip_ptsZ[0]
            ip_pts = C.getAllFields(ip_ptsZC,loc='nodes', api=1)[0]
            ip_pts = Converter.convertArray2Node(ip_pts)
            wallpts = T.projectAllDirs(ip_ptsZC, tb, varsn, oriented=0)
            C._normalize(wallpts,varsn)
            wallpts = C.convertArray2Node(wallpts)
            for var in varsn:
                C._initVars(wallpts,'{%s}=%g*{%s}'%(var,he,var))
            imagepts = T.deform(wallpts, vector=varsn)
            wallpts = C.getFields(Internal.__GridCoordinates__,wallpts, api=1)[0]
            imagepts = C.getFields(Internal.__GridCoordinates__,imagepts, api=1)[0]
            ip_pts = Converter.extractVars(ip_pts,['CoordinateX','CoordinateY','CoordinateZ'])
            wallpts = Converter.extractVars(wallpts,['CoordinateX','CoordinateY','CoordinateZ'])
            imagepts = Converter.extractVars(imagepts,['CoordinateX','CoordinateY','CoordinateZ'])
            # Converter.convertArrays2File(ip_pts,"targetPts_%s.plt"%zname)
            # Converter.convertArrays2File(wallpts,"wallPts_%s.plt"%zname)
            # Converter.convertArrays2File(imagepts,"imagePts_%s.plt"%zname)
            allip_pts.append(ip_pts)
            allimage_pts.append(imagepts)
            allwall_pts.append(wallpts)
        else:

            if C.getMinValue(z,'cellN')==0.:
                z = P.selectCells(z,"{cellN}==1.",strict=1)
                C._rmVars(z,['TurbulentDistance','cellN'])
            else:
                C._rmVars(z,['TurbulentDistance','cellN'])
                z = C.convertArray2Hexa(z)

        parentz[2][noz] = z

    front1=[]
    for zname in frontDict:
        if frontDict[zname] != []:
            front1.append(frontDict[zname])
    if front1 != []: front1 = T.join(front1)
    if check:
        C.convertPyTree2File(front1,'targetFaces_%d.cgns'%rank)

    zones = Internal.getZones(t)
    z = zones[0]
    for noz in range(1,len(zones)):
        z = T.join([z,zones[noz]])
    z[0] = prefzone+'%d'%rank

    # externalBCs:
    for extBC in extBCs:
        extBC[0]='ExtBC_Quad'
        C._addBC2Zone(z,'externalBC',externalBCType,subzone=extBC)
    # recoverIBC familySpecified:IBMWall
    if front1 != []: # can be empty on a proc
        front1[0]='IBMWall'
        C._addBC2Zone(z, 'IBMWall', 'FamilySpecified:IBMWall',subzone=front1)
        for bc in Internal.getNodesFromType(z,'BC_t'):
            famName = Internal.getNodeFromName(bc,'FamilyName')
            if famName is not None:
                if Internal.getValue(famName)=='IBMWall':
                    ibcdataset=Internal.createNode('BCDataSet','BCDataSet_t',parent=bc,value='Null')
                    dnrPts = Internal.createNode("DonorPointCoordinates",'BCData_t',parent=ibcdataset)
                    wallPts = Internal.createNode("WallPointCoordinates",'BCData_t',parent=ibcdataset)
                    #allip_pts=Transform.join(allip_pts)
                    allwall_pts = Transform.join(allwall_pts)
                    allimage_pts = Transform.join(allimage_pts)
                    coordsPI = Converter.extractVars(allimage_pts, ['CoordinateX','CoordinateY','CoordinateZ'])
                    dnrPts[2].append(['CoordinateX',coordsPI[1][0,:], [], 'DataArray_t'])
                    dnrPts[2].append(['CoordinateY',coordsPI[1][1,:], [], 'DataArray_t'])
                    dnrPts[2].append(['CoordinateZ',coordsPI[1][2,:], [], 'DataArray_t'])
                    coordsPW = Converter.extractVars(allwall_pts, ['CoordinateX','CoordinateY','CoordinateZ'])
                    wallPts[2].append(['CoordinateX',coordsPW[1][0,:], [], 'DataArray_t'])
                    wallPts[2].append(['CoordinateY',coordsPW[1][1,:], [], 'DataArray_t'])
                    wallPts[2].append(['CoordinateZ',coordsPW[1][2,:], [], 'DataArray_t'])

        #_addIBDataZSR(z,[allip_pts],[allwall_pts],[allimage_pts], prefix='IBCD_')

    # Add abutting 1to1 on one side only-(smallest proc number chosen)
    datas={}
    for norankopp in dictOfAbuttingSurfaces:
        zf = T.join(dictOfAbuttingSurfaces[norankopp])
        hookVertex = C.createHook(zf,function='nodes')
        idNodes = C.identifyNodes(hookVertex,z)
        C.freeHook(hookVertex)
        PL = idNodes[(idNodes>-1)!=0]
        # create GC
        if PL != [] and rank < norankopp:
            ZGC = Internal.getNodeFromType2(z,'ZoneGridConnectivity_t')
            if ZGC is None:
                ZGC = Internal.createNode('ZoneGridConnectivity', 'ZoneGridConnectivity_t', parent=z)

            matchname = 'match_%d_%d'%(rank,norankopp)

            gcnode = Internal.createNode(matchname,'GridConnectivity_t', parent=ZGC, value='CART_P%d'%norankopp)
            gcnode[2].append(Internal.createNode('GridLocation','GridLocation_t',value='Vertex'))
            lsize = len(PL.shape)
            if lsize == 1: PL = numpy.reshape(PL,(1, PL.shape[0]), order='F')
            gcnode[2].append(Internal.createNode('PointList','IndexArray_t',value=PL))
        if PL != [] and rank > norankopp:
            matchname = 'match_%d_%d'%(norankopp, rank)
            datas[norankopp]=[matchname,PL]

    #transfer PLD
    rcvDataM = Cmpi.sendRecv(datas, graphM)
    if rcvDataM != {}:
        for origProc in rcvDataM:
            matchname=rcvDataM[origProc][0]
            PLD = rcvDataM[origProc][1]
            lsize = len(PLD.shape)
            if lsize == 1: PLD = numpy.reshape(PLD,(1, PLD.shape[0]), order='F')
            gcnode = Internal.getNodeFromName(z,matchname)
            gcnode[2].append(Internal.createNode('PointListDonor','IndexArray_t',value=PLD))
    mergeQuadConn0(z)
    EltRange2FaceCenter(z)
    # add2PyTree
    #Cmpi._setProc(z,rank)
    t[2][1][2]=[z]

    Internal._rmNodesFromName(t,".Solver#Param")
    Internal._renameNode(t,'FlowSolution#Centers','FlisWallDistance')
    # identify elements per processor
    if t_out is not None: Cmpi.convertPyTree2File(t, t_out)

    return t
# IN: zoneS : zone structured with BCs
# IN: zoneU : zone HEXA without BCs - to be recovered as a single QUAD conn + BCs with FaceCenter
# zoneU is modified in place
def _recoverBCsStruct2QUAD(zoneS, zoneU):
    # externalBCs:
    indicesF=[]; f = P.exteriorFaces(zoneU, indices=indicesF)
    indicesF=indicesF[0]
    hook = C.createHook(f, 'elementCenters')
    zbcs=[]; bctypes=[]
    for bc in Internal.getNodesFromType(zoneS,'BC_t'):
        bctype = Internal.getValue(bc)
        if bctype not in bctypes: bctypes.append(bctype)
    for bctype in bctypes:
        zbc = C.extractBCOfType(zoneS,bctype)
        zbc = C.convertArray2Hexa(zbc)
        zbc = T.join(zbc)
        zbcs.append(zbc)

    for zbc in zbcs:
        ids = C.identifyElements(hook, zbc, tol=1e-6)
        ids = ids[ids[:] > -1]
        ids = ids.tolist()
        ids = [ids[i]-1 for i in range(len(ids))]
        zf = T.subzone(f,ids, type='elements')
        _addBC2ZoneLoc(zoneU, zbc[0], bctype, zbc)

    C.freeHook(hook)

    return None

def _addBC2ZoneLoc(z, bndName, bndType, zbc, loc='FaceCenter', zdnrName=None):
    s = bndType.split(':')
    bndType1 = s[0]
    if len(s) > 1: bndType2 = s[1]
    else: bndType2 = ''

    # Analyse zone zbc
    dims = Internal.getZoneDim(zbc)
    neb = dims[2] # nbre d'elts de zbc

    eltType, nf = Internal.eltName2EltNo(dims[3]) # type d'elements de zbc
    # On cherche l'element max dans les connectivites de z
    maxElt = 0
    connects = Internal.getNodesFromType(z, 'Elements_t')
    for cn in connects:
        r = Internal.getNodeFromName1(cn, 'ElementRange')
        m = r[1][1]
        maxElt = max(maxElt, m)

    # on cree un nouveau noeud connectivite dans z1 (avec le nom de la zone z2)
    nebb = neb
    node = Internal.createUniqueChild(z, zbc[0], 'Elements_t', value=[eltType,nebb])
    Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t',
                               value=[maxElt+1,maxElt+neb])
    oldc = Internal.getNodeFromName2(zbc, 'ElementConnectivity')[1]
    newc = numpy.copy(oldc)
    hook = C.createHook(z, 'nodes')
    ids = C.identifyNodes(hook, zbc)
    newc[:] = ids[oldc[:]-1]
    faceList = [i for i in range(neb)]

    Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)

    zoneBC = Internal.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')
    if len(s)==1:
        info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType)
    else: # familyspecified
        info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType1)
        Internal.createUniqueChild(info, 'FamilyName', 'FamilyName_t',
                                   value=bndType2)

    Internal.createUniqueChild(info, 'GridLocation', 'GridLocation_t',
                               value='FaceCenter')
    if isinstance(faceList, numpy.ndarray): r = faceList
    else: r = numpy.array(faceList, dtype=numpy.int32)
    r = r.reshape((1,r.size), order='F')
    info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])
    if bndType == 'Abutting1to1':
        info[2].append(["UserDefined", zdnrName, [], 'UserDefinedData_t'])
    return None

def stretchCart(z, dir=1, N=1, width=0.5, ind=(1,1,1), h=0.):
    i1 = ind[0]; j1 = ind[1]; k1 = ind[2]
    if dir == 1:
        zt = T.subzone(z, (1,j1,k1), (N,j1,k1))
    elif dir == 2:
        zt = T.subzone(z, (i1,1,k1), (i1,N,k1))
        zt = T.reorder(zt, (2,1,3))
    else:
        zt = T.subzone(z, (i1,j1,1), (i1,j1,N))
        zt = T.reorder(zt, (3,2,1))

    l = D.getLength(zt)
    zt = D.getCurvilinearAbscissa(zt)
    distrib = C.cpVars(zt, 's', zt, 'CoordinateX')
    distrib = C.initVars(distrib, 'CoordinateY', 0.)
    distrib = C.initVars(distrib, 'CoordinateZ', 0.)
    distrib = C.rmVars(distrib, 's')

    Nr = int(width*N)
    #print h, l, Nr, i1, j1, k1

    if dir == 1:
        val = C.getValue(zt, 's', i1-1)
        indl = (i1,j1,k1); indp1 = (i1+1,j1,k1); indm1 = (i1-1,j1,k1)
    elif dir == 2:
        val = C.getValue(zt, 's', j1-1)
        indl = (i1,j1,k1); indp1 = (i1,j1+1,k1); indm1 = (i1,j1-1,k1)
    else:
        val = C.getValue(zt, 's', k1-1)
        indl = (i1,j1,k1); indp1 = (i1,j1,k1+1); indm1 = (i1,j1,k1-1)

    if dir == 1:
        if i1 == 1:
            distrib = G.enforcePlusX(distrib, h/l, Nr, Nr)
        elif i1 == N:
            distrib = G.enforceMoinsX(distrib, h/l, Nr, Nr)

    elif dir == 2:
        if j1 == 1:
            distrib = G.enforcePlusX(distrib, h/l, Nr, Nr)
        elif j1 == N:
            distrib = G.enforceMoinsX(distrib, h/l, Nr, Nr)

    elif dir == 3:
        if k1 == 1:
            distrib = G.enforcePlusX(distrib, h/l, Nr, Nr)
        elif k1 == N:
            distrib = G.enforceMoinsX(distrib, h/l, Nr, Nr)

    z = G.map(z, distrib, dir)
    return z

def cartRxLoc(X0, H, N, Nb, depth=0, addCellN=False, addBCMatch=False, rank=None, size=None):
    """Create a set of regular cartesian grids."""
    out = []
    for k in range(Nb[2]):
        for j in range(Nb[1]):
            for i in range(Nb[0]):
                if rank is None or size is None or rank == (i+j*Nb[0]+k*Nb[0]*Nb[1])%size:
                    Xp = [X0[0]+H[0]*(N[0]-1)*i,X0[1]+H[1]*(N[1]-1)*j,X0[2]+H[2]*(N[2]-1)*k]
                    Np = [N[0],N[1],N[2]]
                    if i > 0: Xp[0] -= depth*H[0]; Np[0] += depth
                    if i < Nb[0]-1: Xp[0] += depth*H[0]; Np[0] += depth
                    if j > 0: Xp[1] -= depth*H[1]; Np[1] += depth
                    if j < Nb[1]-1: Xp[1] += depth*H[1]; Np[1] += depth
                    if k > 0: Xp[2] -= depth*H[2]; Np[2] += depth
                    if k < Nb[2]-1: Xp[2] += depth*H[2]; Np[2] += depth
                    z = G.cart(Xp, H, Np); z[0] = 'cart%d.%d.%d'%(i,j,k)
                    if rank is not None:
                        import Converter.Mpi as Cmpi
                        Cmpi._setProc(z, rank)
                    if addCellN:
                        C._initVars(z, 'centers:cellN', 1)
                        cellN = Internal.getNodeFromName2(z, 'cellN')[1]
                        if i > 0: cellN[0:depth,:,:] = 2
                        if i < Nb[0]-1: cellN[Np[0]-depth-1:Np[0]-1,:,:] = 2
                        if j > 0: cellN[:,0:depth,:] = 2
                        if j < Nb[1]-1: cellN[:,Np[1]-depth-1:Np[1]-1,:] = 2
                        if k > 0: cellN[:,:,0:depth] = 2
                        if k < Nb[2]-1: cellN[:,:,Np[2]-depth-1:Np[2]-1] = 2
                    if addBCMatch and depth == 0:
                        if i > 0:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'imin', z, 'imax', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i-1,j,k))
                        if i < Nb[0]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'imax', z, 'imin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i+1,j,k))
                        if j > 0:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'jmin', z, 'jmax', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i,j-1,k))
                        if j < Nb[1]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'jmax', z, 'jmin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i,j+1,k))
                        if k > 0:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'kmin', z, 'kmax', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i,j,k-1))
                        if k < Nb[2]-1:
                            C._addBC2Zone(z, 'match', 'BCMatch', 'kmax', z, 'kmin', [1,2,3])
                            bcs = Internal.getNodesFromType2(z, 'GridConnectivity1to1_t')
                            Internal._setValue(bcs[-1], 'cart%d.%d.%d'%(i,j,k+1))
                    out.append(z)
    return out

# suppose que les connectivites quad sont contigues
def mergeQuadConn(z):
    NBVERTQ = 4
    BCs = Internal.getNodesFromType(z, 'BC_t')

    rminglob = -1; rmaxglob = -1
    rangeMinT={}; rangeMaxT={}
    allElts_t = Internal.getNodesFromType(z, 'Elements_t')
    for elts_t in allElts_t:
        typeEt = Internal.getValue(elts_t)[0]
        if typeEt == 7:
            name = Internal.getName(elts_t)
            eltRange = Internal.getNodeFromName1(elts_t,"ElementRange")
            rangeMinT[name] = Internal.getValue(eltRange)[0]
            rangeMaxT[name] = Internal.getValue(eltRange)[1]

            if rminglob == -1: rminglob = rangeMinT[name]
            else: rminglob = min(rminglob, rangeMinT[name])
            if rmaxglob == -1: rmaxglob = rangeMaxT[name]
            else: rmaxglob = max(rmaxglob, rangeMaxT[name])

    # first quad conn
    for name in rangeMinT:
        if rangeMinT[name]==rminglob:
            rmaxnext = rangeMaxT[name]+1
            EltsT = Internal.getNodeFromName(allElts_t,name)
            rangeMinT[name]=-1; rangeMaxT[name]=-1
            bc2 = Internal.getNodeFromName(BCs, name)
            if bc2 is not None:
                PL = Internal.getNodeFromName(bc2,'PointList')[1]
                nfaces = PL.shape[1]
                for indf in range(nfaces): PL[0,indf] += rminglob

            break

    # searching for next elts
    found = True
    while found:
        found = False
        for name in rangeMinT:
            if rangeMinT[name]==rmaxnext:
                found = True
                name2 = name
                EltsT2 = Internal.getNodeFromName(allElts_t,name2)
                # merge
                rmaxnext = rangeMaxT[name2]+1

                ec1 = Internal.getNodeFromName(EltsT, 'ElementConnectivity')
                ec2 = Internal.getNodeFromName(EltsT2, 'ElementConnectivity')
                eshift = ec1[1].shape[0]//NBVERTQ
                bc2 = Internal.getNodeFromName(BCs, name2)
                if bc2 is not None:
                    PL = Internal.getNodeFromName(bc2,'PointList')[1]
                    nfaces = PL.shape[1]
                    for indf in range(nfaces): PL[0,indf] += eshift+rminglob
                    ec1[1] = numpy.concatenate((ec1[1],ec2[1]))
                    ER = Internal.getNodeFromType(EltsT, 'IndexRange_t')
                    ERVal = Internal.getValue(ER)
                    ERVal[1] = rangeMaxT[name2]
                    Internal.setValue(ER,ERVal)
                    rangeMinT[name2]=-1; rangeMaxT[name2]=-1
                    Internal._rmNode(z,EltsT2)

    EltsT[0]='Quad'
    return None

def mergeQuadConn0(z):
    rangeMinT={}; rangeMaxT={}
    rmaxall = -1; rminall = 100000000

    for elts_t in Internal.getNodesFromType(z, "Elements_t"):
        typeEt = Internal.getValue(elts_t)[0]
        if typeEt == 7:
            name = Internal.getName(elts_t)
            eltRange = Internal.getNodeFromName1(elts_t, "ElementRange")
            rangeMinT[name] = Internal.getValue(eltRange)[0]
            rangeMaxT[name] = Internal.getValue(eltRange)[1]
            rmaxall = max(rangeMaxT[name]+1,rmaxall)
            rminall = min(rangeMinT[name],rminall)
        elif typeEt == 17:
            elts_t[0]='Hexas'
    # init
    rmin = rminall; rmax= -1
    newElts_t=[]; etype = 7
    while rmin < rmaxall:
        #start
        for name in rangeMinT:
            if rangeMinT[name]==rmin:
                rmax = rangeMaxT[name]
                EltsT = Internal.getNodeFromName1(z,name)
                name1 = name
                break
        # looking for following elts_t
        found = 1
        while found == 1:
            found = 0
            for name in rangeMinT:
                if rangeMinT[name]==rmax+1:
                    found = 1
                    name2 = name
                    EltsT2 = Internal.getNodeFromName1(z,name2)
            if found == 0:
                rmax2 = rmax
                #ER = Internal.getNodeFromName(EltsT,'ElementRange')
                #ER[1][1]=rmax
                newElts_t.append(EltsT)
            else:
                rmin2 = rangeMinT[name1]
                rmax2 = rangeMaxT[name2]
                ER = Internal.getNodeFromType(EltsT, 'IndexRange_t')
                ERVal = Internal.getValue(ER)
                ERVal[0] = rmin; ERVal[1] = rmax2
                Internal.setValue(ER,ERVal)
                ec1 = Internal.getNodeFromName(EltsT, 'ElementConnectivity')
                ec2 = Internal.getNodeFromName(EltsT2, 'ElementConnectivity')
                ec1[1] = numpy.concatenate((ec1[1],ec2[1]))
                rangeMinT[name2]=-1; rangeMaxT[name2]=-1
                rangeMinT[name1]=-1; rangeMaxT[name1]=-1

                rmax = rmax2
                if rmax2+1==rmaxall:
                    EltsT[0]='Quads'# a priori, unique name...

                    newElts_t.append(EltsT)
                    found = 0
            rmin = rmax2+1

    for elts_t in Internal.getNodesFromType(z,"Elements_t"):
        typeEt = Internal.getValue(elts_t)[0]
        if typeEt == 7:
            Internal._rmNode(z,elts_t)
    z[2] += newElts_t
    return None

# convert ZoneBC ElementRange to a PointList/FaceCenter according to QUAD connectivity
def EltRange2FaceCenter(z):
    connects = Internal.getNodesFromType(z,'Elements_t')
    connectQUAD=None

    for conn in connects:
        eltType = Internal.getValue(conn)[0]
        if eltType == 7:
            connectQUAD = conn
            break
    ERQUAD = Internal.getValue(Internal.getNodeFromName(connectQUAD,'ElementRange'))
    rminQUAD = ERQUAD[0]; rmaxQUAD = ERQUAD[1]
    for zbc in Internal.getNodesFromType(z,"BC_t"):
        ER = Internal.getNodeFromName(zbc, 'ElementRange')
        ERVal = Internal.getValue(ER)[0]
        rmin = ERVal[0]; rmax = ERVal[1]
        PL = numpy.zeros((1,rmax-rmin+1), dtype=numpy.int32, order='F')
        for i in range(PL.shape[1]):
            PL[0,i] = rmin+i-rminQUAD+1
        zbc[2].append(['PointList', PL, [], 'IndexRange_t'])
        gl = numpy.array([c for c in 'FaceCenter'],'c')
        zbc[2].append(['GridLocation', gl, [], 'GridLocation_t'])
        Internal._rmNodesFromName(zbc,"ElementRange")
    return None

def get1To1Connect(a, b, indicesFacesOrig, bopp, indicesFacesOrigOpp):

    b[0]='extFaces'; bopp[0] = 'extFacesOpp'
    hookExtFaces=C.createHook(bopp,function='elementCenters')

    # identify face centers
    idFaces = C.identifyElements(hookExtFaces,b)

    eltType = Internal.getZoneDim(b)[3]
    if eltType == 'BAR': nfaces = 4
    else: nfaces = 6
    ELTG = []; ELTD = []
    for noface in range(len(idFaces)):
        if idFaces[noface] != -1:
            nofaceopp = idFaces[noface]-1 # idFaces starts at 1
            etg = int(indicesFacesOrig[noface]-1)//nfaces
            etd = int(indicesFacesOrigOpp[nofaceopp]-1)//nfaces
            ELTG.append(etg)
            ELTD.append(etd)
    return [ELTG, ELTD]

def buildQuad4QuadInfo(a, graph=None):
    # quad4quad on same processor
    extFaces = P.exteriorFaces(a)
    Internal._rmNodesFromType(extFaces,'FlowSolution_t')
    Internal._rmNodesFromType(extFaces,'ZoneBC_t')
    Internal._rmNodesFromName(extFaces,'Quads')

    extFaces = T.splitConnexity(extFaces)
    res = []
    for noext in range(len(extFaces)):
        for noext2 in range(len(extFaces)):
            if noext != noext2:
                res+=buildQuad4QuadInfoLocal(a, extFaces[noext], extFaces[noext2])

    # quad4quad global
    if Cmpi.size > 1:# Send info to opposite procs
        rank = Cmpi.rank
        datas={}
        for opprank in graph[rank]:
            if opprank not in datas:
                datas[opprank]=[extFaces]
            else:
                datas[opprank].append(extFaces)

        destDatas=Cmpi.sendRecv(datas,graph)

        for i in destDatas:
            for extFacesOpp in destDatas[i]:
                for noext in range(len(extFaces)):
                    for noext2 in range(len(extFacesOpp)):
                        res+=buildQuad4QuadInfoLocal(a, extFacesOpp[noext2], extFaces[noext])

    erange = numpy.ones(2, dtype=numpy.int32)
    maxRange=-1
    for eltst in Internal.getNodesFromType(a,"Elements_t"):
        ER = Internal.getNodeFromName(eltst,'ElementRange')
        if ER is not None:
            ER = Internal.getValue(ER)
            maxRange = max(maxRange, ER[1]+1)
    erange[0] = maxRange
    erange[1] = maxRange-1 + len(res)//9

    if res == []:
        res = None
        erange[1] = erange[0]-1
    Internal.newElements(name='Quad4Quad', etype=1,
                         econnectivity=res,
                         erange=erange, eboundary=0, parent=a)
    Q4Q = Internal.getNodeFromName(a,"Quad4Quad")

    Internal.createNode('ElementTypeName', 'Descriptor_t', value='Quad4Quad', children=None, parent=Q4Q)
    return None

def buildQuad4QuadInfoLocal(a, b, bopp):
    b[0]='extFaces'; bopp[0] = 'extFacesOpp'
    shiftElt=4; nfaces = 6

    # identify face centers of coarse (local) zone matching with a vertex of bopp (fine)
    hookExtFacesOpp = C.createHook(bopp,function='nodes')
    indicesVEFOpp = C.identifyElements(hookExtFacesOpp,b) #>0: vertex 9 found (= bopp is fine and b is coarse)
    # identify indices of vertices in original zone a with vertices of exterior faces of bopp
    hookA = C.createHook(a,function='nodes')
    indicesVertexA = C.identifyNodes(hookA, bopp)

    cnExtFaceOpp = Internal.getNodeFromType(bopp,'Elements_t')
    cnExtFaceOpp = Internal.getNodeFromName(cnExtFaceOpp,'ElementConnectivity')
    cnExtFaceOpp = Internal.getValue(cnExtFaceOpp)
    sizeCNExtFaceOpp = cnExtFaceOpp.shape[0]

    Q4QArray=[]
    for noEltEF in range(len(indicesVEFOpp)):
        indVertexEF = indicesVEFOpp[noEltEF] # starts at 1
        if indVertexEF !=-1:
            noptr = 0
            GVIndices=[-1]*9 #ghost vertex indices
            # loop on ext faces elts
            GVIndices[8] = indVertexEF
            while noptr < sizeCNExtFaceOpp:
                indV1 = cnExtFaceOpp[noptr]
                indV2 = cnExtFaceOpp[noptr+1]
                indV3 = cnExtFaceOpp[noptr+2]
                indV4 = cnExtFaceOpp[noptr+3]

                if indVertexEF == indV1: #P9P6P3P7
                    GVIndices[5] = indV2
                    GVIndices[2] = indV3
                    GVIndices[6] = indV4
                elif indVertexEF == indV2:#P8P9P7P4
                    GVIndices[7] = indV1
                    GVIndices[6] = indV3
                    GVIndices[3] = indV4

                elif indVertexEF == indV3:#P1P5P9P6
                    GVIndices[0] = indV1
                    GVIndices[4] = indV2
                    GVIndices[5] = indV4

                elif indVertexEF == indV4:#P5P2P7P9
                    GVIndices[4] = indV1
                    GVIndices[1] = indV2
                    GVIndices[6] = indV3

                noptr+=shiftElt

            for indv in GVIndices:
                Q4QArray.append(indicesVertexA[indv-1])
    return Q4QArray

def _addIBCDataSet(bc,correctedPts, wallPts, imagePts):
    coordsPC = Converter.extractVars(correctedPts,['CoordinateX','CoordinateY','CoordinateZ'])[0]
    coordsPW = Converter.extractVars(wallPts, ['CoordinateX','CoordinateY','CoordinateZ'])[0]
    coordsPI = Converter.extractVars(imagePts, ['CoordinateX','CoordinateY','CoordinateZ'])[0]
    bcdataset = Internal.newBCDataSet(name='BCDataSet', gridLocation='FaceCenter', parent=bc)
    targetPointsFieldNode = Internal.newBCData('TargetPointCoordinates',parent=bcdataset)
    targetPointsFieldNode[2].append(['CoordinateX',coordsPC[1][0,:], [], 'DataArray_t'])
    targetPointsFieldNode[2].append(['CoordinateY',coordsPC[1][1,:], [], 'DataArray_t'])
    targetPointsFieldNode[2].append(['CoordinateZ',coordsPC[1][2,:], [], 'DataArray_t'])
    imagePointsFieldNode = Internal.newBCData('ImagePointCoordinates',parent=bcdataset)
    imagePointsFieldNode[2].append(['CoordinateX',coordsPI[1][0,:], [], 'DataArray_t'])
    imagePointsFieldNode[2].append(['CoordinateY',coordsPI[1][1,:], [], 'DataArray_t'])
    imagePointsFieldNode[2].append(['CoordinateZ',coordsPI[1][2,:], [], 'DataArray_t'])
    wallPointsFieldNode = Internal.newBCData('WallPointCoordinates',parent=bcdataset)
    wallPointsFieldNode[2].append(['CoordinateX',coordsPW[1][0,:], [], 'DataArray_t'])
    wallPointsFieldNode[2].append(['CoordinateY',coordsPW[1][1,:], [], 'DataArray_t'])
    wallPointsFieldNode[2].append(['CoordinateZ',coordsPW[1][2,:], [], 'DataArray_t'])
    return None

def _addIBDataZSR(z, correctedPts, wallPts, imagePts=None, prefix='IBCD_'):
    zname = Internal.getName(z)
    nameSubRegion = prefix+zname
    zsr = Internal.getNodeFromName1(z, nameSubRegion)
    if zsr is None:
        Internal._createChild(z, nameSubRegion, 'ZoneSubRegion_t', value=zname)
        zsr = Internal.getNodeFromName1(z, nameSubRegion)

    coordsPC = Converter.extractVars(correctedPts,['CoordinateX','CoordinateY','CoordinateZ'])[0]
    zsr[2].append(['CoordinateX_PC',coordsPC[1][0,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateY_PC',coordsPC[1][1,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateZ_PC',coordsPC[1][2,:], [], 'DataArray_t'])

    coordsPW = Converter.extractVars(wallPts, ['CoordinateX','CoordinateY','CoordinateZ'])[0]
    zsr[2].append(['CoordinateX_PW',coordsPW[1][0,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateY_PW',coordsPW[1][1,:], [], 'DataArray_t'])
    zsr[2].append(['CoordinateZ_PW',coordsPW[1][2,:], [], 'DataArray_t'])

    if imagePts is not None:
        coordsPI = Converter.extractVars(imagePts, ['CoordinateX','CoordinateY','CoordinateZ'])[0]
        zsr[2].append(['CoordinateX_PI',coordsPI[1][0,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateY_PI',coordsPI[1][1,:], [], 'DataArray_t'])
        zsr[2].append(['CoordinateZ_PI',coordsPI[1][2,:], [], 'DataArray_t'])
    return None

#=============================================================================
# Extract info for skin post-processing
# INPUT : numpys of coordinates and fields to be projected onto the surface
# IN/OUT: surface defined by a CGNS/Python tree tb
#=============================================================================
def extractIBMWallFields(XCP, YCP, ZCP, arrayOfFields, tb, variables):
    VARLIST = Converter.getVarNames(arrayOfFields)
    dictOfVarNumber={}
    for var in variables:
        for nov in range(len(VARLIST)):
            if VARLIST[nov]==var:
                dictOfVarNumber[var]=nov
                break

    # 1. Creation of a CGNS zone O-D of cloud points
    zsize = numpy.empty((1,3), numpy.int32, order='F')
    zsize[0,0] = XCP.shape[0]; zsize[0,1] = 0; zsize[0,2] = 0
    z = Internal.newZone(name='IBW_Wall',zsize=zsize,ztype='Unstructured')
    gc = Internal.newGridCoordinates(parent=z)
    coordx = ['CoordinateX',XCP,[],'DataArray_t']
    coordy = ['CoordinateY',YCP,[],'DataArray_t']
    coordz = ['CoordinateZ',ZCP,[],'DataArray_t']
    gc[2] = [coordx,coordy,coordz]
    n = Internal.createChild(z, 'GridElements', 'Elements_t', [2,0])
    Internal.createChild(n, 'ElementRange', 'IndexRange_t', [1,0])
    Internal.createChild(n, 'ElementConnectivity', 'DataArray_t', None)
    FSN = Internal.newFlowSolution(name=Internal.__FlowSolutionNodes__,
                                   gridLocation='Vertex', parent=z)

    for varname in dictOfVarNumber:
        novar = dictOfVarNumber[varname]
        vararrayN = arrayOfFields[1][novar]
        FSN[2].append([varname,vararrayN, [],'DataArray_t'])



    dimPb = Internal.getNodeFromName(tb,'EquationDimension')
    if dimPb is None:
        print('Warning: extractIBMWallFields: pb dimension is set to 3.')
        dimPb = 3
    else:
        dimPb = Internal.getValue(dimPb)
    # Force all the zones to be in a single CGNS base
    td = Internal.copyRef(tb)
    for nob in range(len(td[2])):
        b = td[2][nob]
        if b[3] == 'CGNSBase_t':
            zones = Internal.getNodesFromType1(b, 'Zone_t')
            if zones != []:
                zones = C.convertArray2Tetra(zones)
                zones = T.join(zones); zones = G.close(zones)
                b[2] = [zones]
    for varname in dictOfVarNumber:
        C._initVars(td,varname,0.)

    td = P.projectCloudSolution(z, td, dim=dimPb)
    return td

def extractIBMInfo(t):
    tinfo = C.newPyTree(['Wall','IP','Image'])
    tw = createIBMZones(t,typeOfPoint='Wall')
    tw = C.convertArray2Node(tw)
    tinfo[2][1][2] = Internal.getZones(tw)

    tw = createIBMZones(t,typeOfPoint='IP')
    tw = C.convertArray2Node(tw)
    tinfo[2][2][2] = Internal.getZones(tw)

    tw = createIBMZones(t,typeOfPoint='Image')
    tw = C.convertArray2Node(tw)
    tinfo[2][3][2] = Internal.getZones(tw)

    return tinfo
# --------------------------------------------------------------------------------
# Creation of 0-D zones of name 'Zone#IBCD_*' such that their original zones can be
# retrieved in post processing
#typeOfPoint can be 'Wall','Image','IP'
def createIBMZones(tc,variables=[], typeOfPoint='Wall'):
    if typeOfPoint=='Wall': suff='PW'
    elif typeOfPoint=='Image': suff='PI'
    elif typeOfPoint=='IP': suff='PC'
    tw = C.newPyTree(['IBM_WALL'])
    for z in Internal.getZones(tc):
        ZSR = Internal.getNodesFromType2(z,'ZoneSubRegion_t')
        for IBCD in Internal.getNodesFromName(ZSR,"IBCD_*"):
            xPW = Internal.getNodesFromName(IBCD,"CoordinateX_%s"%suff)[0][1]
            yPW = Internal.getNodesFromName(IBCD,"CoordinateY_%s"%suff)[0][1]
            zPW = Internal.getNodesFromName(IBCD,"CoordinateZ_%s"%suff)[0][1]
            nptsW = xPW.shape[0]
            zw = G.cart((0,0,0),(1,1,1),(nptsW,1,1))
            COORDX = Internal.getNodeFromName2(zw,'CoordinateX'); COORDX[1]=xPW
            COORDY = Internal.getNodeFromName2(zw,'CoordinateY'); COORDY[1]=yPW
            COORDZ = Internal.getNodeFromName2(zw,'CoordinateZ'); COORDZ[1]=zPW
            if variables != [] and variables is not None:
                FSN = Internal.newFlowSolution(parent=zw)
                for varo in variables:
                    fieldV = Internal.getNodeFromName2(IBCD,varo)
                    if fieldV is not None:
                        C._initVars(zw,varo,0.)
                        fieldW = Internal.getNodeFromName2(FSN,varo)
                        fieldW[1] = fieldV[1]

            zw[0]=z[0]+"#"+IBCD[0]+"_%s"%suff
            tw[2][1][2].append(zw)
    return tw

# Create the 1to1 connectivity of two Cartesian zones as a PointList/PointListDonor info
# convertOnly = True : abutting 1to1 gc nodes already exist
def _create1To1Connectivity(t, tskel=None, dim=3, convertOnly=True):
    dictOfDnrZoneDims={}
    if tskel is None or not convertOnly:
        Cmpi._addBXZones(t, depth=3)
        t = X.connectMatch(t, dim=dim)
        # dictionary of donor zone dimensions for each gc
        for z in Internal.getZones(t):
            for gc in Internal.getNodesFromType2(z,'GridConnectivity1to1_t'):
                znamed = Internal.getValue(gc)
                zd = Internal.getNodeFromName2(t,znamed)
                SOD = Internal.getNodeFromName2(zd,'.Solver#ownData')
                if SOD is not None:
                    l2g = Internal.getNodeFromName1(SOD,'loc2glob')
                    if l2g is not None:
                        l2g = Internal.getValue(l2g)
                        dictOfDnrZoneDims[gc[0]]=[l2g[-3],l2g[-2],l2g[-1]]
                else:
                    Internal._rmNodesFromName(z,gc[0])
        Cmpi._rmBXZones(t)
    else:
        for z in Internal.getZones(t):
            for gc in Internal.getNodesFromType2(z,'GridConnectivity1to1_t'):
                znamed = Internal.getValue(gc)
                zd = Internal.getNodeFromName2(tskel,znamed)
                if zd is not None:
                    l2g = Internal.getZoneDim(zd)
                    dictOfDnrZoneDims[gc[0]]=[l2g[-3],l2g[-2],l2g[-1]]

    for z in Internal.getZones(t):
        dimZ = Internal.getZoneDim(z)
        ni = dimZ[1]; nj = dimZ[2]; nk = dimZ[3]
        for gc in Internal.getNodesFromType2(z,'GridConnectivity1to1_t'):
            PR = Internal.getNodeFromName1(gc,'PointRange')
            PRD = Internal.getNodeFromName1(gc,'PointRangeDonor')
            win = Internal.range2Window(PR[1])
            imin = win[0]; imax = win[1]
            jmin = win[2]; jmax = win[3]
            kmin = win[4]; kmax = win[5]
            indicesL = Converter.converter.range2PointList(imin,imax,jmin,jmax,kmin,kmax,
                                                           ni, nj, nk)
            if isinstance(indicesL, numpy.ndarray): r = indicesL
            else:
                r = numpy.array(indicesL, dtype=numpy.int32)
            r = r.reshape((1,r.size), order='F')
            gc[2].append(["PointList", r, [], "IndexArray_t"])
            if gc[0] not in dictOfDnrZoneDims:
                if rank==0:  print("NOT FOUND", z[0], Internal.getValue(gc), gc[0])
            else:
                dimZd = dictOfDnrZoneDims[gc[0]]
                nid = dimZd[0]; njd = dimZd[1]; nkd = dimZd[2]
                win = Internal.range2Window(PRD[1])
                imin = win[0]; imax = win[1]
                jmin = win[2]; jmax = win[3]
                kmin = win[4]; kmax = win[5]
                indicesOppL = Converter.converter.range2PointList(imin, imax, jmin, jmax, kmin, kmax, nid, njd, nkd)

                if isinstance(indicesOppL, numpy.ndarray): r = indicesOppL
                else:
                    r = numpy.array(indicesOppL, dtype=numpy.int32)
                r = r.reshape((1,r.size), order='F')
                gc[2].append(["PointListDonor", r, [], "IndexArray_t"])
                Internal._rmNodesFromName(gc, "PointRange")
                Internal._rmNodesFromName(gc, "PointRangeDonor")
                Internal._rmNodesFromName(gc, 'Transform')
                Internal.setType(gc, 'GridConnectivity_t')
                Internal._createChild(gc, 'GridConnectivityType', 'GridConnectivityType_t', value='Abutting1to1')
                Internal._createChild(gc, 'GridLocation', 'GridLocation_t', value='FaceCenter')
        C._mergeGCs(z)

    return None

def _addExternalBCs(t, bbox, externalBCType='BCFarfield', dimPb=3):
    dirs = [0,1,2,3,4,5]
    rangeDir=['imin','jmin','kmin','imax','jmax','kmax']
    if dimPb == 2: dirs = [0,1,3,4]
    for zp in Internal.getZones(t):
        dimZ = Internal.getZoneDim(zp)
        niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
        indM = niz-1+(njz-1)*niz+(nkz-1)*niz*njz
        x1 = C.getValue(zp,'CoordinateX',0)
        y1 = C.getValue(zp,'CoordinateY',0)
        z1 = C.getValue(zp,'CoordinateZ',0)
        x2 = C.getValue(zp,'CoordinateX',indM)
        y2 = C.getValue(zp,'CoordinateY',indM)
        z2 = C.getValue(zp,'CoordinateZ',indM)
        bbz=[x1,y1,z1,x2,y2,z2]
        for idir in dirs:
            if abs(bbz[idir]-bbox[idir])< 1.e-6:
                C._addBC2Zone(zp, 'external', externalBCType, rangeDir[idir])
    return None

def _setSnear(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'snear', 'DataArray_t', value)
    return None

def _setDfar(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'dfar', 'DataArray_t', value)
    return None

def _setIBCType(z, value):
    zones = Internal.getZones(z)
    for z in zones:
        Internal._createUniqueChild(z, '.Solver#define', 'UserDefinedData_t')
        n = Internal.getNodeFromName1(z, '.Solver#define')
        Internal._createUniqueChild(n, 'ibctype', 'DataArray_t', value)
    return None

def _snearFactor(t, factor=1.):
    zones = Internal.getZones(t)
    for z in zones:
        nodes = Internal.getNodesFromName2(z, 'snear')
        for n in nodes:
            Internal._setValue(n, factor*Internal.getValue(n))
    return None

####################################################
# OBSOLETE
####################################################
# OUT: list [HN_C, HN_F1, HN_F2, HN_F3, HN_F4]: indices of coarse/fine1 to 4 elements, index starts at 0
# indicesFacesOrig start at 1
def getHangingNodesInfoPara(a, b, indicesFacesOrig, bopp, indicesFacesOrigOpp):
    b[0]='extFaces'; bopp[0] = 'extFacesOpp'
    hookExtFacesOpp=C.createHook(bopp,function='nodes')
    # identify face centers of coarse (local) zone matching with a vertex of bopp (fine)
    HN = C.identifyElements(hookExtFacesOpp,b)
    cnExtFaceOpp = Internal.getNodeFromType(bopp,'Elements_t')
    cnExtFaceOpp = Internal.getNodeFromName(cnExtFaceOpp,'ElementConnectivity')
    cnExtFaceOpp = Internal.getValue(cnExtFaceOpp)
    sizeCNExtFaceOpp = cnExtFaceOpp.shape[0]
    eltType = Internal.getZoneDim(b)[3]
    if eltType == 'BAR': shiftElt=2; nfaces = 4
    else: shiftElt=4; nfaces = 6

    # loop on face centers of coarse zone (local b zone)
    HN_COARSE=[]; HN_FINE1=[]; HN_FINE2=[]; HN_FINE3=[]; HN_FINE4=[]
    if eltType=='BAR':
        for noEltEF in range(len(HN)):
            indVertexEF = HN[noEltEF] # indice of vertex of bopp, starts at 1
            if indVertexEF !=-1:
                HN_COARSE.append(int(indicesFacesOrig[noEltEF]-1)//nfaces)
                # looking for opp faces (fine side) with vertex indVertexEF
                noptr = 0; noe = 0
                found = 0; efd = -1; efg = -1
                while noptr < sizeCNExtFaceOpp:
                    indV1 = cnExtFaceOpp[noptr]
                    indV2 = cnExtFaceOpp[noptr+1]
                    if indV1 == indVertexEF:
                        efd = int(indicesFacesOrigOpp[noe]-1)//nfaces
                        found +=1
                    if indV2 == indVertexEF:
                        efg = int(indicesFacesOrigOpp[noe]-1)//nfaces
                        found += 1

                    if found == 2:
                        if efd==efg: efd = -1
                        HN_FINE2.append(efd)
                        HN_FINE1.append(efg)
                        break
                    noptr+=shiftElt
                    noe+=1
        return [HN_COARSE, HN_FINE1, HN_FINE2]

    else:
        for noEltEF in range(len(HN)):
            indVertexEF = HN[noEltEF] # starts at 1
            if indVertexEF !=-1:
                HN_COARSE.append(int(indicesFacesOrig[noEltEF]-1)//nfaces)
                # looking for opp faces (fine side) with vertex indVertexEF
                noptr = 0; noe = 0
                found = 0
                while noptr < sizeCNExtFaceOpp:
                    indV1 = cnExtFaceOpp[noptr]
                    indV2 = cnExtFaceOpp[noptr+1]
                    indV3 = cnExtFaceOpp[noptr+2]
                    indV4 = cnExtFaceOpp[noptr+3]
                    ef1=-1; ef2=-1; ef3=-1; ef4=-1
                    if indVertexEF == indV1:
                        ef4 = int(indicesFacesOrigOpp[noe]-1)//nfaces
                        found +=1

                    elif indVertexEF == indV2:
                        ef3 = int(indicesFacesOrigOpp[noe]-1)//nfaces
                        found += 1

                    elif indVertexEF == indV3:
                        ef1 = int(indicesFacesOrigOpp[noe]-1)//nfaces
                        found += 1

                    elif indVertexEF == indV4:
                        ef2 = int(indicesFacesOrigOpp[noe]-1)//nfaces
                        found += 1

                    if found == 4:
                        if ef2 == ef1: ef2 =-1
                        if ef3 == ef1: ef3 =-1
                        if ef4 == ef1: ef4 =-1
                        if ef2 == ef3: ef3 =-1
                        if ef3 == ef4: ef4 =-1
                        if ef4 == ef2: ef4 =-1
                        HN_FINE1.append(ef1)
                        HN_FINE2.append(ef2)
                        HN_FINE3.append(ef3)
                        HN_FINE4.append(ef4)

                        break
                    noptr+=shiftElt
                    noe+=1
        return [HN_COARSE, HN_FINE1, HN_FINE2, HN_FINE3, HN_FINE4]

def getHangingNodesInfoSeq(a):
    indicesFacesOrig = []
    b = P.exteriorFaces(a,indices=indicesFacesOrig)
    indicesFacesOrig=indicesFacesOrig[0]# index starts at 0
    b[0]='extfaces'
    return getHangingNodesInfoPara(a, b, indicesFacesOrig, b, indicesFacesOrig)

#=================================================================================
# Fonctions patchees de Cassiopee Converter.PyTree disponibles a partir de la v3.4
#=================================================================================
# -- mergeConnectivity
# IN: z1: zone BE
# IN: z2: zone BE (to be merged in z1) avec un seul type d'elements
# si boundary==1, les noeuds de z2 sont identifies dans z1
# si boundary==0, les noeuds de z2 sont merges dans z1 et reidentifies
# IN: boundary: 0 (not a boundary zone), 1 (a boundary zone, add it as a
# boundary connectivity)
def _mergeConnectivityLoc(z1, z2, boundary=0, shared=False):
    # Analyse zone z2
    dims = Internal.getZoneDim(z2)
    neb = dims[2] # nbre d'elts de z2
    eltType, nf = Internal.eltName2EltNo(dims[3]) # type d'elements de z2

    # On cherche l'element max dans les connectivites de z1
    maxElt = 0
    connects = Internal.getNodesFromType(z1, 'Elements_t')
    for cn in connects:
        r = Internal.getNodeFromName1(cn, 'ElementRange')
        m = r[1][1]
        maxElt = max(maxElt, m)

    # connectivite ajoutee=volumique non shared
    if boundary == 0 and not shared:
        # on fusionne les coordonnees
        import Transform.PyTree as T
        zn1 = C.convertArray2Node(z1)
        zn2 = C.convertArray2Node(z2)
        zn = T.join(zn1, zn2)
        # reset Coordinates in z1 with merged zn coordinates
        cont1 = Internal.getNodeFromName1(z1, Internal.__GridCoordinates__)
        contn = Internal.getNodeFromName1(zn, Internal.__GridCoordinates__)
        for name in ['CoordinateX', 'CoordinateY', 'CoordinateZ']:
            p1 = Internal.getNodeFromName1(cont1, name)
            pn = Internal.getNodeFromName1(contn, name)
            p1[1] = pn[1]
        # Recupere le container noeud
        cont1 = Internal.getNodeFromName1(z1, Internal.__FlowSolutionNodes__)
        contn = Internal.getNodeFromName1(zn, Internal.__FlowSolutionNodes__)
        if contn is not None:
            for n in contn[2]:
                if n[3] == 'DataArray_t':
                    p1 = Internal.getNodeFromName1(cont1, n[0])
                    if p1 is not None: p1[1] = n[1]

        # Nouveau nbre de points dans z1
        np = Internal.getZoneDim(zn)[1]
        z1[1] = numpy.copy(z1[1])
        z1[1][0,0] = np

        # Reidentifie les connectivites de z1
        hook = C.createHook(zn, 'nodes')
        ids = C.identifyNodes(hook, z1)
        nodes = Internal.getNodesFromType1(z1, 'Elements_t')
        for n in nodes:
            node = Internal.getNodeFromName1(n, 'ElementConnectivity')
            oldc = node[1]
            newc = numpy.copy(oldc)
            newc[:] = ids[oldc[:]-1]
            node[1] = newc

        # Ajoute les FlowSolutions en centres
        cont1 = Internal.getNodeFromName1(z1, Internal.__FlowSolutionCenters__)
        cont2 = Internal.getNodeFromName1(z2, Internal.__FlowSolutionCenters__)
        if cont2 is not None:
            for n in cont2[2]:
                if n[3] == 'DataArray_t':
                    p1 = Internal.getNodeFromName1(cont1, n[0])
                    if p1 is not None: p1[1] = numpy.concatenate((p1[1],n[1]))

        # Ajoute les connectivites de z2
        nebb = 0
        z1[1][0,1] += neb # nouveau nbre de cellules

        # on identifie les noeuds connectivity de z2 dans zn
        ids = C.identifyNodes(hook, z2)
        # On ajoute les connectivites de z2
        nodes = Internal.getNodesFromType1(z2, 'Elements_t')
        for n in nodes:
            node = Internal.createUniqueChild(z1, n[0]+'-2', 'Elements_t', value=n[1])
            rn = Internal.getNodeFromType1(n, 'IndexRange_t')
            r = Internal.getValue(rn)
            r0 = r[0]+maxElt; r1 = r[1]+maxElt
            Internal.createUniqueChild(node, rn[0], 'IndexRange_t', value=[r0,r1])

            oldc = Internal.getNodeFromName1(n, 'ElementConnectivity')[1]
            newc = numpy.copy(oldc)
            newc[:] = ids[oldc[:]-1]
            Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)

        # decale les ranges des zoneBC de z2
        zbc1 = Internal.getNodeFromType1(z1, 'ZoneBC_t')
        zbc2 = Internal.getNodesFromType2(z2, 'BC_t')
        if zbc1 is None and zbc2 != []: zbc1 = Internal.newZoneBC(parent=z1)
        for bc in zbc2:
            rn = Internal.getNodeFromType1(bc, 'IndexRange_t')
            r = Internal.getValue(rn)
            r0 = r[0,0]+maxElt; r1 = r[0,1]+maxElt
            n = Internal.createUniqueChild(zbc1, bc[0], 'BC_t', value=bc[1])
            Internal.createUniqueChild(n, rn[0], rn[3], value=[[r0,r1]])

        # Tri des connectivites
        Internal._sortNodesInZone(z1)

    elif boundary == 0 and shared: # connectivite ajoutee=volumique shared
        # on cree des nouveaux noeuds connectivites dans z1
        elts = Internal.getNodesFromType2(z2, 'Elements_t')
        z1[1][0,1] += neb # update le nbre d'elements de z1
        nebb = 0
        for e in elts:
            if e[1][1] == 0: # volumique uniquement
                r = Internal.getNodeFromName1(e, 'ElementRange')[1]
                nbe2 = r[1]-r[0]+1
                e2 = Internal.createUniqueChild(z1, e[0], 'Elements_t', value=[eltType,0])
                Internal.createUniqueChild(e2, 'ElementRange', 'IndexRange_t',
                                           value=[maxElt+nebb+1,maxElt+nebb+nbe2])
                newc = Internal.getNodeFromName1(e, 'ElementConnectivity')[1]
                Internal.createUniqueChild(e2, 'ElementConnectivity', 'DataArray_t', value=newc)
                nebb += nbe2

        # Ajoute les FlowSolutions en centres
        cont1 = Internal.getNodeFromName1(z1, Internal.__FlowSolutionCenters__)
        cont2 = Internal.getNodeFromName1(z2, Internal.__FlowSolutionCenters__)
        if cont2 is not None:
            for n in cont2[2]:
                if n[3] == 'DataArray_t':
                    p1 = Internal.getNodeFromName1(cont1, n[0])
                    if p1 is not None: p1[1] = numpy.concatenate((p1[1],n[1]))

        # Tri des connectivites
        Internal._sortNodesInZone(z1)


    else: # connectivite ajoutee=boundary (subzone)
        # on identifie les noeuds de z2 dans z1
        z1ext = P.exteriorFaces(z1)
        hook = C.createHook(z1ext, 'centers')
        ids = C.identifyNodes(hook, z2)
        ids = ids[ids[:]>-1]
        neb = ids.shape[0]

        z1[1] = numpy.copy(z1[1])
        z1[1][0,2] += neb

        # on cree un nouveau noeud connectivite dans z1 (avec le nom de la zone z2)
        nebb = neb
        node = Internal.createUniqueChild(z1, z2[0], 'Elements_t', value=[eltType,nebb])
        Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t',
                                   value=[maxElt+1,maxElt+neb])
        #oldc = Internal.getNodeFromName2(z2, 'ElementConnectivity')[1]
        #newc = numpy.copy(oldc)
        newc = numpy.zeros((neb),dtype=numpy.int32)
        newc[:] = ids[:]-1

        Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)
    return None

# fonction de Converter.PyTree patchee en attendant la v3.4 de Cassiopee
def _addBC2UnstructZoneLoc__(z, bndName, bndType, wrange=[],
                             zoneDonor=[], rangeDonor=[], trirac=[1,2,3],
                             rotationCenter=[], rotationAngle=[], translation=[],
                             faceList=[], pointList=[], elementList=[], elementRange=[], data=None,
                             subzone=None, faceListDonor=None, elementListDonor=None,
                             elementRangeDonor=None, tol=1.e-12, unitAngle=None):
    s = bndType.split(':')
    bndType1 = s[0]
    if len(s) > 1: bndType2 = s[1]
    else: bndType2 = ''
    # si subzone: on cree la connectivite BC, on passe en elementRange
    if subzone is not None and pointList == []:
        bcn = Internal.getNodeFromName1(z, subzone[0])
        if bcn is None:
            _mergeConnectivityLoc(z, subzone, boundary=1)
        bcn = Internal.getNodeFromName1(z, subzone[0])
        bcnr = Internal.getNodeFromName1(bcn, 'ElementRange')
        elementRange = [bcnr[1][0], bcnr[1][1]]

    # si subzone + pointList=True, on identifie la subzone en pointList
    if subzone is not None and pointList == True:
        hook = C.createHook(z, function='nodes')
        pointList = C.identifyNodes(hook, subzone)
        C.freeHook(hook)


    if bndType1 == 'BCMatch' or bndType1 == 'Abutting1to1':
        if (zoneDonor == [] or
                faceListDonor is None and subzone is None and elementListDonor is None and elementRangeDonor is None):
            raise ValueError("addBC2Zone: unstructured match connectivity requires a donor zone and a faceListDonor or a subzone or an elementRangeDonor or an elementListDonor.")
        # si subzone fournie: cree le elementRangeDonor
        if subzone is not None:
            bcn = Internal.getNodeFromName1(zoneDonor, subzone[0])
            if bcn is None:
                C._mergeConnectivity(zoneDonor, subzone, boundary=1)
            bcn = Internal.getNodeFromName1(zoneDonor, subzone[0])
            bcnr = Internal.getNodeFromName1(bcn, 'ElementRange')
            elementRangeDonor = [bcnr[1][0], bcnr[1][1]]

        # Cree le noeud zoneGridConnectivity si besoin
        zoneGC = Internal.createUniqueChild(z, 'ZoneGridConnectivity',
                                            'ZoneGridConnectivity_t')

        if isinstance(zoneDonor, str): v = zoneDonor
        else: v = zoneDonor[0]
        Internal._createChild(zoneGC, bndName, 'GridConnectivity_t', value=v)
        l = len(zoneGC[2])
        info = zoneGC[2][l-1]

        if elementList != []:
            if isinstance(elementList, numpy.ndarray): r = elementList
            else: r = numpy.array(elementList, dtype=numpy.int32)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__ELEMENTLIST__, r, [], 'IndexArray_t'])
        elif elementRange != []:
            r = numpy.empty((1,2), numpy.int32, order='F')
            r[0,0] = elementRange[0]
            r[0,1] = elementRange[1]
            info[2].append([Internal.__ELEMENTRANGE__, r, [], 'IndexRange_t'])
        elif faceList != []:
            Internal._createChild(info, 'GridLocation', 'GridLocation_t', value='FaceCenter')
            if isinstance(faceList, numpy.ndarray): r = faceList
            else: r = numpy.array(faceList, dtype=numpy.int32)
            r = r.reshape((1, r.size), order='F')
            info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

        if elementListDonor is not None:
            if isinstance(elementListDonor, numpy.ndarray):
                r = elementListDonor
            else: r = numpy.array(elementListDonor, dtype=numpy.int32)
            r = r.reshape((1,r.size))
            info[2].append([Internal.__ELEMENTLIST__+'Donor', r, [], 'IndexArray_t'])
        elif elementRangeDonor is not None:
            r = numpy.empty((1,2), numpy.int32, order='F')
            r[0,0] = elementRangeDonor[0]
            r[0,1] = elementRangeDonor[1]
            info[2].append([Internal.__ELEMENTRANGE__+'Donor', r, [], 'IndexRange_t'])
        elif faceListDonor is not None:
            if isinstance(faceListDonor, numpy.ndarray): r = faceList
            else: r = numpy.array(faceListDonor, dtype=numpy.int32)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__FACELIST__+'Donor', r, [], 'IndexArray_t'])

        Internal.createChild(info, 'GridConnectivityType', 'GridConnectivityType_t', 'Abutting1to1')

        # Ajout pour les BC match periodiques
        C._addPeriodicInfoInGC__(info, rotationCenter, rotationAngle, translation, unitAngle=unitAngle)

    elif bndType1 == 'BCOverlap':
        print('Warning: addBC2Zone: BCOverlap not valid for unstructured zones.')

    elif bndType1 == 'BCNearMatch':
        print('Warning: addBC2Zone: BCNearMatch not valid for unstructured zones.')

    elif (bndType1 == 'FamilySpecified' and fnmatch.fnmatch(bndType2, 'BCStage*')) or (bndType1 == 'BCStage'):
        C._addFamilyOfStageGC__(z, bndName, bndType2, typeZone=2, elementRange=elementRange,
                                elementList=elementList, faceList=faceList, zoneDonor=zoneDonor)

    else: # BC classique
        # Cree le noeud zoneBC si besoin
        zoneBC = Internal.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')

        # Cree le noeud de BC
        info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType1)
        if elementList != []:
            if isinstance(elementList, numpy.ndarray): r = elementList
            else: r = numpy.array(elementList, dtype=numpy.int32)
            r = r.reshape((1,r.size), order='F')
            Internal.createChild(info, Internal.__ELEMENTLIST__, 'IndexArray_t', value=r)
        elif elementRange != []:
            n = numpy.empty((1,2), numpy.int32, order='F')
            n[0,0] = elementRange[0]
            n[0,1] = elementRange[1]
            Internal.createUniqueChild(info, Internal.__ELEMENTRANGE__,
                                       'IndexRange_t', value=n)
        elif faceList != []:
            Internal.createUniqueChild(info, 'GridLocation', 'GridLocation_t',
                                       value='FaceCenter')
            if isinstance(faceList, numpy.ndarray): r = faceList
            else: r = numpy.array(faceList, dtype=numpy.int32)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__FACELIST__, r, [], 'IndexArray_t'])

        elif pointList != []:
            Internal.createUniqueChild(info, 'GridLocation', 'GridLocation_t',
                                       value='Vertex')
            if isinstance(pointList, numpy.ndarray): r = pointList
            else: r = numpy.array(pointList, dtype=numpy.int32)
            r = r.reshape((1,r.size), order='F')
            info[2].append([Internal.__POINTLIST__, r, [], 'IndexArray_t'])

        # Ajoute la famille si necessaire
        if bndType1 == 'FamilySpecified':
            Internal.createChild(info, 'FamilyName', 'FamilyName_t', bndType2)

        # Ajoute les Data si necessaire (donnees Dirichlet)
        if data is not None:
            node1 = Internal.createNode('State', 'DataArray_t', value=data)
            node2 = Internal.createNode('DirichletData', 'BCData_t', children=[node1])
            node3 = Internal.createNode('BCDataSet', 'BCDataSet_t', children=[node2])
            info[2].append(node3)
    return None


def _recoverBCs1(a, T, tol=1.e-11):
    """Recover given BCs on a tree.
    Usage: _recoverBCs(a, (BCs, BCNames, BCTypes), tol)"""
    try:import Post.PyTree as P
    except: raise ImportError("_recoverBCs: requires Post module.")
    C._deleteZoneBC__(a)
    zones = Internal.getZones(a)
    (BCs, BCNames, BCTypes) = T
    for z in zones:
        indicesF = []
        try: f = P.exteriorFaces(z, indices=indicesF)
        except: continue
        indicesF = indicesF[0]
        hook = C.createHook(f, 'elementCenters')

        for c in range(len(BCs)):
            b = BCs[c]
            if b == []:
                raise ValueError("_recoverBCs: boundary is probably ill-defined.")
            # Break BC connectivity si necessaire
            elts = Internal.getElementNodes(b)
            size = 0
            for e in elts:
                erange = Internal.getNodeFromName1(e, 'ElementRange')[1]
                size += erange[1]-erange[0]+1
            n = len(elts)
            if n == 1:
                ids = C.identifyElements(hook, b, tol)
            else:
                bb = C.breakConnectivity(b)
                ids = numpy.array([], dtype=numpy.int32)
                for bc in bb:
                    ids = numpy.concatenate([ids, C.identifyElements(hook, bc, tol)])

            # Cree les BCs
            ids0 = ids # keep ids for bcdata
            ids  = ids[ids > -1]
            sizebc = ids.size
            import Transform.PyTree as TP
            zf = TP.subzone(z,ids, type='nodes')
            if sizebc > 0:
                id2 = numpy.empty(sizebc, numpy.int32)
                id2[:] = indicesF[ids[:]-1]
                C._addBC2Zone(z, BCNames[c], BCTypes[c], faceList=id2)
                # recupere la connectivite quad de la face list
                # Recupere BCDataSets
                fsc = Internal.getNodeFromName(b, Internal.__FlowSolutionCenters__)

                if fsc is not None:
                    newNameOfBC = C.getLastBCName(BCNames[c])
                    bcz = Internal.getNodeFromNameAndType(z, newNameOfBC, 'BC_t')

                    ds = Internal.newBCDataSet(name='BCDataSet', value='UserDefined',
                                               gridLocation='FaceCenter', parent=bcz)
                    d = Internal.newBCData('NeumannData', parent=ds)

                    for node in Internal.getChildren(fsc):
                        if Internal.isType(node, 'DataArray_t'):
                            val0 = Internal.getValue(node)
                            val0 = numpy.reshape(val0, val0.size, order='F')
                            val1 = val0[ids0>-1]
                            Internal._createUniqueChild(d, node[0], 'DataArray_t', value=val1)

        C.freeHook(hook)

    return None
