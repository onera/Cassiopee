import Converter.PyTree as C
import Converter.Internal as Internal
import Connector.ToolboxIBM as TIBM
import Transform.PyTree as T
import Converter.Internal as Internal
import Converter
import Generator.PyTree as G
import Post.PyTree as P
import Connector.OversetData as XOD
import Converter.Mpi as Cmpi
import Connector.connector as connector
import numpy
import Dist2Walls.PyTree as DTW
import KCore.test as test
import Connector.PyTree as X
import Connector.Mpi as Xmpi

#==============================================
# IBM prepro for CODA 
#==============================================
def prepare(t_case, t_out, vmin=5, dfarList=[], dfar=10., snears=0.01, NP=0,
            tbox=None, snearsf=None, expand=3, check=False, fileout='octree.cgns'):
    isHanging=False
    symmetry=0
    dfarDir=0
    IBCType=1
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
                else: snearf.append(1.)

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

    o = TIBM.buildOctree(tb, snears=snears, snearFactor=1., dfar=dfar, dfarList=dfarList,
                         to=to, tbox=tbox, snearsf=snearsf,
                         dimPb=dimPb, vmin=vmin, symmetry=symmetry, fileout=None, rank=rank,
                         expand=expand, dfarDir=dfarDir)

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
    # C.convertPyTree2File(p,"octree_%d.cgns"%rank)

    # fill vmin + merge in parallel
    test.printMem(">>> Octree struct [start]")
    res = TIBM.octree2StructLoc__(p, vmin=vmin, ext=-1, optimized=0, parento=parento, sizeMax=1000000)
    del p
    if parento is not None:
        for po in parento: del po
    t = C.newPyTree(['CARTESIAN', res])
    zones = Internal.getZones(t)
    for z in zones: 
        zname = Internal.getName(z)
        zname = zname+'X%d'%rank
        Internal.setName(z,zname)
    Cmpi._setProc(t, rank)

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
    t = P.computeGrad(t, 'TurbulentDistance')
    t = C.center2Node(t,["centers:gradxTurbulentDistance",'centers:gradyTurbulentDistance','centers:gradzTurbulentDistance'])
    #print("After computeGrad : Perform transfers of gradient correctly ????")
    Internal._rmNodesFromName(t,Internal.__FlowSolutionCenters__)
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
    varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
    ibczones = []
    fullzones = []
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
            ibczones.append(z)
            #
            # IBM Wall points
            #
            zname = ip_ptsZ[0]
            ip_pts = C.getAllFields(ip_ptsZC,loc='nodes')
            ip_pts = Converter.convertArray2Node(ip_pts)
            wallpts = T.projectAllDirs(ip_ptsZC, tb, varsn, oriented=0)
            C._normalize(wallpts,varsn)
            wallpts = C.convertArray2Node(wallpts)
            for var in varsn:
                C._initVars(wallpts,'{%s}=%g*{%s}'%(var,he,var))
            imagepts = T.deform(wallpts, vector=varsn)
            wallpts = C.getFields(Internal.__GridCoordinates__,wallpts)
            imagepts = C.getFields(Internal.__GridCoordinates__,imagepts)
            # Converter.convertArrays2File(ip_pts,"targetPts_%s.plt"%zname)
            # Converter.convertArrays2File(wallpts,"wallPts_%s.plt"%zname)
            # Converter.convertArrays2File(imagepts,"imagePts_%s.plt"%zname)
            # indices of elements at border of IP Faces
            ip_ptsZ[0]+='_IPFaces'
 #           C._addBC2Zone(z, 'IBMWall', 'FamilySpecified:IBMWall',subzone=ip_ptsZ)
            VARCOORDS = ['CoordinateX','CoordinateY','CoordinateZ']
            ip_pts = Converter.extractVars(ip_pts,VARCOORDS)
            ip_pts[0]='xc,yc,zc'
            wallpts = Converter.extractVars(wallpts,VARCOORDS) 
            wallpts[0]='xw,yw,zw'
            imagepts = Converter.extractVars(imagepts,VARCOORDS) 
            imagepts[0]='xi,yi,zi'
            IBPts=Converter.addVars([ip_pts,wallpts,imagepts])
            C.setFields([IBPts], ip_ptsZ,loc='centers',writeDim=False)
            # bcs = Internal.getNodesFromType(z,'BC_t')
            # for bc in bcs:
            #     FamName = Internal.getNodeFromName(bc,'FamilyName')
            #     if FamName is not None and not check:
            #         if Internal.getValue(FamName)=='IBMWall':
            #             _addIBCDataSet(bc,ip_pts, wallpts, imagepts)
 
        else:
            if C.getMinValue(z,'cellN')==0.:
                z = P.selectCells(z,"{cellN}==1.",strict=1)
                C._rmVars(z,['TurbulentDistance','cellN'])
            else:
                C._rmVars(z,['TurbulentDistance','cellN'])
                z = C.convertArray2Hexa(z)
            fullzones.append(z)

        parentz[2][noz] = z

    C.convertPyTree2File(ibczones,"ibczones_%d.cgns"%(Cmpi.rank))
    C.convertPyTree2File(fullzones,"fullzones_%d.cgns"%(Cmpi.rank))

    if check:
        front1=[]
        for zname in frontDict:
            if frontDict[zname] != []:
                front1.append(frontDict[zname])
        if front1 != []: front1 = T.join(front1)
        C.convertPyTree2File(front1,'targetFaces_%d.cgns'%(Cmpi.rank))

    zones = Internal.getZones(t)
    z = zones[0]
    for noz in range(1,len(zones)):
        z = T.join([z,zones[noz]])
    z[0] = 'CART_P%d'%rank
    t[2][1][2]=[z]

    # identify hanging nodes
    
    # identify elements per processor
    if t_out is not None:
        Cmpi.convertPyTree2File(t,t_out)

    return t

# OUT: list [HN_C, HN_F1, HN_F2, HN_F3, HN_F4] : indices of coarse/fine1 to 4 elements, index starts at 0
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
    if eltType=='BAR': 
        shiftElt=2; nfaces = 4
    else: 
        shiftElt=4; nfaces = 6

    # loop on face centers of coarse zone (local b zone)
    HN_COARSE=[]; HN_FINE1 = []; HN_FINE2 = []; HN_FINE3 = []; HN_FINE4 = []
    if eltType=='BAR':
        for noEltEF in range(len(HN)):
            indVertexEF = HN[noEltEF] # indice of vertex of bopp, starts at 1
            if indVertexEF !=-1:
                HN_COARSE.append((indicesFacesOrig[noEltEF]-1)//nfaces)
                # looking for opp faces (fine side) with vertex indVertexEF
                noptr = 0; noe = 0
                found = 0
                while noptr < sizeCNExtFaceOpp:
                    indV1 = cnExtFaceOpp[noptr]
                    indV2 = cnExtFaceOpp[noptr+1]
                    if indV1 == indVertexEF: 
                        efd = (indicesFacesOrigOpp[noe]-1)//nfaces
                        HN_FINE2.append(efd)
                        found +=1
                    if indV2 == indVertexEF:
                        efg = (indicesFacesOrigOpp[noe]-1)//nfaces
                        HN_FINE1.append(efg)
                        found += 1

                    if found == 2:
                        break
                    noptr+=shiftElt
                    noe+=1
        return [HN_COARSE, HN_FINE1, HN_FINE2] 

    else:
        for noEltEF in range(len(HN)):
            indVertexEF = HN[noEltEF] # starts at 1
            if indVertexEF !=-1:
                HN_COARSE.append((indicesFacesOrig[noEltEF]-1)//nfaces)
                # looking for opp faces (fine side) with vertex indVertexEF
                noptr = 0; noe = 0
                found = 0
                while noptr < sizeCNExtFaceOpp:
                    indV1 = cnExtFaceOpp[noptr]
                    indV2 = cnExtFaceOpp[noptr+1]
                    indV3 = cnExtFaceOpp[noptr+2]
                    indV4 = cnExtFaceOpp[noptr+3]
                    if indVertexEF == indV1: 
                        efd = (indicesFacesOrigOpp[noe]-1)//nfaces
                        HN_FINE4.append(efd)
                        found +=1

                    elif indVertexEF == indV2:
                        efg = (indicesFacesOrigOpp[noe]-1)//nfaces
                        HN_FINE3.append(efg)
                        found += 1

                    elif indVertexEF == indV3:
                        efg = (indicesFacesOrigOpp[noe]-1)//nfaces
                        HN_FINE1.append(efg)
                        found += 1

                    elif indVertexEF == indV4:
                        efg = (indicesFacesOrigOpp[noe]-1)//nfaces
                        HN_FINE2.append(efg)
                        found += 1

                    if found == 4: break
                    noptr+=shiftElt
                    noe+=1
        return [HN_COARSE, HN_FINE1, HN_FINE2, HN_FINE3, HN_FINE4] 
    
def getHangingNodesInfoSeq(a):
    indicesFacesOrig = []
    b = P.exteriorFaces(a,indices=indicesFacesOrig)
    indicesFacesOrig=indicesFacesOrig[0]# index starts at 0
    b[0]='extfaces'
    return getHangingNodesInfoPara(a, b, indicesFacesOrig, b, indicesFacesOrig)

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
        v = numpy.fromstring(zname, 'c')
        z[2].append([nameSubRegion, v, [],'ZoneSubRegion_t'])
        info = z[2][len(z[2])-1]
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
   
# --------------------------------------------------------------------------------
# Creation of 0-D zones of name 'Zone#IBCD_*' such that their original zones can be
# retrieved in post processing
# Fonction is available in 
def createIBMWZones(tc,variables=[]):
    tw = C.newPyTree(['IBM_WALL'])
    for z in Internal.getZones(tc):
        ZSR = Internal.getNodesFromType2(z,'ZoneSubRegion_t')
        for IBCD in Internal.getNodesFromName(ZSR,"IBCD_*"):
            xPW = Internal.getNodesFromName(IBCD,"CoordinateX_PW")[0][1]
            yPW = Internal.getNodesFromName(IBCD,"CoordinateY_PW")[0][1]
            zPW = Internal.getNodesFromName(IBCD,"CoordinateZ_PW")[0][1]
            nptsW = xPW.shape[0]
            zw = G.cart((0,0,0),(1,1,1),(nptsW,1,1))
            COORDX = Internal.getNodeFromName2(zw,'CoordinateX'); COORDX[1]=xPW
            COORDY = Internal.getNodeFromName2(zw,'CoordinateY'); COORDY[1]=yPW
            COORDZ = Internal.getNodeFromName2(zw,'CoordinateZ'); COORDZ[1]=zPW
            if variables != [] and variables != None:
                FSN = Internal.newFlowSolution(parent=zw)
                for varo in variables:
                    fieldV = Internal.getNodeFromName2(IBCD,varo)
                    if fieldV is not None: 
                        C._initVars(zw,varo,0.)
                        fieldW = Internal.getNodeFromName2(FSN,varo)
                        fieldW[1] = fieldV[1]

            zw[0]=z[0]+"#"+IBCD[0]
            tw[2][1][2].append(zw)
    return tw

#------------------------------------------
# Creation of the subregions for IBM zones
#-----------------------------------------
def _createIB_ZSR(z, facelist, correctedPts, wallPts, imagePts, bctype, loc='faces'):
    varsRM = ["Density","utau","yplus","Pressure","Velocity*"]
    zname = Internal.getName(z)
    nameSubRegion='IBCD_'+zname
    zsr = Internal.getNodesFromName1(z, nameSubRegion)
    # create new subregion for interpolations
    if zsr == []:
        dimZSR = numpy.zeros((1), numpy.int32)
        dimZSR[0] = correctedPts[1].shape[0]
        #v = numpy.fromstring(zname, 'c')
        z[2].append([nameSubRegion, dimZSR, [],'ZoneSubRegion_t'])
        info = z[2][len(z[2])-1]
        info[2].append(['PointList',facelist, [], 'IndexArray_t'])
        Internal.createChild(info,'GridLocation','GridLocation_t','FaceCenter')

    XOD._addIBCCoords__(z, z[0], correctedPts, wallPts, imagePts, bctype)
    # remove some nodes created by addIBCCoord but not useful for CODA
    zsr=Internal.getNodesFromName1(z,nameSubRegion)[0]
    for vrm in varsRM: Internal._rmNodesFromName(zsr,vrm)
    return None

# Extract hanging nodes information of a structured octree mesh with a nearmatch GC (coarse side)
# IN : t : Cartesian mesh with nearmatch info (coarse side)
# OUT : dictionary (key=coarse zone), values = [[znameopp,HC,HF1,HF2,HF3,HF4],...]
#￿where HC : index of element (>=0) on coarse side, HF1 to HF4 element indices on fine side 
#            znameopp : name of opposite (fine) structured zone
# in parallel mode: t must be a skeleton such that dims of opposite zones are defined even if they 
# are not on the same processor
def getHangingNodesFromNearMatch(t_skel, dimPb=3):
    zones = Internal.getZones(t_skel)
    nzones = len(zones)
    dictOfNozOfZones={}
    for noz in range(nzones):
        zname = zones[noz][0]
        dictOfNozOfZones[zname]=noz

    dictOfHN={}
    for noz in range(nzones):
        z = zones[noz]
        if Internal.getNodeFromName1(z,'XZone') is None:
            zname = Internal.getName(z)
            dictOfHN[zname]=[]
            dimZ = Internal.getZoneDim(z)
            niz = dimZ[1]; njz = dimZ[2]; nkz = dimZ[3]
            ni1 = niz-1; nj1 = njz-1; nk1 = nkz-1
            shiftr = 0; shiftd = 0
            for gc in Internal.getNodesFromType(z,'GridConnectivity_t'):
                gctype = Internal.getNodeFromType(gc,'GridConnectivityType_t')
                if gctype is not None:
                    gctype = Internal.getValue(gctype)
                    if gctype=='Abutting':
                        PR = Internal.getNodeFromName(gc,'PointRange')
                        PR = Internal.getValue(PR)
                        PR = Internal.range2Window(PR)
                        PRD = Internal.getNodeFromName(gc,'PointRangeDonor')
                        PRD = Internal.getValue(PRD)
                        PRD = Internal.range2Window(PRD)
                        NMR = Internal.getNodeFromName(gc,'NMRatio')
                        NMR = Internal.getValue(NMR)
                        zdnrname = Internal.getValue(gc)
                        if zdnrname in dictOfNozOfZones:
                            nozd = dictOfNozOfZones[zdnrname]
                            zdnr = zones[nozd]
                            dimZd = Internal.getZoneDim(zdnr)
                            nizd = dimZd[1]; njzd = dimZd[2]; nkzd = dimZd[3]
                            nid1 = nizd-1; njd1 = njzd-1; nkd1 = nkzd-1
                            nid1njd1 = nid1*njd1
                            ni1nj1=ni1*nj1
                            i1 = PR[0]; i2 = PR[1]
                            j1 = PR[2]; j2 = PR[3]
                            k1 = PR[4]; k2 = PR[5]
                            id1 = PRD[0]; id2 = PRD[1]
                            jd1 = PRD[2]; jd2 = PRD[3]
                            kd1 = PRD[4]; kd2 = PRD[5]
                            coarseList=[]
                            fineList1 = []; fineList2 = []
                            fineList3 = []; fineList4 = []

                            if dimPb == 2:
                                if NMR[0]==2:# z is coarse, zd is fine, j=cst
                                    idl = id1-1
                                    if j1 > 1: j1 = j1-1
                                    if jd1 > 1: jd1 = jd1-1
                                    for i in range(i1-1,i2-1):
                                        indr    = i    + (j1-1)*ni1 
                                        indopp1 = idl  + (jd1-1)*nid1
                                        indopp2 = idl+1+ (jd1-1)*nid1
                                        idl = idl+2
                                        coarseList.append(indr+shiftr)
                                        fineList1.append(indopp1+shiftd)
                                        fineList2.append(indopp2+shiftd)

                                elif NMR[1] == 2:
                                    jdl = jd1-1
                                    if i1 > 1: i1 = i1-1
                                    if id1 > 1: id1 = id1-1
                                    for j in range(j1-1,j2-1):
                                        indr    = i1-1  + j*ni1 
                                        indopp1 = id1-1 + jdl*nid1
                                        indopp2 = id1-1 + (jdl+1)*nid1
                                        jdl = jdl+2
                                        coarseList.append(indr+shiftr)
                                        fineList1.append(indopp1+shiftd)
                                        fineList2.append(indopp2+shiftd)


                            else: #3D
                                if NMR[0]==1 and NMR[1]==2 and NMR[2]==2:# i=cst border, z is coarse
                                    jdl = jd1-1; kdl = kd1-1
                                    if i1 > 1: i1 = i1-1
                                    if id1 > 1: id1 = id1-1
                                    id11 = id1-1; i11 = i1-1
                                    for k in range(k1-1,k2-1):
                                        for j in range(j1-1,j2-1):
                                            indr    = i11  + j*ni1+ k*ni1nj1 
                                            indopp1 = id11 +     jdl*nid1 +    kdl*nid1njd1
                                            indopp2 = id11 + (jdl+1)*nid1 +    kdl*nid1njd1
                                            indopp3 = id11 + (jdl+1)*nid1 +(kdl+1)*nid1njd1
                                            indopp4 = id11 +     jdl*nid1 +(kdl+1)*nid1njd1
                                            jdl = jdl+2; kdl = kdl+2
                                            coarseList.append(indr+shiftr)
                                            fineList1.append(indopp1+shiftd)
                                            fineList2.append(indopp2+shiftd)
                                            fineList3.append(indopp3+shiftd)
                                            fineList4.append(indopp4+shiftd)

                                elif NMR[0]==2 and NMR[1]==1 and NMR[2]==2:# j=cst border, z is coarse
                                    idl = id1-1; kdl = kd1-1
                                    if j1 > 1: j1 = j1-1
                                    if jd1 > 1: jd1 = jd1-1
                                    jd11n = (jd1-1)*nid1; j11n = (j1-1)*ni1
                                    for k in range(k1-1,k2-1):
                                        for i in range(i1-1,i2-1):
                                            indr    = i  + j11n + k*ni1nj1 
                                            indopp1 = idl   + jd11n +    kdl*nid1njd1
                                            indopp2 = idl+1 + jd11n +    kdl*nid1njd1
                                            indopp3 = idl+1 + jd11n +(kdl+1)*nid1njd1
                                            indopp4 = idl   + jd11n +(kdl+1)*nid1njd1
                                            idl = idl+2; kdl = kdl+2
                                            coarseList.append(indr+shiftr)
                                            fineList1.append(indopp1+shiftd)
                                            fineList2.append(indopp2+shiftd)
                                            fineList3.append(indopp3+shiftd)
                                            fineList4.append(indopp4+shiftd)

                                elif NMR[0]==2 and NMR[1]==2 and NMR[2]==1:# k=cst border, z is coarse
                                    idl = id1-1; jdl = jd1-1
                                    if k1 > 1: k1 = k1-1
                                    if kd1 > 1: kd1 = kd1-1
                                    kd11n = (kd1-1)*nid1njd1; k11n = (k1-1)*ni1nj1
                                    for j in range(j1-1,j2-1):
                                        for i in range(i1-1,i2-1):
                                            indr    = i  + j*ni1+ k11n
                                            indopp1 = idl   +     jdl*nid1 + kd11n
                                            indopp2 = idl+1 +     jdl*nid1 + kd11n
                                            indopp3 = idl+1 + (jdl+1)*nid1 + kd11n
                                            indopp4 = idl   + (jdl+1)*nid1 + kd11n
                                            idl = idl+2; jdl = jdl+2
                                            coarseList.append(indr+shiftr)
                                            fineList1.append(indopp1+shiftd)
                                            fineList2.append(indopp2+shiftd)
                                            fineList3.append(indopp3+shiftd)
                                            fineList4.append(indopp4+shiftd)
                            if coarseList != []:
                                HNInfoC = [zdnrname, coarseList,
                                           fineList1, fineList2, fineList3, fineList4]
                                dictOfHN[zname].append(HNInfoC)
                           
    return dictOfHN
