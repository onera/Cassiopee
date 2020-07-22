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
import Transform
#==============================================
# IBM prepro for CODA 
#==============================================
def prepare(t_case, t_out, vmin=5, dfarList=[], dfar=10., snears=0.01, NP=0,
            tbox=None, snearsf=None, expand=3, check=False, fileout='octree.cgns'):
    recoverBC=True
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
    C.convertPyTree2File(p,"octree_%d.cgns"%rank)

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
    allip_pts=[]; allwall_pts=[]; allimage_pts=[]

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
            ip_pts = C.getAllFields(ip_ptsZC,loc='nodes')[0]
            ip_pts = Converter.convertArray2Node(ip_pts)
            wallpts = T.projectAllDirs(ip_ptsZC, tb, varsn, oriented=0)
            C._normalize(wallpts,varsn)
            wallpts = C.convertArray2Node(wallpts)
            for var in varsn:
                C._initVars(wallpts,'{%s}=%g*{%s}'%(var,he,var))
            imagepts = T.deform(wallpts, vector=varsn)
            wallpts = C.getFields(Internal.__GridCoordinates__,wallpts)[0]
            imagepts = C.getFields(Internal.__GridCoordinates__,imagepts)[0]
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
        C.convertPyTree2File(front1,'targetFaces_%d.cgns'%(Cmpi.rank))

    zones = Internal.getZones(t)
    z = zones[0]
    for noz in range(1,len(zones)):
        z = T.join([z,zones[noz]])
    z[0] = 'CART_P%d'%rank

    # recoverIBC familySpecified:IBMWall
    if front1 != []: # can be empty on a proc
        front1[0]='IBMWall'    
        C._addBC2Zone(z, 'IBMWall', 'FamilySpecified:IBMWall',subzone=front1)
        allip_pts=Transform.join(allip_pts)
        allwall_pts = Transform.join(allwall_pts)
        allimage_pts = Transform.join(allimage_pts)
        _addIBDataZSR(z,[allip_pts],[allwall_pts],[allimage_pts], prefix='IBCD_')

    # add 2 PyTree
    Cmpi._setProc(z,Cmpi.rank)
    t[2][1][2]=[z]
    
    
    # identify elements per processor
    if t_out is not None: Cmpi.convertPyTree2File(t, t_out)

    return t

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
            if variables != [] and variables != None:
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

#------------------------------------------
# Creation of the subregions for IBM zones
#-----------------------------------------
#def _createIB_ZSR(z, facelist, correctedPts, wallPts, imagePts, bctype, loc='faces'):
#    varsRM = ["Density","utau","yplus","Pressure","Velocity*"]
#    zname = Internal.getName(z)
#    nameSubRegion='IBCD_'+zname
#    zsr = Internal.getNodesFromName1(z, nameSubRegion)
#    # create new subregion for interpolations
#    if zsr == []:
#        dimZSR = numpy.zeros((1), numpy.int32)
#        dimZSR[0] = correctedPts[1].shape[0]
#        #v = numpy.fromstring(zname, 'c')
#        z[2].append([nameSubRegion, dimZSR, [],'ZoneSubRegion_t'])
#        info = z[2][len(z[2])-1]
#        info[2].append(['PointList',facelist, [], 'IndexArray_t'])
#        Internal.createChild(info,'GridLocation','GridLocation_t','FaceCenter')
#
#    XOD._addIBCCoords__(z, z[0], correctedPts, wallPts, imagePts, bctype)
#    # remove some nodes created by addIBCCoord but not useful for CODA
#    zsr=Internal.getNodesFromName1(z,nameSubRegion)[0]
#    for vrm in varsRM: Internal._rmNodesFromName(zsr,vrm)
#    return None
