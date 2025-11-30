"""Parallel Toolbox for CODA."""
import Converter
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Connector.PyTree as X
import Connector.IBM as X_IBM
import Connector.Mpi as Xmpi
import Connector.connector as connector
import Transform
import Transform.PyTree as T
import XCore.PyTree as XC
import Intersector.PyTree as XOR
import Generator.PyTree as G
import Generator.Mpi as Gmpi
import Post.PyTree as P
import Dist2Walls.PyTree as DTW
import KCore.test as test
import Generator.IBM as G_IBM
import Generator.AMR as G_AMR
import Generator.Generator as Generator
import os, sys, time,copy,math
import numpy
from Converter.Internal import E_NpyInt as E_NpyInt
from .QuadratureDG import *
TOL = 1.e-9

OPT = True # distance aux noeuds uniquement - a valider !!!
def prepareAMRData(t_case, t, IBM_parameters=None, check=False, dim=3, localDir='./'):
    sym3D=False; forceAlignment=False;
    Cmpi.trace('AMR prepare IBM...start', master=True)
    t_prep_start = time.perf_counter()
    frontTypeIP = IBM_parameters["integration points"]["front type"]
    if frontTypeIP not in ["1","2"]:
        raise ValueError("FrontTypeIP not implemented: only frontTypeIP==\"1\" and \"2\" are implemented in parallel.")

    frontTypeDP = IBM_parameters["donor points"]["front type"]
    if frontTypeDP not in ["1","2"]:
        raise ValueError("FrontTypeDP not implemented: only frontTypeDP==\"1\" and \"2\" are implemented in parallel.")

    if frontTypeDP == "1":
        depth_DP = IBM_parameters["donor points"]["depth DonorPoints"]
        if depth_DP != 1:
            print("Warning: Only depth_DP=1 is implemented in parallel-AMR. Continuing with depthDP=1.")
            IBM_parameters["donor points"]["depth DonorPoints"] = 1
            depth_DP = 1

    if frontTypeIP == "1":
        depth_IP = IBM_parameters["integration points"]["depth IntegrationPoints"]
        if depth_IP != 0:
            print("Warning: Only depth_IP=0 is implemented in parallel-AMR. Continuing with depthIP=0.")
            IBM_parameters["integration points"]["depth IntegrationPoints"] = 0
            depth_IP = 0

    VPM = False
    if "method" in IBM_parameters["IBM type"].keys():
        if IBM_parameters["IBM type"]["method"] == "VPM":
            VPM = True

    if "use different front for different BCs" not in IBM_parameters["integration points"]:
        different_front_flag = True
    else:
        different_front_flag = IBM_parameters["integration points"]["use different front for different BCs"]

    if IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
        print("Warning: You are using high-order DG/DGSEM spatial discretizations in parallel. This is a development version. For a more validated and robust high-order IBM-preprocessing, switch to serial. ")

        if different_front_flag == False:
            print("Warning:High-order DG/DGSEM spatial discretizations  \n \"use different front for different BCs ==False\" is not implemented.\n Using \"use different front for different BCs==True\" instead.")

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    bbo = Gmpi.bbox(t)
    if dim == 2:
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]
        tb2 = C.initVars(tb, 'CoordinateZ', bb0[2]) # forced
        T._addkplane(tb2)
        T._contract(tb2, (0,0,0), (1,0,0), (0,1,0), dz)
    else:
        tb2 = tb

    tb2_pre = T.join(tb2)
    tb2_pre = G.close(tb2_pre)
    tb2_pre = C.newPyTree(["unstr", tb2_pre])
    t_prep_end = time.perf_counter()
    t_wdist_start = time.perf_counter()
    # Distance to IBCs (all IBCs)
    varnames = C.getVarNames(t,loc="nodes")[0]
    if "TurbulentDistance" not in varnames:
        test.printMem(">>> Wall distance nodes [start]")
        if different_front_flag == True: #True is default
            tb_WD = getBodiesForWallDistanceComputation(tb2)
            DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dim, loc='nodes')
            if not OPT: DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dim, loc='centers')
        else: #False
            DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dim, loc='nodes')
            if not OPT: DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dim, loc='centers')
        test.printMem(">>> Wall distance nodes [end]")
    else:
        if different_front_flag == False: #True is default
            test.printMem(">>> Wall distance nodes [start]")
            Internal._renameNode(t, "TurbulentDistance", "TurbulentDistanceForCFDComputation")
            DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dim, loc='nodes')
            if not OPT: DTW._distance2Walls(t, tb2, type='ortho', signed=0, dim=dim, loc='centers')
            test.printMem(">>> Wall distance nodes [end]")
        else:
            test.printMem(">>> Wall distance nodes : skipped - dist2wall is in input PyTree ")

    t_wdist_end = time.perf_counter()
    t_blank_start = time.perf_counter()
    # blanking
    test.printMem(">>> Blanking [start]")
    # Internal._rmNodesFromName(t,Internal.__FlowSolutionCenters__)
    C._initVars(t,'cellN',1.)

    t = X_IBM.blankByIBCBodies(t, tb2_pre, 'nodes', 3)

    test.printMem(">>> Blanking [end]")
    t_blank_end = time.perf_counter()

    print('Nb of Cartesian grids=%d.'%len(Internal.getZones(t)))
    npts = 0
    NCells = Cmpi.getNCells(t)
    if Cmpi.rank==0: print('Final number of cells=%5.4f millions.'%(NCells*1e-6), flush=True)
    C._initVars(t,'{TurbulentDistance}=-1.*({cellN}<1.)*{TurbulentDistance}+({cellN}>0.)*{TurbulentDistance}')

    t_storebcs_start = time.perf_counter()
    zbcs=[]; bctypes=[]; bcnames=[]
    for bc in Internal.getNodesFromType(t,'BC_t'):
        bctype = Internal.getValue(bc)
        bcname = Internal.getName(bc)
        if bctype not in bctypes :
            bctypes.append(bctype)
            bcnames.append(bcname)

    for bctype in bctypes:
        zbc = C.extractBCOfType(t,bctype)
        Internal._rmNodesByType(zbc,"FlowSolution_t")
        zbc = T.join(zbc)
        zbcs.append(zbc)
    t_storebcs_end = time.perf_counter()
    t_front_start = time.perf_counter()
    print("Extract front faces of IBM target points...")
    frontIP = computeCellNForIBMFronts(t, dim, IBM_parameters, VPM=VPM)
    t_front_end = time.perf_counter()
    t_rem_start = time.perf_counter()
    print(" Removing blanked cells...")

    t = P.selectCells(t,"{cellN}==1.",strict=1)

    Internal._rmNodesFromName(t,"FlowSolution")
    for node in Internal.getNodesFromType(t,"Elements_t"):
        if node[0] != "GridElements":
            Internal._rmNode(t,node)
    Internal._rmNodesFromType(t,"Family_t")
    t_rem_end = time.perf_counter()
    t_frontip_start = time.perf_counter()
    print("Gathering front IP..")
    frontIP = Internal.getZones(frontIP)[0]
    dimfrontIP = numpy.sum(Internal.getValue(frontIP)[0])
    frontIP_gath = Cmpi.allgatherZones(frontIP)
    frontIP_gath = C.newPyTree(["frontIP",frontIP_gath])
    frontIP_gath = T.join(frontIP_gath)
    frontIP_gath = G.close(frontIP_gath)
    for node in Internal.getNodesFromType(frontIP_gath,"Elements_t"):
        if node[0] != "GridElements":
            Internal._rmNode(frontIP_gath,node)
    if Cmpi.rank == 0:
        if check:
            print("Exporting frontIP..")
            C.convertPyTree2File(frontIP_gath,localDir+"frontIP_gath.plt")
            C.convertPyTree2File(frontIP_gath,localDir+"frontIP_gath.cgns")

    t_frontip_end = time.perf_counter()
    t_recbcs_start = time.perf_counter()
    print(" Recovering Boundary Conditions...")
    f_pytree = P.exteriorFaces(t)
    for elt_t in Internal.getNodesFromType(f_pytree, "Elements_t"):
        if not elt_t[0].startswith("GridElements"):
            Internal._rmNode(f_pytree,elt_t)

    _recoverBoundaryConditions(t, f_pytree, zbcs, bctypes, bcnames)
    t_recbcs_end = time.perf_counter()
    t_recbcs2_start = time.perf_counter()
    print(" Cleaning frontIP (IBMWall) per processor...")
    f = Internal.getZones(f_pytree)
    if f != []:
        f = f[0]
        hook = C.createHook(f,"elementCenters")
        ids = C.identifyElements(hook, frontIP_gath, tol=TOL)
        ids = ids[ids[:] > -1]
        ids = ids.tolist()
        ids_IBMWall = [ids[i]-1 for i in range(len(ids))]
        C.freeHook(hook)
        if ids_IBMWall != []:
            frontIP = T.subzone(f,ids_IBMWall, type='elements')
            dimfrontIP = numpy.sum(Internal.getValue(frontIP)[0])
        else:
            frontIP = Internal.newZone(name="frontIP%d"%Cmpi.rank,zsize=[[0,0]],ztype="Unstructured")
            dimfrontIP = 0
    else:
        frontIP = Internal.newZone(name="frontIP%d"%Cmpi.rank,zsize=[[0,0]],ztype="Unstructured")
        dimfrontIP = 0

    t_recbcs2_end = time.perf_counter()
    if VPM == False:
        print(" Extracting front of the donor points...")
        t_frontdp_start = time.perf_counter()
        if frontTypeDP=="1":
            frontDP_gath = extractFrontDP(t, frontIP_gath, dim, sym3D, check, localDir=localDir)
        else:
            frontDP_gath = None

        del frontIP_gath

        t_frontdp_end = time.perf_counter()
        if dimfrontIP>0:
            if IBM_parameters["spatial discretization"]["type"] == "FV":
                t_normals_start = time.perf_counter()
                print(" Computing normals via project ortho..")
                _computeNormalsViaProjectOrtho(frontIP, tb2)
                t_normals_end = time.perf_counter()
                t_ipcenters_start = time.perf_counter()
                print(" Computing frontIP_C...")
                frontIP_C = C.node2Center(frontIP)
                Internal._rmNodesByType(frontIP_C,"Elements_t")
                t_ipcenters_end = time.perf_counter()
            elif IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
                frontIP_C = computeSurfaceQuadraturePoints(t, IBM_parameters, frontIP)
                frontIP_C = computeNormalsForDG(frontIP_C, tb2)
            t_ibpoints_start = time.perf_counter()
            print(" Extracting IBM Points...")
            ip_pts, image_pts, wall_pts = extractIBMPoints(tb2, frontIP, frontIP_C, frontDP_gath, bbo, IBM_parameters, check, dim, forceAlignment, localDir=localDir)
            t_ibpoints_end = time.perf_counter()
            t_dataset_start = time.perf_counter()
            print(" Adding IBCDatasets...")
            _addIBCDatasets(t, f, image_pts, wall_pts, ip_pts, IBM_parameters)
            t_dataset_end = time.perf_counter()
    else:
        if dimfrontIP>0:
            _addIBC2Zone(t, f, frontIP)

    C._rmVars(t,['cellNFront'])

    if IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
        _computeTurbulentDistanceForDG(t, tb2, IBM_parameters)
    for z in Internal.getZones(t):
        Cmpi._setProc(z,Cmpi.rank)

    if IBM_parameters["spatial discretization"]["type"] == "FV":
        if different_front_flag == False: #True is default
            Internal._rmNodesFromName(t,"TurbulentDistance")
            Internal._renameNode(t, "TurbulentDistanceForCFDComputation","TurbulentDistance")

    Internal._renameNode(t,'FlowSolution#Centers','FlisWallDistance')
    Cmpi.trace('AMR prepare IBM...end', master=True)
    return t

def computationDistancesNormals(t, tb, dim=3):
    if dim == 2:
        dz = 0.01
        tb2 = T.addkplane(tb)
        T._contract(tb2, (0,0,0), (1,0,0), (0,1,0), dz)
    else: tb2 = tb

    #if Cmpi.rank==0: C.convertPyTree2File(tb2,"tb2.plt")

    tc = C.node2Center(t)
    tb_WD = getBodiesForWallDistanceComputation(tb2)
    DTW._distance2Walls(t, tb_WD, type='ortho', signed=0, dim=3, loc='centers')
    X._applyBCOverlaps(t, depth=2, loc='centers', val=2, cellNName='cellN')
    C._initVars(t,'{centers:cellNChim}={centers:cellN}')
    Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dim, itype='chimera', order=2, cartesian=False)
    varsn=["gradxTurbulentDistance",'gradyTurbulentDistance','gradzTurbulentDistance']

    # A COMPARER !!
    if OPT: t = P.computeGrad(t, 'TurbulentDistance')
    else: P._computeGrad2(t, 'centers:TurbulentDistance', ghostCells=True, withCellN=False)

    for v in varsn: C._cpVars(t, 'centers:'+v, tc, v)
    C._cpVars(t,'centers:cellNChim',tc,'cellNChim')
    Xmpi._setInterpTransfers(t, tc, variables=varsn, cellNVariable='cellNChim', compact=0, type='ID')
    for v in varsn: t = C.center2Node(t,'centers:'+v)
    if not OPT: t = C.center2Node(t,'centers:TurbulentDistance')
    return t

def extractFrontDP(t, frontIP_gath, dim, sym3D, check, localDir='./'):
    if Cmpi.rank == 0:
        C._deleteEmptyZones(frontIP_gath)
        frontIP_gath = T.join(frontIP_gath)
        frontIP_gath = C.convertArray2Tetra(frontIP_gath)
        frontIP_gath = G.close(frontIP_gath)
        frontIP_gath[0] = "frontIP_gath"
        if dim ==3 and sym3D:
            print("Symmetry of frontIP gathered")
            # symmetry plane xz
            point = (0.0,0.0,0.0)
            vector1 = (1.0,0.0,0.0)
            vector2 = (0.0,0.0,1.0)
            z_sym = T.symetrize(frontIP_gath, point, vector1, vector2)
            z_sym[0] = "sym"
            frontIP_gath = Internal.getZones(frontIP_gath)[0]
            z_sym = Internal.getZones(z_sym)[0]
            frontIP_gath = T.join(frontIP_gath,z_sym)
            frontIP_gath = G.close(frontIP_gath)
    else:
        frontIP_gath = Internal.newZone(name="front",zsize=[[0,0]],ztype="Unstructured")
        gc = Internal.newGridCoordinates(parent=frontIP_gath)
        Internal.newDataArray('CoordinateX', value=numpy.empty(0), parent=gc)
        Internal.newDataArray('CoordinateY', value=numpy.empty(0), parent=gc)
        Internal.newDataArray('CoordinateZ', value=numpy.empty(0), parent=gc)
    frontIP_gath = Cmpi.bcastZone(frontIP_gath, root=0, coord=True)
    frontIP_gath = C.newPyTree(["body",frontIP_gath])
    C._initVars(t,'cellNFront',1.)
    X_IBM._blankByIBCBodies(t, frontIP_gath, 'nodes', 3,cellNName="cellNFront")
    del frontIP_gath
    frontDP = P.frontFaces(t,'cellNFront')
    del t
    frontDP_gath = Cmpi.allgatherZones(frontDP)
    frontDP_gath = T.join(frontDP_gath)
    frontDP_gath = C.newPyTree(["frontDP",frontDP_gath])
    if check:
        print("Exporting front DP...")
        if Cmpi.rank==0:
            C.convertPyTree2File(frontDP_gath,localDir+"frontDP_gath.cgns")
            C.convertPyTree2File(frontDP_gath,localDir+"frontDP_gath.plt")

    return frontDP_gath

def computeCellNForIBMFronts(t, dim, IBM_parameters, VPM=False):
    if VPM == False:
        frontTypeIP = IBM_parameters["integration points"]["front type"]
        print("frontTypeIP=",frontTypeIP)

        if frontTypeIP == "2":
            distance_IP = IBM_parameters["integration points"]["distance IntegrationPoints"]
            C._initVars(t,'distance_IP',distance_IP)
            C._initVars(t,'{cellN}=({TurbulentDistance}>{distance_IP})*{cellN}')
    else:
        snear = IBM_parameters["IBM type"]["size elements body"]
        C._initVars(t,'distance_IP',2*snear)
        C._initVars(t,'{cellN}=1-({TurbulentDistance}<-{distance_IP})')

    frontIP = P.frontFaces(t,'cellN')
    return frontIP

def extractIBMPoints(tb, frontIP, frontIP_C, frontDP, bbo, IBM_parameters, check, dim, forceAlignment=False, localDir='./'):

    frontTypeDP = IBM_parameters["donor points"]["front type"]
    frontTypeIP = IBM_parameters["integration points"]["front type"]
    IBMType = IBM_parameters["IBM type"]["type"]

    if frontTypeIP == "1":
        depth_IP = IBM_parameters["integration points"]["depth IntegrationPoints"]
    elif frontTypeIP == "2":
        distance_IP = IBM_parameters["integration points"]["distance IntegrationPoints"]

    if frontTypeDP == "2":
        distance_DP = IBM_parameters["donor points"]["distance DonorPoints"]
        C._initVars(frontIP_C,'dist',distance_DP)

    ip_pts = C.getAllFields(frontIP_C,loc='nodes', api=1)[0]
    ip_pts = Converter.convertArray2Node(ip_pts)
    ip_pts = [ip_pts]

    # Regrouping of the bodies per BC type
    bodies = []; listOfIBCTypes=[]

    for noz,zone in enumerate(Internal.getZones(tb)):
        body = C.getFields(Internal.__GridCoordinates__, zone, api=1)
        body = Converter.convertArray2Tetra(body)
        body = Transform.join(body)
        bodies.append(body)
        listOfIBCTypes.append("IBMWall%d" %noz)

    varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']

    if frontTypeDP == "2":
        res = connector.getIBMPtsWithoutFront(ip_pts, bodies, varsn, 'dist', 1)
        wallpts = res[0]
        imagepts = res[1]
    elif frontTypeDP == "1":
        frontDP = C.getFields(Internal.__GridCoordinates__, frontDP, api=1)
        frontDP = Converter.convertArray2Tetra(frontDP)
        listOfSnearsLoc=[]
        listOfModelingHeightsLoc = []
        snear = IBM_parameters["IBM type"]["size elements body"]
        if isinstance(snear,list): snear = min(snear)
        listOfSnearsLoc.append(snear)
        if frontTypeIP == "2": listOfModelingHeightsLoc.append(distance_IP)
        else: listOfModelingHeightsLoc.append(0.)
        res = connector.getIBMPtsWithFront(ip_pts, listOfSnearsLoc, listOfModelingHeightsLoc, bodies, frontDP, varsn, 1, 2, False, False)
        wallpts = res[0]
        imagepts = res[1]

        if frontTypeDP == "1": imagepts = moveIBMPoints(ip_pts, imagepts, wallpts, varsn, 1e-8)

    # Check if any of the image points lays outside the bbox. In this case we modify it.
    if isIBMPointInBox(bbo,imagepts)[0] == False:
        print("Warning: At least one donor point lays outside the bbox. The point is being moved closer to the wall...")
        list_ids_outside_box = isIBMPointInBox(bbo, imagepts)[1]
        epsilon = 0.9
        while (epsilon >= 0.1):
            print("Moving the badly located image point at epsilon %.2f %% of the initial distance from the wall." %epsilon)
            imagepts2correct = copy.deepcopy(imagepts)
            imagepts_modified = moveIBMPoints(ip_pts, imagepts2correct, wallpts, varsn, epsilon, list_ids_outside_box, tb)

            if isIBMPointInBox(bbo, imagepts_modified)[0] == False:
                epsilon = epsilon - 0.1
            else:
                imagepts = imagepts_modified
                break
        if isIBMPointInBox(bbo, imagepts)[0] == False:
            raise ValueError("Moving the points has not worked. Exiting..")

    wallpts = Converter.extractVars(wallpts,['CoordinateX','CoordinateY','CoordinateZ'])
    imagepts = Converter.extractVars(imagepts,['CoordinateX','CoordinateY','CoordinateZ'])
    ip_pts = Converter.extractVars(ip_pts,['CoordinateX','CoordinateY','CoordinateZ'])

    array_check = checkRelativeOffset(ip_pts[0][1], wallpts[0][1], imagepts[0][1], forceAlignment)
    if array_check.size != 0 and forceAlignment==True:
        wallpts = moveWallPoints(ip_pts, imagepts, wallpts, array_check, tb)
    checkCoincidentPoints(ip_pts[0][1], imagepts[0][1])

    if check:
        Converter.convertArrays2File(ip_pts  , localDir+"targetPts_proc%s.plt" %Cmpi.rank)
        Converter.convertArrays2File(wallpts , localDir+"wallPts_proc%s.plt" %Cmpi.rank)
        Converter.convertArrays2File(imagepts, localDir+"imagePts_proc%s.plt" %Cmpi.rank)

    dictOfImagePtsByIBCName={}
    dictOfTargetPtsByIBCName={}
    dictOfWallPtsByIBCName={}
    if (len(res)==3 and frontTypeDP=="2") or (len(res)==4 and frontTypeDP=="1"):
        noz = 0 #we always have only one zone
        indicesByTypeForZone = res[2][noz]
        nbTypes = len(indicesByTypeForZone)
        for nob in range(nbTypes):
            ibcTypeL = listOfIBCTypes[nob]
            indicesByTypeL = indicesByTypeForZone[nob]
            if indicesByTypeL.shape[0] > 0:
                ipPtsL = Transform.subzone(ip_pts[noz], indicesByTypeL)
                imagePtsL = Transform.subzone(imagepts[noz], indicesByTypeL)
                wallPtsL = Transform.subzone(wallpts[noz], indicesByTypeL)
            else:
                ipPtsL=[]; imagePtsL = []; wallPtsL = []

            dictOfTargetPtsByIBCName[ibcTypeL] = [ipPtsL]
            dictOfWallPtsByIBCName[ibcTypeL] = [wallPtsL]
            dictOfImagePtsByIBCName[ibcTypeL] = [imagePtsL]
    else:
        raise ValueError("The function connector.getIBMPtsWith/WithoutFront has not worked properly.")

    return dictOfTargetPtsByIBCName, dictOfImagePtsByIBCName,  dictOfWallPtsByIBCName

def isIBMPointInBox(bbox, coords):

    xmin = bbox[0]; ymin = bbox[1]; zmin = bbox[2]
    xmax = bbox[3]; ymax = bbox[4]; zmax = bbox[5]
    coords_x = coords[0][1][0,:]
    ids_x_min = numpy.where(coords_x>xmax)[0]
    ids_x_max = numpy.where(coords_x<xmin)[0]
    coords_y = coords[0][1][1,:]
    ids_y_min = numpy.where(coords_y>ymax)[0]
    ids_y_max = numpy.where(coords_y<ymin)[0]
    coords_z = coords[0][1][2,:]
    ids_z_min = numpy.where(coords_z>zmax)[0]
    ids_z_max = numpy.where(coords_z<zmin)[0]
    out = True
    list_ids_outside_box = numpy.concatenate([ids_x_min,ids_y_min,ids_z_min,ids_x_max,ids_y_max,ids_z_max])
    if list_ids_outside_box.shape[0] != 0: out = False

    return out, list_ids_outside_box

def moveWallPoints(ip_pts, imagepts, wallpts, array_check, tb):
    nb_image_pts = imagepts[0][1][0].size

    zsize = numpy.empty((1,3), E_NpyInt, order='F')
    zsize[0,0] = nb_image_pts; zsize[0,1] = 0; zsize[0,2] = 0
    zone_targetPts = Internal.newZone(name='TargetPoints',zsize=zsize,ztype='Unstructured')
    gc = Internal.newGridCoordinates(parent=zone_targetPts)
    Internal.newDataArray('CoordinateX', value=ip_pts[0][1][0], parent=gc)
    Internal.newDataArray('CoordinateY', value=ip_pts[0][1][1], parent=gc)
    Internal.newDataArray('CoordinateZ', value=ip_pts[0][1][2], parent=gc)

    DTW._distance2Walls(zone_targetPts, tb, type='ortho', signed=0, dim=3, loc='nodes')
    array_turb_dist = Internal.getNodeFromName(zone_targetPts,"TurbulentDistance")
    f_wall = open("wall_misaligned_after_proc%s.dat" %Cmpi.rank,"w")
    for count in array_check:
        dist = array_turb_dist[1][count]
        dirx0 = (imagepts[0][1][0][count]-ip_pts[0][1][0][count])
        diry0 = (imagepts[0][1][1][count]-ip_pts[0][1][1][count])
        dirz0 = (imagepts[0][1][2][count]-ip_pts[0][1][2][count])
        dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5
        dist0 = dist/dirn
        wallpts[0][1][0][count] = ip_pts[0][1][0][count] - dirx0*dist0
        wallpts[0][1][1][count] = ip_pts[0][1][1][count] - diry0*dist0
        wallpts[0][1][2][count] = ip_pts[0][1][2][count] - dirz0*dist0

        f_wall.write("%f %f %f\n" %(wallpts[0][1][0][count],wallpts[0][1][1][count],wallpts[0][1][2][count]))
    f_wall.close()
    return wallpts

def moveIBMPoints(ip_pts, imagepts, wallpts, varsn, epsilon, indices_outside_box=None, tb=None):

    nb_image_pts = imagepts[0][1][0].size
    if indices_outside_box is None:
        dist = epsilon
        for count in range(nb_image_pts):
            dirx0 = (imagepts[0][1][0][count]-wallpts[0][1][0][count])
            diry0 = (imagepts[0][1][1][count]-wallpts[0][1][1][count])
            dirz0 = (imagepts[0][1][2][count]-wallpts[0][1][2][count])
            dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5
            dist0 = dist/dirn
            imagepts[0][1][0][count] = imagepts[0][1][0][count] + dirx0*dist0
            imagepts[0][1][1][count] = imagepts[0][1][1][count] + diry0*dist0
            imagepts[0][1][2][count] = imagepts[0][1][2][count] + dirz0*dist0
    else:

        zsize = numpy.empty((1,3), E_NpyInt, order='F')
        zsize[0,0] = nb_image_pts; zsize[0,1] = 0; zsize[0,2] = 0
        zone_imagePts = Internal.newZone(name='ImagePoints',zsize=zsize,ztype='Unstructured')
        gc = Internal.newGridCoordinates(parent=zone_imagePts)
        Internal.newDataArray('CoordinateX', value=imagepts[0][1][0], parent=gc)
        Internal.newDataArray('CoordinateY', value=imagepts[0][1][1], parent=gc)
        Internal.newDataArray('CoordinateZ', value=imagepts[0][1][2], parent=gc)

        DTW._distance2Walls(zone_imagePts, tb, type='ortho', signed=0, dim=3, loc='nodes')
        array_turb_dist = Internal.getNodeFromName(zone_imagePts,"TurbulentDistance")
        for idx in indices_outside_box:
            dist = epsilon * array_turb_dist[1][idx]
            dist = epsilon * array_turb_dist[1][idx]
            dirx0 = (imagepts[0][1][0][idx]-wallpts[0][1][0][idx])
            diry0 = (imagepts[0][1][1][idx]-wallpts[0][1][1][idx])
            dirz0 = (imagepts[0][1][2][idx]-wallpts[0][1][2][idx])
            dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5
            dist0 = dist/dirn
            imagepts[0][1][0][idx] = wallpts[0][1][0][idx] + dirx0*dist0
            imagepts[0][1][1][idx] = wallpts[0][1][1][idx] + diry0*dist0
            imagepts[0][1][2][idx] = wallpts[0][1][2][idx] + dirz0*dist0

    return imagepts

def _recoverBoundaryConditions(t, f_pytree, zbcs, bctypes, bcnames):
    meshgen = "AMR"
    for z in Internal.getZones(t):
        if z is not None:
            nobc = len(zbcs)
            f = Internal.getZones(f_pytree)[0]
            if Cmpi.rank == 0: print("Performing the 'identifyElements' function (it can be long.)")
            for nobc, zbc in enumerate(zbcs):
                hook = C.createHook(f,"elementCenters")
                # Indices of the elements of f corresponding to the elements of zbc
                ids = C.identifyElements(hook, zbc, tol=TOL)
                len_ids = Internal.getValue(f)[0][1]
                ids = ids[ids[:] > -1] - 1
                ids = ids.tolist()
                if len(ids) > 0 :
                    zf = T.subzone(f,ids, type='elements')
                    if bcnames[nobc] != "QuadNQuad":
                        G_AMR._addBC2Zone__(z, bctypes[nobc], bctypes[nobc], zf)
                    else:
                        G_AMR._addBC2Zone__(z, "QuadNQuad", "FamilySpecified:QuadNQuad", zf)
                    ids_all = list(range(len_ids))
                    ids_new = list(set(ids_all)-set(ids))
                    if len(ids_new) > 0:
                        f = T.subzone(f,ids_new, type='elements')
                elif len(ids) == 0 and bcnames[nobc] == "QuadNQuad":
                    elts = Internal.getNodesFromType(z,"Elements_t")
                    maxElt = Internal.getNodeFromName(elts[-1],"ElementRange")[1][1]
                    CODABCType="QuadNQuad"
                    Internal.newElements(name=CODABCType, etype=7, econnectivity=numpy.empty(0),
                                         erange=[maxElt+1, maxElt], eboundary=1, parent=z)
                    C._addBC2Zone(z,CODABCType,"FamilySpecified:"+CODABCType, elementRange=[maxElt+1,maxElt])
                    zone_bc =  Internal.getNodeFromType(z,"ZoneBC_t")
                    lastbcname = C.getLastBCName(CODABCType)
                    node_bc = Internal.getNodeFromName(zone_bc,lastbcname)
                    node_bc[0] = CODABCType

                C.freeHook(hook)
            z[0] = z[0]+str(Cmpi.rank)
    if meshgen == "AMR": f_pytree[2][1][2] = [f]
    return None

def _addIBCDatasets(t, f, image_pts, wall_pts, ip_pts, IBM_parameters):

    if IBM_parameters["spatial discretization"]["type"]=="FV":
        N_IP_per_face = 1
    else:
        degree = IBM_parameters["spatial discretization"]["degree"]
        if IBM_parameters["spatial discretization"]["type"] == "DG":
            integrationDegree = 2*degree+1
            quadratureType = "GaussLegendre"
        elif IBM_parameters["spatial discretization"]["type"] == "DGSEM":
            integrationDegree = 2*degree-1
            quadratureType = "GaussLobatto"
        N_IP_per_face = GetReferencePointsQuad(integrationDegree, quadratureType)[0]
    list_suffix_datasets = [""]
    list_suffix_datasets.extend(range(1, N_IP_per_face))

    for z in Internal.getZones(t):

        hook = C.createHook(f, 'elementCenters')
        if IBM_parameters["spatial discretization"]["type"] == "FV":

            for nobc,ibc in enumerate(list(ip_pts.values())):
                if ibc!=[[]]:
                    coords_IBC_x = Converter.extractVars(ibc,["CoordinateX"])[0][1][0]
                    coords_IBC_y = Converter.extractVars(ibc,["CoordinateY"])[0][1][0]
                    coords_IBC_z = Converter.extractVars(ibc,["CoordinateZ"])[0][1][0]
                    zibc = Internal.newZone(name="IntegrationPoints",zsize=[[len(coords_IBC_x),0]],ztype="Unstructured")
                    gc = Internal.newGridCoordinates(parent=zibc)
                    Internal.newDataArray('CoordinateX', value=coords_IBC_x, parent=gc)
                    Internal.newDataArray('CoordinateY', value=coords_IBC_y, parent=gc)
                    Internal.newDataArray('CoordinateZ', value=coords_IBC_z, parent=gc)
                    ids = C.identifyNodes(hook, zibc, tol=TOL)
                    ids = ids[ids[:] > -1]
                    ids = ids.tolist()
                    ids = [ids[i]-1 for i in range(len(ids))]
                    zf = T.subzone(f,ids, type='elements')
                    G_AMR._addBC2Zone__(z, "IBMWall%d" %nobc, "FamilySpecified:IBMWall", zf)

        for bc in Internal.getNodesFromType(z,'BC_t'):
            famName = Internal.getNodeFromName(bc,'FamilyName')
            if famName is not None:
                if Internal.getValue(famName)=='IBMWall':
                    namebc = bc[0]
                    print(namebc)
                    ibcdataset=Internal.createNode('BCDataSet','BCDataSet_t',parent=bc,value='Null')
                    for i in range(N_IP_per_face):
                        dnrPts = Internal.createNode("DonorPointCoordinates"+str(list_suffix_datasets[i]),'BCData_t',parent=ibcdataset)
                        wallPts = Internal.createNode("WallPointCoordinates"+str(list_suffix_datasets[i]),'BCData_t',parent=ibcdataset)

                        coordsPD = Converter.extractVars(image_pts[namebc], ['CoordinateX','CoordinateY','CoordinateZ'])
                        coordsPW = Converter.extractVars(wall_pts[namebc], ['CoordinateX','CoordinateY','CoordinateZ'])

                        Internal.newDataArray('CoordinateX', value=coordsPD[0][1][0,:][i::N_IP_per_face], parent=dnrPts)
                        Internal.newDataArray('CoordinateY', value=coordsPD[0][1][1,:][i::N_IP_per_face], parent=dnrPts)
                        Internal.newDataArray('CoordinateZ', value=coordsPD[0][1][2,:][i::N_IP_per_face], parent=dnrPts)

                        Internal.newDataArray('CoordinateX', value=coordsPW[0][1][0,:][i::N_IP_per_face], parent=wallPts)
                        Internal.newDataArray('CoordinateY', value=coordsPW[0][1][1,:][i::N_IP_per_face], parent=wallPts)
                        Internal.newDataArray('CoordinateZ', value=coordsPW[0][1][2,:][i::N_IP_per_face], parent=wallPts)

                        if IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
                            coordsPI = Converter.extractVars(ip_pts[bc[0]], ['CoordinateX','CoordinateY','CoordinateZ'])
                            integrationPts = Internal.createNode("IntegrationPointCoordinates"+str(list_suffix_datasets[i]),'BCData_t',parent=ibcdataset)
                            Internal.newDataArray('CoordinateX', value=coordsPI[0][1][0,:][i::N_IP_per_face], parent=integrationPts)
                            Internal.newDataArray('CoordinateY', value=coordsPI[0][1][1,:][i::N_IP_per_face], parent=integrationPts)
                            Internal.newDataArray('CoordinateZ', value=coordsPI[0][1][2,:][i::N_IP_per_face], parent=integrationPts)

    return None

def getMinimumSpacing(t, dim, snear=1e-1):
    G._getVolumeMap(t)
    vol = Internal.getNodeFromName(t,"vol")[1]
    if dim==2: locsize = (vol/snear)**(1./2.)
    elif dim==3: locsize = (vol)**(1./3.)
    return min(locsize)

def computeDistance_IP_DP_front42_nonAdaptive(t, Reynolds, yplus_target, Lref, dim, snear=1e-2):
    distance_IP = D_IBM.computeModelisationHeight(Re=Reynolds, yplus=yplus_target, L=Lref)
    locsize = getMinimumSpacing(t, dim, snear)
    distance_DP = distance_IP+2*(dim**0.5)*locsize
    distance_DP = min(Cmpi.allgather(distance_DP))
    return distance_IP, distance_DP

def checkRelativeOffset(ip_pts, wallpts, imagepts, forceAlignment=False):

    x_wall = wallpts[0]
    y_wall = wallpts[1]
    z_wall = wallpts[2]

    x_image = imagepts[0]
    y_image = imagepts[1]
    z_image = imagepts[2]

    x_intp = ip_pts[0]
    y_intp = ip_pts[1]
    z_intp = ip_pts[2]

    wallPointFacePointDist = ( (x_wall - x_intp)**2  + (y_wall - y_intp)**2  + (z_wall - z_intp)**2 )**0.5
    donorPointWallPointDist =( (x_wall - x_image)**2 + (y_wall - y_image)**2 + (z_wall - z_image)**2 )**0.5

    x_facePointWallPointVec = (x_intp - x_wall) / wallPointFacePointDist;
    y_facePointWallPointVec = (y_intp - y_wall) / wallPointFacePointDist;
    z_facePointWallPointVec = (z_intp - z_wall) / wallPointFacePointDist;

    x_offsetPointLocation = x_wall + x_facePointWallPointVec * donorPointWallPointDist;
    y_offsetPointLocation = y_wall + y_facePointWallPointVec * donorPointWallPointDist;
    z_offsetPointLocation = z_wall + z_facePointWallPointVec * donorPointWallPointDist;

    epsDonorPointOffset = 1e-6
    offsetcheck = ( (x_offsetPointLocation - x_image)**2 + (y_offsetPointLocation - y_image)**2 + (z_offsetPointLocation - z_image)**2 )**0.5
    array_check =  numpy.where(offsetcheck>(epsDonorPointOffset * donorPointWallPointDist))[0]
    print("Check not aligned points: size array=",array_check.size)
    if array_check.size!=0:
        if forceAlignment:
            print("ATTENTION!!!!!!! Max offset on rank %d = " %Cmpi.rank,numpy.max(offsetcheck))
            f_wall = open("wall_misaligned_before_proc%s.dat" %Cmpi.rank,"w")
            f_target = open("target_misaligned_before_proc%s.dat" %Cmpi.rank,"w")
            f_image = open("image_misaligned_before_proc%s.dat" %Cmpi.rank,"w")
            for i in range(array_check.size):
                f_wall.write("%f %f %f\n" %(x_wall[array_check[i]],y_wall[array_check[i]],z_wall[array_check[i]]))
                f_target.write("%f %f %f\n" %(x_intp[array_check[i]],y_intp[array_check[i]],z_intp[array_check[i]]))
                f_image.write("%f %f %f\n" %(x_image[array_check[i]],y_image[array_check[i]],z_image[array_check[i]]))
            f_wall.close()
            f_target.close()
            f_image.close()
        else:
            raise ValueError("The maximum allowed relative tangential offset was exceeded by one of the donor points.\nMax offset on rank %d = " %Cmpi.rank,numpy.max(offsetcheck))
    return array_check

def checkCoincidentPoints(ip_pts, imagepts):

    x_image = imagepts[0]
    y_image = imagepts[1]
    z_image = imagepts[2]

    x_intp = ip_pts[0]
    y_intp = ip_pts[1]
    z_intp = ip_pts[2]

    donorPointIntegrationPointDist =( (x_intp - x_image)**2 + (y_intp - y_image)**2 + (z_intp - z_image)**2 )**0.5
    epsDonorPointOffset = 1e-9
    array_check =  numpy.where(donorPointIntegrationPointDist < epsDonorPointOffset)[0]
    print("Check coincident points: size array=",array_check.size)
    if array_check.size!=0:
        raise ValueError("ATTENTION!!!!!!! coincident integration and donor points rank %d = " %Cmpi.rank,numpy.min(donorPointIntegrationPointDist))
    return

def getBodiesForWallDistanceComputation(tb2):
    zones_tb = Internal.getZones(tb2)
    zones_tb_WD = []
    for z_tb in zones_tb:
        ibctype = Internal.getNodeFromName(z_tb, "ibctype")
        if ibctype == None:
            ibctype = 0
        else:
            ibctype = ibctype[1][0]
        if ibctype == 0:
            zones_tb_WD.append(z_tb)
        #print("ibctype=",ibctype)
    tb_WD = C.newPyTree(["tbWD", zones_tb_WD])
    return tb_WD

def _computeNormalsViaProjectOrtho(front, tb2):

    varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
    front_centers = C.node2Center(front)
    proj = T.projectOrtho(front_centers, tb2); proj[0] = 'projection'
    x_proj = Internal.getNodeFromName(proj,"CoordinateX")[1]
    y_proj = Internal.getNodeFromName(proj,"CoordinateY")[1]
    z_proj = Internal.getNodeFromName(proj,"CoordinateZ")[1]
    x_front = Internal.getNodeFromName(front_centers,"CoordinateX")[1]
    y_front = Internal.getNodeFromName(front_centers,"CoordinateY")[1]
    z_front = Internal.getNodeFromName(front_centers,"CoordinateZ")[1]
    dirx0 = (x_front-x_proj)
    diry0 = (y_front-y_proj)
    dirz0 = (z_front-z_proj)
    dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5
    dirx0 = dirx0/dirn
    diry0 = diry0/dirn
    dirz0 = dirz0/dirn
    zone = Internal.getZones(front)
    FS = Internal.newFlowSolution(name='FlowSolution#Centers', gridLocation='CellCenter', parent=zone[0])
    Internal.newDataArray(varsn[0], value=dirx0, parent=FS)
    Internal.newDataArray(varsn[1], value=diry0, parent=FS)
    Internal.newDataArray(varsn[2], value=dirz0, parent=FS)
    return None

def _addIBC2Zone(t, f, frontIP):
    for z in Internal.getZones(t):
        hook = C.createHook(f, 'elementCenters')
        ids = C.identifyElements(hook, frontIP, tol=TOL)
        ids = ids[ids[:] > -1]
        ids = ids.tolist()
        ids = [ids[i]-1 for i in range(len(ids))]
        zf = T.subzone(f,ids, type='elements')
        G_AMR._addBC2Zone(z, "IBMWall", "FamilySpecified:IBMWall", zf)
    return None

def computeSurfaceQuadraturePoints(t, IBM_parameters, frontIP):
    zones = Internal.getZones(t)
    f = P.exteriorFaces(zones[0])
    dims_f = Internal.getZoneDim(f)
    #print(dims_f)
    for elt_t in Internal.getNodesFromType(f, "Elements_t"):
        if not elt_t[0].startswith("GridElements"):
            Internal._rmNode(f,elt_t)

    hook = C.createHook(f, 'elementCenters')
    ids = C.identifyElements(hook, frontIP, tol=TOL)
    ids = ids[ids[:] > -1]
    ids = ids.tolist()
    ids_IBMWall = [ids[i]-1 for i in range(len(ids))]
    zf = T.subzone(f,ids_IBMWall, type='elements')
    G_AMR._addBC2Zone__(zones[0], 'IBMWall0', 'FamilySpecified:IBMWall',zf)

    ## Computation of the surface quadrature points
    degree = IBM_parameters["spatial discretization"]["degree"]
    if IBM_parameters["spatial discretization"]["type"] == "DG":
        integrationDegree = 2*degree+1
        quadratureType = "GaussLegendre"
    elif IBM_parameters["spatial discretization"]["type"] == "DGSEM":
        integrationDegree = 2*degree-1
        quadratureType = "GaussLobatto"
    else:
        raise ValueError("Unkown discretization type; options are: FV, DG or DGSEM.")

    coordsX = Internal.getNodeFromName(t,"CoordinateX")[1]
    coordsY = Internal.getNodeFromName(t,"CoordinateY")[1]
    coordsZ = Internal.getNodeFromName(t,"CoordinateZ")[1]
    elts = Internal.getNodesFromType(t,"Elements_t")
    IBMWall_node = Internal.getNodeFromName(elts,"IBMWall0")
    N_IBM_cells = IBMWall_node[1][1]
    IBMWall_EC = (Internal.getNodeFromName(IBMWall_node,"ElementConnectivity")[1]).reshape(N_IBM_cells,4)
    cellType = 4
    N_IP_per_face = GetReferencePointsQuad(integrationDegree, quadratureType)[0]
    N_IP = N_IP_per_face * N_IBM_cells
    tic = time.perf_counter()
    weights, interpolationMatrix = GetReferencePointsData(integrationDegree, quadratureType, cellType)
    quadPoints_surf_location = numpy.empty((N_IP,3))
    for i,j in zip(range(0,N_IP,N_IP_per_face), range(N_IBM_cells)):
        connectivity_IBMFace = IBMWall_EC[j]-1
        pos_first = numpy.argmin(connectivity_IBMFace)
        pos_second_array_min = numpy.array([connectivity_IBMFace[(pos_first+1)%4], connectivity_IBMFace[pos_first-1]])
        value_second = numpy.min(pos_second_array_min)
        pos_second = numpy.argwhere(connectivity_IBMFace==value_second)[0][0]
        index_right_left = numpy.argmin(pos_second_array_min)
        if index_right_left == 0:
            pos_third = (pos_second+1)%4
            pos_fourth = (pos_third+1)%4
        elif index_right_left ==1:
            pos_third = (pos_second-1)
            pos_fourth = (pos_third-1)
        connectivity_IBMFace[[0,1,2,3]] = connectivity_IBMFace[[pos_first, pos_second, pos_third, pos_fourth]]
        nodalData = numpy.hstack([coordsX[connectivity_IBMFace].reshape(4,1),coordsY[connectivity_IBMFace].reshape(4,1),coordsZ[connectivity_IBMFace].reshape(4,1)])
        quadPoints_surf_location[i:i+N_IP_per_face,:] = interpolationMatrix.dot(nodalData)

    toc = time.perf_counter()

    print("Time to compute the surface quadrature points : ", toc-tic)

    N_IP_surf = N_IBM_cells * N_IP_per_face

    z_IP = Internal.newZone(name="SurfaceIntegrationPoints",zsize=[[N_IP_surf,0]],ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=z_IP)
    Internal.newDataArray('CoordinateX', value=quadPoints_surf_location[:,0], parent=gc)
    Internal.newDataArray('CoordinateY', value=quadPoints_surf_location[:,1], parent=gc)
    Internal.newDataArray('CoordinateZ', value=quadPoints_surf_location[:,2], parent=gc)
    return z_IP

def computeNormalsForDG(z_IP, tb):
    wall_Points = T.projectOrtho(z_IP,tb)
    x_wallPoints = Internal.getNodeFromName(wall_Points,"CoordinateX")[1]
    y_wallPoints = Internal.getNodeFromName(wall_Points,"CoordinateY")[1]
    z_wallPoints = Internal.getNodeFromName(wall_Points,"CoordinateZ")[1]

    x_integrationPoints = Internal.getNodeFromName(z_IP,"CoordinateX")[1]
    y_integrationPoints = Internal.getNodeFromName(z_IP,"CoordinateY")[1]
    z_integrationPoints = Internal.getNodeFromName(z_IP,"CoordinateZ")[1]

    dirx0 = (x_wallPoints-x_integrationPoints)
    diry0 = (y_wallPoints-y_integrationPoints)
    dirz0 = (z_wallPoints-z_integrationPoints)

    dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5

    dirx0 = dirx0/dirn
    diry0 = diry0/dirn
    dirz0 = dirz0/dirn
    varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
    #varsn = ["nx","ny","nz"]
    FS = Internal.newFlowSolution(name='FlowSolution', gridLocation='Vertex', parent=z_IP)
    Internal.newDataArray(varsn[0], value=dirx0, parent=FS)
    Internal.newDataArray(varsn[1], value=diry0, parent=FS)
    Internal.newDataArray(varsn[2], value=dirz0, parent=FS)
    return z_IP

def _computeTurbulentDistanceForDG(t, tb, IBM_parameters):

    degree = IBM_parameters["spatial discretization"]["degree"]
    if IBM_parameters["spatial discretization"]["type"] == "DG":
        integrationDegree = 2*degree+1
        quadratureType = "GaussLegendre"
    elif IBM_parameters["spatial discretization"]["type"] == "DGSEM":
        integrationDegree = 2*degree-1
        quadratureType = "GaussLobatto"
    else:
        raise ValueError("Unkown discretization type; options are: FV, DG or DGSEM.")
    coordsX = Internal.getNodeFromName(t,"CoordinateX")[1]
    coordsY = Internal.getNodeFromName(t,"CoordinateY")[1]
    coordsZ = Internal.getNodeFromName(t,"CoordinateZ")[1]
    elts = Internal.getNodesFromType(t,"Elements_t")
    GE_node = Internal.getNodeFromName(elts,"GridElements")
    GE_EC_ravel = Internal.getNodeFromName(GE_node,"ElementConnectivity")[1]
    N_volume_cells = len(GE_EC_ravel)//8
    GE_EC = GE_EC_ravel.reshape(N_volume_cells,8)

    cellType = 8
    N_IP_per_cell = GetReferencePointsHexa(integrationDegree, quadratureType)[0]
    tic = time.perf_counter()
    weights, interpolationMatrix = GetReferencePointsData(integrationDegree, quadratureType, cellType)
    N_IP = N_IP_per_cell * N_volume_cells
    quadPoints_vol_location = numpy.empty((N_IP,3))
    for i,j in zip(range(0,N_IP,N_IP_per_cell), range(N_volume_cells)):
        nodalData = numpy.hstack([coordsX[GE_EC[j]-1].reshape(8,1),coordsY[GE_EC[j]-1].reshape(8,1),coordsZ[GE_EC[j]-1].reshape(8,1)])
        quadPoints_vol_location[i:i+N_IP_per_cell,:] = interpolationMatrix.dot(nodalData)

    toc = time.perf_counter()
    print("Time to compute the position of the volume quadrature points : ", toc-tic)

    N_IP_vol = N_volume_cells * N_IP_per_cell
    z_IP = Internal.newZone(name="VolumeIntegrationPoints",zsize=[[N_IP_vol,0]],ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=z_IP)
    Internal.newDataArray('CoordinateX', value=quadPoints_vol_location[:,0], parent=gc)
    Internal.newDataArray('CoordinateY', value=quadPoints_vol_location[:,1], parent=gc)
    Internal.newDataArray('CoordinateZ', value=quadPoints_vol_location[:,2], parent=gc)
    DTW._distance2Walls(z_IP, tb, type='ortho', signed=0, loc='nodes')

    if N_IP_vol%N_volume_cells != 0: raise ValueError("The division between the number of the integration points divided by the number of IBM faces is not exact. Every face should have the same number of integration points.")

    list_suffix_datasets = [""]
    list_suffix_datasets.extend(range(1, N_IP_per_cell))
    walldistance_volume_ip = Internal.getNodeFromName(z_IP, "TurbulentDistance")[1]

    zones = Internal.getZones(t) #always one single zone
    for i in range(N_IP_per_cell):
        walldistance_dataset=Internal.newFlowSolution('FlisWallDistance'+str(list_suffix_datasets[i]),parent=zones[0],gridLocation="CellCenter")
        Internal.newDataArray('TurbulentDistance',value=walldistance_volume_ip[i::N_IP_per_cell],parent=walldistance_dataset)

    return None

## IMPORTANT NOTE:: this is a template of a python wrapper. Not to be used. It is a placeholder and is very likely to change in subsequent version.
def prepareAMRIBM(tb, levelMax, vmins, dim, IBM_parameters, toffset=None, check=False, opt=False, octreeMode=1,
                  snears=0.01, dfars=10, loadBalancing=False, conformal=False, OutputAMRMesh=False,
                  localDir='./', fileName=None):
    """Generate AMR IBM mesh and prepare AMR IBM data for CODA simulation. 
    Usage: prepareAMRIBM(tb, levelMax, vmins, dim, IBM_parameters, toffset, check, opt, octreeMode,
                  snears, dfars, loadBalancing, conformal, OutputAMRMesh,
                  localDir= fileName)"""

    t_AMR = G_AMR.generateAMRMesh(tb=tb, levelMax=levelMax, vmins=vmins, dim=dim,
                                  toffset=None, check=check, opt=opt, octreeMode=octreeMode, localDir=localDir,
                                  snears=0.01, dfars=10, loadBalancing=False)

    Cmpi.trace('AMR Mesh Dist2Walls...start', master=True)
    if dim == 2: tb2 = T.addkplane(tb)
    DTW._distance2Walls(t_AMR, tb2, type='ortho', signed=0, dim=dim, loc='centers')
    DTW._distance2Walls(t_AMR, tb2, type='ortho', signed=0, dim=dim, loc='nodes')
    del tb2
    Cmpi.trace('AMR Mesh Dist2Walls...end', master=True)

    if OutputAMRMesh: Cmpi.convertPyTree2File(t_AMR, localDir+'tAMRMesh.cgns')
    t_AMR = prepareAMRData(tb, t_AMR, IBM_parameters=IBM_parameters, dim=dim, check=check, localDir=localDir)

    if fileName is not None:
        Cmpi.convertPyTree2File(t_AMR, localDir+fileName)
        return None
    else:
        return t_AMR
