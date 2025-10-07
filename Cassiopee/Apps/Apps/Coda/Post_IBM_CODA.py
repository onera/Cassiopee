import os,sys
#Cassiopee imports
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Connector.PyTree as X
import Intersector.PyTree as XOR
import Connector.IBM as X_IBM
import Connector.connector as connector
import Converter.Internal as Internal
import Post.IBM as P_IBM
import Post.PyTree as P
import Generator.PyTree as G
import Transform.PyTree as T

import numpy
from mpi4py import MPI
#FSDM imports
import FSDM
from FSDataManager import FSClac, FSLog, FSError, FSMesh, FSMeshEnums, FSDataName, FSString, FS_AT_CADGroupID, FSFloatArray, FSStringArray, FSDataSpecArray, FSDataSpec, FSDatasetInfo, FS_AT_GlobalNumber, FSIntArray
from FSDataManager import FSTimer, FSDataLog, FSUnstructVolumeCellTypes, FSUnstructSurfaceCellTypes, FSMeshEnums, FSMeshSelectionOpAttribute, FSMeshSelectionOpCellType, FSMeshSelectionOpDataset, FSQuantityDescArray
from FSDataManager import IsCompiledWith
#FSDM modules import
try:
    from FSZoltan import FSZoltanInterface  # noqa: 401
    WITH_FSZOLTAN = True
except ImportError:
    WITH_FSZOLTAN = False

def interpolationDonorPoints(fsmesh,clac,discParaDict,discSelectionParaDict, marker_IBM, wall_boundary_markers=[]):
    # import from CODA
    from CODA import DiscretizationFactory, TimeIntegrationFactory, FSMeshInterpolationVolumeReconstruction

    fsmesh.RepartitionMeshRCB() or FSError.PrintAndExit()

    if WITH_FSZOLTAN:
        fsmesh.RepartitionMeshZOLTAN(GraphExtraction={"GraphType": 0}) or FSError.PrintAndExit()
    elif IsCompiledWith("PARMETIS"):
        fsmesh.RepartitionMeshPARMETIS(GraphExtraction={"GraphType": 0}) or FSError.PrintAndExit()
    elif clac.GetNProcs() > 1:
        raise NotImplementedError("Your FlowSimulator installation does not have a graph partitioner.")

    if Cmpi.rank==0: print("Interpolation of the flow solution on the donor points...")

    wall = fsmesh.GetUnstructDataset("WallPointCoordinates").GetValues()

    wall_numpy = numpy.array(wall.Buffer(), copy=True)
    nb_vertex = wall.Size(0)
    donor = fsmesh.GetUnstructDataset("DonorPointCoordinates").GetValues()
    donor_numpy = numpy.array(donor.Buffer(), copy=True)
    #Get marker list
    fs_boundary_marker_list = fsmesh.GetCellAttributeValuesWithNames("CADGroupID")
    np_boundary_marker_list = numpy.array(fs_boundary_marker_list.Buffer(), copy=True)

    fs_boundary_markers_celltype = fsmesh.GetCellAttribute("CADGroupID",4)
    np_boundary_markers_celltype = numpy.array(fs_boundary_markers_celltype.Buffer(), copy=True)
    indices = numpy.ravel(numpy.argwhere(np_boundary_markers_celltype==marker_IBM[0]))
    nb_nodes_IBM = len(indices)
    if wall_boundary_markers!=[]:
        wall_numpy_IBM = numpy.concatenate([wall_numpy[indices], wallBF_numpy[indices_wall]])
        donor_numpy_IBM = numpy.concatenate([donor_numpy[indices],wallBF_numpy[indices_wall]])
    else:
        wall_numpy_IBM = wall_numpy[indices]
        donor_numpy_IBM = donor_numpy[indices]

    fsmesh_dp = FSMesh(clac)
    fsmesh_dp.BeginInitialization()

    fsmesh_dp.InitUnstructNodes(nb_nodes_IBM)

    fs_coordinates = FSFloatArray(donor_numpy_IBM.shape[0],3)
    fs_coordinates_wall = FSFloatArray(wall_numpy_IBM.shape[0],3)

    numpy.copyto(numpy.array(fs_coordinates.Buffer(), copy=False), donor_numpy_IBM, casting='unsafe')
    numpy.copyto(numpy.array(fs_coordinates_wall.Buffer(), copy=False), wall_numpy_IBM, casting='unsafe')

    coordNames = FSStringArray(3)
    coordNames[0] = FSDataName.Coordinate().X()
    coordNames[1] = FSDataName.Coordinate().Y()
    coordNames[2] = FSDataName.Coordinate().Z()

    coordSpecs = FSDataSpecArray(3)
    coordSpecs[0].Length()
    coordSpecs[1].Length()
    coordSpecs[2].Length()

    fsmesh_dp.InitUnstructDataset(FSDataName.Coordinates(),
                                  FSDatasetInfo(coordNames, coordSpecs, FSMeshEnums.CT_Node),
                                  fs_coordinates)

    fsmesh_dp.InitUnstructDataset(FSDataName("WallPoints"),
                                  FSDatasetInfo(coordNames, coordSpecs, FSMeshEnums.CT_Node),
                                  fs_coordinates_wall)

    fsmesh_dp.EndInitialization()

    # prepare face-based data for preprocessing
    fsmesh.HasLocalNumbering() or fsmesh.CreateLocalNumbering()


    # create discretization (which includes any preprocessing that may be necessary)
    disc = DiscretizationFactory.GetSingleton().Create(discSelectionParaDict, clac, fsmesh, discParaDict)

    state = disc.CreateZeroFieldVector()
    state.ImportFromFSMesh(disc.GetMeshInterface(), fsmesh, "State") or FSError.PrintAndExit()


    # solution transfer and output
    state.ExportToFSMesh(disc.GetMeshInterface(), fsmesh, "State") or FSError.PrintAndExit()


    fsmesh_dp.RepartitionMeshRCB() or FSError.PrintAndExit()

    if WITH_FSZOLTAN:
        fsmesh_dp.RepartitionMeshZOLTAN(GraphExtraction={"GraphType": 0}) or FSError.PrintAndExit()
    elif IsCompiledWith("PARMETIS"):
        fsmesh_dp.RepartitionMeshPARMETIS(GraphExtraction={"GraphType": 0}) or FSError.PrintAndExit()
    elif clac.GetNProcs() > 1:
        raise NotImplementedError("Your FlowSimulator installation does not have a graph partitioner.")

    # import mesh of donor points
    fsmesh_dp.HasLocalNumbering() or fsmesh_dp.CreateLocalNumbering()


    # prepare for interpolation of quantities from fsmesh to fsmesh_dp
    srcCoords = FSQuantityDescArray()
    fsmesh.DetermineCoordinates(srcCoords) or FSError.PrintAndExit()
    quantities = FSQuantityDescArray()
    dstCoords = FSQuantityDescArray()
    fsmesh_dp.DetermineCoordinates(dstCoords) or FSError.PrintAndExit()
    ipvr = FSMeshInterpolationVolumeReconstruction("sourcePoint")
    ipvr.InterpolateAtNodes(fsmesh, srcCoords, disc, state, fsmesh_dp, dstCoords, quantities, False) or FSError.PrintAndExit()

    fsmesh_dp.ExportMeshTECPLOT(Filename='pointMesh_AugState.dat', FileFormat='ASCII',PrefixDatasetName=True) or FSError.PrintAndExit()
    #fsmesh_dp.ExportMeshHDF5(Filename='pointMesh_AugState.h5') or FSError.PrintAndExit()
    node_coordinates = fsmesh_dp.GetUnstructDataset("Coordinates").GetValues()
    node_wallpoints = fsmesh_dp.GetUnstructDataset("WallPoints").GetValues()
    augState_data = fsmesh_dp.GetUnstructDataset("CODAAugState").GetValues()
    augState_names = fsmesh_dp.GetUnstructDataset("CODAAugState").GetNames()

    augState_names_list = []

    #################### TEST MPI ########################### (temporary)
    nProcs = clac.GetNProcs()
    if clac.GetProcID()<10:
        indices_testMPI = FSFloatArray(2)
        indices_testMPI[0] = clac.GetProcID()
        indices_testMPI[1] = clac.GetProcID()+10
    else:
        indices_testMPI = FSFloatArray(2)
        indices_testMPI.Fill(0.)
    #########################################################

    node_coordinates_numpy = numpy.array(node_coordinates.Buffer(), copy=True)
    node_coordinates_gathered = Cmpi.allgather(node_coordinates_numpy)
    node_wallpoints_numpy = numpy.array(node_wallpoints.Buffer(), copy=True)
    node_wallpoints_gathered = Cmpi.allgather(node_wallpoints_numpy)
    augState_data_numpy = numpy.array(augState_data.Buffer(), copy=True)
    augState_data_gathered = Cmpi.allgather(augState_data_numpy)
    indices_testMPI_numpy = numpy.array(indices_testMPI.Buffer(), copy=True)
    indices_testMPI_gathered = Cmpi.allgather(indices_testMPI_numpy)

    return indices_testMPI_gathered, node_coordinates_gathered,node_wallpoints_gathered,augState_data_gathered,augState_names_list

def createPyTreeForIBMWallFieldsExtraction(donorPoints_numpy,wallPoints_numpy,augState_numpy,discSelectionParaDict):#,augState_names):

    print("Creating pyTree of the quantities from the FSMesh of the donorPoints..")
    nb_node =  donorPoints_numpy.shape[0]
    utau_vector = numpy.zeros(nb_node)
    yplus_vector = numpy.zeros(nb_node)

    zsize = numpy.empty((1,3), numpy.int32, order='F')
    zsize[0,0] = nb_node; zsize[0,1] = 0; zsize[0,2] = 0
    z = Internal.newZone(name='IBMWall',zsize=zsize,ztype='Unstructured')
    pytree = C.newPyTree(["DonorPointsPyTree",z])

    coordx = numpy.ravel(donorPoints_numpy[:,0])
    coordy = numpy.ravel(donorPoints_numpy[:,1])
    coordz = numpy.ravel(donorPoints_numpy[:,2])

    idx_sorted = numpy.lexsort((coordx,coordy,coordz))

    #Create coordinate node
    gc = Internal.newGridCoordinates(parent=z)
    Internal.newDataArray('CoordinateX', value=coordx[idx_sorted], parent=gc)
    Internal.newDataArray('CoordinateY', value=coordy[idx_sorted], parent=gc)
    Internal.newDataArray('CoordinateZ', value=coordz[idx_sorted], parent=gc)
    #Create flow solution node
    FS = Internal.newFlowSolution(name='FlowSolution', gridLocation='CellCenter', parent=z)
    Internal.newDataArray("Density", value=augState_numpy[:,0][idx_sorted], parent=FS)
    if discSelectionParaDict["PDE"]=="RANSSAneg":
        Internal.newDataArray("Pressure", value=augState_numpy[:,6][idx_sorted], parent=FS)
        Internal.newDataArray("ViscosityMolecular", value=augState_numpy[:,9][idx_sorted], parent=FS)
        Internal.newDataArray("Temperature", value=augState_numpy[:,7][idx_sorted], parent=FS)
    elif discSelectionParaDict["PDE"]=="NavierStokes":
        Internal.newDataArray("Pressure", value=augState_numpy[:,5][idx_sorted], parent=FS)
        Internal.newDataArray("ViscosityMolecular", value=augState_numpy[:,8][idx_sorted], parent=FS)
        Internal.newDataArray("Temperature", value=augState_numpy[:,6][idx_sorted], parent=FS)
    elif discSelectionParaDict["PDE"]=="Euler":
        Internal.newDataArray("Pressure", value=augState_numpy[:,5][idx_sorted], parent=FS)
        Internal.newDataArray("ViscosityMolecular", value=0*augState_numpy[:,5][idx_sorted], parent=FS)
        Internal.newDataArray("Temperature", value=augState_numpy[:,6][idx_sorted], parent=FS)
    velX = augState_numpy[:,1][idx_sorted]/(augState_numpy[:,0][idx_sorted])
    Internal.newDataArray("VelocityX", value=velX, parent=FS)
    velY = augState_numpy[:,2][idx_sorted]/(augState_numpy[:,0][idx_sorted])
    Internal.newDataArray("VelocityY", value=velY, parent=FS)
    velZ = augState_numpy[:,3][idx_sorted]/(augState_numpy[:,0][idx_sorted])
    Internal.newDataArray("VelocityZ", value=velZ, parent=FS)
    Internal.newDataArray('CoordinateX_PW', value=numpy.ravel(wallPoints_numpy[:,0])[idx_sorted], parent=FS)
    Internal.newDataArray('CoordinateY_PW', value=numpy.ravel(wallPoints_numpy[:,1])[idx_sorted], parent=FS)
    Internal.newDataArray('CoordinateZ_PW', value=numpy.ravel(wallPoints_numpy[:,2])[idx_sorted], parent=FS)
    Internal.newDataArray('CoordinateX_PI', value=numpy.ravel(donorPoints_numpy[:,0])[idx_sorted], parent=FS)
    Internal.newDataArray('CoordinateY_PI', value=numpy.ravel(donorPoints_numpy[:,1])[idx_sorted], parent=FS)
    Internal.newDataArray('CoordinateZ_PI', value=numpy.ravel(donorPoints_numpy[:,2])[idx_sorted], parent=FS)
    Internal.newDataArray('utau', value=utau_vector, parent=FS)
    Internal.newDataArray('yplus', value=yplus_vector, parent=FS)
    return pytree

def extractIBMWallFields(pytree,tb,discSelectionParaDict,tempViscFlag=False):
    from numpy.linalg import norm

    print("Extracting Wall Fields..")
    z = Internal.getZones(pytree)[0]
    zname = Internal.getName(z)

    X_IBM._computeFrictionVelocity(z)

    ibctype = 3 # 3 = Musker
    nameSubRegion = 'IBCD_%d_%s'%(ibctype,zname)
    FS = Internal.getNodeFromType(z,"FlowSolution_t")

    v = numpy.fromstring(zname, 'c')
    z[2].append([nameSubRegion, v, [],'ZoneSubRegion_t'])
    info = z[2][len(z[2])-1]# les enfants de la subregion
    info[2].append(FS[2][1:])
    info[2] = info[2][0]
    #This is a workaround for the function P_IBM.extractIBMWallFields, which projects the variable "Temperature" and not the "ViscosityMolecular". -> We put then the viscosity into a node called "Temperature".
    if tempViscFlag:
        Internal.getNodeFromName(info,"Temperature")[1] = Internal.getNodeFromName(info,"ViscosityMolecular")[1]
    zw = P_IBM.extractIBMWallFields(z, tb=tb, famZones=[])#, front=1)
    if discSelectionParaDict["PDE"]=="Euler":
        C._initVars(zw,"utau",0)
        C._initVars(zw,"yplus",0)
    return zw

def computeBoundaryQuantities(zw,dictReferenceQuantities,dimPb=3,reorderFlag=False,invertYZ=False):
    if Cmpi.rank==0:
        print("Computing integral coefficients..")
    if dimPb==2:
        zw = C.convertBAR2Struct(zw)
        T._addkplane(zw)

    if reorderFlag==True:
        T._reorder(zw,(-1,))

    if invertYZ == False:
        Sref=dictReferenceQuantities["Sref"]
        alpha=dictReferenceQuantities["alpha"]
        beta=dictReferenceQuantities["beta"]
    else:
        Sref=dictReferenceQuantities["Sref"]
        alpha=dictReferenceQuantities["beta"]
        beta=-dictReferenceQuantities["alpha"]

    zw = C.convertArray2Tetra(zw); zw = G.close(zw)
    zw = C.node2Center(zw,Internal.__FlowSolutionNodes__)
    # add reference state for computation of integrated coefficients
    ref1 = Internal.getNodesFromName(zw,"ReferenceState")
    if ref1!=[]:
        Internal._rmNodesByName(zw,"ReferenceState")
    ref = Internal.newReferenceState("ReferenceState",parent=zw)
    ref[2].append(["VelocityX",dictReferenceQuantities["velX_ref"], [],'DataArray_t'])
    ref[2].append(["VelocityY",dictReferenceQuantities["velY_ref"], [],'DataArray_t'])
    ref[2].append(["VelocityZ",dictReferenceQuantities["velZ_ref"], [],'DataArray_t'])
    ref[2].append(["Density",dictReferenceQuantities["density_ref"], [],'DataArray_t'])
    ref[2].append(["Pressure",dictReferenceQuantities["pressure_ref"], [],'DataArray_t'])
    [res, res2, [clp, cdp, clatp], [clf, cdf, clatf]] = _loads0(zw,Sref=Sref,alpha=alpha,beta=beta,dimPb=dimPb,verbose=True)
    nProcs = Cmpi.size
    cdf = cdf/nProcs
    cdp = cdp/nProcs
    clf = clf/nProcs
    clp = clp/nProcs
    clatf = clatf/nProcs
    clatp = clatp/nProcs
    cd = (cdf+cdp)
    cl = (clf+clp)
    clat = (clatf+clatp)
    return zw, [cd, cl, clat, cdf, cdp, clf, clp, clatf, clatp]


def test_MPI(indices_testMPI_numpy,nProcs):
    ########## TEST MPI ###########   (temporary)
    #print("output concatenate",numpy.concatenate(indices_testMPI_numpy))
    #print(nProcs)
    if nProcs>9:
        test_array_indices = [i for i in range(0,10)]
    else:
        test_array_indices = [i for i in range(0,nProcs)]
    full_test_array_indices = numpy.zeros(nProcs*2)
    counter = 0
    for i in test_array_indices:
        full_test_array_indices[counter] = i
        full_test_array_indices[counter+1] = i+10
        counter = counter+2
    if numpy.equal(numpy.concatenate(indices_testMPI_numpy),full_test_array_indices).all()==False:
        print("Error: full test array indices", full_test_array_indices)
        sys.exit("Problem with MPI. Gather not in the correct order")
    ##############################

    return

def computeExtraVariablesForLocalIBM(ts,dictReferenceQuantities,dimPb=3):

    import math
    import Post.ExtraVariables2 as PE

    print("Computing extra variables..")
    # add reference state for computation of integrated coefficients
    Pref=None; Qref=None;
    if dimPb==2:
        ts = C.convertBAR2Struct(ts)
        T._addkplane(ts)

    ts = C.convertArray2Tetra(ts); ts = G.close(ts)
    ts = C.node2Center(ts,Internal.__FlowSolutionNodes__)

    zones = Internal.getZones(ts)
    zone_IBM = zones[-1] #last zone is IBM
    if zones:
        ref1 = Internal.getNodesFromName(ts,"ReferenceState")
        if ref1!=[]:
            Internal._rmNodesByName(ts,"ReferenceState")
        ref = Internal.newReferenceState("ReferenceState",parent=ts)
        ref[2].append(["VelocityX",dictReferenceQuantities["velX_ref"], [],'DataArray_t'])
        ref[2].append(["VelocityY",dictReferenceQuantities["velY_ref"], [],'DataArray_t'])
        ref[2].append(["VelocityZ",dictReferenceQuantities["velZ_ref"], [],'DataArray_t'])
        ref[2].append(["Density",dictReferenceQuantities["density_ref"], [],'DataArray_t'])
        ref[2].append(["Pressure",dictReferenceQuantities["pressure_ref"], [],'DataArray_t'])
        alpha=dictReferenceQuantities["alpha"];beta=dictReferenceQuantities["beta"]

        RefState = Internal.getNodeFromType(ts,'ReferenceState_t')
        PInf     = Internal.getValue(Internal.getNodeFromName(RefState,"Pressure"))
        RoInf    = Internal.getValue(Internal.getNodeFromName(RefState,"Density"))
        VxInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityX"))
        VyInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityY"))
        VzInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityZ"))
        VInf2    = VxInf*VxInf+VyInf*VyInf+VzInf*VzInf
        VInf     = math.sqrt(VInf2)
        q      = 0.5*RoInf*VInf2
        qinv   = 1./q

        Qref = q
        Pref = PInf
        #only zone IBM!!!


        P_IBM._computeExtraVariables(zone_IBM, Pref, Qref, variables=['Cp','Cf','frictionX','frictionY','frictionZ','frictionMagnitude','ShearStress'])

    return ts

def divideLocalIBMZones(wall_IBM, wall_BF):

    wall_BF = C.convertArray2Tetra(wall_BF)

    CP_BF = Internal.getNodeFromName(wall_BF,"CoefPressure")[1]
    CFX_BF = Internal.getNodeFromName(wall_BF,"CoefSkinFrictionX")[1]
    CFY_BF = Internal.getNodeFromName(wall_BF,"CoefSkinFrictionY")[1]
    CFZ_BF = Internal.getNodeFromName(wall_BF,"CoefSkinFrictionZ")[1]
    CFT_BF = Internal.getNodeFromName(wall_BF,"CoefSkinFrictionTangential")[1]

    z_IBM = Internal.getZones(wall_IBM)[0]

    hook = C.createHook(wall_IBM, 'nodes')
    ids = C.identifyNodes(hook, wall_BF)
    ids_clean = ids[ids[:] > -1]
    ids_clean = ids_clean.tolist()
    zf_clean = T.subzone(wall_IBM,ids_clean, type='nodes')
    hook = C.createHook(wall_IBM, 'elementCenters')
    ids = C.identifyElements(hook, zf_clean)
    ids_clean = ids[ids[:] > -1]-1
    ids_clean = ids_clean.tolist()
    zf_clean_2 = T.subzone(wall_IBM,ids_clean, type='elements')


    list_all_ids_IBM = [i for i in range(z_IBM[1][0][1])]
    ids_ibm_zone = list(set(list_all_ids_IBM)-set(ids_clean))
    zf_ibm_zone = T.subzone(wall_IBM,ids_ibm_zone, type='elements')

    # paste solution on boundaries (Cp/Cf) with new indices
    hook = C.createHook(wall_BF, 'nodes')
    ids_clean_zf2 = C.identifyNodes(hook, zf_clean_2)
    ids_clean_zf2 = ids_clean_zf2[ids_clean_zf2[:] > -1]-1
    ids_clean_zf2 = ids_clean_zf2.tolist()
    FS_BF = Internal.getNodeFromName(zf_clean_2,"FlowSolution")
    FS_BF[2].append(["Cp",CP_BF[ids_clean_zf2],[],'DataArray_t'])
    FS_BF[2].append(["frictionX",CFX_BF[ids_clean_zf2],[],'DataArray_t'])
    FS_BF[2].append(["frictionY",CFY_BF[ids_clean_zf2],[],'DataArray_t'])
    FS_BF[2].append(["frictionZ",CFZ_BF[ids_clean_zf2],[],'DataArray_t'])
    FS_BF[2].append(["Cf",CFT_BF[ids_clean_zf2],[],'DataArray_t'])

    t = C.newPyTree(["New PyTree", [Internal.getZones(zf_clean_2)[0],Internal.getZones(zf_ibm_zone)[0]]])
    return t


def integCfn(teff):
    """Integ tau.n.ds"""
    import Post.Mpi as Pmpi
    retx = Pmpi.integ(teff, 'centers:frictionX')
    rety = Pmpi.integ(teff, 'centers:frictionY')
    retz = Pmpi.integ(teff, 'centers:frictionZ')
    return [retx[0],rety[0],retz[0]]


def _loads0LocalIBM(ts,Sref=1.0,Qref=1.0,alpha=0.0,beta=0.0,dimPb=2,BFtreatment="nodal"):
    dimPb = 2
    import Transform.PyTree as T
    import Generator.PyTree as G
    import math
    import Post.ExtraVariables2 as PE

    T._reorder(ts,(-1,))

    z_BF = Internal.getZones(ts)[0]
    z_IBM = Internal.getZones(ts)[1]


    res_pres_BF = PE.integCp(z_BF)[0]
    res_pres_BF    = [-i/Sref for i in res_pres_BF]
    res_fric_BF = integCfn(z_BF)
    res_fric_BF   = [ i/Sref for i in res_fric_BF]

    res_pres_IBM = PE.integCp(z_IBM)[0]
    res_pres_IBM    = [-i/Sref for i in res_pres_IBM]
    res_fric_IBM = PE.integTaun(z_IBM)
    res_fric_IBM   = [ i/Sref for i in res_fric_IBM]

    alpha  = math.radians(alpha)
    beta   = math.radians(beta)
    calpha = math.cos(alpha); cbeta = math.cos(beta)
    salpha = math.sin(alpha); sbeta = math.sin(beta)

    def computeClCd(res_p,res_f):
        cdp = res_p[0]*calpha*cbeta + res_p[1]*salpha*cbeta - res_p[2]*sbeta
        clp = res_p[1]*calpha       - res_p[0]*salpha
        cdf = res_f[0]*calpha*cbeta + res_f[1]*salpha*cbeta - res_f[2]*sbeta
        clf = res_f[1]*calpha       - res_f[0]*salpha
        return [cdp,clp,cdf,clf]

    cdp_BF,  clp_BF,  cdf_BF,  clf_BF  = computeClCd(res_pres_BF, res_fric_BF)
    cdp_IBM, clp_IBM, cdf_IBM, clf_IBM = computeClCd(res_pres_IBM, res_fric_IBM)


    print(cdp_BF,  clp_BF,  cdf_BF,  clf_BF)
    cd_BF = cdp_BF+cdf_BF
    cl_BF = clp_BF+clf_BF

    cd_IBM = cdp_IBM+cdf_IBM
    cl_IBM = clp_IBM+clf_IBM

    cdf = cdf_BF + cdf_IBM
    cdp = cdp_BF + cdp_IBM
    cd  = cd_BF + cd_IBM

    clf = clf_BF + clf_IBM
    clp = clp_BF + clp_IBM
    cl  = cl_BF + cl_IBM

    FSC = Internal.getNodesFromName(ts,Internal.__FlowSolutionCenters__)
    Internal._rmNodesFromName(FSC, 'ShearStress*')

    return ts, [[cd, cl, cdf, cdp, clf, clp], [cd_BF, cl_BF, cdf_BF, cdp_BF, clf_BF, clp_BF], [cd_IBM, cl_IBM, cdf_IBM, cdp_IBM, clf_IBM, clp_IBM]]

def intersectObstacleMesh(t, surf_detail,list_bcs):
    print("Performing intersection for Local IBM approach..")

    zones = Internal.getZones(t)
    z = zones[0]

    zoneDim = Internal.getZoneDim(z)
    #print(zoneDim[0])
    if zoneDim[0] == "Unstructured":
        nb_vertex = zoneDim[1]
        nb_cells = zoneDim[2]
    elif zoneDim[0] == "Structured":
        index_i, index_j, index_k = zoneDim[1], zoneDim[2], zoneDim[3]
        nb_vertex = index_i * index_j * index_k
        nb_cells = (index_i-1) * (index_j-1) * (index_k-1)

    zbcs=[]; bctypes=[]; bcs=[];
    for bc in Internal.getNodesFromType(z,'BC_t'):
        bctype = Internal.getValue(bc)

        if bctype not in bctypes:
            bctypes.append(bctype)
            bcs.append(bc)

    for i,bctype in enumerate(bctypes):
        zbc = C.extractBCOfType(t,bctype)
        zbc = T.join(zbc)
        zbcs.append(zbc)
        if bctype.startswith("BCWall"):
            nobc = i
            wall_clean = zbc

    #4: intersect surfaces
    surf_detail = T.join(surf_detail)
    surf_detail = C.convertArray2Tetra(surf_detail); surf_detail = G.close(surf_detail)

    wall_clean = Internal.getZones(wall_clean)[0]
    wall_clean = C.convertArray2Tetra(wall_clean); wall_clean = G.close(wall_clean)

    x = XOR.conformUnstr(surf_detail, wall_clean, tol=0., left_or_right=2, itermax=1)
    x = T.splitManifold(x)

    xzs = Internal.getZones(x)
    result = T.join(xzs[0], xzs[2]) # warning : cas-dependant !

    return result

def extractNotIBMWallFields(t,IBM_parameters,list_boundary_values_to_extract,discSelectionParaDict,BFtreatment="nodal"):
    print("Extracting wall fields that are not IBM..")
    C._deleteFlowSolutions__(t)

    zones = Internal.getZones(t)
    z = zones[0]

    zoneDim = Internal.getZoneDim(z)
    print(zoneDim[0])

    for i,bc in enumerate(Internal.getNodesFromType(z,'BC_t')):
        bctype = Internal.getValue(bc)
        if bctype.startswith("BCWall"):
            zbc = C.extractBCOfType(t,bctype)
            break
    zbc = T.join(zbc)
    FS = Internal.getNodesFromName(zbc,"FlowSolution#Centers")
    Internal._rmNodesFromName(FS,"CoefSkinFrictionImmersed")

    names_var = list_boundary_values_to_extract
    FS_1 = []
    for i in range(len(names_var)):
        FS_ = numpy.empty((0,3))
        for j in range(len(FS)):
            FS_ = numpy.append(FS_,FS[j][2][i+1][1])
        FS_1.append(FS_)

    # only one
    pytree_new = C.newPyTree(["PyTree"])
    indicesF=[]; f = P.exteriorFaces(z, indices=indicesF)
    hook = C.createHook(f, 'elementCenters')

    zbc = T.join(zbc)
    ids = C.identifyElements(hook, zbc)
    ids = ids[ids[:] > -1]
    ids = ids.tolist()
    ids = [ids[i]-1 for i in range(len(ids))]
    zf = T.subzone(f,ids, type='elements')

    list_skin_friction_quantities = ["CoefSkinFrictionX","CoefSkinFrictionY","CoefSkinFrictionZ","CoefSkinFrictionTangential"]
    if discSelectionParaDict["PDE"]=="Euler": names_var = names_var+list_skin_friction_quantities
    zf[2] = zf[2][:3]
    GE = Internal.getNodeFromName(zf,"GridElements_HEXA")
    for idx,name_var in enumerate(names_var):
        zf = C.initVars(zf, 'centers:%s' %name_var,0.)
        if discSelectionParaDict["PDE"]!="Euler" or idx==0:
            node_PCFC = Internal.getNodeFromName(zf,name_var)
            node_PCFC[1] = FS_1[idx]

    zf_nodes = C.center2Node(zf,"FlowSolution#Centers")
    Internal._rmNodesByName(zf_nodes,"FlowSolution#Centers")

    return zf_nodes

def _loads0(ts, Sref=None, Pref=None, Qref=None, alpha=0., beta=0., dimPb=3, verbose=False,time=0):

    import Post.ExtraVariables2 as PE
    import math

    zones = Internal.getZones(ts)
    if zones:
        if Sref is None:
            C._initVars(ts, '__ONE__',1.)
            Sref = P.integ(ts, '__ONE__')[0];
            C._rmVars(ts, ['__ONE__', 'centers:vol'])

        RefState = Internal.getNodeFromType(ts,'ReferenceState_t')
        PInf     = Internal.getValue(Internal.getNodeFromName(RefState,"Pressure"))
        RoInf    = Internal.getValue(Internal.getNodeFromName(RefState,"Density"))
        VxInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityX"))
        VyInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityY"))
        VzInf    = Internal.getValue(Internal.getNodeFromName(RefState,"VelocityZ"))
        VInf2    = VxInf*VxInf+VyInf*VyInf+VzInf*VzInf
        VInf     = math.sqrt(VInf2)
        q      = 0.5*RoInf*VInf2

        qinv   = 1./q
        alpha  = math.radians(alpha)
        beta   = math.radians(beta)
        calpha = math.cos(alpha); cbeta = math.cos(beta)
        salpha = math.sin(alpha); sbeta = math.sin(beta)

        if Qref is None: Qref = q
        if Pref is None: Pref = PInf

        P_IBM._computeExtraVariables(ts, Pref, Qref, variables=['Cp','Cf','frictionX','frictionY','frictionZ','frictionMagnitude','ShearStress'])

        #===========================
        # Compute pressure forces
        #===========================
    res = PE.integCp(ts)[0]
    res2 = PE.integTaun(ts)
    clp=0
    cdp=0
    clf=0
    cdf=0
    if zones:
        res    = [-i/Sref for i in res]
        QADIMI = 1./(Qref*Sref)
        res2   = [ i*QADIMI for i in res2]

        calpha = math.cos(alpha); cbeta = math.cos(beta)
        salpha = math.sin(alpha); sbeta = math.sin(beta)
        if dimPb == 3:
            cdp =   res[0]*calpha*cbeta + res[1]*salpha*cbeta - res[2]*sbeta
            clp =   res[1]*calpha       - res[0]*salpha
            clatp = res[0]*calpha*sbeta + res[1]*salpha*sbeta + res[2]*cbeta

            cdf =   res2[0]*calpha*cbeta + res2[1]*salpha*cbeta - res2[2]*sbeta
            clf =   res2[1]*calpha       - res2[0]*salpha
            clatf = res2[0]*calpha*sbeta + res2[1]*salpha*sbeta + res2[2]*cbeta

        else:
            cdp = res[0]*calpha + res[1]*salpha
            clp = res[1]*calpha - res[0]*salpha
            clatp = 0.0
            cdf = res2[0]*calpha + res2[1]*salpha
            clf = res2[1]*calpha - res2[0]*salpha
            clatf = 0.0

        if verbose and Cmpi.rank==0:
            print("Normalized pressure drag = %.4e, lift = %.4e and lateral force = %.4e"%(cdp, clp, clatp))
            print("Vector of pressure loads: (Fx_P,Fy_P,Fz_P)=(%.4e, %.4e, %.4e)"%(res[0],res[1],res[2]))


            if time==0:time=-1
            print("****************************************")
            print("Total Drag (time/nit=%.4e) : %.4e"%(time,(cdp+cdf)))
            print("Total Lift (time/nit=%.4e) : %.4e"%(time,(clp+clf)))
            print("Total Lateral Force (time/nit=%.4e) : %.4e"%(time,(clatp+clatf)))
            print("****************************************")

        FSC = Internal.getNodesFromName(ts,Internal.__FlowSolutionCenters__)
        Internal._rmNodesFromName(FSC, 'ShearStress*')

    return [res, res2, [clp, cdp, clatp], [clf, cdf, clatf]]
