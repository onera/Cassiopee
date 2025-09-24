#Cassiopee imports
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Post.AMR as P_AMR
import os, sys, numpy, math

##IMPORTANT NOTE:: The original Post_IBM_CODA.py is copied in the same directory as this file for refernces & tracking purposes.

#FSDM imports
try:
    from FSDataManager import FSClac, FSError, FSMesh, FSMeshEnums, FSDataName, FSFloatArray, FSStringArray, FSDataSpecArray, FSDatasetInfo, FSQuantityDescArray
    from CODA import DiscretizationFactory, TimeIntegrationFactory, FSMeshInterpolationVolumeReconstruction
except:
    raise ImportError("interpolationDonorPoints: FSDataManager and CODA")


def interpolationDonorPoints(fsmesh, clac, discParaDict, discSelectionParaDict, marker_IBM, wall_boundary_markers=[], localDir='./', check=False):
    
    try:
        from FSZoltan import FSZoltanInterface  # noqa: 401
        WITH_FSZOLTAN = True
    except ImportError:
        WITH_FSZOLTAN = False
    ##Parallel function 

    def _checkRepartitionOption(fsmesh, clac):
        from FSDataManager import IsCompiledWith
        isPARMETIS = IsCompiledWith("PARMETIS")
        fsmesh.RepartitionMeshRCB() or FSError.PrintAndExit()
        if WITH_FSZOLTAN:
            fsmesh.RepartitionMeshZOLTAN(GraphExtraction={"GraphType": 0}) or FSError.PrintAndExit()
        elif isPARMETIS:
            fsmesh.RepartitionMeshPARMETIS(GraphExtraction={"GraphType": 0}) or FSError.PrintAndExit()
        elif clac.GetNProcs() > 1:
            raise NotImplementedError("Your FlowSimulator installation does not have a graph partitioner...FSZoltan & PARMETIS are missing...")
        return None

    
    ## Interpolation of the flow field solution onto the image points....
    ## CODA :: image point --> donor point ... ghost cells --> integration point/face
    _checkRepartitionOption(fsmesh, clac)
    
    if Cmpi.rank==0: print("Interpolation of the flow solution on the donor points...", flush=True)
    
    wall        = fsmesh.GetUnstructDataset("WallPointCoordinates").GetValues()
    donor       = fsmesh.GetUnstructDataset("DonorPointCoordinates").GetValues()
    wall_numpy  = numpy.array(wall.Buffer() , copy=True)
    donor_numpy = numpy.array(donor.Buffer(), copy=True)
    nb_vertex   = wall.Size(0)
    
    #Get marker list
    fs_boundary_marker_list = fsmesh.GetCellAttributeValuesWithNames("CADGroupID")
    np_boundary_marker_list = numpy.array(fs_boundary_marker_list.Buffer(), copy=True)
    
    fs_boundary_markers_celltype = fsmesh.GetCellAttribute("CADGroupID",4)
    np_boundary_markers_celltype = numpy.array(fs_boundary_markers_celltype.Buffer(), copy=True)
    
    indices      = numpy.ravel(numpy.argwhere(np_boundary_markers_celltype==marker_IBM[0]))
    nb_nodes_IBM = len(indices)
    if wall_boundary_markers != []:
      wall_numpy_IBM  = numpy.concatenate([wall_numpy[indices], wallBF_numpy[indices_wall]])
      donor_numpy_IBM = numpy.concatenate([donor_numpy[indices],wallBF_numpy[indices_wall]])
    else:
      wall_numpy_IBM  = wall_numpy[indices]
      donor_numpy_IBM = donor_numpy[indices]
    
    fsmesh_dp = FSMesh(clac)
    fsmesh_dp.BeginInitialization()
    fsmesh_dp.InitUnstructNodes(nb_nodes_IBM)
    
    fs_coordinates_donor = FSFloatArray(donor_numpy_IBM.shape[0],3)
    fs_coordinates_wall  = FSFloatArray(wall_numpy_IBM.shape[0] ,3)
    
    numpy.copyto(numpy.array(fs_coordinates_donor.Buffer(), copy=False), donor_numpy_IBM, casting='unsafe')
    numpy.copyto(numpy.array(fs_coordinates_wall.Buffer() , copy=False), wall_numpy_IBM , casting='unsafe')
    
    coordNames    = FSStringArray(3)
    coordNames[0] = FSDataName.Coordinate().X()
    coordNames[1] = FSDataName.Coordinate().Y()
    coordNames[2] = FSDataName.Coordinate().Z()
    
    coordSpecs = FSDataSpecArray(3)
    coordSpecs[0].Length()
    coordSpecs[1].Length()
    coordSpecs[2].Length()
    
    fsmesh_dp.InitUnstructDataset(FSDataName("DonorPoints"),
                                  FSDatasetInfo(coordNames, coordSpecs, FSMeshEnums.CT_Node),
                                  fs_coordinates_donor)
    
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

    _checkRepartitionOption(fsmesh_dp,clac)
    
    # import mesh of donor points
    fsmesh_dp.HasLocalNumbering() or fsmesh_dp.CreateLocalNumbering()
    
    # prepare for interpolation of quantities from fsmesh to fsmesh_dp
    srcCoords = FSQuantityDescArray()
    fsmesh.DetermineCoordinates(srcCoords) or FSError.PrintAndExit()
    quantities = FSQuantityDescArray()
    dstCoords  = FSQuantityDescArray()
    fsmesh_dp.DetermineCoordinates(dstCoords) or FSError.PrintAndExit()
    ipvr = FSMeshInterpolationVolumeReconstruction("sourcePoint")
    ipvr.InterpolateAtNodes(fsmesh, srcCoords, disc, state, fsmesh_dp, dstCoords, quantities, False) or FSError.PrintAndExit()
    
    if check:
        fsmesh_dp.ExportMeshTECPLOT(Filename=localDir+'pointMesh_AugState.dat', FileFormat='ASCII',PrefixDatasetName=True) or FSError.PrintAndExit()
    node_donorpoints = fsmesh_dp.GetUnstructDataset("DonorPoints").GetValues()
    node_wallpoints  = fsmesh_dp.GetUnstructDataset("WallPoints").GetValues()
    augState_data    = fsmesh_dp.GetUnstructDataset("CODAAugState").GetValues()
    augState_names   = fsmesh_dp.GetUnstructDataset("CODAAugState").GetNames()
    
    node_donorpoints = Cmpi.allgather(numpy.array(node_donorpoints.Buffer(), copy=True))
    node_wallpoints  = Cmpi.allgather(numpy.array(node_wallpoints.Buffer() , copy=True))
    augState_data    = Cmpi.allgather(numpy.array(augState_data.Buffer()   , copy=True))

    
    return node_donorpoints, node_wallpoints, augState_data


def computeSurfValues(fileNameResultIn, tb, CODAInputs, dim=3, fileNameIBMPnts=None, fileNameResultOut=None, localDir='./', check=False, verbose=False):
    ## CODAInputs = [discParaDictAllStages, discSelectionParaDict, IBM_markers, dictReferenceQuantities]
    ## fileNameResultIn = (e.g) output_stage3.h5 (the h5 from the CODA run)
    if isinstance(tb, str): tb = C.convertFile2PyTree(tb)
    else: tb = Internal.copyTree(tb)

    discParaDict            = CODAInputs[0]
    discSelectionParaDict   = CODAInputs[1]
    IBM_markers             = CODAInputs[2]
    dictReferenceQuantities = CODAInputs[3]
    
    clac   = FSClac()
    fsmesh = FSMesh(clac)
    fsmesh.ImportMeshHDF5(Filename=localDir+fileNameResultIn) or FSError.PrintAndExit()
    
    # Reconstruction of the solution at the donor points (in parallel).
    donorPoints_numpy, wallPoints_numpy, augState_numpy = interpolationDonorPoints(fsmesh, clac, discParaDict, discSelectionParaDict, IBM_markers, localDir=localDir, check=check)
    
    # Serial part of the post-processing.
    zw = None
    if Cmpi.master:
        pytree = P_AMR.createPyTreeForIBMWallFieldsExtraction(numpy.concatenate(donorPoints_numpy),
                                                        numpy.concatenate(wallPoints_numpy),
                                                        numpy.concatenate(augState_numpy),
                                                        discSelectionParaDict)
        if fileNameIBMPnts is not None: C.convertPyTree2File(pytree,localDir+fileNameIBMPnts)
        zw = P_AMR.extractIBMWallFields(pytree,tb,discSelectionParaDict)
    
    zw       = Cmpi.bcast(zw,root=0)    
    zw,coefs = P_AMR.computeBoundaryQuantities(zw, dictReferenceQuantities, dim=dim, verbose=verbose)
    
    if Cmpi.master:
        print("\nIntegrated coefficients:")
        print("CD=      %g"%coefs[0])
        print("CDfric=  %g"%coefs[2])
        print("CDpres=  %g"%coefs[3])
        print("CL=      %g"%coefs[1])
        print("CLfric=  %g"%coefs[4])
        print("CLpres=  %g"%coefs[5])
        
    if fileNameResultOut is not None:
        if Cmpi.master: C.convertPyTree2File(zw,localDir+"out_wallNEW.plt")
        return None
    else:
        return zw, coefs

##========================================================================
##========================================================================

def test_MPI_Part1(clac):
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
    indices_testMPI_numpy     = numpy.array(indices_testMPI.Buffer(), copy=True)
    indices_testMPI_gathered  = Cmpi.allgather(indices_testMPI_numpy)
    return indices_testMPI_gathered


def test_MPI_Part2(indices_testMPI_numpy, nProcs):
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

  return None

