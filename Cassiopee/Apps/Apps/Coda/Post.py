#Cassiopee imports
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Post.AMR as P_AMR
import os, sys, numpy, math

##IMPORTANT NOTE:: The original Post_IBM_CODA.py is copied in the same directory as this file for refernces & tracking purposes.

#FSDM imports
try:
    from FSDataManager import FSClac, FSError, FSMesh, FSMeshEnums, FSDataName, \
        FSFloatArray, FSStringArray, FSDataSpecArray, FSDatasetInfo, FSQuantityDescArray
    from CODA import DiscretizationFactory, TimeIntegrationFactory, FSMeshInterpolationVolumeReconstruction
except:
    raise ImportError("interpolationDonorPoints: FSDataManager and CODA")

def outputTime(startTime,functionName='FunctionName'):
    endTime     = time.perf_counter()
    elapsedTime = endTime-startTime
    elapsedTime = Cmpi.allreduce(elapsedTime  ,op=Cmpi.MAX)
    if Cmpi.rank==0: print('Elapsed Time: %s: %g [s] | %g [min] | %g [hr]'%(functionName,elapsedTime,elapsedTime/60,elapsedTime/3600),flush=True)
    return None


def interpolationDonorPoints(fsmesh, clac, discParaDict, discSelectionParaDict, markerIBM, wallBoundaryMarkers=[], check=False):
    """ Interpolate the CODA flow field to the IBM donor points (image points in ghost cell IBM approach). 
    Usage: interpolationDonorPoints(fsmesh, clac, discParaDict, discSelectionParaDict, markerIBM, wallBoundaryMarkers, check)"""
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

    Cmpi.trace('[INTERPOLATIONDONORPOINTS] Interpolation of the flow solution on the donor points..start', master=True)

    wall        = fsmesh.GetUnstructDataset("WallPointCoordinates").GetValues()
    donor       = fsmesh.GetUnstructDataset("DonorPointCoordinates").GetValues()
    wallNumpy  = numpy.array(wall.Buffer() , copy=True)
    donorNumpy = numpy.array(donor.Buffer(), copy=True)

    #Get marker list
    fsBoundaryMarkerList = fsmesh.GetCellAttributeValuesWithNames("CADGroupID")
    npBoundaryMarkerList = numpy.array(fsBoundaryMarkerList.Buffer(), copy=True)

    fsBoundaryMarkersCelltype = fsmesh.GetCellAttribute("CADGroupID",4)
    npBoundaryMarkersCelltype = numpy.array(fsBoundaryMarkersCelltype.Buffer(), copy=True)

    #indices      = numpy.ravel(numpy.argwhere(npBoundaryMarkersCelltype == markerIBM[0]))
    indices = numpy.ravel(numpy.argwhere(numpy.isin(npBoundaryMarkersCelltype, markerIBM)))
    nb_nodes_IBM = len(indices)
    if wallBoundaryMarkers != []:
        ##This makes no sense...i have no indices_wall...taken directory from original toolbox
        wallNumpyIBM  = numpy.concatenate([wallNumpy[indices],  wallBF_numpy[indices_wall]])
        donorNumpyIBM = numpy.concatenate([donorNumpy[indices], wallBF_numpy[indices_wall]])
    else:
        wallNumpyIBM  = wallNumpy[indices]
        donorNumpyIBM = donorNumpy[indices]

    fsmeshDonorPnt = FSMesh(clac)
    fsmeshDonorPnt.BeginInitialization()
    fsmeshDonorPnt.InitUnstructNodes(nb_nodes_IBM)

    fsCoordinatesDonor = FSFloatArray(donorNumpyIBM.shape[0], 3)
    fsCoordinatesWall  = FSFloatArray(wallNumpyIBM.shape[0] , 3)

    numpy.copyto(numpy.array(fsCoordinatesDonor.Buffer(), copy=False), donorNumpyIBM, casting='unsafe')
    numpy.copyto(numpy.array(fsCoordinatesWall.Buffer() , copy=False), wallNumpyIBM , casting='unsafe')

    coordNames    = FSStringArray(3)
    coordNames[0] = FSDataName.Coordinate().X()
    coordNames[1] = FSDataName.Coordinate().Y()
    coordNames[2] = FSDataName.Coordinate().Z()

    coordSpecs = FSDataSpecArray(3)
    coordSpecs[0].Length()
    coordSpecs[1].Length()
    coordSpecs[2].Length()

    fsmeshDonorPnt.InitUnstructDataset(FSDataName("DonorPoints"),
                                       FSDatasetInfo(coordNames, coordSpecs, FSMeshEnums.CT_Node),
                                       fsCoordinatesDonor)

    fsmeshDonorPnt.InitUnstructDataset(FSDataName("WallPoints"),
                                       FSDatasetInfo(coordNames, coordSpecs, FSMeshEnums.CT_Node),
                                       fsCoordinatesWall)

    fsmeshDonorPnt.EndInitialization()

    # prepare face-based data for preprocessing
    fsmesh.HasLocalNumbering() or fsmesh.CreateLocalNumbering()

    # create discretization (which includes any preprocessing that may be necessary)
    disc = DiscretizationFactory.GetSingleton().Create(discSelectionParaDict, clac, fsmesh, discParaDict)

    state = disc.CreateZeroFieldVector()
    state.ImportFromFSMesh(disc.GetMeshInterface(), fsmesh, "State") or FSError.PrintAndExit()

    # solution transfer and output
    state.ExportToFSMesh(disc.GetMeshInterface(), fsmesh, "State") or FSError.PrintAndExit()

    _checkRepartitionOption(fsmeshDonorPnt,clac)

    # import mesh of donor points
    fsmeshDonorPnt.HasLocalNumbering() or fsmeshDonorPnt.CreateLocalNumbering()

    # prepare for interpolation of quantities from fsmesh to fsmeshDonorPnt
    srcCoords = FSQuantityDescArray()
    fsmesh.DetermineCoordinates(srcCoords) or FSError.PrintAndExit()
    quantities = FSQuantityDescArray()
    dstCoords  = FSQuantityDescArray()
    fsmeshDonorPnt.DetermineCoordinates(dstCoords) or FSError.PrintAndExit()
    ipvr = FSMeshInterpolationVolumeReconstruction("sourcePoint")
    ipvr.InterpolateAtNodes(fsmesh, srcCoords, disc, state, fsmeshDonorPnt, dstCoords, quantities, False) or FSError.PrintAndExit()

    if check:
        fsmeshDonorPnt.ExportMeshTECPLOT(Filename='pointMesh_AugState.dat', FileFormat='ASCII', PrefixDatasetName=True) or FSError.PrintAndExit()
    nodeDonorPoints = fsmeshDonorPnt.GetUnstructDataset("DonorPoints").GetValues()
    nodeWallPoints  = fsmeshDonorPnt.GetUnstructDataset("WallPoints").GetValues()
    augStateData    = fsmeshDonorPnt.GetUnstructDataset("CODAAugState").GetValues()
    augStateNames   = fsmeshDonorPnt.GetUnstructDataset("CODAAugState").GetNames()

    nodeDonorPoints = Cmpi.allgather(numpy.array(nodeDonorPoints.Buffer(), copy=True))
    nodeWallPoints  = Cmpi.allgather(numpy.array(nodeWallPoints.Buffer() , copy=True))
    augStateData    = Cmpi.allgather(numpy.array(augStateData.Buffer()   , copy=True))

    Cmpi.trace('[INTERPOLATIONDONORPOINTS] Interpolation of the flow solution on the donor points..end', master=True)
    return nodeDonorPoints, nodeWallPoints, augStateData


def computeSurfValues(fileNameResultIn, tb, CODAInputs, dim=3, fileNameIBMPnts=None, fileNameResultOut=None, check=False, verbose=False, isRevertToOld=False):
    """ Surface quantity post-processing for CODA IBM computation. 
    Usage: computeSurfValues(fileNameResultIn, tb, CODAInputs, dim, fileNameIBMPnts, fileNameResultOut, check, verbose)"""

    ## CODAInputs = [discParaDictAllStages, discSelectionParaDict, IBM_markers, dictReferenceQuantities]
    ## fileNameResultIn = (e.g) output_stage3.h5 (the h5 from the CODA run)
    if isinstance(tb, str): tb = C.convertFile2PyTree(tb)
    else: tb = Internal.copyTree(tb)

    if not isinstance(fileNameResultIn,str):
        ValueError('fileNameResultIn MUST the name of the file. Importing of this file is done in this file through "fsmesh.ImportMeshHDF5".')
    if isRevertToOld and Cmpi.master:
        print("WARNING: Reverting to previous projection approach - extrapolation based...", flush=True)

    discParaDict            = CODAInputs[0]
    discSelectionParaDict   = CODAInputs[1]
    IBMMarkers              = CODAInputs[2]
    dictReferenceQuantities = CODAInputs[3]

    clac   = FSClac()
    fsmesh = FSMesh(clac)
    fsmesh.ImportMeshHDF5(Filename=fileNameResultIn) or FSError.PrintAndExit()

    # Reconstruction of the solution at the donor points (in parallel).
    donorPointsNumpy, wallPointsNumpy, augStateNumpy = interpolationDonorPoints(fsmesh, clac, discParaDict, discSelectionParaDict,
                                                                                IBMMarkers, check=check)

    # Serial part of the post-processing.
    zw = None
    if Cmpi.master:
        pytree = P_AMR.createPyTreeForIBMWallFieldsExtraction(numpy.concatenate(donorPointsNumpy),
                                                              numpy.concatenate(wallPointsNumpy),
                                                              numpy.concatenate(augStateNumpy),
                                                              discSelectionParaDict)
        if fileNameIBMPnts is not None: C.convertPyTree2File(pytree, fileNameIBMPnts)
        zw = P_AMR.extractIBMWallFields(pytree, tb, discSelectionParaDict, isRevertToOld=isRevertToOld)

    zw       = Cmpi.bcast(zw, root=0)
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
        if Cmpi.master: C.convertPyTree2File(zw, fileNameResultOut)
        return None
    else:
        return zw, coefs


def computeSurfValuesFSUI(fileNameResultIn, tb, fileNameRelations, dim=3, fileNameIBMPnts=None, fileNameResultOut=None, fileNameCoefOut='coefLiftDrag.txt', check=False, verbose=False, isRevertToOld=False):
    """ Surface quantity post-processing for CODA IBM computation using FSUI-CODA.
    Usage: computeSurfValues(fileNameResultIn, tb, fileNameRelations, dim, fileNameIBMPnts, fileNameResultOut, check, verbose)"""
    import json

    ## fileNameResultIn  = (e.g) solution.h5    (the h5 from the fsui-coda run)
    ## fileNameRelations = (e.g) relations.json (the json file created after fsui-coda create)
    if isinstance(tb, str): tb = C.convertFile2PyTree(tb)
    else: tb = Internal.copyTree(tb)

    if isRevertToOld and Cmpi.master:
        print("WARNING: Reverting to previous projection approach - extrapolation based...", flush=True)

    if not isinstance(fileNameResultIn,str):
        ValueError('fileNameResultIn MUST be the name of the file. Importing of this file is done in this file through "fsmesh.ImportMeshHDF5".')
    if not isinstance(fileNameRelations,str):
        ValueError('fileNameRelations MUST be the name of the file. Importing of this file is done in this file through "open(RelationsFileName)".')

    ## READ relations.json file
    with open(fileNameRelations, 'r') as f: data = json.load(f)

    discSelectionParaDict = data['CFD']['GeomIBM']['DiscretizationSelection']
    discParaDictTmp       = data['CFD']['GeomIBM']['DiscretizationPara']
    domain                = data['CFD']['GeomIBM']['Domain']
    # Convert domain to boundary treatments
    IBMMarkers = []
    bndy_treat  = []
    for boundary_name, boundary_data in domain.items():
        treatment = {
            "treatment type": boundary_data["FlowSolverBCType"],
            "boundary markers": boundary_data["CADGroupID"]
        }
        if 'Immersed' in boundary_data["FlowSolverBCType"]:
            IBMMarkers.append(boundary_data["CADGroupID"][0])
        # Check if there's a wall model in the CFD section
        if "CFD" in boundary_data and "wall model" in boundary_data["CFD"]:
            treatment["wall model"] = boundary_data["CFD"]["wall model"]
        bndy_treat.append(treatment)
    discParaDict = {
        **discParaDictTmp,
        "boundary treatments": bndy_treat
    }

    alpha    = discParaDict["reference state"]["flow direction specification"]["angle of attack"]
    beta     = discParaDict["reference state"]["flow direction specification"]["angle of sideslip"]
    Mach     = discParaDict["reference state"]["flow speed specification"]["Mach"]
    Reynolds = discParaDict["reference state"]["viscosity specification"]["Reynolds"]
    Lref     = discParaDict["reference state"]["viscosity specification"]["Reynolds_Length"]
    Aref     = discParaDict["boundary integral quantities"]["Coef_Area"]
    ## Post Values
    pressureRef    = 1
    densityRef     = 1
    gammaRef       = 1.4
    aRef           = (gammaRef*pressureRef/densityRef)**0.5
    velooRef       = Mach*aRef

    alpha_rad = numpy.deg2rad(alpha)
    beta_rad  = numpy.deg2rad(beta)
    velXRef = velooRef*numpy.cos(alpha_rad)*numpy.cos(beta_rad)
    velYRef = velooRef*numpy.cos(alpha_rad)*numpy.sin(beta_rad)
    velZRef = velooRef*numpy.sin(alpha_rad)

    dictReferenceQuantities = {
        "Mach_ref" : Mach,
        "Reynolds_ref" : Reynolds,
        "pressure_ref" : pressureRef,
        "density_ref" : densityRef,
        "velX_ref" : velXRef,
        "velY_ref" : velYRef,
        "velZ_ref" : velZRef,
        "alpha" : alpha,
        "beta" : beta,
        "Sref": Aref,
    }

    clac   = FSClac()
    fsmesh = FSMesh(clac)
    fsmesh.ImportMeshHDF5(Filename=fileNameResultIn) or FSError.PrintAndExit()

    # Reconstruction of the solution at the donor points (in parallel).
    donorPointsNumpy, wallPointsNumpy, augStateNumpy = interpolationDonorPoints(fsmesh, clac, discParaDict, discSelectionParaDict,
                                                                                IBMMarkers, check=check)

    # Serial part of the post-processing.
    zw = None
    if Cmpi.master:
        pytree = P_AMR.createPyTreeForIBMWallFieldsExtraction(numpy.concatenate(donorPointsNumpy),
                                                              numpy.concatenate(wallPointsNumpy),
                                                              numpy.concatenate(augStateNumpy),
                                                              discSelectionParaDict)
        if fileNameIBMPnts is not None: C.convertPyTree2File(pytree, fileNameIBMPnts)
        zw = P_AMR.extractIBMWallFields(pytree, tb, discSelectionParaDict, isRevertToOld=isRevertToOld)

    zw       = Cmpi.bcast(zw, root=0)
    zw,coefs = P_AMR.computeBoundaryQuantities(zw, dictReferenceQuantities, dim=dim, verbose=verbose) #ok

    if Cmpi.master:
        with open(fileNameCoefOut, 'w') as f:
            lines = [
                "\nIntegrated coefficients:",
                "\n ==== DRAG ==== :",
                "CD=      %g" % coefs[0],
                "CDfric=  %g" % coefs[2],
                "CDpres=  %g" % coefs[3],
                "\n ==== LIFT ==== :",
                "CL=      %g" % coefs[1],
                "CLfric=  %g" % coefs[4],
                "CLpres=  %g" % coefs[5],
            ]
            output = "\n".join(lines)
            print(output, flush=True)
            f.write(output + "\n") # file

    if fileNameResultOut is not None:
        if Cmpi.master: C.convertPyTree2File(zw, fileNameResultOut)
        return None
    else:
        return zw, coefs

##========================================================================
##========================================================================

def testMPIPart1(clac):
    #################### TEST MPI ########################### (temporary)
    nProcs = clac.GetNProcs()
    if clac.GetProcID()<10:
        indicesTestMPI = FSFloatArray(2)
        indicesTestMPI[0] = clac.GetProcID()
        indicesTestMPI[1] = clac.GetProcID()+10
    else:
        indicesTestMPI = FSFloatArray(2)
        indicesTestMPI.Fill(0.)
    #########################################################
    indicesTestMPINumpy     = numpy.array(indicesTestMPI.Buffer(), copy=True)
    indicesTestMPIGathered  = Cmpi.allgather(indicesTestMPINumpy)
    return indicesTestMPIGathered


def testMPIPart2(indicesTestMPINumpy, nProcs):
    ########## TEST MPI ###########   (temporary)
    if nProcs>9:
        testArrayIndices = [i for i in range(0,10)]
    else:
        testArrayIndices = [i for i in range(0,nProcs)]
    fullTestArrayIndices = numpy.zeros(nProcs*2)
    counter = 0
    for i in testArrayIndices:
        fullTestArrayIndices[counter] = i
        fullTestArrayIndices[counter+1] = i+10
        counter = counter+2
    if numpy.equal(numpy.concatenate(indicesTestMPINumpy),fullTestArrayIndices).all()==False:
        print("Error: full test array indices", fullTestArrayIndices)
        sys.exit("Problem with MPI. Gather not in the correct order")
    ##############################

    return None
