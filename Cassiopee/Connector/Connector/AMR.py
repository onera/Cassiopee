"""Connector for CODA/IBM preprocessing on AMR-based grids"""
import Converter
import Transform
import Converter.PyTree as C
import Converter.Internal as Internal
import Converter.Mpi as Cmpi
import Connector.PyTree as X
import Connector.IBM as X_IBM
import Connector.Mpi as Xmpi
import Connector.connector as connector
import Transform.PyTree as T
import Generator.PyTree as G
import Generator.Mpi as Gmpi
import Post.PyTree as P
import Dist2Walls.PyTree as DTW
import Generator.AMR as G_AMR
import time, copy, numpy
from Converter.Internal import E_NpyInt as E_NpyInt
from .QuadratureDG import *

__TOL__ = 1e-9

def outputTime(startTime,functionName='FunctionName'):
    endTime     = time.perf_counter()
    elapsedTime = endTime-startTime
    elapsedTime = Cmpi.allreduce(elapsedTime  ,op=Cmpi.MAX)
    if Cmpi.rank==0: print('Elapsed Time: %s: %g [s] | %g [min] | %g [hr]'%(functionName,elapsedTime,elapsedTime/60,elapsedTime/3600),flush=True)
    return None

def printMeshInfo(t, message=''):
    Nzones = len(Internal.getZones(t))
    Nzones = Cmpi.allreduce(Nzones, op=Cmpi.SUM)
    NCells = Cmpi.getNCells(t)
    if Cmpi.master:
        header = '' if message == '' else message + ' '
        print('%s Total number of zones=%d.'%(header, Nzones), flush=True)
        print('%s Final number of cells=%5.4f millions.'%(header, NCells*1e-6), flush=True)

    return None

# ===============================================================================================================================
def prepareAMRData(t_case, t, IBM_parameters=None, check=False, dim=3, localDir='./', forceAlignment=False, isFastApproach=True):
    Cmpi.trace('AMR prepare IBM...start', master=True)
    if "method" in IBM_parameters["IBM type"].keys() or IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
        t = prepareAMRDataDG__(t_case, t, IBM_parameters=IBM_parameters, check=check, dim=dim, localDir=localDir, forceAlignment=forceAlignment, isFastApproach=isFastApproach)
    else:
        t = prepareAMRDataFV__(t_case, t, IBM_parameters=IBM_parameters, check=check, dim=dim, localDir=localDir, forceAlignment=forceAlignment, isFastApproach=isFastApproach)
    Cmpi.trace('AMR prepare IBM...end', master=True, cpu=True)
    return t

def prepareAMRDataFV__(t_case, t, IBM_parameters=None, check=False, dim=3, localDir='./', forceAlignment=False, isFastApproach=True):
    VPM = False

    (IBM_parameters, frontTypeIP, frontTypeDP, dir_sym, different_front_flag) = checkInputsIbmParam__(IBM_parameters)

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    ## Note: tb2 is the correct geometry for 3D calcs.
    ##       2D: tb2 = tb with addkplane
    ##       3D: tb2 = tb
    bbo = Gmpi.bbox(t)
    if dim == 2:
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]
        tb2 = C.initVars(tb, 'CoordinateZ', bb0[2]) # forced
        T._addkplane(tb2)
        T._contract(tb2, (0,0,0), (1,0,0), (0,1,0), dz)
    else:
        baseSYM = Internal.getNodesFromName1(tb, "SYM")
        if baseSYM:
            tb2 = Internal.rmNodesByNameAndType(tb, 'SYM', 'CGNSBase_t')
            tb2 = Internal.rmNodesByNameAndType(tb2, '*_sym*', 'Zone_t')
        else:
            tb2 = tb

    #==========================================================
    # STEP 0: Calculate Distance to IBCs (all IBCs) (if needed)
    #==========================================================
    _dist2wallIBM(t, tb2, dim, different_front_flag)

    #==========================================================
    # STEP 1: Blanking by immersed body
    #==========================================================
    Cmpi.trace(">>> Blanking [start]", master=True, cpu=False)
    t = blankingIBM(t, tb2)
    Cmpi.trace(">>> Blanking [end]  ", master=True, cpu=False)

    #=======================================================================
    # STEP 2: Save BC Names & Types - check if QuadNQuad is fully inside IBM
    #=======================================================================
    BCInfo = getBCs(t, tb2, dim)

    #===============================================================================================================================
    # STEP 3: Get Integration Point Front (Recall: IP = integration point (NOT image point) - Target cells in Mittal et al. approach)
    #===============================================================================================================================
    # only done for frontTypeIP=2 - if frontTypeIP=1 it is F1 so blankByIBCBodies is sufficient for the cellN value
    Cmpi.trace("Extract front faces of IBM integration points [start] ", master=True, cpu=False)
    frontIP, maxDistIP = getFrontIP(t, IBM_parameters, localDir, False, VPM=VPM)
    Cmpi.trace("Extract front faces of IBM integration points [end]   ", master=True, cpu=False)

    #====================================================
    # STEP 4: Select Cells - only fluid cells from now on
    #====================================================
    Cmpi.trace(" Removing blanked cells [start]", master=True, cpu=False)
    t = removeBlankedCells(t)
    Cmpi.trace(" Removing blanked cells [end]  ", master=True, cpu=False)

    #===============================
    # STEP 5: Recover BCs - Add IBCs
    #===============================
    Cmpi.trace(" Recovering Boundary Conditions [start]", master=True, cpu=False)
    t, t_exteriorFaces, zibc = recoverBCs(t, frontIP, BCInfo)
    Cmpi.trace(" Recovering Boundary Conditions [end]  ", master=True, cpu=False)

    #===================================================================================
    # STEP 6: Get Donor Point Front (Recall: DP = image point in Mittal et. al approach)
    #         Determine DP points on DPfront
    #         Add IBC Dataset to the PyTree
    #===================================================================================
    Cmpi.trace(" Extracting front of the donor points [start]", master=True, cpu=False)
    frontDP = None
    if frontTypeDP == '1':
        frontDP = getFrontDP(t, tb2, frontIP, dim, dir_sym, check, distIP=maxDistIP, localDir=localDir, isFastApproach=isFastApproach)
    Cmpi.trace(" Extracting front of the donor points [end]  ", master=True, cpu=False)

    #===================================================================================
    # STEP 7: Determine DP points on DPfront & Add IBC Dataset to the PyTree
    #===================================================================================
    Cmpi.trace(" Getting DP IBM points and adding IBC Dataset [start]", master=True, cpu=False)
    _setIBCData(t, t_exteriorFaces, tb2, zibc, frontDP, bbo, IBM_parameters, check, forceAlignment, localDir, different_front_flag)
    Cmpi.trace(" Getting DP IBM points and adding IBC Dataset [end]  ", master=True, cpu=False)

    return t

def prepareAMRDataDG__(t_case, t, IBM_parameters=None, check=False, dim=3, localDir='./', forceAlignment=False, isFastApproach=True):
    VPM = False

    (IBM_parameters, frontTypeIP, frontTypeDP, dir_sym, different_front_flag) = checkInputsIbmParam__(IBM_parameters)

    # VPM = volume penalization method - primarily used by Universidad Politécnica de Madrid (UPM)
    if "method" in IBM_parameters["IBM type"].keys():
        if IBM_parameters["IBM type"]["method"] == "VPM":
            VPM = True

    if IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
        if Cmpi.master: print("Warning: You are using high-order DG/DGSEM spatial discretizations in parallel. This is a development version. For a more validated and robust high-order IBM-preprocessing, switch to serial. ", flush=True)

        if different_front_flag == False:
            if Cmpi.master: print("Warning:High-order DG/DGSEM spatial discretizations  \n \"use different front for different BCs ==False\" is not implemented.\n Using \"use different front for different BCs==True\" instead.", flush=True)

    if isinstance(t_case, str): tb = C.convertFile2PyTree(t_case)
    else: tb = t_case

    ## Note: tb2 is the correct geometry for 3D calcs.
    ##       2D: tb2 = tb with addkplane
    ##       3D: tb2 = tb
    bbo = Gmpi.bbox(t)
    if dim == 2:
        z0 = Internal.getNodeFromType2(t, "Zone_t")
        bb0 = G.bbox(z0); dz = bb0[5]-bb0[2]
        tb2 = C.initVars(tb, 'CoordinateZ', bb0[2]) # forced
        T._addkplane(tb2)
        T._contract(tb2, (0,0,0), (1,0,0), (0,1,0), dz)
    else:
        baseSYM = Internal.getNodesFromName1(tb, "SYM")
        if baseSYM:
            tb2 = Internal.rmNodesByNameAndType(tb, 'SYM', 'CGNSBase_t')
            tb2 = Internal.rmNodesByNameAndType(tb2, '*_sym*', 'Zone_t')
        else:
            tb2 = tb

    #==========================================================
    # STEP 0: Calculate Distance to IBCs (all IBCs) (if needed)
    #==========================================================
    _dist2wallIBM(t, tb2, dim, different_front_flag)

    #==========================================================
    # STEP 1: Blanking by immersed body
    #==========================================================
    Cmpi.trace(">>> Blanking [start]", master=True, cpu=False)
    t = blankingIBM(t, tb2)
    Cmpi.trace(">>> Blanking [end]  ", master=True, cpu=False)

    #=======================================================================
    # STEP 2: Save BC Names & Types - check if QuadNQuad is fully inside IBM
    #=======================================================================
    BCInfo = getBCs(t, tb2, dim)

    #===============================================================================================================================
    # STEP 3: Get Integration Point Front (Recall: IP = integration point (NOT image point) - Target cells in Mittal et al. approach)
    #===============================================================================================================================
    # only done for frontTypeIP=2 - if frontTypeIP=1 it is F1 so blankByIBCBodies is sufficient for the cellN value
    Cmpi.trace("Extract front faces of IBM integration points [start] ", master=True, cpu=False)
    frontIP, maxDistIP = getFrontIP(t, IBM_parameters, localDir, False, VPM=VPM)
    Cmpi.trace("Extract front faces of IBM integration points [end]   ", master=True, cpu=False)

    #====================================================
    # STEP 4: Select Cells - only fluid cells from now on
    #====================================================
    Cmpi.trace(" Removing blanked cells [start]", master=True, cpu=False)
    t = removeBlankedCells(t)
    Cmpi.trace(" Removing blanked cells [end]  ", master=True, cpu=False)

    #===============================
    # STEP 5: Recover BCs - Add IBCs
    #===============================
    Cmpi.trace(" Recovering Boundary Conditions [start]", master=True, cpu=False)
    t, t_exteriorFaces, zibc = recoverBCs(t, frontIP, BCInfo)
    Cmpi.trace(" Recovering Boundary Conditions [end]  ", master=True, cpu=False)

    #===================================================================================
    # STEP 6: Get Donor Point Front (Recall: DP = image point in Mittal et. al approach)
    #         Determine DP points on DPfront
    #         Add IBC Dataset to the PyTree
    #===================================================================================
    if VPM == False:
        Cmpi.trace(" Extracting front of the donor points [start]", master=True, cpu=False)
        frontDP = None
        if frontTypeDP == '1':
            frontDP = getFrontDP(t, tb2, frontIP, dim, dir_sym, check, distIP=maxDistIP, localDir=localDir, isFastApproach=isFastApproach)
        Cmpi.trace(" Extracting front of the donor points [end]  ", master=True, cpu=False)

        #===================================================================================
        # STEP 7: Determine DP points on DPfront & Add IBC Dataset to the PyTree
        #===================================================================================
        Cmpi.trace(" Getting DP IBM points and adding IBC Dataset [start]", master=True, cpu=False)
        _setIBCDataDG(t, t_exteriorFaces, tb2, zibc, frontDP, bbo, IBM_parameters, check, forceAlignment, localDir, different_front_flag)
        Cmpi.trace(" Getting DP IBM points and adding IBC Dataset [end]  ", master=True, cpu=False)
    else:
        #===================================================================================
        # STEP 7: Determine DP points on DPfront & Add IBC Dataset to the PyTree
        #===================================================================================
        ncellsIP = C.getNCells(zibc)
        if ncellsIP>0: _addIBC2Zone__(t, t_exteriorFaces, zibc)


    return t

# ===============================================================================================================================
def prepareAMRIBM(tb, vmins, dim, IBM_parameters, levelMax=0, toffset=None, check=False, opt=False, octreeMode=1,
                  snears=0.01, dfars=10, loadBalancing=False, OutputAMRMesh=False,
                  localDir='./', fileName=None, tbox=None, vminsTbox=5, forceAlignment=False,
                  tIn=None, isFastApproach=True, **kwargs):
    """Generate AMR IBM mesh and prepare AMR IBM data for CODA simulation. 
    Usage: prepareAMRIBM(tb, levelMax, vmins, dim, IBM_parameters, toffset, check, opt, octreeMode,
                         snears, dfars, loadBalancing, OutputAMRMesh, localDir, fileName, tbox, vminsTbox, tbv2, forceAlignment)"""

    import gc

    # debug parameters
    tbv2 = kwargs.get('tbv2', None)

    ## =========================
    ## ==== Mesh Generation ====
    ## =========================
    t_AMR = G_AMR.generateAMRMesh(tb=tb, levelMax=levelMax, vmins=vmins, dim=dim,
                                  toffset=toffset, check=check, opt=opt, octreeMode=octreeMode, localDir=localDir,
                                  snears=snears, dfars=dfars, loadBalancing=loadBalancing,
                                  tbox=tbox, vminsTbox=vminsTbox, tbv2=tbv2,
                                  tIn=tIn)

    if OutputAMRMesh: Cmpi.convertPyTree2File(t_AMR, localDir+'tAMRMesh.cgns')

    printMeshInfo(t_AMR, message='[MESH GEN.]')

    ### Clear memory
    Cmpi.trace("AMR Memory clean & memory check...start", master=True)
    gc.collect()
    Cmpi.trace("AMR Memory clean & memory check...end", master=True)
    Cmpi.barrier()

    ## ==================
    ## ==== IBM Prep ====
    ## ==================
    t_AMR = prepareAMRData(tb, t_AMR, IBM_parameters=IBM_parameters, dim=dim, check=check, localDir=localDir,
                           forceAlignment=forceAlignment, isFastApproach=isFastApproach)

    printMeshInfo(t_AMR, message='[IBM PREP.]')

    if fileName is not None:
        Cmpi.convertPyTree2File(t_AMR, localDir+fileName)
        return None
    else:
        return t_AMR

# ===============================================================================================================================
def checkInputsIbmParam__(IBM_parametersIn):
    IBM_parameters = copy.deepcopy(IBM_parametersIn)

    frontTypeIP = IBM_parameters["integration points"]["front type"]
    if frontTypeIP not in ["1","2"]:
        raise ValueError("FrontTypeIP not implemented: only frontTypeIP==\"1\" and \"2\" are implemented in parallel.")
        Cmpi.abort(errorcode=1)

    frontTypeDP = IBM_parameters["donor points"]["front type"]
    if frontTypeDP not in ["1","2"]:
        raise ValueError("FrontTypeDP not implemented: only frontTypeDP==\"1\" and \"2\" are implemented in parallel.")
        Cmpi.abort(errorcode=1)

    if frontTypeIP == "1":
        depth_IP = IBM_parameters["integration points"]["depth IntegrationPoints"]
        if depth_IP != 0:
            if Cmpi.master: print("Warning: Only depth_IP=0 is implemented in parallel-AMR. Continuing with depthIP=0.", flush=True)
            IBM_parameters["integration points"]["depth IntegrationPoints"] = 0

    if frontTypeDP == "1":
        depth_DP = IBM_parameters["donor points"]["depth DonorPoints"]
        if depth_DP != 1:
            if Cmpi.master: print("Warning: Only depth_DP=1 is implemented in parallel-AMR. Continuing with depthDP=1.", flush=True)
            IBM_parameters["donor points"]["depth DonorPoints"] = 1

    # symmetry & direction
    dir_sym = 0
    if "symmetryPlane" in IBM_parameters["IBM type"].keys():
        dir_sym = int(IBM_parameters["IBM type"]["symmetryPlane"])
        if IBM_parameters["IBM type"]["symmetryPlane"] < 0 or IBM_parameters["IBM type"]["symmetryPlane"] > 3:
            if Cmpi.master: print("Warning: Symmetry plane direction can only be : 1 (x-direction), 2 (y-direction) or 3 (z-direction)... exiting", flush=True)
            raise ValueError("Choose a valid symmetry plane direction. Exiting..")
            Cmpi.abort(errorcode=1)

    # Important Note: this use of the flag is still ambiguous - related to local IBMs??
    different_front_flag = True
    if "use different front for different BCs" in IBM_parameters["integration points"]:
        different_front_flag = IBM_parameters["integration points"]["use different front for different BCs"]

    return (IBM_parameters, frontTypeIP, frontTypeDP, dir_sym, different_front_flag)

# ===============================================================================================================================
def removeBlankedCells(t):
    t = P.selectCells(t, '{cellN}==1.', strict=1)
    # Make sure that the only node of type Elements_t is 'GridElements'
    for node in Internal.getNodesFromType(t, 'Elements_t'):
        if node[0] != 'GridElements':
            Internal._rmNode(t, node)

    Internal._rmNodesFromName(t, Internal.__FlowSolutionNodes__)
    Internal._rmNodesFromType(t, 'Family_t')

    for z in Internal.getZones(t): Cmpi._setProc(z, Cmpi.rank)

    return t

# ===============================================================================================================================
def getBCs(t, tb2, dim):
    # Identity the BCTypes & BCNames in t
    zbcs=[]; bctypes=[]; bcnames=[]
    for bc in Internal.getNodesFromType(t, 'BC_t'):
        bctype = Internal.getValue(bc)
        bcname = Internal.getName(bc)
        if bctype not in bctypes:
            bctypes.append(bctype)
            bcnames.append(bcname)

    # Save the boundary conditions for later use
    for bctype in bctypes:
        zbc = C.extractBCOfType(t, bctype)
        Internal._rmNodesByType(zbc, "FlowSolution_t")
        zbc = T.join(zbc)
        zbcs.append(zbc)

    # Blanking to check if the QuadNQuad BC is inside the geometry
    # Needed to avoid wrong BCs
    # Get largest length of the bases
    tTmp = G.getVolumeMap(t)
    hminGlobal = (C.getMinValue(tTmp,"centers:vol"))**(1/dim)
    hminGlobal = Cmpi.allreduce(hminGlobal, op=Cmpi.MIN)
    del tTmp
    L1 = 0.0
    for bodyLocal in Internal.getBases(tb2):
        bb1 = G.bbox(bodyLocal)
        L1 = max(L1, bb1[3]-bb1[0])
        L1 = max(L1, bb1[4]-bb1[1])
        if dim == 3: L1 = max(L1, bb1[5]-bb1[2])

    for nobc, zbc in enumerate(zbcs):
        if bcnames[nobc] == "QuadNQuad":
            XRAYDIM1 = int(L1/hminGlobal) + 10;
            XRAYDIM1 = max(1500, min(15000, XRAYDIM1)) #x3
            bodies = [Internal.getBases(tb2)]; nbodies = len(Internal.getBases(tb2))
            BM = numpy.ones((1, nbodies), dtype=Internal.E_NpyInt)
            zbcTemp = C.newPyTree(["BASE", Internal.getZones(zbc)])
            # Check if BCs is inside the geometry
            zbcTemp = X.blankCells(zbcTemp, bodies, BM, blankingType='center_in', dim=dim, XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM1)
            maxBlankVal = C.getMaxValue(zbcTemp, 'centers:cellN')
            # If BC is entirely inside the geometry - change its name
            if maxBlankVal < 1: bcnames[nobc] = 'QuadNQuad_Empty'
            del zbcTemp
            del bodies

    return (zbcs, bctypes, bcnames)

# ===============================================================================================================================
def recoverBCs(t, frontIP, BCInfo):
    tp = Internal.copyRef(t)
    t_exteriorFaces = P.exteriorFaces(tp)

    zbcs, bctypes, bcnames = BCInfo
    for elt_t in Internal.getNodesFromType(t_exteriorFaces, "Elements_t"):
        if not elt_t[0].startswith("GridElements"):
            Internal._rmNode(t_exteriorFaces, elt_t)
    _recoverBCs__(tp, t_exteriorFaces, zbcs, bctypes, bcnames)
    zibc, t_exteriorFaces = recoverIBCs__(t_exteriorFaces, frontIP)

    return tp, t_exteriorFaces, zibc

def _recoverBCs__(t, t_exteriorFaces, zbcs, bctypes, bcnames):
    meshgen = "AMR"
    f = None
    for z in Internal.getZones(t):
        if z is not None:
            nobc = len(zbcs)
            f = Internal.getZones(t_exteriorFaces)[0]
            if Cmpi.master: print("Performing the 'identifyElements' function (it can be long.)", flush=True)
            for nobc, zbc in enumerate(zbcs):
                hook = C.createHook(f, "elementCenters")
                # Indices of the elements of f corresponding to the elements of zbc
                # Note: zbc is before selectCells (mesh generated with G_AMR)
                #       f (hook) is after selectCells
                #       Check which elements in f correspond to the zbc
                ids = C.identifyElements(hook, zbc, tol=__TOL__)
                len_ids = Internal.getValue(f)[0][1]
                ids = ids[ids[:] > -1] - 1 # consider the positive numbers only & index starts at 0
                ids = ids.tolist()
                #C.freeHook(hook)
                if len(ids) > 0:
                    # subzone of elements that are in zbc & f
                    # zf: elements of that BC that have to be conserved
                    zf = T.subzone(f, ids, type='elements')
                    if bcnames[nobc] != "QuadNQuad":
                        G_AMR._addBC2Zone__(z, bctypes[nobc], bctypes[nobc], zf)
                    else:
                        G_AMR._addBC2Zone__(z, "QuadNQuad", "FamilySpecified:QuadNQuad", zf)
                    # ids_all: all (match & unmatched) ids for that BC
                    # ids_new: consider the ids of the unmatched elements for that BC
                    ids_all = list(range(len_ids))
                    ids_new = list(set(ids_all)-set(ids))
                    if len(ids_new) > 0:
                        f = T.subzone(f, ids_new, type='elements')
                elif len(ids) == 0 and bcnames[nobc] == "QuadNQuad":
                    elts = Internal.getNodesFromType1(z, "Elements_t")
                    maxElt = Internal.getNodeFromName(elts[-1], "ElementRange")[1][1]
                    CODABCType = "QuadNQuad"
                    Internal.newElements(name=CODABCType, etype=7, econnectivity=numpy.empty(0),
                                         erange=[maxElt+1, maxElt], eboundary=1, parent=z)
                    C._addBC2Zone(z, CODABCType, "FamilySpecified:"+CODABCType, elementRange=[maxElt+1,maxElt])
                    zone_bc = Internal.getNodeFromType1(z, 'ZoneBC_t')
                    lastbcname = C.getLastBCName(CODABCType)
                    node_bc = Internal.getNodeFromName(zone_bc, lastbcname)
                    node_bc[0] = CODABCType
                C.freeHook(hook)

            z[0] = z[0]+str(Cmpi.rank)
    # t_exteriorFaces - should only contain the integration front BC condition
    #                   CODA needs this as it applies the IBC on this front.
    #                   This is the unmatched BC that will become the IBC
    if meshgen == "AMR" and f is not None: t_exteriorFaces[2][1][2] = [f]
    return None

def recoverIBCs__(t_exteriorFaces, frontIP):
    # Here:
    # t_exteriorFaces - is ONLY the exteriorFaces of the integration front on which CODA applies the IBCs
    Cmpi.trace(" Adding the IBC BC tag (per processor) for CFD solver [start]", master=True, cpu=False)
    if Cmpi.master: print("Performing the 'identifyElements' function (it can be long.)", flush=True)
    startTime = time.perf_counter()
    f = Internal.getZones(t_exteriorFaces)
    if f != []:
        f = f[0]
        hook = C.createHook(f,"elementCenters")
        ids = C.identifyElements(hook, frontIP, tol=__TOL__)
        ids = ids[ids[:] > -1]
        ids = ids.tolist()
        ids_IBMWall = [ids[i]-1 for i in range(len(ids))]
        C.freeHook(hook)
        if ids_IBMWall != []:
            zibc = T.subzone(f, ids_IBMWall, type='elements')
        else:
            # Needed for MPI all gather
            zibc = Internal.newZone(name="frontIP%d"%Cmpi.rank, zsize=[[0,0]], ztype="Unstructured")
    else:
        # Needed for MPI all gather
        zibc = Internal.newZone(name="frontIP%d"%Cmpi.rank, zsize=[[0,0]], ztype="Unstructured")
    outputTime(startTime,functionName='identifyElementsPrt2')
    Cmpi.trace(" Adding the IBC BC tag (per processor) for CFD solver [end]", master=True, cpu=False)

    return zibc, f

# ===============================================================================================================================
def blankingIBM(t, tb):
    C._initVars(t, 'cellN', 1.)
    t = X_IBM.blankByIBCBodies(t, tb, 'nodes', 3)
    C._initVars(t,'{TurbulentDistance}=-1.*({cellN}<1.)*{TurbulentDistance}+({cellN}>0.)*{TurbulentDistance}')

    return t

# ===============================================================================================================================
def _dist2wallIBM(t, tb, dim, different_front_flag):
    if different_front_flag: # True is default
        tbLocal = getBodiesDist2wall__(tb) # for local IBMs
    else:
        tbLocal = Internal.copyRef(tb)

    varnames = C.getVarNames(t, loc="nodes")[0]
    if "TurbulentDistance" not in varnames:
        Cmpi.trace(">>> Wall distance nodes [start]", master=True, cpu=False)
        DTW._distance2Walls(t, tbLocal, type='ortho', signed=0, dim=dim, loc='nodes')
        Cmpi.trace(">>> Wall distance nodes [end]  ", master=True, cpu=False)
    else:
        Cmpi.trace(">>> Wall distance nodes : skipped - dist2wall is in input PyTree ", master=True, cpu=False)

    varnames = C.getVarNames(t, loc="centers")[0]
    if "TurbulentDistance" not in varnames:
        Cmpi.trace(">>> Wall distance centers [start]", master=True, cpu=False)
        DTW._distance2Walls(t, tbLocal, type='ortho', signed=0, dim=dim, loc='centers')
        Cmpi.trace(">>> Wall distance centers [end]  ", master=True, cpu=False)
    else:
        Cmpi.trace(">>> Wall distance centers : skipped - dist2wall is in input PyTree ", master=True, cpu=False)

    return None

def getBodiesDist2wall__(tb):
    zones_tb = Internal.getZones(tb)
    zones_tb_WD = []
    for z_tb in zones_tb:
        ibctype = Internal.getNodeFromName(z_tb, "ibctype")
        if ibctype == None:
            ibctype = 0
        elif isinstance(Internal.getValue(ibctype), str):
            print(Internal.getValue(ibctype))
            ibctype = 0
        else:
            ibctype = ibctype[1][0]
        if ibctype == 0:
            zones_tb_WD.append(z_tb)
    tb_WD = C.newPyTree(["tbWD", zones_tb_WD])
    return tb_WD

# ===============================================================================================================================
def getFrontIP(t, IBM_parameters, localDir, check=False, VPM=False):
    frontIP = extractFrontIP__(t, IBM_parameters, VPM=VPM)

    maxDistIP = 0.0
    turbDistanceTmp = Internal.getNodeFromName(frontIP, 'TurbulentDistance')[1]
    if len(turbDistanceTmp)>0: maxDistIP = C.getMaxValue(frontIP, 'TurbulentDistance')
    maxDistIP = Cmpi.allreduce(maxDistIP, op=Cmpi.MAX)
    if Cmpi.master: print('extractFrontDP__ Info: maxDistanceFrontIP=%g'%maxDistIP, flush=True)

    frontIP = gatherFrontIP__(frontIP, localDir, check)

    return frontIP, maxDistIP

def extractFrontIP__(t, IBM_parameters, VPM=False):
    if VPM == False:
        frontTypeIP = IBM_parameters["integration points"]["front type"]
        print("Rank: %d :: frontTypeIP=%d"%(Cmpi.rank, int(frontTypeIP)), flush=True)

        if frontTypeIP == "2":
            distance_IP = IBM_parameters["integration points"]["distance IntegrationPoints"]
            C._initVars(t, 'distance_IP', distance_IP)
            C._initVars(t, '{cellN}=({TurbulentDistance}>{distance_IP})*{cellN}')
    else:
        snear = IBM_parameters["IBM type"]["size elements body"]
        C._initVars(t, 'distance_IP', 2*snear)
        C._initVars(t, '{cellN}=1-({TurbulentDistance}<-{distance_IP})')

    frontIP = P.frontFaces(t, 'cellN')
    return frontIP

def gatherFrontIP__(frontIP, localDir, check):
    Cmpi.trace("Gathering front IP [start]", master=True, cpu=False)
    frontIP = Internal.getZones(frontIP)[0]
    frontIP_gath = Cmpi.allgatherZones(frontIP)
    frontIP_gath = C.newPyTree(["frontIP", frontIP_gath])
    frontIP_gath = T.join(frontIP_gath)
    frontIP_gath = G.close(frontIP_gath)
    for node in Internal.getNodesFromType(frontIP_gath, "Elements_t"):
        if node[0] != "GridElements": Internal._rmNode(frontIP_gath, node)
    Cmpi.trace("Gathering front IP [end]  ", master=True, cpu=False)

    if Cmpi.master and check:
        print("Exporting frontIP..", flush=True)
        C.convertPyTree2File(frontIP_gath, localDir+"frontIP_gath.plt")
        C.convertPyTree2File(frontIP_gath, localDir+"frontIP_gath.cgns")

    return frontIP_gath

# ===============================================================================================================================
def getFrontDP(t, tb, frontIP, dim, dir_sym, check, distIP, localDir='./', isFastApproach=True):
    frontDP = extractFrontDP__(t, tb, frontIP, dim, dir_sym, distIP, isFastApproach)
    frontDP = gatherFrontDP__(frontDP, localDir, check, isFastApproach)

    return frontDP

def extractFrontDP__(t, tb, frontIP_gath, dim, dir_sym, distIP, isFastApproach=True):
    import Geom.IBM as D_IBM
    if dim == 2 and not isFastApproach:
        isFastApproach = True
        if Cmpi.master: print("extractFrontDP__: for 2D test cases... Robust approach == Fast approach.", flush=True)
    if isFastApproach:
        startTimeExtract = time.perf_counter()
        if Cmpi.master: print("extractFrontDP__ - using Fast approach based on the integration points. This approach may yield unsatisfactory results for small resolutions", flush=True)
        ##Orig Approach - Based on Integration points front (frontIP_gath)
        ##                fast approach but can lead to errors in the IBM points - encountered when running CODA
        ##                Recall: mushroom clouds on CRM
        if Cmpi.master:
            C._deleteEmptyZones(frontIP_gath)
            frontIP_gath = T.join(frontIP_gath)
            frontIP_gath = C.convertArray2Tetra(frontIP_gath)
            frontIP_gath = G.close(frontIP_gath)
            frontIP_gath[0] = "frontIP_gath"
            if dim == 3 and dir_sym > 0:
                print("Symmetry of frontIP: Sym. Plane: %d"%dir_sym)
                frontIP_gath = C.newPyTree(["Base", frontIP_gath])
                frontIP_gath= Internal.getNodeFromName(frontIP_gath, 'Base')
                minval = C.getMinValue(frontIP_gath, ['CoordinateX', 'CoordinateY','CoordinateZ'])
                minval = minval[dir_sym-1]
                if dir_sym == 1: symPlane=(minval,0,0)
                elif dir_sym == 2: symPlane=(0,minval,0)
                else: symPlane=(0,0,minval)
                D_IBM._symmetrizeBody(frontIP_gath, dir_sym=dir_sym, symPlane=symPlane) # expect base input
                frontIP_gath = G.close(frontIP_gath)
                frontIP_gath = T.join(frontIP_gath)
        frontIP_gath = Cmpi.bcastZone(frontIP_gath)
        frontIP_gath = C.newPyTree(["Base", frontIP_gath])
        C._initVars(t, 'cellNFront', 1.)
        X_IBM._blankByIBCBodies(t, frontIP_gath, 'nodes', 3, cellNName="cellNFront")
        del frontIP_gath
        frontDP = P.frontFaces(t, 'cellNFront')
        del t
        outputTime(startTimeExtract,functionName='extractFrontDP__ - Fast Approach')
    else:
        startTimeExtract = time.perf_counter()
        if Cmpi.master:
            print("extractFrontDP__ - using Robust approach based on the offsets & dist2wall. This approach can take some time.", flush=True)
            if dir_sym > 0: print("Symmetry of frontIP: Sym. Plane: %d"%dir_sym, flush=True)
        ## Robust - based on tb (input geomtery), offset, selectcells, & dist2wall approach
        ##          more expensive but proven to be more robust

        # Get snear
        G._getVolumeMap(t)
        hminTmp = (C.getMinValue(t,"centers:vol"))**(1/dim)
        hminTmp = Cmpi.allreduce(hminTmp, op=Cmpi.MIN)
        if Cmpi.master: print('extractFrontDP__ Info: Smallest cell size (snear): %g'%hminTmp, flush=True)
        # Generate Offset - scaled tb
        frontIP_gathScale = localOffset__(tb, dim=dim, dir_sym=dir_sym, minSnear=hminTmp, distIP=distIP)
        #Cmpi.convertPyTree2File(frontIP_gathScale, 'check_frontIP_gathScale.cgns') # Keep for now - debugging

        # blankcells - what is inside the offset
        C._initVars(t, 'cellNTmp', 1.)
        bodiesTmp = [Internal.getZones(frontIP_gathScale)]
        nbodies = len(bodiesTmp)
        nbases = len(Internal.getBases(t))
        XRAYDIM1 = 2500
        BM = numpy.ones((nbases,nbodies),dtype=Internal.E_NpyInt)
        t = X.blankCells(t, bodiesTmp, BM, blankingType='node_in', XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM1,
                         dim=dim, cellNName='cellNTmp')

        # select cells - small region to do dist2wall
        tTmp = P.selectCells(t, "{cellNTmp}<1", strict=0)
        C._rmVars(t,['cellNTmp','vol'])

        # dist2wall on region of select cells
        if len(Internal.getZones(tTmp))>0:
            DTW._distance2Walls(tTmp, frontIP_gath, type='ortho', signed=0, dim=dim, loc='nodes')
            # new cellNFront
            C._initVars(tTmp, '{cellNFront}=({TurbulentDistance}>0.9*%g)'%hminTmp)
            #C.convertPyTree2File(tTmp, 'check_tTmp_final_proc%d.cgns'%Cmpi.rank)
            del frontIP_gath
            del t
            frontDP = P.frontFaces(tTmp, 'cellNFront')
            #C.convertPyTree2File(frontDP, 'check_frontDP_proc%d.cgns'%Cmpi.rank)
            del frontIP_gathScale
            del tTmp
        else:
            frontDP = Internal.newZone(name="front", zsize=[[0,0]], ztype="Unstructured")
            gc = Internal.newGridCoordinates(parent=frontDP)
            Internal.newDataArray('CoordinateX', value=numpy.empty(0), parent=gc)
            Internal.newDataArray('CoordinateY', value=numpy.empty(0), parent=gc)
            Internal.newDataArray('CoordinateZ', value=numpy.empty(0), parent=gc)
        outputTime(startTimeExtract,functionName='extractFrontDP__ - Robust Approach')
    ## Continue - same as orig.

    return frontDP

def gatherFrontDP__(frontDP, localDir, check, isFastApproach=True):
    Cmpi.trace("Gathering front DP [start]", master=True, cpu=False)
    frontDP_gath = Cmpi.allgatherZones(frontDP)
    C._deleteEmptyZones(frontDP_gath)
    frontDP_gath = T.join(frontDP_gath)
    frontDP_gath = C.newPyTree(["frontDP", frontDP_gath])
    if Cmpi.master and check:
        print("Exporting front DP...", flush=True)
        if isFastApproach:
            C.convertPyTree2File(frontDP_gath, localDir+"frontDP_gath_FastApproach.cgns")
            C.convertPyTree2File(frontDP_gath, localDir+"frontDP_gath_FastApproach.plt")
        else:
            C.convertPyTree2File(frontDP_gath, localDir+"frontDP_gath_RobustApproach.cgns")
            C.convertPyTree2File(frontDP_gath, localDir+"frontDP_gath_RobustApproach.plt")

    return frontDP_gath

# ===============================================================================================================================
# To be replaced with Generator.AMR counterpart
# ===============================================================================================================================
def localOffset__(tbLocal, dim, dir_sym, minSnear, distIP):
    # A lot of redundancies with Generator/AMR.py - TODO: can some parts be generalized
    import Geom.IBM as D_IBM

    distOffset = distIP + 7*minSnear # 5 (real) & 2 for security

    # [Connector/AMR.py specific]
    # Copy tb about the symmetry plane
    tbSym = Internal.copyTree(tbLocal)
    if dir_sym > 0:
        D_IBM._setSnear(tbSym, minSnear)
        D_IBM._setDfar(tbSym, 10)
        D_IBM._symmetrizePb(tbSym, 'Base', snear_sym=minSnear, dir_sym=dir_sym)
        baseSYM = Internal.getNodesFromName1(tbSym, "SYM")
        if baseSYM is not None: tbSym=Internal.rmNodesByNameAndType(tbSym, 'SYM', 'CGNSBase_t')
        C._rmVars(tbSym,['centers:cellN'])

    # [Connector/AMR.py specific]
    # coarsening the tb offset - like in Generator/AMR.py
    for nob in range(len(tbSym[2])):
        if Internal.getType(tbSym[2][nob]) == 'CGNSBase_t':
            z = Internal.getZones(tbSym[2][nob])
            z = C.convertArray2Tetra(z)
            z = T.join(z)
            bbz = G.bbox(z)
            # [TODO] this needs to be related to the hmin
            # hausd is a length and must be adapted to the dimensions of each case
            hausd = max(bbz[3]-bbz[0], bbz[4]-bbz[1], bbz[5]-bbz[2])/10000.
            hmax = hausd*1000
            if Cmpi.master:
                print('Remeshing (tbSym) surface mesh --> Maximum chordal deviation between final and initial mesh::%g || Maximum mesh step in final mesh::%g'%(hausd, hmax), flush=True)
            # exteriorFaces currently crashes if the surface is closed
            try:
                fixedConstraints = P.exteriorFaces(z)
            except:
                fixedConstraints = []
            z = G.mmgs(z, hausd=hausd, hmax=hmax, fixedConstraints=fixedConstraints)
            tbSym[2][nob][2] = Internal.getZones(z)

    # Background mesh for offset - only around orig. tb & not its symmetry one
    BB = G.bbox(tbLocal)
    ni = 150; nj = 150; nk = 150
    XRAYDIM1 = 3*ni; XRAYDIM2 = 3*nj

    # CARTRX
    delta2 = max(BB[3]-BB[0], BB[4]-BB[1], BB[5]-BB[2])*0.02 # 2% seems enough for the external cases already tested
    xmin_core = BB[0]-delta2
    ymin_core = BB[1]-delta2
    zmin_core = BB[2]-delta2
    xmax_core = BB[3]+delta2
    ymax_core = BB[4]+delta2
    zmax_core = BB[5]+delta2

    # [Connector/AMR.py specific]
    xmin = BB[0]-2*delta2; ymin = BB[1]-2*delta2; zmin = BB[2]-2*delta2
    xmax = BB[3]+2*delta2; ymax = BB[4]+2*delta2; zmax = BB[5]+2*delta2

    ## Get factor of lengths - cartRX core is rectilinear [Connector/AMR.py specific]
    lenX = xmax_core-xmin_core; lenY = ymax_core-ymin_core; lenZ = zmax_core-zmin_core
    minLen = min(lenX, lenY)
    minLen = min(minLen, lenZ)
    factorX = int(lenX/minLen); factorY = int(lenY/minLen); factorZ = 1
    factorZ = int(lenZ/minLen)
    factorX = min(factorX, 2); factorY = min(factorY, 2); factorZ = min(factorZ, 2) # can cause some wrinkles on the surface

    ni_core = 61; nj_core = 61; nk_core = 61
    hi_core = (xmax_core-xmin_core)/(ni_core-1)
    hj_core = (ymax_core-ymin_core)/(nj_core-1)
    hk_core = (zmax_core-zmin_core)/(nk_core-1)
    h_core = min(hi_core, hj_core)
    h_core = min(h_core, hk_core)
    h_core = min(h_core, 4.*minSnear)

    # Do not extend the CartCore beyond the symmetry plane (symClose)
    if dir_sym > 0:
        if   dir_sym == 1: xmin_core += delta2
        elif dir_sym == 2: ymin_core += delta2
        elif dir_sym == 3: zmin_core += delta2
    # smaller and finer Cartesian core, bigger geometric factor
    XC0 = (xmin_core, ymin_core, zmin_core); XF0 = (xmin, ymin, zmin)
    XC1 = (xmax_core, ymax_core, zmax_core); XF1 = (xmax, ymax, zmax)
    b = G.cartRx3(XC0, XC1, (factorX*h_core, factorY* h_core, factorZ*h_core), XF0, XF1, (1.3, 1.3, 1.3), dim=dim, rank=Cmpi.rank, size=Cmpi.size)

    t0 = time.perf_counter()
    DTW._distance2Walls(b, tbSym, type='ortho', loc='nodes', signed=0)
    tElapse = time.perf_counter()-t0
    tElapse = Cmpi.allreduce(tElapse, op=Cmpi.MAX)
    if Cmpi.master: print("Generate offset frontDP: dist2wall: %.2fs"%tElapse, flush=True)

    C._initVars(b,"cellN",1.)
    # merging of symmetrical bodies in the original blanking bodies
    # required for blankCells as a closed set of surfaces
    bodies = [Internal.getZones(tbSym)]; nbodies = len(bodies)
    BM = numpy.ones((1, nbodies), dtype=numpy.int32)
    t = C.newPyTree(["BASE", Internal.getZones(b)])
    X._blankCells(t, bodies, BM, blankingType='node_in', dim=dim, XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM1)
    C._initVars(t, '{TurbulentDistance}={TurbulentDistance}*({cellN}>0.)-{TurbulentDistance}*({cellN}<1.)')
    iso = P.isoSurfMC(t, 'TurbulentDistance', distOffset)
    iso = Cmpi.allgatherZones(iso)
    iso = C.convertArray2Tetra(iso)
    iso = T.join(iso)
    return iso

# ===============================================================================================================================
def _setIBCData(t, t_exteriorFaces, tb, zibc, frontDP, bbo, IBM_parameters, check, forceAlignment, localDir, different_front_flag):
    ncellsIP = C.getNCells(zibc)

    if ncellsIP > 0:
        Cmpi.trace(" Computing normals via project ortho [start]", master=False, cpu=False)
        _computeIBCNormals__(zibc, tb)
        Cmpi.trace(" Computing normals via project ortho [end]  ", master=False, cpu=False)

        zcibc = C.node2Center(zibc)
        Internal._rmNodesByType(zcibc, "Elements_t")

        Cmpi.trace(" Extracting IBM Points [start]", master=False, cpu=False)
        integrationPts, donorPts, wallPts = getAllIBMPoints__(tb, zcibc, frontDP, bbo, IBM_parameters, check, forceAlignment, localDir)
        Cmpi.trace(" Extracting IBM Points [end]"  , master=False, cpu=False)

        Cmpi.trace(" Adding IBCDatasets [start]", master=False, cpu=False)
        _addIBCData__(t, t_exteriorFaces, donorPts, wallPts, integrationPts, IBM_parameters)
        Cmpi.trace(" Adding IBCDatasets [end]  ", master=False, cpu=False)

    C._rmVars(t, ['cellNFront'])

    if different_front_flag == False: #True is default
        Internal._rmNodesFromName(t, "TurbulentDistance")
        Internal._renameNode(t, "TurbulentDistanceForCFDComputation", "TurbulentDistance")

    Internal._renameNode(t, 'FlowSolution#Centers', 'FlisWallDistance')

    return None

def _setIBCDataDG(t, t_exteriorFaces, tb, zibc, frontDP, bbo, IBM_parameters, check, forceAlignment, localDir, different_front_flag):
    ncellsIP = C.getNCells(zibc)

    if ncellsIP > 0:
        if IBM_parameters["spatial discretization"]["type"] == "FV":
            Cmpi.trace(" Computing normals via project ortho [start]", master=False, cpu=False)
            _computeIBCNormals__(zibc, tb)
            Cmpi.trace(" Computing normals via project ortho [end]  ", master=False, cpu=False)

            zcibc = C.node2Center(zibc)
            Internal._rmNodesByType(zcibc, "Elements_t")
        elif IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
            zcibc = computeSurfaceQuadraturePoints__(t, IBM_parameters, zibc)
            zcibc = computeNormalsForDG__(zcibc, tb)

        Cmpi.trace(" Extracting IBM Points [start]", master=False, cpu=False)
        integrationPts, donorPts, wallPts = getAllIBMPoints__(tb, zcibc, frontDP, bbo, IBM_parameters, check, forceAlignment, localDir)
        Cmpi.trace(" Extracting IBM Points [end]"  , master=False, cpu=False)

        Cmpi.trace(" Adding IBCDatasets [start]", master=False, cpu=False)
        _addIBCData__(t, t_exteriorFaces, donorPts, wallPts, integrationPts, IBM_parameters)
        Cmpi.trace(" Adding IBCDatasets [end]  ", master=False, cpu=False)

    C._rmVars(t, ['cellNFront'])

    if IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
        _computeTurbulentDistanceForDG__(t, tb, IBM_parameters)

    if IBM_parameters["spatial discretization"]["type"] == "FV":
        if different_front_flag == False: #True is default
            Internal._rmNodesFromName(t, "TurbulentDistance")
            Internal._renameNode(t, "TurbulentDistanceForCFDComputation","TurbulentDistance")

    Internal._renameNode(t, 'FlowSolution#Centers', 'FlisWallDistance')

    return None

def _computeIBCNormals__(front, tb):

    varsn = ['gradxTurbulentDistance','gradyTurbulentDistance','gradzTurbulentDistance']
    front_centers = C.node2Center(front)
    proj = T.projectOrtho(front_centers, tb); proj[0] = 'projection'
    x_proj = Internal.getNodeFromName(proj, "CoordinateX")[1]
    y_proj = Internal.getNodeFromName(proj, "CoordinateY")[1]
    z_proj = Internal.getNodeFromName(proj, "CoordinateZ")[1]
    x_front = Internal.getNodeFromName(front_centers, "CoordinateX")[1]
    y_front = Internal.getNodeFromName(front_centers, "CoordinateY")[1]
    z_front = Internal.getNodeFromName(front_centers, "CoordinateZ")[1]
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

def getAllIBMPoints__(tb, zcibc, frontDP, bbo, IBM_parameters, check, forceAlignment=False, localDir='./'):
    projAlgo = 0

    frontTypeDP = IBM_parameters["donor points"]["front type"]
    frontTypeIP = IBM_parameters["integration points"]["front type"]
    IBMType = IBM_parameters["IBM type"]["type"]

    if frontTypeIP == "1":
        depth_IP = IBM_parameters["integration points"]["depth IntegrationPoints"]
    elif frontTypeIP == "2":
        distance_IP = IBM_parameters["integration points"]["distance IntegrationPoints"]

    if frontTypeDP == "2":
        distance_DP = IBM_parameters["donor points"]["distance DonorPoints"]
        C._initVars(zcibc, 'dist', distance_DP)

    integrationPts = C.getAllFields(zcibc, loc='nodes', api=1)[0]
    integrationPts = Converter.convertArray2Node(integrationPts)
    integrationPts = [integrationPts]

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
        res = connector.getIBMPtsWithoutFront(integrationPts, bodies, varsn, 'dist', 1)
        wallPts = res[0]
        donorPts = res[1]
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
        res = connector.getIBMPtsWithFront(integrationPts, listOfSnearsLoc, listOfModelingHeightsLoc, bodies, frontDP, varsn, 1, 2, projAlgo, 0, 0)
        wallPts = res[0]
        donorPts = res[1]

        ## Ouput the IBM points that have a type 3 and type 4 projection
        if len(res) > 3:
            allWallPts = res[0]
            allWallPts = Converter.extractVars(allWallPts, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])

            allInterpPts = res[1]
            allInterpPts = Converter.extractVars(allInterpPts, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])

            allCorrectedPts = Converter.extractVars(integrationPts, ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
            nzonesR         = len(allInterpPts)

            nameZone = ['IBM', 'Wall', 'Donor']
            tLocal3 = C.newPyTree(nameZone)
            tLocal4 = C.newPyTree(nameZone)
            isWrite3 = 0
            isWrite4 = 0
            allProjectPts = res[3]
            allProjectPts = Converter.extractVars(allProjectPts, ['ProjectionType'])
            outputProjection3 = [[],[],[],[],[],[],[],[],[]]
            outputProjection4 = [[],[],[],[],[],[],[],[],[]]
            for noz in range(nzonesR):
                arrayLocal = allProjectPts[noz][1][0]
                type_3 = numpy.count_nonzero(arrayLocal==3)
                type_4 = numpy.count_nonzero(arrayLocal==4)

                if type_3 > 0: X_IBM._prepOutputProject__(outputProjection3, 3, arrayLocal, allCorrectedPts[noz][1], allWallPts[noz][1], allInterpPts[noz][1])
                if type_4 > 0: X_IBM._prepOutputProject__(outputProjection4, 4, arrayLocal, allCorrectedPts[noz][1], allWallPts[noz][1], allInterpPts[noz][1])

            if outputProjection3[0] and check:
                tLocal3  = X_IBM._writeOutputProject__(outputProjection3, tLocal3)
                isWrite3 = 1
            if outputProjection4[0] and check:
                tLocal4  = X_IBM._writeOutputProject__(outputProjection4, tLocal4)
                isWrite4 = 1

            if check:
                print("Rank: %d :: Writing projection files..."%Cmpi.rank, flush=True)
                if isWrite3 > 0:
                    print("Rank: %d :: projection 3 file..."%Cmpi.rank, flush=True)
                    C.convertPyTree2File(tLocal3, localDir+'projection3_Proc_%d.cgns'%Cmpi.rank)
                if isWrite4 > 0:
                    print("Rank: %d :: projection 4 file..."%Cmpi.rank, flush=True)
                    C.convertPyTree2File(tLocal4, localDir+'projection4_Proc_%d.cgns'%Cmpi.rank)
                #if Cmpi.allreduce(isWrite3, op=Cmpi.MAX) > 0: Cmpi.convertPyTree2File(tLocal3, localDir+'projection3.cgns')
                #if Cmpi.allreduce(isWrite4, op=Cmpi.MAX) > 0: Cmpi.convertPyTree2File(tLocal4, localDir+'projection4.cgns')

            del tLocal3
            del tLocal4

        donorPts = projectDPPoints__(integrationPts, donorPts, wallPts, varsn, 1e-8)
    # Check if any of the donor points lays outside the bbox. In this case we modify it.
    if isDPinDomain__(bbo,donorPts)[0] == False:
        print("Rank: %d :: Warning: At least one donor point lays outside the bbox. The point is being moved closer to the wall..."%Cmpi.rank, flush=True)
        list_ids_outside_box = isDPinDomain__(bbo, donorPts)[1]
        epsilon = 0.9
        while (epsilon >= 0.1):
            print("Rank: %d :: Moving the badly located donor point at epsilon %.2f %% of the initial distance from the wall."%(Cmpi.rank, epsilon), flush=True)
            donorPts2correct = copy.deepcopy(donorPts)
            donorPts_modified = projectDPPoints__(integrationPts, donorPts2correct, wallPts, varsn, epsilon, list_ids_outside_box, tb)

            if isDPinDomain__(bbo, donorPts_modified)[0] == False:
                epsilon = epsilon - 0.1
            else:
                donorPts = donorPts_modified
                break
        if isDPinDomain__(bbo, donorPts)[0] == False:
            raise ValueError("Moving the points has not worked. Exiting..")
            Cmpi.abort(errorcode=1)

    wallPts  = Converter.extractVars(wallPts,  ['CoordinateX','CoordinateY','CoordinateZ'])
    donorPts = Converter.extractVars(donorPts, ['CoordinateX','CoordinateY','CoordinateZ'])
    integrationPts   = Converter.extractVars(integrationPts,   ['CoordinateX','CoordinateY','CoordinateZ'])
    array_check = checkMisalignedWallPoints__(integrationPts[0][1], wallPts[0][1], donorPts[0][1], forceAlignment, localDir=localDir)
    if array_check.size != 0 and forceAlignment==True:
        wallPts = projectMisalignedWallPoints__(integrationPts, donorPts, wallPts, array_check, tb, localDir=localDir)
    array_check = checkDPtoIPDistance__(integrationPts[0][1], donorPts[0][1])

    if check:
        print("Rank: %d :: Writing IBM tecplot files..."%Cmpi.rank, flush=True)
        Converter.convertArrays2File(integrationPts  , localDir+"integrationPts_proc%s.plt" %Cmpi.rank)
        Converter.convertArrays2File(wallPts , localDir+"wallPts_proc%s.plt" %Cmpi.rank)
        Converter.convertArrays2File(donorPts, localDir+"donorPts_proc%s.plt" %Cmpi.rank)

    dictOfDonorPtsByIBCName={}
    dictOfIntegrationPtsByIBCName={}
    dictOfWallPtsByIBCName={}
    if (len(res) == 3 and frontTypeDP == "2") or (len(res) == 4 and frontTypeDP == "1"):
        noz = 0 #we always have only one zone
        indicesByTypeForZone = res[2][noz]
        nbTypes = len(indicesByTypeForZone)
        for nob in range(nbTypes):
            ibcTypeL = listOfIBCTypes[nob]
            indicesByTypeL = indicesByTypeForZone[nob]
            if indicesByTypeL.shape[0] > 0:
                ipPtsL = Transform.subzone(integrationPts[noz], indicesByTypeL)
                donorPtsL = Transform.subzone(donorPts[noz], indicesByTypeL)
                wallPtsL = Transform.subzone(wallPts[noz], indicesByTypeL)
            else:
                ipPtsL=[]; donorPtsL = []; wallPtsL = []

            dictOfIntegrationPtsByIBCName[ibcTypeL] = [ipPtsL]
            dictOfWallPtsByIBCName[ibcTypeL] = [wallPtsL]
            dictOfDonorPtsByIBCName[ibcTypeL] = [donorPtsL]
    else:
        raise ValueError("The function connector.getIBMPtsWith/WithoutFront has not worked properly.")
        Cmpi.abort(errorcode=1)
    return dictOfIntegrationPtsByIBCName, dictOfDonorPtsByIBCName,  dictOfWallPtsByIBCName

def projectDPPoints__(integrationPts, donorPts, wallPts, varsn, epsilon, indices_outside_box=None, tb=None):

    nb_donor_pts = donorPts[0][1][0].size
    if indices_outside_box is None:
        dist = epsilon
        for count in range(nb_donor_pts):
            dirx0 = (donorPts[0][1][0][count]-wallPts[0][1][0][count])
            diry0 = (donorPts[0][1][1][count]-wallPts[0][1][1][count])
            dirz0 = (donorPts[0][1][2][count]-wallPts[0][1][2][count])
            dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5
            dist0 = dist/dirn
            donorPts[0][1][0][count] = donorPts[0][1][0][count] + dirx0*dist0
            donorPts[0][1][1][count] = donorPts[0][1][1][count] + diry0*dist0
            donorPts[0][1][2][count] = donorPts[0][1][2][count] + dirz0*dist0
    else:

        zsize = numpy.empty((1,3), E_NpyInt, order='F')
        zsize[0,0] = nb_donor_pts; zsize[0,1] = 0; zsize[0,2] = 0
        zone_donorPts = Internal.newZone(name='DonorPoints', zsize=zsize, ztype='Unstructured')
        gc = Internal.newGridCoordinates(parent=zone_donorPts)
        Internal.newDataArray('CoordinateX', value=donorPts[0][1][0], parent=gc)
        Internal.newDataArray('CoordinateY', value=donorPts[0][1][1], parent=gc)
        Internal.newDataArray('CoordinateZ', value=donorPts[0][1][2], parent=gc)

        DTW._distance2Walls(zone_donorPts, tb, type='ortho', signed=0, dim=3, loc='nodes')
        array_turb_dist = Internal.getNodeFromName(zone_donorPts, "TurbulentDistance")
        for idx in indices_outside_box:
            dist = epsilon * array_turb_dist[1][idx]
            dist = epsilon * array_turb_dist[1][idx]
            dirx0 = (donorPts[0][1][0][idx]-wallPts[0][1][0][idx])
            diry0 = (donorPts[0][1][1][idx]-wallPts[0][1][1][idx])
            dirz0 = (donorPts[0][1][2][idx]-wallPts[0][1][2][idx])
            dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5
            dist0 = dist/dirn
            donorPts[0][1][0][idx] = wallPts[0][1][0][idx] + dirx0*dist0
            donorPts[0][1][1][idx] = wallPts[0][1][1][idx] + diry0*dist0
            donorPts[0][1][2][idx] = wallPts[0][1][2][idx] + dirz0*dist0

    return donorPts

def isDPinDomain__(bbox, coords):

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
    list_ids_outside_box = numpy.concatenate([ids_x_min, ids_y_min, ids_z_min, ids_x_max, ids_y_max, ids_z_max])
    if list_ids_outside_box.shape[0] != 0: out = False

    return out, list_ids_outside_box

def _addIBCData__(t, f, donorPts, wallPts, integrationPts, IBM_parameters):

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
            for nobc,ibc in enumerate(list(integrationPts.values())):
                if ibc!=[[]]:
                    coords_IBC_x = Converter.extractVars(ibc,["CoordinateX"])[0][1][0]
                    coords_IBC_y = Converter.extractVars(ibc,["CoordinateY"])[0][1][0]
                    coords_IBC_z = Converter.extractVars(ibc,["CoordinateZ"])[0][1][0]
                    zibc = Internal.newZone(name="IntegrationPoints",zsize=[[len(coords_IBC_x),0]],ztype="Unstructured")
                    gc = Internal.newGridCoordinates(parent=zibc)
                    Internal.newDataArray('CoordinateX', value=coords_IBC_x, parent=gc)
                    Internal.newDataArray('CoordinateY', value=coords_IBC_y, parent=gc)
                    Internal.newDataArray('CoordinateZ', value=coords_IBC_z, parent=gc)
                    ids = C.identifyNodes(hook, zibc, tol=__TOL__)
                    ids = ids[ids[:] > -1]
                    ids = ids.tolist()
                    ids = [ids[i]-1 for i in range(len(ids))]
                    zf = T.subzone(f,ids, type='elements')
                    G_AMR._addBC2Zone__(z, "IBMWall%d" %nobc, "FamilySpecified:IBMWall", zf)

        #C.freeHook(hook)
        for bc in Internal.getNodesFromType2(z, 'BC_t'):
            famName = Internal.getNodeFromName(bc, 'FamilyName')
            if famName is not None:
                if Internal.getValue(famName) == 'IBMWall':
                    namebc = bc[0]
                    ibcdataset = Internal.createNode('BCDataSet','BCDataSet_t', parent=bc,value='Null')
                    for i in range(N_IP_per_face):
                        dnrPts = Internal.createNode("DonorPointCoordinates"+str(list_suffix_datasets[i]), 'BCData_t', parent=ibcdataset)
                        wallPtsTmp = Internal.createNode("WallPointCoordinates"+str(list_suffix_datasets[i]), 'BCData_t', parent=ibcdataset)

                        coordsPD = Converter.extractVars(donorPts[namebc], ['CoordinateX', 'CoordinateY', 'CoordinateZ'])
                        coordsPW = Converter.extractVars(wallPts[namebc], ['CoordinateX', 'CoordinateY', 'CoordinateZ'])

                        Internal.newDataArray('CoordinateX', value=coordsPD[0][1][0,:][i::N_IP_per_face], parent=dnrPts)
                        Internal.newDataArray('CoordinateY', value=coordsPD[0][1][1,:][i::N_IP_per_face], parent=dnrPts)
                        Internal.newDataArray('CoordinateZ', value=coordsPD[0][1][2,:][i::N_IP_per_face], parent=dnrPts)

                        Internal.newDataArray('CoordinateX', value=coordsPW[0][1][0,:][i::N_IP_per_face], parent=wallPtsTmp)
                        Internal.newDataArray('CoordinateY', value=coordsPW[0][1][1,:][i::N_IP_per_face], parent=wallPtsTmp)
                        Internal.newDataArray('CoordinateZ', value=coordsPW[0][1][2,:][i::N_IP_per_face], parent=wallPtsTmp)

                        if IBM_parameters["spatial discretization"]["type"] in ["DG", "DGSEM"]:
                            coordsPI = Converter.extractVars(integrationPts[bc[0]], ['CoordinateX','CoordinateY','CoordinateZ'])
                            integrationPtsTmp = Internal.createNode("IntegrationPointCoordinates"+str(list_suffix_datasets[i]),'BCData_t', parent=ibcdataset)
                            Internal.newDataArray('CoordinateX', value=coordsPI[0][1][0,:][i::N_IP_per_face], parent=integrationPtsTmp)
                            Internal.newDataArray('CoordinateY', value=coordsPI[0][1][1,:][i::N_IP_per_face], parent=integrationPtsTmp)
                            Internal.newDataArray('CoordinateZ', value=coordsPI[0][1][2,:][i::N_IP_per_face], parent=integrationPtsTmp)

    return None

def checkMisalignedWallPoints__(integrationPts, wallPts, donorPts, forceAlignment=False, localDir="./"):
    x_wall = wallPts[0]
    y_wall = wallPts[1]
    z_wall = wallPts[2]

    x_donor = donorPts[0]
    y_donor = donorPts[1]
    z_donor = donorPts[2]

    x_intp = integrationPts[0]
    y_intp = integrationPts[1]
    z_intp = integrationPts[2]

    wallPointFacePointDist = ( (x_wall - x_intp)**2  + (y_wall - y_intp)**2  + (z_wall - z_intp)**2 )**0.5
    donorPointWallPointDist =( (x_wall - x_donor)**2 + (y_wall - y_donor)**2 + (z_wall - z_donor)**2 )**0.5

    x_facePointWallPointVec = (x_intp - x_wall) / wallPointFacePointDist;
    y_facePointWallPointVec = (y_intp - y_wall) / wallPointFacePointDist;
    z_facePointWallPointVec = (z_intp - z_wall) / wallPointFacePointDist;

    x_offsetPointLocation = x_wall + x_facePointWallPointVec * donorPointWallPointDist;
    y_offsetPointLocation = y_wall + y_facePointWallPointVec * donorPointWallPointDist;
    z_offsetPointLocation = z_wall + z_facePointWallPointVec * donorPointWallPointDist;

    epsDonorPointOffset = 1e-6
    offsetcheck = ( (x_offsetPointLocation - x_donor)**2 + (y_offsetPointLocation - y_donor)**2 + (z_offsetPointLocation - z_donor)**2 )**0.5
    array_check =  numpy.where(offsetcheck > (epsDonorPointOffset * donorPointWallPointDist))[0]
    print("Rank: %d :: Check not aligned points: size array=%d"%(Cmpi.rank,array_check.size), flush=True)
    if array_check.size != 0:
        #if forceAlignment:
        print("Rank: %d :: ATTENTION!!!!!!! Max offset on rank = %g"%(Cmpi.rank, numpy.max(offsetcheck)), flush=True)
        f_wall   = open(localDir+"wall_misaligned_before_proc%s.dat" %Cmpi.rank, "w")
        f_integration = open(localDir+"integration_misaligned_before_proc%s.dat" %Cmpi.rank, "w")
        f_donor  = open(localDir+"donor_misaligned_before_proc%s.dat" %Cmpi.rank, "w")
        for i in range(array_check.size):
            f_wall.write("%f %f %f\n" %(x_wall[array_check[i]], y_wall[array_check[i]], z_wall[array_check[i]]))
            f_integration.write("%f %f %f\n" %(x_intp[array_check[i]], y_intp[array_check[i]], z_intp[array_check[i]]))
            f_donor.write("%f %f %f\n" %(x_donor[array_check[i]], y_donor[array_check[i]], z_donor[array_check[i]]))
        f_wall.close()
        f_integration.close()
        f_donor.close()

        if not forceAlignment:
            raise ValueError("The maximum allowed relative tangential offset was exceeded by one of the donor points. Max offset on rank %d = %g"%(Cmpi.rank, numpy.max(offsetcheck)))
            Cmpi.abort(errorcode=1)
    return array_check

def projectMisalignedWallPoints__(integrationPts, donorPts, wallPts, array_check, tb, localDir='./'):
    nb_donorPts = donorPts[0][1][0].size

    zsize = numpy.empty((1,3), E_NpyInt, order='F')
    zsize[0,0] = nb_donorPts; zsize[0,1] = 0; zsize[0,2] = 0
    zone_integrationPts = Internal.newZone(name='IntegrationPoints', zsize=zsize, ztype='Unstructured')
    gc = Internal.newGridCoordinates(parent=zone_integrationPts)
    Internal.newDataArray('CoordinateX', value=integrationPts[0][1][0], parent=gc)
    Internal.newDataArray('CoordinateY', value=integrationPts[0][1][1], parent=gc)
    Internal.newDataArray('CoordinateZ', value=integrationPts[0][1][2], parent=gc)

    DTW._distance2Walls(zone_integrationPts, tb, type='ortho', signed=0, dim=3, loc='nodes')
    array_turb_dist = Internal.getNodeFromName(zone_integrationPts, "TurbulentDistance")
    f_wall = open(localDir+"wall_misaligned_after_proc%s.dat" %Cmpi.rank,"w")
    for count in array_check:
        dist = array_turb_dist[1][count]
        dirx0 = (donorPts[0][1][0][count]-integrationPts[0][1][0][count])
        diry0 = (donorPts[0][1][1][count]-integrationPts[0][1][1][count])
        dirz0 = (donorPts[0][1][2][count]-integrationPts[0][1][2][count])
        dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5
        dist0 = dist/dirn
        wallPts[0][1][0][count] = integrationPts[0][1][0][count] - dirx0*dist0
        wallPts[0][1][1][count] = integrationPts[0][1][1][count] - diry0*dist0
        wallPts[0][1][2][count] = integrationPts[0][1][2][count] - dirz0*dist0

        f_wall.write("%f %f %f\n" %(wallPts[0][1][0][count], wallPts[0][1][1][count], wallPts[0][1][2][count]))
    f_wall.close()
    return wallPts

def checkDPtoIPDistance__(integrationPts, donorPts):

    x_donor = donorPts[0]
    y_donor = donorPts[1]
    z_donor = donorPts[2]

    x_intp = integrationPts[0]
    y_intp = integrationPts[1]
    z_intp = integrationPts[2]

    donorPointIntegrationPointDist =( (x_intp - x_donor)**2 + (y_intp - y_donor)**2 + (z_intp - z_donor)**2 )**0.5
    epsDonorPointOffset = 1e-9
    array_check =  numpy.where(donorPointIntegrationPointDist < epsDonorPointOffset)[0]
    print("Rank: %d :: Check not aligned points: size array=%d"%(Cmpi.rank, array_check.size), flush=True)

    if array_check.size != 0:
        raise ValueError("WARNING: Coincident integration and donor points on rank %d = " %Cmpi.rank, numpy.min(donorPointIntegrationPointDist))
        Cmpi.abort(errorcode=1)
    return array_check

# ===============================================================================================================================
# DG-related functions
# ===============================================================================================================================
def _addIBC2Zone__(t, f, frontIP):
    for z in Internal.getZones(t):
        hook = C.createHook(f, 'elementCenters')
        ids = C.identifyElements(hook, frontIP, tol=__TOL__)
        ids = ids[ids[:] > -1]
        ids = ids.tolist()
        ids = [ids[i]-1 for i in range(len(ids))]
        #C.freeHook(hook)
        zf = T.subzone(f, ids, type='elements')
        G_AMR._addBC2Zone(z, "IBMWall", "FamilySpecified:IBMWall", zf)
    return None

def computeSurfaceQuadraturePoints__(t, IBM_parameters, frontIP):
    zones = Internal.getZones(t)
    f = P.exteriorFaces(zones[0])
    dims_f = Internal.getZoneDim(f)
    for elt_t in Internal.getNodesFromType(f, "Elements_t"):
        if not elt_t[0].startswith("GridElements"):
            Internal._rmNode(f, elt_t)

    hook = C.createHook(f, 'elementCenters')
    ids = C.identifyElements(hook, frontIP, tol=__TOL__)
    ids = ids[ids[:] > -1]
    ids = ids.tolist()
    ids_IBMWall = [ids[i]-1 for i in range(len(ids))]
    #C.freeHook(hook)
    zf = T.subzone(f, ids_IBMWall, type='elements')
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
        Cmpi.abort(errorcode=1)

    coordsX = Internal.getNodeFromName(t, "CoordinateX")[1]
    coordsY = Internal.getNodeFromName(t, "CoordinateY")[1]
    coordsZ = Internal.getNodeFromName(t, "CoordinateZ")[1]
    elts = Internal.getNodesFromType(t, "Elements_t")
    IBMWall_node = Internal.getNodeFromName(elts, "IBMWall0")
    N_IBM_cells = IBMWall_node[1][1]
    IBMWall_EC = (Internal.getNodeFromName(IBMWall_node, "ElementConnectivity")[1]).reshape(N_IBM_cells, 4)
    cellType = 4
    N_IP_per_face = GetReferencePointsQuad(integrationDegree, quadratureType)[0]
    N_IP = N_IP_per_face * N_IBM_cells
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
        nodalData = numpy.hstack([coordsX[connectivity_IBMFace].reshape(4,1), coordsY[connectivity_IBMFace].reshape(4,1), coordsZ[connectivity_IBMFace].reshape(4,1)])
        quadPoints_surf_location[i:i+N_IP_per_face,:] = interpolationMatrix.dot(nodalData)

    N_IP_surf = N_IBM_cells * N_IP_per_face

    z_IP = Internal.newZone(name="SurfaceIntegrationPoints", zsize=[[N_IP_surf,0]], ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=z_IP)
    Internal.newDataArray('CoordinateX', value=quadPoints_surf_location[:,0], parent=gc)
    Internal.newDataArray('CoordinateY', value=quadPoints_surf_location[:,1], parent=gc)
    Internal.newDataArray('CoordinateZ', value=quadPoints_surf_location[:,2], parent=gc)
    return z_IP

def computeNormalsForDG__(z_IP, tb):
    wall_Points = T.projectOrtho(z_IP, tb)
    x_wallPoints = Internal.getNodeFromName(wall_Points, "CoordinateX")[1]
    y_wallPoints = Internal.getNodeFromName(wall_Points, "CoordinateY")[1]
    z_wallPoints = Internal.getNodeFromName(wall_Points, "CoordinateZ")[1]

    x_integrationPoints = Internal.getNodeFromName(z_IP, "CoordinateX")[1]
    y_integrationPoints = Internal.getNodeFromName(z_IP, "CoordinateY")[1]
    z_integrationPoints = Internal.getNodeFromName(z_IP, "CoordinateZ")[1]

    dirx0 = (x_wallPoints-x_integrationPoints)
    diry0 = (y_wallPoints-y_integrationPoints)
    dirz0 = (z_wallPoints-z_integrationPoints)

    dirn = (dirx0*dirx0+diry0*diry0+dirz0*dirz0)**0.5

    dirx0 = dirx0/dirn
    diry0 = diry0/dirn
    dirz0 = dirz0/dirn
    varsn = ['gradxTurbulentDistance', 'gradyTurbulentDistance', 'gradzTurbulentDistance']
    #varsn = ["nx","ny","nz"]
    FS = Internal.newFlowSolution(name='FlowSolution', gridLocation='Vertex', parent=z_IP)
    Internal.newDataArray(varsn[0], value=dirx0, parent=FS)
    Internal.newDataArray(varsn[1], value=diry0, parent=FS)
    Internal.newDataArray(varsn[2], value=dirz0, parent=FS)
    return z_IP

def _computeTurbulentDistanceForDG__(t, tb, IBM_parameters):

    degree = IBM_parameters["spatial discretization"]["degree"]
    if IBM_parameters["spatial discretization"]["type"] == "DG":
        integrationDegree = 2*degree+1
        quadratureType = "GaussLegendre"
    elif IBM_parameters["spatial discretization"]["type"] == "DGSEM":
        integrationDegree = 2*degree-1
        quadratureType = "GaussLobatto"
    else:
        raise ValueError("Unkown discretization type; options are: FV, DG or DGSEM.")
    coordsX = Internal.getNodeFromName(t, "CoordinateX")[1]
    coordsY = Internal.getNodeFromName(t, "CoordinateY")[1]
    coordsZ = Internal.getNodeFromName(t, "CoordinateZ")[1]
    elts = Internal.getNodesFromType(t, "Elements_t")
    GE_node = Internal.getNodeFromName(elts, "GridElements")
    GE_EC_ravel = Internal.getNodeFromName(GE_node, "ElementConnectivity")[1]
    N_volume_cells = len(GE_EC_ravel)//8
    GE_EC = GE_EC_ravel.reshape(N_volume_cells, 8)

    cellType = 8
    N_IP_per_cell = GetReferencePointsHexa(integrationDegree, quadratureType)[0]
    weights, interpolationMatrix = GetReferencePointsData(integrationDegree, quadratureType, cellType)
    N_IP = N_IP_per_cell * N_volume_cells
    quadPoints_vol_location = numpy.empty((N_IP,3))
    for i,j in zip(range(0,N_IP,N_IP_per_cell), range(N_volume_cells)):
        nodalData = numpy.hstack([coordsX[GE_EC[j]-1].reshape(8,1), coordsY[GE_EC[j]-1].reshape(8,1), coordsZ[GE_EC[j]-1].reshape(8,1)])
        quadPoints_vol_location[i:i+N_IP_per_cell,:] = interpolationMatrix.dot(nodalData)


    N_IP_vol = N_volume_cells * N_IP_per_cell
    z_IP = Internal.newZone(name="VolumeIntegrationPoints", zsize=[[N_IP_vol,0]], ztype="Unstructured")
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
        walldistance_dataset=Internal.newFlowSolution('FlisWallDistance'+str(list_suffix_datasets[i]), parent=zones[0], gridLocation="CellCenter")
        Internal.newDataArray('TurbulentDistance', value=walldistance_volume_ip[i::N_IP_per_cell], parent=walldistance_dataset)

    return None

# ========================================= CURRENTLY NOT USED!! ================================================================
# def computationDistancesNormals(t, tb, dim=3):
#     if dim == 2:
#         dz = 0.01
#         tb2 = T.addkplane(tb)
#         T._contract(tb2, (0,0,0), (1,0,0), (0,1,0), dz)
#     else: tb2 = tb

#     #if Cmpi.rank==0: C.convertPyTree2File(tb2,"tb2.plt")

#     tc = C.node2Center(t)
#     tb_WD = getBodiesDist2wall__(tb2)
#     DTW._distance2Walls(t, tb_WD, type='ortho', signed=0, dim=3, loc='centers')
#     X._applyBCOverlaps(t, depth=2, loc='centers', val=2, cellNName='cellN')
#     C._initVars(t,'{centers:cellNChim}={centers:cellN}')
#     Xmpi._setInterpData(t, tc, nature=1, loc='centers', storage='inverse', sameName=1, sameBase=1, dim=dim, itype='chimera', order=2, cartesian=False)
#     varsn=["gradxTurbulentDistance", 'gradyTurbulentDistance', 'gradzTurbulentDistance']

#     # A COMPARER !!
#     if OPT: t = P.computeGrad(t, 'TurbulentDistance')
#     else: P._computeGrad2(t, 'centers:TurbulentDistance', ghostCells=True, withCellN=False)

#     for v in varsn: C._cpVars(t, 'centers:'+v, tc, v)
#     C._cpVars(t, 'centers:cellNChim', tc, 'cellNChim')
#     Xmpi._setInterpTransfers(t, tc, variables=varsn, cellNVariable='cellNChim', compact=0, type='ID')
#     for v in varsn: t = C.center2Node(t, 'centers:'+v)
#     if not OPT: t = C.center2Node(t, 'centers:TurbulentDistance')
#     return t

# def getMinimumSpacing__(t, dim, snear=1e-1):
#     G._getVolumeMap(t)
#     vol = Internal.getNodeFromName(t, "vol")[1]
#     locsize = (vol/snear)**(1./float(dim))
#     return min(locsize)

# def computeDistance_IP_DP_front42_nonAdaptive__(t, Reynolds, yplus_target, Lref, dim, snear=1e-2):
#     # not used currently - not sure what it does... need to look into it
#     import Geom.IBM as D_IBM
#     distance_IP = D_IBM.computeModelisationHeight(Re=Reynolds, yplus=yplus_target, L=Lref)
#     locsize = getMinimumSpacing__(t, dim, snear)
#     distance_DP = distance_IP+2*(dim**0.5)*locsize
#     distance_DP = min(Cmpi.allgather(distance_DP))
#     return distance_IP, distance_DP
