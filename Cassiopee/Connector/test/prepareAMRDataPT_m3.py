# kpython -n 12 -t 4 test.py
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Converter.Internal as Internal
import Generator.PyTree as G
import Generator.Mpi as Gmpi

import Connector.AMR as X_AMR
import Generator.IBMmodelHeight as G_IBM_Height
import Transform.PyTree as T
import KCore.test as test

__TOL__ = 1e-9
LOCAL = test.getLocal()+'/'

# ===============================================================================================================================
def prepareAMRDataFV__(t_case, t, IBM_parameters=None, check=False, dim=3, localDir='./', forceAlignment=False, isFastApproach=True):
    VPM = False

    (IBM_parameters, frontTypeIP, frontTypeDP, dir_sym, different_front_flag) = X_AMR.checkInputsIbmParam__(IBM_parameters)

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
    X_AMR._dist2wallIBM(t, tb2, dim, different_front_flag)
    if Cmpi.master: test.testT(t,1)

    #==========================================================
    # STEP 1: Blanking by immersed body
    #==========================================================
    Cmpi.trace(">>> Blanking [start]", master=True, cpu=False)
    t = X_AMR.blankingIBM(t, tb2)
    Cmpi.trace(">>> Blanking [end]  ", master=True, cpu=False)
    if Cmpi.master: test.testT(t,10)

    #=======================================================================
    # STEP 2: Save BC Names & Types - check if QuadNQuad is fully inside IBM
    #=======================================================================
    BCInfo = X_AMR.getBCs(t, tb2, dim)
    if Cmpi.master: test.testO(BCInfo,20)

    #===============================================================================================================================
    # STEP 3: Get Integration Point Front (Recall: IP = integration point (NOT image point) - Target cells in Mittal et al. approach)
    #===============================================================================================================================
    # only done for frontTypeIP=2 - if frontTypeIP=1 it is F1 so blankByIBCBodies is sufficient for the cellN value
    Cmpi.trace("Extract front faces of IBM integration points [start] ", master=True, cpu=False)
    frontIP, maxDistIP = X_AMR.getFrontIP(t, IBM_parameters, localDir, check=check, VPM=VPM)
    Cmpi.trace("Extract front faces of IBM integration points [end]   ", master=True, cpu=False)
    if Cmpi.master: test.testT(frontIP, 30)

    #====================================================
    # STEP 4: Select Cells - only fluid cells from now on
    #====================================================
    Cmpi.trace(" Removing blanked cells [start]", master=True, cpu=False)
    t = X_AMR.removeBlankedCells(t)
    Cmpi.trace(" Removing blanked cells [end]  ", master=True, cpu=False)
    if Cmpi.master: test.testT(t, 40)
    
    #===============================
    # STEP 5: Recover BCs - Add IBCs
    #===============================
    Cmpi.trace(" Recovering Boundary Conditions [start]", master=True, cpu=False)
    t, zibc = X_AMR.recoverBCs(t, frontIP, BCInfo)
    Cmpi.trace(" Recovering Boundary Conditions [end]  ", master=True, cpu=False)
    if Cmpi.master:
        test.testT(t, 50)
        test.testT(zibc, 52)

    #===================================================================================
    # STEP 6: Get Donor Point Front (Recall: DP = image point in Mittal et. al approach)
    #         Determine DP points on DPfront
    #         Add IBC Dataset to the PyTree
    #===================================================================================
    Cmpi.trace(" Extracting front of the donor points [start]", master=True, cpu=False)
    frontDP = None
    if frontTypeDP == '1':
        frontDP = X_AMR.getFrontDP(t, tb2, frontIP, dim, dir_sym, check, distIP=maxDistIP, localDir=localDir, isFastApproach=isFastApproach)
    Cmpi.trace(" Extracting front of the donor points [end]  ", master=True, cpu=False)
    if Cmpi.master: test.testT(frontDP, 60)

    #===================================================================================
    # STEP 7: Determine DP points on DPfront & Add IBC Dataset to the PyTree
    #===================================================================================
    Cmpi.trace(" Getting DP IBM points and adding IBC Dataset [start]", master=True, cpu=False)
    X_AMR._setIBCData(t, frontIP, tb2, zibc, frontDP, bbo, IBM_parameters, check, forceAlignment, localDir, different_front_flag)
    Cmpi.trace(" Getting DP IBM points and adding IBC Dataset [end]  ", master=True, cpu=False)
    if Cmpi.master: test.testT(t, 70)

    return t

## =====================================================================================================
t = Cmpi.convertFile2PyTree('tAMR2DMpi.cgns', proc=Cmpi.rank)
tb = C.convertFile2PyTree('tb_NACA0012_singleSnear.cgns')

##Mesh
snearReal  = 1.e-04
snear      = 10.0*snearReal #target:: 1.0e-4
dFar       = 15
levelMax   = 20
vmins      = [25,24,24,20,21,5,4]
dim        = 2

##IBM
frontType    = 2
Reynolds     = 6.0e06
yplus_target = 100 #100*2.5
Lref         = 1.0
Aref         = 1.0
## =====================================================================================================


depth_IP     = 1
distance_IP  = G_IBM_Height.computeModelisationHeight(Re=Reynolds, yplus=yplus_target, L=Lref)
if Cmpi.master and frontType==2: print("distance_IP=", distance_IP,flush=True)
################### IBM pre-processing #########################
if frontType==1:
    integrationPoints= {
        "front type": "1",
        "depth IntegrationPoints": depth_IP,  # for type 1
    }
else:
    integrationPoints= {
        "front type": "2",
        "distance IntegrationPoints": distance_IP,  # for type 2
    }

IBM_parameters = {
    "donor points": {
        "front type": "1",
        "depth DonorPoints": 1,  # for type 1
        # "distance DonorPoints": distance_DP,  # for type 2
    },
    "integration points": integrationPoints,
    "IBM type": {
        "type": "global",
        "size elements body": snear,
    },
    "spatial discretization": {
        "type": "FV",
    },
}

t = prepareAMRDataFV__(tb, t, IBM_parameters=IBM_parameters, dim=dim, check=False, localDir=LOCAL)

Cmpi.convertPyTree2File(t, 'check_t.cgns')
