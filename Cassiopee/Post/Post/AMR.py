#Cassiopee imports
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Intersector.PyTree as XOR
import Connector.IBM as X_IBM
import Converter.Internal as Internal
import Post.IBM as P_IBM
import Post.PyTree as P
import Generator.PyTree as G
import Transform.PyTree as T
import Post.ExtraVariables2 as PE
import os, sys, numpy, math

def createPyTreeForIBMWallFieldsExtraction(donorPointsNumpy, wallPointsNumpy, augStateNumpy, discSelectionParaDict):
    """ Create a CGNS IBM PyTree from IBM numpy data. 
    Usage: createPyTreeForIBMWallFieldsExtraction(donorPointsNumpy, wallPointsNumpy, augStateNumpy, discSelectionParaDict)"""

    ##Serial functions
    print("Creating pyTree of the quantities from the FSMesh of the donorPoints..", flush=True)
    nb_node     = donorPointsNumpy.shape[0]
    utauVector  = numpy.zeros(nb_node)
    yplusVector = numpy.zeros(nb_node)

    zsize = numpy.zeros((1,3), numpy.int32, order='F'); zsize[0,0] = nb_node
    z      = Internal.newZone(name='IBMWall', zsize=zsize, ztype='Unstructured')
    pytree = C.newPyTree(["DonorPointsPyTree",z])
    coordx = numpy.ravel(donorPointsNumpy[:,0])
    coordy = numpy.ravel(donorPointsNumpy[:,1])
    coordz = numpy.ravel(donorPointsNumpy[:,2])

    idxSorted = numpy.lexsort((coordx, coordy, coordz))

    #Create coordinate node
    gc = Internal.newGridCoordinates(parent=z)
    Internal.newDataArray('CoordinateX', value=coordx[idxSorted], parent=gc)
    Internal.newDataArray('CoordinateY', value=coordy[idxSorted], parent=gc)
    Internal.newDataArray('CoordinateZ', value=coordz[idxSorted], parent=gc)

    #Create flow solution node
    pdeList       = ["Euler", "NavierStokes", "RANSSAneg"]
    listVariables = ["X", "Y", "Z"]
    FS = Internal.newFlowSolution(name='FlowSolution', gridLocation='CellCenter', parent=z)
    Internal.newDataArray("Density"           , value=augStateNumpy[:,0][idxSorted], parent=FS)
    Internal.newDataArray("Pressure"          , value=augStateNumpy[:,5+pdeList.index(discSelectionParaDict["PDE"])//2][idxSorted], parent=FS)
    Internal.newDataArray("Temperature"       , value=augStateNumpy[:,6+pdeList.index(discSelectionParaDict["PDE"])//2][idxSorted], parent=FS)
    if discSelectionParaDict["PDE"] != 'Euler':
        Internal.newDataArray("ViscosityMolecular", value=augStateNumpy[:,8+pdeList.index(discSelectionParaDict["PDE"])//2][idxSorted], parent=FS)

    for i, varName in enumerate(listVariables): Internal.newDataArray('Velocity'  +varName      , value=augStateNumpy[:,i+1][idxSorted]/(augStateNumpy[:,0][idxSorted]), parent=FS)
    for i, varName in enumerate(listVariables): Internal.newDataArray('Coordinate'+varName+'_PW', value=numpy.ravel(wallPointsNumpy[:,i])[idxSorted] , parent=FS)
    for i, varName in enumerate(listVariables): Internal.newDataArray('Coordinate'+varName+'_PI', value=numpy.ravel(donorPointsNumpy[:,i])[idxSorted], parent=FS)

    Internal.newDataArray('utau' , value=utauVector , parent=FS)
    Internal.newDataArray('yplus', value=yplusVector, parent=FS)

    return pytree


def extractIBMWallFields(pytree, tb, discSelectionParaDict, ibctype=3):
    """ Extract IBM wall data. Interface between the CODA approach and the FastIBM function.
    Usage: extractIBMWallFields(pytree, tb, discSelectionParaDict, ibctype)"""
    #ibctype = 3 --> Musker
    from numpy.linalg import norm

    print("Extracting Wall Fields..")
    z = Internal.getZones(pytree)[0]

    X_IBM._computeFrictionVelocity(z)

    nameSubRegion = 'IBCD_%d_%s'%(ibctype, z[0])
    FS = Internal.getNodeFromType(z, "FlowSolution_t")

    v = numpy.fromstring(z[0], 'c')
    z[2].append([nameSubRegion, v, [], 'ZoneSubRegion_t'])
    info = z[2][len(z[2])-1]
    info[2].append(FS[2][1:])
    info[2] = info[2][0]

    zw = P_IBM.extractIBMWallFields(z, tb=tb, famZones=[])

    return zw


def computeBoundaryQuantities(zw, dictReferenceQuantities, dim=3, reorderFlag=False, invertYZ=False, verbose=False, time=-1):
    """  Computes the aerodynamic loads at the wall.
    Usage: computeBoundaryQuantities(zw, dictReferenceQuantities, dim, reorderFlag, invertYZ, verbose, time)"""
    if Cmpi.master: print("Computing integral coefficients..")
    if dim == 2:
        zw = C.convertBAR2Struct(zw)
        T._addkplane(zw)
    if reorderFlag == True: T._reorder(zw, (-1,))

    Sref = dictReferenceQuantities["Sref"]
    if invertYZ == False:
        alpha= dictReferenceQuantities["alpha"]
        beta = dictReferenceQuantities["beta"]
    else:
        alpha= dictReferenceQuantities["beta"]
        beta =-dictReferenceQuantities["alpha"]

    zw = C.convertArray2Tetra(zw); zw = G.close(zw)
    zw = C.node2Center(zw, Internal.__FlowSolutionNodes__)
    # add reference state for computation of integrated coefficients
    ref1 = Internal.getNodesFromName(zw, "ReferenceState")
    if ref1 != []: Internal._rmNodesByName(zw, "ReferenceState")
    ref = Internal.newReferenceState("ReferenceState", parent=zw)
    ref[2].append(["VelocityX", dictReferenceQuantities["velX_ref"]    , [], 'DataArray_t'])
    ref[2].append(["VelocityY", dictReferenceQuantities["velY_ref"]    , [], 'DataArray_t'])
    ref[2].append(["VelocityZ", dictReferenceQuantities["velZ_ref"]    , [], 'DataArray_t'])
    ref[2].append(["Density"  , dictReferenceQuantities["density_ref"] , [], 'DataArray_t'])
    ref[2].append(["Pressure" , dictReferenceQuantities["pressure_ref"], [], 'DataArray_t'])
    [res, res2, PressureCoef, FrictionCoef] = P_IBM._loads0(zw, Sref=Sref, alpha=alpha, beta=beta, dimPb=dim, verbose=False)
    ## The division with Cmpi.size is taken 'as is' from the original Post_IBM_CODA.py
    for lst in [res, res2, PressureCoef, FrictionCoef]: lst[:] = [x / Cmpi.size for x in lst]
    [clp, cdp] = PressureCoef
    [clf, cdf] = FrictionCoef
    if verbose and Cmpi.master:
        print("Normalized pressure drag = %.4e and lift = %.4e"%(cdp, clp))
        print("Vector of pressure loads: (Fx_P,Fy_P,Fz_P) = (%.4e, %.4e, %.4e)"%(res[0],res[1],res[2]))

        print("Normalized skin friction drag = %.4e and lift = %.4e"%(cdf, clf))
        print("Vector of skin friction loads: (Fx_f,Fy_f,Fz_f) = (%.4e,%.4e,%.4e)"%(res2[0], res2[1], res2[2]))

        infoTime = ' (time = %.4e)'%time if time >= 0 else ''
        print("******************************************")
        print("Total Drag%s: %.12e"%(infoTime,(cdp+cdf)))
        print("Total Lift%s: %.12e"%(infoTime,(clp+clf)))
        print("******************************************")
    #skin surface, C_(D,total), C_(L,total), C_(D,friction), C_(D,pressure), C_(L,friction), C_(L,pressure)
    return zw, [cdp + cdf, clp + clf, cdf, cdp, clf, clp]


##========================================================================
##========================================================================
## I do not know what the functions below are used for. They are not used in the post scripts obtained from the DLR gitlab.
## Their names suggest they are for the local IBMs but without a test case I do not recommend refactoring them (just yet).

def computeExtraVariablesForLocalIBM(ts, dictReferenceQuantities, dimPb=3):
    print("Computing extra variables..")
    # add reference state for computation of integrated coefficients
    Pref=None; Qref=None;
    if dimPb==2:
        ts = C.convertBAR2Struct(ts)
        T._addkplane(ts)

    ts = C.convertArray2Tetra(ts); ts = G.close(ts)
    ts = C.node2Center(ts, Internal.__FlowSolutionNodes__)

    zones    = Internal.getZones(ts)
    zone_IBM = zones[-1] #last zone is IBM
    if zones:
        ref1 = Internal.getNodesFromName(ts, "ReferenceState")
        if ref1!=[]: Internal._rmNodesByName(ts, "ReferenceState")
        ref = Internal.newReferenceState("ReferenceState", parent=ts)
        ref[2].append(["VelocityX",dictReferenceQuantities["velX_ref"]    , [],'DataArray_t'])
        ref[2].append(["VelocityY",dictReferenceQuantities["velY_ref"]    , [],'DataArray_t'])
        ref[2].append(["VelocityZ",dictReferenceQuantities["velZ_ref"]    , [],'DataArray_t'])
        ref[2].append(["Density"  ,dictReferenceQuantities["density_ref"] , [],'DataArray_t'])
        ref[2].append(["Pressure" ,dictReferenceQuantities["pressure_ref"], [],'DataArray_t'])
        alpha = dictReferenceQuantities["alpha"]; beta = dictReferenceQuantities["beta"]

        RefState = Internal.getNodeFromType(ts, 'ReferenceState_t')
        PInf     = Internal.getValue(Internal.getNodeFromName(RefState, "Pressure"))
        RoInf    = Internal.getValue(Internal.getNodeFromName(RefState, "Density"))
        VxInf    = Internal.getValue(Internal.getNodeFromName(RefState, "VelocityX"))
        VyInf    = Internal.getValue(Internal.getNodeFromName(RefState, "VelocityY"))
        VzInf    = Internal.getValue(Internal.getNodeFromName(RefState, "VelocityZ"))
        VInf2    = VxInf*VxInf+VyInf*VyInf+VzInf*VzInf
        VInf     = math.sqrt(VInf2)
        q        = 0.5*RoInf*VInf2
        qinv     = 1./q
        Qref     = q
        Pref     = PInf
        P_IBM._computeExtraVariables(zone_IBM, Pref, Qref, variables=['Cp','Cf','frictionX','frictionY','frictionZ','frictionMagnitude','ShearStress'])

    return ts


def divideLocalIBMZones(wall_IBM, wall_BF):

    wall_BF = C.convertArray2Tetra(wall_BF)

    CP_BF  = Internal.getNodeFromName(wall_BF, "CoefPressure")[1]
    CFX_BF = Internal.getNodeFromName(wall_BF, "CoefSkinFrictionX")[1]
    CFY_BF = Internal.getNodeFromName(wall_BF, "CoefSkinFrictionY")[1]
    CFZ_BF = Internal.getNodeFromName(wall_BF, "CoefSkinFrictionZ")[1]
    CFT_BF = Internal.getNodeFromName(wall_BF, "CoefSkinFrictionTangential")[1]

    z_IBM = Internal.getZones(wall_IBM)[0]

    hook       = C.createHook(wall_IBM, 'nodes')
    ids        = C.identifyNodes(hook, wall_BF)
    ids_clean  = ids[ids[:] > -1]
    ids_clean  = ids_clean.tolist()
    zf_clean   = T.subzone(wall_IBM, ids_clean, type='nodes')
    hook       = C.createHook(wall_IBM, 'elementCenters')
    ids        = C.identifyElements(hook, zf_clean)
    ids_clean  = ids[ids[:] > -1]-1
    ids_clean  = ids_clean.tolist()
    zf_clean_2 = T.subzone(wall_IBM, ids_clean, type='elements')


    list_all_ids_IBM = [i for i in range(z_IBM[1][0][1])]
    ids_ibm_zone     = numpy.setdiff1d(list_all_ids_IBM, ids_clean).tolist()
    zf_ibm_zone      = T.subzone(wall_IBM, ids_ibm_zone, type='elements')

    # paste solution on boundaries (Cp/Cf) with new indices
    hook = C.createHook(wall_BF, 'nodes')
    ids_clean_zf2 = C.identifyNodes(hook, zf_clean_2)
    ids_clean_zf2 = ids_clean_zf2[ids_clean_zf2[:] > -1]-1
    ids_clean_zf2 = ids_clean_zf2.tolist()
    FS_BF = Internal.getNodeFromName(zf_clean_2,"FlowSolution")
    FS_BF[2].append(["Cp"       , CP_BF[ids_clean_zf2] , [], 'DataArray_t'])
    FS_BF[2].append(["frictionX", CFX_BF[ids_clean_zf2], [], 'DataArray_t'])
    FS_BF[2].append(["frictionY", CFY_BF[ids_clean_zf2], [], 'DataArray_t'])
    FS_BF[2].append(["frictionZ", CFZ_BF[ids_clean_zf2], [], 'DataArray_t'])
    FS_BF[2].append(["Cf"       , CFT_BF[ids_clean_zf2], [], 'DataArray_t'])

    t = C.newPyTree(["New PyTree", [Internal.getZones(zf_clean_2)[0], Internal.getZones(zf_ibm_zone)[0]]])
    return t


def integCfn(teff):
    """Integ tau.n.ds"""
    import Post.Mpi as Pmpi
    retx = Pmpi.integ(teff, 'centers:frictionX')
    rety = Pmpi.integ(teff, 'centers:frictionY')
    retz = Pmpi.integ(teff, 'centers:frictionZ')
    return [retx[0], rety[0], retz[0]]


def _loads0LocalIBM(ts, Sref=1.0, Qref=1.0, alpha=0.0, beta=0.0, dimPb=2, BFtreatment="nodal"):
    dimPb = 2
    T._reorder(ts, (-1,))

    z_BF  = Internal.getZones(ts)[0]
    z_IBM = Internal.getZones(ts)[1]


    res_pres_BF = PE.integCp(z_BF)[0]
    res_pres_BF = [-i/Sref for i in res_pres_BF]
    res_fric_BF = integCfn(z_BF)
    res_fric_BF = [ i/Sref for i in res_fric_BF]

    res_pres_IBM = PE.integCp(z_IBM)[0]
    res_pres_IBM = [-i/Sref for i in res_pres_IBM]
    res_fric_IBM = PE.integTaun(z_IBM)
    res_fric_IBM = [ i/Sref for i in res_fric_IBM]

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

    FSC = Internal.getNodesFromName(ts, Internal.__FlowSolutionCenters__)
    Internal._rmNodesFromName(FSC, 'ShearStress*')

    return ts, [[cd, cl, cdf, cdp, clf, clp], [cd_BF, cl_BF, cdf_BF, cdp_BF, clf_BF, clp_BF], [cd_IBM, cl_IBM, cdf_IBM, cdp_IBM, clf_IBM, clp_IBM]]


def intersectObstacleMesh(t, surf_detail,list_bcs):
    print("Performing intersection for Local IBM approach..")

    zones = Internal.getZones(t)
    z     = zones[0]

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

    xzs    = Internal.getZones(x)
    result = T.join(xzs[0], xzs[2]) # warning : cas-dependant !

    return result


def extractNotIBMWallFields(t, IBM_parameters, list_boundary_values_to_extract, discSelectionParaDict, BFtreatment="nodal"):
    print("Extracting wall fields that are not IBM..")
    C._deleteFlowSolutions__(t)

    zones = Internal.getZones(t)
    z     = zones[0]

    zoneDim = Internal.getZoneDim(z)
    print(zoneDim[0])

    for i,bc in enumerate(Internal.getNodesFromType(z, 'BC_t')):
        bctype = Internal.getValue(bc)
        if bctype.startswith("BCWall"):
            zbc = C.extractBCOfType(t, bctype)
            break
    zbc = T.join(zbc)
    FS  = Internal.getNodesFromName(zbc, "FlowSolution#Centers")
    Internal._rmNodesFromName(FS, "CoefSkinFrictionImmersed")

    names_var = list_boundary_values_to_extract
    FS_1 = []
    for i in range(len(names_var)):
        FS_ = numpy.empty((0,3))
        for j in range(len(FS)):
            FS_ = numpy.append(FS_, FS[j][2][i+1][1])
        FS_1.append(FS_)

    # only one
    pytree_new = C.newPyTree(["PyTree"])
    indicesF=[]; f = P.exteriorFaces(z, indices=indicesF)
    hook = C.createHook(f, 'elementCenters')

    zbc = T.join(zbc)
    ids = C.identifyElements(hook, zbc)
    ids = ids[ids[:] > -1]
    ids = (ids-1).tolist()
    zf  = T.subzone(f, ids, type='elements')

    list_skin_friction_quantities = ["CoefSkinFrictionX","CoefSkinFrictionY","CoefSkinFrictionZ","CoefSkinFrictionTangential"]
    if discSelectionParaDict["PDE"]=="Euler": names_var = names_var+list_skin_friction_quantities
    zf[2] = zf[2][:3]
    GE = Internal.getNodeFromName(zf, "GridElements_HEXA")
    for idx,name_var in enumerate(names_var):
        zf = C.initVars(zf, 'centers:%s' %name_var,0.)
        if discSelectionParaDict["PDE"]!="Euler" or idx==0:
            node_PCFC = Internal.getNodeFromName(zf, name_var)
            node_PCFC[1] = FS_1[idx]

    zf_nodes = C.center2Node(zf, "FlowSolution#Centers")
    Internal._rmNodesByName(zf_nodes, "FlowSolution#Centers")

    return zf_nodes
