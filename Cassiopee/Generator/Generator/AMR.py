# IMPROVEMENTS
# cartRX in offsets ?
# snear to adjust as in FAST/IBM
import Converter.PyTree as C
import Converter.Mpi as Cmpi
import Transform.PyTree as T
import Converter.Internal as Internal
import Connector.PyTree as X
import Dist2Walls.PyTree as DTW
import Post.PyTree as P
import XCore.PyTree as XC
from . import PyTree as G
import os, numpy, math, time

__TOL__ = 1e-9

def vminsInputCheck__(vminsIN, numBaseTMP, levelMaxTMP):
    import copy
    vminsTMP = copy.deepcopy(vminsIN)
    # list of vminsTMP
    if isinstance(vminsTMP,list):
        # vmin = [] & vmin =[[],[]] for 3 bases
        if not isinstance(vminsTMP[0],list): vminsTMP = [vminsTMP] # vmin = [] --> vmin = [[]]

        # vmin = [[],[]] --> vmin = [[],[],[]] for 3 bases
        while len(vminsTMP) < numBaseTMP: vminsTMP.append(vminsTMP[-1])
    else:
        # vmin = 10
        vminsTMP = numpy.ones((numBaseTMP,levelMaxTMP))*vminsTMP
        vminsTMP = vminsTMP.tolist() # needed
    return vminsTMP

# Important note: This function (__addItemDict) must be moved into adaptMesh__ if only used in that function.
#                 Leaving it here temporarily - no certainty that it wont be used elsewhere.
#                 If by June 2026 it is not used elsewhere make definite move to adaptMesh__
def _addItemDict__(d, key, value):
    if key in d:
        d[key].append(value)
    else:
        d[key] = [value]
    return None

def getListSnear__(tb,snears):
    # List of snears
    if not isinstance(snears, list):
        snearsList = []
        for baseTmp in Internal.getBases(tb):
            bodies    = Internal.getZones(baseTmp)
            snearsTmp = [snears*1.]*len(bodies)
            for c, z in enumerate(bodies):
                n = Internal.getNodeFromName2(z, 'snear')
                if n is not None: snearsTmp[c] = Internal.getValue(n)*1.
            snearsList.append(snearsTmp)
    else:
        bodies = Internal.getZones(tb)
        if len(bodies) != len(snears): raise ValueError('generateAMRMesh (generateSkeletonMesh__): Number of bodies is not equal to the size of snears.')
    snears = snearsList
    tbTMP = Internal.copyTree(tb)
    baseSYM = Internal.getNodesFromName1(tbTMP, "SYM")
    if baseSYM: tbTMP = Internal.rmNodesByNameAndType(tbTMP, 'SYM', 'CGNSBase_t')
    numBase = len(Internal.getBases(tbTMP))
    return snears, numBase

def checkBaseNames__(tb,tbox):
    # check consistency of base names in tb and tbox
    # check that tb does not have base names that repeat
    # check that the tb and tbox do not have the same names
    baseNameListLocal = []
    for b in Internal.getBases(tb):
        if b[0] in baseNameListLocal:
            raise ValueError('Tb must have bases with distinct names. Exiting...')
            Cmpi.abort(errorcode=1)
        else:
            baseNameListLocal.append(b[0])
    if tbox is not None:
        for b in Internal.getBases(tbox):
            if b[0] in baseNameListLocal:
                raise ValueError('Tb and Tbox have a base with the same name. Tbox must have base names that are distinct to those found in tb. Exiting...')
                Cmpi.abort(errorcode=1)
            else:
                baseNameListLocal.append(b[0])
    return None

# Generation of the list of offset surfaces starting from tb
# IN: offsetValues : list of float values defining the offset distance to tb
# if opt: mmgs is used to coarsen the tb surface to optimize distance field
# returns a tree toffset
def generateListOfOffsets__(tb, snears, offsetValues=[], dim=3, opt=False, numTbox=0, tbv2=None):
    import Geom.IBM as D_IBM
    if offsetValues == []: return []
    numBase = len(Internal.getBases(tb))
    if Cmpi.master: print('Generating list of offsets...start', flush=True)

    minSnear=1e10
    for snearLocal in snears: minSnear = min(minSnear, snearLocal[0])

    dir_sym = getSymmetryPlaneInfo__(tb, dim=dim)
    baseSYM = Internal.getNodesFromName1(tb, "SYM")
    listShiftBase = []
    if baseSYM is not None:
        # Keep only the real closed tb
        # Remove SYM Base & Zones
        tb = Internal.rmNodesByNameAndType(tb, 'SYM', 'CGNSBase_t')
        for b in Internal.getBases(tb):
            node = Internal.getNodesFromNameAndType(b, '*_sym', 'Zone_t')
            if node: listShiftBase.append(b[0])
        tb = Internal.rmNodesByNameAndType(tb, '*_sym*', 'Zone_t')
    numBase = len(Internal.getBases(tb))
    if opt and dim == 3:
        for nob in range(len(tb[2])):
            if Internal.getType(tb[2][nob]) == 'CGNSBase_t':
                z = Internal.getZones(tb[2][nob])
                z = C.convertArray2Tetra(z)
                z = T.join(z)
                bbz = G.bbox(z)
                # [TODO] this needs to be related to the hmin
                # hausd is a length and must be adapted to the dimensions of each case
                hausd = max(bbz[3]-bbz[0], bbz[4]-bbz[1], bbz[5]-bbz[2])/10000.
                hmax = hausd*1000
                if Cmpi.master:
                    print('Remeshing (tb) surface mesh --> Maximum chordal deviation between final and initial mesh::%g || Maximum mesh step in final mesh::%g'%(hausd, hmax), flush=True)
                # exteriorFaces currently crashes if the surface is closed
                try:
                    fixedConstraints = P.exteriorFaces(z)
                except:
                    fixedConstraints = []
                z = G.mmgs(z, hausd=hausd, hmax=hmax, fixedConstraints=fixedConstraints)
                tb[2][nob][2] = Internal.getZones(z)

        #[Temp. Patch] - This need to be generalized.
        # tbv2: tb version 2 - it is tb with an small extrusion outward of the calc. domain to guarantee a correct offset for the near body offsets.
        # This is a (very) temporary fix & is intended to be replaced with a more suitable approach after extensive testing.
        # Current use: 3D symmetric cases & needed for the CRM case 1 of the HLPW5.
        if tbv2 is not None:
            for nob in range(len(tbv2[2])):
                if Internal.getType(tbv2[2][nob]) == 'CGNSBase_t':
                    z = Internal.getZones(tbv2[2][nob])
                    z = C.convertArray2Tetra(z)
                    z = T.join(z)
                    bbz = G.bbox(z)
                    # hausd is a length and must be adapted to the dimensions of each case
                    hausd = max(bbz[3]-bbz[0], bbz[4]-bbz[1], bbz[5]-bbz[2])/10000.
                    hmax = hausd*1000
                    if Cmpi.master:
                        print('Remeshing aux. temp. (tbv2) surface mesh --> Maximum chordal deviation between final and initial mesh::%g || Maximum mesh step in final mesh::%g'%(hausd, hmax), flush=True)
                    # exteriorFaces currently crashes if the surface is closed
                    try:
                        fixedConstraints = P.exteriorFaces(z)
                    except:
                        fixedConstraints = []
                    z = G.mmgs(z, hausd=hausd, hmax=hmax, fixedConstraints=fixedConstraints)
                    tbv2[2][nob][2] = Internal.getZones(z)

    tbTmp = Internal.copyTree(tb)
    if tbv2 is not None: tbTmp = tbv2

    t_offset = C.newPyTree()
    no_offsetGlobal = 0
    for nBase, tbLocal in enumerate(Internal.getBases(tb)):
        BB = G.bbox(tbLocal)
        ni = 150; nj = 150; nk = 150
        XRAYDIM1 = 3*ni; XRAYDIM2 = 3*nj

        offsetValMin = min(offsetValues[nBase])
        offsetValMax = max(offsetValues[nBase])
        alpha = 1.1
        delta = alpha*offsetValMax
        xmin = BB[0]-delta; ymin = BB[1]-delta; zmin = BB[2]-delta
        xmax = BB[3]+delta; ymax = BB[4]+delta; zmax = BB[5]+delta

        # CARTRX
        # delta2 = sum(offsetValues)/len(offsetValues)n
        delta2 = max(BB[3]-BB[0], BB[4]-BB[1], BB[5]-BB[2])*0.02 # 2% seems enough for the external cases already tested
        xmin_core = BB[0]-delta2
        ymin_core = BB[1]-delta2
        zmin_core = BB[2]-delta2
        xmax_core = BB[3]+delta2
        ymax_core = BB[4]+delta2
        zmax_core = BB[5]+delta2

        if dim == 2: ni_core = 101; nj_core = 101; nk_core = 101
        else: ni_core = 61; nj_core = 61; nk_core = 61
        hi_core = (xmax_core-xmin_core)/(ni_core-1)
        hj_core = (ymax_core-ymin_core)/(nj_core-1)
        hk_core = (zmax_core-zmin_core)/(nk_core-1)
        h_core = min(hi_core, hj_core)
        if dim == 3:
            h_core = min(h_core, hk_core)
            h_core = min(h_core, 8.*minSnear)
        if dim == 2:
            zmin = 0; zmax = 0
            zmin_core = 0.; zmax_core = 0.
            hk_core = 0.
            # Pull request note: h_core may cause regressions in the mesh generation
            h_core = min(h_core, 16.*minSnear)

        # Do not extend the CartCore beyond the symmetry plane (symClose)
        if dir_sym > 0 and tbLocal[0] in listShiftBase:
            if   dir_sym == 1: xmin_core += delta2
            elif dir_sym == 2: ymin_core += delta2
            elif dir_sym == 3: zmin_core += delta2
        # smaller and finer Cartesian core, bigger geometric factor
        XC0 = (xmin_core, ymin_core, zmin_core); XF0 = (xmin, ymin, zmin)
        XC1 = (xmax_core, ymax_core, zmax_core); XF1 = (xmax, ymax, zmax)
        b = G.cartRx3(XC0, XC1, (h_core, h_core, h_core), XF0, XF1, (1.3, 1.3, 1.3), dim=dim, rank=Cmpi.rank, size=Cmpi.size)

        tbLocalTmp = Internal.getNodeFromNameAndType(tbTmp, tbLocal[0], 'CGNSBase_t')
        if tbLocalTmp is None: tbLocalTmp = tbLocal
        C._initVars(tbLocalTmp, "cellN", 1.)

        t0 = time.perf_counter()
        DTW._distance2Walls(b, tbLocalTmp, type='ortho', loc='nodes', signed=0)
        tElapse = time.perf_counter()-t0
        tElapse = Cmpi.allreduce(tElapse, op=Cmpi.MAX)
        if Cmpi.master: print("Generate list of offsets: Base %s Num. %d:dist2wall: %.2fs"%(tbLocal[0], nBase, tElapse), flush=True)

        C._initVars(b,"cellN",1.)
        # merging of symmetrical bodies in the original blanking bodies
        # required for blankCells as a closed set of surfaces
        bodies = [tbLocalTmp]; nbodies = 1
        BM = numpy.ones((1, nbodies), dtype=numpy.int32)
        t = C.newPyTree(["BASE", Internal.getZones(b)])
        X._blankCells(t, bodies, BM, blankingType='node_in', dim=dim, XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2)
        C._initVars(t, '{TurbulentDistance}={TurbulentDistance}*({cellN}>0.)-{TurbulentDistance}*({cellN}<1.)')
        #Cmpi.convertPyTree2File(b, 'meshForOffsetBase%d.cgns'%nBase) # DEBUG ONLY

        # all body offsets are prefaced by 'z_offsetBase' - only the zone name
        # all tbox offsets are prefaced by 'Tbox_offsetBase' - only the zone name
        preffixLocal = 'z_offsetBase'
        if nBase >= numBase-numTbox: preffixLocal = 'Tbox_offsetBase'
        for no_offset, offsetval in enumerate(offsetValues[nBase]):
            if Cmpi.master: print("Offset %d - value: %g - snear: %g"%(no_offset, offsetval, snears[nBase][0]*2**no_offset), flush=True)
            iso = P.isoSurfMC(t, 'TurbulentDistance', offsetval)
            iso = Cmpi.allgatherZones(iso)
            iso = C.convertArray2Tetra(iso)
            iso = T.join(iso)
            iso = G.close(iso, tol=1.e-6)
            iso = T.smooth(iso)
            iso[0]='%s%d_%d'%(preffixLocal, nBase, no_offset)
            D_IBM._setSnear(iso, snears[nBase][0]*2**no_offset)
            C._addBase2PyTree(t_offset, 'OFFSETBase%d_%d'%(nBase, no_offset))
            t_offset[2][no_offsetGlobal+1][2]=[iso]
            no_offsetGlobal += 1
    return t_offset

# Generates an isotropic skeleton mesh to be adapted then by AMR
def generateSkeletonMesh__(tb, snears, dfars=10., dim=3, levelSkel=7, octreeMode=0):
    surfaces=[]; dfarList=[]; snearsList=[]; levelSkelList=[]
    # This clips the upper limit on the number of offset level to the input value. Important for tests & devs.
    # This is an expert expert parameter & should be used with (a lot of) caution.
    # This is needed to bypass G.adaptOctree that can be very expensive when we need a fine background (outside the offset levels) grid.
    forceUpperLimitOffset = False

    # list of dfars
    bodies = Internal.getZones(tb)
    if not isinstance(dfars, list):
        dfars = [dfars*1.]*len(bodies)
        for c, z in enumerate(bodies):
            n = Internal.getNodeFromName2(z, 'dfar')
            if n is not None: dfars[c] = Internal.getValue(n)*1.
    else:
        if len(bodies) != len(dfars): raise ValueError('generateAMRMesh (generateSkeletonMesh__): Number of bodies is not equal to the size of dfars.')
    levelSkelInput = levelSkel
    for c, z in enumerate(bodies):
        levelSkel = levelSkelInput
        if dfars[c] > -1: #body snear is only considered if dfar_loc > -1
            surfaces.append(z)
            # Pull request note: levelSkelLoc causes regressions in the mesh generation
            levelSkelLoc = int(math.log2(dfars[c]/snears[c])) # as the dfar is fixed we do not need a fraction of the dfar to get the levelSkelLoc
            #levelSkelLoc = int(math.log2(0.2*dfars[c]/snears[c])) # Old levelSkelLoc. Stays here in case it is needed in the future
            if not forceUpperLimitOffset: levelSkel = max(levelSkel, levelSkelLoc) # security so that levelSkel is not too small
            dfarloc = dfars[c]
            snearloc = 2**levelSkel*snears[c]
            while snearloc > dfarloc/2: # security so that levelSkel is not too big
                snearloc  /= 2.
                levelSkel -= 1
            levelSkelList.append(levelSkel)
            snearsList.append(snearloc)
            dfarList.append(dfarloc)

    o = G.octree(surfaces, snearList=snearsList, dfarList=dfarList, balancing=1, octreeMode=octreeMode)
    levelSkel = max(levelSkelList)
    #
    # SYMMETRY - select only cells from one side
    #
    # determine where the symmetry plane is
    dir_sym = getSymmetryPlaneInfo__(tb, dim=dim)
    [xmin, ymin, zmin, xmax, ymax, zmax] = G.bbox(o)
    if dir_sym == 1:
        coordsym = 'CoordinateX'
        valsym = 0.5*(xmin+xmax)
    elif dir_sym == 2:
        coordsym = 'CoordinateY'
        valsym = 0.5*(ymin+ymax)
    elif dir_sym == 3:
        coordsym = 'CoordinateZ'
        valsym = 0.5*(zmin+zmax)
    if dir_sym > 0: o = P.selectCells(o,'{%s}>%g'%(coordsym, valsym-__TOL__), strict=1)

    if forceUpperLimitOffset:
        box = G.bbox(o)
        dxdydz_local      = min(snearsList)
        nCellsCartesian_x = int((box[3]-box[0])/dxdydz_local)
        nCellsCartesian_y = int((box[4]-box[1])/dxdydz_local)
        nCellsCartesian_z = 1
        dz_local          = box[5]-box[2]
        if dim == 3:
            nCellsCartesian_z = int((box[5]-box[2])/dxdydz_local)
            dz_local          = dxdydz_local
        o = G.cart((box[0],box[1],box[2]), (dxdydz_local, dxdydz_local, dz_local), (nCellsCartesian_x+1, nCellsCartesian_y+1,nCellsCartesian_z+1))
    else:
        # adapt the mesh to get a single refinement level - uniform grid
        refined = True
        G._getVolumeMap(o)
        volminAll = C.getMinValue(o, "centers:vol")
        tol_vol = 1e-2*volminAll
        while refined:
            C._initVars(o, '{centers:indicator}=({centers:vol}>%g)'%(volminAll+tol_vol))
            if C.getMaxValue(o, "centers:indicator") > 1.-__TOL__:
                o = G.adaptOctree(o, 'centers:indicator', balancing=1)
                G._getVolumeMap(o)
            else:
                break
        C._rmVars(o, ['centers:indicator', 'centers:vol'])

    while(C.getNCells(o) < Cmpi.size):
        # when the background grid has less # of cells than the Cmpi.size the load and split during the adaptMesh
        # yields procs with 0 cells. This is a safeguard to avoid this.
        if Cmpi.master:
            print("=======================================================================", flush=True)
            print("========================       WARNING!        ========================", flush=True)
            print("======================== Background Mesh Mod.! ========================", flush=True)
            print("======= The background has been divided due to Nmpi > Ncells!!  =======", flush=True)
            print("=== The background has been divided due to Nmpi [%d]> Ncells [%d]!! ==="%(Cmpi.size, C.getNCells(o)), flush=True)
            print("=======================================================================", flush=True)
        G._getVolumeMap(o)
        volminAll = C.getMinValue(o, "centers:vol")/2.
        tol_vol = 1e-2*volminAll
        C._initVars(o, '{centers:indicator}=({centers:vol}>%g)'%(volminAll+tol_vol))
        o = G.adaptOctree(o, 'centers:indicator', balancing=1)
        levelSkel -= 1

    if dim == 2: T._addkplane(o)
    o = C.convertArray2NGon(o)
    _addPhysicalBCs__(o, tb, dim=dim)
    Internal._adaptNGon32NGon4(o)
    if Cmpi.master: print('Generating skeleton mesh...end', flush=True)
    return o, levelSkel

def tagOutsideBody__(o, tbTMP, dim=3, h_target=-1., opt=False, noffsets=None, coarseXray=False):
    # To avoid adapting inside the bodies when the bodies and the tbox intersect we have this function.
    # It tags the inside of the bodies as cellN=0 and then multuiplies the indicator. i.e. the parts inside the body will be zero.

    to = C.newPyTree(["OCTREE"]); to[2][1][2]=Internal.getZones(o)
    C._initVars(to, 'centers:indicatorTmp', 0.)

    if dim == 2: tbTMP = T.addkplane(tbTMP)
    bodies1 = [Internal.getZones(tbTMP)]
    bb1 = G.bbox(tbTMP)
    L1 = max(bb1[3]-bb1[0], bb1[4]-bb1[1])
    if dim == 3: L1 = max(L1, bb1[5]-bb1[2])

    BM = numpy.ones((1,1), dtype=Internal.E_NpyInt)

    # ideally we should use blankCellsTri to avoid XRAYDIM but currently not safe
    XRAYDIM1 = int(L1/h_target)+10;
    # [Temp. Patch] - This need to be generalized. Works for all test cases currently considered but the sample size is limited to 2 in 3D (M6 & CRM Case1 HLPW5)
    if (noffsets is None) or (noffsets > 2) or coarseXray: XRAYDIM1 = max(500, min(5000, XRAYDIM1))  #x1
    elif noffsets == 2: XRAYDIM1 = max(1500, min(15000, XRAYDIM1));                  #x3
    elif noffsets == 1: XRAYDIM1 = max(2500, min(25000, XRAYDIM1));                  #x5
    elif noffsets == 0: XRAYDIM1 = max(5000, min(50000, XRAYDIM1));                  #x10
    # XRAYDIM1 = max(500, min(5000, XRAYDIM1)); # XRAYDIM1 = max(5000, min(50000, XRAYDIM1)); is too expensive need to find another solution
    C._initVars(to, "cellNIn", 1.)

    to = X.blankCells(to, bodies1, BM, blankingType='node_in',
                      XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM1, dim=dim,
                      cellNName='cellN')
    if opt:
        ## needed - see comment in tagInsideOffset
        depthLocal = 0
        if noffsets == 0:   depthLocal = -2
        elif noffsets == 1: depthLocal = -1
        if abs(depthLocal) > 0: to = X.setHoleInterpolatedPoints(to, depth=depthLocal, cellNName='cellN', loc='nodes')
    to = C.node2Center(to,["cellN"])
    C._initVars(to, "{centers:indicatorTmp}=({centers:cellN}>0)")
    C._rmVars(to, ["cellN","centers:cellN"])
    o = Internal.getZones(to)[0]
    return o

def tagInsideOffset__(o, offset1=None, offset2=None, dim=3, h_target=-1., opt=False, noffsets=None, coarseXray=False):
    to = C.newPyTree(["OCTREE"]); to[2][1][2] = Internal.getZones(o)
    C._initVars(to, 'centers:indicatorTmp', 0.)

    isTbox = False
    if offset1 is None:
        offset1 = Internal.copyTree(offset2)
        isTbox = True

    if dim==2:
        offset1 = T.addkplane(offset1)
        offset2 = T.addkplane(offset2)
    bodies1 = [Internal.getZones(offset1)] # tag for refinement outside of offset1
    bodies2 = [Internal.getZones(offset2)] # tag for refinement inside of offset2

    BM = numpy.ones((1,1), dtype=Internal.E_NpyInt)

    bb1 = G.bbox(offset1)
    bb2 = G.bbox(offset2)
    L1 = max(bb1[3]-bb1[0], bb1[4]-bb1[1])
    L2 = max(bb2[3]-bb2[0], bb2[4]-bb2[1])
    if dim == 3:
        L1 = max(L1, bb1[5]-bb1[2])
        L2 = max(L1, bb2[5]-bb2[2])

    # ideally we should use blankCellsTri to avoid XRAYDIM but currently not safe
    XRAYDIM1 = int(L1/h_target)+10; XRAYDIM2 = int(L2/h_target)+10
    # [Temp. Patch] - This need to be generalized. Works for all test cases currently considered but the sample size is limited to 2 in 3D (M6 & CRM Case1 HLPW5)
    # XRAYDIM1 = min(5000, XRAYDIM1); XRAYDIM2 = min(5000, XRAYDIM2)
    # XRAYDIM1 = max(500, XRAYDIM1); XRAYDIM2 = max(500, XRAYDIM2)
    # Pull request note: causes regressions in the mesh generation
    if (noffsets is None) or (noffsets > 2) or coarseXray:
        XRAYDIM1 = max(500, min(5000, XRAYDIM1))  #x1
        XRAYDIM2 = max(500, min(5000, XRAYDIM2))  #x1
    elif noffsets == 2:
        XRAYDIM1 = max(1500, min(15000, XRAYDIM1)); #x3
        XRAYDIM2 = max(1500, min(15000, XRAYDIM2))  #x3
    elif noffsets == 1:
        XRAYDIM1 = max(2500, min(25000, XRAYDIM1)); #x5
        XRAYDIM2 = max(2500, min(25000, XRAYDIM2))  #x5
    elif noffsets == 0:
        XRAYDIM1 = max(5000, min(50000, XRAYDIM1)); #x10
        XRAYDIM2 = max(5000, min(50000, XRAYDIM2))  #x10

    C._initVars(to, "cellNOut", 1.)
    C._initVars(to, "cellNIn", 1.)

    to = X.blankCells(to, bodies2, BM, blankingType='node_in',
                      XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=dim,
                      cellNName='cellNIn')
    if isTbox: C._initVars(to, '{cellN}=({cellNIn}<1)')
    else:
        to = X.blankCells(to, bodies1, BM, blankingType='node_in',
                          XRaydim1=XRAYDIM1, XRaydim2=XRAYDIM2, dim=dim,
                          cellNName='cellNOut')
        # artificial shift the location of the boundary by depthLocal cells inwards (inside the body)
        # needed for the offset closest to the body as opt=True might coarsen certain critical regions (observed with CRM case 1 of the HLPW5)
        depthLocal = -1
        if opt:
            if noffsets == 0:   depthLocal = -3
            elif noffsets == 1: depthLocal = -2

        to = X.setHoleInterpolatedPoints(to, depth=depthLocal, cellNName='cellNOut', loc='nodes')
        C._initVars(to, '{cellN}=({cellNIn}<1)*({cellNOut}>0.)')

    to = C.node2Center(to, ["cellN"])
    C._initVars(to, "{centers:indicatorTmp}=({centers:cellN}>0.)")
    #
    G._getVolumeMap(to)
    vol_target = h_target**dim * 1.01 # add a tolerance
    C._initVars(to, "{centers:indicatorTmp}={centers:indicatorTmp}*({centers:vol}>%g)"%(vol_target))
    #
    C._rmVars(to, ["cellN", "cellNIn", "cellNOut", "centers:cellN", "centers:vol", "centers:h"])
    o = Internal.getZones(to)[0]
    return o

def createQuadSurfaceFromNgonPointListBigFace__(a, cranges, indices_owners=[], dimPb=3):
    a = Internal.getZones(a)[0] # get first zone
    # Recuperation des donnees de connectivite
    faces = Internal.getNGonNode(a)
    vol_cells = Internal.getNFaceNode(a)

    EC_faces = Internal.getNodeFromName(faces, "ElementConnectivity")[1]
    EC_volcells = Internal.getNodeFromName(vol_cells, "ElementConnectivity")[1]
    offset_faces = Internal.getNodeFromName(faces, "ElementStartOffset")[1]
    offset_volcells = Internal.getNodeFromName(vol_cells, "ElementStartOffset")[1]
    length_faces = offset_faces[1:] - offset_faces[:-1]

    coords_x = Internal.getNodeFromName(a, "CoordinateX")[1]
    coords_y = Internal.getNodeFromName(a, "CoordinateY")[1]
    coords_z = Internal.getNodeFromName(a, "CoordinateZ")[1]
    nb_vertices = len(coords_x)

    # Initialisation
    if indices_owners is None:
        indices_owners = range(len(offset_volcells) - 1)
    n_smallfaces = 4 if dimPb == 3 else 2

    big_faces = []
    # Traitement des volumes
    for idx_vol in indices_owners:
        stride = cranges[idx_vol]
        stride_offset = numpy.cumsum(stride)
        start = offset_volcells[idx_vol]

        # Faces candidates à fusionner
        for side, count in enumerate(stride):
            if count == n_smallfaces:
                start_idx = stride_offset[side - 1] if side > 0 else 0
                indices_faces = EC_volcells[start + start_idx : start + stride_offset[side]]

                conn_Nfaces = numpy.zeros((n_smallfaces, 4), dtype=Internal.E_NpyInt)

                indices_faces = numpy.array(indices_faces)
                face_lengths = length_faces[indices_faces - 1]
                quad_mask = face_lengths == 4
                non_quad_mask = ~quad_mask

                quad_indices = indices_faces[quad_mask]
                non_quad_indices = indices_faces[non_quad_mask]

                # On récupère directement les 4 nœuds pour chaque face
                n_quad_indices = len(quad_indices)
                n_non_quad_indices = len(non_quad_indices)

                # Création du tableau pour stocker les connectivités quad
                conn_Nfaces_quad = numpy.zeros((n_quad_indices, 4), dtype=Internal.E_NpyInt)
                conn_Nfaces_non_quad = numpy.zeros((n_non_quad_indices, 4), dtype=Internal.E_NpyInt)

                # quad faces
                f_start = offset_faces[quad_indices - 1]
                for i in range(4):
                    conn_Nfaces_quad[:, i] = EC_faces[f_start + i]

                # non quad faces
                i=0
                for face_idx in non_quad_indices:
                    f_start = offset_faces[face_idx - 1]
                    f_len = length_faces[face_idx - 1]
                    face_nodes = EC_faces[f_start : f_start + f_len]

                    xyz = numpy.vstack((coords_x[face_nodes - 1], coords_y[face_nodes - 1], coords_z[face_nodes - 1]))
                    center = (xyz.max(axis=1) + xyz.min(axis=1)) / 2
                    dist = numpy.sqrt(numpy.sum((xyz - center[:, numpy.newaxis])**2, axis=0))
                    mask = dist >= 0.8 * dist.max()
                    filtered_nodes = face_nodes[mask]
                    conn_Nfaces_non_quad[i, :len(filtered_nodes)] = filtered_nodes
                    i+=1
                conn_Nfaces = numpy.concatenate((conn_Nfaces_quad, conn_Nfaces_non_quad))

                # Reconstitution du gros quadrangle
                reorder = reorderNodesInCanonicalOrderForBigFace3D if dimPb == 3 else reorderNodesInCanonicalOrderForBigFace2D
                big_faces.append(reorder(conn_Nfaces))
    # Flatten et finalisation
    flattened_faces = numpy.array([pt for quad in big_faces for pt in quad], dtype=Internal.E_NpyInt)
    n_faces = len(big_faces)

    zone = Internal.newZone(name="Zone", zsize=[[nb_vertices, n_faces, 0]], ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=zone)
    Internal.newDataArray("CoordinateX", coords_x, parent=gc)
    Internal.newDataArray("CoordinateY", coords_y, parent=gc)
    Internal.newDataArray("CoordinateZ", coords_z, parent=gc)
    Internal.newElements("GridElements", etype=7, econnectivity=flattened_faces, erange=[1, n_faces], eboundary=0, parent=zone)
    return zone

def createQuadSurfaceFromNgonPointListSmallFace__(a, PL):
    a = Internal.getZones(a)[0] # get first zone
    faces = Internal.getNGonNode(a)
    #vol_cells = Internal.getNFaceNode(a)

    EC_faces = Internal.getNodeFromName(faces, "ElementConnectivity")[1]
    offset_faces = Internal.getNodeFromName(faces, "ElementStartOffset")[1]
    length_faces = offset_faces[1:] - offset_faces[:-1]

    coords_x = Internal.getNodeFromName(a, "CoordinateX")[1]
    coords_y = Internal.getNodeFromName(a, "CoordinateY")[1]
    coords_z = Internal.getNodeFromName(a, "CoordinateZ")[1]
    nb_vertices = len(coords_x)

    estimated_size = len(PL) * 4
    new_EC_faces = numpy.zeros(estimated_size, dtype=Internal.E_NpyInt)
    idx_new = 0

    for idx in PL:
        start = offset_faces[idx - 1]
        n_nodes = length_faces[idx - 1]
        face_nodes = EC_faces[start : start + n_nodes]

        if n_nodes == 4:
            new_EC_faces[idx_new : idx_new + 4] = face_nodes
            idx_new += 4
        else:
            nodes_idx = face_nodes - 1
            x, y, z = coords_x[nodes_idx], coords_y[nodes_idx], coords_z[nodes_idx]

            cx = (x.max() + x.min()) / 2
            cy = (y.max() + y.min()) / 2
            cz = (z.max() + z.min()) / 2
            dist = numpy.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
            mask = dist >= 0.8 * dist.max()

            kept = face_nodes[mask]
            count = len(kept)
            new_EC_faces[idx_new : idx_new + count] = kept
            idx_new += count

    new_EC_faces = new_EC_faces[:idx_new]

    # Création de la zone
    zone = Internal.newZone(name="Zone", zsize=[[nb_vertices, len(PL), 0]], ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=zone)
    Internal.newDataArray("CoordinateX", value=coords_x, parent=gc)
    Internal.newDataArray("CoordinateY", value=coords_y, parent=gc)
    Internal.newDataArray("CoordinateZ", value=coords_z, parent=gc)

    Internal.newElements(
        name="GridElements",
        etype=7,
        econnectivity=new_EC_faces,
        erange=[1, len(PL)],
        eboundary=0,
        parent=zone)

    return zone

def findFifthNode__(idx, offset_faces, length_faces, EC_faces, coordsx_a, coordsy_a, coordsz_a):
    idx_face_init = offset_faces[idx-1]
    len_face = length_faces[idx-1]
    EC_ngon_faces = EC_faces[idx_face_init:idx_face_init+len_face]
    x_5 = coordsx_a[EC_ngon_faces-1]
    y_5 = coordsy_a[EC_ngon_faces-1]
    z_5 = coordsz_a[EC_ngon_faces-1]

    x_center = (numpy.max(x_5)+numpy.min(x_5))/2
    y_center = (numpy.max(y_5)+numpy.min(y_5))/2
    z_center = (numpy.max(z_5)+numpy.min(z_5))/2

    dists_from_center = ((x_5-x_center)**2+(y_5-y_center)**2+(z_5-z_center)**2)**0.5
    max_dist = numpy.max(dists_from_center)
    point5_to_remove = numpy.where(numpy.abs(dists_from_center)<0.8*max_dist)[0]
    return point5_to_remove

def createPseudoBCQuadNQuadIntra__(a, owners, levels, halo_levels, neighbours, cranges, dimPb=3):
    point_list_small = []
    point_list_big = []
    indices_owners = []

    for PL in range(len(owners)):
        own = owners[PL]
        nei = neighbours[PL]
        own_lev = levels[own]
        nei_lev = levels[nei]
        if nei > -1:
            if own_lev != nei_lev:
                point_list_small.append(PL+1)
                point_list_big.append(PL+1)
                indices_owners.append(own)
                indices_owners.append(nei)
    if len(point_list_small) > 0:
        zsnc_small_internal = createQuadSurfaceFromNgonPointListSmallFace__(a, point_list_small)
    else:
        zsnc_small_internal = createEmptyQuadZone__()
    if len(point_list_big) > 0:
        zsnc_big_internal = createQuadSurfaceFromNgonPointListBigFace__(a, cranges, indices_owners, dimPb=dimPb)
    else:
        zsnc_big_internal = createEmptyQuadZone__()
    zone_nonconformal_in_proc = T.join(zsnc_big_internal, zsnc_small_internal)
    return zone_nonconformal_in_proc

def createPseudoBCQuadNQuadInter__(a, owners, levels, halo_levels, neighbours, cranges, dimPb=3):
    point_list_commpatch = []
    zgc_nodes = Internal.getNodesFromType(a, "GridConnectivity1to1_t")
    for zgc_node in zgc_nodes:
        PL = Internal.getNodeFromName(zgc_node, "PointList")[1]
        point_list_commpatch.append(PL)

    point_list_small = []
    point_list_big = []
    for nei_proc in range(len(zgc_nodes)):
        PL_nei_proc = point_list_commpatch[nei_proc]
        for i in range(len(PL_nei_proc)):
            PL = PL_nei_proc[i]
            own = owners[PL-1]
            own_lev = levels[own]
            nei_lev = halo_levels[nei_proc][i]
            if own_lev < nei_lev:
                point_list_big.append(PL)
            elif own_lev > nei_lev:
                point_list_small.append(PL)
    if len(point_list_small)>0:
        zsnc_small_MPI = createQuadSurfaceFromNgonPointListSmallFace__(a, point_list_small)
    else:
        zsnc_small_MPI = createEmptyQuadZone__()

    if len(point_list_big)>0:
        point_list_big = numpy.array(point_list_big)
        indices_owners = owners[point_list_big-1]
        zsnc_big_MPI = createQuadSurfaceFromNgonPointListBigFace__(a, cranges, indices_owners, dimPb)
    else:
        zsnc_big_MPI = createEmptyQuadZone__()

    zone_nonconformal_between_procs = T.join(zsnc_big_MPI, zsnc_small_MPI)
    return zone_nonconformal_between_procs

def reorderNodesInCanonicalOrderForBigFace2D(conn_2faces):
    point0=None; point1=None; point2=None; point3=None
    unique, counts = numpy.unique(conn_2faces, return_counts=True)
    indices_count_2 = numpy.argwhere(counts==2)
    indices_count_2 = indices_count_2.reshape(2)
    indices_count_1 = numpy.argwhere(counts==1)
    indices_count_1 = indices_count_1.reshape(4)
    points45 = unique[indices_count_2]
    corners_sort = unique[indices_count_1]
    point0 = corners_sort[0] # small index
    if point0 in conn_2faces[0]:
        face0 = 0
        face1 = 1
    else:
        face0 = 1
        face1 = 0
    conn_face0 = conn_2faces[face0]
    pos_point0 = numpy.where(conn_face0==point0)[0][0]

    if conn_face0[(pos_point0+1)%4] in points45:
        point4 = conn_face0[(pos_point0+1)%4]
        point5 = conn_face0[(pos_point0+2)%4]
        point3 = conn_face0[(pos_point0+3)%4]
    else:
        point4 = conn_face0[(pos_point0-1)%4]
        point5 = conn_face0[(pos_point0-2)%4]
        point3 = conn_face0[(pos_point0-3)%4]

    conn_face1 = conn_2faces[face1]
    pos_point4 = numpy.where(conn_face1==point4)[0][0]

    if conn_face1[(pos_point4+1)%4] in points45:
        point2 = conn_face1[(pos_point4+2)%4]
        point1 = conn_face1[(pos_point4+3)%4]
    else:
        point2 = conn_face1[(pos_point4-2)%4]
        point1 = conn_face1[(pos_point4+1)%4]

    if point1>point3:
        tmp = point1
        point1 = point3
        point3 = tmp
    return [point0, point1, point2, point3]

def reorderNodesInCanonicalOrderForBigFace3D(conn_4faces):
    point0 = None; point1 = None; point2 = None; point3 = None
    unique, counts = numpy.unique(conn_4faces, return_counts=True)
    index_count_4 = numpy.argwhere(counts==4)[0][0]
    indices_count_1 = numpy.argwhere(counts==1)
    indices_count_1 = indices_count_1.reshape(4)
    point8 = unique[index_count_4]
    corners = unique[indices_count_1]

    corner_points = []
    pos_corner_points = []
    for EC in conn_4faces:
        pos8 = numpy.where(EC==point8)[0][0]
        pos_corner_point = (pos8 + 2) % 4
        pos_corner_points.append(pos_corner_point)
        corner_points.append(EC[pos_corner_point])

    elt0 = numpy.argmin(corner_points)
    point0 = corner_points[elt0]
    pos_corner_in_elt0 = pos_corner_points[elt0]

    point4 = conn_4faces[elt0][(pos_corner_in_elt0+1)%4]
    point7 = conn_4faces[elt0][(pos_corner_in_elt0-1)]

    elt1 = numpy.where(conn_4faces==point4)[0]
    elt1 = elt1[elt1!=elt0][0]

    point1 = corner_points[elt1]
    #pos_corner_in_elt1 = pos_corner_points[elt1]

    point5 = list(set(conn_4faces[elt1])-set([point4, point8, point1]))[0]

    elt2 = numpy.where(conn_4faces==point5)[0]
    elt2 = elt2[elt2!=elt1][0]

    point2 = corner_points[elt2]
    #pos_corner_in_elt2 = pos_corner_points[elt2]

    point6 = list(set(conn_4faces[elt2])-set([point5,point8,point2]))[0]
    elt3 = list(set([0, 1, 2, 3])-set([elt0, elt1, elt2]))[0]
    point3 = list(set(conn_4faces[elt3])-set([point6, point7, point8]))[0]
    return [point0, point1, point2, point3]

def createEmptyQuadZone__():
    zone = Internal.newZone(name="empty", zsize=[[0,0]], ztype="Unstructured")
    gc = Internal.newGridCoordinates(parent=zone)
    Internal.newDataArray('CoordinateX', value=numpy.empty(0), parent=gc)
    Internal.newDataArray('CoordinateY', value=numpy.empty(0), parent=gc)
    Internal.newDataArray('CoordinateZ', value=numpy.empty(0), parent=gc)
    elt = Internal.newElements(name="GridElements", etype=7, econnectivity=numpy.empty(0), erange=[1, 0], eboundary=0, parent=zone)
    return zone

def _createQuadConnectivityFromNgonPointList__(a_hexa, a, PL, bcname, bctype):
    a = Internal.getZones(a)[0] # get first zone
    faces = Internal.getNGonNode(a)
    #vol_cells = Internal.getNFaceNode(a)

    EC_faces = Internal.getNodeFromName(faces, "ElementConnectivity")[1]
    offset_faces = Internal.getNodeFromName(faces, "ElementStartOffset")[1]
    length_faces = offset_faces[1:] - offset_faces[:-1]

    coords_x = Internal.getNodeFromName(a, "CoordinateX")[1]
    coords_y = Internal.getNodeFromName(a, "CoordinateY")[1]
    coords_z = Internal.getNodeFromName(a, "CoordinateZ")[1]

    # --- Initialisation ---
    estimated_len = len(PL) * 4
    new_connectivity = numpy.zeros(estimated_len, dtype=Internal.E_NpyInt)
    idx_new = 0

    # --- Boucle optimisée ---
    # Face indices à traiter
    PL2 = numpy.array(PL)

    # Longueurs des faces sélectionnées
    face_lengths = length_faces[PL2 - 1]

    # Masque des faces à 4 noeuds
    quad_mask = face_lengths == 4
    non_quad_mask = ~quad_mask

    # Indices quads et non quads
    quad_indices = PL2[quad_mask]
    non_quad_indices = PL2[non_quad_mask]

    # Préallocation (optimiste)
    new_connectivity = numpy.zeros(len(PL2) * 4, dtype=Internal.E_NpyInt)
    idx_new = 0

    # --- 1. Traitement des faces quad
    for idx in quad_indices:
        start = offset_faces[idx - 1]
        new_connectivity[idx_new : idx_new + 4] = EC_faces[start : start + 4]
        idx_new += 4

    # --- 2. Traitement des faces non quad ---
    for idx in non_quad_indices:
        start = offset_faces[idx - 1]
        n_nodes = length_faces[idx - 1]
        nodes = EC_faces[start : start + n_nodes]

        xyz = numpy.vstack((coords_x[nodes - 1], coords_y[nodes - 1], coords_z[nodes - 1]))
        center = (xyz.max(axis=1) + xyz.min(axis=1)) / 2
        dist = numpy.sqrt(numpy.sum((xyz - center[:, numpy.newaxis])**2, axis=0))
        keep_mask = dist >= 0.8 * dist.max()

        kept = nodes[keep_mask]
        count = len(kept)
        new_connectivity[idx_new : idx_new + count] = kept
        idx_new += count

    # --- Ajout des nouveaux éléments à la zone ---
    zone = Internal.getZones(a_hexa)[0]
    elt_nodes = Internal.getNodesFromType(zone, "Elements_t")
    last_id = Internal.getNodeFromName(elt_nodes[-1], "ElementRange")[1][1]

    Internal.newElements(
        name=bcname,
        etype=7,
        econnectivity=new_connectivity[:idx_new],
        erange=[last_id + 1, last_id + len(PL)],
        eboundary=0,
        parent=zone)
    C._addBC2Zone(zone, bcname, bctype, elementRange=[last_id + 1, last_id + len(PL)])
    return None

def _createBCNearMatch__(a_hexa, bcnearmatch):
    z = Internal.getZones(a_hexa)
    if Internal.getValue(bcnearmatch)[0][1] > 0:
        f = P.exteriorFaces(z[0])
        hook = C.createHook(f, "elementCenters")
        ids = C.identifyElements(hook, bcnearmatch, tol=__TOL__)

        ids = ids[ids[:] > -1]
        ids = ids.tolist()
        ids = [ids[i]-1 for i in range(len(ids))]
        if len(ids)>0:
            zf = T.subzone(f,ids, type='elements')
            _addBC2Zone__(z[0], "QuadNQuad", "FamilySpecified:QuadNQuad", zf)
        C.freeHook(hook)
    return None

def _createBCStandard__(a_hexa, a):
    nodes_bcs = Internal.getNodesFromType(a, "BC_t")
    for node_bc in nodes_bcs:
        bctype = Internal.getValue(node_bc)
        bcname = Internal.getName(node_bc)
        PL = Internal.getNodeFromName(node_bc, "PointList")[1]
        _createQuadConnectivityFromNgonPointList__(a_hexa, a, PL, bcname, bctype)
    return None

def adaptMesh__(fileSkeleton, hmin, tb, bbo, toffset=None, dim=3, loadBalancing=False, opt=False, numTbox=0):
    from mpi4py import MPI # for MPI_Init
    coarseXray = False
    bbtb = G.bbox(tb)
    lenMax = 0.0
    for i in range(dim): lenMax = max(bbo[i+3]-bbo[i], lenMax)
    if lenMax/hmin < 100: coarseXray = True

    numBase = len(Internal.getZones(tb))
    o, res = XC.loadAndSplitNGon(fileSkeleton)
    Cmpi.barrier()
    gcells = res[5]
    gfaces = res[6]
    comm = res[1]
    if dim == 3: normal2D=None
    else: normal2D = numpy.array([0.0, 0.0, 1.0])
    hookAM = XC.AdaptMesh_Init(o, normal2D, comm=comm, gcells=gcells, gfaces=gfaces)
    offset_zones = Internal.getZones(toffset)
    offset_inside = [Internal.getZones(tb)]
    noffsetBase = []

    sortDicOffsetIBM  = {}
    sortDicOffsetTbox = {}
    newOffsetsTbox    = []
    # combine zones with the same snear
    for i in offset_zones:
        # dict={snear1:[OffsetzoneName1,OffsetzoneName2], snear2:[OffsetzoneName3,OffsetzoneName4]}
        if 'Tbox' in i[0]:
            hminLocal = Internal.getValue(Internal.getNodeFromName2(i, 'snear'))
            _addItemDict__(sortDicOffsetTbox, hminLocal, i[0])
        else:
            hminLocal = Internal.getValue(Internal.getNodeFromName2(i, 'snear'))
            _addItemDict__(sortDicOffsetIBM, hminLocal, i[0])

    # newOffsetsIBM=[[Offsetzone1,Offsetzone2], [Offsetzone3,Offsetzone4]]
    # corresponding for tbox also
    newOffsetsIBM = []
    for snearLocal in sortDicOffsetIBM:
        tmpOffset = []
        for tt in sortDicOffsetIBM[snearLocal]:tmpOffset.append(Internal.getNodeFromName(offset_zones, tt))
        newOffsetsIBM.append(tmpOffset)

    if numTbox > 0:
        for snearLocal in sortDicOffsetTbox:
            tmpOffset = []
            for tt in sortDicOffsetTbox[snearLocal]:tmpOffset.append(Internal.getNodeFromName(offset_zones, tt))
            newOffsetsTbox.append(tmpOffset)
        offset_inside.append(None)

    noffsetBase = [len(newOffsetsIBM), len(newOffsetsTbox)]
    offset_zonesNew = [newOffsetsIBM, newOffsetsTbox]
    offset_name = ['IBM body', 'tbox'] # for I/O purposes only
    noffsets = max(noffsetBase) # max offset for all the bases
    for i in range(noffsets-1, -1,-1):
        if Cmpi.master: print('\n------------------------> Adapt Offset level %d ... start'%i, flush=True)
        for nBase in range(1+min(1,numTbox)): #nBase = 0(IBM), 1(tbox - if it exists))
            # if i (offset) is greater than the num. offset for nBase - we continue to the next offset number
            if i > noffsetBase[nBase]-1: continue
            if Cmpi.master: print("~~~~~~~~~~Base %s AdaptMesh...start"%offset_name[nBase], flush=True)
            offsetloc = offset_zonesNew[nBase][i][0] # as all offset at this level have the same snear we just take the first one in the list
            hminLocal = Internal.getValue(Internal.getNodeFromName2(offsetloc, 'snear'))
            hx        = hminLocal# * 2**i
            adaptPass = 0
            adapting  = True
            while adapting:
                Ncells = Cmpi.getNCells(o)
                C._initVars(o,'centers:indicator',0.)
                # loop through the offsets in the list for base=nBase and offset level=i
                for numOffTmp, offsetlocTmp in enumerate(offset_zonesNew[nBase][i]):
                    # body offset: tag the region between the body & the offset
                    # tbox offset: tage the region enclosed by the offset
                    o = tagInsideOffset__(o,  offset1=offset_inside[nBase], offset2=offsetlocTmp, dim=dim, h_target=hx, opt=opt, noffsets=i, coarseXray=coarseXray)
                    C._initVars(o,"{centers:indicator}={centers:indicator}+{centers:indicatorTmp}")
                    C._rmVars(o, ["centers:indicatorTmp"])
                # Pull request note: tagOutsideBody causes regressions in the mesh generation for test cases: Connector/prepAMRFull_*.py
                for offsetlocTmp in offset_inside[0]:
                    # tag cellN=0 the region enclosed inside the body - need to avoid adapting inside the body when the tbox cuts the body
                    o = tagOutsideBody__(o, tbTMP=offsetlocTmp, dim=dim, h_target=hx, opt=opt, noffsets=i, coarseXray=coarseXray)
                    C._initVars(o,"{centers:indicator}={centers:indicator}*{centers:indicatorTmp}")
                    C._rmVars(o, ["centers:indicatorTmp"])

                indicMax = C.getMaxValue(o,"centers:indicator")
                indicMax = Cmpi.allgather(indicMax)
                indicMax = max(indicMax)
                if indicMax<1. or (i == 0 and adaptPass > 0):
                    adapting=False
                    C._rmVars(o,["centers:indicator"])
                    break
                else:
                    if Cmpi.master: print("......Recursive AdaptMesh:: level %d...start"%adaptPass, flush=True)
                    f = Internal.getNodeFromName(o, 'indicator')[1]
                    REF = f.astype(dtype=Internal.E_NpyInt)
                    XC.AdaptMesh_AssignRefData(hookAM, REF)
                    # [TODO] check if 4 cells is a global limit
                    # it was found that at least 4 cells are required if not the PT-Scotch partitioning can have partitions with 0 cells. This is empirical and subject to change as the sample
                    # size of the test cases increases
                    if loadBalancing and Ncells//Cmpi.size>3: XC.AdaptMesh_LoadBalance(hookAM)
                    XC.AdaptMesh_Adapt(hookAM)
                    o = XC.AdaptMesh_ExtractMesh(hookAM, conformize=1)
                    o = Internal.getZones(o)[0]
                    # Cmpi.convertPyTree2File(o,'check_loadBalancingFix%d_pass%d.cgns'%(i,adaptPass)) # Leave here for now. very useful for debugging
                    if Cmpi.master: print("......Recursive AdaptMesh:: level %d...end"%adaptPass, flush=True)
                    adaptPass+=1
            if Cmpi.master: print("~~~~~~~~~~Base %s AdaptMesh...end"%offset_name[nBase], flush=True)
    if Cmpi.master: print('------------------------> Adapt Offset level %d ... end'%i, flush=True)

    o = XC.AdaptMesh_ExtractMesh(hookAM, conformize=1) #ok - base per proc
    o = Internal.getZones(o)[0]
    owners = XC.AdaptMesh_ExtractOwners(hookAM)
    levels = XC.AdaptMesh_ExtractCellLevels(hookAM)
    halo_levels = XC.AdaptMesh_ExtractHaloCellLevels(hookAM)
    neighbours = XC.AdaptMesh_ExtractNeighbours(hookAM)
    cranges = XC.AdaptMesh_ExtractCellRanges(hookAM)
    cart_hexa = XC.AdaptMesh_ExtractMesh(hookAM, conformize=0)
    XC.AdaptMesh_Exit(hookAM)

    # BCs
    zone_nonconformal_inter = createPseudoBCQuadNQuadInter__(o, owners, levels, halo_levels, neighbours, cranges, dimPb=dim)
    zone_nonconformal_intra = createPseudoBCQuadNQuadIntra__(o, owners, levels, halo_levels, neighbours, cranges, dimPb=dim)
    zone_nonconformal = T.join(zone_nonconformal_inter, zone_nonconformal_intra)
    _createBCNearMatch__(cart_hexa, zone_nonconformal)
    _createBCStandard__(cart_hexa, o)
    del o

    Cmpi._setProc(cart_hexa, Cmpi.rank)
    return cart_hexa

def _addBC2Zone__(z, bndName, bndType, zbc, loc='FaceCenter', zdnrName=None):
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
    node = Internal.createUniqueChild(z, bndName, 'Elements_t', value=[eltType, nebb])
    Internal.createUniqueChild(node, 'ElementRange', 'IndexRange_t', value=[maxElt+1, maxElt+neb])
    oldc = Internal.getNodeFromName2(zbc, 'ElementConnectivity')[1]
    newc = numpy.copy(oldc)
    hook = C.createHook(z, 'nodes')
    ids = C.identifyNodes(hook, zbc)
    newc[:] = ids[oldc[:]-1]
    Internal.createUniqueChild(node, 'ElementConnectivity', 'DataArray_t', value=newc)

    zoneBC = Internal.createUniqueChild(z, 'ZoneBC', 'ZoneBC_t')
    if len(s)==1:
        info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType)
    else: # familyspecified
        info = Internal.createChild(zoneBC, bndName, 'BC_t', value=bndType1)
        Internal.createUniqueChild(info, 'FamilyName', 'FamilyName_t',value=bndType2)

    Internal.createUniqueChild(info, 'GridLocation', 'GridLocation_t', value='FaceCenter')
    Internal.createUniqueChild(info, 'ElementRange', 'IndexRange_t', value=numpy.array([[maxElt+1, maxElt+neb]]))
    return None

# determine where the symmetry plane is by dir
# 1 : xmin, 2 : ymin, 3 = zmin
def getSymmetryPlaneInfo__(tb, dim=3):
    baseSYM = Internal.getNodesFromName1(tb, "SYM")
    dir_sym = 0
    if baseSYM:
        symplane = []
        for zsym in Internal.getZones(baseSYM):
            if C.getMaxValue(zsym,'centers:cellN')>0.: symplane.append(zsym)
        [xmin, ymin, zmin, xmax, ymax, zmax] = G.bbox(symplane)
        if abs(xmax-xmin) < __TOL__:
            dir_sym=1;
        elif abs(ymax-ymin) < __TOL__:
            dir_sym=2;
        elif abs(zmax-zmin) < __TOL__ and dim == 3:
            dir_sym=3;
    return dir_sym

#==================================================================
# Addition of Physical BCs on the skeleton NGON Mesh
#==================================================================
def _addPhysicalBCs__(z_ngon, tb, dim=3):
    bbo = G.bbox(z_ngon)

    dir_sym = getSymmetryPlaneInfo__(tb,dim=dim)

    xmin = bbo[0]; ymin = bbo[1]; zmin = bbo[2]
    xmax = bbo[3]; ymax = bbo[4]; zmax = bbo[5]

    extFaces = P.exteriorFaces(z_ngon)
    extFaces = T.splitSharpEdges(extFaces, 80.)

    zbc_xmin = []; zbc_xmax = [];
    zbc_ymin = []; zbc_ymax = [];
    zbc_zmin = []; zbc_zmax = []
    for e in extFaces:
        bbe = G.bbox(e)
        xmine = bbe[0]; ymine = bbe[1]; zmine = bbe[2]
        xmaxe = bbe[3]; ymaxe = bbe[4]; zmaxe = bbe[5]
        dx = xmaxe-xmine; dy = ymaxe-ymine; dz = zmaxe-zmine
        if dx < __TOL__ and abs(xmine-xmin) < __TOL__: zbc_xmin.append(e)
        elif dx < __TOL__ and abs(xmaxe-xmax) < __TOL__: zbc_xmax.append(e)
        elif dy < __TOL__ and abs(ymine-ymin) < __TOL__: zbc_ymin.append(e)
        elif dy < __TOL__ and abs(ymaxe-ymax) < __TOL__: zbc_ymax.append(e)
        elif dz < __TOL__ and abs(zmine-zmin) < __TOL__: zbc_zmin.append(e)
        elif dz < __TOL__ and abs(zmaxe-zmax) < __TOL__: zbc_zmax.append(e)

    zbc_xmin = T.join(zbc_xmin); zbc_xmin[0] = 'xmin'
    zbc_xmax = T.join(zbc_xmax); zbc_xmax[0] = 'xmax'
    zbc_ymin = T.join(zbc_ymin); zbc_ymin[0] = 'ymin'
    zbc_ymax = T.join(zbc_ymax); zbc_ymax[0] = 'ymax'
    zbc_zmin = T.join(zbc_zmin); zbc_zmin[0] = 'zmin'
    zbc_zmax = T.join(zbc_zmax); zbc_zmax[0] = 'zmax'
    if dir_sym == 1:
        zbc_xmin[0] = C.getZoneName("BCSymmetryPlane")
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane", subzone=zbc_xmin)
        C._addBC2Zone(z_ngon, 'BCFarfield'     , "BCFarfield"    , subzone=zbc_xmax)
    else:
        zbc_xmin[0] = C.getZoneName("BCFarfield")
        C._addBC2Zone(z_ngon, 'BCFarfield', "BCFarfield", subzone=zbc_xmin)
        C._addBC2Zone(z_ngon, 'BCFarfield', "BCFarfield", subzone=zbc_xmax)

    if dir_sym == 2:
        zbc_ymin[0] = C.getZoneName("BCSymmetryPlane")
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane',"BCSymmetryPlane", subzone=zbc_ymin)
        C._addBC2Zone(z_ngon, 'BCFarfield'     , "BCFarfield"    , subzone=zbc_ymax)

    else:
        zbc_ymin[0] = C.getZoneName("BCFarfield")
        C._addBC2Zone(z_ngon, 'BCFarfield', "BCFarfield", subzone=zbc_ymin)
        C._addBC2Zone(z_ngon, 'BCFarfield', "BCFarfield", subzone=zbc_ymax)

    if dim == 2:
        zbc_zmin[0] = C.getZoneName("BCSymmetryPlane")
        zbc_zmax[0] = C.getZoneName("BCSymmetryPlane")
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane', "BCSymmetryPlane", subzone=zbc_zmin)
        C._addBC2Zone(z_ngon, 'BCSymmetryPlane', "BCSymmetryPlane", subzone=zbc_zmax)
    else:
        if dir_sym == 3:
            zbc_zmin[0] = C.getZoneName("BCSymmetryPlane")
            C._addBC2Zone(z_ngon, 'BCSymmetryPlane', "BCSymmetryPlane", subzone=zbc_zmin)
            C._addBC2Zone(z_ngon, 'BCFarfield'     , "BCFarfield"     , subzone=zbc_zmax)

        else:
            zbc_zmin[0] = C.getZoneName("BCFarfield")
            C._addBC2Zone(z_ngon, 'BCFarfield', "BCFarfield", subzone=zbc_zmin)
            C._addBC2Zone(z_ngon, 'BCFarfield', "BCFarfield", subzone=zbc_zmax)

    return None

#==================================================================
# Generation of the AMR mesh for IBMs
# MAIN FUNCTION
# opt = True: for offset surface generation if it takes too long (depending on the resolution of tb)
#==================================================================
def generateAMRMesh(tb, toffset=None, levelMax=7, vmins=11, snears=0.01, dfars=10, dim=3, check=False,
                    opt=False, loadBalancing=False, octreeMode=0, localDir='./', tbox=None, vminsTbox=None, tbv2=None):
    NumMinDxLarge=1
    Cmpi.trace('AMR Mesh Generation...start', master=True)
    fileSkeleton = 'skeleton.cgns'
    pathSkeleton = os.path.join(localDir, fileSkeleton)

    # NumMinDxLarge : min. number of cells at the largest Deltax between the first refinement & the domain boundary conditions
    # e.g. NumMinDxLarge = 1
    #      |   |   |
    #      |   |---|
    #      |   |   |
    # e.g. NumMinDxLarge = 2
    #      |   |   |   |
    #      |   |   |---|
    #      |   |   |   |
    if NumMinDxLarge<1:
        if Cmpi.master: print("NumMinDxLarge - the min. num. of cells in any direction for the max. Dx - cannot be less than 1. Setting NumMinDxLarge to at least 1.", flush=True)
        NumMinDxLarge=1
    NumMinDxLarge+=1 # we add one as the check later on is on nodes & not cells.

    # check to see if the sym bases have snear?
    baseSYM    = Internal.getNodesFromName1(tb,"SYM")
    isSymLocal = False
    if baseSYM:
        zoneTmp = Internal.getZones(baseSYM)
        if Internal.getNodeFromName(baseSYM, 'snear'):
            isSymLocal = True

    checkBaseNames__(tb,tbox)

    # snears & numBase for the input tb only.
    # If isSymLocal snear includes the 6 for the symb base & a copy of all the oringal zones.
    # NumBase is only the num. of 'real' bodies
    snears, numBase = getListSnear__(tb, snears)
    snearsFlat = [item for sub in snears for item in sub] # needed for G.octree
    if tbv2 is not None: tbv2 = C.convertFile2PyTree(tbv2)

    # list of vmins
    # This section::
    # ================== SECTION START ==================
    # Checks that the vmin input is correct & if needed corrects it (if info. is missing it will apply default values & default copies)
    vmins=vminsInputCheck__(vmins, numBase, levelMax)

    snearsTbox = []
    numTbox    = 0
    tb_tbox    = Internal.copyTree(tb)
    if tbox:
        numTbox    = len(Internal.getBases(tbox))
        tb_tbox[2]+=Internal.getBases(tbox)

        if vminsTbox is None:
            vminsTbox = numpy.ones((numTbox,levelMax))*5
            vminsTbox = vminsTbox.tolist()
        else:
            vminsTbox=vminsInputCheck__(vminsTbox, numTbox, levelMax)

        Cmpi.barrier()
        snearsTbox, tmpRMV = getListSnear__(tbox, 1)
        vmins.extend(vminsTbox)

    # add the tbox snear parameters
    snearEnd = len(snears)
    if isSymLocal: snearEnd = -1
    snearsTbTbox = snears[0:snearEnd]
    if snearsTbox: snearsTbTbox.extend(snearsTbox)

    numBase = numBase + numTbox
    vminsLocal = numpy.ones((numBase, levelMax))
    for nBase in range(numBase):
        if not isinstance(vmins[nBase],list):
            vminsLocal[nBase][:] = vmins[nBase][:]
        elif len(vmins[nBase]) < levelMax:
            vminsLocal[nBase][:len(vmins[nBase])] = vmins[nBase][:]
            vminsLocal[nBase][len(vmins[nBase]):] = vmins[nBase][-1]
        elif len(vmins[nBase]) > levelMax:
            vminsLocal[nBase][:] = vmins[nBase][:levelMax]
        else:
            vminsLocal[nBase][:] = vmins[nBase][:]
    vmins=[]
    for nBase in range(numBase):
        vmins.append(list(vminsLocal[nBase]))
        vmins[nBase] = [max(5,v) for v in vmins[nBase]] # vmin values should not be inferior to a given threshold
    # ================== SECTION END  ==================

    # levelSkel: initial refinement level of the skeleton octree
    # might be tuned
    # Input: tb, no tbox; snears of tb only in flat (list of zone snears), dfars - integer value is ok
    o, newLevelMax = generateSkeletonMesh__(tb, snears=snearsFlat, dfars=dfars, dim=dim, levelSkel=levelMax, octreeMode=octreeMode)
    if newLevelMax != levelMax:
        if Cmpi.master: print('Warning: modified number of AMR Levels. Old levelMax = %d || New levelMax = %d'%(levelMax, newLevelMax), flush=True)
        for nBase in range(numBase):
            while len(vmins[nBase]) < newLevelMax: vmins[nBase].append(vmins[nBase][-1]) # if newLevelMax > levelMax
            vmins[nBase]    = vmins[nBase][:newLevelMax] # if newLevelMax < levelMax
            levelMax = newLevelMax
    G._getVolumeMap(o)
    hmin_skel = (C.getMinValue(o,"centers:vol"))**(1/dim)
    hmin = hmin_skel * 2 ** (-levelMax)
    if Cmpi.master: print(" Minimum spacing = ", hmin, hmin_skel, flush=True)

    # Leaving here for now just in case we decided to reactive this option in the future. The chosen option is to reduce the cell size (increase the # of cells) of the octree skeleton grid.
    #Ncells = C.getNCells(o)
    #if Ncells < Cmpi.size:
    #    raise ValueError('There are more MPI processes (Nmpi) [%d] than number of cells in the background skeleton/octree mesh (Ncells) [%d]. Note: Nmpi ≤ Ncells. Exiting...'%(Cmpi.size, Ncells))
    #    Cmpi.abort(errorcode=1)

    # Modifies the snears such that they are always a multiple of the smallest snear
    # Needed to quarantee a smooth transition of the flow field
    # If octree=0, this guarantees that the min snear previous obtained is propagated to the other snears
    minSnearsOrig = min(snearsFlat)
    if abs(hmin-minSnearsOrig) > __TOL__:
        for nBase in range(numBase):
            # we only check the first zone of each base as CODA (currently)
            # cant handle snear changes in the same base - no change in res at the surface of the IB (immersed boundary)
            snearMult = snearsTbTbox[nBase][0]//minSnearsOrig
            snearsTbTbox[nBase] = [snearMult*hmin]

    # mandatory save file for loadAndSplit for adaptation
    if Cmpi.master: C.convertPyTree2File(o, pathSkeleton)
    Cmpi.barrier()

    bbo = G.bbox(o)

    dir_sym      = getSymmetryPlaneInfo__(tb, dim=dim)
    # adaptation of the mesh wrt to the bodies (finest level) and offsets
    # only a part is returned per processor
    tb_tboxLocal = Internal.copyTree(tb_tbox)
    if baseSYM:
        # Remove SYM Base & Zones - keep real closed tb & tbox
        tb_tboxLocal = Internal.rmNodesByNameAndType(tb_tboxLocal, 'SYM', 'CGNSBase_t')
        tb_tboxLocal = Internal.rmNodesByNameAndType(tb_tboxLocal, '*_sym*', 'Zone_t')

    # New calc. of dfar max. This is done based on the distance between the IBM & the domain edges.
    # The min of this becomes the dfar max. Need to consider the tbox also as cannot have hanging nodes @ BC
    # is per base so there are no cross-contamination constraints (e.g. tbox constaining tb)
    dfarmaxLocal = []
    for nBase, tbLocal in enumerate(Internal.getBases(tb_tboxLocal)):
        bbTbLocal = G.bbox(tbLocal)
        dfarmax   = 1e10
        for i in range(dim):
            if i+1 == dir_sym: dfarmaxTmp = abs(bbTbLocal[i+3]-bbo[i+3])
            else:              dfarmaxTmp = min(abs(bbTbLocal[i]-bbo[i]), abs(bbTbLocal[i+3]-bbo[i+3]))
            dfarmax    = min(dfarmaxTmp, dfarmax)
        dfarmaxLocal.append(dfarmax-NumMinDxLarge*hmin_skel)

    # Old dfar max calc. Stays here in case it is needed in the future
    dfarmax = min(bbo[3]-bbo[0], bbo[4]-bbo[1])
    if dim==3: dfarmax = min(dfarmax, bbo[5]-bbo[2])

    if toffset==None:
        offsetValues = []
        for nBase in range(numBase):
            offsetprev       = 0.
            offsetValuesBase = []
            for no_adapt in range(len(vmins[nBase])):
                hminLocal = snearsTbTbox[nBase][0]
                offsetloc = offsetprev + hminLocal*(2**no_adapt)*vmins[nBase][no_adapt]
                # Pull request note: the line below causes regressions in the mesh generation
                if offsetloc < 0.99*dfarmaxLocal[nBase]:
                    #if offsetloc < 0.99*dfarmax: # for old dfar max cal. Stays here in case it is needed in the future
                    offsetValuesBase.append(offsetloc)
                    offsetprev=offsetloc
            if not offsetValuesBase:
                no_adapt = 0
                offsetloc = offsetprev + hminLocal*(2**no_adapt)*vmins[nBase][no_adapt]
                raise ValueError('Base #%d has no offset values. The first offset (closest to the body) is at a distance of %g which is larger than the max allowable distance of %g. Exiting...'%(nBase, offsetloc, 0.99*dfarmaxLocal[nBase]))
                Cmpi.abort(errorcode=1)
            offsetValues.append(offsetValuesBase)

        # generate list of offsets
        # tb & tbox
        if Cmpi.master: print("Generate list of offsets for rank ", Cmpi.rank, flush=True)
        toffset = generateListOfOffsets__(tb_tbox, snearsTbTbox, offsetValues=offsetValues, dim=dim, opt=opt, numTbox=numTbox, tbv2=tbv2)
        if check and Cmpi.master: C.convertPyTree2File(toffset, os.path.join(localDir, "offset.cgns"))
    Cmpi.barrier()

    # adaptation of the mesh wrt to the bodies (finest level) and offsets
    # only a part is returned per processor
    if baseSYM:
        # Remove SYM Base & Zones - keep real closed tb
        tb = Internal.rmNodesByNameAndType(tb, 'SYM', 'CGNSBase_t')
        tb = Internal.rmNodesByNameAndType(tb, '*_sym*', 'Zone_t')

    Cmpi.barrier()
    # only tb --> for blanking & tagging inside the geometry
    o = adaptMesh__(pathSkeleton, hmin, tb, bbo, toffset=toffset, dim=dim, loadBalancing=loadBalancing, opt=opt, numTbox=numTbox)
    Cmpi.trace('AMR Mesh Generation...end', master=True)
    return o # requirement for X_AMR (one zone per base, one base per proc)



# ORIG F.Basile files - kept 'just in case'
#def createQuadSurfaceFromNgonPointListBigFaceOrig__(a, cranges, indices_owners=[], dimPb=3):
#    a = Internal.getZones(a)[0] # get first zone
#    faces = Internal.getNGonNode(a)
#    vol_cells = Internal.getNFaceNode(a)
#
#    EC_volcells = Internal.getNodeFromName(vol_cells, "ElementConnectivity")[1]
#    EC_faces = Internal.getNodeFromName(faces, "ElementConnectivity")[1]
#
#    ER_volcells = Internal.getNodeFromName(vol_cells, "ElementRange")[1]
#    ER_faces = Internal.getNodeFromName(faces, "ElementRange")[1]
#
#    offset_volcells = Internal.getNodeFromName(vol_cells, "ElementStartOffset")[1]
#    offset_faces = Internal.getNodeFromName(faces, "ElementStartOffset")[1]
#
#    length_volcells = offset_volcells[1:] - offset_volcells[:-1]
#    length_faces = offset_faces[1:] - offset_faces[:-1]
#    lens_vol = numpy.unique(length_volcells)
#    lens_faces = numpy.unique(length_faces)
#    coordsx_a = Internal.getNodeFromName(a, "CoordinateX")[1]
#    coordsy_a = Internal.getNodeFromName(a, "CoordinateY")[1]
#    coordsz_a = Internal.getNodeFromName(a, "CoordinateZ")[1]
#    nb_vertices_a = len(coordsx_a)
#
#    if len(indices_owners) == 0:
#        indices_owners = range(len(length_volcells))
#
#    if dimPb == 3: n_smallfaces = 4
#    else: n_smallfaces = 2
#
#    list_bigface = []
#    for index_volume in indices_owners:
#        stride = cranges[index_volume]
#        stride_offset = numpy.cumsum(stride)
#        idx_vol_init = offset_volcells[index_volume]
#
#        idx_sides = []
#        for i in range(6):
#            if stride[i] == n_smallfaces:
#                idx_sides.append(i)
#        if idx_sides != []:
#            for idx_side in idx_sides:
#                indices_faces = []
#                start_idx = stride_offset[idx_side-1] if idx_side>0 else 0
#                for i in range(start_idx,stride_offset[idx_side]):
#                    face = EC_volcells[idx_vol_init+i]
#                    indices_faces.append(face)
#                conn_Nfaces = numpy.zeros((n_smallfaces,4), dtype=Internal.E_NpyInt)
#                for i,idx_face in enumerate(indices_faces):
#                    idx_face_init = offset_faces[idx_face-1]
#                    len_face = length_faces[idx_face-1]
#
#                    if len_face != 4:
#                        indices_face_5nodes = findFifthNode__(idx_face, offset_faces, length_faces, EC_faces, coordsx_a, coordsy_a, coordsz_a)
#                    k = 0
#                    for j in range(len_face):
#                        if len_face == 4 or (len_face != 4 and not j in indices_face_5nodes):
#                            conn_Nfaces[i][k] = EC_faces[idx_face_init+j]
#                            k+=1
#                if dimPb == 3: [point0, point1, point2, point3] = reorderNodesInCanonicalOrderForBigFace3D(conn_Nfaces)
#                else: [point0, point1, point2, point3] = reorderNodesInCanonicalOrderForBigFace2D(conn_Nfaces)
#                list_bigface.append([point0, point1, point2, point3])
#
#    len_new_faces = len(list_bigface)
#    flat_list_bigfaces = numpy.array([x for xs in list_bigface for x in xs])
#
#    zsnc_big_internal = Internal.newZone(name="Zone",zsize=[[nb_vertices_a, len_new_faces,0]], ztype="Unstructured")
#    gcBI = Internal.newGridCoordinates(parent=zsnc_big_internal)
#    Internal.newDataArray('CoordinateX', value=coordsx_a, parent=gcBI)
#    Internal.newDataArray('CoordinateY', value=coordsy_a, parent=gcBI)
#    Internal.newDataArray('CoordinateZ', value=coordsz_a, parent=gcBI)
#    Internal.newElements(name="GridElements", etype=7, econnectivity=flat_list_bigfaces, erange=[1, len_new_faces], eboundary=0, parent=zsnc_big_internal)
#    return zsnc_big_internal
#
#def createQuadSurfaceFromNgonPointListSmallFaceOrig__(a, PL):
#    a = Internal.getZones(a)[0] # get first zone
#    faces = Internal.getNGonNode(a)
#    vol_cells = Internal.getNFaceNode(a)
#
#    #EC_volcells = Internal.getNodeFromName(vol_cells, "ElementConnectivity")[1]
#    EC_faces = Internal.getNodeFromName(faces, "ElementConnectivity")[1]
#
#    #ER_volcells = Internal.getNodeFromName(vol_cells, "ElementRange")[1]
#    #ER_faces = Internal.getNodeFromName(faces, "ElementRange")[1]
#
#    offset_volcells = Internal.getNodeFromName(vol_cells, "ElementStartOffset")[1]
#    offset_faces = Internal.getNodeFromName(faces, "ElementStartOffset")[1]
#
#    length_volcells = offset_volcells[1:] - offset_volcells[:-1]
#    length_faces = offset_faces[1:] - offset_faces[:-1]
#    #lens_vol = numpy.unique(length_volcells)
#    #lens_faces = numpy.unique(length_faces)
#
#    coordsx_a = Internal.getNodeFromName(a, "CoordinateX")[1]
#    coordsy_a = Internal.getNodeFromName(a, "CoordinateY")[1]
#    coordsz_a = Internal.getNodeFromName(a, "CoordinateZ")[1]
#    nb_vertices_a = len(coordsx_a)
#
#    len_new_faces = len(PL)
#    len_EC_faces = len(PL)*4
#    new_EC_ngon_faces = numpy.zeros(len_EC_faces, dtype=Internal.E_NpyInt)
#    idx_new = 0
#    for idx in PL:
#        idx_face_init = offset_faces[idx-1]
#        len_face = length_faces[idx-1]
#        if len_face !=4:
#            indices_face_5nodes = findFifthNode__(idx, offset_faces, length_faces, EC_faces, coordsx_a, coordsy_a, coordsz_a)
#        for i in range(len_face):
#            if len_face == 4 or (len_face != 4 and not i in indices_face_5nodes):
#                new_EC_ngon_faces[idx_new] = EC_faces[idx_face_init+i]
#                idx_new += 1
#
#    zone = Internal.newZone(name="Zone",zsize=[[nb_vertices_a, len_new_faces, 0]], ztype="Unstructured")
#    gc = Internal.newGridCoordinates(parent=zone)
#    Internal.newDataArray('CoordinateX', value=coordsx_a, parent=gc)
#    Internal.newDataArray('CoordinateY', value=coordsy_a, parent=gc)
#    Internal.newDataArray('CoordinateZ', value=coordsz_a, parent=gc)
#    elt = Internal.newElements(name="GridElements", etype=7, econnectivity=new_EC_ngon_faces, erange=[1, len_new_faces], eboundary=0, parent=zone)
#    return zone
#
#def _createQuadConnectivityFromNgonPointListOrig__(a_hexa, a, PL, bcname, bctype):
#    a = Internal.getZones(a)[0] # get first zone
#    faces = Internal.getNGonNode(a)
#    vol_cells = Internal.getNFaceNode(a)
#
#    EC_volcells = Internal.getNodeFromName(vol_cells, "ElementConnectivity")[1]
#    EC_faces = Internal.getNodeFromName(faces, "ElementConnectivity")[1]
#
#    ER_volcells = Internal.getNodeFromName(vol_cells, "ElementRange")[1]
#    ER_faces = Internal.getNodeFromName(faces, "ElementRange")[1]
#
#    offset_volcells = Internal.getNodeFromName(vol_cells, "ElementStartOffset")[1]
#    offset_faces = Internal.getNodeFromName(faces, "ElementStartOffset")[1]
#
#    length_volcells = offset_volcells[1:] - offset_volcells[:-1]
#    length_faces = offset_faces[1:] - offset_faces[:-1]
#    lens_vol = numpy.unique(length_volcells)
#    lens_faces = numpy.unique(length_faces)
#
#    coordsx_a = Internal.getNodeFromName(a, "CoordinateX")[1]
#    coordsy_a = Internal.getNodeFromName(a, "CoordinateY")[1]
#    coordsz_a = Internal.getNodeFromName(a, "CoordinateZ")[1]
#    nb_vertices_a = len(coordsx_a)
#
#    len_new_faces = len(PL)
#    len_EC_faces = len(PL) * 4
#    new_EC_ngon_faces = numpy.zeros(len_EC_faces, dtype=Internal.E_NpyInt)
#    idx_new = 0
#    for idx in PL:
#        idx_face_init = offset_faces[idx - 1]
#        len_face = length_faces[idx - 1]
#        if len_face != 4:
#            indices_face_5nodes = findFifthNode__(idx, offset_faces, length_faces, EC_faces, coordsx_a, coordsy_a, coordsz_a)
#        for i in range(len_face):
#            if len_face == 4 or (len_face != 4 and not i in indices_face_5nodes):
#                new_EC_ngon_faces[idx_new] = EC_faces[idx_face_init + i]
#                idx_new += 1
#
#    elts = Internal.getNodesFromType(a_hexa, "Elements_t")
#    maxElt = Internal.getNodeFromName(elts[-1], "ElementRange")[1][1]
#    zones = Internal.getZones(a_hexa)
#    Internal.newElements(
#        name=bcname,
#        etype=7,
#        econnectivity=new_EC_ngon_faces,
#        erange=[maxElt + 1, maxElt + len_new_faces],
#        eboundary=0,
#        parent=zones[0]
#    )
#    C._addBC2Zone(zones[0], bcname, bctype, elementRange=[maxElt + 1, maxElt + len_new_faces])
#    return None
#
#
